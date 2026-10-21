"""Partition compiled element batches while retaining a replicated serial host.

Every rank holds the mesh and old state, gathers R/K/trial state in original
element order, and follows the same Newton/line-search decisions. A root-only
linear solve broadcasts its update. This is distinct from the distributed
matrix and ghost-vector host in :mod:`petsc_mpi`. PETSc is loaded lazily;
mpi4py is not required by either host.
"""

import importlib.util
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
import time

import numpy as np

from .compiled_element import CompiledAbaqusElement, _signature, build_compiled_uel
from ._collective import collective_status, gather
from .petsc_backend import _petsc, make_petsc_solver


def build_module(weakform, build_root, verify=True, comm=None):
    """Build once on rank zero and load the same extension from a shared folder.

    ``build_root`` must be visible to every rank. The public weak-form builder
    supplies a unique, neutral module name; no stale build cache is reused.
    """
    PETSc = _petsc()
    comm = PETSc.COMM_WORLD if comm is None else comm
    module, error = None, None
    if comm.getRank() == 0:
        try:
            module = build_compiled_uel(weakform, workdir=build_root, verify=verify)
        except Exception as exc:
            error = exc
    collective_status(error, comm=comm)
    name = bytes(gather(list(module.__name__.encode('ascii')) if module is not None else [], comm)
                 .astype(np.uint8)).decode('ascii')
    try:
        if module is None:
            folder = Path(build_root).resolve()
            extension = next(folder / (name + suffix) for suffix in EXTENSION_SUFFIXES
                             if (folder / (name + suffix)).is_file())
            spec = importlib.util.spec_from_file_location(name, str(extension))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module._abaqus_ufl_build_directory = None
        module._abaqus_ufl_signature = _signature(weakform)
    except Exception as exc:
        error = exc
    collective_status(error, comm=comm)
    return module


class MPIBatchedElement(CompiledAbaqusElement):
    """Replicated full-batch adapter with unanimous cache-hit decisions."""

    def __init__(self, *args, comm=None, **kwargs):
        super().__init__(*args, **kwargs)
        PETSc = _petsc()
        self.comm = PETSc.COMM_WORLD if comm is None else comm
        self.rank, self.size = self.comm.getRank(), self.comm.getSize()
        self.batch_calls = 0
        self.batch_elapsed_s = 0.0

    def _agree_batch_cache_hit(self, hit, error):
        collective_status(error, comm=self.comm)
        return bool(np.all(gather([int(hit)], self.comm) == 1))

    def _call_batch(self, coords_all, U_all, DU_all, sv_all):
        started = time.monotonic()
        ne = len(U_all)
        counts = gather([ne], self.comm)
        if np.any(counts != ne):
            raise ValueError("MPI batch element counts must agree on all ranks")
        begin, end = ne * self.rank // self.size, ne * (self.rank + 1) // self.size
        owned = slice(begin, end)
        previous_ids = self.elem_ids
        error, pn = None, 1.0
        packet = np.empty((end - begin, self.ndofel + self.ndofel**2 + self.svars_size))
        try:
            ids = self._batch_ids(ne)
            if (coords_all.shape != (ne, 4, 2) or U_all.shape != (ne, self.ndofel)
                    or DU_all.shape != U_all.shape or sv_all.shape != (ne, self.svars_size)):
                raise ValueError("MPI batch input shapes disagree with the Quad4 layout")
            if not all(np.isfinite(a).all() for a in (coords_all, U_all, DU_all, sv_all)):
                raise ValueError("MPI batch inputs must be finite")
            self.elem_ids = ids[owned]
            if end > begin:
                result = super()._call_batch(coords_all[owned], U_all[owned], DU_all[owned], sv_all[owned])
                self._check_result(result)
                rhs, k, sv, pn = result
                packet[:] = np.column_stack((rhs.T, k.transpose(2, 0, 1).reshape(end - begin, -1), sv.T))
        except Exception as exc:
            error = exc
            pn = getattr(exc, 'pnewdt', 1.0)
        finally:
            self.elem_ids = previous_ids
        collective_status(error, pn, self.comm)
        data = gather(packet, self.comm).reshape(ne, packet.shape[1])
        n = self.ndofel
        self.batch_calls += 1
        self.batch_elapsed_s += time.monotonic() - started
        return (data[:, :n].T, data[:, n:n + n*n].reshape(ne, n, n).transpose(1, 2, 0),
                data[:, n + n*n:].T, 1.0)

    def root_linear_solver(self, **options):
        """Return a collective CSR solver that solves only on rank zero."""
        PETSc = _petsc()
        solve, error = None, None
        if self.rank == 0:
            try:
                solve = make_petsc_solver(comm=PETSc.COMM_SELF, **options)
            except Exception as exc:
                error = exc
        collective_status(error, comm=self.comm)

        def linear(K, b):
            update, error = np.empty(0), None
            if self.rank == 0:
                try:
                    update = solve(K, b)
                except Exception as exc:
                    error = exc
            collective_status(error, comm=self.comm)
            return gather(update, self.comm)
        return linear
