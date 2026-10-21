"""Small PETSc-only collectives shared by the optional MPI hosts."""

import sys

import numpy as np

from .compiled_element import IncrementTooLarge
from .petsc_backend import _petsc


def gather(values, comm=None):
    """Concatenate rank-owned real arrays on every rank, including empty ones."""
    PETSc = _petsc()
    comm = PETSc.COMM_WORLD if comm is None else comm
    values = np.ascontiguousarray(values, dtype=PETSc.ScalarType).ravel()
    vec = PETSc.Vec().createMPI((len(values), PETSc.DECIDE), comm=comm)
    try:
        vec.getArray()[:] = values
        scatter, full = PETSc.Scatter.toAll(vec)
        try:
            scatter.scatter(vec, full, addv=PETSc.InsertMode.INSERT_VALUES,
                            mode=PETSc.ScatterMode.FORWARD)
            return full.getArray().copy()
        finally:
            scatter.destroy()
            full.destroy()
    finally:
        vec.destroy()


def collective_status(error, pnewdt=1.0, comm=None):
    """All ranks finish this vote before any raise or subsequent data gather."""
    PETSc = _petsc()
    comm = PETSc.COMM_WORLD if comm is None else comm
    code = 1 if isinstance(error, IncrementTooLarge) else 2 if error is not None else 0
    if error is not None:
        print("rank {}: {}: {}".format(comm.getRank(), type(error).__name__, error),
              file=sys.stderr, flush=True)
    pn = float(pnewdt) if np.isfinite(pnewdt) else 0.0
    status = gather([code, pn], comm).reshape(-1, 2)
    if np.any(status[:, 0] == 2):
        ranks = np.flatnonzero(status[:, 0] == 2).tolist()
        raise RuntimeError("MPI operation failed on ranks {}; see rank diagnostic".format(ranks)) from error
    if np.any(status[:, 0] == 1) or np.min(status[:, 1]) < 1.0:
        failure = IncrementTooLarge("an MPI partition requested an increment cutback")
        failure.pnewdt = float(np.min(status[:, 1]))
        raise failure from error
