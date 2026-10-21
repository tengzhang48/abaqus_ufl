"""Generate and call a standard Quad4 UEL through a packaged f2py wrapper."""

import importlib.util
from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile
import uuid

import numpy as np

from abaqus_ufl.core.fields import VectorField
from abaqus_ufl.core.verify import verify_before_generation
from abaqus_ufl.generators.uel_gen import generate_uel
from abaqus_ufl.generators.umat_gen import _nstate_per_gp, _state_var_info

_DRIVE_UEL = Path(__file__).with_name("drive_uel.f90")


class IncrementTooLarge(RuntimeError):
    """A UEL rejected a trial increment through PNEWDT."""


def _signature(weakform):
    if weakform.ndim != 2 or any(field.degree != 1 for field in weakform.fields.values()):
        raise ValueError("serial FE supports 2D Quad4 degree-one nodal fields only")
    if weakform.local_fields:
        raise ValueError("compiled FE excludes condensed local fields")
    return (type(weakform), type(weakform._mat), tuple(weakform._mat.props_names),
            tuple((name, type(field), field.degree) for name, field in weakform.fields.items()),
            tuple((name, entry['shape']) for name, entry in
                  _state_var_info(weakform._mat).items()))


def build_compiled_uel(weakform, workdir=None, verify=True):
    """Generate a standard Quad4 UEL and return its compiled f2py module.

    The wrapper ships in the wheel. No tests directory or external solver is
    needed. Without ``workdir``, the module owns a temporary build directory;
    a supplied directory retains the source and compiler log for inspection.
    Each build has a unique module name to prevent stale extension reuse.
    ``verify`` runs ``weakform.verify()`` before anything is built: True
    (default state), a dict of ``verify()`` arguments for a model-specific
    state, or False.
    """
    signature = _signature(weakform)
    verify_before_generation(weakform, verify)
    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required to build a compiled UEL")
    temporary = tempfile.TemporaryDirectory(prefix="abaqus_ufl_fe_") if workdir is None else None
    build = Path(temporary.name if temporary else workdir).resolve()
    build.mkdir(parents=True, exist_ok=True)
    source = build / "element.for"
    generate_uel(weakform, str(source), element="Quad4", formulation="standard",
                 verify=False)
    wrapper = _DRIVE_UEL
    module_name = "_abaqus_ufl_fe_" + uuid.uuid4().hex
    command = [sys.executable, "-m", "numpy.f2py", "-c", str(wrapper), str(source),
               "-m", module_name, "only:", "drive_uel", "drive_uel_batch", ":"]
    if np.lib.NumpyVersion(np.__version__) >= "1.26.0":
        command.extend(["--backend", "meson"])
    result = subprocess.run(command, cwd=str(build), capture_output=True,
                            text=True, encoding="utf-8", errors="replace",
                            env={**os.environ, "PYTHONUTF8": "1"}, check=False)
    log = result.stdout + "\n" + result.stderr
    (build / "compile.log").write_text(log, encoding="utf-8")
    if result.returncode:
        raise RuntimeError("f2py UEL build failed:\n" + log)
    from importlib.machinery import EXTENSION_SUFFIXES
    extension = next((build / (module_name + suffix) for suffix in EXTENSION_SUFFIXES
                      if (build / (module_name + suffix)).is_file()), None)
    if extension is None:
        raise RuntimeError("f2py did not produce the requested extension:\n" + log)
    spec = importlib.util.spec_from_file_location(module_name, str(extension))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module._abaqus_ufl_signature = signature
    module._abaqus_ufl_build_directory = temporary
    return module


class CompiledAbaqusElement:
    """Standard Quad4 adapter: R=-RHS and K=AMATRX=dR/dU.

    ``DU=U-U_previous`` carries nodal field history. The adapter uses a normal
    implicit residual-and-tangent request, rejects cutbacks/nonfinite output,
    and reads property overrides from the material instance. Material SVARS
    require ``n_elem``; evaluations are trials until explicitly committed.
    ``cache_batches=True`` reuses the last exact batch, including trial state.
    Use it only for kernels determined by their arguments without required
    external side effects. Scalar calls remain uncached.
    """

    def __init__(self, module, weakform, dt=1.0, time=None, *, n_elem=None,
                 cache_batches=False):
        if getattr(module, "_abaqus_ufl_signature", None) != _signature(weakform):
            raise ValueError("compile this weak-form type with build_compiled_uel first")
        self.module = module
        self.dof_per_node = sum(2 if isinstance(field, VectorField) else 1
                                for field in weakform.fields.values())
        self.props = np.array([getattr(weakform._mat, name)
                              for name in weakform._mat.props_names], dtype=float)
        if not np.isfinite(self.props).all():
            raise ValueError("material properties must be finite")
        self.ndofel = 4 * self.dof_per_node
        nstate = 4 * _nstate_per_gp(_state_var_info(weakform._mat))
        self.svars_size = max(1, nstate)
        if n_elem is not None and (isinstance(n_elem, bool)
                or not isinstance(n_elem, (int, np.integer)) or n_elem < 0):
            raise ValueError("n_elem must be a nonnegative integer")
        if nstate and n_elem is None:
            raise ValueError("stateful elements require n_elem for committed SVARS")
        self.svars = np.zeros((n_elem, nstate)) if nstate else None
        self.elem_ids = None
        self._jprops = np.array([0], dtype=np.int32)
        self._lflags = np.array([1, 0, 1, 0, 0, 0], dtype=np.int32)
        self._params = np.zeros(3)
        self.cache_batches = bool(cache_batches)
        self._batch_cache = None
        self.batch_evaluations = 0
        self.batch_cache_hits = 0
        self.dt = float(dt)
        self.begin_increment(1, time_begin=0.0 if time is None else time)

    def begin_increment(self, kinc, time_begin=None, dt=None):
        """Set beginning-of-increment TIME and KINC for one analysis step.

        For uniform increments, time_begin defaults to (kinc-1)*dt. Supply it
        explicitly for nonuniform increments. This runtime uses KSTEP=1.
        """
        if isinstance(kinc, bool) or not isinstance(kinc, (int, np.integer)) or kinc < 1:
            raise ValueError("kinc must be a positive integer")
        increment = self.dt if dt is None else float(dt)
        beginning = (kinc - 1) * increment if time_begin is None else float(time_begin)
        if not np.isfinite([beginning, increment]).all() or increment <= 0 or beginning < 0:
            raise ValueError("time_begin must be nonnegative and dt positive, both finite")
        self.time, self.dt, self.kinc = beginning, increment, int(kinc)
        self.clear_batch_cache()

    def clear_batch_cache(self):
        """Discard the saved trial; evaluation/hit counters are cumulative."""
        self._batch_cache = None

    def reject_step(self):
        """Discard a rejected trial without modifying committed SVARS."""
        self.clear_batch_cache()

    def _batch_ids(self, nelem):
        if self.elem_ids is None:
            return np.arange(1, nelem + 1, dtype=np.int32)
        ids = np.asarray(self.elem_ids)
        if (ids.shape != (nelem,) or not np.issubdtype(ids.dtype, np.integer)
                or np.any(ids < 1) or np.any(ids > np.iinfo(np.int32).max)):
            raise ValueError("elem_ids must contain one positive int32 global ID per element")
        return ids.astype(np.int32, copy=False)

    def _batch_key(self, coords, U, DU, sv):
        # Exact pre-call snapshots include shapes, dtypes and bit patterns.
        arrays = (coords, U, DU, sv, self.props, self._jprops,
                  np.array([self.time, self.time]), self._lflags,
                  self._params, self._batch_ids(len(U)), np.array([self.dt]))
        return (self.kinc, tuple((a.shape, a.dtype.str, a.tobytes())
                                for a in map(np.asarray, arrays)))

    def _check_batch_error(self, error):
        """MPI overrides vote on preparation errors before any rank raises."""
        if error is not None:
            raise error

    def _agree_batch_cache_hit(self, hit, error):
        """MPI overrides must enter a collective even for disabled caches."""
        self._check_batch_error(error)
        return hit

    def _evaluate_batch(self, coords, U, DU):
        key, error, hit = None, None, False
        entry = self._batch_cache
        try:
            # Preparation can fail on only one MPI rank. Include it in the
            # same error vote as key creation, before any cache/evaluation vote.
            coords, U, DU = (np.asarray(a, dtype=float) for a in (coords, U, DU))
            if U.ndim != 2:
                raise ValueError("compiled Quad4 batch U must be a two-dimensional array")
            ne = len(U)
            sv = self.svars.copy() if self.svars is not None else np.zeros((ne, 1))
            if coords.shape != (ne, 4, 2) or U.shape != (ne, self.ndofel) or DU.shape != U.shape:
                raise ValueError("compiled Quad4 batch input shapes disagree with the nodal field layout")
            if sv.shape != (ne, self.svars_size):
                raise ValueError("committed SVARS shape disagrees with the element batch")
            if not all(np.isfinite(a).all() for a in (coords, U, DU, sv)):
                raise ValueError("element inputs must be finite")
            if self.cache_batches:
                key = self._batch_key(coords, U, DU, sv)
                hit = entry is not None and entry[0] is self.module and entry[1] == key
        except Exception as exc:
            error = exc
        try:
            if self._agree_batch_cache_hit(hit, error):
                self.batch_cache_hits += 1
                return (*(a.copy() for a in entry[2][:3]), entry[2][3])
            self.clear_batch_cache()
            self.batch_evaluations += 1
            result = self._call_batch(coords, U, DU, sv)
            self._check_result(result)
            ne = len(U)
            if (np.shape(result[0]) != (self.ndofel, ne)
                    or np.shape(result[1]) != (self.ndofel, self.ndofel, ne)
                    or np.shape(result[2]) != (self.svars_size, ne)):
                raise RuntimeError("compiled UEL returned invalid batch output shapes")
            if self.cache_batches:
                saved = tuple(np.array(a, copy=True, order="K") for a in result[:3])
                for a in saved:
                    a.flags.writeable = False
                self._batch_cache = (self.module, key, (*saved, float(result[3])))
            return result
        except Exception:
            self.clear_batch_cache()
            raise

    @staticmethod
    def _check_result(result):
        rhs, tangent, state, pnewdt = result
        if not np.isfinite(pnewdt) or pnewdt < 1.0:
            error = IncrementTooLarge("generated UEL requested a cutback (PNEWDT={})".format(pnewdt))
            error.pnewdt = float(pnewdt)
            raise error
        if not all(np.isfinite(value).all() for value in (rhs, tangent, state)):
            raise RuntimeError("generated UEL returned nonfinite values")

    def _call(self, ei, coords, U_e, DU_e, sv):
        jelem = ei + 1 if self.elem_ids is None else int(self.elem_ids[ei])
        return self.module.drive_uel(
            sv, np.asfortranarray(coords.T), U_e, DU_e, self.props,
            self._jprops, np.array([self.time, self.time]), self.dt,
            1.0, self._lflags, self._params, 1, 1, self.kinc, jelem, 1.0)

    def element_rk(self, ei, coords, U_e, DU_e):
        coords = np.asarray(coords, dtype=float)
        U_e, DU_e = np.asarray(U_e, dtype=float), np.asarray(DU_e, dtype=float)
        ndof = 4 * self.dof_per_node
        if coords.shape != (4, 2) or U_e.shape != (ndof,) or DU_e.shape != (ndof,):
            raise ValueError("compiled Quad4 input shapes disagree with the nodal field layout")
        if not all(np.isfinite(value).all() for value in (coords, U_e, DU_e)):
            raise ValueError("element inputs must be finite")
        sv = self.svars[ei].copy() if self.svars is not None else np.zeros(1)
        try:
            result = self._call(ei, coords, U_e, DU_e, sv)
            self._check_result(result)
        except Exception:
            self.clear_batch_cache()
            raise
        return -np.asarray(result[0], dtype=float), np.ascontiguousarray(result[1], dtype=float)

    def _call_batch(self, coords, U, DU, sv):
        ne = len(U)
        ids = self._batch_ids(ne)
        if not ne:
            return np.empty((self.ndofel, 0)), np.empty((self.ndofel, self.ndofel, 0)), sv.T.copy(), 1.0
        return self.module.drive_uel_batch(
            np.asfortranarray(sv.T), np.asfortranarray(coords.transpose(2, 1, 0)),
            np.asfortranarray(U.T), np.asfortranarray(DU.T), self.props,
            self._jprops, np.array([self.time, self.time]), self.dt, 1.0,
            self._lflags, self._params, 1, 1, self.kinc, 1.0, ids)

    def element_rk_batch(self, coords_all, U_all, DU_all):
        """Trial R/K in element-major order; committed state remains unchanged."""
        rhs, tangent, _, _ = self._evaluate_batch(coords_all, U_all, DU_all)
        return -np.asarray(rhs, float).T, np.ascontiguousarray(tangent.transpose(2, 0, 1))

    def commit_step(self, nodes, elems, dof_per_node, U, U_previous):
        """Commit the accepted trial from old state, reusing an exact cached batch."""
        if self.svars is None:
            self.clear_batch_cache()
            return
        try:
            error, inputs = None, None
            try:
                if dof_per_node != self.dof_per_node:
                    raise ValueError("commit DOF layout disagrees with the compiled element")
                elems, nodes = np.asarray(elems), np.asarray(nodes, float)
                gm = (elems[:, :, None] * dof_per_node + np.arange(dof_per_node)).reshape(len(elems), -1)
                inputs = nodes[elems], np.asarray(U)[gm], (np.asarray(U) - U_previous)[gm]
            except Exception as exc:
                error = exc
            self._check_batch_error(error)
            self.commit_group(*inputs)
        finally:
            self.clear_batch_cache()

    def commit_group(self, coords, U, DU):
        """Commit one adapter's full batch; never commit a rejected evaluation."""
        try:
            if self.svars is not None:
                _, _, state, _ = self._evaluate_batch(coords, U, DU)
                self.svars = np.ascontiguousarray(state.T, dtype=float)
        finally:
            self.clear_batch_cache()

    def commit_element(self, ei, coords, U, DU):
        """Commit one local element, for partitioned hosts."""
        try:
            if self.svars is not None:
                args = [np.asarray(a, float) for a in (coords, U, DU)]
                result = self._call(ei, *args, self.svars[ei].copy())
                self._check_result(result)
                self.svars[ei] = result[2]
        finally:
            self.clear_batch_cache()
