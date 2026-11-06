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


class NonfiniteElementOutput(RuntimeError):
    """An element returned a nonfinite residual, tangent or state.

    Solvers that backtrack treat such a trial as an infinite residual; any
    other exception from an element is not a numerical rejection and
    propagates. It subclasses RuntimeError for backward compatibility.
    """


def _signature(weakform):
    if weakform.ndim != 2 or any(field.degree != 1 for field in weakform.fields.values()):
        raise ValueError("serial FE supports 2D Quad4 degree-one nodal fields only")
    if weakform.local_fields:
        raise ValueError("compiled FE excludes condensed local fields")
    return (type(weakform), type(weakform._mat), tuple(weakform._mat.props_names),
            tuple((name, type(field), field.degree) for name, field in weakform.fields.items()),
            tuple((name, entry['shape']) for name, entry in
                  _state_var_info(weakform._mat).items()))


def _build_directory(workdir, prefix):
    """A retained ``workdir`` or a temporary directory owned by the module."""
    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required to build a compiled UEL")
    temporary = tempfile.TemporaryDirectory(prefix=prefix) if workdir is None else None
    build = Path(temporary.name if temporary else workdir).resolve()
    build.mkdir(parents=True, exist_ok=True)
    return temporary, build


def _link_with_wrapper(build, inputs, module_name, temporary, prior_log=""):
    """Link Fortran sources or objects with the packaged wrapper and import it."""
    command = [sys.executable, "-m", "numpy.f2py", "-c", str(_DRIVE_UEL),
               *(str(item) for item in inputs),
               "-m", module_name, "only:", "drive_uel", "drive_uel_batch", ":"]
    if np.lib.NumpyVersion(np.__version__) >= "1.26.0":
        command.extend(["--backend", "meson"])
    result = subprocess.run(command, cwd=str(build), capture_output=True,
                            text=True, encoding="utf-8", errors="replace",
                            env={**os.environ, "PYTHONUTF8": "1"}, check=False)
    log = prior_log + result.stdout + "\n" + result.stderr
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
    module._abaqus_ufl_build_directory = temporary
    return module


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
    temporary, build = _build_directory(workdir, "abaqus_ufl_fe_")
    source = build / "element.for"
    generate_uel(weakform, str(source), element="Quad4", formulation="standard",
                 verify=False)
    module = _link_with_wrapper(build, [source], "_abaqus_ufl_fe_" + uuid.uuid4().hex,
                                temporary)
    module._abaqus_ufl_signature = signature
    return module


def build_external_uel(sources, workdir=None, *, compile_flags=(), module_name=None):
    """Compile a standard Quad4 UEL written outside abaqus_ufl; return its module.

    ``sources`` is one Fortran file or a sequence of them, compiled in the
    given order by ``gfortran -c -O2 -fPIC -J <build directory>`` followed by
    ``compile_flags`` (later flags take precedence, so ``("-O3",)`` or
    ``("-ffixed-line-length-none",)`` adjust the defaults).  The objects are
    linked with the packaged ``drive_uel.f90`` through numpy.f2py, so the
    subroutine must be named ``UEL`` and have the standard Abaqus argument
    list.  ``workdir`` behaves as in :func:`build_compiled_uel`; the
    compiler log is ``compile.log``.  ``module_name`` is a readable prefix:
    every build appends a unique suffix, so a rebuild (even into the same
    ``workdir``) never reuses native code already loaded under an earlier
    name or path, and never replaces another entry of ``sys.modules``.  Use
    the module with :meth:`CompiledAbaqusElement.from_module`.
    """
    if isinstance(sources, (str, os.PathLike)):
        sources = [sources]
    sources = [Path(item).resolve() for item in sources]
    if not sources:
        raise ValueError("sources must name at least one Fortran file")
    missing = [str(item) for item in sources if not item.is_file()]
    if missing:
        raise FileNotFoundError("UEL source files not found: " + ", ".join(missing))
    if isinstance(compile_flags, str) or not all(isinstance(flag, str) for flag in compile_flags):
        raise ValueError("compile_flags must be a sequence of strings")
    if module_name is None:
        module_name = "_abaqus_ufl_ext"
    elif not (isinstance(module_name, str) and module_name.isidentifier()):
        raise ValueError("module_name must be a valid Python identifier")
    module_name = module_name + "_" + uuid.uuid4().hex     # unique per build
    temporary, build = _build_directory(workdir, "abaqus_ufl_ext_")
    objects, log = [], ""
    for index, source in enumerate(sources):
        target = build / "external_{}.o".format(index)
        command = ["gfortran", "-c", "-O2", "-fPIC", "-J", str(build), *compile_flags,
                   str(source), "-o", str(target)]
        result = subprocess.run(command, cwd=str(source.parent), capture_output=True,
                                text=True, encoding="utf-8", errors="replace", check=False)
        log += " ".join(command) + "\n" + result.stdout + result.stderr + "\n"
        if result.returncode:
            (build / "compile.log").write_text(log, encoding="utf-8")
            raise RuntimeError("gfortran failed on {}:\n{}".format(source, result.stderr))
        objects.append(target)
    return _link_with_wrapper(build, objects, module_name, temporary, prior_log=log)


def _is_count(value):
    """A nonnegative Python or NumPy integer, not a bool."""
    return (not isinstance(value, (bool, np.bool_))
            and isinstance(value, (int, np.integer)) and value >= 0)


def _int32_vector(values, name):
    values = np.asarray(values)
    if values.size == 0 and values.ndim == 1:
        return np.zeros(0, dtype=np.int32)
    if (values.ndim != 1 or not np.issubdtype(values.dtype, np.integer)
            or values.min() < np.iinfo(np.int32).min or values.max() > np.iinfo(np.int32).max):
        raise ValueError(name + " must be a one-dimensional sequence of int32 integers")
    return values.astype(np.int32)


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
        self._setup(
            module,
            dof_per_node=sum(2 if isinstance(field, VectorField) else 1
                             for field in weakform.fields.values()),
            props=np.array([getattr(weakform._mat, name)
                            for name in weakform._mat.props_names], dtype=float),
            nstate=4 * _nstate_per_gp(_state_var_info(weakform._mat)),
            n_elem=n_elem, svars_initial=None,
            jprops=np.array([0], dtype=np.int32),
            lflags=np.array([1, 0, 1, 0, 0, 0], dtype=np.int32),
            dt=dt, time=time, cache_batches=cache_batches)

    @classmethod
    def from_module(cls, module, *, dof_per_node, props, n_elem=None, nsvars=0,
                    svars_initial=None, jprops=(), lflags=(1, 0, 1, 0, 0, 0),
                    dt=1.0, time=None, cache_batches=False):
        """Adapter for a standard Quad4 UEL compiled outside abaqus_ufl.

        ``module`` must expose ``drive_uel`` and ``drive_uel_batch`` of the
        packaged ``drive_uel.f90``; :func:`build_external_uel` builds one.
        ``dof_per_node`` gives NDOFEL = 4 * dof_per_node in interleaved
        node/component order.  ``props`` (PROPS) and ``jprops`` (JPROPS,
        int32) are passed unchanged; ``lflags`` are six integers (the default
        requests residual and tangent).  ``nsvars`` is NSVARS per element:
        with ``nsvars > 0`` committed SVARS need ``n_elem`` and start from
        ``svars_initial`` (zeros by default) in every element; with 0 the
        element is stateless.  ``dt``, ``time`` and ``cache_batches`` are as
        in the constructor.  The adapter passes beginning-of-increment TIME;
        a subclass may override :meth:`begin_increment` for another
        convention, and ``from_module`` then returns that subclass.
        """
        for name in ("drive_uel", "drive_uel_batch"):
            if not callable(getattr(module, name, None)):
                raise ValueError("module must expose drive_uel and drive_uel_batch "
                                 "(build it with build_external_uel)")
        if not _is_count(dof_per_node) or dof_per_node < 1:
            raise ValueError("dof_per_node must be a positive integer")
        props = np.array(props, dtype=float)
        if props.ndim != 1 or not np.isfinite(props).all():
            raise ValueError("props must be a finite one-dimensional sequence")
        if not _is_count(nsvars):
            raise ValueError("nsvars must be a nonnegative integer")
        if svars_initial is not None:
            svars_initial = np.array(svars_initial, dtype=float)
            if svars_initial.shape != (nsvars,) or not np.isfinite(svars_initial).all():
                raise ValueError("svars_initial must hold nsvars finite values")
        jprops = _int32_vector(jprops, "jprops")
        lflags = _int32_vector(lflags, "lflags")
        if lflags.shape != (6,):
            raise ValueError("lflags must hold six integers")
        if (isinstance(dt, (bool, np.bool_)) or not isinstance(dt, (int, float, np.number))
                or not np.isfinite(dt) or dt <= 0):
            raise ValueError("dt must be positive and finite")
        element = object.__new__(cls)
        element._setup(module, dof_per_node=int(dof_per_node), props=props,
                       nstate=int(nsvars), n_elem=n_elem, svars_initial=svars_initial,
                       jprops=jprops, lflags=lflags, dt=dt, time=time,
                       cache_batches=cache_batches)
        return element

    def _setup(self, module, *, dof_per_node, props, nstate, n_elem, svars_initial,
               jprops, lflags, dt, time, cache_batches):
        """Attribute setup shared by the constructor and :meth:`from_module`."""
        self.module = module
        self.dof_per_node = dof_per_node
        self.props = props
        if not np.isfinite(self.props).all():
            raise ValueError("material properties must be finite")
        self.ndofel = 4 * self.dof_per_node
        self.svars_size = max(1, nstate)
        if n_elem is not None and not _is_count(n_elem):
            raise ValueError("n_elem must be a nonnegative integer")
        if nstate and n_elem is None:
            raise ValueError("stateful elements require n_elem for committed SVARS")
        if not nstate:
            self.svars = None
        elif svars_initial is None:
            self.svars = np.zeros((n_elem, nstate))
        else:
            self.svars = np.tile(svars_initial, (n_elem, 1))
        self.elem_ids = None
        self._jprops = jprops
        self._lflags = lflags
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
            raise NonfiniteElementOutput("generated UEL returned nonfinite values")

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
