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
from abaqus_ufl.generators.uel_gen import generate_uel


def _signature(weakform):
    if weakform.ndim != 2 or any(field.degree != 1 for field in weakform.fields.values()):
        raise ValueError("serial FE supports 2D Quad4 degree-one nodal fields only")
    if weakform.local_fields or getattr(weakform._mat, "state_vars", {}):
        raise ValueError("serial FE excludes material state variables and local fields")
    return (type(weakform), type(weakform._mat), tuple(weakform._mat.props_names),
            tuple((name, type(field), field.degree) for name, field in weakform.fields.items()))


def build_compiled_uel(weakform, workdir=None):
    """Generate a standard Quad4 UEL and return its compiled f2py module.

    The wrapper ships in the wheel. No tests directory or external solver is
    needed. Without ``workdir``, the module owns a temporary build directory;
    a supplied directory retains the source and compiler log for inspection.
    Each build has a unique module name to prevent stale extension reuse.
    """
    signature = _signature(weakform)
    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required to build a compiled UEL")
    temporary = tempfile.TemporaryDirectory(prefix="abaqus_ufl_fe_") if workdir is None else None
    build = Path(temporary.name if temporary else workdir).resolve()
    build.mkdir(parents=True, exist_ok=True)
    source = build / "element.for"
    generate_uel(weakform, str(source), element="Quad4", formulation="standard")
    wrapper = Path(__file__).with_name("drive_uel.f90")
    module_name = "_abaqus_ufl_fe_" + uuid.uuid4().hex
    command = [sys.executable, "-m", "numpy.f2py", "-c", str(wrapper), str(source),
               "-m", module_name, "only:", "drive_uel", ":"]
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
    """Stateless standard Quad4 adapter: R=-RHS and K=AMATRX=dR/dU.

    ``DU=U-U_previous`` carries nodal field history. The adapter uses a normal
    implicit residual-and-tangent request, rejects cutbacks/nonfinite output,
    and reads property overrides from the material instance.
    """

    def __init__(self, module, weakform, dt=1.0, time=None):
        if getattr(module, "_abaqus_ufl_signature", None) != _signature(weakform):
            raise ValueError("compile this weak-form type with build_compiled_uel first")
        self.module = module
        self.dof_per_node = sum(2 if isinstance(field, VectorField) else 1
                                for field in weakform.fields.values())
        self.props = np.array([getattr(weakform._mat, name)
                              for name in weakform._mat.props_names], dtype=float)
        if not np.isfinite(self.props).all():
            raise ValueError("material properties must be finite")
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

    def element_rk(self, ei, coords, U_e, DU_e):
        coords = np.asarray(coords, dtype=float)
        U_e, DU_e = np.asarray(U_e, dtype=float), np.asarray(DU_e, dtype=float)
        ndof = 4 * self.dof_per_node
        if coords.shape != (4, 2) or U_e.shape != (ndof,) or DU_e.shape != (ndof,):
            raise ValueError("compiled Quad4 input shapes disagree with the nodal field layout")
        if not all(np.isfinite(value).all() for value in (coords, U_e, DU_e)):
            raise ValueError("element inputs must be finite")
        rhs, tangent, state, pnewdt = self.module.drive_uel(
            np.zeros(1), np.asfortranarray(coords.T), U_e, DU_e, self.props,
            np.array([0], dtype=np.int32), np.array([self.time, self.time]), self.dt,
            1.0, np.array([1, 0, 1, 0, 0, 0], dtype=np.int32), np.zeros(3),
            1, 1, self.kinc, ei + 1, 1.0)
        if not np.isfinite(pnewdt) or pnewdt < 1.0:
            raise RuntimeError("generated UEL requested a cutback (PNEWDT={})".format(pnewdt))
        if not all(np.isfinite(value).all() for value in (rhs, tangent, state)):
            raise RuntimeError("generated UEL returned nonfinite values")
        return -np.asarray(rhs, dtype=float), np.asarray(tangent, dtype=float)
