"""Regenerate, compile, and directly execute the generated coupled UEL.

Gates: deterministic regeneration against the committed source, gfortran
compile, and direct f2py element calls compared against the Python
reference assembly (which check_assembled.py has already gated against
closed forms and finite differences):

- RHS and AMATRX parity at the reference state, at the uniform-T
  thermal-stress state, and at a generic deformed/heated state with a
  nonzero increment;
- an independent finite-difference check of the COMPILED tangent,
  AMATRX = -dRHS/dU, at the generic state.
"""

from pathlib import Path
import importlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from abaqus_ufl.core.reference_assembly import assemble_element
from abaqus_ufl.generators.uel_gen import generate_uel

from check_assembled import (
    T_DOFS,
    distorted_quad,
    props_array,
    set_T,
    unit_square,
)


def run_checked(command, cwd):
    result = subprocess.run(
        [str(item) for item in command],
        cwd=str(cwd),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Command failed: {}\nstdout:\n{}\nstderr:\n{}".format(
                " ".join(str(item) for item in command),
                result.stdout,
                result.stderr,
            )
        )


def call_uel(module, coords, u, du, dtime=0.1, lflags3=1, lflags1=72):
    svars = np.zeros(1)
    jprops = np.array([0], dtype=np.int32)
    time = np.array([0.0, 0.0])
    # LFLAGS(1)=72 transient coupled thermal-stress procedure (the
    # procedure the demo decks request), LFLAGS(2)=1 nlgeom, LFLAGS(3)
    # per request type (1 = normal residual+stiffness evaluation).
    lflags = np.array([lflags1, 1, lflags3, 0, 0, 0], dtype=np.int32)
    params = np.zeros(3)
    rhs, amatrx, svars_out, pnewdt = module.drive_uel(
        svars, np.asarray(coords, dtype=float), np.asarray(u, dtype=float),
        np.asarray(du, dtype=float), props_array(), jprops, time, dtime,
        1.0, lflags, params, 1, 1, 1, 1, 0.0,
    )
    return rhs, amatrx, pnewdt


def check():
    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required for the compiled pipeline gate")

    committed = HERE / "scalar_diffusion_uel.for"
    if not committed.exists():
        raise RuntimeError("Missing committed source; run `python build.py`")

    from build import HeatDiffusionProblem

    with tempfile.TemporaryDirectory(prefix="abaqus_ufl_sdiff_") as tmp:
        work = Path(tmp)
        generated = work / committed.name
        generate_uel(HeatDiffusionProblem(), str(generated),
                     element="Quad4", formulation="standard")

        if generated.read_bytes() != committed.read_bytes():
            raise AssertionError(
                "Committed Fortran is stale; regenerate with `python build.py`"
            )
        print("[PASS] deterministic generated source")

        run_checked(
            [
                "gfortran",
                "-c",
                "-ffixed-form",
                "-ffixed-line-length-none",
                generated,
                "-o",
                work / "scalar_diffusion_uel.o",
            ],
            work,
        )
        print("[PASS] gfortran compile")

        module_name = "_abaqus_ufl_scalar_diffusion_uel"
        run_checked(
            [
                sys.executable,
                "-m",
                "numpy.f2py",
                "-c",
                HERE / "f2py" / "drive_uel.f90",
                generated,
                "-m",
                module_name,
                "only:",
                "drive_uel",
                ":",
            ],
            work,
        )

        sys.path.insert(0, str(work))
        try:
            sys.modules.pop(module_name, None)
            module = importlib.import_module(module_name)
        finally:
            sys.path.pop(0)

        problem = HeatDiffusionProblem()
        generic_U = np.array([
            0.013, -0.008, 1.2, 0.021, 0.004, 0.7, -0.011, 0.017, -0.4,
            0.006, -0.015, 2.1,
        ])
        states = (
            ("reference", unit_square(), np.zeros(12), np.zeros(12)),
            ("uniform-T thermal stress", unit_square(),
             set_T(np.zeros(12), [3.0] * 4), np.zeros(12)),
            ("generic coupled", distorted_quad(), generic_U, 0.3 * generic_U),
        )
        rhs_errors, amx_errors = [], []
        for name, coords, U, DU in states:
            rhs_f, amx_f, pnewdt = call_uel(module, coords, U, DU)
            if not (np.isfinite(rhs_f).all() and np.isfinite(amx_f).all()):
                raise AssertionError("non-finite compiled output at " + name)
            rhs_p, amx_p = assemble_element(
                problem, coords, U, DU, 0.1, props_array(), element="quad4")
            scale_r = max(1.0e-12, float(np.max(np.abs(rhs_p))))
            scale_a = float(np.max(np.abs(amx_p)))
            err_r = float(np.max(np.abs(rhs_f - rhs_p)))
            err_a = float(np.max(np.abs(amx_f - amx_p))) / scale_a
            if err_r > 1e-9 * max(1.0, scale_r):
                raise AssertionError(
                    "RHS parity failed at {}: max abs err {}".format(
                        name, err_r))
            if err_a > 1e-9:
                raise AssertionError(
                    "AMATRX parity failed at {}: rel err {}".format(
                        name, err_a))
            if pnewdt <= 0.0:
                raise AssertionError("invalid PNEWDT at " + name)
            rhs_errors.append(err_r)
            amx_errors.append(err_a)
            print("[PASS] f2py RHS/AMATRX parity at {} state".format(name))

        # Independent FD check of the compiled tangent at the generic state.
        coords, U, DU = distorted_quad(), generic_U, 0.3 * generic_U
        rhs0, amx0, _ = call_uel(module, coords, U, DU)
        eps = 1e-7
        fd = np.zeros((12, 12))
        for j in range(12):
            Up = U.copy()
            Up[j] += eps
            DUp = DU.copy()
            DUp[j] += eps
            rhs_p, _, _ = call_uel(module, coords, Up, DUp)
            fd[:, j] = (rhs_p - rhs0) / eps
        rel = np.max(np.abs(amx0 - (-fd))) / (np.max(np.abs(amx0)) + 1e-30)
        if not np.isfinite(rel):
            raise AssertionError("non-finite tangent comparison")
        if rel > 5e-5:
            raise AssertionError(
                "compiled AMATRX vs -dRHS/dU failed: rel err {:.3e}".format(
                    rel))
        print("[PASS] compiled AMATRX vs -dRHS/dU (rel err {:.1e})".format(
            rel))

        # Contract: only request types 1, 2, 5 are supported; everything
        # else (mass, damping, initial acceleration, or any unknown code)
        # must return ZEROED arrays, not the static residual and stiffness.
        for req in (3, 4, 6, 100):
            rhs_m, amx_m, _ = call_uel(module, coords, U, DU, lflags3=req)
            if np.max(np.abs(rhs_m)) != 0.0 or np.max(np.abs(amx_m)) != 0.0:
                raise AssertionError(
                    "LFLAGS(3)={} returned nonzero arrays".format(req))
        rhs_p, amx_p2, _ = call_uel(module, coords, U, DU, lflags1=99)
        if np.max(np.abs(rhs_p)) != 0.0 or np.max(np.abs(amx_p2)) != 0.0:
            raise AssertionError(
                "unsupported LFLAGS(1)=99 returned nonzero arrays")
        print("[PASS] unsupported request types and procedures return "
              "zeroed arrays")

        # Contract: an invalid deformation state (det F <= 0) must request a
        # cutback instead of returning finite but meaningless matrices.
        U_bad = np.array(U, dtype=float).copy()
        for node in range(4):
            U_bad[node * 3] = -2.0 * coords[0][node]
        _, _, pnewdt_bad = call_uel(module, coords, U_bad, DU)
        if not pnewdt_bad < 1.0:
            raise AssertionError(
                "det(F)<=0 state did not request a cutback "
                "(PNEWDT={})".format(pnewdt_bad))
        print("[PASS] nonpositive det(F) requests a cutback")

        return {
            "max_rhs_abs_error": max(rhs_errors),
            "max_amatrx_rel_error": max(amx_errors),
            "compiled_fd_tangent_rel_error": rel,
        }


if __name__ == "__main__":
    metrics = check()
    print("  max RHS parity abs err    = {:.3e}".format(
        metrics["max_rhs_abs_error"]))
    print("  max AMATRX parity rel err = {:.3e}".format(
        metrics["max_amatrx_rel_error"]))
    print("  compiled FD tangent rel   = {:.3e}".format(
        metrics["compiled_fd_tangent_rel_error"]))
