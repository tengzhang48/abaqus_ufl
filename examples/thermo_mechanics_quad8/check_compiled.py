"""Regenerate, compile, and directly execute the mixed-order Quad8 UEL.

Gates: deterministic regeneration, gfortran compile, f2py RHS/AMATRX
parity against the Python reference assembly (itself gated by
check_assembled.py) at the reference state, the exact-stretch pull-back
state, and a generic coupled state on a distorted element, plus an
independent finite-difference check of the compiled 20x20 tangent.
The parity states exercise the node-dependent mixed-order DOF maps in
the exact Fortran Abaqus would compile.
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
    NDOFEL,
    build_U,
    distorted_quad,
    props_array,
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

    committed = HERE / "thermo_mechanics_quad8_uel.for"
    if not committed.exists():
        raise RuntimeError("Missing committed source; run `python build.py`")

    from build import ThermoMechanicalQuad8

    with tempfile.TemporaryDirectory(prefix="abaqus_ufl_tm8_") as tmp:
        work = Path(tmp)
        generated = work / committed.name
        generate_uel(ThermoMechanicalQuad8(), str(generated),
                     element="Quad8", formulation="standard")

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
                work / "thermo_mechanics_quad8_uel.o",
            ],
            work,
        )
        print("[PASS] gfortran compile")

        module_name = "_abaqus_ufl_thermo_quad8_uel"
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

        problem = ThermoMechanicalQuad8()
        U_pull = build_U(unit_square(), 1.3, lambda x, y: 2.0 * x)
        U_gen = build_U(distorted_quad(), 1.05,
                        lambda x, y: 1.5 * x + 0.4 * y)
        U_gen = U_gen + np.linspace(-0.02, 0.02, NDOFEL)
        states = (
            ("reference", unit_square(), np.zeros(NDOFEL), np.zeros(NDOFEL)),
            ("stretch + pull-back flux", unit_square(), U_pull,
             np.zeros(NDOFEL)),
            ("generic coupled", distorted_quad(), U_gen, 0.3 * U_gen),
        )
        rhs_errors, amx_errors = [], []
        for name, coords, U, DU in states:
            rhs_f, amx_f, pnewdt = call_uel(module, coords, U, DU)
            if not (np.isfinite(rhs_f).all() and np.isfinite(amx_f).all()):
                raise AssertionError("non-finite compiled output at " + name)
            rhs_p, amx_p = assemble_element(
                problem, coords, U, DU, 0.1, props_array(), element="quad8")
            scale_r = max(1.0e-12, float(np.max(np.abs(rhs_p))))
            scale_a = float(np.max(np.abs(amx_p)))
            err_r = float(np.max(np.abs(rhs_f - rhs_p)))
            err_a = float(np.max(np.abs(amx_f - amx_p))) / scale_a
            if err_r > 1e-9 * max(1.0, scale_r) or err_a > 1e-9:
                raise AssertionError(
                    "parity failed at {}: rhs {} amatrx {}".format(
                        name, err_r, err_a))
            if pnewdt <= 0.0:
                raise AssertionError("invalid PNEWDT at " + name)
            rhs_errors.append(err_r)
            amx_errors.append(err_a)
            print("[PASS] f2py RHS/AMATRX parity at {} state".format(name))

        coords, U, DU = distorted_quad(), U_gen, 0.3 * U_gen
        rhs0, amx0, _ = call_uel(module, coords, U, DU)
        eps = 1e-7
        fd = np.zeros((NDOFEL, NDOFEL))
        for j in range(NDOFEL):
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
        for node in range(8):
            ux = 3 * node if node < 4 else 12 + 2 * (node - 4)
            U_bad[ux] = -2.0 * coords[0][node]
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
