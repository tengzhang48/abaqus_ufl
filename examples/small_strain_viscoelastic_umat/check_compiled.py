"""Regenerate, compile, and directly execute the generated SLS UMAT.

The compiled drive repeats the step-shear relaxation of
``check_reference.py`` through the f2py one-point driver and checks three
things the Python oracle cannot: the ``DSTRAN`` engineering-shear
boundary, the relaxation history produced by the compiled state update,
and the COLUMN-MAJOR ``STATEV(1..9)`` layout of the tensor state
``eps_v`` (symmetric shear slots populated, diagonal empty, every slot
matching the closed form).
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

import abaqus_ufl as au

from build import SmallStrainViscoelastic
from check_reference import (
    DT,
    G_INF,
    G_V,
    K,
    NINC,
    SHEAR,
    TAU,
    discrete_relaxation,
    viscous_shear,
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


def call_umat(module, stress, statev, dstran, stran, time_now):
    cmname = b"SLS" + b" " * 77
    return module.drive_umat(
        stress,
        statev,
        np.eye(3),
        np.eye(3),
        dstran,
        stran,
        np.array([K, G_INF, G_V, TAU]),
        np.array([time_now, time_now]),
        DT,
        0.0,
        0.0,
        np.array([0.0]),
        np.array([0.0]),
        1.0,
        np.zeros(3),
        np.eye(3),
        3,
        3,
        cmname,
    )


def check():
    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required for the compiled pipeline gate")

    committed = HERE / "small_strain_viscoelastic.for"
    if not committed.exists():
        raise RuntimeError("Missing committed source; run `python build.py`")

    with tempfile.TemporaryDirectory(prefix="abaqus_ufl_sls_") as tmp:
        work = Path(tmp)
        generated = work / committed.name
        au.generate_small_strain_umat(SmallStrainViscoelastic(), str(generated))

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
                work / "small_strain_viscoelastic.o",
            ],
            work,
        )
        print("[PASS] gfortran compile")

        module_name = "_abaqus_ufl_small_strain_viscoelastic"
        run_checked(
            [
                sys.executable,
                "-m",
                "numpy.f2py",
                "-c",
                HERE / "f2py" / "drive_umat.f90",
                generated,
                "-m",
                module_name,
                "only:",
                "drive_umat",
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

        # Step engineering shear gamma = 2*SHEAR in increment 1, then hold.
        step = np.array([0.0, 0.0, 0.0, 2.0 * SHEAR, 0.0, 0.0])
        hold = np.zeros(6)

        stress = np.zeros(6)
        statev = np.zeros(9)
        stran = np.zeros(6)
        tau_errors, statev_errors = [], []
        for n in range(1, NINC + 1):
            dstran = step if n == 1 else hold
            stress, statev, ddsdde, pnewdt = call_umat(
                module, stress, statev, dstran, stran, (n - 1) * DT
            )
            stran = stran + dstran

            tau_ref = discrete_relaxation(n)
            tau_errors.append(abs(float(stress[3]) - tau_ref))
            if tau_errors[-1] > 1e-12:
                raise AssertionError(
                    "increment {}: compiled shear stress {} != closed form "
                    "{}".format(n, stress[3], tau_ref))

            # Column-major STATEV layout of the 3x3 tensor eps_v:
            # slot k = 3*(j-1) + i for component (i, j), 1-based.
            # STATEV stores the USER-side compression-positive tensor, so a
            # positive Abaqus DSTRAN shear yields a NEGATIVE stored eps_v
            # (the generated UMAT flips sign only at the STRESS/DSTRAN
            # boundary, never inside the state).
            expected = np.zeros(9)
            ev = -viscous_shear(n)
            expected[3] = ev   # (1,2) -> slot 4
            expected[1] = ev   # (2,1) -> slot 2
            statev_errors.append(float(np.max(np.abs(statev - expected))))
            if statev_errors[-1] > 1e-12:
                raise AssertionError(
                    "increment {}: STATEV layout mismatch:\nobserved={}\n"
                    "expected={}".format(n, statev, expected))

            if not (np.isfinite(stress).all() and np.isfinite(statev).all()
                    and np.isfinite(ddsdde).all()):
                raise AssertionError("non-finite UMAT output")
            if pnewdt <= 0.0:
                raise AssertionError("invalid PNEWDT")
            # Exact backward-Euler algorithmic tangent of the SLS in shear:
            # DDSDDE(4,4) = G_inf + G_v/(1+dt/tau), independent of history.
            expected_tan = G_INF + G_V / (1.0 + DT / TAU)
            if abs(ddsdde[3, 3] - expected_tan) > 1e-8 * expected_tan:
                raise AssertionError(
                    "DDSDDE(4,4)={} != closed-form SLS tangent {}".format(
                        ddsdde[3, 3], expected_tan))

        print("[PASS] f2py relaxation history vs closed form at every increment")
        print("[PASS] tensor STATEV column-major layout (symmetric shear "
              "slots, empty diagonal)")
        print("[PASS] exact SLS algorithmic tangent at every increment")
        return {
            "max_tau_abs_error": max(tau_errors),
            "max_statev_abs_error": max(statev_errors),
        }


if __name__ == "__main__":
    metrics = check()
    print("  max tau abs error    = {:.3e}".format(
        metrics["max_tau_abs_error"]
    ))
    print("  max statev abs error = {:.3e}".format(
        metrics["max_statev_abs_error"]
    ))
