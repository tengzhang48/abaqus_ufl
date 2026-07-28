"""Regenerate, compile, and directly execute the generated J2 UMAT.

The compiled drive repeats the monotonic pure-shear path of
``check_reference.py`` through the f2py one-point driver. Because Abaqus
passes ENGINEERING shear in ``DSTRAN(4)`` (gamma = 2 eps_xy), this gate
also verifies the generated Voigt/tensor conversion at the subroutine
boundary. ``STATEV(1)`` is compared with the closed-form equivalent
plastic strain after every increment, which checks the state round trip
through the ``TIME(2)==0`` initialization branch and the reread path.
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

from build import SmallStrainJ2
from check_reference import (
    FINAL_SHEAR,
    G,
    H,
    LAM,
    NINC,
    SIGMA_Y,
    YIELD_SHEAR,
    closed_form,
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


def call_umat(module, stress, statev, dstran, stran, time_now, dtime):
    cmname = b"J2" + b" " * 78
    return module.drive_umat(
        stress,
        statev,
        np.eye(3),
        np.eye(3),
        dstran,
        stran,
        np.array([G, LAM, SIGMA_Y, H]),
        np.array([time_now, time_now]),
        dtime,
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

    committed = HERE / "small_strain_j2.for"
    if not committed.exists():
        raise RuntimeError("Missing committed source; run `python build.py`")

    with tempfile.TemporaryDirectory(prefix="abaqus_ufl_j2_") as tmp:
        work = Path(tmp)
        generated = work / committed.name
        au.generate_small_strain_umat(SmallStrainJ2(), str(generated))

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
                work / "small_strain_j2.o",
            ],
            work,
        )
        print("[PASS] gfortran compile")

        module_name = "_abaqus_ufl_small_strain_j2"
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

        # Monotonic pure shear in Abaqus convention: DSTRAN(4) carries the
        # ENGINEERING shear increment gamma = 2 * (tensor shear increment).
        de = FINAL_SHEAR / NINC
        dstran = np.array([0.0, 0.0, 0.0, 2.0 * de, 0.0, 0.0])
        dtime = 1.0 / NINC

        stress = np.zeros(6)
        statev = np.zeros(1)
        elastic_checked = False
        ep_errors = []
        for inc in range(NINC):
            stran = inc * dstran
            stress, statev, ddsdde, pnewdt = call_umat(
                module, stress, statev, dstran, stran, inc * dtime, dtime
            )
            if not (np.isfinite(stress).all() and np.isfinite(statev).all()
                    and np.isfinite(ddsdde).all()):
                raise AssertionError("non-finite UMAT output at increment "
                                     + str(inc + 1))
            e_now = (inc + 1) * de
            tau_ref, ep_ref = closed_form(e_now)
            if abs(stress[3] - tau_ref) > 1e-10:
                raise AssertionError(
                    "increment {}: compiled shear stress {} != closed form "
                    "{}".format(inc + 1, stress[3], tau_ref))
            ep_errors.append(abs(float(statev[0]) - ep_ref))
            if ep_errors[-1] > 1e-10:
                raise AssertionError(
                    "increment {}: compiled STATEV(1)={} != closed-form "
                    "ep={}".format(inc + 1, statev[0], ep_ref))
            if not np.isfinite(ddsdde).all():
                raise AssertionError("non-finite DDSDDE")
            if pnewdt <= 0.0:
                raise AssertionError("invalid PNEWDT")
            if not elastic_checked and e_now < YIELD_SHEAR:
                if abs(statev[0]) > 1e-14:
                    raise AssertionError("plasticity before yield")
                # Elastic engineering-shear tangent must be exactly G.
                if abs(ddsdde[3, 3] - G) > 1e-8:
                    raise AssertionError(
                        "elastic DDSDDE(4,4)={} != G".format(ddsdde[3, 3]))
                elastic_checked = True

        if not elastic_checked:
            raise AssertionError("path never sampled the elastic branch")

        # Exact algorithmic consistent tangent for radial return with
        # linear isotropic hardening in monotonic pure shear:
        # engineering-shear modulus DDSDDE(4,4) = G H / (3 G + H).
        expected_plastic = G * H / (3.0 * G + H)
        if abs(ddsdde[3, 3] - expected_plastic) > 1e-6 * expected_plastic:
            raise AssertionError(
                "plastic DDSDDE(4,4)={} != closed-form consistent tangent "
                "{}".format(ddsdde[3, 3], expected_plastic))

        tau_ref, ep_ref = closed_form(FINAL_SHEAR)
        print("[PASS] f2py shear path vs closed form at every increment")
        print("[PASS] STATEV(1) round trip matches closed-form ep")
        print("[PASS] elastic DDSDDE(4,4)=G; plastic DDSDDE(4,4)="
              "GH/(3G+H) exact")
        return {
            "final_tau_abs_error": abs(float(stress[3]) - tau_ref),
            "max_ep_abs_error": max(ep_errors),
        }


if __name__ == "__main__":
    metrics = check()
    print("  final tau abs error = {:.3e}".format(
        metrics["final_tau_abs_error"]
    ))
    print("  max ep abs error    = {:.3e}".format(metrics["max_ep_abs_error"]))
