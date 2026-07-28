"""Regenerate, compile, and directly execute the generated UMAT with f2py."""

from pathlib import Path
import importlib
import math
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

from build import TemplateElastic
from check_reference import (
    FINAL_STRETCH_X,
    G,
    LAM,
    NINC,
    analytical_response,
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


def call_umat(module, stress, statev, dstran, stran):
    cmname = b"TEMPLATE" + b" " * 72
    return module.drive_umat(
        stress,
        statev,
        np.eye(3),
        np.eye(3),
        dstran,
        stran,
        np.array([G, LAM]),
        np.array([0.0, 0.0]),
        1.0 / NINC,
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

    committed = HERE / "template_umat.for"
    if not committed.exists():
        raise RuntimeError("Missing committed source; run `python build.py`")

    with tempfile.TemporaryDirectory(prefix="abaqus_ufl_template_") as tmp:
        work = Path(tmp)
        generated = work / committed.name
        au.generate_small_strain_umat(TemplateElastic(), str(generated))

        if generated.read_bytes() != committed.read_bytes():
            raise AssertionError(
                "Committed Fortran is stale; regenerate with `python build.py`"
            )
        print("[PASS] deterministic generated source")

        obj = work / "template_umat.o"
        run_checked(
            [
                "gfortran",
                "-c",
                "-ffixed-form",
                "-ffixed-line-length-none",
                generated,
                "-o",
                obj,
            ],
            work,
        )
        print("[PASS] gfortran compile")

        module_name = "_abaqus_ufl_template_umat"
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

        stress = np.zeros(6)
        statev = np.zeros(1)
        total_strain = math.log(FINAL_STRETCH_X)
        dstran = np.array([total_strain / NINC, 0.0, 0.0, 0.0, 0.0, 0.0])
        ddsdde = None
        for increment in range(NINC):
            stran = increment * dstran
            stress, statev, ddsdde, pnewdt = call_umat(
                module, stress, statev, dstran, stran
            )

        expected = analytical_response()
        if not np.allclose(
            stress, expected["stress_abq"], rtol=1.0e-11, atol=1.0e-11
        ):
            raise AssertionError(
                "Compiled stress mismatch:\nobserved={}\nexpected={}".format(
                    stress, expected["stress_abq"]
                )
            )

        expected_tangent = np.array(
            [
                [2.0 * G + LAM, LAM, LAM, 0.0, 0.0, 0.0],
                [LAM, 2.0 * G + LAM, LAM, 0.0, 0.0, 0.0],
                [LAM, LAM, 2.0 * G + LAM, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, G, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, G, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, G],
            ]
        )
        if not np.allclose(
            ddsdde, expected_tangent, rtol=1.0e-11, atol=1.0e-11
        ):
            raise AssertionError(
                "Compiled DDSDDE mismatch:\nobserved={}\nexpected={}".format(
                    ddsdde, expected_tangent
                )
            )
        if not np.isfinite(stress).all() or not np.isfinite(ddsdde).all():
            raise AssertionError("Compiled UMAT returned non-finite output")
        if pnewdt <= 0.0:
            raise AssertionError("Compiled UMAT returned invalid PNEWDT")

        print("[PASS] f2py UMAT call against closed form")
        return {
            "stress_max_abs_error": float(
                np.max(np.abs(stress - np.asarray(expected["stress_abq"])))
            ),
            "tangent_max_abs_error": float(
                np.max(np.abs(ddsdde - expected_tangent))
            ),
        }


if __name__ == "__main__":
    metrics = check()
    print("  stress max abs error  = {:.3e}".format(
        metrics["stress_max_abs_error"]
    ))
    print("  tangent max abs error = {:.3e}".format(
        metrics["tangent_max_abs_error"]
    ))
