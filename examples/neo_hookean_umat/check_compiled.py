"""Regenerate, compile, and directly execute the generated neo-Hookean UMAT.

Gates, in order: deterministic regeneration against the committed source,
gfortran compile, and a direct one-point f2py call whose ``STRESS`` is
compared with the hand-derived closed forms at a uniaxial-stretch state and
a simple-shear state. ``DDSDDE`` is checked against a Jaumann-corrected
finite-difference tangent built from repeated compiled stress calls, so the
tangent evidence is independent of the complex-step machinery that
generated it.
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

from build import NeoHookean
from check_reference import (
    G,
    K,
    shear_F,
    shear_cauchy_closed_form,
    uniaxial_F,
    uniaxial_cauchy_closed_form,
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


def call_umat(module, F1):
    stress = np.zeros(6)
    statev = np.zeros(1)
    dstran = np.zeros(6)
    stran = np.zeros(6)
    cmname = b"NEOHOOKEAN" + b" " * 70
    stress, statev, ddsdde, pnewdt = module.drive_umat(
        stress,
        statev,
        np.asarray(F1, dtype=float),
        np.eye(3),
        dstran,
        stran,
        np.array([G, K]),
        np.array([0.0, 0.0]),
        1.0,
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
    return stress, ddsdde, pnewdt


def fd_jaumann_tangent(module, F1, eps=1.0e-6):
    """Jaumann-corrected finite-difference DDSDDE from compiled calls.

    Column j perturbs F multiplicatively with a small symmetric strain,
    F' = (I + de) F, and corrects with the tr(de) sigma term of the
    Jaumann rate of Kirchhoff stress divided by J.
    """
    stress0, _, _ = call_umat(module, F1)
    tangent = np.zeros((6, 6))
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    for j, (a, b) in enumerate(pairs):
        de = np.zeros((3, 3))
        if a == b:
            de[a, a] = eps
            tr_de = eps
        else:
            de[a, b] = eps / 2.0
            de[b, a] = eps / 2.0
            tr_de = 0.0
        stress_p, _, _ = call_umat(module, (np.eye(3) + de) @ F1)
        tangent[:, j] = (stress_p - stress0) / eps + (tr_de / eps) * stress0
    return stress0, tangent


def check():
    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required for the compiled pipeline gate")

    committed = HERE / "neo_hookean_umat.for"
    if not committed.exists():
        raise RuntimeError("Missing committed source; run `python build.py`")

    with tempfile.TemporaryDirectory(prefix="abaqus_ufl_neo_") as tmp:
        work = Path(tmp)
        generated = work / committed.name
        au.generate_umat(NeoHookean(), str(generated))

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
                work / "neo_hookean_umat.o",
            ],
            work,
        )
        print("[PASS] gfortran compile")

        module_name = "_abaqus_ufl_neo_hookean_umat"
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

        stress_errors = []
        for name, F1, expected in (
            ("uniaxial", uniaxial_F(), uniaxial_cauchy_closed_form()),
            ("simple shear", shear_F(), shear_cauchy_closed_form()),
        ):
            stress, ddsdde, pnewdt = call_umat(module, F1)
            if not np.allclose(stress, expected, rtol=1e-9, atol=1e-11):
                raise AssertionError(
                    "Compiled {} stress mismatch:\n"
                    "observed={}\nexpected={}".format(name, stress, expected)
                )
            if not np.isfinite(ddsdde).all():
                raise AssertionError("Compiled DDSDDE contains non-finite values")
            if pnewdt <= 0.0:
                raise AssertionError("Compiled UMAT returned invalid PNEWDT")
            stress_errors.append(float(np.max(np.abs(stress - expected))))
            print("[PASS] f2py {} stress vs closed form".format(name))

        _, ddsdde, _ = call_umat(module, uniaxial_F())
        _, fd = fd_jaumann_tangent(module, uniaxial_F())
        denom = max(1.0, float(np.max(np.abs(ddsdde))))
        tangent_error = float(np.max(np.abs(ddsdde - fd))) / denom
        if not np.isfinite(tangent_error):
            raise AssertionError("non-finite tangent comparison")
        if tangent_error > 5.0e-5:
            raise AssertionError(
                "Compiled DDSDDE disagrees with Jaumann-corrected FD tangent "
                "(relative error {:.3e})".format(tangent_error)
            )
        print("[PASS] DDSDDE vs Jaumann-corrected FD tangent")
        return {
            "stress_max_abs_error": max(stress_errors),
            "tangent_rel_error": tangent_error,
        }


if __name__ == "__main__":
    metrics = check()
    print("  stress max abs error = {:.3e}".format(
        metrics["stress_max_abs_error"]
    ))
    print("  tangent rel error    = {:.3e}".format(
        metrics["tangent_rel_error"]
    ))
