"""Regenerate, compile, and directly execute the generated Ogden UMAT.

The compiled gates target the spectral path specifically:

- stress parity against the Python model at a generic rotated state, at a
  pair-repeated uniaxial state, and at a triple-repeated dilation;
- an eig-free closed form driven THROUGH the compiled module with
  ``alpha = 2`` runtime properties (Ogden degenerates to isochoric
  neo-Hookean);
- ``DDSDDE`` versus a Jaumann-corrected finite-difference tangent at both
  a distinct-spectrum state and the REPEATED-spectrum uniaxial state.
  The repeated-state tangent runs the complex-step perturbation through
  the rotation-safe ``eig33z`` fallback; the pre-correction guard
  (``V = I``) silently corrupted exactly this case.
"""

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

from build import OgdenOneTerm
from check_reference import (
    ALPHA,
    DILATION,
    K,
    MU,
    UNIAXIAL,
    cauchy_from_model,
    generic_F,
    neo_hookean_closed_form,
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


def voigt(sigma):
    return np.array([
        sigma[0, 0], sigma[1, 1], sigma[2, 2],
        sigma[0, 1], sigma[0, 2], sigma[1, 2],
    ])


def call_umat(module, F1, props):
    stress = np.zeros(6)
    statev = np.zeros(1)
    dstran = np.zeros(6)
    stran = np.zeros(6)
    cmname = b"OGDEN" + b" " * 75
    stress, statev, ddsdde, pnewdt = module.drive_umat(
        stress,
        statev,
        np.asarray(F1, dtype=float),
        np.eye(3),
        dstran,
        stran,
        np.asarray(props, dtype=float),
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


def fd_jaumann_tangent(module, F1, props, eps=1.0e-6):
    stress0, _, _ = call_umat(module, F1, props)
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
        stress_p, _, _ = call_umat(module, (np.eye(3) + de) @ F1, props)
        tangent[:, j] = (stress_p - stress0) / eps + (tr_de / eps) * stress0
    return stress0, tangent


def check():
    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required for the compiled pipeline gate")

    committed = HERE / "ogden_umat.for"
    if not committed.exists():
        raise RuntimeError("Missing committed source; run `python build.py`")

    with tempfile.TemporaryDirectory(prefix="abaqus_ufl_ogden_") as tmp:
        work = Path(tmp)
        generated = work / committed.name
        au.generate_umat(OgdenOneTerm(), str(generated))

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
                work / "ogden_umat.o",
            ],
            work,
        )
        print("[PASS] gfortran compile")

        module_name = "_abaqus_ufl_ogden_umat"
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

        props = [MU, ALPHA, K]
        model = OgdenOneTerm()
        uniaxial_F = np.diag([
            UNIAXIAL,
            1.0 / math.sqrt(UNIAXIAL),
            1.0 / math.sqrt(UNIAXIAL),
        ])
        states = (
            ("generic rotated", generic_F()),
            ("pair-repeated uniaxial", uniaxial_F),
            ("triple-repeated dilation", DILATION * np.eye(3)),
        )
        parity_errors = []
        for name, F in states:
            stress, ddsdde, pnewdt = call_umat(module, F, props)
            expected = voigt(cauchy_from_model(model, F))
            if not np.allclose(stress, expected, rtol=1e-10, atol=1e-12):
                raise AssertionError(
                    "compiled {} stress parity failed:\nobserved={}\n"
                    "expected={}".format(name, stress, expected))
            if not np.isfinite(ddsdde).all():
                raise AssertionError("non-finite DDSDDE at " + name)
            if pnewdt <= 0.0:
                raise AssertionError("invalid PNEWDT at " + name)
            parity_errors.append(float(np.max(np.abs(stress - expected))))
            print("[PASS] f2py stress parity at {} state".format(name))

        # Runtime-property degeneracy: alpha = 2 through the compiled module
        # against the eig-free closed form.
        F = generic_F()
        stress, _, _ = call_umat(module, F, [MU, 2.0, K])
        expected = voigt(neo_hookean_closed_form(F, MU, K))
        if not np.allclose(stress, expected, rtol=1e-9, atol=1e-11):
            raise AssertionError(
                "compiled alpha=2 closed form failed:\nobserved={}\n"
                "expected={}".format(stress, expected))
        closed_form_error = float(np.max(np.abs(stress - expected)))
        print("[PASS] compiled alpha=2 vs eig-free closed form at generic F")

        # Tangent gates: distinct AND repeated spectra. The repeated case is
        # the regression the corrected eig33z fallback exists for.
        tangent_errors = {}
        for name, F in (
            ("distinct", generic_F()),
            ("repeated", uniaxial_F),
        ):
            ddsdde = call_umat(module, F, props)[1]
            _, fd = fd_jaumann_tangent(module, F, props)
            denom = max(1.0, float(np.max(np.abs(ddsdde))))
            rel = float(np.max(np.abs(ddsdde - fd))) / denom
            if not np.isfinite(rel):
                raise AssertionError(
                    "non-finite tangent comparison at " + name)
            if rel > 5.0e-5:
                raise AssertionError(
                    "DDSDDE vs FD-Jaumann failed at {} spectrum "
                    "(rel err {:.3e})".format(name, rel))
            tangent_errors[name] = rel
            print("[PASS] DDSDDE vs FD-Jaumann tangent at {} spectrum".format(
                name))

        return {
            "stress_parity_max_abs_error": max(parity_errors),
            "alpha2_closed_form_max_abs_error": closed_form_error,
            "tangent_rel_error_distinct": tangent_errors["distinct"],
            "tangent_rel_error_repeated": tangent_errors["repeated"],
        }


if __name__ == "__main__":
    metrics = check()
    print("  stress parity max abs error = {:.3e}".format(
        metrics["stress_parity_max_abs_error"]))
    print("  alpha=2 closed form err     = {:.3e}".format(
        metrics["alpha2_closed_form_max_abs_error"]))
    print("  tangent rel err (distinct)  = {:.3e}".format(
        metrics["tangent_rel_error_distinct"]))
    print("  tangent rel err (repeated)  = {:.3e}".format(
        metrics["tangent_rel_error_repeated"]))
