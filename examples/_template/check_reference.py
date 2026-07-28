"""Check the Python model against a closed-form constrained-strain solution.

This is a physics oracle, not a code-vs-itself check.  It independently
evaluates isotropic elasticity for the same 2% compression path used by the
optional Abaqus deck.
"""

from pathlib import Path
import math
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from build import TemplateElastic
from tools.reference_utils import small_strain_path


G = 10.0
LAM = 20.0
FINAL_STRETCH_X = 0.98
NINC = 10


def analytical_response():
    """Return the closed-form response in Abaqus output conventions."""
    eps_abq = math.log(FINAL_STRETCH_X)
    eps_cp = -eps_abq
    return {
        "stress_abq": [
            -(2.0 * G + LAM) * eps_cp,
            -LAM * eps_cp,
            -LAM * eps_cp,
            0.0,
            0.0,
            0.0,
        ],
        "mises": 2.0 * G * eps_cp,
        "strain_abq": [eps_abq, 0.0, 0.0, 0.0, 0.0, 0.0],
    }


def check():
    observed = small_strain_path(
        TemplateElastic,
        final_stretch_x=FINAL_STRETCH_X,
        ninc=NINC,
    )
    expected = analytical_response()

    for name in ("stress_abq", "strain_abq"):
        if not np.allclose(
            observed[name], expected[name], rtol=1.0e-12, atol=1.0e-12
        ):
            raise AssertionError(
                "{} mismatch:\nobserved={}\nexpected={}".format(
                    name, observed[name], expected[name]
                )
            )
    if not math.isclose(
        observed["mises"], expected["mises"], rel_tol=1.0e-12, abs_tol=1.0e-12
    ):
        raise AssertionError(
            "mises mismatch: observed={} expected={}".format(
                observed["mises"], expected["mises"]
            )
        )
    return observed


if __name__ == "__main__":
    result = check()
    print("[PASS] closed-form constrained-strain reference")
    print("  stress = {}".format(result["stress_abq"]))
    print("  mises  = {:.12e}".format(result["mises"]))
