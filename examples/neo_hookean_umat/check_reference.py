"""Independent physics oracle for the neo-Hookean example.

The expected Cauchy stresses below are hand-derived closed forms; they are
evaluated without calling the implemented ``stress_PK1``.

Uniaxial stretch, F = diag(lam, 1, 1), J = lam:

    sigma_xx = G (lam - 1/lam) + K ln(lam) / lam
    sigma_yy = sigma_zz = K ln(lam) / lam

Simple shear, F = I + gamma e_x otimes e_y, J = 1:

    sigma = G [[gamma^2, gamma, 0],
               [gamma,   0,     0],
               [0,       0,     0]]

which shows the classical neo-Hookean normal-stress effect in shear.
"""

from pathlib import Path
import math
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

G = 0.5
K = 50.0
STRETCH = 1.2
GAMMA = 0.3


def uniaxial_F():
    return np.diag([STRETCH, 1.0, 1.0])


def shear_F():
    F = np.eye(3)
    F[0, 1] = GAMMA
    return F


def uniaxial_cauchy_closed_form():
    lam = STRETCH
    sxx = G * (lam - 1.0 / lam) + K * math.log(lam) / lam
    syy = K * math.log(lam) / lam
    # Abaqus Voigt order (11, 22, 33, 12, 13, 23)
    return np.array([sxx, syy, syy, 0.0, 0.0, 0.0])


def shear_cauchy_closed_form():
    g = GAMMA
    return np.array([G * g * g, 0.0, 0.0, G * g, 0.0, 0.0])


def python_cauchy(model, F):
    """Cauchy stress from the implemented model, sigma = P F^T / J."""
    P = model.stress_PK1(F)
    J = np.linalg.det(F)
    sigma = P @ F.T / J
    return np.array([
        sigma[0, 0], sigma[1, 1], sigma[2, 2],
        sigma[0, 1], sigma[0, 2], sigma[1, 2],
    ])


def check():
    from build import NeoHookean

    model = NeoHookean()

    for name, F, expected in (
        ("uniaxial", uniaxial_F(), uniaxial_cauchy_closed_form()),
        ("simple shear", shear_F(), shear_cauchy_closed_form()),
    ):
        observed = python_cauchy(model, F)
        if not np.allclose(observed, expected, rtol=1e-12, atol=1e-12):
            raise AssertionError(
                "{} mismatch:\nobserved={}\nexpected={}".format(
                    name, observed, expected
                )
            )
        print("[PASS] closed-form {} state".format(name))
        print("  sigma = {}".format(np.array2string(expected, precision=6)))

    # Broken control: the oracle must reject a deliberately wrong response
    # (volumetric term dropped). If this "defect" passed, the check would
    # be too loose to gate anything.
    F = uniaxial_F()
    P_wrong = G * (F - np.linalg.inv(F).T)
    sigma_wrong = P_wrong @ F.T / np.linalg.det(F)
    wrong = np.array([
        sigma_wrong[0, 0], sigma_wrong[1, 1], sigma_wrong[2, 2],
        sigma_wrong[0, 1], sigma_wrong[0, 2], sigma_wrong[1, 2],
    ])
    if np.allclose(wrong, uniaxial_cauchy_closed_form(), rtol=1e-6, atol=1e-8):
        raise AssertionError(
            "Broken control failed: oracle accepted the K-term-free response"
        )
    print("[PASS] broken control (missing volumetric term is rejected)")


if __name__ == "__main__":
    check()
