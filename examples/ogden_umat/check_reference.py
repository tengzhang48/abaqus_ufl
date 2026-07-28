"""Independent physics oracle for the one-term Ogden example.

Three hand-derived checks that never call the implemented spectral
reconstruction, plus a broken control:

1. Pure dilation, F = c I: the isochoric stretches are 1, so
   ``sigma = K ln(c^3) / c^3 I`` exactly, for any alpha.

2. Isochoric uniaxial stretch, F = diag(l, 1/sqrt(l), 1/sqrt(l)), J = 1:
   the principal Cauchy stresses are the classical Ogden values
   ``sigma_i = (2 mu / alpha) (l_i^alpha - (1/3) sum_j l_j^alpha)``,
   evaluated directly from the diagonal stretches. The transverse
   eigenvalue is REPEATED, so this state also exercises the repeated
   spectrum of ``eig`` on the value path.

3. alpha = 2 degeneracy: the one-term Ogden energy with alpha = 2 equals
   the isochoric neo-Hookean energy, so at ANY deformation
   ``sigma = [ mu dev(bbar) + K ln(J) I ] / J`` with
   ``bbar = J^{-2/3} F F^T``. This closed tensor formula is evaluated in
   plain numpy (no eigendecomposition) at a rotated, non-isochoric F, so
   it checks the spectral reconstruction against an eig-free route on a
   fully generic state.

Broken control: dropping the isochoric split (using lambda_i instead of
``J^{-1/3} lambda_i``) must be rejected by check 3 at J != 1.
"""

from pathlib import Path
import math
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MU = 1.0
ALPHA = 3.5
K = 100.0

DILATION = 1.2
UNIAXIAL = 1.6


def rotation(angles=(0.4, -0.3, 0.7)):
    ax, ay, az = angles
    Rx = np.array([
        [1.0, 0.0, 0.0],
        [0.0, math.cos(ax), -math.sin(ax)],
        [0.0, math.sin(ax), math.cos(ax)],
    ])
    Ry = np.array([
        [math.cos(ay), 0.0, math.sin(ay)],
        [0.0, 1.0, 0.0],
        [-math.sin(ay), 0.0, math.cos(ay)],
    ])
    Rz = np.array([
        [math.cos(az), -math.sin(az), 0.0],
        [math.sin(az), math.cos(az), 0.0],
        [0.0, 0.0, 1.0],
    ])
    return Rz @ Ry @ Rx


def generic_F():
    """Rotated, sheared, non-isochoric deformation gradient."""
    U = np.array([
        [1.15, 0.08, 0.00],
        [0.08, 0.95, 0.05],
        [0.00, 0.05, 1.06],
    ])
    return rotation() @ U


def cauchy_from_model(model, F):
    P = np.asarray(model.stress_PK1(F), dtype=complex).real
    J = np.linalg.det(F)
    return P @ F.T / J


def neo_hookean_closed_form(F, mu, k_mod):
    """alpha = 2 Ogden equals isochoric neo-Hookean; eig-free formula."""
    J = np.linalg.det(F)
    bbar = J ** (-2.0 / 3.0) * (F @ F.T)
    dev_bbar = bbar - (np.trace(bbar) / 3.0) * np.eye(3)
    return (mu * dev_bbar + k_mod * math.log(J) * np.eye(3)) / J


def check():
    from build import OgdenOneTerm

    model = OgdenOneTerm()

    # 1. Pure dilation (triple-repeated eigenvalue on the value path).
    c = DILATION
    sigma = cauchy_from_model(model, c * np.eye(3))
    expected = K * math.log(c ** 3) / c ** 3 * np.eye(3)
    if not np.allclose(sigma, expected, rtol=1e-11, atol=1e-11):
        raise AssertionError(
            "dilation mismatch:\nobserved={}\nexpected={}".format(
                sigma, expected))
    print("[PASS] pure dilation closed form (triple-repeated spectrum)")

    # 2. Isochoric uniaxial (pair-repeated eigenvalue on the value path).
    l = UNIAXIAL
    stretches = np.array([l, 1.0 / math.sqrt(l), 1.0 / math.sqrt(l)])
    F = np.diag(stretches)
    powers = stretches ** ALPHA
    expected_diag = (2.0 * MU / ALPHA) * (powers - powers.sum() / 3.0)
    sigma = cauchy_from_model(model, F)
    if not np.allclose(np.diag(sigma), expected_diag, rtol=1e-11, atol=1e-11):
        raise AssertionError(
            "uniaxial mismatch:\nobserved={}\nexpected={}".format(
                np.diag(sigma), expected_diag))
    off = sigma - np.diag(np.diag(sigma))
    if np.max(np.abs(off)) > 1e-11:
        raise AssertionError("uniaxial state produced shear stress")
    print("[PASS] isochoric uniaxial closed form (pair-repeated spectrum)")

    # 3. alpha = 2 equals isochoric neo-Hookean at a generic rotated F.
    model2 = OgdenOneTerm(mu=MU, alpha=2.0, K=K)
    F = generic_F()
    sigma = cauchy_from_model(model2, F)
    expected = neo_hookean_closed_form(F, MU, K)
    if not np.allclose(sigma, expected, rtol=1e-10, atol=1e-12):
        raise AssertionError(
            "alpha=2 degeneracy mismatch:\nobserved={}\nexpected={}".format(
                sigma, expected))
    print("[PASS] alpha=2 degeneracy vs eig-free neo-Hookean at generic F")

    # Broken control: dropping the isochoric split must fail check 3.
    J = np.linalg.det(F)
    b = F @ F.T
    dev_b = b - (np.trace(b) / 3.0) * np.eye(3)
    sigma_wrong = (MU * dev_b + K * math.log(J) * np.eye(3)) / J
    if np.allclose(sigma_wrong, expected, rtol=1e-6, atol=1e-8):
        raise AssertionError(
            "broken control failed: missing J^(-1/3) split accepted")
    print("[PASS] broken control (missing isochoric split is rejected)")


if __name__ == "__main__":
    check()
