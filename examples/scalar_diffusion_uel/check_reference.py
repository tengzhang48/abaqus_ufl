"""Independent physics oracle for the thermo-mechanical Quad4 UEL.

Material-point closed forms evaluated without the framework:

- Fourier flux: for grad_T = g, flux = -k g exactly.
- Storage: for a temperature step dT over dt, storage = rho_cp dT/dt;
  for T = T_old, storage = 0 exactly.
- Thermal stress: at F = I, P = -K alpha T I exactly (the neo-Hookean
  part vanishes and F^{-T} = I).

Plus the WeakForm tangent-consistency gate at a state with nonzero
temperature, gradient, and deformation, and a broken control that
rejects a sign-flipped Fourier law.
"""

from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

G = 1.0
K = 10.0
ALPHA = 1e-3
COND = 0.5
RHO_CP = 1.0


def check():
    from build import (
        HeatDiffusionMaterial,
        HeatDiffusionProblem,
        verification_state,
    )

    problem = HeatDiffusionProblem()
    if not problem.verify(state=verification_state(), verbose=False):
        raise AssertionError("WeakForm verification failed")
    print("[PASS] WeakForm CS-vs-FD tangent consistency")

    mat = HeatDiffusionMaterial()

    grad_T = np.array([2.0, -1.0, 0.5])
    flux = np.asarray(mat.solvent_flux(np.eye(3), 10.0, grad_T), dtype=float)
    if not np.allclose(flux, -COND * grad_T, rtol=0, atol=1e-14):
        raise AssertionError("Fourier flux mismatch: {}".format(flux))
    print("[PASS] Fourier flux closed form")

    s0 = mat.solvent_storage(np.eye(3), np.eye(3), 7.0, 7.0, 0.1)
    s1 = mat.solvent_storage(np.eye(3), np.eye(3), 8.0, 5.0, 0.25)
    if abs(float(np.real(s0))) > 1e-14:
        raise AssertionError("storage nonzero for T = T_old")
    if abs(float(np.real(s1)) - RHO_CP * 3.0 / 0.25) > 1e-12:
        raise AssertionError("storage step mismatch: {}".format(s1))
    print("[PASS] storage closed forms (zero and step)")

    T = 25.0
    P = np.asarray(mat.stress_PK1(np.eye(3), T), dtype=complex).real
    expected = -K * ALPHA * T * np.eye(3)
    if not np.allclose(P, expected, rtol=0, atol=1e-12):
        raise AssertionError(
            "thermal stress at F=I mismatch:\n{}\nexpected\n{}".format(
                P, expected))
    print("[PASS] thermal stress closed form at F = I")

    # Broken control: a sign-flipped Fourier law must be rejected.
    wrong = +COND * grad_T
    if np.allclose(wrong, -COND * grad_T, rtol=1e-6, atol=1e-9):
        raise AssertionError("broken control failed: flux sign undetectable")
    print("[PASS] broken control (sign-flipped Fourier law is rejected)")


if __name__ == "__main__":
    check()
