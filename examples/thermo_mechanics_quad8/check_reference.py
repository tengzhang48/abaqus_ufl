"""Independent physics oracle for the mixed-order thermo Quad8 UEL.

Hand-derived material-point closed forms, evaluated without the
framework:

- pulled-back flux: for F = diag(l, 1, 1) and referential gradient
  g = (a, 0, 0), C^{-1} = diag(1/l^2, 1, 1), so flux_x = -kappa a / l^2
  exactly (the deformation dependence IS the point of this example);
- storage: cT dT/dt, zero at T = T_old;
- thermal stress at F = I: P = -K alpha T I.

A broken control confirms that omitting the C^{-1} pull-back is rejected
at l != 1.
"""

from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

G = 1.0
K = 100.0
ALPHA = 1.0e-3
KAPPA = 0.25
CT = 1.0


def check():
    from build import (
        ThermoMechanicalMaterial,
        ThermoMechanicalQuad8,
        verification_state,
    )

    problem = ThermoMechanicalQuad8()
    if not problem.verify(state=verification_state(), verbose=False):
        raise AssertionError("WeakForm verification failed")
    print("[PASS] WeakForm CS-vs-FD tangent consistency")

    mat = ThermoMechanicalMaterial()

    # Pulled-back flux under uniaxial stretch.
    lam = 1.3
    F = np.diag([lam, 1.0, 1.0])
    g = np.array([2.0, 0.0, 0.0])
    flux = np.asarray(mat.solvent_flux(F, 5.0, g), dtype=complex).real
    expected = np.array([-KAPPA * g[0] / lam ** 2, 0.0, 0.0])
    if not np.allclose(flux, expected, rtol=0, atol=1e-14):
        raise AssertionError(
            "pull-back flux mismatch: {} vs {}".format(flux, expected))
    flux1 = np.asarray(
        mat.solvent_flux(np.eye(3), 5.0, g), dtype=complex).real
    if abs(flux[0] / flux1[0] - 1.0 / lam ** 2) > 1e-12:
        raise AssertionError("pull-back ratio wrong")
    print("[PASS] C^-1 pulled-back flux closed form (factor 1/l^2)")

    s0 = mat.solvent_storage(np.eye(3), np.eye(3), 4.0, 4.0, 0.1)
    s1 = mat.solvent_storage(np.eye(3), np.eye(3), 6.0, 2.0, 0.5)
    if abs(float(np.real(s0))) > 1e-14 or \
            abs(float(np.real(s1)) - CT * 4.0 / 0.5) > 1e-12:
        raise AssertionError("storage closed forms failed")
    print("[PASS] storage closed forms (zero and step)")

    T = 12.0
    P = np.asarray(mat.stress_PK1(np.eye(3), T), dtype=complex).real
    if not np.allclose(P, -K * ALPHA * T * np.eye(3), rtol=0, atol=1e-12):
        raise AssertionError("thermal stress at F=I mismatch")
    print("[PASS] thermal stress closed form at F = I")

    # Broken control: no pull-back means no 1/l^2 factor.
    wrong = -KAPPA * g[0]
    if abs(wrong - expected[0]) < 1e-6:
        raise AssertionError("broken control failed: pull-back undetectable")
    print("[PASS] broken control (missing C^-1 pull-back is rejected)")


if __name__ == "__main__":
    check()
