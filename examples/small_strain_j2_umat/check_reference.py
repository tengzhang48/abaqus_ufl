"""Independent physics oracle for the small-strain J2 example.

Monotonic pure shear admits an exact closed form for J2 with linear
isotropic hardening, because radial return is exact under proportional
deviatoric loading. With tensor shear strain ``e`` (component
``eps_xy = e``), shear stress ``tau = sigma_xy``, and ``n_xy = sqrt(3)/2``,

```text
elastic:  tau = 2 G e                    while sqrt(3) tau <= sigma_y
yield:    e_y = sigma_y / (2 sqrt(3) G)
plastic:  ep  = (2 sqrt(3) G e - sigma_y) / (3 G + H)
          tau = (sigma_y + H ep) / sqrt(3)
```

These expressions are hand-derived from the consistency condition
``sqrt(3) tau = sigma_y + H ep`` and the split
``e = tau / (2G) + (sqrt(3)/2) ep``; they never call the implemented
update. The oracle drives the Python model incrementally through both
branches and compares stress and state against the closed form.
"""

from pathlib import Path
import math
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

G = 10.0
LAM = 20.0
SIGMA_Y = 0.1
H = 5.0

FINAL_SHEAR = 0.02      # tensor shear strain, well past yield
NINC = 40

YIELD_SHEAR = SIGMA_Y / (2.0 * math.sqrt(3.0) * G)


def closed_form(e):
    """Return (tau, ep) for monotonic tensor shear strain e >= 0."""
    if math.sqrt(3.0) * 2.0 * G * e <= SIGMA_Y:
        return 2.0 * G * e, 0.0
    ep = (2.0 * math.sqrt(3.0) * G * e - SIGMA_Y) / (3.0 * G + H)
    tau = (SIGMA_Y + H * ep) / math.sqrt(3.0)
    return tau, ep


def drive_python(model, e_total, ninc):
    """Incremental pure-shear drive of the Python stress_update."""
    de = e_total / ninc
    dstrain = np.zeros((3, 3))
    dstrain[0, 1] = dstrain[1, 0] = de
    sigma = np.zeros((3, 3))
    strain = np.zeros((3, 3))
    ep = 0.0
    for _ in range(ninc):
        sigma, state = model.stress_update(
            sigma, strain, dstrain, ep, 1.0 / ninc)
        sigma = np.real(np.asarray(sigma, dtype=complex))
        ep = float(np.real(state["ep"]))
        strain = strain + dstrain
    return sigma, ep


def check():
    from build import SmallStrainJ2

    model = SmallStrainJ2()

    # Elastic branch: stop safely below yield, no plastic strain.
    e_el = 0.5 * YIELD_SHEAR
    sigma, ep = drive_python(model, e_el, 10)
    tau_expected, ep_expected = closed_form(e_el)
    assert abs(ep - ep_expected) < 1e-14, "elastic branch produced plasticity"
    assert abs(sigma[0, 1] - tau_expected) < 1e-12
    print("[PASS] elastic branch (tau = 2 G e, ep = 0)")

    # Plastic branch: monotonic loading past yield matches the closed form.
    sigma, ep = drive_python(model, FINAL_SHEAR, NINC)
    tau_expected, ep_expected = closed_form(FINAL_SHEAR)
    if abs(sigma[0, 1] - tau_expected) > 1e-10:
        raise AssertionError(
            "plastic shear stress mismatch: observed {} expected {}".format(
                sigma[0, 1], tau_expected))
    if abs(ep - ep_expected) > 1e-10:
        raise AssertionError(
            "equivalent plastic strain mismatch: observed {} expected {}".format(
                ep, ep_expected))
    print("[PASS] plastic branch (consistency closed form)")
    print("  tau = {:.10f}, ep = {:.10f}".format(tau_expected, ep_expected))

    # Regime coverage: the plastic run must actually have yielded.
    assert ep_expected > 5.0 * YIELD_SHEAR, "test path barely yields"

    # Broken control: perfect-plasticity mistake (H dropped from the return
    # denominator) must be rejected by the closed form.
    ep_wrong = (2.0 * math.sqrt(3.0) * G * FINAL_SHEAR - SIGMA_Y) / (3.0 * G)
    if abs(ep_wrong - ep_expected) < 1e-6:
        raise AssertionError("broken control failed: H-free return accepted")
    print("[PASS] broken control (H-free return mapping is rejected)")


if __name__ == "__main__":
    check()
