"""Independent physics oracle for the SLS viscoelastic example.

Step-shear relaxation admits an exact closed form for the DISCRETE
backward-Euler update, hand-derived from the recursion. With a step
tensor shear strain ``e`` applied in increment 1 and then held, and
``r = dt / tau``,

```text
dev_e - eps_v^n = dev_e / (1 + r)^n
tau_n = 2 e [ G_inf + G_v / (1 + r)^n ]      (shear stress after n increments)
```

so the stress relaxes geometrically from the instantaneous limit
``2 e (G_inf + G_v)`` toward the equilibrium limit ``2 e G_inf``. As
``r -> 0`` the factor ``(1+r)^{-n}`` converges to the continuous
``exp(-t/tau)``. These formulas never call the implemented update.
"""

from pathlib import Path
import math
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

K = 30.0
G_INF = 10.0
G_V = 5.0
TAU = 0.5

SHEAR = 0.01          # tensor shear strain of the held step
DT = 0.1              # increment size; r = DT/TAU = 0.2
NINC = 30             # step increment + 29 holds


def discrete_relaxation(n):
    """Exact shear stress of the backward-Euler SLS after n increments."""
    r = DT / TAU
    return 2.0 * SHEAR * (G_INF + G_V / (1.0 + r) ** n)


def viscous_shear(n):
    """Exact viscous shear-strain component after n increments."""
    r = DT / TAU
    return SHEAR * (1.0 - 1.0 / (1.0 + r) ** n)


def drive_python(model, ninc):
    """Step shear then hold; return per-increment (tau, eps_v) lists."""
    dstrain = np.zeros((3, 3))
    dstrain[0, 1] = dstrain[1, 0] = SHEAR
    hold = np.zeros((3, 3))
    sigma = np.zeros((3, 3))
    strain = np.zeros((3, 3))
    eps_v = np.zeros((3, 3))
    taus, eps_vs = [], []
    for n in range(ninc):
        de = dstrain if n == 0 else hold
        sigma, state = model.stress_update(sigma, strain, de, eps_v, DT)
        sigma = np.real(np.asarray(sigma, dtype=complex))
        eps_v = np.real(np.asarray(state["eps_v"], dtype=complex))
        strain = strain + de
        taus.append(float(sigma[0, 1]))
        eps_vs.append(eps_v.copy())
    return taus, eps_vs


def check():
    from build import SmallStrainViscoelastic

    model = SmallStrainViscoelastic()
    taus, eps_vs = drive_python(model, NINC)

    for n in range(1, NINC + 1):
        expected = discrete_relaxation(n)
        if abs(taus[n - 1] - expected) > 1e-13:
            raise AssertionError(
                "increment {}: tau {} != closed form {}".format(
                    n, taus[n - 1], expected))
        if abs(eps_vs[n - 1][0, 1] - viscous_shear(n)) > 1e-13:
            raise AssertionError(
                "increment {}: eps_v {} != closed form {}".format(
                    n, eps_vs[n - 1][0, 1], viscous_shear(n)))
    print("[PASS] discrete relaxation closed form at every increment")

    # Physical limits bracket the response.
    instantaneous = 2.0 * SHEAR * (G_INF + G_V)
    equilibrium = 2.0 * SHEAR * G_INF
    assert taus[0] < instantaneous, "first increment must already relax"
    assert taus[-1] > equilibrium, "stress must stay above equilibrium"
    decayed = (instantaneous - taus[-1]) / (instantaneous - equilibrium)
    assert decayed > 0.95, "relaxation path barely decayed; extend NINC"
    print("[PASS] instantaneous/equilibrium limits bracket the response")

    # Continuous-limit consistency: halving dt moves the discrete factor
    # toward exp(-t/tau) (first-order backward Euler).
    t_total = 1.0
    errors = []
    for dt in (0.1, 0.05):
        n = int(round(t_total / dt))
        discrete = 1.0 / (1.0 + dt / TAU) ** n
        errors.append(abs(discrete - math.exp(-t_total / TAU)))
    assert errors[1] < 0.6 * errors[0], "discrete factor not converging"
    print("[PASS] backward-Euler factor converges to exp(-t/tau)")

    # Broken control: dropping the implicit denominator (forward-Euler-like
    # dashpot) must be rejected by the discrete closed form.
    r = DT / TAU
    eps_v_wrong = 0.0
    for _ in range(NINC):
        eps_v_wrong = eps_v_wrong + r * (SHEAR - eps_v_wrong)
    tau_wrong = 2.0 * SHEAR * G_INF + 2.0 * G_V * (SHEAR - eps_v_wrong)
    if abs(tau_wrong - discrete_relaxation(NINC)) < 1e-8:
        raise AssertionError(
            "broken control failed: explicit dashpot accepted")
    print("[PASS] broken control (explicit-dashpot update is rejected)")


if __name__ == "__main__":
    check()
