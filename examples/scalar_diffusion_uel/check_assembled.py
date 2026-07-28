"""Assembled-element checks for the thermo-mechanical Quad4 UEL.

``problem.verify()`` gates material-point derivatives only, so this
script exercises the ASSEMBLED residual and tangent through the Python
reference assembly, which mirrors the generated Fortran element:

1. equilibrium: undeformed element, uniform T = T_old -> RHS = 0;
2. invariance: rigid translation adds no mechanical residual;
3. heat balance (partition of unity): summing the thermal residual rows
   kills the flux term exactly, so on a DISTORTED element with a linear
   steady temperature field the thermal-row sum is zero, and with a
   uniform temperature step on the unit square it equals
   rho_cp dT/dt * area, both hand-derivable;
4. assembled tangent: AMATRX vs -dRHS/dU by finite differences over all
   12 DOFs.

DOF layout (generated header): per node (u1, u2, T), four corner nodes,
NDOFEL = 12; T DOFs are 2, 5, 8, 11 (0-based).
"""

from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from abaqus_ufl.core.reference_assembly import assemble_element

from check_reference import RHO_CP

T_DOFS = [2, 5, 8, 11]
U_DOFS = [0, 1, 3, 4, 6, 7, 9, 10]


def props_array():
    from build import HeatDiffusionMaterial
    mat = HeatDiffusionMaterial()
    return np.array([getattr(mat, name) for name in mat.props])


def unit_square():
    return np.array([
        [0.0, 1.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 1.0],
    ])


def distorted_quad():
    return np.array([
        [0.0, 1.1, 1.3, -0.2],
        [0.0, -0.1, 1.2, 0.9],
    ])


def assemble(problem, coords, U, DU, dtime=0.1):
    return assemble_element(
        problem, coords, U, DU, dtime, props_array(), element="quad4")


def set_T(U, values):
    for dof, val in zip(T_DOFS, values):
        U[dof] = val
    return U


def check():
    from build import HeatDiffusionProblem

    problem = HeatDiffusionProblem()

    # 1a. Undeformed, T = 0 everywhere: zero residual.
    RHS0, _ = assemble(problem, unit_square(), np.zeros(12), np.zeros(12))
    if np.max(np.abs(RHS0)) > 1e-12:
        raise AssertionError("nonzero RHS at reference state: {}".format(RHS0))
    print("[PASS] zero residual at the reference state")

    # 1b. Undeformed, uniform steady T: the element carries the uniform
    # thermal stress P = -K alpha T I, whose nodal forces on the unit
    # square are hand-derivable from the shape-gradient integrals
    # int grad(N_a) dV = (+-1/2, +-1/2). With RHS = -f_int,
    #   RHS_mech = K alpha T * [-.5,-.5, .5,-.5, .5,.5, -.5,.5],
    # and the thermal rows stay zero (steady, uniform field).
    Tval = 3.0
    U = set_T(np.zeros(12), [Tval] * 4)
    RHS, AMATRX = assemble(problem, unit_square(), U, np.zeros(12))
    from check_reference import ALPHA, K
    g = np.array([-0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5])
    expected_mech = K * ALPHA * Tval * g
    if not np.allclose(RHS[U_DOFS], expected_mech, rtol=0, atol=1e-12):
        raise AssertionError(
            "thermal-stress nodal forces mismatch:\nobserved={}\n"
            "expected={}".format(RHS[U_DOFS], expected_mech))
    if np.max(np.abs(RHS[T_DOFS])) > 1e-10:
        raise AssertionError("thermal rows nonzero at steady state")
    print("[PASS] closed-form thermal-stress nodal forces (uniform T)")

    # 2. Rigid translation: the residual is UNCHANGED (the thermal-stress
    # forces of check 1b persist; a translation must add nothing).
    U2 = U.copy()
    for dof in U_DOFS[0::2]:
        U2[dof] += 0.37          # uniform x-translation
    RHS2, _ = assemble(problem, unit_square(), U2, np.zeros(12))
    if not np.allclose(RHS2, RHS, rtol=0, atol=1e-10):
        raise AssertionError(
            "rigid translation changed the residual:\n{}\nvs\n{}".format(
                RHS2, RHS))
    print("[PASS] rigid-translation invariance of the residual")

    # 3a. Distorted element, linear steady field: thermal-row sum is zero.
    coords = distorted_quad()
    a, b = 2.0, -0.7
    T_lin = [a * coords[0, n] + b * coords[1, n] for n in range(4)]
    U3 = set_T(np.zeros(12), T_lin)
    RHS3, _ = assemble(problem, coords, U3, np.zeros(12))
    thermal_sum = float(np.sum(RHS3[T_DOFS]))
    if abs(thermal_sum) > 1e-10:
        raise AssertionError(
            "flux term leaked into the heat balance: {}".format(thermal_sum))
    print("[PASS] heat balance on a distorted element (steady linear field)")

    # 3b. Uniform step dT on the unit square. The storage term enters the
    # residual R positively and RHS = -R, so the SIGNED thermal-row sum is
    #   sum(RHS[T]) = -rho_cp dT/dt * area  (negative for heating).
    dT, dt = 2.5, 0.2
    U4 = set_T(np.zeros(12), [dT] * 4)
    DU4 = set_T(np.zeros(12), [dT] * 4)
    RHS4, _ = assemble(problem, unit_square(), U4, DU4, dtime=dt)
    expected = -RHO_CP * dT / dt * 1.0
    observed = float(np.sum(RHS4[T_DOFS]))
    if abs(observed - expected) > 1e-9 * abs(expected):
        raise AssertionError(
            "signed storage balance mismatch: observed {} expected {}".format(
                observed, expected))
    print("[PASS] heat balance with uniform temperature step")

    # 4. Assembled tangent vs finite differences on a generic coupled state.
    rng_state = np.array([
        0.013, -0.008, 1.2, 0.021, 0.004, 0.7, -0.011, 0.017, -0.4,
        0.006, -0.015, 2.1,
    ])
    DU5 = 0.3 * rng_state
    RHS5, AMATRX5 = assemble(problem, distorted_quad(), rng_state, DU5)
    eps = 1e-7
    fd = np.zeros((12, 12))
    for j in range(12):
        Up = rng_state.copy()
        Up[j] += eps
        DUp = DU5.copy()
        DUp[j] += eps          # perturbing U at fixed old state
        RHSp, _ = assemble(problem, distorted_quad(), Up, DUp)
        fd[:, j] = (RHSp - RHS5) / eps
    rel = np.max(np.abs(AMATRX5 - (-fd))) / (np.max(np.abs(AMATRX5)) + 1e-30)
    if not np.isfinite(rel):
        raise AssertionError("non-finite tangent comparison")
    if rel > 5e-5:
        raise AssertionError(
            "assembled AMATRX vs -dRHS/dU failed: rel err {:.3e}".format(rel))
    print("[PASS] assembled AMATRX vs -dRHS/dU (rel err {:.1e})".format(rel))


if __name__ == "__main__":
    check()
