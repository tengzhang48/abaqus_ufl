"""Assembled-element checks for the mixed-order thermo Quad8 UEL.

The distinct capability here is the mixed-order layout: quadratic
displacement on all eight nodes, bilinear temperature on the four
corners, per-node interleaved (corners (u1, u2, T) -> DOFs 1..12,
midsides (u1, u2) -> DOFs 13..20). Checks:

1. zero residual at the reference state; residual invariance under a
   rigid translation of a heated state;
2. pulled-back flux element oracle: on the unit square, a homogeneous
   x-stretch l (exactly representable by the quadratic basis) with a
   steady corner-linear temperature field T = a X gives constant
   referential flux -kappa a / l^2 e_x, so each thermal row equals that
   flux contracted with the bilinear shape-gradient integrals
   (+-1/2, +-1/2); the l-run and the undeformed run must differ by the
   factor 1/l^2 EXACTLY, which discriminates the C^{-1} pull-back;
3. partition-of-unity heat balance on a distorted Quad8 (zero flux-term
   sum for a steady linear field; exact cT dT/dt x area for a uniform
   step);
4. assembled AMATRX vs -dRHS/dU over all 20 DOFs.
"""

from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from abaqus_ufl.core.reference_assembly import assemble_element

from check_reference import CT, KAPPA

NDOFEL = 20
T_DOFS = [2, 5, 8, 11]
UX_DOFS = [0, 3, 6, 9, 12, 14, 16, 18]


def props_array():
    from build import ThermoMechanicalMaterial
    mat = ThermoMechanicalMaterial()
    return np.array([getattr(mat, name) for name in mat.props])


def quad8_coords(corners):
    """Corners (2,4) -> full Quad8 coords (2,8) with edge midpoints."""
    corners = np.asarray(corners, dtype=float)
    mids = np.column_stack([
        0.5 * (corners[:, 0] + corners[:, 1]),
        0.5 * (corners[:, 1] + corners[:, 2]),
        0.5 * (corners[:, 2] + corners[:, 3]),
        0.5 * (corners[:, 3] + corners[:, 0]),
    ])
    return np.column_stack([corners, mids])


def unit_square():
    return quad8_coords([[0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0]])


def distorted_quad():
    return quad8_coords([[0.0, 1.1, 1.3, -0.2], [0.0, -0.1, 1.2, 0.9]])


def dof_indices(node):
    """(start, has_T) for 0-based node index in the interleaved layout."""
    if node < 4:
        return 3 * node, True
    return 12 + 2 * (node - 4), False


def build_U(coords, stretch_x=1.0, T_of_x=None):
    """Homogeneous x-stretch plus a corner-supported temperature field."""
    U = np.zeros(NDOFEL)
    for node in range(8):
        start, has_T = dof_indices(node)
        U[start] = (stretch_x - 1.0) * coords[0, node]
        if has_T and T_of_x is not None:
            U[start + 2] = T_of_x(coords[0, node], coords[1, node])
    return U


def assemble(problem, coords, U, DU, dtime=0.1):
    return assemble_element(
        problem, coords, U, DU, dtime, props_array(), element="quad8")


def check():
    from build import ThermoMechanicalQuad8

    problem = ThermoMechanicalQuad8()

    # 1. Reference state and rigid-translation invariance.
    RHS0, _ = assemble(problem, unit_square(), np.zeros(NDOFEL),
                       np.zeros(NDOFEL))
    if np.max(np.abs(RHS0)) > 1e-12:
        raise AssertionError("nonzero RHS at reference state")
    U_hot = build_U(unit_square(), 1.0, lambda x, y: 3.0)
    RHS_hot, _ = assemble(problem, unit_square(), U_hot, np.zeros(NDOFEL))
    U_tr = U_hot.copy()
    for dof in UX_DOFS:
        U_tr[dof] += 0.41
    RHS_tr, _ = assemble(problem, unit_square(), U_tr, np.zeros(NDOFEL))
    if not np.allclose(RHS_tr, RHS_hot, rtol=0, atol=1e-9):
        raise AssertionError("rigid translation changed the residual")
    print("[PASS] reference state + rigid-translation invariance")

    # 2. Pulled-back flux element oracle. Steady linear T = a X, exact
    # homogeneous stretch l. Thermal rows carry the constant referential
    # flux contracted with the bilinear gradient integrals g_a.
    a = 2.0
    g_ax = np.array([-0.5, 0.5, 0.5, -0.5])
    thermal_rows = {}
    for lam in (1.0, 1.3):
        U = build_U(unit_square(), lam, lambda x, y: a * x)
        RHS, _ = assemble(problem, unit_square(), U, np.zeros(NDOFEL))
        thermal_rows[lam] = RHS[T_DOFS].copy()
        # RHS = -R and R_T,a = -int grad(theta_a) . flux with the constant
        # referential flux flux_x = -kappa a / l^2, so
        #   RHS_T,a = -(kappa a / l^2) g_ax.
        expected_rows = -KAPPA * a / lam ** 2 * g_ax
        if not np.allclose(thermal_rows[lam], expected_rows,
                           rtol=1e-10, atol=1e-12):
            raise AssertionError(
                "thermal rows at l={} are {} (expected {})".format(
                    lam, thermal_rows[lam], expected_rows))
        if not np.allclose(np.abs(np.sum(thermal_rows[lam])), 0.0,
                           atol=1e-10):
            raise AssertionError("flux term leaked into heat balance")
    ratio = thermal_rows[1.3] / thermal_rows[1.0]
    if not np.allclose(ratio, 1.0 / 1.3 ** 2, rtol=1e-9):
        raise AssertionError(
            "pull-back factor missing in assembled flux: ratio {}".format(
                ratio))
    print("[PASS] pulled-back flux element oracle (ratio 1/l^2 exact)")

    # 3. Heat balance: steady linear field on a distorted element sums to
    # zero; a uniform step on the unit square sums to cT dT/dt x area.
    U3 = build_U(distorted_quad(), 1.0,
                 lambda x, y: 2.0 * x - 0.7 * y)
    RHS3, _ = assemble(problem, distorted_quad(), U3, np.zeros(NDOFEL))
    if abs(float(np.sum(RHS3[T_DOFS]))) > 1e-10:
        raise AssertionError("distorted-element heat balance failed")
    # Signed balance: the storage term enters R positively and RHS = -R,
    # so sum(RHS[T]) = -cT dT/dt * area (negative for heating).
    dT, dt = 2.5, 0.2
    U4 = build_U(unit_square(), 1.0, lambda x, y: dT)
    DU4 = U4.copy()
    RHS4, _ = assemble(problem, unit_square(), U4, DU4, dtime=dt)
    expected = -CT * dT / dt
    observed = float(np.sum(RHS4[T_DOFS]))
    if abs(observed - expected) > 1e-9 * abs(expected):
        raise AssertionError(
            "signed storage balance mismatch: observed {} expected {}".format(
                observed, expected))
    print("[PASS] heat balance (distorted steady + uniform step)")

    # 4. Assembled tangent vs finite differences over all 20 DOFs.
    rng = np.linspace(-0.02, 0.02, NDOFEL)
    U5 = build_U(distorted_quad(), 1.05, lambda x, y: 1.5 * x + 0.4 * y)
    U5 = U5 + rng
    DU5 = 0.3 * U5
    RHS5, AMATRX5 = assemble(problem, distorted_quad(), U5, DU5)
    eps = 1e-7
    fd = np.zeros((NDOFEL, NDOFEL))
    for j in range(NDOFEL):
        Up = U5.copy()
        Up[j] += eps
        DUp = DU5.copy()
        DUp[j] += eps
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
