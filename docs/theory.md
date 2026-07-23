# Mathematical Framework

## 1. General Structure

All problems share this pattern: a deformable body with one or more
scalar fields transported through it.

**Solid mechanics (momentum balance):**

    Div(P) + b_0 = 0    in Omega_0

**Scalar field equations (transport/constraint):**

    storage_i + Div(j_i) = s_i    in Omega_0,    i = 1, ..., n_fields

where P is the first Piola-Kirchhoff stress, j_i are referential
fluxes, s_i are source/sink terms.

**Formulation choices:**
- Total Lagrangian (reference configuration for all integrals)
- PK-1 stress: P = d(psi_R)/dF — one differentiation from the energy
- No push-forward during element assembly, no Jaumann correction
- Mixed u-p formulation for near-incompressibility

## 2. Weak Form and Discretization

For each field equation, the weak form is obtained by multiplying by a
test function and integrating by parts.

**Momentum (test function v):**

    integral P : Grad(v) dV = integral t . v dS + integral b . v dV

**Algebraic constraint (test function q):**

    integral r_p * q dV = 0

**Transport (test function w):**

    integral c_dot * w dV - integral j_R . Grad(w) dV = -integral j_n * w dS

The user provides: P, r_p, c_dot, j_R (constitutive functions).
The framework provides: shape functions, quadrature, assembly.

## 3. Complex-Step Tangent Computation

All tangent blocks are computed via:

    d(f)/d(x) = Im[ f(x + ih) ] / h,    h = 1e-10

For the element Jacobian, this means perturbing each DOF value (or its
gradient) by ih at each Gauss point and extracting the imaginary part of
the residual contribution.

**Cost per Gauss point** (three-field gel model):

| Block | Perturbations | Evaluations |
|-------|:---:|:---:|
| dP/dF | 9 (3x3 components) | 9 stress |
| dP/dp | 1 | 1 stress |
| dP/dmu | 0 (known zero) | 0 |
| dr_p/dF | 9 | 9 constraint |
| dr_p/dp | 1 | 1 constraint |
| dr_p/dmu | 1 | 1 constraint |
| dj_R/dF | 9 | 9 flux |
| dj_R/dp | 1 | 1 flux |
| dj_R/dmu | 1 | 1 flux |
| dj_R/d(grad_mu) | 3 (in 3D) | 3 flux |
| dc_dot/dF | 9 | 9 storage |
| dc_dot/dp | 1 | 1 storage |
| **Total** | **45** | **45** |

## 4. Mixed Interpolation

For nearly incompressible materials, equal-order interpolation of all
fields leads to volumetric locking or pressure oscillations. The
framework uses mixed-order interpolation:

| Field | Interpolation | Nodes (Quad8) | Rationale |
|-------|:---:|:---:|---|
| u (displacement) | Quadratic (degree 2) | All 8 | Standard for accuracy |
| p (pressure) | Linear (degree 1) | 4 corners | Inf-sup stability |
| mu (chemical potential) | Quadratic (degree 2) | All 8 | Coupled to u |

This satisfies the LBB (inf-sup) condition for the mixed u-p formulation.

## 5. Application: Elastomeric Gel

**Fields:** displacement u, elastic pressure p, chemical potential mu

**Kinematics:**

    F = I + Grad(u),    J = det(F)
    Je = exp(p/K),      phi = phi0 * Je / J

**PK-1 stress:**

    P = G(F - F^{-T}) + p F^{-T}

**Pressure constraint:**

    mu = mu0 + RT[ln(1-phi) + phi + chi*phi^2] - (phi/phi0)*Omega*p

**Referential flux:**

    j_R = -(D*cR/RT) C^{-1} Grad(mu)

**Mass balance (exact backward difference):**

    c_dot = (J/Je - J_old/Je_old) / (Omega * dt)

## 6. References

1. Chester, S.A., Di Leo, C.V., Anand, L. (2015). Int. J. Solids Struct.
2. Mao, Y., Anand, L. (2018). J. Mech. Phys. Solids.
3. Datta, B., Nguyen, T.D. (2024). JHU Technical Report.
4. Korelc, J. (2002). Comput. Mech. (AceGen)
5. Alnaes, M.S., et al. (2014). ACM Trans. Math. Softw. (UFL)
