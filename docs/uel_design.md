# UEL Code Generator Design

## Overview

`uel_gen.py` takes a `WeakForm` instance and produces a single self-contained
`.for` file that works as an Abaqus UEL subroutine. The user writes only Python.
No Fortran knowledge required.

```python
problem = GelProblem()
problem.verify()
import abaqus_ufl as au
au.generate_uel(problem, 'gel_uel.for', element='Quad8')
```

## What the generator produces

```
┌─────────────────────────────────────────────────┐
│ Section 1: Shipped Fortran templates            │
│   tensor_ops.for     (det33z, inv33z, ...)      │
│   shape_quad8.for    (shape functions)           │
│   shape_quad4.for    (for degree=1 fields)      │
│   gauss_rules.for    (3x3 volume, 3pt edge)     │
│   isoparametric.for  (map_grad_2d + apply_jinv) │
│   edge_quad8.for     (surface integration)      │
├─────────────────────────────────────────────────┤
│ Section 2: Material subroutines (COMPLEX*16)    │
│   AST-translated from Python Material methods   │
│   Helper methods (_phi etc.) as separate subs   │
├─────────────────────────────────────────────────┤
│ Section 3: CS tangent engine                    │
│   Perturbation loops for all tangent blocks     │
├─────────────────────────────────────────────────┤
│ Section 4: UEL subroutine                       │
│   DOF parsing (current + old from U and DU)     │
│   Gauss loop with correct sign assembly         │
│   SVARS management                              │
├─────────────────────────────────────────────────┤
│ Section 5: UVARM subroutine (optional)          │
└─────────────────────────────────────────────────┘
```

---

## Section 1: Shipped templates

Static Fortran files included verbatim. One notable detail:

**isoparametric.for** returns `Jinv` from `map_grad_2d` and provides an
`apply_jinv` routine. This is needed because degree=1 fields (pressure) use
Quad4 shape functions but must be mapped with the Quad8 geometry Jacobian.

```fortran
      SUBROUTINE map_grad_2d(dshxi, coords, nNode, dsh, detJ,
     &                       Jinv_out, stat)
C     Returns Jinv(2,2) for reuse with lower-degree fields
      ...
      END SUBROUTINE

      SUBROUTINE apply_jinv(dshxi, Jinv, nNode, dsh)
C     Map dshxi -> dsh using a pre-computed Jinv
      ...
      END SUBROUTINE
```

---

## Section 2: Material subroutine translation

The AST translator converts each Python Material method to a COMPLEX*16
Fortran subroutine.

### Supported operations (tensor vocabulary)
```
Python                  Fortran
──────                  ───────
det(A)                  det33z(Az)
inv(A)                  CALL inv33z(Az, Az_inv)
A.T                     CALL transpose33z(Az, Az_T)
A @ B                   CALL matmul33z(Az, Bz, Cz)
log(x)                  LOG(xz)
exp(x)                  EXP(xz)
sqrt(x)                 SQRT(xz)
trace(A)                Az(1,1)+Az(2,2)+Az(3,3)
self.G                  DCMPLX(props(1), 0.0d0)
+, -, *, /, **          same operators
```

### Method signatures → Fortran subroutine signatures

```
Python argument    Fortran type              Notes
───────────────    ────────────              ─────
F                  DOUBLE COMPLEX :: F(3,3)  deformation gradient
p                  DOUBLE COMPLEX :: p       pressure scalar
mu                 DOUBLE COMPLEX :: mu      chemical potential scalar
grad_mu            DOUBLE COMPLEX :: grad_mu(3)  gradient vector
F_old              DOUBLE COMPLEX :: F_old(3,3)  history (previous step)
p_old              DOUBLE COMPLEX :: p_old   history scalar
dt                 DOUBLE PRECISION :: dt    time increment (always real)
```

Note: `dt` is DOUBLE PRECISION (not complex) because we never differentiate
with respect to it.

### Return types → Fortran output arguments

```
Method              Returns         Fortran output
──────              ───────         ──────────────
stress_PK1          P (3x3 tensor)  DOUBLE COMPLEX :: P_out(3,3)
pressure_resid      r_p (scalar)    DOUBLE COMPLEX :: rp_out
solvent_flux        j_R (3-vector)  DOUBLE COMPLEX :: jR_out(3)
solvent_storage     c_dot (scalar)  DOUBLE COMPLEX :: cdot_out
```

### State variable initialization

The `state_vars` dict declared at class level is used by the generator to
determine SVARS layout. However, if a prop like `phi0` controls the initial
value, the class-level `state_vars` may be out of sync. Always set
`self.state_vars` in `__init__` to match the actual props:

```python
def __init__(self, **kwargs):
    super().__init__(**kwargs)
    self.state_vars = {'phi': float(self.phi0)}
```

### Scaling coupled equations

For coupled diffusion-deformation at SI scale, the mechanical and chemical
stiffness blocks can differ by >10 orders. The solver fails. Add a `scale`
property and multiply transport outputs:

```python
props = dict(G=1e5, K=1e7, D=5e-11, scale=1e9)

def solvent_flux(self, F, mu, grad_mu, ...):
    return -m * Cinv @ grad_mu * self.scale

def solvent_storage(self, F, F_old, mu, ...):
    return c_dot * self.scale
```

This row-scales the Jacobian without changing the physical solution.

### Helper methods

The translator scans the Material class for all methods called by the
registered material methods:

1. For each registered method body, collect all `self._xxx(...)` calls via AST
2. Translate each helper as a separate Fortran subroutine
3. Emit helpers BEFORE the methods that call them
4. `self._phi(F, p)` becomes `CALL mat_phi(Fz, pz, props, phi_out)`

Return type is inferred by calling the helper with test inputs. If inference
fails, a warning is emitted and scalar is assumed.

---

## Section 3: CS tangent engine

Generated from the WeakForm's tangent block list. For the three-field gel:

```
Block                  Perturb   Output shape   Perturbations
─────                  ───────   ────────────   ─────────────
dP/dF                  F(k,l)    (3,3,3,3)      9
dP/dp                  p         (3,3)          1
dP/dmu                 mu        (3,3)          1
drp/dF                 F(k,l)    (3,3)          9
drp/dp                 p         scalar         1
drp/dmu                mu        scalar         1
djR/dF                 F(k,l)    (3,3,3)        9
djR/dp                 p         (3,)           1
djR/dmu                mu        (3,)           1
djR/dgrad_mu           grad_mu(k)(3,3)          3
dcdot/dF               F(k,l)    (3,3)          9
dcdot/dp               p         scalar         1
                                         Total: 46
```

Key properties:
- CS_H = 1e-10 (verified identical to 1e-30 across all blocks; avoids
  subnormal issues in Fortran intermediate computations)
- Full state reset before every perturbation block (prevents cross-contamination)
- Uses AIMAG (portable), DCMPLX, DBLE only — no non-standard intrinsics
- Tangent blocks sorted by stable key for deterministic argument order

---

## Section 4: UEL subroutine

### 4.0 Array initialization and safety guards

**Zero RHS and AMATRX**: Abaqus documentation warns not to assume
arrays are pre-initialized.

**DTIME guard**: Abaqus calls UEL with DTIME=0 for initial
stiffness evaluation. Division by zero in storage rate would crash.

```fortran
C     Zero RHS and AMATRX
      DO i = 1, NDOFEL
        RHS(i, 1) = 0.0d0
        DO j = 1, NDOFEL
          AMATRX(i, j) = 0.0d0
        END DO
      END DO

C     DTIME guard
      dt_safe = DTIME
      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14
```

### 4.1 DOF parsing — current AND old state + element-local DOF maps

**F_old from U - DU**: Abaqus does not pass previous-step fields.
The UEL receives `U(NDOFEL)` (total DOFs at end of increment) and
`DU(MLVARX,1)` (incremental DOFs). Old-step values are `U - DU`.

**edof mapping arrays**: Assembly loops need element-local DOF
indices (1 to NDOFEL) to index into RHS and AMATRX. The `edof` arrays
are built during DOF extraction — the `idx` counter tracks the 1-based position.

```fortran
C     Parse current + old DOFs, build element-local DOF maps
      idx = 0
      DO a = 1, NNODE
C       --- displacement (degree=2, all nodes) ---
        DO i = 1, ndim
          idx = idx + 1
          u_node(i, a) = U(idx)
          u_old(i, a)  = U(idx) - DU(idx, 1)
          edof_u(i, a) = idx
        END DO
C       --- pressure (degree=1, corner nodes only) ---
        IF (a .LE. NCORNER) THEN
          idx = idx + 1
          p_node(a) = U(idx)
          p_old_node(a) = U(idx) - DU(idx, 1)
          edof_p(a) = idx
        END IF
C       --- chemical potential (degree=2, all nodes) ---
        idx = idx + 1
        mu_node(a) = U(idx)
        mu_old_node(a) = U(idx) - DU(idx, 1)
        edof_mu(a) = idx
      END DO
```

**Key:** This loop is GENERATED from `WeakForm._field_order`. No hardcoded
field names in the generator itself.

### 4.2 Gauss loop structure

**Degenerate Jacobian trap**: If detJ ≤ 0 (inverted element),
set PNEWDT < 1 and return immediately.

**Shared Jinv**: Compute Jacobian once from 8-node geometry,
apply to both Quad8 and Quad4 shape function derivatives.

```fortran
      CALL gauss_2d_3x3(xi_gp, w_gp, ngp)

      DO kk = 1, ngp

C       Shape functions + geometry Jacobian
        CALL shape_quad8(xi_gp(kk,1), xi_gp(kk,2), sh8, dshxi8)
        CALL map_grad_2d(dshxi8, coords, 8, dsh8, detJ, Jinv, stat)

C       Trap degenerate Jacobian
        IF (stat .EQ. 0) THEN
          PNEWDT = 0.25d0
          RETURN
        END IF

        CALL shape_quad4(xi_gp(kk,1), xi_gp(kk,2), sh4, dshxi4)
        CALL apply_jinv(dshxi4, Jinv, 4, dsh4)

        wdetJ = detJ * w_gp(kk)

C       Field interpolation ...
C       Material evaluation ...
C       CS tangent engine ...
C       RHS assembly ...
C       AMATRX assembly ...
      END DO
```

### 4.3 Field interpolation

#### F₃₃ = 1 for plane strain

Initialize F and F_old as 3×3 identity before adding the 2D displacement
gradient. Otherwise F(3,3) = 0 and det(F) = 0, crashing inv(F).

```fortran
      CALL eye33d(F)
      CALL eye33d(F_old)
      DO a = 1, 8
        DO i = 1, ndim      ! ndim = 2
          DO J = 1, ndim
            F(i,J)     = F(i,J)     + u_node(i,a) * dsh8(a,J)
            F_old(i,J) = F_old(i,J) + u_old(i,a)  * dsh8(a,J)
          END DO
        END DO
      END DO
```

#### Scalar fields

Degree determines shape function set:
- degree=2 → `sh8`, `dsh8` (Quad8, 8 nodes)
- degree=1 → `sh4`, `dsh4` (Quad4, 4 corner nodes, mapped with shared Jinv)

Gradients are computed only for fields that appear as `grad_X` in equation
signatures.

### 4.4 Residual (RHS) assembly — with sign tracking

**Sign convention**: Each assembly term has a sign that enters RHS.
Transport equation returns (c_dot, j_R) with different signs.

```
Term type          RHS formula                           RHS sign
─────────          ───────────                           ────────
stress (P)         RHS(row,1) -= P_iJ * dN_a/dX_J * wdetJ   -1
constraint (r_p)   RHS(row,1) -= r_p * N_a * wdetJ           -1
storage (c_dot)    RHS(row,1) -= c_dot * N_a * wdetJ         -1
flux (j_R)         RHS(row,1) += j_R_J * dN_a/dX_J * wdetJ  +1
```

Note: RHS uses 2D indexing `RHS(row, 1)` per Abaqus convention.

### 4.5 Tangent (AMATRX) assembly — with sign tracking

**AMATRX = -d(RHS)/d(U)**: Since RHS terms have different signs,
the AMATRX contribution sign is the NEGATIVE of the RHS sign.

Complete sign table:

```
Tangent block    From term     RHS sign  AMATRX sign  Pattern
─────────────    ─────────     ────────  ───────────  ───────
dP/dF            stress          -1         +1        +grad-grad
dP/dp            stress          -1         +1        +grad-value
dP/dmu           stress          -1         +1        +grad-value
drp/dF           constraint      -1         +1        +value-grad
drp/dp           constraint      -1         +1        +value-value
drp/dmu          constraint      -1         +1        +value-value
dcdot/dF         storage         -1         +1        +value-grad
dcdot/dp         storage         -1         +1        +value-value
djR/dF           flux            +1         -1        -grad-grad
djR/dp           flux            +1         -1        -grad-value
djR/dmu          flux            +1         -1        -grad-value
djR/dgrad_mu     flux            +1         -1        -grad-grad
```

The 4 assembly patterns:

```
Pattern          Formula
───────          ─────────────────────────────────────────────
grad-grad        sign * sum_{J,L} T_iJkL * dN_a/dX_J * dN_b/dX_L * wdetJ
grad-value       sign * sum_J     T_iJ   * dN_a/dX_J * N_b       * wdetJ
value-grad       sign * sum_L     T_kL   * N_a       * dN_b/dX_L * wdetJ
value-value      sign * T        * N_a   * N_b       * wdetJ
```

The generator determines the pattern from two pieces of metadata:
- Row side: equation's assembly type (`grad` or `value`)
- Column side: differentiated variable's type (`matrix`→grad, `scalar`→value, `vector`→grad)

### 4.6 Fortran compatibility

**No BLOCK construct:** All complex temporaries declared at subroutine top
(F90-compatible, no Fortran 2008 features required).

**Column 72 wrapping:** All generated lines are wrapped at column 72 using
`_wrap_fortran()`. The wrapper searches backwards for break characters
(`' '`, `','`, `'+'`, `'-'`, `'*'`, `'('`, `')'`) and inserts `&` continuation
markers. Comment detection checks column 1 only (not after indentation — a
bug was found and fixed where variables starting with `c` like `cR0` were
misidentified as comments).

---

## Summary of what the generator needs from WeakForm

```python
weakform.ndim              # 2 or 3
weakform.field_names       # ['u', 'p', 'mu']
weakform.fields            # {name: VectorField/ScalarField}
weakform.ndofel            # 28
weakform._dof_map          # (node, field) -> [global DOF indices]
weakform._field_shape      # {field: 'quad8' or 'quad4'}
weakform._field_nodes      # {field: [node indices]}
weakform.equations         # {eq_name: {test_field, field_vars, ...}}
weakform.tangent_blocks()  # OrderedDict{(eq_name, var): {assembly, sign, ...}}
weakform._mat              # Material instance (for AST translation + props)
```

---

## Verification

The Python reference assembly (`reference_assembly.py`) mirrors the intended
UEL operations — DOF parsing, field interpolation, CS tangents, and assembly
patterns. Its FD tangent check confirms:

```
AMATRX ≈ -dRHS/dU    (max relative error: 5.12e-07)
```

This is Python assembled-residual/tangent evidence only. It neither executes
the generated Fortran nor proves the weak form matches the theory. Complete the
code-generation gate by calling the generated UEL through f2py or another
checked compiled runtime and comparing `RHS`, `AMATRX`, `SVARS`, and `PNEWDT`
with this assembly. Keep an independent quantitative physics oracle as a
separate gate.

---

## F-bar for coupled multi-physics

F-bar stabilization modifies the deformation gradient and adds a rank-1
correction to the tangent. For pure mechanical problems it eliminates
volumetric locking. For coupled problems at SI scale, it is no longer the
preferred production route.

**When it works:** Non-SI parameters (G~1, K~100, h~1) — F-bar converges
normally.

**When it fails:** SI-scale coupled problems (G~1e5, K~1e7, h~1e-4) — the
F-bar correction adds entries ~K/h² that dominate the standard stiffness,
producing condition numbers >1e15 and solver failure.

**Current recommendation:**
- Pure mechanics (any scale): F-bar is safe and recommended
- Coupled swelling/diffusion on Quad4/Hex8: use local pressure condensation
  when the pressure equation is algebraic. Store one pressure value per
  element in `SVARS`, solve it by local Newton, and statically condense it
  with all coupled tangent blocks included.
- Coupled high-order/reference studies: global three-field (u, p, μ) remains
  useful, especially for comparing against the mixed formulation.
- Disabling F-bar (`fbar=False`) is only a fallback for diagnosis because it
  converges but locks when K/G is large.
- The generated F-bar tangent has been FD-verified mathematically; the
  failure is a solver-conditioning issue, not a tangent bug

The local-pressure generator path (`abaqus_ufl/generators/uel_local_pressure.py`,
exposed as `au.generate_uel_local_pressure`) emits a UEL that stores `SVARS(1)`
as the condensed pressure state, can mirror `phi` and pressure into extra
`SVARS` slots, and emits `UVARM1/UVARM2` for Abaqus visualization.

### Why the coupled F-bar path is fragile

The current evidence points to the chemical-mechanical coupling as the source
of the F-bar difficulty, not to mechanical F-bar itself.

In pure mechanics, F-bar only modifies the volumetric part of the mechanical
deformation gradient. The residual and tangent depend on the corrected
deformation gradient through stress, so the chain rule is contained inside the
mechanical stiffness block. This has been verified by finite-difference checks,
and by a near-incompressible Neo-Hookean Quad4 cantilever (`K/G=1000`): the
F-bar element produces much larger load-control tip deflections than the
standard Quad4, as expected for a locking-relieved element.

In a swelling gel, the same volumetric quantity also controls non-mechanical
physics:

- polymer volume fraction `phi`,
- chemical-potential residual,
- transient solvent storage,
- mobility and flux,
- old-state terms through `F_old` and `p_old`.

Therefore a coupled F-bar formulation has two unattractive choices:

1. Apply F-bar only to the mechanical stress. This relieves mechanical locking
   but gives an inconsistent coupled tangent because transport still sees the
   physical `F` and `J`.
2. Apply F-bar consistently to every quantity that depends on volume. This
   requires chain-rule terms through `Jbar`, `phi`, storage, flux, and all
   `u-mu` tangent blocks, which is complicated and poorly conditioned at SI
   scale.

For this package, the design decision is:

- keep F-bar as the low-order mechanical formulation,
- keep coupled F-bar as an experimental/comparison branch,
- use local pressure condensation or global mixed `u,p,mu` formulations for
  production coupled swelling/diffusion.

## Deferred to follow-up

- Surface BC integration (NDLOAD/JDLTYP face mapping)
- General SVARS read/write for arbitrary post-processing variables
- Time-dependent boundary conditions
- Full original-package-style dummy-mesh UVARM/global-state visualization
- LFLAGS check for initial stiffness call
