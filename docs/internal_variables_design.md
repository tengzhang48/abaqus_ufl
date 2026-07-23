# Internal Variables Design

## Architecture Decision

Internal variables (plastic strain, damage, hardening parameters) live
inside the material function, which the UEL calls directly at each Gauss
point. The UEL handles multi-physics coupling (weak form assembly). The
material function handles constitutive law + internal variable evolution.

**The UEL never calls a UMAT.** A UMAT is a separate output target
(for use with standard Abaqus elements) that wraps the same material
function with push-forward → Jaumann → Voigt → DDSDDE. The UEL
works in total Lagrangian with PK1 stress and dP/dF directly — no
push-forward, no DDSDDE.

```
UEL path:
  stress_PK1(F, p, mu, state_old) → CS engine → dP/dF, dP/dp, ...
                                                → assemble AMATRX

UMAT path (separate, for standard Abaqus elements):
  stress_PK1(F, state_old) → CS engine → dP/dF → push-forward
                                        → Jaumann → Voigt → DDSDDE
```

Same material function, different wrappers.

```
┌──────────────────────────────────────────────┐
│  UEL subroutine (Abaqus entry point)         │
│  - DOF parsing, field interpolation          │
│  - Gauss loop                                │
│  - SVARS management (read/write state)       │
│  - RHS + AMATRX assembly                     │
│                                              │
│  At each Gauss point:                        │
│    ┌──────────────────────────────────────┐   │
│    │  Material function (NOT a UMAT)      │   │
│    │  - Takes F (complex), state_old (real)│   │
│    │  - Return mapping / state update     │   │
│    │  - Returns P (complex), state_new    │   │
│    └──────────────────────────────────────┘   │
│                                              │
│  CS tangent engine:                          │
│    - Perturbs F, p, mu, grad_mu (complex)    │
│    - Calls same material function            │
│    - state_old is REAL (not perturbed)       │
│    - Extracts AIMAG → all tangent blocks     │
│    - Captures algorithmic tangent (not       │
│      continuum tangent)                      │
└──────────────────────────────────────────────┘
```

## Why This Works

### CS captures the algorithmic tangent

For plasticity with return mapping, the stress P depends on F through
the return mapping algorithm:

```
1. Trial elastic state: F_e_trial = F * (F_p_old)^{-1}
2. Check yield: f = vonMises(F_e_trial) - sigma_y
3. If f > 0: return mapping (Newton iteration for Δγ)
4. Update: F_p_new = exp(Δγ * n) * F_p_old
5. Stress: P = function(F * (F_p_new)^{-1})
```

When CS perturbs F → F + i·h·e_kl, the entire return mapping runs
in complex arithmetic. The Newton iteration converges on DBLE(residual)
(branching on real part is CS-safe). The resulting P is complex, and
AIMAG(P)/h gives the exact algorithmic tangent — consistent with the
return mapping, not the continuum tangent. This is a known advantage
of CS over analytical derivation.

### State variables are always real inputs

The key rule: **state_old comes from SVARS (real). It is never perturbed.**

During the CS perturbation for dP/dF:
- F is complex (perturbed)
- state_old is real (from SVARS, same for all perturbations)
- The return mapping produces complex state_new (because it depends on F)
- AIMAG(state_new)/h gives d(state)/dF (not needed for assembly, but
  available if needed)
- Only DBLE(state_new) from the REAL evaluation gets written to SVARS

```
Real evaluation:     F (real)    + state_old (real)  → P (real),    state_new (real) → SVARS
CS perturbation:     F (complex) + state_old (real)  → P (complex)  → AIMAG for tangent
                                                        state_new NOT written to SVARS
```

## Python Interface

```python
class J2Plasticity(au.Material):
    props = dict(E=200e3, nu=0.3, sigma_y=250.0, H=1000.0)
    
    # Declare state variables with initial values
    state_vars = dict(
        ep=0.0,                    # equivalent plastic strain (scalar)
        Fp=np.eye(3),              # plastic deformation gradient (tensor)
    )
    
    def stress_PK1(self, F, ep_old, Fp_old, dt):
        # Elastic trial
        Fp_old_inv = inv(Fp_old)
        Fe_trial = F @ Fp_old_inv
        
        # ... return mapping (all in complex arithmetic) ...
        
        # Return stress AND updated state
        return P, {'ep': ep_new, 'Fp': Fp_new}
```

### Key differences from current Material class

1. `state_vars` dict declares the internal variables and their initial
   values (scalar or tensor)
2. Material methods receive `xxx_old` arguments for each state variable
3. Material methods return a tuple: `(output, state_dict)`
4. `state_old` arguments are always REAL (DOUBLE PRECISION in Fortran)
5. `dt` is REAL (not differentiated)

### The WeakForm doesn't change

```python
class MechProblem(au.WeakForm):
    material = J2Plasticity
    
    def define_fields(self):
        self.u = au.VectorField('u', degree=1)
    
    def momentum_equation(self, v, F):
        P, state_new = self.material.stress_PK1(F)
        return P
```

The WeakForm only sees P (the stress). The state variable management is
entirely inside the material subroutine and the UEL's SVARS handling.

## Generated Fortran Structure

### SVARS layout

```
SVARS = [GP1_state, GP2_state, ..., GPn_state]

GP_state = [ep, Fp_11, Fp_12, Fp_13, Fp_21, Fp_22, Fp_23, 
            Fp_31, Fp_32, Fp_33]
```

For the J2 example: 1 scalar + 9 tensor components = 10 values per GP.
For 4 GPs (Quad4): NSVARS = 40. For 9 GPs (Quad8): NSVARS = 90.

### Generated SVARS code in UEL

```fortran
C     Read state at Gauss point kk
      LOC = (kk - 1) * NSTATE_PER_GP
      ep_old = SVARS(LOC + 1)
      DO i = 1, 3
        DO j = 1, 3
          Fp_old(i,j) = SVARS(LOC + 1 + (i-1)*3 + j)
        END DO
      END DO

C     Material evaluation (real, for RHS + SVARS update)
      CALL mat_stress_PK1(Fbar_z, ep_old_z, Fp_old_z,
     &  dt_safe, PROPS, Pz_eval, ep_new_z, Fp_new_z)
      P_bar = DBLE(Pz_eval)

C     Write updated state to SVARS (real part only)
      SVARS(LOC + 1) = DBLE(ep_new_z)
      DO i = 1, 3
        DO j = 1, 3
          SVARS(LOC + 1 + (i-1)*3 + j) = DBLE(Fp_new_z(i,j))
        END DO
      END DO

C     CS tangent (state_old is real, NOT perturbed)
      DO k = 1, 3
        DO l = 1, 3
          CALL real2complex33(Fbar, Fbar_z)
          ep_old_z = DCMPLX(ep_old, 0.0d0)
          CALL real2complex33(Fp_old, Fp_old_z)
          Fbar_z(k,l) = Fbar_z(k,l) + DCMPLX(0.0d0, CS_H)
          CALL mat_stress_PK1(Fbar_z, ep_old_z, Fp_old_z,
     &      dt_safe, PROPS, Pz, ep_z, Fp_z)
          dPdF(:,:,k,l) = AIMAG(Pz) / CS_H
C         Note: ep_z and Fp_z contain d(state)/dF in imaginary part
C         but we don't need them for assembly — just discard
        END DO
      END DO
```

### Generated material subroutine signature

```fortran
      SUBROUTINE mat_stress_PK1(F, ep_old, Fp_old, dt,
     &  props, P_out, ep_new, Fp_new)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: F(3,3)
      DOUBLE COMPLEX, INTENT(IN)  :: ep_old       ! state: real input
      DOUBLE COMPLEX, INTENT(IN)  :: Fp_old(3,3)  ! state: real input
      DOUBLE PRECISION, INTENT(IN) :: dt
      DOUBLE PRECISION, INTENT(IN) :: props(4)
      DOUBLE COMPLEX, INTENT(OUT) :: P_out(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: ep_new       ! state: complex output
      DOUBLE COMPLEX, INTENT(OUT) :: Fp_new(3,3)  ! state: complex output
```

Note: `ep_old` and `Fp_old` are declared DOUBLE COMPLEX but are passed
real values (zero imaginary part). This is because during CS perturbation,
F is complex and the return mapping produces complex state updates. The
interface must be uniform — all field-dependent arguments are COMPLEX.

The `dt` stays DOUBLE PRECISION (never differentiated).

## Standalone UMAT Generation (Separate Target)

The same material function can generate a standard Abaqus UMAT for use
with built-in elements. This is a **separate output**, not called by the UEL.

```python
model = J2Plasticity(E=200e3, nu=0.3, sigma_y=250, H=1000)
model.verify()
generate_umat(model, 'j2_umat.for')      # UMAT for built-in Abaqus elements
generate_uel(problem, 'coupled.for')      # UEL for multi-physics (calls material directly)
```

The UMAT wrapper adds:
- CS engine for dP/dF (9 perturbations)
- Push-forward: dP/dF → spatial tangent
- Jaumann correction: σ⊗I + I⊗σ terms
- Voigt packing: 4th-order tensor → 6×6 DDSDDE
- STATEV read/write (same layout as UEL's SVARS)

The UEL does NOT use any of this — it works directly with P and dP/dF
in the reference configuration. No push-forward, no Jaumann, no Voigt.

## What Changes in the Framework

### 1. Material class — return signature

Materials with internal variables return a tuple `(output, state_dict)`.
Materials without internal variables return just `output` (backwards
compatible).

```python
# Without state (existing, unchanged):
def stress_PK1(self, F, p, mu):
    return P

# With state (new):
def stress_PK1(self, F, ep_old, Fp_old, dt):
    ...
    return P, {'ep': ep_new, 'Fp': Fp_new}
```

Detection: if `hasattr(mat, 'state_vars')` and `len(mat.state_vars) > 0`,
the material has internal variables. Otherwise, backwards-compatible path.

### 2. WeakForm — unpack tuple returns

The WeakForm equation methods call `self.material.stress_PK1(...)` and
expect just `P`. For materials with state, the WeakForm wraps the call:

```python
def momentum_equation(self, v, F):
    result = self.material.stress_PK1(F)
    if isinstance(result, tuple):
        P, state_new = result
    else:
        P = result
    return P
```

Or more cleanly, the framework handles this in the tangent engine — the
material method is called, the first element of the tuple is used for
tangent computation, and the state_dict is ignored at the WeakForm level.

### 3. Fortran subroutine argument ordering

Consistent ordering across all generated subroutines:

```
1. Field variables (COMPLEX):    F, p, mu, grad_mu
2. State variables old (COMPLEX): ep_old, Fp_old, ...
3. Real parameters:              props, dt
4. Output (COMPLEX):             P_out, rp_out, jR_out, ...
5. State variables new (COMPLEX): ep_new, Fp_new, ...
```

State variables are DOUBLE COMPLEX (not DOUBLE PRECISION) because during
CS perturbation, the return mapping produces complex state updates even
though the inputs are real (zero imaginary part).

### 4. CS engine extension

The CS engine must pass state_old alongside field variables. State_old
is NOT perturbed — same real value for all perturbations:

```fortran
C     CS perturbation for dP/dF
      DO k = 1, 3
        DO l = 1, 3
C         Reset ALL complex variables
          CALL real2complex33(F, Fz)
          ep_old_z = DCMPLX(ep_old, 0.0d0)     ! real, not perturbed
          CALL real2complex33(Fp_old, Fp_old_z)  ! real, not perturbed
C         Perturb F only
          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)
          CALL mat_stress_PK1(Fz, ep_old_z, Fp_old_z,
     &      dt_safe, PROPS, Pz, ep_z, Fp_z)
          dPdF(:,:,k,l) = AIMAG(Pz) / CS_H
        END DO
      END DO
```

### 5. SVARS layout and size calculation

The generator computes `NSTATE_PER_GP` from `state_vars`:

| Python type | Size | Storage order |
|-------------|:----:|---------------|
| `float` (scalar) | 1 | single value |
| `np.zeros(3)` (vector) | 3 | components 1,2,3 |
| `np.eye(3)` (3×3 tensor) | 9 | column-major (Fortran native) |

Total SVARS = NSTATE_PER_GP × number_of_Gauss_points.

Column-major for tensors: `Fp(i,j)` stored at offset `(j-1)*3 + i`.
This matches Fortran's native memory layout and avoids transposition bugs.

### 6. verify() with nonzero initial state

`verify()` should test with both zero and nonzero initial state to catch
return mapping bugs:

```python
model = J2Plasticity(E=200e3, nu=0.3, sigma_y=250, H=1000)

# Default: verifies at initial state (ep=0, Fp=I)
model.verify()

# With nonzero state: verifies algorithmic tangent from return mapping
model.verify(state={'ep': 0.05, 'Fp': some_Fp_matrix})
```

If `state` is not provided, `verify()` uses the `state_vars` initial
values. If `state_vars` is not defined, the existing behavior is unchanged.

### 7. Backwards compatibility

Materials without `state_vars` (NeoHookean, GelMaterial) work exactly
as before. The generator checks:

```python
has_state = hasattr(mat, 'state_vars') and len(mat.state_vars) > 0
```

If `has_state` is False:
- No SVARS read/write generated
- Material subroutine has no state arguments
- CS engine has no state arguments
- NSVARS = 0 (or whatever the user sets for post-processing)

### Summary of changes

| Component | With state_vars | Without state_vars |
|-----------|----------------|-------------------|
| Material return | `(P, state_dict)` | `P` (unchanged) |
| WeakForm | Unpacks tuple | Unchanged |
| verify() | Tests at nonzero state | Unchanged |
| Fortran signature | F, state_old, props, dt, P_out, state_new | F, props, P_out (unchanged) |
| CS engine | Passes state_old (real) | Unchanged |
| UEL | SVARS read/write per GP | No SVARS code |
| UMAT | STATEV read/write | No STATEV code |

## CS-Safety Rules for Internal Variables

From `complex_step_patterns.md`:

1. **State_old is always real** — read from SVARS with zero imaginary part
2. **Never force state back to real during CS** — `ep_new = DBLE(ep_new)` 
   destroys the algorithmic tangent
3. **Converge return mapping on DBLE(residual)** — branch on real part only
4. **State_new is complex** — its imaginary part carries d(state)/d(input)
5. **Only write DBLE(state_new) to SVARS** — during the real evaluation only

## When to Use UEL vs UMAT

| Problem | Generate | Why |
|---------|----------|-----|
| Single-field mechanics (NeoHookean, J2, Ogden) | UMAT | Abaqus built-in elements sufficient |
| Single-field + F-bar (nearly incompressible) | UEL (Quad4) | F-bar needs custom element |
| Two-field (u + μ, u + damage) | UEL | Multi-physics coupling |
| Three-field (u + p + μ) | UEL | Mixed interpolation |
| Any of above + internal variables | Same choice | Internal vars live in material function, same in both |
| Thermo-mechanical (2 fields) | UEL or UMAT+UMATHT | UEL is more general |

**Key point:** The material function (COMPLEX*16, with state variables)
is identical regardless of the target. The choice of UEL vs UMAT is about
the element technology and field coupling, not about the constitutive law.
