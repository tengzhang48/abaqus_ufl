# UINTER Generator Design Document

## 1. Overview

`uinter_gen.py` generates complete Fortran UINTER subroutines from Python
surface interaction definitions. The generated `.for` file implements
user-defined interfacial constitutive behavior for Abaqus/Standard
contact pairs.

The architecture mirrors `umat_gen.py`:

```
Python SurfaceInteraction class
  → _generate_traction_subroutine()    [COMPLEX*16 traction law]
  → _generate_cs_engine()              [CS tangent for all Jacobian blocks]
  → _generate_uinter_wrapper()         [Abaqus UINTER interface]
  → generate_uinter()                  [assemble .for file]
```

### What the User Defines

```python
class MyCohesive(au.SurfaceInteraction):
    props = dict(K_n=1e6, K_t=1e5, sigma_c=1e3)

    state_vars = dict(d=0.0, H=0.0)  # optional

    def traction(self, rdisp, d_old, H_old, dt):
        """
        Args:
            rdisp: relative displacement vector (NDIR components)
                   rdisp[0] = normal (positive = penetration)
                   rdisp[1:] = tangential
            d_old, H_old: state from previous increment
            dt: time increment

        Returns:
            stress: traction vector (NDIR components)
                    positive normal = compression
            state_new: dict of updated state variables
        """
        ...
        return stress, {'d': d_new, 'H': H_new}

    def flux(self, rdisp, temp, d_old, dt):  # optional
        """
        Args:
            rdisp: relative displacement
            temp: [T_slave, T_master]
            d_old: state from previous increment
            dt: time increment

        Returns:
            flux: [q_slave, q_master]
                  positive = heat into surface
        """
        ...
        return flux
```

### What Gets Generated

A self-contained `.for` file containing:

1. **UINTER wrapper** — standard Abaqus signature, STATEV read/write,
   calls the CS engine, packs outputs
2. **Traction subroutine** — COMPLEX\*16, translated from Python
3. **Flux subroutine** — COMPLEX\*16, translated from Python (if defined)
4. **CS tangent engine** — computes DDSDDR, DDFDDT, DDSDDT, DDFDDR
5. **tensor\_ops.for** — shared math utilities

---

## 2. Abaqus UINTER Interface

### 2.1 Fortran Signature

```fortran
SUBROUTINE UINTER(STRESS,DDSDDR,DVISCOUS,DSTRUCTURAL,
     1   FLUX,DDFDDT,DDSDDT,DDFDDR,
     2   STATEV,SED,SFD,SPD,SVD,SCD,PNEWDT,
     3   RDISP,DRDISP,TEMP,DTEMP,PREDEF,DPRED,
     4   TIME,DTIME,FREQR,CINAME,SLNAME,MSNAME,
     5   PROPS,COORDS,ALOCALDIR,DROT,AREA,CHRLNGTH,
     6   NODE,NDIR,NSTATV,NPRED,NPROPS,MCRD,
     7   KSTEP,KINC,KIT,LINPER,LOPENCLOSE,LSTATE,
     8   LSDI,LPRINT)
```

### 2.2 Variables to Define

| Variable | Dimension | Description |
|----------|-----------|-------------|
| STRESS | (NDIR) | Interface traction. Positive normal = compression. Passed in as old value, must be updated to end-of-increment value. |
| DDSDDR | (NDIR,NDIR) | $\partial \text{STRESS}(i) / \partial \text{RDISP}(j)$ — interface stiffness matrix |
| FLUX | (2) | Heat flux into [slave, master] surfaces. Positive = into surface. |
| DDFDDT | (2,2) | $-\partial \text{FLUX}(i) / \partial \text{TEMP}(j)$ — note the negative sign |
| DDSDDT | (NDIR,2) | $\partial \text{STRESS}(i) / \partial \text{TEMP}(j)$ — thermo-mechanical coupling |
| DDFDDR | (2,NDIR) | $\partial \text{FLUX}(i) / \partial \text{RDISP}(j)$ — mechano-thermal coupling |

### 2.3 Variables That Can Be Updated

| Variable | Description |
|----------|-------------|
| STATEV(NSTATV) | Solution-dependent state variables |
| SED, SFD, SPD, SVD, SCD | Energy dissipation quantities (for output only) |
| PNEWDT | Suggested time step ratio (< 1 forces cutback) |
| LOPENCLOSE | Contact status flag: 0=open, 1=closed. Changes trigger SDI. |
| LSTATE | General state flag (user-defined integer) |
| LSDI | Set to 1 to force severe discontinuity iteration |

### 2.4 Variables Passed In for Information

| Variable | Dimension | Description |
|----------|-----------|-------------|
| RDISP | (NDIR) | Current relative position. RDISP(1) > 0 = penetration, < 0 = open gap. |
| DRDISP | (NDIR) | Increment in relative position |
| TEMP | (2) | Temperature at [slave, master] surface points |
| DTEMP | (2) | Temperature increment |
| PREDEF | (2,NPRED) | Predefined field variables at [slave, master] |
| DPRED | (2,NPRED) | Increments in predefined fields |
| TIME | (2) | [step time, total time] |
| DTIME | | Current time increment |
| PROPS | (NPROPS) | User-defined property values |
| COORDS | (MCRD) | Current coordinates of the contact point |
| ALOCALDIR | (3,3) | Local coordinate system (column 1 = normal, 2,3 = tangential) |
| DROT | (2,2) | Rotation increment matrix (for 3D rigid surface contact) |
| AREA | | Surface area associated with contact point |
| CHRLNGTH | | Characteristic contact surface face dimension |
| NODE | | Slave node number |
| NDIR | | Number of traction components (2 in 2D, 3 in 3D) |
| KSTEP, KINC | | Step and increment numbers |
| KIT | | Iteration number |
| LINPER | | Linear perturbation flag |
| LPRINT | | 1 if detailed contact printout requested |

### 2.5 Key Conventions

- **Sign convention for STRESS**: positive normal stress = compression
  (opposite to usual continuum mechanics convention). This matches
  Abaqus's contact pressure convention (CPRESS > 0 in compression).
- **Sign convention for RDISP(1)**: positive = penetration into master
  surface. Negative = open gap. For open points with no master pairing,
  RDISP(1) = $-10^{36}$.
- **Sign convention for DDFDDT**: the **negative** of $\partial$FLUX/$\partial$TEMP.
  This is Abaqus's convention for thermal conductance matrices.
- **NDIR**: 2 for 2D/axisymmetric (normal + 1 tangential), 3 for 3D
  (normal + 2 tangential).
- **STATEV initialization**: use the `TIME(2) == 0` init pattern (the
  same one used in `umat_gen.py`).
- **LOPENCLOSE**: set to $-1$ at analysis start. Setting to 0 (open) or
  1 (closed) enables contact status tracking. Changes trigger SDI
  automatically.

---

## 3. Complex-Step Tangent Engine for UINTER

### 3.1 What Needs to Be Differentiated

The UINTER must return four Jacobian blocks:

1. **DDSDDR(NDIR,NDIR)**: $\frac{\partial \text{STRESS}_i}{\partial \text{RDISP}_j}$
   — perturb each RDISP component, read imaginary part of STRESS

2. **DDFDDT(2,2)**: $-\frac{\partial \text{FLUX}_i}{\partial \text{TEMP}_j}$
   — perturb each TEMP component, read imaginary part of FLUX,
   negate (Abaqus convention)

3. **DDSDDT(NDIR,2)**: $\frac{\partial \text{STRESS}_i}{\partial \text{TEMP}_j}$
   — from the same TEMP perturbation, read imaginary part of STRESS

4. **DDFDDR(2,NDIR)**: $\frac{\partial \text{FLUX}_i}{\partial \text{RDISP}_j}$
   — from the same RDISP perturbation, read imaginary part of FLUX

### 3.2 Perturbation Strategy

**Phase 1: Perturb RDISP** ($j = 1, \ldots, \text{NDIR}$)

```
For j = 1 to NDIR:
    rdisp_z = DCMPLX(RDISP, 0)          ! reset to real
    state_old_z = DCMPLX(state_old, 0)   ! reset state (never perturbed)
    rdisp_z(j) = rdisp_z(j) + ih         ! perturb component j
    
    call traction(rdisp_z, state_old_z, dt, props, stress_z, state_new_z)
    call flux(rdisp_z, temp_z_real, state_old_z, dt, props, flux_z)  ! if flux defined
    
    DDSDDR(:,j) = AIMAG(stress_z) / h
    DDFDDR(:,j) = AIMAG(flux_z) / h       ! if flux defined
```

**Phase 2: Perturb TEMP** ($j = 1, 2$) — only if flux is defined

```
For j = 1, 2:
    rdisp_z = DCMPLX(RDISP, 0)          ! reset
    temp_z = DCMPLX(TEMP, 0)            ! reset
    state_old_z = DCMPLX(state_old, 0)
    temp_z(j) = temp_z(j) + ih          ! perturb slave (j=1) or master (j=2)
    
    call traction(rdisp_z, state_old_z, dt, props, stress_z, state_new_z)
    call flux(rdisp_z, temp_z, state_old_z, dt, props, flux_z)
    
    DDSDDT(:,j) = AIMAG(stress_z) / h    ! if traction depends on temp
    DDFDDT(:,j) = -AIMAG(flux_z) / h     ! note: NEGATIVE (Abaqus convention)
```

### 3.3 State Variable Handling in CS

Same rules as `umat_gen.py` and `uel_gen.py`:

- `state_old` is DOUBLE PRECISION, never perturbed
- Before each perturbation, convert to COMPLEX with zero imaginary part
- `state_new` from CS calls is COMPLEX but only used for its effect on
  STRESS/FLUX — the imaginary part carries the derivative information
- The actual state update uses the REAL (unperturbed) call only
- `dt` is DOUBLE PRECISION, never differentiated

### 3.4 Efficiency

For a 3D problem (NDIR=3) with both traction and flux:
- Phase 1: 3 perturbations → 3 traction + 3 flux calls
- Phase 2: 2 perturbations → 2 traction + 2 flux calls
- Total: **5 complex-arithmetic evaluations** per contact point

For 2D (NDIR=2): 4 evaluations. Without thermal coupling: NDIR evaluations.

This is trivially cheap — a single contact point evaluation costs
microseconds. The entire UINTER call (including CS) is negligible
compared to one element stiffness assembly.

---

## 4. Generated File Structure

```
┌─────────────────────────────────────────────────┐
│  Header (material name, props, state vars)       │
├─────────────────────────────────────────────────┤
│  Section 1: UINTER wrapper                       │
│    - STATEV read (TIME(2)==0 init pattern)       │
│    - Call CS engine → DDSDDR, DDFDDT, etc.       │
│    - Call real traction → STRESS, state_new      │
│    - STATEV write (DBLE of state_new)            │
│    - LOPENCLOSE / LSDI flag management           │
│    - PNEWDT safety                               │
├─────────────────────────────────────────────────┤
│  Section 2: Traction subroutine (COMPLEX*16)     │
│    - Translated from Python traction()           │
│    - Signature: (rdisp_z, [state_old_z], dt,     │
│                  props, stress_z, [state_new_z]) │
├─────────────────────────────────────────────────┤
│  Section 3: Flux subroutine (COMPLEX*16)         │
│    - Translated from Python flux() [if defined]  │
│    - Signature: (rdisp_z, temp_z, [state_old_z], │
│                  dt, props, flux_z)              │
├─────────────────────────────────────────────────┤
│  Section 4: CS tangent engine                    │
│    - Phase 1: perturb RDISP → DDSDDR, DDFDDR    │
│    - Phase 2: perturb TEMP → DDSDDT, DDFDDT     │
│    - All state_old reset per perturbation        │
├─────────────────────────────────────────────────┤
│  Section 5: Helper subroutines (if any)          │
│    - Translated from Python _helper() methods    │
├─────────────────────────────────────────────────┤
│  Section 6: tensor_ops.for (shared utilities)    │
└─────────────────────────────────────────────────┘
```

---

## 5. UINTER Wrapper Logic

### 5.1 Pseudocode

```fortran
SUBROUTINE UINTER(STRESS, DDSDDR, ..., STATEV, ..., RDISP, ...)

C     --- Read state variables ---
      IF (TIME(2) .EQ. 0.0d0) THEN
          d_old     = 0.0d0      ! initial value from state_vars dict
          H_old     = 0.0d0
      ELSE
          d_old     = STATEV(1)
          H_old     = STATEV(2)
      END IF

C     --- Compute Jacobian blocks via complex-step ---
      CALL prefix_cs_tangent(RDISP, TEMP,
     &     d_old, H_old, DTIME,
     &     PROPS, NDIR,
     &     DDSDDR, DDFDDT, DDSDDT, DDFDDR)

C     --- Compute real traction and state update ---
      CALL prefix_traction_real(RDISP,
     &     d_old, H_old, DTIME,
     &     PROPS,
     &     STRESS, d_new, H_new)

C     --- Compute real flux (if defined) ---
      CALL prefix_flux_real(RDISP, TEMP,
     &     d_old, DTIME, PROPS, FLUX)

C     --- Write state variables ---
      STATEV(1) = d_new
      STATEV(2) = H_new

C     --- Contact status flags ---
      IF (d_new .GT. 0.99d0 .AND. RDISP(1) .LT. 0.0d0) THEN
          IF (LOPENCLOSE .NE. 0) LSDI = 1   ! status changed
          LOPENCLOSE = 0                      ! open
      ELSE
          IF (LOPENCLOSE .NE. 1) LSDI = 1   ! status changed
          LOPENCLOSE = 1                      ! closed
      END IF

      RETURN
      END
```

### 5.2 Key Design Decisions

**Traction is computed twice**: once with complex arithmetic (inside CS
engine, for tangent) and once with real arithmetic (for actual STRESS
and state update). This follows the `umat_gen.py` pattern. Only the
real call updates state variables.

**Alternatively**, the real call can be the "unperturbed" CS call. In
`umat_gen.py`, the UMAT wrapper calls the stress subroutine directly
(real version) and separately calls the CS engine. For UINTER we can
do the same: the CS engine only produces Jacobians, and a separate
real-arithmetic wrapper computes STRESS and state updates.

**The real traction wrapper** calls the same translated Fortran logic
but with DOUBLE PRECISION instead of DOUBLE COMPLEX. This requires
generating two versions of the traction subroutine (complex and real),
or generating one complex version and a real wrapper that converts.
Following `umat_gen.py`'s pattern, we generate a single COMPLEX\*16
subroutine and call it with purely real inputs for the actual stress
evaluation, then extract the real part with `DBLE()`.

**LOPENCLOSE management**: the wrapper automatically sets the open/close
flag based on the damage state variable and the normal gap. The user
can override this by returning a `lopenclose` value from `traction()`.
If not returned, the default logic applies.

---

## 6. Python Class Interface: `SurfaceInteraction`

### 6.1 Base Class

```python
class SurfaceInteraction:
    """Base class for surface interaction definitions."""

    props = {}          # dict of property_name: default_value
    state_vars = {}     # dict of state_name: initial_value

    def traction(self, rdisp, ..., dt):
        """
        Compute interface traction.
        Must return: stress_vector, state_dict (if state_vars)
                 or: stress_vector (if no state_vars)
        """
        raise NotImplementedError

    def flux(self, rdisp, temp, ..., dt):
        """
        Compute heat flux (optional).
        Must return: [q_slave, q_master]
        """
        return None  # not defined by default
```

### 6.2 Method Signature Detection

The generator inspects `traction()` and `flux()` signatures to determine:

- Which arguments are "displacement-like" (`rdisp`) vs "temperature-like"
  (`temp`) — this determines which CS perturbation phase they belong to
- Which state variables are used (from the `_old` suffix convention)
- Whether `dt` is present (for rate-dependent models)
- Whether `flux()` is defined (for thermal coupling)

The recognized argument names:

| Argument | Type | Perturbation phase |
|----------|------|-------------------|
| `rdisp` | vector(NDIR) | Phase 1 (RDISP perturbation) |
| `temp` | vector(2) | Phase 2 (TEMP perturbation) |
| `*_old` | state variable | Not perturbed (real, reset per CS call) |
| `dt` | scalar | Not perturbed |

### 6.3 Dimension Handling (NDIR)

The user writes dimension-independent code using `rdisp[0]` for normal
and `rdisp[1:]` for tangential. The generator produces Fortran that
works for any NDIR (passed as argument). The traction subroutine
receives `rdisp(NDIR)` and returns `stress(NDIR)`.

For the CS engine, the perturbation loop runs from 1 to NDIR.

---

## 7. Comparison with umat_gen.py

| Aspect | umat_gen.py | uinter_gen.py |
|--------|-------------|---------------|
| Primary input | F (3x3 tensor) | RDISP (NDIR vector) |
| Primary output | PK1 stress (3x3) | Traction (NDIR vector) |
| Secondary input | — | TEMP (2-vector) |
| Secondary output | — | FLUX (2-vector) |
| Tangent | dPdF (3x3x3x3) | DDSDDR (NDIR×NDIR) + 3 cross blocks |
| Post-processing | PK1→Cauchy push-forward + Jaumann | None (already in local frame) |
| CS perturbations | 9 (3×3 components of F) | NDIR + 2 (RDISP + TEMP) |
| Spatial | Per integration point | Per contact point |
| Shape functions | Not involved (UMAT level) | Not involved (point-wise) |
| State vars | STATEV on element | STATEV on contact point |
| Special flags | — | LOPENCLOSE, LSTATE, LSDI |
| Time step control | PNEWDT | PNEWDT |

The UINTER generator is simpler than UMAT because:
- No push-forward / Jaumann correction needed
- Fewer perturbations (NDIR + 2 vs 9)
- No tensor index gymnastics (vectors, not tensors)
- No DDSDDE Voigt packing

---

## 8. vector_ops.for — Vector and Scalar Utilities

### 8.1 Rationale

`tensor_ops.for` handles 3×3 matrix operations (det33, inv33,
matmul33, etc.) and is shared by UMAT and UEL generators. UINTER
operates on variable-length vectors (NDIR = 2 or 3) and scalars,
not 3×3 tensors. A separate `vector_ops.for` keeps the two domains
clean:

- `tensor_ops.for` — fixed 3×3 tensors (UMAT, UEL)
- `vector_ops.for` — variable-length vectors, scalars (UINTER)

Both provide real and complex versions of each routine.

### 8.2 Required Subroutines and Functions

**Scalar operations (COMPLEX*16):**

| Routine | Signature | Description |
|---------|-----------|-------------|
| `smooth_ramp_z` | `FUNCTION smooth_ramp_z(x, eps)` | $S_\epsilon(x) = \frac{1}{2}(x + \sqrt{x^2 + \epsilon^2})$ |
| `smooth_abs_z` | `FUNCTION smooth_abs_z(x, eps)` | $\sqrt{x^2 + \epsilon^2}$ |
| `smooth_max_z` | `FUNCTION smooth_max_z(a, b, eps)` | $\frac{1}{2}(a + b + \sqrt{(a-b)^2 + \epsilon^2})$ |

**Vector operations (COMPLEX*16, variable length):**

| Routine | Signature | Description |
|---------|-----------|-------------|
| `dot_z` | `FUNCTION dot_z(a, b, n)` | $\sum_{i=1}^{n} a_i b_i$ |
| `norm_z` | `FUNCTION norm_z(a, n)` | $\sqrt{\sum_{i=1}^{n} a_i^2}$ |
| `norm_sq_z` | `FUNCTION norm_sq_z(a, n)` | $\sum_{i=1}^{n} a_i^2$ (avoids sqrt for energy) |

**Type conversion:**

| Routine | Signature | Description |
|---------|-----------|-------------|
| `real2complex_vec` | `SUBROUTINE(a_real, a_z, n)` | Convert real vector to complex with zero imaginary |
| `complex2real_vec` | `SUBROUTINE(a_z, a_real, n)` | Extract real part of complex vector |

**DOUBLE PRECISION versions** (same names without `_z` suffix) for
the real-arithmetic calls in the wrapper.

### 8.3 Design Notes

- All vector routines take `n` (length) as an argument — no fixed
  dimension. This handles NDIR=2 (2D) and NDIR=3 (3D) with the
  same code.
- `smooth_ramp_z` is the primary tool for CS-safe activation
  functions. It replaces all uses of `max(x, 0)` in constitutive
  laws where the argument is CS-perturbed.
- `smooth_max_z(a, b, eps)` handles `max(a, b)` for two
  CS-perturbed arguments (e.g., history variable update
  $\mathcal{H} = \max(\mathcal{H}_{\text{old}}, D)$).
- The template file is included in the generated `.for` file
  alongside `tensor_ops.for` (both are needed if the model calls
  tensor functions, e.g., for `det` in an advanced model).
  For most UINTER models, only `vector_ops.for` is needed.

### 8.4 File Location

```
abaqus_ufl/
  generators/
    templates/
      tensor_ops.for    ← existing (3×3 tensors)
      vector_ops.for    ← new (variable-length vectors + scalars)
```

---

## 9. AST Translator Extension for Vectors

### 9.1 Current State

`FortranTranslator` in `umat_gen.py` recognizes:
- **Tensors**: `F`, state variables with shape (3,3) — translated
  with `(i,j)` indexing
- **Scalars**: properties, scalar state variables, temporaries
- **Tensor functions**: `det(A)`, `inv(A)`, `trace(A)`, etc.

It does NOT recognize:
- **Vectors**: `rdisp`, `temp` — variable-length, indexed with `[i]`
- **Vector functions**: `dot()`, `norm()`
- **Smooth activation**: `smooth_ramp()`

### 9.2 Required Extensions

**New variable kind: `'vector'`**

The translator's `_var_kinds` dict gains a third kind alongside
`'scalar'` and `'tensor'`:

```python
_var_kinds = {
    'F': 'tensor',           # existing
    'rdisp': 'vector',       # new
    'temp': 'vector',        # new
    'd_old': 'scalar',       # existing pattern
}
```

**Subscript translation:**

| Python | Fortran | Rule |
|--------|---------|------|
| `rdisp[0]` | `rdisp_z(1)` | 0-indexed → 1-indexed |
| `rdisp[1]` | `rdisp_z(2)` | |
| `rdisp[2]` | `rdisp_z(3)` | Only valid when NDIR ≥ 3 |
| `temp[0]` | `temp_z(1)` | Slave surface |
| `temp[1]` | `temp_z(2)` | Master surface |

The translator converts `ast.Subscript` nodes for vector variables
by adding 1 to the index. For literal integer indices, this is
straightforward. For variable indices (`rdisp[i]`), it becomes
`rdisp_z(i+1)` — but this pattern is unlikely in typical UINTER
models where users access named components.

**NDIR-conditional indexing:**

If the user writes `rdisp[2]` (the second tangential component,
3D only), the generated code should guard:

```fortran
IF (ndir .GE. 3) THEN
    delta_t2 = rdisp_z(3)
END IF
```

Implementation: the translator tracks the maximum index used for
each vector variable. If `max_index >= 3` for `rdisp`, a warning
is emitted and a guard is generated. For Phase 1 (2D only), we
can require that `rdisp[2]` is not used and raise an error if it
appears.

**Vector return (list → stress_z array):**

The `return [stress_n, stress_t]` statement maps to:

```fortran
stress_z(1) = stress_n
stress_z(2) = stress_t
```

The translator detects `ast.Return` with an `ast.List` value and
generates component-wise assignment. For 3D:
`return [stress_n, stress_t1, stress_t2]` → three assignments.

**New function mappings:**

| Python | Fortran | Notes |
|--------|---------|-------|
| `dot(a, b)` | `dot_z(a, b, ndir)` | From `vector_ops.for` |
| `norm(a)` | `norm_z(a, ndir)` | From `vector_ops.for` |
| `norm_sq(a)` | `norm_sq_z(a, ndir)` | Avoids sqrt |
| `smooth_ramp(x)` | `smooth_ramp_z(x, eps_prop)` | $\epsilon$ from props |
| `smooth_max(a, b)` | `smooth_max_z(a, b, eps_prop)` | For history update |

The `smooth_ramp` function is special: if the user defines it as a
helper method (`self._smooth_ramp(x)`), the translator emits a call
to the generated helper subroutine (existing pattern). If it's a
recognized built-in name, the translator maps directly to the
`vector_ops.for` function. Both approaches work; the built-in path
is more convenient but the helper method path already works with
the existing translator.

**Warning for bare `max()` with CS-perturbed arguments:**

When the translator encounters `max(expr1, expr2)` where either
argument involves a CS-perturbed variable (`rdisp`, `temp`, or any
intermediate computed from them), it should emit a warning:

```
WARNING: max() is not differentiable at the transition point.
Consider using smooth_ramp() or smooth_max() for CS compatibility.
```

This is a translator-level diagnostic, not an error — bare `max()`
will still generate `MAX(DBLE(x), 0.0d0)` which works but gives
zero derivative at the kink.

### 9.3 Implementation Scope

Estimated additions to `FortranTranslator`:

| Feature | Lines | Difficulty |
|---------|-------|-----------|
| `'vector'` kind in `_var_kinds` | ~10 | Easy |
| Subscript `[i]` → `(i+1)` for vectors | ~20 | Easy |
| List return → component assignment | ~15 | Easy |
| `dot()`, `norm()`, `norm_sq()` mapping | ~20 | Easy |
| `smooth_ramp()`, `smooth_max()` mapping | ~15 | Easy |
| `max()` CS-safety warning | ~10 | Easy |
| NDIR-conditional guard for `rdisp[2]` | ~20 | Medium |
| **Total** | **~110** | |

This is a focused extension, not a rewrite. The existing translator
handles the hard parts (AST walking, expression translation, scope
management, helper method detection). The vector extension adds a
new kind and a handful of function mappings.

### 9.4 Prerequisite: Test Existing Translator First

The translator extension is validated only after `umat_gen.py` and
`uel_gen.py` are tested in Abaqus. Any bugs in the base translator
should be found and fixed on the simpler tensor case before extending
to vectors. The recommended sequence:

1. Test `umat_gen.py` with NeoHookean in Abaqus → fix translator
   bugs
2. Test `uel_gen.py` with a single element → fix UEL bugs
3. Extend translator for vectors → implement in `uinter_gen.py`
4. Test `uinter_gen.py` with linear cohesive law

---

## 10. User-Controlled Contact Flags and Time Step

### 10.1 Current Behavior

The UINTER wrapper currently sets LOPENCLOSE based on a hard-coded
rule: open if `RDISP(1) < -CHRLNGTH*1e-3`, closed otherwise. LSDI
is set to 1 whenever the status changes. PNEWDT is never modified.

This is insufficient for models where the contact status depends on
the constitutive state (e.g., keep contact "closed" during healing
even if the gap is slightly open).

### 10.2 Proposed Mechanism

The user's `traction()` method can return additional control values
in the state dict:

```python
return [stress_n, stress_t], {
    'd': d_new,
    'H': H_new,
    # Optional control overrides:
    'lopenclose': 1,       # force closed (e.g., during healing)
    'lsdi': 0,             # no SDI this iteration
    'pnewdt': 0.5,         # request time step cutback
}
```

### 10.3 Wrapper Logic

The generated UINTER wrapper checks for these keys after the real
traction call:

```fortran
C     --- Contact status ---
C     Check if user returned overrides
      IF (has_lopenclose_override) THEN
          IF (NINT(lopenclose_val) .NE. LOPENCLOSE) LSDI = 1
          LOPENCLOSE = NINT(lopenclose_val)
      ELSE
C         Default: based on normal gap
          IF (RDISP(1) .LT. -CHRLNGTH*1.0d-3) THEN
              IF (LOPENCLOSE .EQ. 1) LSDI = 1
              LOPENCLOSE = 0
          ELSE
              IF (LOPENCLOSE .EQ. 0) LSDI = 1
              LOPENCLOSE = 1
          END IF
      END IF

C     --- PNEWDT ---
      IF (has_pnewdt_override) THEN
          PNEWDT = pnewdt_val
      END IF
```

### 10.4 Implementation

The control overrides are **not** CS-perturbed — they are extracted
from the real (unperturbed) traction call only. They are stored as
regular DOUBLE PRECISION output arguments of the traction subroutine,
separate from the COMPLEX*16 stress and state outputs.

The generator detects these keys by checking whether the return dict
contains `'lopenclose'`, `'lsdi'`, or `'pnewdt'` during AST
analysis of the `traction()` method. If present, the wrapper
generates the override logic; if absent, the default behavior is
used.

---

## 11. DVISCOUS and DSTRUCTURAL

These arrays are part of the Abaqus UINTER signature and control
viscous and structural damping for dynamic analysis procedures
(direct steady-state, mode-based transient, complex eigenvalue
extraction, etc.).

For the quasi-static contact problems targeted by this generator,
these arrays should be **explicitly zeroed** in the wrapper:

```fortran
      DVISCOUS = 0.0d0
      DSTRUCTURAL = 0.0d0
```

If a future application requires dynamic damping at the interface,
the generator can be extended to accept a `damping()` method that
returns DVISCOUS and DSTRUCTURAL. This is not planned for the
initial implementation.

---

## 12. Shared Utilities: codegen_utils.py

### 12.1 Motivation

Currently, `uinter_gen.py` imports utility functions directly from
`umat_gen.py`. This creates a fragile dependency: `umat_gen.py`
owns functions that are logically shared across all generators.
The same issue affects `uel_gen.py`.

### 12.2 Proposed Module

Extract shared utilities into `codegen_utils.py`:

```
abaqus_ufl/
  generators/
    codegen_utils.py     ← shared utilities
    umat_gen.py          ← imports from codegen_utils
    uel_gen.py           ← imports from codegen_utils
    uinter_gen.py        ← imports from codegen_utils
    umatht_gen.py        ← future, imports from codegen_utils
    templates/
      tensor_ops.for
      vector_ops.for
```

**Functions to extract:**

| Function | Current location | Used by |
|----------|-----------------|---------|
| `_state_var_info()` | umat_gen.py | all |
| `_nstate_per_gp()` | umat_gen.py | all |
| `_fortran_dims()` | umat_gen.py | all |
| `_fortran_call()` | umat_gen.py | all |
| `_wrap_fortran_source()` | umat_gen.py | all |
| `FortranTranslator` | umat_gen.py | all |
| `Material` base class | currently in core | umat, uel |
| `SurfaceInteraction` base | uinter_gen.py | uinter |

### 12.3 Timing

This refactoring should happen **after** the existing generators
are tested and before the translator is extended for vectors.
The extraction is mechanical (move functions, update imports) and
reduces the risk of import breakage when adding new generators.

---

## 13. Translation Rules

The AST translator from `umat_gen.py` is reused with extensions
described in Section 9. Summary of all mappings:

### 13.1 Variable Access

| Python | Fortran | Variable kind |
|--------|---------|---------------|
| `rdisp[0]` | `rdisp_z(1)` | vector |
| `rdisp[1]` | `rdisp_z(2)` | vector |
| `rdisp[2]` | `rdisp_z(3)` | vector (3D only) |
| `temp[0]` | `temp_z(1)` | vector |
| `temp[1]` | `temp_z(2)` | vector |
| `self.K_n` | `props(1)` | property (by index) |
| `d_old` | `d_old_z` | scalar state (complex in CS) |
| `F[i,j]` | `F(i,j)` | tensor (existing) |

### 13.2 Function Mappings

| Python | Fortran | Source |
|--------|---------|--------|
| `dot(a, b)` | `dot_z(a, b, ndir)` | vector_ops.for |
| `norm(a)` | `norm_z(a, ndir)` | vector_ops.for |
| `norm_sq(a)` | `norm_sq_z(a, ndir)` | vector_ops.for |
| `smooth_ramp(x)` | `smooth_ramp_z(x, eps)` | vector_ops.for |
| `smooth_max(a, b)` | `smooth_max_z(a, b, eps)` | vector_ops.for |
| `sqrt(x)` | `CDSQRT(x)` | intrinsic (existing) |
| `exp(x)` | `CDEXP(x)` | intrinsic (existing) |
| `log(x)` | `CDLOG(x)` | intrinsic (existing) |
| `abs(x)` | `CDABS(x)` | intrinsic (existing) |
| `det(A)` | `det33z(A)` | tensor_ops.for (existing) |
| `inv(A)` | via `CALL inv33z` | tensor_ops.for (existing) |
| `max(a, b)` | `MAX(DBLE(a), DBLE(b))` | intrinsic (**CS warning**) |

### 13.3 Return Statement Mapping

| Python return | Fortran output | Context |
|---------------|----------------|---------|
| `return P` | `P_out(ii,jj) = P(ii,jj)` | UMAT (tensor) |
| `return [s_n, s_t]` | `stress_z(1) = s_n; stress_z(2) = s_t` | UINTER (vector) |
| `return P, {'ep': ep_new}` | `P_out = P; ep_new = ep_new` | UMAT + state |
| `return [s_n, s_t], {'d': d_new}` | `stress_z(1:2) + state` | UINTER + state |
| `return [q_s, q_m]` | `flux_z(1) = q_s; flux_z(2) = q_m` | UINTER flux |

---

## 14. Target Application: Healing Cohesive Law

The primary application motivating `uinter_gen.py` is an
Anand-framework-inspired cohesive law with independent strength control
and healing. The Python definition would be:

```python
class HealingCohesive(au.SurfaceInteraction):
    props = dict(
        K_n=1e6,           # normal stiffness
        K_t=1e5,           # tangential stiffness
        sigma_n=1e3,       # normal adhesion strength (independent)
        sigma_t=5e2,       # tangential adhesion strength (independent)
        zeta_n=1e4,        # normal damage growth rate
        zeta_t=1e4,        # tangential damage growth rate
        psi_star=1e3,      # fracture resistance energy density
        eta=1.0,           # viscous regularization
        k_h=0.1,           # healing rate
        delta_heal=0.01,   # healing activation gap
        K_pen=1e7,         # contact penalty for compression
        eps=1e-4,          # smooth ramp parameter
    )

    state_vars = dict(
        d=0.0,             # damage variable
        H=0.0,             # history variable (max driving force)
        t_contact=0.0,     # accumulated contact time
    )

    def _smooth_ramp(self, x):
        """Smooth approximation of max(x, 0)."""
        return 0.5 * (x + sqrt(x**2 + self.eps**2))

    def traction(self, rdisp, d_old, H_old, t_contact_old, dt):
        delta_n = rdisp[0]       # positive = penetration
        delta_t_sq = rdisp[1]**2  # 2D; extend for 3D

        # Undamaged cohesive energies (tension only for normal)
        psi_n = 0.5 * self.K_n * self._smooth_ramp(-delta_n)**2
        psi_t = 0.5 * self.K_t * delta_t_sq

        # Driving forces with independent thresholds
        psi_cr_n = self.sigma_n**2 / (2.0 * self.K_n)
        psi_cr_t = self.sigma_t**2 / (2.0 * self.K_t)
        D_n = self.zeta_n * self._smooth_ramp(psi_n / psi_cr_n - 1.0)
        D_t = self.zeta_t * self._smooth_ramp(psi_t / psi_cr_t - 1.0)
        D = D_n + D_t

        # History update (max over time)
        H_new = max(H_old, D)  # Note: max of two scalars, OK for CS

        # Healing activation
        h = self._smooth_ramp(-delta_n - self.delta_heal)

        # Damage evolution (backward Euler, implicit in d)
        # eta * (d - d_old) / dt = 2(1-d)*H - 2*psi_star*d - k_h*d*h
        # Rearrange: d * (eta/dt + 2*H + 2*psi_star + k_h*h) 
        #          = d_old * eta/dt + 2*H
        numer = d_old * self.eta / dt + 2.0 * H_new
        denom = self.eta / dt + 2.0 * H_new + 2.0 * self.psi_star + self.k_h * h
        d_new = numer / denom
        # Note: bounds [0,1] handled by the structure of numer/denom
        #       when H >= 0, psi_star > 0, k_h >= 0, h >= 0

        # Degraded traction
        g = (1.0 - d_new)**2

        # Adhesive traction (tension only for normal)
        T_n_adhesion = -g * self.K_n * self._smooth_ramp(-delta_n)
        T_t_adhesion = g * self.K_t * rdisp[1]

        # Contact repulsion (compression, independent of damage)
        T_n_contact = self.K_pen * self._smooth_ramp(delta_n)

        stress_n = T_n_contact + T_n_adhesion
        stress_t = T_t_adhesion

        # Contact time tracking
        t_contact_new = t_contact_old + dt if delta_n > -self.delta_heal else 0.0

        return [stress_n, stress_t], {
            'd': d_new, 'H': H_new, 't_contact': t_contact_new
        }
```

Key features of this model:

1. **Independent strength and energy**: $\sigma_n$ and $\sigma_t$ control
   initiation; $\psi^\star$ controls fracture resistance. Inspired by
   Anand (2026).

2. **Implicit damage update**: The backward Euler solve for $d$ is
   algebraic (linear in $d$ after rearrangement). This ensures stability
   and avoids negative $d$ when the denominator is positive.

3. **All smooth**: `_smooth_ramp()` replaces `max()` everywhere, ensuring
   $C^\infty$ response and quadratic Newton convergence via CS tangent.

4. **Healing**: activated by the $h(\delta_n)$ term when surfaces are
   in proximity. Rate controlled by $k_h$.

5. **History variable**: $\mathcal{H} = \max(H_{\text{old}}, D)$ for
   irreversible damage direction; healing operates independently.

---

## 15. Implementation Plan

### Prerequisites (must complete before uinter_gen work)

1. **Test `umat_gen.py` with Abaqus** — run NeoHookean UMAT, verify
   stress and DDSDDE against analytical. Fix translator bugs.
2. **Test `uel_gen.py`** — single Quad8 element, verify
   residual and tangent. Fix UEL bugs.
3. **Extract `codegen_utils.py`** — move shared utilities, update
   imports in all generators. Verify all existing tests still pass.

### Phase 1: vector_ops.for + AST extension

- Write `vector_ops.for` template (Section 8)
- Extend `FortranTranslator` for `'vector'` kind (Section 9)
- Add `smooth_ramp`, `dot`, `norm` function mappings
- Add list-return translation for vector outputs
- Add `max()` CS-safety warning

### Phase 2: Complete traction translation

- Implement `_generate_traction_subroutine()` with full AST
  translation
- Handle `_smooth_ramp` and other helper methods (existing pattern)
- Test with stateless linear cohesive law:
    - `T_n = K_n * rdisp[0]`, `T_t = K_t * rdisp[1]`
    - Verify DDSDDR = `[[K_n, 0], [0, K_t]]` from CS matches
      analytical

### Phase 3: State variables + healing model

- Test with stateful damage model (implicit damage update)
- Verify CS tangent against finite differences
- Implement user-controlled flags (Section 10)
- Implement DVISCOUS/DSTRUCTURAL zeroing (Section 11)

### Phase 4: Flux (thermal coupling)

- Implement `_generate_flux_subroutine()` with AST translation
- CS engine Phase 2: TEMP perturbation
- Cross-coupling tangent blocks (DDSDDT, DDFDDR)
- Test with gap conductance model

### Phase 5: 3D support

- Handle NDIR=3 (two tangential components)
- NDIR-conditional guards for `rdisp[2]`
- Test with 3D contact patch

### Phase 6: Validation against Abaqus

- Run healing cohesive model (Section 14) in Abaqus
- Compare Newton iteration counts with built-in cohesive behavior
- Verify re-adhesion behavior that built-in cannot reproduce

---

## 16. Files

| File | Description |
|------|-------------|
| `uinter_gen.py` | Main generator module (`generate_uinter`, `SurfaceInteraction`) |
| `generators/templates/tensor_ops.for` | Shared 3×3 math utilities |
| `vector_ops.for` | Vector/scalar template (proposed; see Section 8) |
| `codegen_utils.py` | Shared utilities extracted from `umat_gen.py` (proposed; see Section 12) |

---

## 17. References

- Abaqus 2017 Subroutines Reference Manual, UINTER section
- Anand L (2026). Fracture of rock-like materials: A gradient-damage
  theory. *Int J Solids Struct* 325:113739.
- Chester SA. `umat_simpleFeFp.for` — reference implementation for
  STATEV handling patterns.
- Wu JY (2017). A unified phase-field theory for the mechanics of
  damage and quasi-brittle failure. *J Mech Phys Solids* 103:72–99.
- Wriggers P, Schröder J, Schwarz A (2013). A finite element method
  for contact using a third medium. *Comput Mech* 52:837–847.
