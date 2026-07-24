# Writing Complex-Step Safe Constitutive Models

The complex-step method computes exact first derivatives by evaluating
your material function with a small imaginary perturbation:

    df/dx = Im[ f(x + ih) ] / h,    h = 1e-10

This works for any function that uses only **holomorphic** (complex-analytic)
operations. Most standard math operations are holomorphic. A few common
patterns are not, and will silently produce wrong tangents.

## Rules

### Rule 1: No `abs`, `max`, `min` on complex arguments

These discard or scramble the imaginary part.

```python
# WRONG
sigma_eq = sqrt(abs(s11**2 + s22**2 - s11*s22 + 3*s12**2))

# CORRECT — abs is not needed if the argument is always positive
sigma_eq = sqrt(s11**2 + s22**2 - s11*s22 + 3*s12**2)
```

### Rule 2: Branch on the real part only

`if` statements must compare real values, not complex values.

```python
# WRONG — comparing complex numbers is undefined
if f_trial > 0:
    # plastic correction

# CORRECT — branch on real part, arithmetic stays complex
if f_trial.real > 0:     # Python
# IF (DBLE(f_trial) .GT. 0.0d0) THEN   ! Fortran
    # plastic correction (all math here uses complex variables)
```

The imaginary perturbation is ~1e-10, so `x.real > 0` gives the same
branch as the unperturbed evaluation. The complex arithmetic inside
the branch propagates the derivative information correctly.

### Rule 3: State variables — previous step is real, current step is complex

State variables from the **previous converged step** (read from SVARS)
are always real. The **current step** update must be computed in complex
arithmetic to capture the consistent (algorithmic) tangent.

State arrives as `<name>_old` arguments (always real), and the method returns
`(stress, {name: new_value})`; the generator writes the dict back to `STATEV`.

```python
class J2Plasticity(au.Material):
    state_vars = dict(ep=0.0)   # equivalent plastic strain

    def stress_PK1(self, F, ep_old, dt):
        # ep_old is REAL (previous converged step); current-step quantities are
        # COMPLEX during tangent computation.

        # Trial elastic state (complex)
        Fe_trial = F   # simplified; real model has multiplicative split
        sigma_trial = self.elastic_stress(Fe_trial)
        s_trial = dev(sigma_trial)
        f_trial = von_mises(s_trial) - (self.sigma_y + self.H * ep_old)

        # Plastic correction — branch on the real part
        if f_trial.real > 0:
            # Return mapping (Newton on real part, arithmetic on complex)
            dgamma = f_trial / (3*self.G + self.H)   # COMPLEX
            # ... radial return ...
            ep_new = ep_old + dgamma   # COMPLEX; captures d(ep)/d(strain)
        else:
            ep_new = ep_old            # elastic step

        # Return (stress, updated-state dict); ep_new is COMPLEX during the
        # tangent pass and real during the residual pass.
        return P, {"ep": ep_new}
```

**Critical:** Do NOT cast the returned state back to real during tangent
computation. Writing `ep_new = float(ep_new)` (or `DBLE(...)` in Fortran) loses
the algorithmic tangent and falls back to the continuum tangent, destroying
quadratic Newton convergence.

### Rule 4: Iterative algorithms — converge on real, compute on complex

For Newton iterations inside the material (e.g., return mapping):

```python
# Return mapping Newton loop
dgamma = complex(0.0, 0.0)     # complex accumulator
for iteration in range(20):
    # Compute residual (complex)
    f = von_mises(s_trial - 2*G*dgamma*n_trial) - (sigma_y + H*(ep_old + dgamma))

    # Check convergence on REAL PART only
    if abs(f.real) < tol:
        break

    # Newton update (complex — preserves derivative info)
    df_ddgamma = -2*G - H
    dgamma = dgamma - f / df_ddgamma
```

### Rule 5: Eigenvalues work naturally

When you compute eigenvalues of a real symmetric matrix that has been
perturbed by `ih` in one component, the perturbed matrix is no longer
symmetric (it's complex non-Hermitian). Its eigenvalues are distinct
complex numbers even when the real eigenvalues are repeated.

This means `eig(C)` works without any special handling:

```python
def stress_PK1(self, F):
    C = F.T @ F               # complex 3x3
    lam, N = eig(C)           # complex eigenvalues, always distinct
    lam_stretch = sqrt(lam)   # complex principal stretches
    # ... Ogden energy, Hencky strain, etc.
```

No branch-free algorithms are needed to obtain a well-defined *tangent*: the
imaginary complex-step perturbation splits otherwise-repeated eigenvalues, so
`eig(C)` differentiates cleanly. This does **not** make the `eig`-based
matrix-function backend robust at near-degenerate states — its near-diagonal
guards can zero shear-block derivatives — which is why the generator's default
matrix backend is `iterative`, not `eig`.

### Rule 6: Allowed operations

| Operation | Complex-safe? | Notes |
|-----------|:---:|-------|
| `+`, `-`, `*`, `/` | Yes | |
| `**n` (integer power) | Yes | |
| `exp`, `log`, `sqrt` | Yes | Use complex intrinsics (CDEXP, CDLOG, CDSQRT) |
| `sin`, `cos`, `tan` | Yes | Rarely needed in constitutive models |
| `det(A)` | Yes | Explicit formula, no branching |
| `inv(A)` | Yes | Explicit cofactor formula |
| `A @ B` (matmul) | Yes | |
| `A.T` (transpose) | Yes | |
| `eig(A)` / `eigh(A)` | Yes | Self-contained CS-safe eigensolver, no LAPACK; `eigh` is an alias for the general `eig` |
| `trace(A)` | Yes | |
| `0.5*(A + A.T)`, `dev(A)` | Yes | Symmetric/deviatoric parts; the translator has no unary `sym` |
| `abs(x)` | **No** | Use `x` directly or `sqrt(x*x)` if needed |
| `max(a, b)` | **No** | Use `if a.real > b.real` branching |
| `min(a, b)` | **No** | Same |
| `sign(x)` | **No** | Branch on real part instead |
| `np.clip` | **No** | Branch on real part |

### Rule 7: SVARS layout in generated Fortran

State variables are stored in the flat SVARS array, carved up per Gauss point:

```fortran
C     Read state at Gauss point KINTK
      LOC = (KINTK - 1) * nStatePer
      ep_old = SVARS(LOC + 1)                 ! scalar state (1 slot)
      Fp_old = SVARS(LOC + 2 : LOC + 10)      ! (3,3) tensor state (9 slots)

C     ... compute new state (complex during tangent) ...

C     Write state back (real part only — complex part is derivative info,
C     not stored between increments)
      SVARS(LOC + 1)            = DBLE(ep_new)
      SVARS(LOC + 2 : LOC + 10) = DBLE(Fp_new)
```

When computing the tangent, `ep_new` and `Fp_new` are complex. Their real
parts are the physical state update; their imaginary parts carry the derivative
information used by the tangent engine. Only the real parts are written to SVARS.

Supported state-variable shapes are scalar (1 slot), vector `(3,)` (3 slots),
and tensor `(3,3)` (9 slots).
