# Matrix Functions Design

This document describes two self-contained, complex-step-compatible ways to
evaluate matrix functions (`sqrtm`, `logm`, `expm`) on 3×3 tensors, and when to
use each. The iterative backend is the default in both the Python oracle and the
generated Fortran (`tensor.py`, `tensor_ops.for`); the eigendecomposition
backend is available as an explicit generator option.

## Motivation

Constitutive models based on logarithmic strain (Hencky), exponential maps
(plasticity return mapping), and matrix square roots need matrix functions
like `logm(A)`, `expm(A)`, `sqrtm(A)`.

Two approaches are available, both self-contained (no LAPACK, no external
dependencies) and CS-compatible.

## Approach A: Eigendecomposition

For a 3×3 symmetric matrix A with eigenvalues λ_i and eigenvectors V:

```
sqrtm(A) = V @ diag(sqrt(λ_i)) @ V^{-1}
logm(A)  = V @ diag(log(λ_i))  @ V^{-1}
expm(A)  = V @ diag(exp(λ_i))  @ V^{-1}
```

### Why this works for complex-step but NOT for symbolic AD

In symbolic AD (FEniCSx, JAX), the derivative of eigendecomposition
involves `1/(λ_i − λ_j)`, which is singular when eigenvalues are repeated
(common at F = I and near isotropic states). This is why FEniCSx does not
provide eigenvalue-based `logm`.

In our complex-step framework, this problem does not exist. The tangent is:
```
d(logm(A))/dA_kl = AIMAG(logm(A + i·h·e_kl)) / h
```
The eigenvalues of `A + i·h·e_kl` (a complex non-Hermitian matrix) are
always distinct — the complex perturbation breaks degeneracy. No derivative
of eigendecomposition is ever computed explicitly. The CS derivative flows
through `eig → log → reconstruct` naturally.

### Components (all self-contained, no external dependencies)

**`eig33z` — eigenvalues via Cardano's cubic formula (~50 lines)**

The characteristic polynomial of a 3×3 matrix is a cubic. Cardano's formula
gives the three roots analytically using only `+`, `-`, `*`, `/`, `sqrt`,
`cbrt`, and trig functions — all available as Fortran COMPLEX*16 intrinsics.
No LAPACK, no iteration.

```fortran
      SUBROUTINE eig33z(A, lam, V)
C     Eigenvalues and eigenvectors of 3x3 symmetric matrix
C     Input:  A(3,3) — COMPLEX*16
C     Output: lam(3) — eigenvalues, V(3,3) — eigenvector columns
C     Uses Cardano's formula for the cubic characteristic polynomial
      ...
      END SUBROUTINE
```

**`sqrtm33z`, `logm33z`, `expm33z` — ~10 lines each on top of `eig33z`**

For a symmetric matrix with orthonormal eigenvectors, reconstruction can use
`V^T`, as the illustrative code below does. The generated implementation instead
returns unnormalized cross-product eigenvectors and uses `inv(V)` for
reconstruction, which does not require normalizing the eigenvectors.

```fortran
      SUBROUTINE logm33z(A, L)
C     Matrix logarithm of symmetric A:
C       L = V * diag(log(lam)) * V^T
C     (V^T used because eigenvectors are orthonormal for symmetric A)
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: L(3,3)
      DOUBLE COMPLEX :: lam(3), V(3,3), VT(3,3)
      DOUBLE COMPLEX :: D(3,3), T(3,3)

      CALL eig33z(A, lam, V)
      CALL transpose33z(V, VT)

C     D = diag(log(lam))
      D = DCMPLX(0.0d0, 0.0d0)
      D(1,1) = LOG(lam(1))
      D(2,2) = LOG(lam(2))
      D(3,3) = LOG(lam(3))

C     L = V * D * V^T
      CALL matmul33z(V, D, T)
      CALL matmul33z(T, VT, L)

      END SUBROUTINE
```

Same pattern for `sqrtm33z` (replace `LOG` with `SQRT`) and `expm33z`
(replace `LOG` with `EXP`).

**Note on `eig33z` implementation:** The eigenvector computation must
produce orthonormal vectors. For a symmetric matrix with distinct
eigenvalues (guaranteed by CS perturbation), eigenvectors of different
eigenvalues are automatically orthogonal. Each vector must be normalized
to unit length. The computation uses cross products of rows of `(A - λ_i I)`
followed by normalization — approximately 20 lines of Fortran.

### Advantages
- Simple: ~10 lines per matrix function on top of `eig33z`
- Fast: one cubic solve + two matrix multiplies
- General: works for any scalar function f(λ) applied to a matrix
- No iteration, no convergence parameters
- Zero external dependencies

### Limitations and cautions

**Branch cut trap in Cardano's formula:** When the input is complex
(`A + i·h·E_kl`), the cubic formula involves complex cube roots and
(in the trigonometric form) complex arccosine, both of which have branch
cuts. If the perturbation pushes a value across a branch cut, eigenvalues
can swap or jump, destroying the CS derivative.

Mitigation:
- Use the algebraic form of Cardano (not trigonometric)
- Sort eigenvalues by real part to maintain consistent ordering
- Test rigorously at A = I (triply degenerate) and near-degenerate states
- Verify CS derivatives against FD at these states

**Eigenvector computation:** The design shows eigenvalues but the
reconstruction `V @ diag(f(λ)) @ V^{-1}` also needs eigenvectors V.
For symmetric matrices with distinct eigenvalues (guaranteed by CS
perturbation), eigenvectors can be computed via cross products of rows
of `(A - λ_i I)` — approximately 20 extra lines of Fortran, no
external dependencies.

**Intended for symmetric matrices:** Cardano gives eigenvalues for any
3×3 matrix, but eigenvector reconstruction via `V^{-1}` can be unstable
if V is nearly singular (ill-conditioned eigenvectors for non-symmetric
matrices). In continuum mechanics, the tensors we apply matrix functions
to (C, B, U, σ, S) are always symmetric, so this is not a practical
limitation. For non-symmetric cases, use Approach B.

---

## Approach B: Iterative (alternative)

Uses only `matmul` and `inv` — no eigendecomposition at all. Originally
developed for symbolic AD frameworks where eigendecomposition derivatives
are problematic. Useful as a fallback or for non-symmetric matrices.

### Matrix Square Root: Denman-Beavers Iteration

```
Y_0 = A,  Z_0 = I

Y_{k+1} = (Y_k + Z_k^{-1}) / 2
Z_{k+1} = (Z_k + Y_k^{-1}) / 2

Converges: Y_k → A^{1/2}, Z_k → A^{-1/2}
```

Quadratic convergence. For 3×3 SPD matrices, 15-20 iterations give
machine precision.

```fortran
      SUBROUTINE sqrtm33z_iter(A, S, Sinv)
C     Matrix square root via Denman-Beavers iteration
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: S(3,3), Sinv(3,3)
      DOUBLE COMPLEX :: Y(3,3), Z(3,3), Yinv(3,3), Zinv(3,3)
      INTEGER :: iter
      INTEGER, PARAMETER :: MAXITER = 20

      Y = A
      CALL eye33z(Z)
      DO iter = 1, MAXITER
        CALL inv33z(Z, Zinv)
        CALL inv33z(Y, Yinv)
        Y = (Y + Zinv) * DCMPLX(0.5d0, 0.0d0)
        Z = (Z + Yinv) * DCMPLX(0.5d0, 0.0d0)
      END DO
      S = Y
      Sinv = Z
      END SUBROUTINE
```

### Matrix Logarithm: Inverse Scaling & Squaring

```
Step 1: B = A
        for k = 1 to NSCALE:        (fixed count, e.g. 10)
            B = sqrtm(B)
        B is now A^{1/1024} ≈ I

Step 2: X = (B - I)(B + I)^{-1}     (Cayley transform)
        log(B) = 2 * [X + X³/3 + X⁵/5 + ...]

Step 3: log(A) = 2^NSCALE * log(B)
```

Using a fixed iteration count (NSCALE=10) eliminates all branching on
matrix values, making the algorithm perfectly CS-compatible.

```fortran
      SUBROUTINE logm33z_iter(A, L)
C     Matrix logarithm via inverse scaling & squaring
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: L(3,3)
      DOUBLE COMPLEX :: B(3,3), Bsqrt(3,3), Binv(3,3)
      DOUBLE COMPLEX :: X(3,3), Xpow(3,3), T(3,3)
      DOUBLE COMPLEX :: BpI(3,3), BpI_inv(3,3), II(3,3)
      INTEGER :: k
      INTEGER, PARAMETER :: NSCALE = 10
      INTEGER, PARAMETER :: NSERIES = 8

C     Step 1: Fixed scaling — B = A^{1/2^NSCALE}
      B = A
      DO k = 1, NSCALE
        CALL sqrtm33z_iter(B, Bsqrt, Binv)
        B = Bsqrt
      END DO

C     Step 2: Gregory series for log(B) where B ≈ I
      CALL eye33z(II)
      BpI = B + II
      CALL inv33z(BpI, BpI_inv)
      CALL matmul33z(B - II, BpI_inv, X)

      L = X
      CALL matmul33z(X, X, T)
      Xpow = X
      DO k = 3, 2*NSERIES+1, 2
        CALL matmul33z(Xpow, T, Xpow)
        L = L + Xpow / DCMPLX(DBLE(k), 0.0d0)
      END DO
      L = L * DCMPLX(2.0d0, 0.0d0)

C     Step 3: Undo scaling
      L = L * DCMPLX(DBLE(2**NSCALE), 0.0d0)

      END SUBROUTINE
```

### Matrix Exponential: Scaling & Squaring + Padé

Same idea in reverse: scale A down by a fixed factor, apply Padé
approximant for exp, then square up. Fixed scaling count, no branching.

### Advantages
- No eigendecomposition needed — avoids branch cut issues entirely
- Works for non-symmetric matrices
- Denman-Beavers gives both `sqrtm` and `sqrtm_inv` simultaneously
- Historically the most robust choice for plasticity (exponential map)

### Limitations and cautions

**SPD assumption:** Denman-Beavers converges to the principal square root
only for symmetric positive definite (SPD) matrices. For non-SPD inputs
(which should not arise in normal mechanics — C and B are always SPD),
the iteration may converge to a non-principal branch.

**Fixed vs dynamic iteration count:** Two options, both CS-safe:

Option 1 (fixed): Always do NSCALE=10 iterations. Simple, predictable,
but wasteful when A ≈ I.

Option 2 (dynamic, recommended): Branch on the REAL part of the error.
Because the imaginary perturbation (h = 1e-10) is infinitesimal, it does
not change the convergence trajectory. The complex part "rides along" for
exactly as many iterations as the real physics dictates:

```fortran
      Y_old = Y
      DO iter = 1, 100
C       ... Denman-Beavers update ...
C       Check convergence on REAL part only (CS-safe!)
        err = 0.0d0
        DO i = 1, 3
          DO j = 1, 3
            err = MAX(err, ABS(DBLE(Y(i,j))-DBLE(Y_old(i,j))))
          END DO
        END DO
        IF (err .LT. 1.0d-14) EXIT
        Y_old = Y
      END DO
```

This is perfectly CS-compatible because the branch condition depends only
on the real part — the same code path is taken for perturbed and
unperturbed evaluations.

### Algorithm selection by use case

| Use case | Recommended approach | Why |
|----------|---------------------|-----|
| Hencky strain (logm of C) | A (eigendecomposition) | C is SPD, need principal stretches |
| Ogden model | A (eigendecomposition) | Need individual λ_i explicitly |
| J2 plasticity (expm) | B (iterative) | Plastic flow direction often degenerate |
| Viscoelasticity (expm) | B (iterative) | More robust for general matrices |
| Damage (principal stresses) | A (eigendecomposition) | Need individual eigenvalues |

---

## Comparison

| | Approach A (eigendecomposition) | Approach B (iterative) |
|---|---|---|
| Dependencies | `eig33z` (~70 lines) | None beyond matmul/inv |
| Lines per function | ~10 | ~30-40 |
| Speed | Fast (cubic + 2 matmuls) | Moderate (adaptive iteration) |
| CS-compatible | ✅ (with branch cut care) | ✅ (branch on DBLE only) |
| Non-symmetric matrices | ❌ (symmetric only) | ✅ (SPD for sqrtm) |
| Symbolic AD compatible | ❌ (1/(λ_i−λ_j) singularity) | ✅ |
| External dependencies | None | None |
| Degenerate matrices | Requires testing at A≈I | Robust (no eigenvalues) |

## Recommendation

Implement both. Each has a natural use case:

- **Approach A** — high-performance default for hyperelasticity (Hencky,
  Ogden, damage). The tensors involved (C, B, σ) are always symmetric
  and well-conditioned.
- **Approach B** — robust fallback for plasticity (exponential map) and
  any case where the matrix may be degenerate or non-symmetric.

In `tensor_ops.for`:
- `eig33z`, `sqrtm33z`, `logm33z`, `expm33z` — eigendecomposition-based
- `sqrtm33z_iter`, `logm33z_iter`, `expm33z_iter` — iterative

Current `tensor.py` uses the iterative matrix functions directly. Current UMAT
generation defaults to the iterative Fortran backend; the eig-based backend is
available only by explicit generator option.

## What These Enable

| Model | Needs | Recommended approach |
|-------|-------|---------------------|
| Hencky strain | ε = log(U) = ½ log(C) | A (`logm33z`) |
| Ogden | ψ(λ_1, λ_2, λ_3) | A (`eig33z` directly) |
| J2 plasticity (exp map) | F^p_{n+1} = exp(Δε^p) F^p_n | B (`expm33z_iter`) |
| Viscoelasticity | exp(-Δt/τ · A) | B (`expm33z_iter`) |
| Damage (principal stress) | max(σ_1, σ_2, σ_3) | A (`eig33z` directly) |
| Phase-field fracture | tension-compression split | A (`eig33z` directly) |

## CS Compatibility

Both approaches use only holomorphic operations (matmul, inv/transpose,
scalar functions). The complex-step derivative flows through naturally.

Approach A: CS breaks eigenvalue degeneracy, so `eig33z` always returns
distinct eigenvalues for perturbed input. The reconstruction
`V @ diag(f(λ)) @ Vᵀ` is smooth. **Caution:** branch cuts in Cardano's
formula require sorting eigenvalues by real part and testing at degenerate
states (see Approach A limitations).

Approach B: Dynamic convergence on `DBLE(err)` ensures the same iteration
count for real and perturbed evaluations. Alternatively, fixed iteration
counts guarantee identical code paths unconditionally.

## Fallback Logic

If the user requests `method='eig'` but the input matrix is detected as
non-symmetric (e.g., `max|A - Aᵀ| > tol`), automatically switch to
`method='iter'` with a warning. This prevents silent numerical errors
from ill-conditioned eigenvectors.

## Matrix Exponential: Iterative (Scaling & Squaring + Padé)

For completeness, here is the iterative expm sketch:

```fortran
      SUBROUTINE expm33z_iter(A, E)
C     Matrix exponential via scaling & squaring + [2/2] Padé
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: E(3,3)
      DOUBLE COMPLEX :: B(3,3), II(3,3)
      DOUBLE COMPLEX :: N2(3,3), D2(3,3), T1(3,3), T2(3,3)
      DOUBLE COMPLEX :: D2inv(3,3)
      INTEGER :: s, k
      INTEGER, PARAMETER :: NSCALE = 10

      CALL eye33z(II)

C     Step 1: Scale down — B = A / 2^NSCALE
      B = A / DCMPLX(DBLE(2**NSCALE), 0.0d0)

C     Step 2: [2/2] Padé approximant
C     exp(B) ≈ (I - B/2 + B²/12)^{-1} (I + B/2 + B²/12)
      CALL matmul33z(B, B, T1)       ! T1 = B²
      N2 = II + B*DCMPLX(0.5d0,0d0)
     &     + T1*DCMPLX(1d0/12d0,0d0) ! numerator
      D2 = II - B*DCMPLX(0.5d0,0d0)
     &     + T1*DCMPLX(1d0/12d0,0d0) ! denominator
      CALL inv33z(D2, D2inv)
      CALL matmul33z(D2inv, N2, E)   ! E = D2^{-1} * N2

C     Step 3: Square up — E = E^{2^NSCALE}
      DO k = 1, NSCALE
        CALL matmul33z(E, E, T1)
        E = T1
      END DO

      END SUBROUTINE
```

The [2/2] Padé gives ~6 digits of accuracy for small B. With NSCALE=10,
B = A/1024 is very small, so the approximation is excellent. For machine
precision, use [6/6] Padé (more terms but same structure).

## Testing Plan

Essential tests before using matrix functions in production:

1. **CS vs FD at F = I** (triply degenerate eigenvalues):
   Compute `d(logm(C))/dC_kl` via CS and FD. Must agree to ~1e-6.
   This is the hardest test for branch cut handling.

2. **CS vs FD near F = I** (nearly degenerate):
   F = I + 0.001 * random. Tests that eigenvalue sorting is consistent
   under tiny perturbations.

3. **CS vs FD at large deformation** (well-separated eigenvalues):
   F = diag(2.0, 0.7, 1.0). Easy case — should pass trivially.

4. **Approach A vs Approach B consistency:**
   For the same SPD input, `logm33z(A)` and `logm33z_iter(A)` must
   agree to machine precision. Any discrepancy indicates a bug in one
   of the implementations.

5. **Known analytical results:**
   `logm(diag(a,b,c))` = `diag(log(a), log(b), log(c))`.
   `expm(diag(a,b,c))` = `diag(exp(a), exp(b), exp(c))`.
   `sqrtm(I)` = `I`.

## Implementation Priority

1. `eig33z` — Cardano's formula + orthonormal eigenvectors
2. `logm33z` — eigendecomposition-based (Hencky strain)
3. `expm33z` — eigendecomposition-based (plasticity)
4. `sqrtm33z` — eigendecomposition-based (stretch tensor U)
5. `sqrtm33z_iter` — Denman-Beavers (fallback)
6. `logm33z_iter` — inverse scaling & squaring (fallback)
7. `expm33z_iter` — scaling & squaring + Padé (fallback)
