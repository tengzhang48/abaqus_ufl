# Matrix Functions Design: `abaqus_ufl/core/tensor.py`

Three-by-three matrix functions for constitutive models.  Every
implementation is complex-step (CS) safe: fixed iteration counts, no
branching on imaginary parts, no `abs`/`max`/`sign` on complex
variables.

---

## Current Status

There are two active backends:

| Context | `eig(A)` | `sqrtm/logm/expm/polar` | Status |
|---------|----------|-------------------------|--------|
| Python material oracle (`abaqus_ufl/core/tensor.py`) | Trigonometric cubic + guards | Iterative fixed-count algorithms | Active |
| Generated Fortran default (`matrix_backend="iterative"`) | Same `eig33z` available when user calls `eig` | Iterative helpers `sqrtm33z_iter/logm33z_iter/expm33z_iter/polar33z_iter` | **Active generator default** |
| Generated Fortran alternate (`matrix_backend="eig"`) | Trigonometric cubic + guards | Eigendecomposition helpers `sqrtm33z/logm33z/expm33z/polar33z` | Backward compatibility |

This split is important.  The Python path is the reference oracle used by
`verify()` and material-point tests.  The generated Fortran **default is now
iterative** because the eig-based backend is not CS-safe at diagonal-F states:
`eig33z`'s absolute near-diagonal guard (`ABS(p1) < 1e-18`) triggers for the
tiny imaginary off-diagonals produced by complex-step perturbation, returning
V=I and zeroing the shear block of the tangent.  The iterative backend avoids
eigendecomposition entirely inside matrix functions and is therefore CS-safe.

`eig` (and its alias `eigh`) map to the same general eigensolver primitive;
`eigh` is not a separate symmetric solver. The `eig` matrix-function backend
remains available for backward compatibility:

```python
au.generate_umat(model, "model_eig.for", matrix_backend="eig")
```

## Algorithm table

| Function | Algorithm | Reference | Iterations / order |
|----------|-----------|-----------|-------------------:|
| `eig(A)` | Trigonometric depressed-cubic + diagonal guard | §3.5, Higham *FM* | 1 shot |
| Python `sqrtm(A)` | Denman–Beavers | Denman & Beavers (1976) | 20 |
| Python `logm(A)` | Inverse scaling-and-squaring + Gregory series | §11.3, Higham *FM* | 6 sqrts, 8 terms |
| Python `expm(A)` | Scaling-and-squaring + [6/6] Padé | §10.3, Higham *FM* | s = 10 |
| Fortran `sqrtm33z/logm33z/expm33z` | `eig33z` + reconstruction with `inv(V)` | local implementation | 1 eigensolve |
| `polar(F)` | `C = F^T @ F`, `U = sqrtm(C)`, `R = F @ inv(U)` | — | inherits `sqrtm` |

Higham *FM* = Higham, *Functions of Matrices* (2nd ed.).

---

## Rationale per function

### `eig(A)` / `eig33z` — trigonometric cubic, not Cardano, not Jacobi

Cardano's algebraic form involves complex cube roots with branch cuts
on the negative real axis.  A CS perturbation can push an intermediate
value across the cut, causing an eigenvalue to jump.  The
trigonometric form `λ_k = q + 2p·cos(φ + 2πk/3)` avoids cube roots
entirely; `arccos` and `cos` on complex arguments are analytic away
from their own cuts (which the spread guard handles).

Jacobi rotations were investigated and abandoned: computing the
rotation angle via `τ = (aqq−app)/(2·apq)` divides by a CS-perturbed
`apq` that can be nearly zero imaginary, producing division-by-zero or
catastrophic eigenvector error.  Direct 2×2 eigenstructure computation
avoids the `arctan` cancellation but still yields poor CS accuracy
(~100 % relative error on random matrices).

**Diagonal guard:** When `|p1|/|p2| < 1e−14` or `|p1| < 1e−18`, the
matrix is nearly diagonal and the trigonometric formula is unstable
(two eigenvalues nearly equal, `|r| ≈ 1`). We return the diagonal
elements directly. The absolute guard catches pure complex-step shear
perturbations where `p1` is `O(CS_H^2)` and the relative guard can be
too strict near the undeformed state.

The Fortran template mirrors the same guard and clamps the real part of
the cubic invariant `r` before complex arccos. These protections were
added after generated finite-strain UMATs produced singular eigenvectors
and NaN cascades for diagonal/near-diagonal matrices in `logm33z` and in
small symmetric tensile tangents.

Important boundary: this guard prevents a numerical NaN cascade; it is not a
general replacement for consistent spectral tangents at repeated eigenvalues.
For matrix functions such as `logm`, `sqrtm`, and `polar`, repeated eigenvalues
can often be handled with eigenspace-invariant formulas or iterative backends.
For constitutive laws that explicitly use principal directions, tensile
spectral splits, or slip-system construction, the real material update may
legitimately require eigenvectors. The tangent must then be formulated with
active-set or projector-aware logic, or with a deterministic basis selection
tied to the real state, because raw eigenvectors inside a repeated eigenspace
are not unique.

### Python `sqrtm(A)` — Denman–Beavers, always

The old hybrid code used eigendecomposition for well-separated
eigenvalues and fell back to Denman–Beavers only near degeneracy.
Denman–Beavers has excellent CS accuracy (~1e−15 relative) and needs
no `eig` call, so we use it as the sole path.  20 iterations is
overkill for value accuracy (convergence is quadratic) but the cost is
negligible for 3×3 matrices and the fixed count removes all branching.

### Python `logm(A)` — inverse scaling-and-squaring

Repeated square roots bring `B` toward `I`; the Gregory series
`log(B) = 2·(X + X³/3 + X⁵/5 + …)` with `X = (B−I)(B+I)⁻¹` converges
rapidly.  After 6 square roots, `B` is close enough that 8 terms give
function-value accuracy ~1e−14.

The cost is CS noise: each sqrtm step adds ~1e−15 to 1e−13 imaginary
garbage, and after 6 steps the accumulated error is visible in
derivatives.  See **Tuning** below.

### Python `expm(A)` — scaling-and-squaring + [6/6] Padé

Scale `A` by `2⁻¹⁰`, evaluate a [6/6] Padé approximant, then square
back 10 times.  The [6/6] approximant is accurate to machine precision
for `‖A‖ ≤ 0.01`, which the `s = 10` scaling guarantees for the input
ranges seen in plasticity (`‖dt·D_p‖ ≤ 10`).

CS accuracy is ~1e−05 relative — adequate for the tolerance budget
but the noisiest of the three matrix functions.  If margin erodes
during plasticity development, `expm` is the second place to look
after `logm`.

---

## Tuning: `logm` nscale 8 → 6

| nscale | function error | CS rel. error | % of 1e−6 tol |
|--------|---------------:|--------------:|--------------:|
| 6      | ~4e−14         | 8.72e−08      | 8.7 %         |
| 7      | ~1e−13         | 5.27e−07      | 52.7 %        |
| 8      | ~6e−13         | 1.11e−06      | 111 % **FAIL**|
| 10     | ~1e−12         | ~1e−05        | ~1000 %       |

A finite-strain UMAT verify at nscale=6 lands at 1.00e−07 → **10.0 %** of
budget.

**What the trade bought and what it cost:**

The old hybrid code (eigensolver for well-separated eigenvalues,
iterative fallback for degeneracy) gave ~2.1e−09 relative error on the
same test — about **0.2 %** of the 1e−6 tolerance.  The pure iterative
refactor (nscale=6) gives **8–10 %** of the tolerance.  That's a
~50× loss of CS safety margin, exchanged for a decomposition-free
implementation that never hits the eigenvalue-degeneracy
misclassification bug that motivated the refactor.

The mathematically optimal nscale for *value* accuracy is 8–10
(Higham).  The safe nscale for *derivative* accuracy in the CS
framework is 6.  The sweet spot is below the textbook optimum because
textbook optima optimize value accuracy, not derivative accuracy.

---

## Abandonment-of-decomposition: a design principle

When eigendecomposition-based matrix functions fail at degenerate
states, the response is to **abandon the paradigm** rather than patch
within it.

The specific failure mode that triggers the abandonment: the old
hybrid `eig` → `reconstruct(V, f(λ))` path produced spurious
eigenvalue splitting when two eigenvalues were nearly equal, causing
`sqrtm` and `logm` to return wrong values that misclassified material
mechanisms (brittle vs plastic).  Adding guards inside the
reconstruction path proved fragile; the robust fix was to remove the
eigendecomposition dependency entirely from `sqrtm`, `logm`, and
`expm` in the Python oracle.

The generated UMAT default has now completed that migration for matrix
functions: `generate_umat(..., matrix_backend="iterative")` is the default.
The eigendecomposition-based backend remains available by explicit request
(`matrix_backend="eig"`) for backward compatibility and debugging, but it is
not the recommended path for models whose tangents pass through diagonal or
near-degenerate matrix-function states.

## Explored Methods And Outcome

| Method | Where considered | Outcome |
|--------|------------------|---------|
| Algebraic Cardano cubic | Early `eig33z` design and scalar cubic-root equations | Kept only for scalar cubic roots; avoided for matrix eig due to complex cube-root branch cuts |
| Trigonometric depressed cubic | Python `eig`, Fortran `eig33z` | Current eigensolver |
| Jacobi rotations | Matrix eigensolver alternative | Rejected; CS perturbations around small off-diagonal terms produced unstable rotation angles |
| Direct 2x2 eigensystem fallback | Matrix eigensolver alternative | Rejected; poor CS derivative accuracy on random tests |
| Eigendecomposition matrix functions | Explicit generated Fortran option (`matrix_backend="eig"`) | Backward-compatible, but riskier near repeated eigenvalues |
| Iterative matrix functions | Current Python oracle and generated UMAT default | Preferred robust direction for matrix functions when repeated eigenvalues are a concern |

## Validation Coverage

The matrix functions are exercised by:

- Python value identities, `logm(F.T@F)` CS-vs-FD at degenerate,
  near-degenerate, and large-deformation states, `expm` CS-vs-FD, and eig
  reconstruction.
- Generated-Fortran-vs-Python-oracle comparison of finite-strain UMATs at
  nontrivial finite strain, including undeformed, tensile, small symmetric
  tension, and compressive-plastic states. The small symmetric-tension check
  confirms the generated shear tangent block remains nonzero and matches the
  Python reference.
- Known caveat: near plastic-yield thresholds with repeated principal
  stresses, Python and Fortran can choose different valid eigenvector bases
  and cross branch conditions at slightly different strains. Tests should use
  states safely inside the intended branch unless they are explicitly testing
  threshold behavior.

---

## What to watch when adding plasticity

Slip-system construction adds:
- Eigenvalue sorting (already present, but more frequent)
- Six dyad products `M = d ⊗ d`
- Six rate-law evaluations with exponent `1/m` at small `m` (e.g. `m = 0.02`)

Each operation amplifies the imaginary part.  If the CS budget
consumption climbs from the current 8–10 % toward 30–50 %, the
diagnostic order is:

1. **Suspect `logm` first.** Its inner `sqrtm` is the dominant CS
   noise source.  Check whether the regression appears on models that
   call `logm` frequently vs models that call `expm` frequently.
2. **Suspect `expm` second.** The [6/6] Padé is the least accurate of
   the three matrix functions in CS derivatives.
3. **Suspect the new plasticity code third.** Only after the matrix
   functions are ruled out should the slip-system construction be
   traced for branching-on-imaginary-part bugs.

Do **not** retune `nscale` as a first response.  If the matrix
functions were stable before the plasticity work and unstable after,
the cause is in the new code, not in the old knobs.

---

## Process discipline: reporting CS verification status

When reporting CS verification status, always give all four numbers:

1. Raw `err` (max absolute difference)
2. Raw `norm` (max absolute reference)
3. Tolerance value
4. Ratio `err / norm`

Restating only the ratio (or only the absolute error) makes it
impossible to catch unit mismatches across runs.
