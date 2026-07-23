# Lessons: Complex-Step Differentiation & CS-Safe Matrix Functions

Complex-step (CS) differentiation computes a derivative from a single
function evaluation on a complex-perturbed input:

```
f'(x) ≈ Im[ f(x + i·h) ] / h
```

Because the perturbation lives entirely in the imaginary axis, there is no
subtractive cancellation — the result is accurate to machine precision for
any sufficiently small `h`, unlike a finite difference. This makes CS an
excellent way to obtain consistent constitutive tangents (`dP/dF`, material
Jacobians) directly from a scalar-valued material model, provided every
operation the perturbation flows through is **holomorphic** (complex
analytic). The entire discipline below is about keeping that path
holomorphic and keeping the imaginary part clean.

Foundational references for the technique: Squire & Trapp (1998), *Using
complex variables to estimate derivatives of real functions*; Martins,
Sturdza & Alonso (2003), *The complex-step derivative approximation*.

---

## Complex-step best practices

### Perturbation size (`CS_H`)

`CS_H = 1e-10` is a robust default. Because CS has no subtractive
cancellation, any small `h` gives machine-precision tangents — verify this
by confirming that `1e-10` and `1e-30` produce identical tangents across all
tangent blocks. `1e-10` is preferred over smaller values to avoid subnormal
issues in intermediate computations.

**Rule:** use `CS_H = 1e-10`, and document the choice so reviewers do not
mistake it for a finite-difference step size.

### Reset all complex state before each perturbation

**The trap:** if complex working variables are reused across perturbation
blocks without being reset, the imaginary part left over from a previous
perturbation leaks into the next block and produces cross-contaminated
tangents.

**Rule:** at the start of *every* perturbation, rebuild all complex
variables from their real counterparts, then perturb exactly one component:

```fortran
      DO k = 1, 3
        DO l = 1, 3
          CALL real2complex33(F, Fz)               ! fresh copy every time
          pz = DCMPLX(p, 0.0d0)                     ! fresh
          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)   ! perturb ONE component
          CALL material(Fz, pz, ..., Pz)
          dPdF(:,:,k,l) = AIMAG(Pz) / CS_H
        END DO
      END DO
```

### Use `AIMAG`, not `DIMAG`

`DIMAG` is a compiler extension. `AIMAG` is standard Fortran and works
correctly with `DOUBLE COMPLEX`. Use `AIMAG` to extract the imaginary part.

### Keep `dt` in double precision

The time increment `dt` is never differentiated with respect to. Passing it
as `DOUBLE COMPLEX` wastes the imaginary part and risks a subtle bug if that
imaginary part leaks into a rate computation.

**Rule:** declare `dt` (and other genuinely-real scalars) as
`DOUBLE PRECISION` in material-subroutine signatures, not `DOUBLE COMPLEX`.

### CS-safe branching rules

Never branch on the imaginary part of a complex variable. Branch on the
**real part only** so that the perturbed and unperturbed evaluations follow
the same code path.

| Operation | CS-safe? | Alternative |
|-----------|:--------:|-------------|
| `+`, `-`, `*`, `/`, `**` | yes | |
| `LOG`, `EXP`, `SQRT` | yes | Standard intrinsics |
| `SIN`, `COS` | yes | |
| `det(A)`, `inv(A)` | yes | Explicit formula |
| `A @ B`, `A.T` | yes | |
| `ABS(x)` | no | Use `x` directly, or `SQRT(x*x)` |
| `MAX(a,b)`, `MIN(a,b)` | no | Branch on `DBLE(a) > DBLE(b)` |
| `SIGN(x)` | no | Branch on `DBLE(x) > 0` |
| `IF (x > threshold)` with complex `x` | no | `IF (DBLE(x) > threshold)` |

`ABS`, `MAX`, `MIN`, and `SIGN` are non-holomorphic because they encode a
branch on magnitude or sign — exactly the kind of comparison that must be
taken on the real part only.

### Defensive `.real` casts are a CS leak

**The trap:** when a complex value must feed a real-only operation (a
comparison, a `max`, a Macaulay bracket), the natural defensive move is to
cast to real early. This silently breaks the chain rule:

```python
tau = float((s.T @ M @ m).real)   # kills the derivative here
if tau > 0.0:
    gamma = (tau / C) ** (1.0 / m)
```

The imaginary part of `s.T @ M @ m` carries the derivative of the resolved
quantity with respect to the perturbation. Stripping it at the source means
everything downstream never sees the perturbation.

**Correct pattern:** branch on `.real`, but compute on the complex value:

```python
tau = s.T @ M @ m                    # keep complex
if tau.real > 0.0:                   # branch on real part only
    gamma = (tau / C) ** (1.0 / m)   # gamma is complex, derivative flows
else:
    gamma = 0.0 * tau                # preserve complex type
```

The same rule applies to *every* intermediate that participates in a branch
(back-stress, denominators, rate quantities). If any one of them is cast to
real before the branch, a derivative leak is introduced at that point and
propagates silently.

**Diagnostic:** if removing a `.real` cast does not change the verification
budget number, the test states are simply not exercising the sensitive path.
That is useful information about test coverage — but it does not make the
cast harmless. The cast is still wrong; the test is just not catching it.

**Rules:**
- Treat every `.real` cast as a bug until proven otherwise. The proof is
  that the downstream operation genuinely requires a real scalar and has no
  physical derivative through this variable.
- Substepping counts, convergence tolerances, iteration limits, array sizes,
  and print statements are legitimate real-only uses. Constitutive
  quantities (stress, strain, slip rate, damage) are not.

### Iterative methods compound CS error linearly

**Rule:** every iteration step in a complex-step evaluation injects roughly
`1e-15` to `1e-13` of spurious imaginary noise. After `N` iterations that
error accumulates roughly linearly with `N`. Function-value accuracy (more
iterations is better) and CS-derivative accuracy (fewer iterations is
better) pull in opposite directions, so the CS sweet spot is often *below*
the textbook iteration count.

**Concrete cost — `logm` via inverse scaling-and-squaring**, as a function of
the squaring count `nscale`, against a `1e-6` tolerance:

| nscale | function error | CS rel. error | fraction of tol |
|--------|---------------:|--------------:|----------------:|
| 6      | ~4e-14         | 8.72e-08      | 8.7 %           |
| 7      | ~1e-13         | 5.27e-07      | 52.7 %          |
| 8      | ~6e-13         | 1.11e-06 FAIL | 111 %           |
| 10     | ~1e-12         | ~1e-05        | ~1000 %         |

A second, independently demanding state at `nscale = 6` lands near
`1.00e-07`, about 10 % of the budget. More function-value accuracy (higher
`nscale`) directly *degrades* the derivative.

**The robustness/accuracy trade.** A hybrid implementation — eigensolver for
well-separated eigenvalues, iterative fallback only for degeneracy — can
reach `~2e-09` relative error (about 0.2 % of a `1e-6` tolerance). A pure
fixed-count iterative implementation (`nscale = 6`) sits at 8–10 % of the
tolerance. That is roughly a 50× loss of CS safety margin, bought in
exchange for a decomposition-free implementation that never hits the
eigenvalue-degeneracy misclassification bug that eigensolver backends carry
(see the matrix-functions half). Both halves of the trade are real: the lost
margin is real, and the robustness gain is real. Choose deliberately and
record the reasoning.

**Margin-watch discipline.** Operations that amplify the imaginary part —
eigenvalue sorting, chains of dyad products, rate laws with small exponents
(e.g. `x ** (1/m)` at `m ≈ 0.02`) — each erode the CS margin. If budget
consumption climbs (say from ~10 % toward 30–50 %) as a model grows, that is
a signal to **root-cause where the new amplification comes from**, not to
tweak the iteration count. The first response to margin erosion is
diagnosis, not knob-turning.

**Reporting discipline — always state four numbers.** When reporting CS
verification status, give all of:

1. Raw `err` (max absolute difference),
2. Raw `norm` (max absolute reference),
3. the tolerance value, and
4. the ratio `err / norm`.

Restating only the ratio (or only the absolute error) makes it impossible to
catch unit mismatches, and hides the difference between a small absolute
error at a small reference and a genuine failure.

**Rules:**
- Document the iteration-count budget in the routine's docstring.
- Test CS derivatives at the *most demanding* state the model will
  encounter, not just the reference state (`F = I`).
- Periodically tighten the test-suite tolerance (e.g. to `1e-10`) so margin
  erosion is caught before it becomes a hard failure.

### Verification can surface constitutive-law ambiguities, not just code bugs

**The trap:** when a verification check fails at a particular state, the
default assumption is that a numerical bug has been found — a CS leak, an
overflow, a wrong sign. But a verification failure can also mean the
constitutive law itself is ambiguous at that state.

Two effects are easy to confuse:

1. **A backend/code bug** — reproducible, state-independent, fixable with a
   single change. Example: at diagonal or nearly-diagonal states, a matrix-log
   routine that calls an eigendecomposition helper hits the helper's
   near-diagonal guard, receives `V = I`, and zeroes the shear tangent (see
   the matrix-functions half). Switching the matrix functions to an iterative
   backend restores agreement with the reference to the full tolerance.

2. **A genuine constitutive ambiguity** — state-dependent, tied to a physical
   regime, with no single-line fix. Example: at a yielded, nearly
   axisymmetric state where two principal stresses are nearly equal, a model
   may construct discrete slip systems from the principal directions of a
   stress measure. In the degenerate plane the eigenvectors are *arbitrary* —
   any orthonormal basis is valid — so different valid bases yield different
   slip systems and different plastic responses unless the model defines an
   active-set or tie-break rule. CS and FD can then disagree not because of a
   code bug but because the response is genuinely multivalued:

   | Component | CS derivative | FD derivative | Relative error |
   |-----------|--------------:|--------------:|---------------:|
   | axial     | 1.960×10¹⁰    | 1.960×10¹⁰    | 1.1×10⁻⁶       |
   | in-plane  | 5.338×10⁹     | 5.277×10⁹     | 1.2×10⁻²       |

   The axial component matches to tolerance while the in-plane component does
   not — the signature of a degenerate-eigenspace ambiguity, not a leak.

**Rules:**
- When verification fails, first decide whether the error is a backend/code
  bug (reproducible, single fix) or a constitutive ambiguity (state-dependent,
  physical regime, no single-line fix). Eliminate backend leaks (such as
  eig-based matrix-function guards) *first*, then retest the sensitive state.
- For a genuine ambiguity, document the regime, quantify the error, and defer
  reformulation until a full simulation shows the ambiguity causes a visible
  problem.
- Do not relax verification tolerances to make the failure disappear. The
  verification did its job; the model has a known limitation.

### Generator-inserted constructs can break CS even when the source is safe

**The trap:** a Python material model can have perfect CS hygiene — no
`.real` casts on constitutive quantities, no branching on complex values, no
`abs()` of complex intermediates — and *still* produce a broken tangent,
because a code generator or a called backend routine inserts a CS-breaking
construct that does not exist in the source.

A CS leak can enter at three levels:

1. **The user's source** — `.real` casts, `abs()`/`max()` on complex values.
2. **The generator's translation** — type mismatches, missing declarations,
   defensive casts inserted around comparisons or loop bounds.
3. **The backend subroutines** called by the generated code — guards,
   branching on tiny values, convergence tolerances hidden inside a utility.

Category 3 is the hardest to catch, because the leak is invisible at the call
site: the generated code reads `CALL logm(...)`, and the developer assumes
the utility is correct. The canonical example: a matrix-log routine that
internally calls an eigendecomposition helper whose near-diagonal guard
returns the identity for a CS-perturbed diagonal input, zeroing the shear
block of the tangent (detailed in the matrix-functions half).

**Rules:**
- Treat every backend subroutine that contains a conditional on a complex
  value as potentially CS-unsafe until proven otherwise. The proof is that
  the conditional is invariant under the complex-step perturbation, or that
  it is a conditional on a genuinely real quantity (iteration count, array
  size).
- For matrix functions, prefer fixed-count iterative backends over
  eigendecomposition backends in generated code, because eigenvector guards
  are inherently perturbation-scale-dependent.
- Add regression tests that check tangent *correctness* — not just finiteness
  — at states where the unperturbed matrix is diagonal or nearly diagonal.

### Reject unbounded `while` loops in generated code

**The trap:** iterative solves (Newton, secant) are naturally written as
`while` loops. In a Fortran UMAT/UEL an unbounded `while` that fails to
converge infinite-loops and hard-crashes the host FE solver with no error
message.

**Design rule:** a CS-safe code generator should reject unbounded `while`
loops and only accept bounded iteration — `for i in range(N)` with an early
`break`. This forces every loop to terminate. Rewrite:

```python
# Rejected: unbounded
while abs(f) > tol:
    df = (f1 - f0) / (x1 - x0)
    x_new = x1 - f1 / df
    x0, x1 = x1, x_new
    f0, f1 = f1, f(x_new)

# Accepted: bounded, terminates even if the step wanders
for _ in range(20):
    df = (f1 - f0) / (x1 - x0)
    x_new = x1 - f1 / df
    x0, x1 = x1, x_new
    f0, f1 = f1, f(x_new)
    if abs(f1) < tol:
        break
```

The bounded form is behaviorally identical for well-posed problems but
cannot hang.

### IF-comparison right-hand sides must be plain real

**The trap:** a generator that wraps every Python float constant in
`DCMPLX(value, 0.0d0)` produces type-mismatched comparisons:
`IF (DBLE(f) .GT. DCMPLX(0.0d0, 0.0d0))` — a `DOUBLE PRECISION` left side
against a `DOUBLE COMPLEX` right side.

**Fix:** when the comparator is a constant, emit a plain real literal instead
of wrapping it, giving `IF (DBLE(f) .GT. 0.0d0) THEN`. Comparisons must be
real-against-real.

### Integer loop bounds must not be complex-wrapped

**The trap:** blindly complexifying every literal turns `for i in range(20)`
into `DO i = 1, INT(DBLE(DCMPLX(20.0d0, 0.0d0)))` — three nested conversions
for a plain integer.

**Fix:** when a `range()` argument is a constant, emit the integer directly
(`DO i = 1, 20`). Reserve `INT(DBLE(...))` for variable bounds.

### The declaration scanner must skip loop constructs

**The trap:** a declaration scanner that discovers user variables by looking
for assignment (`=`, not `CALL`) will match `DO iteration = 1, 20` and emit a
spurious `DOUBLE COMPLEX :: DO iteration`.

**Fix:** extend the scanner's skip list to include loop keywords
(`DO`, `END DO`, `EXIT`) alongside the conditional keywords.

### Loop variables must be `INTEGER`, not `DOUBLE COMPLEX`

**The trap:** a Fortran `DO` variable must be `INTEGER`. If the scanner does
not handle it, the loop variable is either declared `DOUBLE COMPLEX` (type
error) or left undeclared (compile error).

**Fix:** after the main declaration block, scan for `DO varname = ...`
patterns and emit `INTEGER :: varname` for each. Loop counters are integers,
never part of the complex-differentiated state.

### State-output variables: no duplicate declarations, no self-copies

**The trap (duplicate declaration):** when a model returns an updated state
variable (e.g. `{'Fp': Fp_new}`), the translator sees `Fp_new = ...` in the
body and declares it `DOUBLE COMPLEX`, but `Fp_new` is also an `INTENT(OUT)`
argument — a duplicate declaration that fails to compile.

**Fix:** build the set of state-output names *before* the declaration loop
and skip any name already declared as a subroutine argument.

**The trap (no-op self-copy):** when the state-return map is
`{'ep': ep_new}`, a naive generator emits `ep_new = ep_new` — harmless but
confusing.

**Fix:** in the state-copy section, skip the copy when the value expression
is exactly `<name>_new`.

---

## CS-safe matrix functions

Constitutive models routinely need matrix functions of a symmetric tensor:
`logm`, `expm`, `sqrtm`, the polar decomposition, and eigenvalues. Each has a
CS trap around eigenvectors, branch cuts, or degeneracy. The unifying safe
strategy is: use holomorphic reconstruction where possible, and fall back to
**fixed-count iterative** methods (matmul and `inv` only) at degeneracy.

### Eigenvalues are CS-compatible; eigenvector guards are the trap

Eigen*values* are well-behaved under CS. In fact a complex perturbation often
*breaks* a real degeneracy: even when the real matrix has repeated
eigenvalues (e.g. `C = I` at the reference state), the perturbed matrix
`C + i·h·e_kl` can carry the derivative information CS needs. This is why
eigendecomposition can be useful inside a CS framework, whereas a symbolic
tangent containing `1/(λ_i − λ_j)` blows up exactly at repeated eigenvalues.

The trap is the **eigenvectors**. This works only if the eigensolver actually
propagates the perturbation. A near-diagonal guard inside a `3×3`
eigendecomposition helper — for instance `ABS(p1) < 1.0d-18` returning the
identity for inputs that "look diagonal" — silently destroys the derivative:
for a `1e-10` complex perturbation of a diagonal matrix, the off-diagonal is
`~1e-20`, below the guard, so the helper returns `V = I`. The eigenvector
derivative is zero, and any downstream matrix function (`logm`, `sqrtm`,
`expm`, `polar`) loses the shear tangent (`C44 = C55 = C66 = 0`).

**Distinction that matters:** an explicit `eig()` call inside a model may
still be legitimately required — for principal-stress logic or slip-system
construction — and those calls need active-set/eigenspace care near repeated
eigenvalues. Matrix *functions* are different: they should not depend on
arbitrary eigenvectors at diagonal states. Prefer fixed-count iterative
matrix-function helpers as the generated default so this backend leak cannot
occur.

### Eigenvector normalization is non-holomorphic

**The trap:** the standard normalization `v / sqrt(v·v)` breaks CS. `sqrt`
carries a branch cut, and dividing by a complex function of the eigenvector
components is not analytic.

**Fix:** return *raw, unnormalized* eigenvectors from the decomposition and
reconstruct with `inv(V)` rather than `Vᵀ`:

```
f(A) = V @ diag(f(λ)) @ inv(V)
```

`inv(V)` is holomorphic (cofactor formula) and the arbitrary column scaling
cancels between `V` and `inv(V)`.

**Corollary:** `Vᵀ = V⁻¹` holds for real symmetric matrices with normalized
eigenvectors and *approximately* holds for CS-perturbed matrices — but
relying on that approximation is fragile. Use `inv(V)` explicitly.

### Conjugation is not holomorphic

**The trap:** a Hermitian inner product (e.g. `np.vdot(v, v)`) conjugates its
first argument. Conjugation is not complex-differentiable, so any CS path
through it destroys the tangent.

**Fix:** use the **bilinear** inner product without conjugation —
`v[0]*v[0] + v[1]*v[1] + v[2]*v[2]` — not the sesquilinear (Hermitian) form.
This returns a complex number rather than a real one, which is exactly what
CS requires.

### Cardano's algebraic form has branch-cut issues; use the trigonometric form

**The trap:** for a real symmetric `3×3` matrix all three eigenvalues are
real, but Cardano's *algebraic* formula routes through complex cube roots
that cancel to a real result. The complex cube root `exp(log(z)/3)` has a
branch cut on the negative real axis, so an eigenvalue can jump when the CS
perturbation pushes an intermediate value across the cut.

**Fix:** use the *trigonometric* form of the cubic solution:

```
λ_k = q + 2·p·cos(φ + 2πk/3),   φ = arccos(r) / 3,   k = 0, 1, 2
```

This avoids complex cube roots entirely. On complex arguments during CS,
`arccos` and `cos` are analytic (away from the `arccos` branch cuts, which a
degeneracy/spread fallback handles).

### Eigenvalue degeneracy needs an iterative fallback

**The trap:** at `F = I`, the right Cauchy–Green tensor `C = I` has triply
degenerate eigenvalues. Eigenvector computation via cross products of the
rows of `A − λI` yields zero vectors, and reconstruction fails.

**Fix:** measure the *relative eigenvalue spread*. If the spread is below a
threshold, fall back to fixed-count iterative matrix functions that use only
`matmul` and `inv` — no eigendecomposition — and are perfectly CS-safe
because they never branch on matrix values:

- `sqrtm`: Denman–Beavers iteration (fixed count, e.g. 30 iterations);
- `logm`: inverse scaling-and-squaring plus a Gregory series;
- `expm`: scaling-and-squaring plus a `[2/2]` Padé approximant.

**Threshold guidance:** a looser threshold (`~1e-6`) is appropriate where the
iterative methods themselves accumulate CS noise over many iterations; a
tighter threshold (`~1e-10`) is appropriate for more efficient
implementations. Note the tension with the iteration-count budget above:
more iterations improve the function value but erode the CS derivative.

### Polar decomposition must use `inv(U)`, not eigenvector reconstruction

**The trap:** computing `U_inv` via `reconstruct(V, 1/sqrt(λ))` re-decomposes
the tensor — it calls `eig(C)` a second time and uses the eigenvectors
directly, bypassing the degeneracy safeguard *even when* `sqrtm(C)` correctly
fell back to the iterative method. The polar rotation then re-inherits the
degenerate-eigenvector trap that `sqrtm` was careful to avoid.

**Fix:** reuse the already-safe `U = sqrtm(C)` and form

```
R = F @ inv(U)
```

Do not re-decompose. `inv(U)` is holomorphic and preserves whatever
degeneracy handling `sqrtm` already applied.

### Mixed absolute/relative tolerance for tangent verification

**The trap:** verifying a tangent with a *relative* error `err / norm` alone
fails when the true tangent is zero. The reference `norm` is `~0`, so the
relative error diverges even though the absolute error is `~1e-15`.

**Fix:** use a mixed tolerance. If `norm > atol` (e.g. `atol = 1e-12`),
compare relatively; otherwise compare the absolute error against `atol`. This
is the companion to the four-number reporting discipline: keep both the
absolute and relative views so a genuinely-zero tangent is not mistaken for a
failure, and a small-reference near-miss is not mistaken for success.
