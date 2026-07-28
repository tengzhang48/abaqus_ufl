# Lessons: Code Generation, Testing & DSL Design

Hard-won lessons on generating correct Fortran from a Python DSL and verifying
it end to end. Grouped by where the leverage is: emitting the code, proving it
correct, and designing the language the models are written in.

## Code generation

### Generate the entry point first in the generated file

Someone opening a generated `.for` file wants to see the solver interface (the
UEL subroutine) first, not hundreds of lines of utility functions. Fortran
allows calling subroutines defined later in the same file, so order the file
for the reader:

```
Section 1: UEL subroutine (solver entry point)
Section 2: Material subroutines (problem-specific)
Section 3: Complex-step tangent engine (problem-specific)
Section 4: Shipped templates (generic utilities)
```

### Sort tangent blocks deterministically

Python dict iteration is insertion-order, but tangent blocks built from set
operations or dict merges can vary in order between runs. Wrap them in
`OrderedDict(sorted(...))` so the generated tangent engine has the same
argument order every time. Non-deterministic output defeats diffing and
caching.

### Detect helper methods via the AST

When a material class defines helper methods (`_phi(self, F, p)`), the
generator must find and translate them. Scan each registered method's AST for
`self._xxx(...)` calls and emit the helpers before their callers.

Pitfall: inferring a helper's return type by calling it with test inputs fails
when the helper needs a specific physical state. Emit a `warnings.warn()` on
that fallback so the user sees which helper failed type inference.

### Node count comes from the element type, not the weak form

A weak form may infer node count from field degree (degree=2 → 8-node quad).
But if the user declares degree=1 and passes an explicit element (`Quad4`), the
weak form's node count can still be hardcoded for the higher-order topology.
Compute node count from the element-type parameter, not from the weak form.

### Single-letter Fortran names are landmines

Fortran is case-insensitive: `J` (a determinant) and `j` (a loop index) are the
same symbol, and a loop variable `a`/`b` collides with a UEL dummy argument
`A(NDOFEL)`. Python-level tests pass; the compiler rejects it. Never emit
single uppercase letters (I, J, K, L, M, N, A, U, V) as generated variable
names — use descriptive names (`Jdet`, `detJ`, `JdetInv`). Add a compile smoke
test that runs gfortran on generated `.for` files; structural Python tests
cannot catch compiler-level name clashes.

### Fortran array slices are inclusive — off-by-one is silent

In Fortran, `A(i:j)` includes both ends: the count is `j - i + 1`. So
`A(iof:iof+ns)` passes `ns+1` elements and silently reads/writes one element of
the adjacent block; the correct slice for `ns` elements is `A(iof:iof+ns-1)`.
Such a bug hides for years when examples use zero state variables — with `ns`
small or dummy the overflow is harmless — and only corrupts adjacent state once
a real stateful element runs. When you fix one instance, grep for the pattern
and fix every copy: the same off-by-one tends to be copy-pasted across many
sites in several files.

## Testing & verification

### The finite-difference tangent check is the primary automated check

If `AMATRX ≈ -dRHS/dU` within FD accuracy (~1e-6), the tangent is consistent
with the residual across the pipeline: DOF parsing, field interpolation,
material evaluation, tangent assembly, and the sign of `AMATRX` relative to
`RHS`. It verifies the tangent is consistent with the residual, but it does
**not** prove the physics, the residual, or the PDE sign convention is correct —
a residual can be internally tangent-consistent yet encode the wrong equation.
Those need an independent quantitative benchmark. An optional Abaqus
integration check supplies production-interface, convention, and output-bridge
evidence; it is not automatically a physics oracle (see the pipeline below).

```python
for col in range(NDOFEL):
    U_p = U.copy(); U_p[col] += eps
    U_m = U.copy(); U_m[col] -= eps
    RHS_p, _ = assemble(U_p, ...)
    RHS_m, _ = assemble(U_m, ...)
    fd_col = -(RHS_p - RHS_m) / (2*eps)
    assert np.allclose(fd_col, AMATRX[:, col], rtol=1e-4)
```

### Test at non-trivial deformation states

Testing at F = I hides bugs because many tangent terms vanish there. Always
test at a random deformation state with every field nonzero.

### Analytical benchmarks for simple cases

Uniaxial tension on a single element has a closed-form reaction force. Compare
both the Python assembly and the directly executed generated element against
it. The first checks the formulation/assembly; the second adds code-generation
and compiled-interface evidence.

### Build a Python reference assembly

Write a Python function that reimplements the intended element assembly — DOF
parsing, shape functions, and assembly loops — and use it for fast residual,
tangent, and bookkeeping checks during development. It does not execute the
generated source, so code-generation parity still requires a compiled UEL
comparison. Because both paths implement the same formulation, neither is an
independent physics oracle; use a closed form, hand calculation, published
benchmark, or independently formulated comparison for that gate.

### Determine matrix symmetry from the formulation

Do not infer symmetry from the number or names of fields. A variational
saddle-point formulation can be symmetric indefinite, while time-dependent,
history-dependent, stabilized, or non-associated coupling can be genuinely
nonsymmetric. Inspect the off-diagonal blocks and measure
`||K-K.T||/||K||` at representative states. Request unsymmetric tangent
treatment only when the actual operator requires it. A solver symmetry mismatch
can corrupt Newton corrections even while a per-entry stiffness check passes,
so compare the declared treatment with the measured operator before changing
the physics or iteration limits.

### A passing stiffness check does not guarantee Newton convergence

A per-entry stiffness check confirms ∂R/∂u is correct at one equilibrium point.
It can pass while Newton still fails, for three reasons:

1. **Solver symmetry** — the check tests entries independently and never sees
   whether the global solve treats the matrix as symmetric (see above).
2. **Uninitialized time step** — rate-dependent tangents scale as 1/dt. If the
   time increment is unset during the stiffness check, both the analytical
   tangent and the FD derivative blow up identically (~1e13) and agree, but the
   values are meaningless. Set a nonzero time increment before the stiffness
   check for rate-dependent problems.
3. **State-update bugs** — the check perturbs DOFs at a single instant; it
   never advances time or commits state variables, so state read/write and
   reset bugs only surface during a real Newton solve.

Always follow a passing stiffness check with an actual Newton solve — that is
the real test. When the check passes but Newton fails, the bug is in the solver
infrastructure (symmetry, BCs, time stepping), not in the element tangent.

### Add one feature at a time — test simpler before complex

Jumping from a working single-field element straight to a full multi-field
formulation leaves several candidate causes suspected at once, and debugging
them simultaneously is combinatorial. Add one feature at a time and test at
each step, e.g.:

1. Single-field displacement (Quad4/Quad8) — basic assembly
2. Single-field displacement (Hex8/Hex20) — 3D extension
3. Two-field (u, p) mixed — multi-field assembly
4. Three-field (u, p, μ), mechanical only — pressure + coupling
5. Three-field with diffusion — full formulation

When step N breaks, the bug is in what changed between N-1 and N. N separate
tests cost linearly; debugging N features at once costs combinatorially.

### Combine systematic search with domain judgment

Different bug classes need different tools. Mechanical, systematic bugs —
off-by-one slices copy-pasted across files, case-insensitive name collisions,
platform-specific defaults, duplicate declarations — are found fastest by
automated search: grep across the tree, compile tests, AST audits. Conceptual
bugs — an unsymmetric-solver flag, an integration-scheme choice, a
mixed-element DOF layout, which test sequence to run — are found fastest by
domain expertise. Use systematic audits to eliminate the mechanical failures so
the remaining problem is the conceptual one. The most expensive debugging hours
go to problems that domain instinct settles in seconds, but only after the
mechanical causes have been ruled out.

## DSL design

### Workarounds in model code usually mask framework gaps

When model code starts post-processing generated `.for` files, inlining helpers
to dodge a limitation, or sprinkling projections in unexpected places, ask "is
this a model-side rewrite or a framework gap?" — not "how do I make the
workaround robust?" Workarounds in the build script calcify and spread
silently; framework fixes are additive and reusable across all materials.
Heuristic: if the workaround would have to be re-applied to every new model in
the same class, the fix belongs in the framework. Gaps a workaround can hide:

- A regex pass that de-duplicates declarations masks a missing "already
  declared" check in the main renderer that the helper renderers already have.
- Inlining tuple-returning helpers masks a translator that only accepts
  specially-tagged helpers; plain Python helpers returning tuples should
  translate the same way.
- Rewriting `if …: return x` ladders masks a generator that silently overwrites
  earlier returns (see below).

### Stop and surface a missing primitive

When the translator can't express a piece of math, the right response is to
stop and report it, not to improvise a workaround. Missing primitives resolve
into small, additive, model-independent extensions (`dyad(a, b) → tensor(3,3)`,
`sym3(...) → tensor(3,3)`, each ~20 lines) that expand what every future
material author can express. A hand-written Fortran sidecar, by contrast, hides
a reusable primitive (e.g. a spectral solve) as tribal knowledge for every
model that later needs the same math. One stop-and-surface hand-off is far
cheaper than a model-side workaround that lives forever.

### Silent code drops are the worst failure mode

A translator that tracks a single return-expression field overwrites it on each
`return` it visits, so

```python
if x.real > 0:
    return a
return b
```

translates as if only `return b` existed — a real branch silently vanishes,
caught only when results diverge numerically. Rule: when the translator can't
represent a construct, it must raise. Silently dropping code is far worse than
an ugly warning or a refusal to compile. The fix is a one-pass AST check that
rejects `return` nested inside control flow, with a clear "assign-and-
return-once" remediation message.

### "I want a 6×6" usually means "I want a tensor equation"

A 6×6 Voigt system is often a 3×3 tensor equation in disguise. For example
`sym(F·ΔεP + ΔεP·F) = Y` is a Sylvester equation solvable spectrally with no
Voigt vectors and no subscript assignment:

```python
lam, V = eigh(F)
Vinv = inv(V)              # eig(F) does not promise an orthonormal basis:
                           # use inv(V), not V.T
Y_eig = Vinv @ Y @ Vinv.T
ones3, _ = eigh(eye(3))
M = 0.5 * (dyad(lam, ones3) + dyad(ones3, lam))
X_eig = Y_eig / M
deps_p = V @ X_eig @ V.T
```

The 3×3 spectral form is shorter than the Voigt form and translates through the
DSL with no framework change once `dyad` exists. Mind the `inv(V)` vs `V.T`
trap: a `V.T` form can appear to work for diagonal or near-diagonal tensors
because the eigenvectors happen to be near-orthonormal there, and fails once
the tensor rotates. Always test spectral linear algebra on a rotated,
non-diagonal tensor and check the residual of the tensor equation, not just the
downstream stress curve.

Idiom map: `np.einsum('ij,ij->', a, b)` → `trace(a.T @ b)`;
`np.einsum('ij,jk->ik', a, b)` → `a @ b`; `np.zeros((3,3))` → `0.0 * eye(3)`; a
constant symmetric tensor → `sym3(...)`.

### Multi-output helpers replace state dataclasses

A DSL with a narrow kind vocabulary (`scalar`, `vector(3)`, `tensor(3,3)`)
can't translate helpers that take or return a dataclass. Replace the dataclass
with the helper's calling convention: flat tensor/scalar arguments in, a tuple
of tensors/scalars out.

```python
def _step(self, eps_e_old, eps_r_old, F_old, cg_old, deps, dt):
    ...
    return eps_e_new, eps_r_new, F_new, cg_new, sigma_new
```

The translator emits a Fortran subroutine with the inputs as `INTENT(IN)` and
the outputs as trailing `INTENT(OUT)` arguments; the caller unpacks with tuple
assignment. The math stays readable in Python while staying inside the DSL's
type system.

### Validation envelopes must cover the branch being claimed

A model can have many green tests and a generated routine that compiles and
runs, yet leave the claimed physics untested. When the tested histories all
stay in an elastic or near-diagonal envelope, the defects live in the branch
that never fired: rotated fabric, nonzero residual strain, split yield
intervals. A common cause is the default state — if the initializer seeds state
at a value that drives the interesting branch to zero, the verify path only
exercises the elastic wrapper. Rule: every claim needs a matching validation
state. For a "full plastic" claim, include at least one verification path that
activates the plastic branch, seeds the nonzero internal variables so their
terms are live, uses a non-diagonal rotated tensor so spectral-basis
assumptions are tested, and exercises multi-interval logic if the model
supports it. Tests outside that envelope are still useful, but report them as
basic-path validation, not full-model validation.

### State initialization needs one source of truth

When a script-side dataclass and the generated routine seed state differently —
one seeds a variable with a tiny isotropic value to avoid a residual-stiffness
singularity, the other leaves it zero — the two paths start from different
mathematical states and diverge (one guards the zero away and never leaves it;
the other produces NaNs from the same zero). Any state seed required for domain
validity belongs in the generator-visible state layout, and the script-side
constructor should derive from that same source. Tests must round-trip the
initial state through the generated state layout before validating physics.

### Independent review is most valuable attacking the validation envelope

Multi-agent review pays off less as extra code and more as adversarial review
of the validation envelope. The high-value question is not "are the formulas
right?" but "which branches have the acceptance tests never exercised?" — that
is what exposes defects in regimes the green tests miss. A workable division of
labour:

- one agent implements against the current scope;
- a second agent reviews the validation envelope and asks which branches are
  unreached;
- the implementing agent fixes only confirmed defects and adds adversarial
  regressions;
- a review agent re-checks whether any remaining mismatch is a model error, a
  codegen/reference parity gap, a calibration issue, or a solver issue.

The reviewer should design counterexamples, not just read formulas:
non-diagonal tensors, nonzero seeded state, split intervals, non-homogeneous
gradients, near-singular limits, and loading paths that activate the claimed
physics. The payoff comes when the models disagree for specific, testable
reasons — and the honest output is a sharper status ("these defects fixed; full
parity still open"), not a premature all-clear.
