---
name: abaqus-ufl
description: Use when developing, reviewing, validating, or debugging abaqus_ufl UMAT/UEL examples, generated Abaqus Fortran, f2py checks, Abaqus validation decks, or paper-to-code model implementations. Guides agents through target selection, generator-friendly Python, validation ladders, and known Abaqus/Fortran pitfalls.
metadata:
  short-description: Develop and validate abaqus_ufl UMAT/UEL models
---

# abaqus_ufl Model Development

This skill is a compact operational guide distilled from the repository's
tests, examples, and lessons. Use it before implementing or reviewing any
`abaqus_ufl` UMAT, UEL, f2py, or Abaqus-validation work.

## First Moves

1. Check the worktree:
   ```bash
   git status --short
   ```
   Do not overwrite collaborator edits, generated files, or Abaqus-validation
   artifacts without understanding them.

2. Read the live API guide in the repo:
   - `docs/API_USAGE.md`

3. Pick the target deliberately:
   - **finite-strain UMAT**: local constitutive law, built-in Abaqus elements
     are enough, model returns `stress_PK1(F)`.
   - **small-strain UMAT**: incremental material update with local `STATEV`;
     use `SmallStrainMaterial.stress_update(...)`.
   - **UEL**: any extra solved nodal field, gradient term, mixed interpolation,
     phase field, pore pressure, diffusion, temperature, concentration, or
     chemical potential.
   - **f2py**: point/element oracle for generated Fortran before solver runs.
   - **Abaqus validation**: final convention and production-solver check.

## Core Rule

Do not start from Fortran. Start from the theory, write generator-friendly
Python, verify the Python model, generate Fortran, compile, then climb the
validation ladder.

The validation ladder is:

1. Python reference checks.
2. `model.verify()` or `problem.verify()`.
3. Material-point or element-level regime checks for every claimed branch.
4. Generated Fortran compile test.
5. f2py single-point or single-element test when the model is nontrivial.
6. Abaqus one-element validation deck.
7. Published-figure reproduction only after the previous rungs are stable.

## Verify, Do Not Guess

When a build fails, a solve diverges, or you suspect a framework limitation,
**test the hypothesis before acting on it**, and check the theory, the API, and
the existing examples first. Most "the framework can't do X" or "it must be Y"
beliefs are wrong:

- **Check the examples before claiming a limitation.** A local Newton / implicit
  return map *inside* the element is fully supported (branch on `.real`, a fixed
  `for k in range(N)` loop, analytic Jacobian; the complex step rides through the
  converged iterate and yields the consistent tangent). For STIFF exponential rate
  laws (`exp(x/C)`, small C) at implicit dt, solve in log-slip variables (exact
  backward-Euler of the smoothed law, no active-set branching, early-exit `break`
  after two consecutive machine-zero residuals so the CS imaginary part settles).
  `grep` the examples before concluding something is impossible.
- **Localize a failure with a test, not a guess.** Reproduce the exact failing
  state; cross-check the compiled element against `reference_assembly`; confirm a
  suspected culprit actually misbehaves *at that state* before fixing it. In one
  session eig and the elastic/plastic branch were both wrongly blamed: eig was
  finite at the degenerate states (Python *and* Fortran), and the branch is
  `.real`-invariant under the imaginary complex-step perturbation. The real cause
  was a formulation choice (the wrong unknown in a return map — see pitfalls).
- **Follow the theory.** A non-smooth/divergent local solve is usually the wrong
  *unknown* or formulation, not a solver-tuning problem. Reach for the standard
  computational-plasticity construction, not an ad-hoc damper or step limiter.
- **A NaN is a bug to find, not a number to compare.** `NaN >= tol` is False, so a
  NaN slips silently through a convergence check; guard residual checks with
  `isfinite`. (When diagnosing, note an f2py recompile is ~60 s — do not read a
  wall-clock timeout as non-convergence; compile once per debugging process.)
- **Match the linear solver to the MATRIX TYPE, then VERIFY it converged.** Before
  choosing a KSP/PC, determine the operator's symmetry from the physics — a single
  field elliptic problem (Poisson, diffusion, linear elasticity) is **SPD** → `cg` +
  AMG (`gamg`/`hypre`); a **coupled multi-field tangent** (u-eₚ plasticity, u-φ/u-c-φ
  fracture, plastic tangents) is **NON-SYMMETRIC** → `gmres` (never `cg` — it stalls
  silently) + `pc="fieldsplit"` to isolate the blocks. **Plain AMG DIVERGES on a
  coupled non-symmetric matrix** (measured: gamg 166–193 KSP iters → cap; FieldSplit
  with a strong per-block sub-solve → 14). If unsure of symmetry, check
  `‖K−Kᵀ‖/‖K‖`. Then **verify**: a KSP at its iteration cap (`info["ksp_its"]`,
  `info["ksp_diverged"]`, or `-ksp_converged_reason`) means the **preconditioner is
  wrong for the operator — change the PC, do not raise `max_it`.**

## Generator-Friendly Python

Use:

- `abaqus_ufl.core.tensor` helpers: `eye`, `trace`, `det`, `inv`,
  `dev`, `eig`/`eigh`, `logm`, `expm`, `sqrtm`, `dyad`, `sym3`. For the
  symmetric part write `0.5*(A+A.T)` — the translator has no unary `sym`.
  `eigh` is an alias for the general `eig` (not a separate symmetric solver).
- 3x3 tensor form until the Abaqus boundary.
- bounded `for k in range(N)` loops with optional `break`.
- explicit scalar/tensor state variables in `state_vars`.
- explicit three-item compare/swap logic instead of `np.argsort`.
- smooth regularization where complex-step derivatives matter.

Avoid in generated methods:

- `while`, dynamic `append`, list comprehensions, dataclasses as generated
  state carriers, `np.einsum`, hidden file I/O, general dynamic arrays, and
  undocumented parameter guesses.

If generation fails, prefer rewriting the model into supported 3x3 tensor
operations. Add generator support only when the pattern will recur. Use
`@au.fortran_helper` sidecars only when the DSL path is not practical.

### Finite-Strain Tensor-Return Rule

For finite-strain UMAT/UEL methods returning a tensor, prefer explicit scalar
and tensor temporaries around matrix-to-scalar calls:

```python
J = det(F)
FinvT = inv(F).T
return G * (F - FinvT) + K * log(J) * FinvT
```

This is clearer than a one-line return such as
`G*(F-inv(F).T) + K*log(det(F))*inv(F).T`. The direct form is now regression-
tested, but it previously exposed a generator bug that emitted
`det33z(F(ii,jj))`, which compiled and failed only at runtime with NaNs. After
generating finite-strain Fortran, grep for impossible component-indexed helper
calls before blaming physics:

```bash
rg "det33z\\([^)]*\\(ii,jj\\)" path/to/generated.for
```

If any matrix helper receives `(...(ii,jj))`, treat it as a code-generation
shape bug, not a constitutive or solver issue.

## API Patterns

Finite-strain UMAT:

```python
import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log

class Model(au.Material):
    props = dict(...)

    def stress_PK1(self, F):
        ...

model = Model()
model.verify()
au.generate_umat(model, "model.for")  # default matrix_backend="iterative"
```

Small-strain UMAT:

```python
from abaqus_ufl.core.tensor import eye

class Model(au.SmallStrainMaterial):
    props = dict(...)
    state_vars = dict(alpha=0.0, eps_p=0.0 * eye(3))
    stress_convention = "compression_positive"

    def stress_update(self, sigma_old, strain_old, dstrain,
                      alpha_old, eps_p_old, dt):
        ...
        return sigma_new, {"alpha": alpha_new, "eps_p": eps_p_new}

model.verify()
au.generate_small_strain_umat(model, "model.for")
```

The values in `state_vars` are both the initial values and the type/layout
declarations. Scalars consume one `STATEV` slot; 3x3 tensors consume nine.

UEL:

```python
class Problem(au.WeakForm):
    material = Model
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)
        self.c = au.ScalarField("c", degree=1, test="zeta")

    def momentum_equation(self, v, F, c):
        return self.material.stress_PK1(F, c)

    def species_transport_equation(self, zeta, F, c, grad_c, c_old, dt):
        return storage, flux

problem.verify()
au.generate_uel(problem, "model_uel.for", element="quad4",
                formulation="standard")
```

Scalar UEL sign rule: `return storage, flux` means the generated weak residual
is `storage * test - flux . grad(test)`. For ordinary diffusion with
`storage = c_dot`, return the physical flux `-D * grad_c`. For AT2
damage/fracture with
`storage = Gc/ell*d - 2*(1-d)*H`, return `-Gc*ell*grad_d` so the weak form has
`+Gc*ell*grad(d).grad(test)`. For phase-field models, do not reason from
diffusion flux examples. Write the positive weak-form gradient term first, then
return the adapter required by the generator convention. Do not copy a flux sign
from another phase-field example by field name; derive it from the exact
`storage` quantity returned by that equation.

For branch-dependent gradient damage, the flux method must receive the same
state needed to choose the branch as the storage method. If `phase_storage`
uses `psi_star_ductile` in compression and `phase_flux` hard-codes the brittle
`psi_star`, CS-vs-FD still passes because both sides differentiate the same
wrong residual. Add a material-point gate that measures the effective
regularization length per branch and fails if it differs from the declared
`ell`.

Major symmetry: `verify()` checks dP/dF symmetry for stateless materials by
default. A material whose stress is DELIBERATELY not an energy derivative —
e.g. VMS-stabilized mixed formulations like the Scovazzi u-theta element,
where `P = F S(E(u, theta))` approximates the momentum weak form and the
theta equation carries the missing coupling — declares
`symmetric_tangent = False` on the Material class to skip the check (with a
code comment citing why). Do NOT use this flag to silence a symmetry failure
in a material that should be hyperelastic: there it is a bug signal, which
is why the default stays strict.

For UELs with more than one scalar old field, such as `phi_old` and `c_old`,
inspect the generated `*_cs_tangents` subroutine and its call site. The
history arguments must use the same field-declaration order in both places.
A residual-only smoke can pass even when the tangent receives swapped old
values.

Use the standard UEL path for coupled fields. Use F-bar/local-pressure
formulations only when the element theory justifies them; they are not default
stabilization switches for gels, diffusion, phase fields, or transport.

Local-pressure UELs:

- Supported as a prototype `formulation="local_pressure"` path for Quad4 and
  Hex8 `u,mu` gel-style elements with one condensed element-local pressure.
- `SVARS(1)` is the solver pressure state.
- If `Variables >= 1 + 2*NGP`, diagnostics are mirrored as
  `SVARS(1+i)=phi_i` and `SVARS(1+NGP+i)=p_i`.
- Abaqus gets the same diagnostics through generated `UVARM1=phi` and
  `UVARM2=p`.

## When To Read References

- For choosing a prior example: read `references/example-map.md`.
- For test and validation expectations: read `references/validation-ladder.md`.
- For known traps and debugging: read `references/pitfalls.md`.

These references are short summaries. For authoritative detail, open the repo
docs linked in each reference.

## Done Definition

A new example is not done until it has:

- scoped equations and non-scope documented;
- properties, units, sign convention, and `STATEV`/`SVARS` layout documented;
- Python reference checks or material-point tests;
- explicit regime-entry checks for each claimed branch;
- generated Fortran compile coverage;
- f2py or Abaqus validation appropriate to the model;
- **at least one INDEPENDENT, QUANTITATIVE oracle** — a digitized
  paper-figure point, a closed-form limit, or a hand computation — that
  is NOT derived from the code under test (see below);
- clear status labels: "code implementation complete", "solver stabilization
  open", or "paper reproduction complete".

Do not collapse these statuses. A model can be code-complete without being a
paper-figure reproduction.

**Why the independent oracle is non-negotiable.** Four hard models this
project shipped with green suites all had real bugs (wrong signs, wrong
exponents, wrong push-forward `Fp`, 100x-off parameters). The suites
could not see them because they were built from: (1) consistency checks
(CS-vs-FD / `verify()` / compile differentiate the *wrong* code
faithfully); (2) code-vs-itself oracles (f2py comparing generated
Fortran to the Python reference — both share the bug, agree to ~0, both
wrong vs the paper); (3) qualitative criteria ("damage grows", "peak ~
X") that pass through quantitative errors; and (4) a benign loading path
(axisymmetric / elastic-only / `grad=0` / `Fp=I`) that never exercises
the suspect term. A green suite with these four properties certifies
almost nothing about constitutive correctness. The one gate that breaks
all four is a quantitative comparison to something the code did not
produce. Add it.

## Independent Review Pattern

For hard research models, use independent AI review as a validation-envelope
attack, not just a style review. The method that found every bug: read the
paper equation, re-implement that one term independently, evaluate the code at
a state that actually exercises the term, and compare numerically. Do not trust
the green suite, your own paper reading, or a sub-helper's transcription until
the numbers agree (a helper mis-stated a matrix factor order this session; the
numerical check caught it). A reviewer should ask:

- Which branch did `model.verify()` actually exercise? (A default state at
  `e=e_max`, `grad=0`, or `Fp=I` can make the whole nonlinear core elastic /
  inert and uncertified.)
- Did the tests reach plastic flow, damage, residual strain, fabric rotation,
  split yield intervals, non-homogeneous gradients, near-singular limits, or
  NON-COAXIAL loading? (Axisymmetric/diagonal loading masks every
  commute-dependent algebra bug.)
- Are readable Python and generator-facing paths mathematically isomorphic, or
  did a codegen rewrite change data structures? Do they AGREE numerically in
  the regime where they could diverge?
- Is there an independent quantitative oracle, or only consistency / code-vs-
  itself / qualitative checks?
- Does the status claim match the hardest path that has actually passed?
- Every regression should ship a deliberately-broken control: reintroduce the
  bug and confirm the test fails. A formula asserted against itself tests
  nothing.

The recurring generation failure modes to hunt (see `references/pitfalls.md`):
pattern-matched term vs derived
dual; coaxial-only matrix algebra; constraint-used-as-definition;
parameter-table cross-contamination; one symbol reused for two roles; a second
code path drifting from the validated reference; finite-strain tensor-return
mis-tensorization when `det(F)` is nested inside a direct tensor return instead
of assigned to a temporary.

A multi-agent implement/review/fix loop is the reference pattern: one agent
implements, an independent reviewer finds untested codegen defects, a second
agent fixes and adds adversarial regressions, and a re-review localizes any
remaining mismatch instead of hiding it. Keep this loop explicit for any model
that is more than a clean Level-1 UMAT.
