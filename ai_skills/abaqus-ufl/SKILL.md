---
name: abaqus-ufl
description: Use when developing, reviewing, verifying, or debugging abaqus_ufl UMAT/UEL examples, generated Abaqus Fortran, f2py checks, Abaqus integration decks, or paper-to-code model implementations. Guides agents through target selection, generator-friendly Python, the example pipeline, and known Abaqus/Fortran pitfalls.
metadata:
  short-description: Develop and verify abaqus_ufl UMAT/UEL models
---

# abaqus_ufl Model Development

This skill is a compact operational guide distilled from the repository's
tests, examples, and lessons. Use it before implementing or reviewing any
`abaqus_ufl` UMAT, UEL, f2py, or Abaqus-integration work.

## First Moves

1. Check the worktree:
   ```bash
   git status --short
   ```
   Do not overwrite collaborator edits, generated files, or Abaqus-integration
   artifacts without understanding them.

2. Read the live API guide in the repo:
   - `docs/API_USAGE.md`
   - `HOWTO_ADD_AN_EXAMPLE.md`

3. Pick the target deliberately:
   - **finite-strain UMAT**: local constitutive law, built-in Abaqus elements
     are enough, model returns `stress_PK1(F)`.
   - **small-strain UMAT**: incremental material update with local `STATEV`;
     use `SmallStrainMaterial.stress_update(...)`.
   - **UEL**: any extra solved nodal field, gradient term, mixed interpolation,
     phase field, pore pressure, diffusion, temperature, concentration, or
     chemical potential.
   - **f2py**: point/element check for generated Fortran before solver runs.
   - **Abaqus integration check**: production-solver, convention, and
     output-bridge evidence; not automatically an independent physics oracle.

## Core Rule

Do not start from Fortran. Start from the theory, write generator-friendly
Python, verify the Python model, generate Fortran, compile, then execute the
example pipeline.

The existing example pipeline is:

1. Document theory, scope, conventions, fields/DOFs, and state layout.
2. Add Python reference checks and an independent quantitative oracle.
3. Run `model.verify()` or `problem.verify()` for implemented-method tangent
   consistency, plus regime checks for every claimed branch.
4. For a UEL, check assembled `RHS`/`AMATRX`, DOF/state layout, and an
   appropriate patch or invariant; `problem.verify()` alone is not this gate.
5. Generate deterministically and compile every generated Fortran source.
6. Call nontrivial generated code through f2py or another checked compiled
   runtime.
7. Add a small solver run only when it provides useful evidence.
8. Validate the output bridge before using solver results in comparisons or
   plots.
9. Attempt a published-figure reproduction only after the smaller rungs are
   stable.

## Verify, Do Not Guess

When a build fails, a solve diverges, or you suspect a framework limitation,
**test the hypothesis before acting on it**, and check the theory, the API, and
the existing examples first. Most "the framework can't do X" or "it must be Y"
beliefs are wrong:

- **Check available evidence before claiming a limitation.** Bounded local
  Newton/implicit updates inside a generated method are supported when written
  with generator-safe operations and explicit convergence/finite checks.
  Inspect the public manifest, present example folders, package source, and
  linked tests or evidence before concluding something is impossible. Do not
  infer public capabilities from an unpublished model.
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
- **Keep solver setup separate from package correctness.** Determine whether a
  tangent is symmetric and declare the Abaqus interface consistently, but do
  not make a user- or machine-specific solver configuration part of the
  package contract.

## Generator-Friendly Python

Use:

- `abaqus_ufl.core.tensor` helpers: `eye`, `trace`, `det`, `inv`,
  `dev`, `eig`, `logm`, `expm`, `sqrtm`, `dyad`, `sym3`. For a matrix's symmetric
  part write `0.5*(A+A.T)`: the translator has no unary `sym` (it exists in Python
  for the oracle but does not translate to Fortran). `eigh` works inside a
  translated method as an alias for `eig`, but is not importable in Python.
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

For a branch-dependent gradient model, the flux method must receive the same
state needed to choose the branch as the storage/reaction method. CS-vs-FD can
still pass when both sides differentiate the same wrong branch. Add a
branch-specific invariant that checks the declared reaction/gradient
relationship.

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
- assembled residual/tangent plus patch/invariant checks for a UEL;
- generated Fortran compile coverage;
- direct compiled execution appropriate to the model;
- output-bridge mappings, identity, coverage, conventions, and status for every
  claimed solver result;
- **at least one INDEPENDENT, QUANTITATIVE oracle** — a digitized
  paper-figure point, a closed-form limit, or a hand computation — that
  is NOT derived from the code under test (see below);
- clear status labels: "code implementation complete", "solver stabilization
  open", or "paper reproduction complete".

Do not collapse these statuses. A model can be code-complete without being a
paper-figure reproduction.

**Why the independent oracle is non-negotiable.** Reviews of difficult
development examples have repeatedly found real equation, convention,
state-update, and parameter-transcription bugs behind green suites. The suites
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
attack, not just a style review. A repeatedly effective method is to read the
paper equation, re-implement that one term independently, evaluate the code at
a state that actually exercises the term, and compare numerically. Do not trust
a green suite or a transcription until the numbers agree. A reviewer should
ask:

- Which branch did `model.verify()` actually exercise? (A default state at
  a homogeneous, reference, or identity configuration can leave the nonlinear
  core inert and uncertified.)
- Did the tests reach every claimed branch and state evolution, plus
  non-homogeneous, near-singular, and non-coaxial states where applicable?
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
that is more than a minimal clean UMAT.

## Public Capability Boundary

The public AI example map is the positive release manifest, not a history of
everything developed in the lab. Reference only example folders that are
actually shipped and listed in `examples/README.md`. Do not disclose,
reconstruct, or advertise internal, unpublished, license-unclear, or deferred
models as v1 capabilities.
