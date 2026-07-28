# Field Notes: General Lessons

Reusable engineering lessons for anyone building finite-element solvers,
constitutive code generators, and coupled-multiphysics runtimes. Distilled from
a development log; the project-specific chronology is stripped out. Each note is
meant to transfer to a new model, element, or solver.

## Verification: what a green suite actually certifies

- **Consistency is not correctness.** Complex-step-vs-finite-difference tangent
  checks, compile gates, and smoke runs certify that the code computes *what it
  says* — they differentiate the (possibly wrong) code faithfully. A sign flip,
  a wrong exponent, or a wrong parameter passes all of them. Worse are
  code-vs-itself oracles: a generated kernel checked against the reference it was
  generated from shares every bug, so the two agree to round-off while both
  deviate from the truth. Every model needs at least one *independent,
  quantitative* oracle — a digitized figure point, a closed-form limit, or a
  hand computation not derived from the code under test. And every regression
  ships a *broken control*: reintroduce the bug and assert the test fails. A test
  that asserts a formula against itself — or one generated copy against another —
  cannot fail and tests nothing; it must assert an *independent* property
  (definiteness, a dispersion sign, a manufactured-solution residual).

- **Verify the regime, not one state.** A single benign verification state
  (smooth, coaxial, homogeneous, gradient-free)
  hides whole branches, non-coaxial terms,
  the anti-diffusive term that vanishes at homogeneous states. "All tests pass"
  is a claim about the *validated envelope*, not the model — so state the
  envelope explicitly and drive coverage states that actually enter the claimed
  regime. A material-point gate harness that sweeps a `step(F, state, dt)` API
  across regimes and asserts (coverage) the branch was entered, (convention) the
  measured regularization length equals the declared one, and (robustness)
  finite stress/state/tangent across a softening path, catches what a single
  state cannot.

- **Operator-level and invariant gates see structural facts consistency checks
  cannot.** Block definiteness for storage/flux sign pairings, coupled
  dispersion for mixed pairs, kinematic consistency, major symmetry, and
  per-field block-scale balance are all invisible to a differentiate-the-code
  check but decide whether the physics is right. When you can, *assemble the
  block and test it directly* rather than inferring from a heuristic — assembly
  self-calibrates and covers sibling sub-cases (a screen tuned for
  reaction-diffusion can be silently wrong for pure diffusion). Code-generation
  string bugs are equally invisible to any weak-form check: gate them at the
  generated-code level (assert the bad emission pattern is absent) plus a
  finite-tangent run.

- **Know the characteristic failure modes of fast constitutive ports.** Green,
  carefully-written code fails in predictable ways: pattern-matching a term
  against a convention or a sibling instead of deriving it; matrix algebra that
  is exact when everything commutes (diagonal, axisymmetric) and wrong off-axis;
  conflating a constraint (`≥ 0`, "must hold") with a definition; parameter
  cross-contamination between sibling tables/mechanisms; reusing one symbol for
  two physical roles; and a second code path drifting from the validated
  reference. When a quantity and its variational dual both appear (stress ↔
  flux from the same energy, residual ↔ tangent), derive *both* from the weak
  form — getting one right is evidence about diligence, not about the other.

- **Check the already-validated cases before treating a discrepancy as a bug.**
  Before escalating, establish what is *already* validated: the validation
  ladder, the oracle it was checked against, and the exact quantity. Distinguish
  "the model is validated" (constitutive, equilibrium, against an oracle) from
  "this new comparison matches one external run" — the latter failing is not a
  model-validation failure and is often lower-stakes. Reading the validated
  setup first surfaces cheap mismatches (a temperature offset, a stale
  parameter) you would otherwise chase for days; going to the authoritative
  source (the paper) both confirms the validated quantity and characterizes the
  open one in minutes.

- **A pattern is not an oracle.** Some outcomes are genuinely non-unique:
  a degenerate bifurcation (e.g. conjugate shear bands at ±45°) admits both a
  crossed pattern and a single band as legitimate solutions of the same problem,
  selected by perturbation source and boundary kinematics. Never pass or fail a
  reproduction on an image (band orientation/multiplicity) alone; gate on
  pattern-insensitive observables (force-displacement curve, drop location, band
  angle family).

## Solvers and performance

- **Diagnose the bottleneck before rewriting; name the real lever.** Profile
  first — measurement routinely refutes assumptions. A framework switch is not
  an algorithm switch: the same direct solver behind a new interface buys
  nothing; the value is access to a *different algorithm* (iterative + multigrid
  + field-split). Language is not the lever when the cost is in compiled kernels
  (a hot assembly at 5% and a compiled solve at 92% means a rewrite is justified
  for driver flexibility and algorithm access, not because the new language is
  faster). And measure, don't recite the textbook: a symmetric factorization can
  be slower than unsymmetric LU, and adding parallel ranks to a tiny per-rank
  workload can degrade.

- **Fixing a bottleneck moves it — re-profile to find the new one.** Delegating
  the linear solve to optimal multigrid can move the cost to the assembly loop;
  removing per-element marshalling can expose the sparse-scatter cost. Also
  profile the *representative* workload state, not a trivial one — an element
  that does almost nothing at `U = 0` is dominated by call overhead there and
  compute-bound in the real plastic regime, so a speedup projected from the
  cheap state evaporates. And profile before optimizing: the budget is often in
  matrix functions (logm/expm/polar/eig repeated for the tangent), not the local
  Newton iteration count you were about to cut. And verifying *correctness* is
  not testing *performance*: manufactured-solution rates and benchmark matches
  say nothing about scaling, so run the performance benchmark explicitly — it is
  what closes the loop on the *premise* of a rewrite.

- **Field-wise convergence is the whole game for coupled problems.** A coupled
  residual that mixes a large block (momentum ~ modulus) with a tiny-coefficient
  block (transport) is dominated by the large field in a single global `‖R‖`, so
  a global gate stops while the weak field is still iterations short. A small
  residual on a near-singular block is a *large solution error* (`err ≈
  R/coeff`), so gate each field on its own residual relative to its own
  characteristic scale, and balance the line search the same way. Diagnostic
  signature: an error that *compounds over time* (early match, late drift)
  points at convergence gating, never at the physics or BCs (those diverge from
  t=0); a healthy 3-5 Newton iters/step means the tangent is fine and the gate
  is the bug, while 15-40 iters implicates the tangent. Field-wise convergence
  fixes the *gate* but not the matrix *conditioning*: huge per-field coefficient
  ratios (SI-unit artifacts) still wreck any iterative solve, so
  non-dimensionalize the weak form (length, potential, time, stress) — a
  two-sided row *and* column scaling, not a magic constant in the driver.

- **Field-split choice follows coupling strength; diagnose with a direct block
  sub-solve.** Set each block's sub-preconditioner to an exact solve: if the
  coupled system then converges in O(1) iterations, the block structure and
  coupling are excellent and the difficulty is a *single block's*
  preconditioner. Strong two-way coupling wants a Schur split; weak or
  semi-implicit coupling is fine with additive/multiplicative. A degraded
  stiffness block whose coefficient spans many orders (e.g. `g(d)·C` with `g →
  0`) defeats AMG coarsening — that within-block *contrast*, not inter-block
  magnitude, is the obstacle, so equation/diagonal scaling does not help and can
  hurt smoothed aggregation; a robust AMG (hypre) or a block-direct solve is the
  reliable route. Beware: a weak field-split preconditioner combined with
  inexact-Newton forcing can report *false convergence* — pair forcing only with
  a strong split or an exact block solve. And before claiming "monolithic stalls,
  staggered rescues," check whether the model's time discretization already
  decoupled the stiff block (a semi-implicit driving force frozen at the step
  start removes the worst coupling indefiniteness, half-staging the monolithic);
  when a coupled solve walls out, check *where* the residual lives — if it is
  entirely in one field's own kernel, the fix is in that constitutive law, not a
  fancier coupling scheme.

## Distributed and MPI

- **More ranks running slower is usually thread oversubscription, not "the
  solver doesn't scale."** With `OMP_NUM_THREADS` unset, every MPI rank's
  threaded BLAS/direct-solver spawns threads across all cores — N ranks × all
  cores thrashes, worse with each rank. Pin one thread per rank
  (`OMP_NUM_THREADS=1`, and the BLAS-specific variants). A monotonic ~N×
  slowdown is an *artifact* signature; real "direct doesn't strong-scale"
  plateaus or degrades gently. Tell-tale: assembly scales while the solve
  anti-scales.

- **Container/cloud MPI has two traps that pass at 2 ranks and break above.**
  Prefer a fork-based launcher (MPICH/Hydra) over a daemon-mesh launcher
  (OpenMPI/PRRTE) for single-node multi-rank — the daemon mesh probes network
  interfaces and hangs when virtual/bridge interfaces are present. Keep your MPI
  Python bindings ABI-consistent (mixing builds linked against different MPI
  corrupts collectives above 2 ranks); if in doubt, use one binding's MPI and
  don't import the other. Always validate distributed correctness *at >2 ranks*
  (identical result across 1/2/4/8/16) — a 2-rank match exercises none of the
  reduction-tree, off-rank-scatter, or preconditioner-blocking general cases.
  And clean up orphaned launcher daemons between runs; they poison the next
  launch.

- **A Fortran-order array fed to a row-major consumer silently transposes your
  matrix.** An `f2py`/column-major element tangent handed to a row-major API
  (most C libraries) assembles Kᵀ. Element-level value parity does *not* catch
  this — the values are identical, only the storage order differs. A symmetric
  operator hides it entirely, so it is a latent trap that first bites on a
  non-symmetric (e.g. plasticity) tangent, showing up as "correct serially,
  wrong in the path that hands the raw buffer across the boundary." Enforce
  C-contiguity at the boundary.

- **Bounding a field in a staggered driver: commit-time projection beats
  in-sweep masking.** Clamping after every field update and masking the clamped
  DOFs out of the convergence residual is correct on paper and wrong in
  practice: the active set cycles between sweeps and the loop grinds, and a
  rank-local guard around collective calls is a latent parallel deadlock.
  Instead let each step converge exactly as the proven unconstrained scheme,
  then clip the bounded field's owned DOFs once before the state commit — a pure
  local op, rank-independent by construction. Diagnostic: all rank counts
  stalling at the *same* step is a deterministic algorithmic grind; stalls at
  *different* steps suggest a parallel race.

- **A scaling result is a claim about one problem class, not the code.**
  Separate *architectural* scaling (the pattern has no serial bottleneck —
  transferable) from *achieved* scaling (iteration count, memory, wall time at N
  DOFs — problem-class-specific). An AMG iteration count flat for scalar Poisson
  can grow for an indefinite coupled operator; "validated to N on Poisson" does
  not license an N-on-the-coupled-problem number. Quote scaling only for the
  class you measured.

## Code generation and elements

- **Read the paper PDF, not an AI summary, before writing equations.** Summaries
  confidently report equations that are nowhere in the paper; implementing from
  one produces a confidently-wrong model that still passes consistency checks.
  Use search to find and orient, then read the saved PDF for the verbatim
  numbered equations, and cite the exact equation number at the code site.

- **Text-level code-generation passes produce code that compiles but fails at
  runtime; guard the emission.** Two recurring traps. (1) *Case-insensitive
  namespace*: a scalar field `d`, a state variable `d_old`, and a local `D`
  collide, and `_` is not a valid name — use descriptive field names (`damage`,
  `phase`, `concentration`), reserve short physics symbols for comments, and make
  any post-AST text scan reason about the *base name* of an indexed left-hand
  side (strip the `(...)` suffix) so it doesn't re-declare indexed locations as
  new variables. (2) *Expression shape*: a component-wise tensor-return path can
  index a tensor inside a scalar function call (`det(F)` → `det33z(F(ii,jj))`) —
  valid Fortran with no explicit interface, so it compiles but returns NaN at
  runtime. Assign scalar and tensor subexpressions to temporaries before the
  tensor return (`J = det(F)`, `FinvT = inv(F).T`); and emit `dx*dx`, not `z**2`
  (which lowers to complex `exp(2·log z)` and is NaN at `z=0`), for squares of
  possibly-zero quantities.

- **Branchy, rate-independent laws can be made complex-step-safe.** The fixed
  recipe: smooth viscoplastic flow instead of a return-map branch, smooth
  Macaulay `⟨x⟩ = ½(x + √(x²+reg²))` for every bracket, `√(·+reg)` for every
  norm, `½(F+Fᵀ)` for `sym`, an `exp`-built `tanh`. A *local Newton inside the
  element* is also CS-safe: branch on `.real`, a fixed-count loop, an analytic
  Jacobian — the complex-step perturbation rides through the converged iterate
  and yields the consistent tangent for free (an `if f > 0` elastic/plastic
  branch is CS-invariant because the perturbation is purely imaginary).

- **Adding an element type is template + config only.** A new element goes in as
  a shape-function template plus config entries, reusing dimension-generic
  Jacobian helpers, with no generator-code changes. The decisive gate is a
  *distorted-mesh patch test driven through the compiled element* (affine
  Dirichlet everywhere; interior fields exact to round-off) — it catches
  node-ordering, quadrature-weight, and mapping bugs in one shot.

- **Every runtime adapter must derive its layout from the element config, and
  seed/preserve state.** State-variable count, DOF map, coordinate dimension,
  and Gauss count must be read from the element's active layout, never
  hardcoded to one topology's default — a generated element that passes in one
  runtime does not certify a *different* adapter, which infers its own layout.
  Adapters that zero-fill state must seed nontrivial initial conditions (`Fp=I`,
  `S=S₀`, condensed pressure) or the first call hits `inv(0)`, and must persist
  field-level condensed state across steps that the host solver preserved
  automatically.

- **Convert a sign convention once, at the seam, and finite-difference the
  assembled tangent.** Elements that return `AMATRX = −dRHS/dU` need one adapter
  that flips to `K = dR/dU`; scatter the conversion through the solver and Newton
  silently goes uphill (an FD check gives `fd == −K`, which merely looks like
  ill-conditioning). Always one-column-FD an assembled tangent before trusting a
  hand-rolled solver — and for a tangent taken at a *fixed previous state*,
  perturb the total and increment together, or the check gives a spurious
  mismatch on correct code.

- **Know your generator's defaults, and don't flip a flag to mask a bug.** A
  displacement-only element may silently auto-select a stabilized formulation
  (e.g. F-bar) — right for production, wrong as a *locking baseline*; check the
  summary output. A locking-control formulation is not a generic "make low-order
  better" switch you can force onto a coupled weak form without a real
  derivation. And a symmetry opt-out (`symmetric_tangent = False`) is legitimate
  only where the formulation is asymmetric by design (stabilized methods where a
  companion equation carries the missing coupling) — never reach for it to
  silence asymmetry in a material that *should* be hyperelastic; there it is a
  bug signal.

## Constitutive and numerical formulation

- **Solve in a smooth, well-scaled local unknown.** A Newton on the obvious
  increment can be singular or extremely stiff even when the physical state is
  benign. Reformulate the consistency equation using a standard
  computational-plasticity unknown with a finite derivative. This is a
  constitutive-integration choice, not a reason to add arbitrary solver
  damping.

- **Integrate internal-variable ODEs with their exact regime solutions — the
  closed form becomes a free oracle.** Forward-Euler on a stiff hardening/
  softening ODE needs tiny steps; using the exact regime integral (a saturating
  hardening integral, exponential relaxation) is unconditionally stable and exact in each
  pure regime. Bonus: the material-point test can then check the accumulated
  state against the paper's integrated form to relative round-off — an
  independent oracle that costs nothing. When an exponent makes units
  ambiguous, cross-check a hardening modulus against a *digitized figure*, not
  just the table.

- **An explicit trial-based return map fails at the resolution a sharp feature
  needs; the fix is constitutive, not a solver swap.** A smooth-explicit slip
  update chosen for CS-safety overshoots at a resolved stress concentration and
  diverges once the mesh is fine enough to see it — the linear solver is a red
  herring (a robust direct solve merely gives the apples-to-apples comparison
  that *locates* the element as the culprit). Swap to an implicit return map.
  Separately, guard convergence tests against NaN (`NaN ≥ tol` is `False`, so a
  diverged step reports "converged"), and never cut `dt` at a
  storage-stiffness stall — backward-Euler storage scales as `1/dt`, so cutting
  `dt` makes the block *stiffer* and the residual floor *rise*.

## Reproducing papers and debugging discipline

- **Benchmark reproduction includes the loading model.** A result systematically
  offset by a *uniform* amount across element types is the fingerprint of a
  BC/load-model difference, not an element bug — suspect it first. Classify
  every load as dead (reference-frame) vs follower (pressure, deformed-face
  cof(F) map) before running; at large end rotations they differ substantially.
  The small-load linear limit against an analytic solution (Euler-Bernoulli,
  MMS) is the cheapest possible oracle and belongs as a standard rung.

- **Tuned-for-convergence is not reproduced.** A softening solve that localizes
  is a *solver* result; a *figure* reproduction needs the paper's exact
  constants, geometry, BCs, and load rate. Relatedly, "it stalls at X" needs a
  run that actually *reaches* X — a converged run that stops short of a regime
  says nothing about that regime — and "it converged" is not "it reproduced the
  physics" (a coarse run can converge to a stress-concentration artifact at the
  BC corners). Match the observable to the claim.

- **Localize before concluding — symptoms mislead.** When consistency checks and
  compile pass but the solver NaNs, the bug is an input the test states don't
  cover — most often the *deck* (positional property arrays silently misalign
  when a model's parameter count changes, worst when the misaligned slot is a
  denominator). Isolate it by probing the generated subroutine directly from a
  tiny driver: probe-finite-vs-solver-NaN at the *same* state points straight at
  the one argument the probe didn't match, in minutes instead of spelunking the
  generator. Apply the same discipline elsewhere: grep the notes and code before
  asserting "we can't do X because we lack Y" (the capability, or the reason it
  isn't needed, is often already written down); check the ranks × threads vs
  cores confound before blaming an algorithm; and treat a consistency-check
  failure at an *extreme* state as suspect FD, not a tangent bug — sweep the FD
  `eps` first, and if FD approaches the analytic value from above as `eps` grows,
  the analytic tangent is right.

- **Operational hygiene for long runs.** Run long jobs unbuffered with flushed
  per-step progress (stdout to a file is block-buffered; an empty log is not
  "hung"). Once a fast solver is validated, use it for production — if
  convenience features live only in a slow path, add them to the fast path
  rather than run the slow one. And remember a recompile can eat a wall-clock
  budget: don't read a timeout as non-convergence when the process may still be
  building. The first real end-to-end run *is* the test — a one-line fix in a
  code path nothing exercises is feature work, so wire new machinery into a real
  run to flush out dead-code rot.

- **Verify a generated mesh geometrically, with an oracle.** A mesh that imports
  with positive Jacobians can still be the wrong shape. Assert element purity
  (any stray triangle among quads is a hard failure, not cosmetic), reorient
  every element to a positive Jacobian independently (shoelace signed area, then
  an independent per-Gauss-point `detJ > 0` check — not the assembly's own map,
  which is circular), and compare total mesh area to the analytic CAD area. Cut
  geometry (curved arcs) can otherwise come back as a staircase.

## Reuse and architecture

- **Reuse the kernel; don't reimplement it.** When a new runtime drives the same
  generated element as the reference/production path, the physics bug surface
  does not grow — only the driver and linear algebra are new. Keep the slow
  reference path precisely as the *oracle* that cross-validates the fast path;
  don't delete it once the fast path works. When weighing rewrite vs reuse,
  count how much *correct, tested* machinery you can keep — reusing the elements
  and the linear-algebra library and writing only the driver is what makes a
  solver tractable in a handful of steps.

- **Get the state invariant right, and confirm which of two homes a value has.**
  For internal state variables: evaluate every element against a *copy* of the
  committed old state during Newton iterations, commit new→old only after the
  step converges, and recover the new state by *reusing the assembler's own
  Gauss-point loop* — a separate hand-written recover can subtly diverge from
  the residual evaluation. And whenever a value can live in two places (a
  constructor-set instance attribute vs a class-level defaults dict), confirm
  which one the runtime path actually reads — and that your *other* path reads
  the same one; the symptom of a mismatch is agreement at the state where the
  value barely matters and blow-up where it does.
