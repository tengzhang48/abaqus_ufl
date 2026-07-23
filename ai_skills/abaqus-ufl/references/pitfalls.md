# abaqus_ufl Pitfalls

This file lists failure modes repeatedly found in tests, examples, and notes.
Use it when debugging generated UMAT/UEL behavior.

## Abaqus Interface

- `RHS` is two-dimensional: use `RHS(i,1)`, not `RHS(i)`.
- UEL tangent convention is `AMATRX = -d(RHS)/d(U)`.
- Abaqus does not pass old UEL field values. Reconstruct old values as
  `old = U - DU`.
- `DTIME = 0` can occur during initial stiffness calls. Guard divisions by
  `DTIME`.
- Explicitly zero `RHS` and `AMATRX`; do not trust incoming memory.
- Use `PNEWDT` for inverted elements, large increments, or detected invalid
  states. Set `PNEWDT` to a fraction such as `0.25` to request a cutback and
  to `1.0` when the state is acceptable. Larger growth requests are
  solver-dependent and should be used deliberately.
- Use `UNSYMM` for coupled UELs and genuinely nonsymmetric UMAT tangents unless
  symmetry has been proven. If the tangent is unsymmetric and `UNSYMM` is
  omitted, Abaqus may use the wrong tangent treatment and convergence can
  degrade or fail. If `UNSYMM` is declared for a symmetric tangent, Abaqus pays
  the cost of an unsymmetric solver. `UNSYMM` does not fix bad physics or large
  steps.

## DOF and State Layout

- Generated UELs are node-major:
  ```text
  node 1: ux, uy, scalar_1, scalar_2, ...
  node 2: ux, uy, scalar_1, scalar_2, ...
  ```
- Abaqus `*User Element` active DOF labels must match generated ordering.
- Mixed-degree fields must define which nodes own which scalar DOFs.
- Tensor UMAT `STATEV` values consume 9 slots each.
- Document every `STATEV`/`SVARS` slot before writing an Abaqus deck.
- For generated local-pressure UELs, do not overwrite `SVARS(1)`: it is the
  condensed pressure state. Diagnostics require extra slots:
  `SVARS(1+i)=phi_i` and `SVARS(1+NGP+i)=p_i`.
- Multi-history UEL tangents need explicit caller/callee order checks. A
  2026 audit found generated CS tangent subroutines declaring old scalar field
  histories alphabetically while the caller passed field-declaration order.
  Example bad pattern: callee `F_old, c_old, phi_old` and caller
  `F_old_gp, phi_old_gp, c_old_gp`. The residual can still be correct, so this
  shows up as Newton robustness/cutback trouble rather than a wrong-looking
  one-step result. Add or run a regression where old scalar fields have
  different values.

## Material / WeakForm Authoring

- Equation and material method names are **matched, not validated**. Only the
  canonical names in `core/defs.py:METHOD_INFO` and
  `core/weakform.py:EQUATION_INFO` are detected (`stress_PK1`, `pressure_resid`,
  `solvent_storage`, `solvent_flux`, `momentum_equation`, `pressure_equation`,
  `transport_equation`, ...). A method with a different name (e.g. a paper-spec
  sketch's `pp_potential`, `constraint_resid`) is **silently ignored** — no
  tangent block, no Fortran. When porting a spec, reconcile method names against
  `defs.py`/`weakform.py` first.
- Test-function-to-field mapping is by convention (`v->u`, `q->p`, `w->mu`;
  `weakform.py:TEST_TO_FIELD`). Map physics onto these slots: the value-assembly
  `pressure_equation` is the natural home for an algebraic constraint; the mixed
  `transport_equation` for a conserved field's storage+flux.
- `verify()`'s stress-free check is hardwired to a **zero-field** reference
  (`F=I, p=0, mu=0`; `core/verify.py:_check_stress_free`). A material that is
  stress-free at a **nonzero** natural field value (e.g. a volume fraction
  `Phi=Phi_ref`) fails this check spuriously while every tangent block passes.
  Drive the stress-free check at the true reference yourself; keep the
  framework's CS-vs-FD tangent and major-symmetry checks.
- Anything referenced as `self.<x>` inside a generated method must be a declared
  `prop` (it becomes `PROPS(i)`). Derived attributes computed in `__init__` do
  not translate.
- A `self._helper()` is translated to its own Fortran subroutine (UMAT and UEL
  paths). The helper body must use translator-supported idioms: branch with
  **statement-form** `if x.real > 0.0:` (maps to `IF (DBLE(x) .GT. 0)`), not a
  `hasattr`/ternary; Python floats already carry `.real`. A shared helper is
  emitted once. (Before mid-2026 the UEL helper path was broken; if you see a
  `field_args` NameError or "Helper method dispatch requires a material
  instance", the generator is out of date.)
- **State-variable update conventions (UEL `stress_PK1`)** — three rules, each a
  real compile/codegen failure in a coupled u-c-phi element:
  - Return the new state as a **dict literal in the `return`**:
    `return P, {'ep': ep_new, 'Fp': Fp_new}`. Assigning it first
    (`state = {...}; return P, state`) makes the translator hit a bare `Dict`
    statement -> "Unsupported expression: Dict".
  - `stress_PK1` must list **every** declared state var's `_old` in its Python
    signature — even one it recomputes and never reads (e.g. `hydro_old` for a
    `sigma_h` state). The codegen adds all state `_old` args to the subroutine; a
    missing one is passed REAL into a COMPLEX arg -> "Type mismatch ... passed
    REAL(8) to COMPLEX(8)" at compile, not at generate.
  - A **tensor** state var (`eps_p`, `Fp`) is typed as a tensor only inside
    `stress_PK1`; in any other method it is mis-inferred as **scalar**. A
    `self._helper()` called from both `stress_PK1` and a scalar-equation method
    with that tensor arg fails: "called with multiple argument kind signatures".
    Work around by computing the needed scalar (e.g. `sigma_h`) in `stress_PK1`
    and passing it to the scalar equation as a **scalar state var**. This also
    keeps the J2 return map out of the scalar equation.

## Fortran Generation

- Fixed-format column 72 still matters. Use the repo line wrapper.
- No `BLOCK`, variable-length arrays, or modern Fortran features in Abaqus
  fixed-form output.
- Use standard complex intrinsics: `LOG`, `EXP`, `SQRT`, not GNU-specific names.
- Fortran is case-insensitive. Avoid single-letter uppercase variables that
  collide with loop indices — **and with field/param names in coupled models**:
  a local `D` (a driving force) collides with a damage **field** `d` ->
  "Symbol 'd' already has basic type of COMPLEX". Rename shadowing locals
  (`Dval`; `triax` for a triaxiality `T` that would clash with `TIME`).
- A throwaway `_` in tuple unpacking (`a, _, c = f()`) becomes a Fortran variable
  named `_`, which is invalid. Name every unpacked value (`a, ep_tmp, c = ...`).
- Not every `core/tensor` function is translator-supported. Supported: `det, inv,
  log, exp, sqrt, trace, dev, eye, dyad, sym3, sqrtm, logm, expm, tanh`. **`sym`
  is NOT** — write `0.5*(F + F.T)`. `tanh` **IS** supported (emitted as
  `cs_tanh`); if you prefer you can also build it from `exp`
  (`(exp(2x)-1)/(exp(2x)+1)`). An unsupported name raises "Unsupported function:
  <name>" at generate time (verify() is unaffected).
- Plane strain still uses 3x3 tensors with `F33 = 1`.
- **`z**2` (integer power) codegens to complex `z ** DCMPLX(2.0)` = `exp(2*log z)`,
  which is **NaN at `z = 0`** (`log 0 = -inf`). For a square of any quantity that
  can be zero — e.g. `(X[0] - R0)**2` where a Gauss point may sit at `R0` — write
  the explicit product `dx*dx`, not `dx**2`. `verify()` and a direct probe at a
  generic state will NOT catch this; only an input where the base is exactly 0
  (a centroid Gauss point, a symmetric coordinate) trips it.
- `X` (Gauss-point reference position, a vector) and `time` (total time) are
  injectable UEL params (`PARAM_ARGS`): add them to a material/equation method
  signature and the generator threads `X_gp = sum_a N_a X_a` and `TIME(2)` in.
  Use them for spatiotemporal drivers `f(X)`, `C(t)`.
- UEL material generation must preserve `.real` variables as real variables.
  If a method does `x_r = x.real` and then branches on `x_r`, the generated
  declaration must be `DOUBLE PRECISION :: x_r`, not `DOUBLE COMPLEX`. UMAT
  already had this tracking; an F-bar UEL exposed the same requirement.
- Stateful material subroutine calls in UEL F-bar helpers must match the
  generated signature: all old state variables first, then `dt`. Do not blindly
  replay the Python method order when the generated subroutine groups state
  arguments differently.

## F-bar Routing

- `formulation="fbar_mechanics"` is for pure displacement UELs only. Do not
  force it onto coupled `u + damage`, `u + mu`, or thermal
  weak forms.
- If a low-order coupled model needs locking control, first add a
  pure-mechanics F-bar version to isolate the mechanical BVP. A coupled F-bar
  variant needs a separate derivation of how `Jbar` enters every
  volume-dependent storage/source/gradient term and the corresponding tangent
  chain rules.
- For a low-order coupled `u + damage` model, keep the pure-mechanics F-bar
  element and the `u + damage` gradient-damage UEL as separate formulations;
  they are not interchangeable.
- F-bar fixes the low-order volumetric element route; it does not fix a
  constitutive tangent that is already nonfinite in a single Gauss-point
  compression probe.

## Complex-Step and Matrix Functions

- Complex-step tangent checks catch residual/tangent consistency, not physical
  correctness.
- Branch on `.real` when a branch is unavoidable. Do not branch on the
  imaginary perturbation.
- A branch on `.real` is **invariant** under the imaginary complex-step
  perturbation, so an elastic/plastic (yield) switch is *not* a complex-step
  discontinuity — do not blame the branch for a tangent NaN.
- Repeated eigenvalues are dangerous for eig-based matrix functions — but
  *test* whether eig actually misbehaves at the failing state before blaming
  degeneracy. In one debugging session the degenerate states were finite in
  both the Python `eig` and the compiled `eig33z`; the NaN was elsewhere
  (a return-map formulation issue). Confirm the culprit at the exact state.
- Finite-strain UMAT generation defaults to `matrix_backend="iterative"` for
  `logm`, `sqrtm`, `expm`, and `polar` because eig-based guards can erase
  complex-step perturbations near diagonal states.
- One specific leak was an `eig33z` near-diagonal guard:
  tiny complex-step off-diagonals caused the eig backend to return `V = I`,
  zeroing shear-block derivatives in `logm33z`. Lowering the guard threshold
  only moves the failure; the robust default is the iterative backend.
- Use `matrix_backend="eig"` only deliberately, usually for compatibility or
  controlled performance comparison.
- **A CS-vs-FD failure at an extreme state can be FD noise, not a CS bug.**
  `verify()`'s FD eps is hardwired (1e-7); a heavy element evaluation
  (eig + iterative logm/expm + inner Newton) has a numerical noise floor that
  central FD divides by eps. Diagnostic: sweep the FD eps — if FD CONVERGES TO
  the CS value as eps GROWS (e.g. 9e-6 rel at eps=1e-5
  vs 1.2e-3 at 1e-7), the CS tangent is right; probe that state directly with
  an appropriate eps instead of "fixing" the model. Instrument the local
  Newton's real AND imaginary residuals before blaming it.

## UEL Weak-Form Signs

For scalar UEL equations the generated tuple convention is:

```text
RHS += -storage * N + flux . grad(N)
```

Here `flux` is a code-generation name for the coefficient of
`grad(test)`, not necessarily a physical flux. Do not add a separate
physical-flux mapping layer unless the paper itself is written as a
conservation law. The shortest safe pipeline is:

1. write the weak form term you want;
2. match it to `storage * test - flux . grad(test)`;
3. return that coefficient directly.

If the paper gives only a strong form, derive the weak form once and then use
the same direct matching. Avoid translating through "physical flux" and then
translating again into the tuple API; that extra mapping is where most sign
mistakes have entered.

Do not infer flux sign from field name. Phase and species equations can
require different storage/flux sign choices.

AT2 damage/fracture is a common trap. If
`storage = Gc/ell*d - 2*(1-d)*H`, the desired weak form contains
`+ Gc*ell*grad(d).grad(test)`, so the generated material method must return:

```python
def phase_flux(self, grad_d):
    return -self.Gc * self.ell * grad_d
```

The old `phasefield_fracture_uel` pattern with `+Gc*ell*grad_d` is a stale bad
template, not a rule to copy. Homogeneous smoke tests and CS tangent checks can
pass with the wrong sign because they only prove residual/tangent consistency,
not PDE correctness. Add a non-homogeneous gradient-sign check for every new
damage or phase-field UEL.

Do not design phase-field equations by analogy to diffusion fluxes. Diffusion
is a conserved balance with a physical flux. AT2 damage and many order-parameter
models are variational gradient equations; write the positive weak-form
gradient coefficient first and then adapt its sign to the tuple convention.

The same trap appears in a pure rate-plus-gradient equation.
For a generated equation with

```text
storage = (theta - theta_old) / dt - drive
paper: theta_dot = drive + (K / eta) * grad^2(theta)
```

the generated method must return:

```python
def phase_flux(self, grad_theta):
    return -(self.K_F / self.eta_D0) * grad_theta
```

and similarly for order-parameter gradients:

```python
def species_flux(self, grad_S):
    return -(self.K_E / self.eta_S) * grad_S
```

This is easy to miss because a dt-dominated transient equation can look stable
in homogeneous tests and can even defeat derivative-sign heuristics. Use an
assembled operator check at the resolving scale, e.g. `K = s*M - Df*L` at
`h = sqrt(|Df|/|s|)/10`, and include a deliberately flipped-flux control case.

Scope note: the built-in `OperatorSignWarning` is a screening tool, not a proof.
It covers the current isotropic scalar-gradient examples. For anisotropic
gradient energy or cross-gradient coupling, inspect the full `dflux/dgrad`
tensor spectrum or assemble the full element block directly.

## Multiplicative-Plasticity and Finite-Strain Kinematics

For `F = Fe Fp` plasticity, a recurring class of bug is exact under
coaxial/axisymmetric loading and wrong under rotation — so it passes
every diagonal-`F` test and only shows up off-axis. Check these directly,
on a deliberately NON-coaxial state (rotated principal directions, a
non-diagonal `Fp_old` that does not commute with `expm(dt*Dp)`):

- **Push-forward must use the UPDATED `Fp`.** The PK1
  `P = Re @ Me @ inv(Ue) @ inv(Fp).T` (or any `... Fp^{-T}`) must use
  `Fp_new`, not `Fp_old`. Using `Fp_old` after computing `Me` on the
  updated `Fe` is internally inconsistent and is a real stress error
  even for axisymmetric loading (it does not need rotation to bite — it
  cost 2.6% on a compression path). The Mandel stress
  and the pull-back map must reference the same plastic state. Dilatant
  plasticity (`tr Dp != 0`, `Jp != 1`) also needs the `Jp` factor the
  small-strain form drops.
- **`Fe` reconstruction factor order.** The updated elastic gradient is
  `Fe = F @ inv(Fp_new)` (equivalently `Fe_tr @ inv(expm(dt*Dp))`).
  Beware rearrangements like `Fe_tr @ inv(Fp_new) @ Fp_old` — they equal
  the correct value ONLY when `Fp_old` commutes with `expm(dt*Dp)`
  (coaxial), and silently violate `F = Fe Fp_new` otherwise. Assert
  `F == Fe @ Fp_new` numerically at a non-coaxial state.
- **Eigenvector-based slip systems / spectral reconstructions.** `V.T`
  vs `inv(V)`: the framework `eig` returns UNNORMALIZED eigenvectors
  (use `inv(V)`, not `V.T`, for reconstruction). Wrong here is exact at
  diagonal inputs and ~100% wrong at rotated inputs. At repeated
  eigenvalues the eigenvectors are
  non-unique — a genuine constitutive ambiguity, not just a numerical
  one.

General rule: treat every matrix-algebra rearrangement as unproven until
you have checked it on a non-commuting example. "Exact on a diagonal
case" is not evidence.

## Finite-Strain Tensor-Return Shape Trap

When a finite-strain UMAT/UEL material method returns a tensor, prefer assigning
`det(F)` to a scalar temporary before using it in the returned tensor expression:

```python
# Historically risky in finite-strain UMAT/UEL generation
return G * (F - inv(F).T) + K * log(det(F)) * inv(F).T
```

This formula is mathematically correct. It historically exposed a
component-wise tensor-return codegen bug that emitted invalid Fortran like:

```fortran
LOG(det33z(F(ii,jj)))
```

`det33z` needs the full `3x3` tensor, not a scalar component. The runtime symptom
is often an all-NaN element residual/tangent even for a simple undeformed element.

Use explicit scalar/tensor temporaries for readability and easier generated-code
inspection:

```python
J = det(F)
FinvT = inv(F).T
return G * (F - FinvT) + K * log(J) * FinvT
```

Regression tests now cover the direct-return `det(F)` case, but after
generation it is still worth grepping the `.for` file for impossible patterns
such as `det33z(F(ii,jj))`. If any matrix helper receives a component-indexed
tensor argument, treat that as a code-generation shape bug, not a mechanics or
solver-convergence issue.

## Local Return Maps and Nonlinear Hardening

A local Newton / implicit return map *inside* the element is supported and is
often the correct choice — do not replace it with an explicit/smooth-viscoplastic
shortcut to avoid a branch. The CS-safe idiom:
branch the elastic/plastic decision on `.real`, run a **fixed-count** `for k in
range(N)` loop (no `while`, no convergence `break` needed), update with an analytic
Jacobian; the complex step differentiates through the converged iterate and yields
the consistent tangent automatically. An explicit trial-stress update is only
conditionally stable and overshoots at stress concentrations (notch roots).

- **Solve for the right unknown.** For power-law / Ludwik hardening
  `sigma_y = sigma_y0 + K*eps_p^n` with `n < 1`, the forward slope
  `d(sigma_y)/d(eps_p) = K n eps_p^(n-1) -> inf` at `eps_p = 0`. A Newton on the
  plastic increment `dgamma` is therefore **non-smooth at the onset of yield** —
  not a solver bug, the wrong unknown. No reparametrization of `dgamma` fixes it
  (linear space overshoots to `(neg)^n` = NaN; log space overflows `exp` = NaN).
- **Fix: solve the consistency for the FLOW STRESS** via the inverse hardening
  `eps_p(sigma) = ((sigma - sigma_y0)/K)^(1/n)`, whose slope `1/n >= 1` is smooth at
  the onset (`g'(sigma_y0) = -1/(3G)`, "simple at the onset"). One CS-safe Newton
  then covers the hardening-dominated onset and the elastic-dominated high-overstress
  branch; use a root-adjacent init `sigma0 = sigma_y0 + K*(eps_p_old + f_trial/3G)^n`
  so the high-overstress case converges fast. This is standard computational-
  plasticity practice; recognise it is a *formulation* choice, not solver tuning.
- **Case-insensitive Fortran bites local-Newton temporaries.** A residual named `g`
  collides with shear modulus `G`; `dg` with anything `DG`. Name them `gres`/`dgres`.
- **Stiff exponential rate laws (`gdot = gdot0 exp(x/C) y(x)`, small C): solve in
  LOG-SLIP variables.** A Newton on the slip increment (or the explicit
  trial-rate update) explodes at implicit dt because a few percent
  overstress multiplies the rate by orders of magnitude. Change unknowns to
  `z = ln(dgam/(dt*gdot0))`; the residual `R = x + C*ln y(x) - C*z` is the exact
  backward-Euler smoothed law, needs NO active/inactive branching (inactive
  systems sit at z ~ -100, contributing exp(-100) ~ 0), and stays exp-overflow
  free with an elastic-return initial guess plus |dz| and z clamps. Early-exit
  on TWO consecutive machine-zero residuals so the CS imaginary part (a linear
  problem once the real part converges) settles.
- **Evaluate stress-state-dependent criteria (fracture loci, nucleation
  thresholds) on the RETURNED stress, not the elastic trial.** During flow the
  trial deviator overshoots the flow stress by the elastic predictor while the
  mean stress is untouched, so trial triaxiality is biased LOW, dt-dependently
  — and a Hosford-Coulomb-type locus amplifies it by the 1/n power (1/n = 10
  turned a benign-looking bias into a 3.5% locus error). If a gate is needed
  BEFORE the return map runs (dilatancy on/off), carry the previous
  increment's locus as a state variable. Also freeze (`.real`) any locus path
  through `arccos` — simple tension sits EXACTLY at xi = +1 where the
  derivative is unbounded; stress RATIOS carry no g(d)/scale tangent
  information, so freezing loses nothing that matters.

## Transcribing From The Paper

These transcription errors recurred and no consistency test catches
them. When porting equations and tables:

- **A constraint is not a definition.** A paper inequality / restriction
  (something that "must hold", `>= 0`, a dissipation condition) is not
  the defining formula. One recurring error used a dilatancy restriction
  `[tau - beta*sigma] > 0` (a positivity constraint) as the slip-rate
  numerator, when the flow rule's numerator is `tau` alone. Take
  numerators / yield functions from the DEFINITION equation; cite its
  number at the code site.
- **Cross-check every parameter against its specific table, digit and
  unit.** When a paper has sibling tables (brittle vs ductile; material
  A vs B), values migrate from the wrong table or carry over from the
  previously-coded sibling. In one case brittle `psi*`/`eta` values were
  used for the ductile mechanism; a `zeta_d=1e8` was actually correct, and
  a false review claim that it should be `1e6` was caught only by rereading
  the source table. Annotate each prop with its source table.
- **One symbol, one role.** Do not reuse a variable (e.g. an exponent)
  for two distinct physical meanings. One recurring case reused a dilatancy
  exponent `p` as the ductile-damage exponent (which the paper fixes at
  2).
- **A quantity and its variational dual must be derived separately.**
  Getting the stress sign right (derived from the energy) says nothing
  about the dual gradient-flux sign — one case got the stress right and the
  flux wrong. Derive both; check consistency.

## Coupled Dummy Materials In Abaqus

For a coupled temperature-displacement dummy element used alongside a UEL,
Abaqus may require density and specific heat in addition to conductivity:

```abaqus
*Material, name=extra_material
*Elastic
1.0e-20, 0.3
*Density
1.0,
*Conductivity
0.25,
*Specific Heat
1.0,
```

Use `au.ensure_coupled_dummy_material(...)` when rebuilding coupled UEL decks
through the input scaffold path.

## Non-Abaqus Regression Runtimes

- A non-Abaqus regression runtime is not Abaqus, and its nonlinear driver is
  typically simpler.
- Legacy UELs with approximate tangents can fail in a simpler runtime even if
  Abaqus can stabilize them.
- MUMPS/MPI-MUMPS helps large generated UEL/UMAT runs, but solver success is
  not Abaqus validation.
- Parse VTK connectivity; do not assume VTK point order equals input node-label
  order.
- Template a same-NDOFEL single-element deck and set the initial DOF at the
  model's true equilibrium so "did it solve?" is unambiguous. Use an unsymmetric
  solver for coupled/transient tangents.
- An **independent** tangent check (a runtime that FDs the residual vs the
  generated `AMATRX`) is a second confirmation beyond the Python complex-step
  `verify()`. When such a check is printed **twice** — element-level *before* the
  solve (authoritative) and *after* with BCs applied (constrained-DOF columns
  zeroed, FD noise) — read the first block, and **parse the numbers in code** —
  do not eyeball a multi-thousand-line output; you will read the wrong block or
  the wrong column.
- A non-Abaqus runtime does not call Abaqus `UVARM`. If a generated UEL exposes
  Abaqus diagnostics through `UVARM`, mirror the same quantities into `SVARS` and
  request them as field variables.
- **The deck's `PROPERTIES` list is positional and silent.** When you change a
  model's `props` (count OR order), UPDATE the deck — the runtime passes whatever
  is in the deck, so a stale list reinterprets every prop by position. This
  bites hard when a now-misaligned prop is a **denominator** (a rate constant,
  a time scale): a stale `0` there gives `0/0 -> NaN` in the residual, which
  surfaces as a tangent-check NaN, not an obvious "wrong value".
- **Diagnosis order when verify()/compile pass but the runtime NaNs:** the
  residual is NaN but the generated code is fine for *some* inputs. (1)
  Instrument the UEL kernel to print the interpolated `F, p, mu, X_gp, TIME(2)`
  per Gauss point. (2) Directly probe the generated material subroutine (compile
  the `.for`, call it from a tiny Fortran main) at those *exact* inputs. (3) If
  the probe is **finite but the runtime NaNs at identical `F/p/mu/X/time`**, the
  only remaining difference is **`PROPS`** — the deck. This isolates stale-deck
  bugs in minutes vs. spelunking the codegen. (Globals like `TIME`/`DTIME` are
  now default-initialized so the pre-solve tangent check passes a valid time.)
- A static Tecplot/state output that applies the projection solve is the
  normalized field check. VTK may show the raw integrated projection unless the
  lumped projection path is active; use VTK first as a finite/non-NaN
  visualization check, not as the only scalar-value oracle.
- For UELs, generated identifiers live in Fortran's case-insensitive namespace.
  Do not name a scalar field `d` if the material also has `d_old` state or local
  variables named `D`. Prefer explicit field names like `damage` and mirror into
  compact STATEV labels only in documentation. Also watch generated temporaries
  from tensor column slices (`V[:,0]`): declaration logic must treat `COL(ii)`
  as an indexed use of `COL`, not a new symbol.

## Driving a Generated Element in Python (reference_assembly)

`core/reference_assembly.assemble_element(weakform, coords, U, DU, DTIME, PROPS,
element=...)` runs the Python weak form through real element assembly. It is the
right tool for a physics BVP / size-effect benchmark (it exercises the element +
B matrix + quadrature — though the Python material, not the generated Fortran).

- `coords` is shape `(ndim, n_nodes)` — transpose your `(n_node, 2)` node list.
- Pass `element="Quad4"` (or the intended config). The default layout config is
  **Quad8 geometry** (degree-1 fields give 12 DOFs but 8 geometric nodes), so a
  4-node `coords`/`U` mismatches and indexes out of bounds.
- DOF order is node-major, field-interleaved (`[ux, uy, p, ...]` per node).
- **AMATRX is returned in the Abaqus convention `AMATRX = -d(RHS)/dU`.** A global
  Newton step is therefore `K dU = +R`, NOT `-R`. An FD check on the assembled
  tangent gives `fd == -K` exactly (relative error 2.0 on every column); the
  wrong sign diverges and *looks like* a conditioning problem. Confirm the sign
  with a one-column FD before trusting any hand-rolled solver.
- A stiff coupled solve from rest (viscoplastic, `e_p/eps_eq`, damage) overshoots
  a near-singular soft mode at the full Newton step; add a backtracking line
  search (the descent direction is `solve(K, +R)`). Match `dt` to the viscoplastic
  rate (`dt ~ gamma/(n_steps*edot0)`) or the overstress `(seq/sigflow)^(1/m)`
  blows up.

## Hard Paper Models

- Do not jump to a final published figure. Build framework rungs first.
- Missing parameters are common. Do not invent them from memory.
- Preserve physics, not necessarily the paper's software implementation frame.
- "Eulerian" in a paper may mean current-state variables and spatial fluxes,
  not an Eulerian mesh. A total-Lagrangian UEL can pull current quantities back
  with `J`, `F^{-1}`, and `C^{-1}`.
- Use lower-order pressure for poromechanics until stability is demonstrated.
- For swelling/plasticity papers, check whether "damage", "degradation", or
  "swelling" is a real scalar damage law, eigenstrain, residual plastic strain,
  or simply output terminology.

## Debugging Order

When something fails:

1. Check generated `.for` compiles.
2. Check Python reference at the exact state.
3. Check generated Fortran through f2py if possible.
4. Check sign convention and state layout.
5. Check active DOF declaration, `UNSYMM`, `*Depvar`, property count, and dummy
   thermal material blocks.
6. Check time step, boundary conditions, hourglass control, and solver controls.
7. Only then alter the constitutive theory.

Do not skip ahead by guessing the cause. For each suspected culprit (eig
degeneracy, a `.real` branch, the linear solver, the compiled-vs-Python path),
**reproduce the exact failing state and test that hypothesis** before fixing it —
cross-check the compiled element against `reference_assembly`, and isolate which
sub-computation goes non-finite. A wall-clock timeout is not evidence of
non-convergence (an f2py recompile alone is ~60 s; compile once per debug run).
When the local solve itself diverges or NaNs, suspect the *formulation* (the wrong
unknown — see "Local Return Maps and Nonlinear Hardening") before solver tuning.

## Verification: the State and the Regime, Not Just the Method

`verify()` (CS-vs-FD), generate+compile, and code-vs-itself f2py are
**consistency** checks — they cannot see convention errors or untested branches.
One audit found a real gradient-damage bug that all of these
passed. Add a **material-point regime sweep**:

- **Verify every branch, not one state.** `verify()` runs a single *tensile* state;
  a model with a tensile (brittle) and a compressive (ductile) branch needs BOTH
  exercised. The bug and the nonfinite tangent lived only in the unverified branch.
- **CS-vs-FD passes "wrong-but-consistent".** A mismatched `psi_star`/`ell` in a
  gradient term, a degraded-vs-undamaged drive, or a scrambled call-arg order agree
  on both sides. Confirm generator arg order by reading the emitted `.for`.
- **Assert the test entered its regime.** `assert_regime_entered` — a qualitative
  `gamma_bar_p > 1e-3` passes even when the path stayed elastic (uniaxial *strain*
  over-confines → no yield). A vacuous pass is worse than a red.
- **Regularization self-consistency** (gradient damage): measure `ell_eff` from the
  gradient/reaction stiffness ratio in each branch, assert `== ell`. `phase_flux`
  must use the same branch `psi_star` (and `ell`) as `phase_storage` — the paper
  gives each mechanism its own ψ\*, ℓ, η.
- **Robustness is a gate.** Tangent finite across softening; dt-refinement converges
  (distinguishes an explicit-overshoot artifact from a constitutive defect — the
  stiff `1/m` slip flow overshoots to `d=1`+NaN at large per-step strain, but the
  peak stress converges under dt-refinement, so the fix is an implicit/sub-stepped
  return map, not the equations).
- For an open known bug, pin the CORRECT expected value with `xfail(strict)`.
  Remove the marker as soon as the fix turns it into XPASS.
