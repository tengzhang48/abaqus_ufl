# Lessons: Coupled Multi-Physics

General lessons for building and validating coupled multi-physics user
elements (UELs): diffusion–deformation (swelling gels), phase-field and
gradient-damage, and thermo-hygro-mechanical problems. The recurring theme is
that a coupled UEL can pass every local (single-element, tangent-consistency)
check and still be wrong at the system scale — because of conditioning, locking,
an operator sign, or an output path that never gets exercised by a smoke test.

## Scaling & conditioning

### State variables must match initial material properties

A common trap is hardcoding an initial state variable in the material class
(`state_vars = dict(phi=0.99)`) while the constitutive law also accepts an
initial value as a property (`phi0`). When a user passes a different `phi0`, the
state variable stays at the hardcoded default, producing a non-zero residual
stress at the reference configuration:

```
P = (K/phi) * log(phi/phi0) * I   →   ~1 MPa for K=10 MPa, phi=0.99, phi0=0.9
```

That residual stress drives spurious deformation and causes Newton divergence
from the first step. Sync state variables with properties in `__init__`:

```python
class GelMaterial(au.Material):
    props = dict(G=1.0, K=100.0, phi0=0.99)
    state_vars = dict(phi=0.99)  # fallback only

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.state_vars = {'phi': float(self.phi0)}
```

### Scale the transport equation to match mechanical stiffness

In SI units the mechanical stiffness block (u–u) is enormous compared with the
chemical/transport block (μ–μ). For diffusion–deformation the mismatch can span
~15 orders of magnitude: mechanical stress scales with `G ~ 1e5 Pa`, while a
transport flux scaling with `D/(R·T)` for `D ~ 1e-11 m²/s` is O(1e-17). The
stiffness ratio is roughly `G·R·T/D ≈ 1e15`. A direct solver cannot resolve the
coupled system; NaN/Inf propagate and assembly fails.

Fix: give the material a `scale` property and multiply every transport output
(flux and storage) by it. This is equivalent to row-scaling the Jacobian for the
transport equation.

```python
class GelMaterial(au.Material):
    props = dict(G=1e5, K=1e7, D=5e-11, scale=1e9)

    def solvent_flux(self, F, mu, grad_mu, phi_old, dt):
        ...
        return -m * Cinv @ grad_mu * self.scale

    def solvent_storage(self, F, F_old, mu, phi_old, dt):
        ...
        return phi_dot / (J * Vmol * phi**2) * self.scale
```

Choose `scale` so the μ–μ block is comparable to the u–u block; for elastomeric
gels at SI scale a value near `1e9` is typical. This changes the *units of the
residual*, not the solution: the solver sees a well-conditioned matrix and the
DOFs converge to the same physical values (μ still in J/mol).

### Use one consistent unit system and verify every conversion

Pick a single consistent unit system (for example MPa–mm–s) and stay in it. A
constant that looks wrong is often correct after conversion, not a bug: e.g. a
coefficient of `35.3` in MPa–mm–s equals a published `35.3 × 10⁶ J/m³`. Before
running a benchmark, tabulate each parameter as *(code value, source value,
conversion)* and confirm they match.

The dangerous cases are placeholder values that are dimensionally plausible but
off by orders of magnitude — a Young's modulus left at `190` instead of `190e3`,
or a mobility/relaxation constant left at `1.0`. These do not error; they
silently produce the wrong response. Audit them explicitly before trusting any
comparison curve.

### D = 0 does not mean "no transport"

Setting the diffusivity to zero for an inert layer (e.g. an elastic backing)
does **not** stop the material routine from evolving the internal field. The
chemical-potential equation is still solved, a local Newton step can still move
`phi` away from `phi_old`, and that produces a non-zero storage rate and
spurious flux. Branch on `D == 0` in *both* the flux and storage methods:

```python
def solvent_flux(self, F, mu, grad_mu, phi_old, dt):
    if self.D == 0:
        return au.VectorField.zeros(3)
    # ... normal flux computation

def solvent_storage(self, F, F_old, mu, phi_old, dt):
    if self.D == 0:
        return 0.0 * self.scale
    # ... normal storage computation
```

Note: material-method code is translated by a restricted AST, so `np.zeros()`
and helper-function calls are not available inside these methods. Use the
framework's zero constructors (`au.VectorField.zeros(3)`) or inline literals.

## Locking & element choice

### F-bar can destabilize coupled problems at SI scale

F-bar cures volumetric locking in *pure mechanical* low-order elements by adding
a rank-1 correction to the tangent. In an SI-scale coupled problem
(`G ~ 1e5`, `K ~ 1e7`, `h ~ 1e-4 m`) that correction scales like `K·h⁻² ≈ 1e15`,
dominating the standard stiffness and driving condition numbers that crash the
solver.

- Non-SI parameters (`G=1`, `K=100`, `h~1`): F-bar converges in 1–2 iterations.
- SI parameters (`G=1e5`, `K=1e7`, `h~1e-4`): F-bar causes immediate assembly
  failure.

For coupled diffusion–deformation at SI scale, prefer instead:

1. A three-field (u, p, μ) formulation with a proper mixed discretization
   (e.g. lower-order pressure with continuous displacement and μ).
2. Higher-order displacement elements, which lock less severely.
3. Disabling F-bar only as a last resort — the run converges but shows locking
   artifacts (artificially stiff response).

This does **not** invalidate F-bar for pure mechanics: a near-incompressible
Neo-Hookean low-order cantilever with `K/G = 1000` runs cleanly with an F-bar
element and gives the expected larger, locking-relieved tip deflection. Keep
F-bar as the mechanical low-order technology; do not make coupled
swelling/diffusion production depend on it.

### Volumetric locking in two-field (u, μ) low-order elements

A two-field formulation (displacement + a scalar transport field) on a linear
element cannot represent volume-preserving shear when `K/G` is large. The
element locks.

Symptoms:

- Converges in one iteration (the problem is effectively linear within a step).
- Response is too stiff — displacements smaller than expected.
- Stress ratio `σ_yy/σ_xx ≈ 1` in a bending test (should be near 0).

Fix: add a pressure field (three-field u–p–μ) or use higher-order displacement
elements. The three-field split decouples volumetric and deviatoric response and
removes the locking without F-bar.

### Local-pressure condensation is the preferred low-order coupled path

For coupled swelling/diffusion on low-order elements, treat the pressure-like
volumetric variable as an *element-local* scalar when its equation is algebraic.
Solve it inside the element with a scalar Newton iteration, store it in the
element state, and statically condense it out of the global system:

```
K_cond = K_xx - K_xp * inv(K_pp) * K_px
R_cond = R_x  - K_xp * inv(K_pp) * R_p
```

This leaves a global system in the mechanical and transported nodal variables
(e.g. `u, μ`) while keeping a pressure variable that relieves volumetric
locking.

Critical detail: the condensed tangent must include *every* coupling to the
local pressure — not only the mechanical blocks (`K_up`, `K_pu`) but also the
transport–pressure blocks (`K_μp`, `K_pμ`), because the volumetric field,
storage, and mobility all depend on pressure through the constitutive law.

Do **not** condense diffusive variables. A field with gradients, boundary
conditions, or a global conservation requirement (chemical potential,
temperature, concentration) must stay a global nodal field. Local condensation
is only for genuinely algebraic element-local variables.

Recommended low-order coupled path: linear elements with local pressure
condensation for production, keeping a global three-field higher-order path as a
robust reference and high-order option.

## Operator signs in coupled weak forms

### Derive the weak form; tangent checks cannot prove the PDE sign

Scalar transport/phase equations are assembled from a fixed tuple convention.
The generated residual is:

```
RHS += -storage * N + flux . grad(N)
```

equivalently

```
∫ storage * η dV  -  ∫ flux . grad(η) dV = 0
```

Finite-difference / complex-step tangent checks confirm that the generated
tangent is consistent with the generated residual — but they **cannot** prove
the residual has the intended PDE sign. A residual with a flipped gradient term
is internally tangent-consistent and passes `verify()`, yet the PDE is
anti-diffusive/anti-elliptic: refining the mesh makes it *worse*, and stronger
bound penalties only mask the instability.

For every new equation template, start from the weak-form term you want and map
it directly onto `storage * test - flux . grad(test)`. If the paper gives only a
strong form, derive the weak form once. Do not introduce an intermediate
"physical flux" mapping unless the equation is a genuine conserved balance —
that mapping is exactly where sign errors hide.

### The flux sign depends on the storage convention, not the field name

Two equations can share the same gradient operator yet require opposite flux
signs, purely because their `*_storage()` methods return different quantities.

For a phase equation whose storage returns `-φ̇/L - ∂ψ/∂φ`, the gradient
contribution keeps the **positive** conductivity sign:

```python
def phase_flux(self, F, phi, grad_phi):
    return self.kappa * grad_phi
```

For a species equation whose storage returns `ċ`, the same-looking diffusion
term must be returned **negated**, because the weak form places the physical
diffusion vector on the opposite side:

```python
def species_flux(self, F, phi, grad_phi, grad_c):
    return -self.Dm * (grad_c - dh * (cSe - cLe) * grad_phi)
```

Derive the sign from the exact quantity each `*_storage()` method returns, never
from the field's name.

### Gradient-damage / AT2 has the opposite flux sign from a phase equation

Do not treat any one phase equation as the template for all phase-field or
damage models. For a common AT2 damage equation with history field `H`,

```
∫ [Gc/ℓ · d - 2(1-d)H] η dV  +  ∫ Gc·ℓ · grad(d) . grad(η) dV = 0
```

the material methods under the tuple convention are:

```python
def phase_storage(self, d, H):
    return self.Gc / self.ell * d - 2.0 * (1.0 - d) * H

def phase_flux(self, grad_d):
    return -self.Gc * self.ell * grad_d
```

Returning `+Gc·ℓ·grad(d)` gives the generated residual the wrong gradient sign.
It stays tangent-consistent — `verify()` and homogeneous single-element tests
pass — but the PDE is anti-elliptic in the damage gradient. For phase-field and
damage models, do not think of the return value as a diffusion flux: write the
positive weak-form gradient term first, then return whatever sign the tuple API
requires.

**Test rule:** add at least one *non-homogeneous* scalar-field check for every
new phase/damage element. A homogeneous smoke test cannot see the gradient sign.
Good discriminators: a 1-D diffusion/damage operator-sign test, a
positive-definiteness check on the scalar gradient block after applying the
generator convention, or a manufactured profile with a known Laplacian sign.

### Match caller/callee ordering for multi-scalar history arguments

In a multi-scalar coupled element (`u-φ-c`, `u-φ-μ`, …) a subtle, silent bug is
an ordering mismatch between how the tangent subroutine *declares* its scalar
history arguments and how the caller *passes* them. If one side is alphabetical
and the other follows field-declaration order, a routine expecting
`(c_old, phi_old)` can be handed `(phi_old, c_old)`. The residual path still
runs, so it surfaces as poor Newton convergence rather than a compile error.

Derive scalar history arguments from a single shared ordered helper tied to the
declared scalar-field order, and make the regression test deliberately use a
non-alphabetical order (e.g. `phi_old` before `c_old`), checking both the
generated subroutine signature and the call site. Relatedly, represent known
`*_old` field/state arguments symbolically in generated tangents rather than
baking their default initial values into the expression.

## Reproducing multiphysics papers

### Test simpler before complex

When bringing up a new coupled simulation, climb a ladder of increasing
difficulty rather than jumping to the full problem:

1. **Single coupled element**, homogeneous response — validates the local
   (e.g. `phi`) solve.
2. **Single block, coarse mesh**, free response — validates transport and time
   stepping.
3. **Bilayer with a soft second material** (`G₂ ≈ G₁`) — validates the
   interface.
4. **Bilayer with a stiff second material** (`G₂ = 500 × G₁`) — the real
   problem.

Skipping to step 4 produces combinatorial debugging where mesh, interface,
scaling, locking, and inert-layer transport are all suspect at once. A coarse
free-response test at step 2 will catch scaling and state-variable bugs in a
handful of steps; the stiff bilayer would be nearly impossible to debug without
those fixes already in place.

### Hard multiphysics needs framework rungs before parameter tuning

A hard multiphysics paper is often not a local material-subroutine problem: the
essential unknowns may be displacement, temperature, liquid and vapor
concentrations, and gas pressure simultaneously. Treating it as a UMAT would
miss the physics. The right abstraction is a mixed UEL with lower-order
pressure, exposed-surface boundary residuals, and time-dependent process
schedules.

Build it in rungs — reduced fields, then field splits, then boundary elements,
then released mechanics, then the low-order pressure, then comparison scripts.
Each rung makes the next gap visible. Only once the element, boundary terms,
tangents, a full solve, and postprocessing all work does the open question
change from "can the framework express this model?" to "what geometry, surface
scaling, source strength, outlet condition, and parameter set reproduce the
paper?" Do not tune parameters before the framework rungs are in place.

### "Eulerian" physics can live in a total-Lagrangian UEL

The word "Eulerian" in a paper does not require an Eulerian mesh. It usually
means the variables are current-state quantities: current porosity,
current-volume concentrations, spatial gradients, current fluxes. Preserve those
physics in a total-Lagrangian element by pulling them back:

```
dv       = J dV
grad_x a = F^{-T} grad_X a
flux_X   = J C^{-1} grad_X a
```

This is a formulation choice, not a physics compromise. The mistake is copying a
solver-specific implementation literally instead of reformulating the same
balances into the element framework.

### Make the framework observable before production

Run the model through a lightweight reference FE driver before the production
solver, and use it to expose *every* layer: the mixed DOF layout, generated
tangent consistency, boundary-element dispatch, each solved field, the released
mechanics, the low-order pressure, and paper-style postprocessing. When those
are all observable and working, a later disagreement with the paper is
interpretable as model fidelity and calibration — not a hidden sign convention
or assembly bug. Treat this reference gate as mandatory for hard UEL papers, not
optional.

### Solver completion is not output validation

A run can compile, link, converge, and complete the full load schedule while the
quantities needed for comparison are unusable — state fields written as
non-finite, or a derived history stuck at zero. That is an execution pass *and*
an output-validation failure, not a reproduction. It is especially common when a
UEL relies on duplicate visualization elements and a dummy-material bridge to get
its internal fields into the output database.

Validate the bridge explicitly, as separate checks:

1. **Solver execution** — the job compiles, links, and reaches the target step.
2. **UEL state validity** — selected integration-point variables are finite
   *inside* the UEL before being written to visualization storage.
3. **Bridge validity** — the dummy material copies the same finite values into
   its state array.
4. **Extractor validity** — the postprocessor reads the intended fields and
   frames.
5. **Response validation** — curves and contours match the reference / paper.

Use a short *serial* debug run before any threaded or long run, and compare one
known integration point on both sides of the bridge (UEL-side visualization
value vs. the value the dummy material writes into its state array). Module-level
scratch arrays and global scalars used for visualization are often unsafe under
threaded execution — an N-core request can run as one rank × N threads — so
thread safety matters even when the job is launched through an MPI-style script.

Also, do not read a zero or scaled visualization field as proof that the local
solve failed. A missing or lumped projection path can show zeros while the
element state updated correctly. Confirm state validity at the source before
blaming the physics.

### Separate scalar-transfer failures from state-history failures

Coupled-UEL bring-up often hides *two* independent failure modes that look the
same from the output:

- **Scalar transfer / DOF activation.** When a UEL uses coupled scalar DOFs
  (e.g. temperature-like DOFs beyond the mechanical ones), the analysis code may
  not pass useful scalar values into the UEL until a native coupled dummy element
  and matching initial-condition block are added. Until then the scalar output
  path is zero or unusable — a *setup* issue, not a UEL bug.
- **Uncommitted state history.** After scalar transfer is fixed, a separate bug
  can remain: an early nonzero-time call sees zeroed element state before the
  first state update is committed. A term like `-φ̇ / x_old` with `φ̇ = 0` and
  `x_old = 0` then evaluates `0/0 → NaN`, even though the DOFs transfer
  correctly.

For property-initialized scalar state variables that are physically positive,
guard the state read and fall back to the mapped property when the stored value
is NaN, huge, or non-positive. This is independent of tangent verification:
complex-step tangent checks can pass while the analysis-code call sequencing
still exposes an uncommitted-state path.

### A short output smoke before full-mesh production

After a one-element fix passes, the full mesh still needs its own short run
before a production run. The correct intermediate gate is a shortened run on the
production mesh that shows *finite, nonzero* visualization fields at a large
fraction of integration points — proving the big mesh receives the scalar DOFs
and the visualization bridge writes useful fields — without yet claiming the
physical growth/response is reproduced.

Keep three distinct milestones and do not collapse them:

1. one-element scalar-transfer pass,
2. full-mesh finite-output smoke pass,
3. full-duration paper-comparison pass.

A derived quantity that stays at finite zero in the short smoke run is a
production/reproduction target, not evidence from the smoke test.

### Long-run convergence is a separate validation gate

A short run with finite output proves that scalar transfer and visualization
work; it does **not** prove the full nonlinear schedule is stable. A corrected
deck can reach far deeper into the load schedule than the earlier broken-output
runs and still fail late — typically on mechanical convergence
(displacement-correction warnings, negative eigenvalues, too many attempts at an
increment), a different failure mode from the earlier scalar-NaN path.

Treat this as its own gate. When debugging late-time failure, change one control
at a time: reduce the maximum increment size, review boundary/load release near
the failing time, and compare the convergence pattern against a reference
implementation.
