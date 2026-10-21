# Changes of approach

A pivot records a trigger, evidence and replacement. Explanations belong in
the linked lessons; operational mistakes belong in the
[pitfall catalog](../ai_skills/abaqus-ufl/references/pitfalls.md).
Return to the [four-kind index](lessons/README.md) to choose an item by task.

## P1 In-plane complex step

- **Trigger:** 2D UEL tangents spent constitutive calls on columns never assembled.
- **Evidence:** compiled comparisons cover standard Quad4, mixed Quad8, F-bar,
  stateful updates and plane-strain UMATs; original 3D source hashes remain fixed.
  The public Ogden timing measures a kernel speedup, not global scaling.
- **Replacement:** keep full 3x3 kinematics but perturb active in-plane entries;
  omit provably unused state-output calculations in pure straight-line methods.
- **Check and limits:** [regressions](../tests/test_2d_tangents.py) compare with
  full perturbations. Opaque/control-flow updates retain the full method.
  See [the lesson](lessons/performance-and-parallel.md#perturb-only-active-columns).

## P2 Compatible interpolation for stress-gradient coupling

- **Trigger:** a refined equal-order bilinear solve converges in Newton yet
  retains a stress-gradient or swelling-driven damage bias.
- **Evidence:** independent mixed-order comparisons in an internal benchmark
  distinguished interpolation error from quadrature and convergence error.
- **Replacement:** derive a continuous auxiliary stress-like field, use
  displacement one order higher where appropriate, or derive condensed local
  enhanced modes when topology is fixed.
- **Check and limits:** compare integrals and fields under mesh/time refinement
  with an independent higher-order solution. These are formulation choices,
  not interchangeable generator switches or new public enhanced-mode support.
  See [the lesson](lessons/finite-element.md#gradient-of-stress-consistency).

## P3 One agent with an independent oracle

- **Trigger:** the skill called a multi-agent implement/review/fix loop the
  reference pattern regardless of task.
- **Evidence:** an internal benchmark found higher cost without better element
  accuracy from that workflow; the independent reference exposed the useful gaps.
- **Replacement:** default to one agent and seek an independently derived
  oracle/reference solution. Delegate only when requested and useful.
- **Check and limits:** require branch-reaching counterexamples and quantitative
  agreement. A second agent or a second evaluation of the same formulas does
  not establish scientific independence. See [the lesson](lessons/agent-workflow.md#default-to-one-agent-and-an-independent-oracle).

## P4 One canonical copy

- **Trigger:** divergent packages, installed skills and lesson copies hide fixes
  and let agents follow stale guidance.
- **Evidence:** the public skill repeated pitfalls and review advice; field
  notes duplicated themed lessons and retained superseded recommendations.
- **Replacement:** one task index, one explanation per topic, a concise routing
  skill and links from compatibility pages. Develop the public line through
  reviewed feature branches.
- **Check and limits:** check link targets and the executing revision; preserve
  dated evidence separately. This draft does not establish that every external
  package/skill installation has been synchronized. See [the lesson](lessons/agent-workflow.md#keep-one-canonical-implementation-and-index).

## P5 Exact batch reuse

- **Trigger:** line-search checks and commit can repeat a successful full batch.
- **Evidence:** cache-on/off regressions compare R/K and state exactly; MPI
  probes exercise one-rank misses, disabled caches and rejected commits.
- **Replacement:** opt-in exact input snapshots and immutable saved outputs,
  including trial state, with a unanimous MPI hit vote.
- **Check and limits:** [cache tests](../tests/test_fe_batch_cache.py) and
  [MPI probes](../tests/mpi_fe_contract.py). Kernels with required external side
  effects must keep per-call execution. See [the lesson](lessons/performance-and-parallel.md#exact-reuse-needs-complete-inputs-and-state-ownership).

## P6 Clipping is not constrained equilibrium

- **Trigger:** older guidance described commit-time clipping as the preferred
  remedy for active-set cycling.
- **Evidence:** a clipped field can violate the equilibrium residual even if
  the preceding unconstrained Newton solve converged.
- **Replacement:** distinguish a numerical safeguard from a bounded or
  irreversible formulation certified by its own projected/KKT residual.
- **Check and limits:** recompute the appropriate residual after any projection
  and verify rollback of genuinely mutated state. No bounded/staggered driver
  is added by this update. See [the lesson](lessons/constitutive-integration.md#protect-accumulator-arithmetic-and-extreme-driving-forces).
