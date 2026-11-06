# Start here for lessons, pitfalls, pivots and skills

Choose by the question you are answering. Each item has one canonical home;
other pages link to it. The package is a **formulation framework**.

| Kind | Question it answers | Canonical home |
|---|---|---|
| Lessons | What did we learn, and why does it matter? | The themed pages in [this directory](#lesson-homes) |
| Pitfalls | What specific mistake causes this symptom, and which check catches it? | [Pitfall catalog](../../ai_skills/abaqus-ufl/references/pitfalls.md): one row per symptom, linking to the entry that explains it |
| Pivots | What triggered a change of approach, what was the evidence, and what replaced it? | [Pivot records](../pivots.md) |
| Skills | When should an agent use which resource or workflow? | [Skill entry point](../../ai_skills/abaqus-ufl/SKILL.md) and its [references](../../ai_skills/abaqus-ufl/references/validation-ladder.md) |

The [API guide](../API_USAGE.md), [host API](../../abaqus_ufl/fe/README.md) and
[example pipeline](../../HOWTO_ADD_AN_EXAMPLE.md) own executable usage.
Lessons explain decisions; they do not extend the public capability boundary.
An internal benchmark observation is motivation unless linked public evidence
supports a more specific claim.

## Task-to-item map

| Task or symptom | Read this item first | Find the check or changed approach |
|---|---|---|
| Writing a 2D UEL or plane-strain UMAT | [Active dimensions](finite-element.md#active-dimensions) | [In-plane pivot](../pivots.md#p1-in-plane-complex-step), [compiled parity](../../tests/test_2d_tangents.py) |
| A swelling model fails the stress-free check | [Reference assumptions](verification.md#state-the-reference-assumed-by-a-gate) | [Pitfalls](../../ai_skills/abaqus-ufl/references/pitfalls.md#verification-the-state-and-the-regime-not-just-the-method) |
| Choosing element order for a stress-gradient term | [Gradient consistency](finite-element.md#gradient-of-stress-consistency) | [Interpolation pivot](../pivots.md#p2-compatible-interpolation-for-stress-gradient-coupling) |
| Choosing element order when a scalar drives swelling or thermal strain | [FE1 Eigenstrain on equal-order elements](finite-element.md#fe1-a-scalar-driven-eigenstrain-on-equal-order-low-order-elements-can-introduce-parasitic-stress) | [Interpolation pivot](../pivots.md#p2-compatible-interpolation-for-stress-gradient-coupling), [mixed-order example](../../examples/thermo_mechanics_quad8/) |
| Stress, damage or history too high on coarse meshes although Newton converges | [FE1 Eigenstrain on equal-order elements](finite-element.md#fe1-a-scalar-driven-eigenstrain-on-equal-order-low-order-elements-can-introduce-parasitic-stress) | [V1 Error size, not only its trend](verification.md#v1-a-falling-error-establishes-neither-the-limit-nor-an-acceptable-error) |
| A global Newton step goes uphill | [AMATRX sign](abaqus-and-fortran.md#amatrx-sign-convention) | [Host/interface pitfalls](../../ai_skills/abaqus-ufl/references/pitfalls.md#abaqus-interface) |
| Cache safety, hits or state commit | [Exact reuse](performance-and-parallel.md#exact-reuse-needs-complete-inputs-and-state-ownership) | [Cache pivot](../pivots.md#p5-exact-batch-reuse), [cache tests](../../tests/test_fe_batch_cache.py) |
| MPI disagrees with serial or hangs | [Distribution and error handling](performance-and-parallel.md#keep-replicated-batching-distinct-from-distributed-algebra) | [MPI probes](../../tests/mpi_fe_contract.py), [runtime command](../../abaqus_ufl/fe/README.md#petsc-and-mpi) |
| More ranks are slower | [Representative profiling](performance-and-parallel.md#measure-the-representative-bottleneck) | [MPI environment](performance-and-parallel.md#diagnose-the-mpi-environment-before-the-algorithm) |
| Tangent checks pass but the physics is wrong | [What FD certifies](codegen-and-testing.md#the-finite-difference-tangent-check-is-the-primary-automated-check) | [Model/setup/result distinction](verification.md#separate-model-setup-and-reported-quantity) |
| A CS tangent fails at an extreme or repeated spectrum | [Matrix-function lessons](complex-step-and-matrix-functions.md) | [CS pitfalls](../../ai_skills/abaqus-ufl/references/pitfalls.md#complex-step-and-matrix-functions) |
| Rotated plastic states fail while diagonal states pass | [Objectivity and isotropy](constitutive-integration.md#objectivity-and-material-isotropy-are-different-tests) | [Kinematic pitfalls](../../ai_skills/abaqus-ufl/references/pitfalls.md#multiplicative-plasticity-and-finite-strain-kinematics) |
| Damage drifts over long cycles | [History time level](constitutive-integration.md#lagging-history-changes-temporal-accuracy) | [Long-history verification](verification.md#validate-long-histories-and-nonunique-patterns-honestly) |
| A fatigue floor or accumulator loses precision | [Accumulator arithmetic](constitutive-integration.md#protect-accumulator-arithmetic-and-extreme-driving-forces) | [Clipping pivot](../pivots.md#p6-clipping-is-not-constrained-equilibrium) |
| A mesh-order conversion changes the solution | [Topology-dependent sets](finite-element.md#mesh-order-changes-require-new-set-membership-checks) | [DOF pitfalls](../../ai_skills/abaqus-ufl/references/pitfalls.md#dof-and-state-layout) |
| Solver output or a figure seems inconsistent | [Output bridge](coupled-multiphysics.md#solver-completion-is-not-output-validation) | [Retained artifact evidence](verification.md#a-claim-belongs-to-its-retained-artifact) |
| Starting a paper reproduction | [Build smaller rungs](coupled-multiphysics.md#hard-multiphysics-needs-framework-rungs-before-parameter-tuning) | [Validation ladder](../../ai_skills/abaqus-ufl/references/validation-ladder.md) |
| Choosing an agent workflow or finding the right copy | [One agent and an oracle](agent-workflow.md) | [Workflow/canonical-copy pivots](../pivots.md#p3-one-agent-with-an-independent-oracle) |

## Entry format

New and rewritten entries use this form; other entries move to it when they
are next edited. Aim for 80 to 200 words, but never drop a limit to save space.

```markdown
## <ID> <the rule or the failure, as one sentence>

Also: <other words an agent may search for>

- **When:** the problem class and assumptions in which it applies.
- **Symptom:** what an agent or user observes, including results that look like success.
- **Do:** the action. **Why:** the reason it works.
- **Check:** the decisive check: quantity, reference, conditions and acceptance.
- **Evidence and limits:** a named test, example or derivation; what was shown and what was not.
- **Checked:** the commit, and whether the source was inspected or a test was executed.
```

An ID is a page prefix and a number (FE1 in finite-element.md, V1 in
verification.md) and is never reused. When a heading changes, keep its old
anchor with `<a id="..."></a>` so existing links still work. A pitfall row
names the mistake and its symptom and links to the entry; it does not repeat
the explanation. Merge two entries only when they guard the same decision
(the same *when* and the same *do*), not because they share words.
`python tools/check_guidance.py` checks links, anchors and named tests, and
the test suite runs it.

## Lesson homes

- [Abaqus and Fortran](abaqus-and-fortran.md): interface conventions, fixed form and deck contracts.
- [Finite elements](finite-element.md): kinematics, mapping, interpolation, F-bar and mesh order.
- [Complex step and matrix functions](complex-step-and-matrix-functions.md): derivative-safe operations, spectra and generator details.
- [Code generation and testing](codegen-and-testing.md): authoritative generation, consistency checks and DSL design.
- [Coupled multiphysics](coupled-multiphysics.md): scaling, locking, weak-form signs and output bridges.
- [Constitutive integration](constitutive-integration.md): local updates, rotated state, history time levels and precision.
- [Performance and parallel hosts](performance-and-parallel.md): profiling, active work, exact reuse and MPI contracts.
- [Verification and evidence](verification.md): assumptions, falsifiable gates, observables and retained evidence.
- [Agent workflow](agent-workflow.md): task scope, independent oracles, canonical copies and reproducible handoffs.

[Field notes](field-notes.md) is a compatibility map for former aggregate
headings; it contains no second copy of the lessons. For additions, choose a
home, add a task-map entry if needed, and cross-link the symptom/check or pivot
rather than reproducing its explanation.
