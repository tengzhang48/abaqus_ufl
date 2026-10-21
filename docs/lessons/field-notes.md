# Field notes compatibility map

The former aggregate repeated themed lessons and retained conflicting advice.
Read the [task index](README.md) first. These headings preserve existing links;
the linked pages own the explanations.

## Verification: what a green suite actually certifies

See [verification and evidence](verification.md) and
[finite-difference consistency](codegen-and-testing.md#the-finite-difference-tangent-check-is-the-primary-automated-check).

## Solvers and performance

See [performance and parallel hosts](performance-and-parallel.md), including
representative profiling, field-wise scales and preconditioning.

## Distributed and MPI

See [distribution and environment](performance-and-parallel.md#keep-replicated-batching-distinct-from-distributed-algebra).
A machine-specific launcher observation is not a universal vendor rule.
The former commit-time-clipping recommendation is superseded by
[pivot P6](../pivots.md#p6-clipping-is-not-constrained-equilibrium).

## Code generation and elements

See [code generation](codegen-and-testing.md), [finite elements](finite-element.md)
and the [specific mistake/check catalog](../../ai_skills/abaqus-ufl/references/pitfalls.md).
Integer powers and supported DSL operations are specified by the
[API guide](../API_USAGE.md#tensor-dsl-rules); the former warning that every
square lowers through a complex logarithm is obsolete.

## Constitutive and numerical formulation

See [constitutive integration](constitutive-integration.md).

## Reproducing papers and debugging discipline

See [paper reproduction](coupled-multiphysics.md#reproducing-multiphysics-papers)
and [verification evidence](verification.md). Symptoms suggest hypotheses;
they do not uniquely identify a cause.

## Reuse and architecture

See [agent workflow](agent-workflow.md#keep-one-canonical-implementation-and-index)
and [state ownership](performance-and-parallel.md#exact-reuse-needs-complete-inputs-and-state-ownership).
The [host API](../../abaqus_ufl/fe/README.md) owns supported reuse and commit calls.
