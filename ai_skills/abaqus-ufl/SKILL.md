---
name: abaqus-ufl
description: Develop, review, verify or debug abaqus_ufl UMAT/UEL declarations, generated Fortran, f2py checks and Abaqus integration examples. Route formulation, tangent, state, host and validation questions to the repository's canonical lessons and checks.
metadata:
  short-description: Develop and verify abaqus_ufl UMAT/UEL models
---

# abaqus_ufl model development

Use this skill for this formulation framework's material and element work.
The skill is routing metadata; explanations and checks have separate homes.

## Read first

Read the [task-to-item index](../../docs/lessons/README.md) first. It distinguishes
lessons (what and why), pitfalls (mistake, symptom, check), pivots (changed
approach and evidence), and skills (when to use which resource).

Then read the nearby example README/tests and only the relevant references:

| Need | Resource |
|---|---|
| Live signatures, supported targets and generator options | [API guide](../../docs/API_USAGE.md) |
| Add an example or change its declared scope | [Example pipeline](../../HOWTO_ADD_AN_EXAMPLE.md) |
| Choose a public starting example | [Example map](references/example-map.md) |
| Choose the next executable check | [Validation ladder](references/validation-ladder.md) |
| Diagnose a specific failure | [Pitfalls](references/pitfalls.md) |
| Understand a formulation or workflow change | [Pivots](../../docs/pivots.md) |
| Run compiled elements, exact batches or PETSc/MPI | [Host API and limits](../../abaqus_ufl/fe/README.md) |

## Select the target

- Use a finite-strain UMAT for a local `stress_PK1(F)` law with built-in elements.
  `ndim=2` is an explicit plane-strain specialization; the default is full 3D.
- Use a small-strain UMAT for `SmallStrainMaterial.stress_update` and local state.
- Use a UEL for extra solved nodal fields, gradient terms or mixed interpolation.
  Choose the formulation explicitly when comparing element behavior.
- Use f2py for point/element interface checks, a small FE host for mesh/solve
  checks, and genuine Abaqus for its production call environment.

## Work from the declared theory

Start with equations, conventions and generator-friendly Python. Verify,
regenerate, compile and execute the smallest useful check before enlarging the
problem. The [API guide](../../docs/API_USAGE.md#tensor-dsl-rules) owns supported
DSL operations; the [pitfalls](references/pitfalls.md) index known traps.
Inspect unsupported syntax or a missing primitive before rewriting physics to
make generation succeed. Update the authoritative declaration/generator,
then regenerate; preserve dated submitted artifacts as evidence.

Default to one agent. Independence means a separately derived oracle or
reference solution; a team review is not a substitute. The rationale and
conditional delegation advice live in [agent workflow lessons](../../docs/lessons/agent-workflow.md).

Use a nondegenerate state that enters the claimed branch. Keep tangent
consistency, independent physical evidence, compiled-interface behavior,
solver convergence and output validation distinct. The
[validation ladder](references/validation-ladder.md) routes those gates.

## Preserve the task and public boundary

Check the worktree and the user's scope before editing. Preserve existing
API, state layout, property order, theory and unrelated changes unless the
task authorizes a change. Do not infer permission to run expensive experiments,
disclose material, send messages, merge or publish from repository access.
Existing explicit task authorization controls those actions.

The [public example map](references/example-map.md) is the capability boundary.
Use examples actually shipped in this checkout. Do not disclose private
project names, paths, models, data or provenance, or advertise an unpublished
model as public support. Preserve the paper citation and existing credits.
Report the executing revision, checks run, measured results, unrun checks and
remaining evidence limits.
