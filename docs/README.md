# Documentation

Start with the top-level [README](../README.md) for install + quickstart, and
[HOWTO_ADD_AN_EXAMPLE](../HOWTO_ADD_AN_EXAMPLE.md) for the complete
theory-to-compiled-code-to-output example pipeline.

## Usage

- [API_USAGE.md](API_USAGE.md) — import patterns, the generator entry points
  (`generate_umat` / `generate_uel`), the tensor DSL rules, and the
  verification ladder.

## Theory

- [theory.md](theory.md) — the mathematical framework: weak form and
  discretization, complex-step tangents, mixed interpolation.
- [complex_step_patterns.md](complex_step_patterns.md) — how to write
  complex-step-safe constitutive code (branching, `.real` leaks, iterative maps).
- [JAUMANN_RESOLUTION.md](JAUMANN_RESOLUTION.md) — the Jaumann-rate tangent used
  for the UMAT `DDSDDE`, with derivation and numerical verification.

## Design

- [API_EXTENSIBILITY_PLAN.md](API_EXTENSIBILITY_PLAN.md) — the DSL design
  philosophy and the extensibility model.
- [uel_design.md](uel_design.md) — the UEL code-generator design (templates,
  material translation, CS tangent engine, UEL assembly).
- [SYMBOLIC_TANGENT_DESIGN.md](SYMBOLIC_TANGENT_DESIGN.md) — the symbolic (SymPy)
  tangent engine that emits pure-real Fortran tangents.
- [matrix_functions_design.md](matrix_functions_design.md) — **complex-step-safe
  matrix log / exp / sqrt and eigenvalue routines** (iterative vs.
  eigendecomposition backends, Cardano, degeneracy handling). These underpin
  finite-strain and plasticity models.
- [internal_variables_design.md](internal_variables_design.md) — how state /
  internal variables are declared, ordered, and threaded through the generator.

## Lessons Learned

General, hard-won lessons for Abaqus UMAT/UEL development and code generation —
see the [lessons/](lessons/) index (Abaqus/Fortran traps, finite-element
pitfalls, complex-step and matrix functions, code generation and testing,
coupled multiphysics, and field notes).

For an agent-oriented operational guide, see [ai_skills/](../ai_skills/).
