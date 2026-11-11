# Changelog

## Unreleased

- Fix comparisons in material methods: both operands are now emitted as real
  values. A property or a constant on the left (`if self.flag > 2.5:`), and any
  real expression on either side (`x.real > self.G * 2`), used to be emitted as
  `DCMPLX(...)` and the generated file did not compile. A complex variable on
  the right (`x.real > y`) now raises the same error as on the left.
- Fix vector state variables in `stress_PK1` of a UEL: they were registered as
  scalars there, so a helper taking `<name>_old` could not be called from both
  `stress_PK1` and a flux method. A helper may now also return a vector
  expression (`return a + c * b`), which was emitted as `(expr)(ii)`.
- Add `build_external_uel` and `CompiledAbaqusElement.from_module` for a
  standard Quad4 `UEL` compiled outside abaqus_ufl: validated arguments, the
  same adapter and drivers, and a documented TIME convention that a subclass
  may override through `begin_increment`.
- Add serial stepping drivers (`solve_steps`, `solve_steps_staggered`,
  `solve_steps_multimaterial`, `assemble_groups`) and a bounded KKT solver
  (`solve_steps_bounded`, `BoxSpec`, `CompiledElementState`) with transactional
  cutbacks; `solve_steps_staggered` also takes restart dictionaries. A step is
  committed only below `tol`, and bounded routes keep prescribed values exactly.
  Compiled elements receive beginning-of-increment `TIME`, each attempt's
  `DTIME`, and `KINC` that advances on acceptance. Nonfinite element output
  raises `NonfiniteElementOutput`, a `RuntimeError` subclass.
- Add opt-in exact compiled batch caching, batched assembly and explicit
  material-state commits for the standard Quad4 host. Cached trial outputs
  cannot alias callers, and commit/rejection invalidates reuse.
- Add optional PETSc CSR/FieldSplit solvers, distributed MPIAIJ/ghost-vector
  hosts and a replicated MPI batch adapter with collective rejection and
  unanimous cache decisions. The package imports without petsc4py or mpi4py.
- Add a two-rank public mechanics/heat example checked against serial,
  an affine patch and an independent discrete heat eigenmode.
- Add opt-in `generate_umat(..., ndim=2)` for plane strain, and automatically
  restrict standard/F-bar 2D UEL perturbations to active entries. Pure,
  straight-line state-output work may be omitted from tangent-only calls;
  stress dependencies and the full accepted-state update are preserved.
  Default 3D generated files remain byte-identical in six pinned cases.
- Measure the public Ogden Quad4 at 56.55 versus 31.59 microseconds per element
  for nine versus four perturbations (1.79x; identical R/K, compilation
  excluded). `tools/benchmark_2d_tangents.py` reproduces the measurement.
- Fix the stress-free gate to retain supplied scalar/history values or use
  an explicit/material-declared `reference_state`; report assumptions on
  failure. Public heat examples declare their cold reference separately.
- Refresh current 2D generated exports while retaining submitted snapshots.
- Establish one task index for lessons, pitfalls, pivots and skills; replace
  duplicated field notes with compatibility links and make the skill routing
  metadata. Default to one agent and an independently derived oracle. Retire
  unscoped MPI-vendor advice, obsolete power warnings and clipping-as-equilibrium
  claims.

The [host guide](abaqus_ufl/fe/README.md), [API guide](docs/API_USAGE.md) and
[task index](docs/lessons/README.md) describe supported scope and checks. Tiny
MPI examples and kernel timings do not establish production scaling or new
paper-scale Abaqus validation.
