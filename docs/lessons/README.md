# Lessons Learned

Hard-won, general lessons for developing Abaqus **UMAT / UEL** subroutines and
generating them from a Python DSL. Distilled from building and validating a range
of constitutive and coupled-field models — kept general so they transfer to any
Abaqus user-subroutine project.

- [abaqus-and-fortran.md](abaqus-and-fortran.md) — Abaqus UEL/UMAT interface
  traps, fixed-form Fortran compilation, and Abaqus input-deck compatibility.
- [finite-element.md](finite-element.md) — plane-strain `F₃₃`, Jacobian mapping
  for mixed-degree fields, element-local DOF indexing, F-bar, and the
  Jaumann-rate tangent trap.
- [complex-step-and-matrix-functions.md](complex-step-and-matrix-functions.md) —
  complex-step differentiation, and **CS-safe matrix log / exp / sqrt and
  eigenvalue routines** (non-holomorphic eigenvectors, Cardano branch cuts,
  degeneracy fallback, polar decomposition via `inv(U)`).
- [codegen-and-testing.md](codegen-and-testing.md) — generating correct Fortran,
  verifying it (the finite-difference tangent check), and DSL design principles.
- [coupled-multiphysics.md](coupled-multiphysics.md) — scaling and conditioning,
  volumetric locking and element choice, operator signs in coupled weak forms,
  and reproducing multiphysics papers.
- [field-notes.md](field-notes.md) — shorter general lessons on solvers,
  performance, distributed/MPI runs, verification discipline, and reproduction.
