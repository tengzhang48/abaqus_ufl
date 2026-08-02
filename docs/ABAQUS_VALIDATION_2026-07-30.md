# Public Abaqus Validation, 2026-07-30

This report records a fresh-clone validation of the public repository. The
clone started at revision `0f52533` on `main`; the corrections found during
the pass and this report were committed in `c5d2a59`. Solver output was
written to temporary compute-node storage and is not retained in Git.

## Environment

- Abaqus/Standard 2022
- Intel Fortran 19.1.1.217, selected by the Abaqus installation
- Python 3.11.15
- NumPy 2.4.2
- SymPy 1.14.0
- pytest 9.0.3

The conda environment required the standard-library distutils implementation
for its NumPy f2py checks:

```bash
SETUPTOOLS_USE_DISTUTILS=stdlib python -m pytest -q
```

Without that environment setting, the Python and reference checks still
passed, but f2py skipped after a setuptools/distutils constructor conflict.
This was an environment compatibility issue rather than a generated-Fortran
failure.

## Python and Compiled Checks

The complete public test suite passed with f2py enabled: `137 passed`. It
covered all six released examples through deterministic generation, reference
and assembled oracles, gfortran compilation, and direct calls to the compiled
UMATs or UELs.

All four manuscript generation entry points were also exercised:

| Package | Generation result |
|---|---|
| Phase-field corrosion | Material/tangent verification passed; main and diagonal sources regenerated |
| Stabilized Tet4 | Three verification states passed; source regenerated |
| Mixed-order gel bilayer | Material/tangent verification and UEL generation passed; deck conversion skipped because the separately distributed mesh seed was absent |
| Local-pressure Hex8 morphing | Material/tangent verification passed; current source regenerated |

All freshly generated manuscript Fortran sources passed a fixed-form gfortran
syntax compilation.

## Abaqus Results

Small cases were run to completion. Large or expensive paper cases were
limited to source compilation and Abaqus datacheck.

| Case | Scope | Result |
|---|---|---|
| Working UMAT template | Full one-element solve | 10 increments, zero cutbacks, successful completion; all 11 ODB/reference checks passed |
| Scalar-diffusion Quad4 UEL | Full one-element solve | 10 increments, zero cutbacks, successful completion |
| Mixed-order thermo-mechanical Quad8 UEL | Full one-element solve | 10 increments, zero cutbacks, successful completion |
| Mixed-order gel Quad8 | Full one-element pressure-equilibrium smoke solve | 1 increment, zero cutbacks, successful completion |
| Mixed-order gel Quad8 | Full two-element chemical-potential-gradient smoke solve | 1 increment, zero cutbacks, successful completion |
| Stabilized Tet4 | Generated `n=2` full smoke solve | 86 nodes, 282 mixed Tet4 elements, 24 increments, zero cutbacks, successful completion |
| Stabilized Tet4 Figure 14 | Tracked `n=16` deck datacheck | Complete, 17 `.dat` warnings, zero errors |
| Cui corrosion Figure 3 | Tracked full deck with freshly generated diagonal UEL, datacheck only | Complete, 25 `.dat` and 1 `.msg` warnings, zero errors |
| Pasta morphing Figure 6 | Exact tracked production deck and archived submitted UEL, datacheck only | Complete, 13 `.dat` and 1 `.msg` warnings, zero errors |

The completed small analyses reported zero numerical-problem messages and zero
negative-eigenvalue messages. Their warnings were expected integration-deck
messages such as unsupported direct element output for UELs, the remote
thermal activation element forming a second unconnected region, default
absolute-zero handling, and zero force or heat flux in deliberately
equilibrated smoke states.

## Defects Found and Corrected

### Corrosion generation layout

The corrosion postprocessor searched for the old unguarded `SVARS` write
layout. The current generator commits state only for normal
`LFLAGS(3)=1` calls, so the exact marker no longer matched. Generation stopped
after otherwise successful material and tangent checks.

The bridge insertion now matches the final state assignment itself and remains
inside the current state-commit guard. The default main and diagonal outputs
now go to the public `generated/` directory; stale copying to the removed
lab-only `abaqus_test_from_cui/` directory was removed.

### Gel deck regeneration boundary

The gel builder successfully regenerated the UEL and property include, then
failed while opening the non-redistributed Chester-Anand mesh seed. The
included accepted Abaqus deck does not require that seed.

The builder now exits successfully after UEL/property generation when the seed
is absent, reports exactly what was skipped, and points to the accepted deck.
If the seed is supplied at the documented location, deck conversion still
runs.

Regression tests cover both corrections.

## Scope Boundaries

- The full corrosion, `n=16` Tet4, gel-bilayer, and pasta analyses were not
  repeated in this pass.
- The current Hex8 generator source is not interchangeable with the archived
  pasta source because the production dummy-element labels use a different
  UVARM offset. The exact archived source was therefore used for the production
  deck datacheck.
- A successful datacheck verifies input processing, user-subroutine
  compilation/linking, procedure compatibility, and model setup. It does not
  establish nonlinear convergence or reproduce a published result.
- Abaqus result files and transient compiler/solver artifacts remain outside
  the repository.
