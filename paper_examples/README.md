# Paper example packages

These four packages collect materials associated with the examples in the
manuscript *Making coupled-field Abaqus user elements simple*. Across the
packages, those materials include Python declarations, generation entry
points, generated or submitted Fortran, Abaqus decks, reduced reference data,
and figure inputs. They are evidence bundles, not collectively clean-archive
end-to-end reproduction packages.

| Package | Manuscript section / figure | Pipeline command |
|---|---|---|
| `phasefield_corrosion/` | Sec. 3.1, Fig. 3 | `python phasefield_corrosion_cui.py` (main), `python -c "import phasefield_corrosion_cui as m; m.generate_diag()"` (diagonal comparison variant) |
| `stabilized_tet4/` | Sec. 3.2, Fig. 4 | `python scovazzi_block.py --element tet4` |
| `gel_bilayer/` | Sec. 3.3, Fig. 5 | `python build.py` |
| `morphing_hex8/` | Sec. 3.4, Fig. 6 | `python pipeline_hex8.py` |

The commands above are the intended generation entry points and should be run
inside each package directory. Current portability limitations are recorded
below.

## Current portability and provenance status

The July 30 fresh-clone validation found:

- all four declaration-to-Fortran entry points now complete in the public
  layout; when the non-redistributed gel mesh seed is absent, `build.py`
  regenerates the UEL and properties and explicitly skips only deck
  conversion;
- all freshly generated manuscript Fortran sources compile with gfortran;
- the current corrosion source passes an Abaqus/Standard 2022 datacheck with
  the tracked full deck, the current Tet4 source completes a generated `n=2`
  Abaqus smoke solve, and the current mixed-order gel source completes both
  tracked one- and two-element smoke solves;
- the exact archived pasta source and production deck pass a fresh Abaqus 2022
  datacheck;
- all four figure scripts still retain assumptions from the original
  manuscript-production layout and were not validated as standalone
  public-archive entry points; and
- the included corrosion raw/reduced inputs still do not reconstruct the
  complete Figure 3 comparison.

See `../docs/ABAQUS_VALIDATION_2026-07-30.md` for the environment, exact
scope, solver outcomes, expected warnings, and remaining boundaries.

The exact submitted sources and retained reduced results remain useful
provenance. Do not describe the four packages collectively as cleanly
reproducible until these package-specific gaps are closed and tested from a
fresh archive.

## Exact artifacts vs current pipeline

Where a package contains an `archived_submitted/` folder, those files
are the **exact sources submitted with the completed Abaqus runs**,
preserved byte-for-byte for provenance. The current pipeline uses today's
generator, which has since gained procedure-dispatch and invalid-state guards,
to emit sources for the same declared formulations;
regenerated sources are therefore not byte-identical to the archived
ones, and both are kept.

## Naming note

The archived morphing source and some symbols carry a historical
internal prefix derived from Chester, Di Leo, and Anand, whose gel
theory and finite-element treatment the model builds on. The
formulation used here is the pressure-based adaptation described in the
paper (and in its Supplementary Material), **not** the original
Chester--Di Leo--Anand UEL; the current pipeline uses the neutral
`pressuregel` prefix.

## Third-party material

The corrosion package compares against the reference implementation of
Cui, Ma, and Martinez-Paneda (JMPS 2021). The authors' UEL header says the
code is distributed under a BSD license, but the audited distributions do not
include the exact BSD text or copyright notice. Our comparison deck has mesh
lineage from that distribution (with attribution); the reference UEL and deck
themselves are **not** redistributed here. Obtain them from the authors'
repository and verify against the SHA-256 hashes and open provenance boundary
recorded in the package README.

## Figure kits

Each `figure/` folder holds the retained manuscript figure script and selected
image/data inputs. The scripts were written for the manuscript repository
layout and still require portability repair. The corrosion package also needs
additional source data before its complete published comparison can be
reconstructed from this checkout.

These packages are excluded from the Python source distribution; they live in
the repository only. A public tag and archive DOI remain pending.
