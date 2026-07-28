# Paper example packages

These four packages collect materials associated with the examples in the
manuscript *Making coupled-field Abaqus user elements simple*. Across the
packages, those materials include Python declarations, generation entry
points, generated or submitted Fortran, Abaqus decks, reduced reference data,
and figure inputs. At public revision `e0985d9`, they are evidence bundles,
not yet clean-archive end-to-end reproduction packages.

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

The July 28 clean-archive audit found:

- all four figure scripts retain paths from the manuscript/lab layout and do
  not run unmodified from the public archive;
- the corrosion pipeline generates source but then fails on a missing lab
  destination, and the included raw/reduced inputs do not reconstruct the
  complete Figure 3 comparison;
- the included accepted bilayer deck is this project's deck for the
  swell-induced bending problem of Chester, Di Leo, and Anand, and its mesh
  discretization follows their supplemental example; only rerunning
  `build.py` requires their separately obtained input as the mesh seed;
- at revision `a1b5825`, the morphing run instructions referenced generated
  source outside the stated working directory; the current package README
  corrects that path; and
- the public suite at `e0985d9` reports `135 passed` and includes the complete
  current-generator Quad4/Hex8 local-pressure contract and condensed-Jacobian
  tests, but it does not execute the `paper_examples/` pipelines, figure
  scripts, or exact archived production sources.

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
Cui, Ma, and Martinez-Paneda (JMPS 2021), distributed by its authors
under a BSD license. Our comparison decks derive their mesh from that
distribution (with attribution); the reference UEL and deck themselves
are **not** redistributed here — obtain them from the authors'
repository and verify against the SHA-256 hashes recorded in the
package README.

## Figure kits

Each `figure/` folder holds the retained manuscript figure script and selected
image/data inputs. The scripts were written for the manuscript repository
layout and still require portability repair. The corrosion package also needs
additional source data before its complete published comparison can be
reconstructed from this checkout.

These packages are excluded from the Python source distribution; they
live in the repository only. A public tag and archive DOI were still pending
at revision `e0985d9`.
