# Paper examples — exact reproduction packages

These four packages reproduce the examples of the manuscript *Making
coupled-field Abaqus user elements simple*. Each package contains the
**pipeline** (the Python declaration and the command that generates the
UEL), the **generated Fortran**, the **Abaqus input decks**, the
**reduced reference data**, and the **figure kit** (the manuscript
figure script with its exact inputs).

| Package | Manuscript section / figure | Pipeline command |
|---|---|---|
| `phasefield_corrosion/` | Sec. 3.1, Fig. 3 | `python phasefield_corrosion_cui.py` (main), `python -c "import phasefield_corrosion_cui as m; m.generate_diag()"` (diagonal comparison variant) |
| `stabilized_tet4/` | Sec. 3.2, Fig. 4 | `python scovazzi_block.py --element tet4` |
| `gel_bilayer/` | Sec. 3.3, Fig. 5 | `python build.py` |
| `morphing_hex8/` | Sec. 3.4, Fig. 6 | `python pipeline_hex8.py` |

Run each pipeline from inside its package directory; scripts locate the
`abaqus_ufl` package two levels up (this repository), so no
installation is required beyond the repository checkout and its Python
dependencies.

## Exact artifacts vs current pipeline

Where a package contains an `archived_submitted/` folder, those files
are the **exact sources submitted with the completed Abaqus runs**,
preserved byte-for-byte for provenance. The current pipeline
regenerates functionally equivalent elements with today's generator
(which has since gained procedure-dispatch and invalid-state guards);
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

Each `figure/` folder holds the manuscript figure script and its exact
image/data inputs. The scripts were written for the manuscript
repository layout; their input files are included here so the mapping
from data to published figure is complete and auditable, and any path
adjustments needed to run them elsewhere are limited to the constants
at the top of each script.

These packages are excluded from the Python source distribution; they
live in the repository (and its tagged, DOI-archived releases) only.
