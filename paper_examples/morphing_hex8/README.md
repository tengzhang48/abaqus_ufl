# Three-dimensional diffusion-driven morphing (Sec. 3.4, Fig. 6)

Grooved-sheet morphing with the Hex8 element-local-pressure gel UEL
(global u and mu; one condensed pressure per element).

- `pipeline_hex8.py` — the CURRENT declaration-to-Fortran pipeline
  (verify + generate with the neutral `pressuregel` prefix).
- `archived_submitted/chester_anand_local_pressure_hex8_pasta.for` —
  the EXACT generated UEL submitted with the completed Abaqus/Standard
  2022 production run, preserved unmodified (historical internal
  prefix; see the top-level naming note).
- `abaqus/` — the exact production deck, the C3D8 visualization-mesh
  include, the material properties, and the run README copied
  unmodified from the completed run.
- `reference_data/` — the reduced ODB record of the completed run
  (frame times, final fields, and maximum displacement).
- `figure/` — the Fig. 6 script and the six Abaqus/CAE exports.

The condensed-tangent verification of this element formulation (black
box, pressure re-solved per perturbation, relative error 6.5e-10 with a
frozen-Schur negative control at 0.76) ships with the repository test
suite.
