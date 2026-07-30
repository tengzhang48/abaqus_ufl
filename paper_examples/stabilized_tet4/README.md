# Stabilized mixed u-theta tetrahedron (Sec. 3.2, Fig. 4)

Kinematically stabilized equal-order Tet4 of Scovazzi, Zorrilla, and
Rossi (2023), block-compression benchmark.

- `scovazzi_block.py` — declaration and generation pipeline
  (`python scovazzi_block.py --element tet4`).
- `abaqus/` — the n=16 quarter-block deck, its deck generator
  (`generate_inp.py`), the ODB extraction and render scripts, and the
  reduced results (center-displacement history and summary JSON) from
  the completed Abaqus run.
- `figure/` — the Fig. 4 script with the exact Abaqus renderings
  (`Tet4_u_mag.png`, 0-0.6962 mm native legend; `Tet4_NT11.png`,
  theta-tilde in [-7.5, +1.0]e-4).

The published-curve comparison value is u3 = -0.696244 mm at full
follower pressure, within about 1% of the published n=16 curve.

## Fresh-clone Abaqus check

On 2026-07-30, the current source regenerated and compiled. The tracked `n=16`
deck completed an Abaqus/Standard 2022 datacheck with zero errors. A generated
`n=2` smoke model (86 nodes and 282 mixed Tet4 elements) then completed all 24
increments with zero cutbacks. The full `n=16` analysis was not repeated.
