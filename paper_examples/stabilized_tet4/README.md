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
