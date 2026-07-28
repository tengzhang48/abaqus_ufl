# Mixed-order gel bilayer (Sec. 3.3, Fig. 5)

Plane-strain swelling of a gel-rubber bilayer with the global
u-p-mu mixed-order Quad8 element (quadratic displacement and chemical
potential, bilinear corner pressure; 28 DOFs; node-dependent Abaqus DOF
maps).

- `build.py` — declaration, generation, and deterministic deck builder
  (`python build.py` regenerates the UEL and the bilayer deck).
- `abaqus/` — the accepted bilayer deck (2160 increments to 21600 s,
  final right-edge straightness 2.3e-9 m), the one- and two-element
  smoke decks, and the shared property include.
- `figure/` — the Fig. 5 script with the exact Abaqus/CAE exports and
  legends at t = 0, 30 min, 1 h, 6 h.

This example is an execution demonstration (see the paper's Table S1);
the formulation is the pressure-based adaptation described in the
paper's Supplementary Material.

## Regenerating the deck

`build.py` derives the bilayer mesh from the original Chester--Di
Leo--Anand supplementary input (IJSS 2011 supplementary material),
which is not redistributed here. To regenerate the deck, obtain their
`gelUEL` supplementary folder and place it at
`../gelUEL/gelUEL/code/` relative to this package; the UEL generation
itself and the included accepted deck do not require it.
