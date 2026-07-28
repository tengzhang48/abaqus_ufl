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

This example is an execution demonstration; its evidence level is stated in
the manuscript Supplement. The formulation is the pressure-based adaptation
described there.

## Regenerating the deck

The accepted input deck is our deck for the swell-induced bending problem of
Chester, Di Leo, and Anand. Its mesh discretization follows their supplemental
example, with attribution; their original supplemental files are not
redistributed here. To regenerate the deck with `build.py`, obtain their
`gelUEL` supplemental folder and place it at `../gelUEL/gelUEL/code/` relative
to this package. The UEL generation itself and the included accepted deck do
not require that external mesh seed.
