# Mixed-order gel bilayer (Sec. 3.3, Fig. 5)

Plane-strain swelling of a gel-rubber bilayer with the global
u-p-mu mixed-order Quad8 element (quadratic displacement and chemical
potential, bilinear corner pressure; 28 DOFs; node-dependent Abaqus DOF
maps).

- `build.py` — declaration, generation, and deterministic deck builder
  (`python build.py` always regenerates the UEL and regenerates the bilayer
  deck when the separately distributed mesh seed is available).
- `abaqus/` — the accepted bilayer deck, configured for a final time of
  21600 s and a 50000-increment ceiling, plus the one- and two-element
  smoke decks and shared property include.
- `figure/` — the Fig. 5 script with the exact Abaqus/CAE exports and
  legends at t = 0, 30 min, 1 h, 6 h.

The historical project record reported final right-edge straightness of
2.3e-9 m. This checkout does not retain the full-run solver log or a
machine-readable result record supporting that value, so it is not presented
as retained public evidence.

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
not require that external mesh seed. Without it, `python build.py` verifies and
regenerates the UEL and property include, reports that deck conversion was
skipped, and exits successfully.

## Fresh-clone Abaqus check

On 2026-07-30, the current source regenerated and compiled. The included
one-element pressure-equilibrium and two-element chemical-potential-gradient
decks both completed in Abaqus/Standard 2022 with zero cutbacks, errors,
numerical-problem warnings, or negative-eigenvalue warnings. The full bilayer
analysis was not repeated.
