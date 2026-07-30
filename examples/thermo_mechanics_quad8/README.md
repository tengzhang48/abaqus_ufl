# Mixed-Order Thermo-Mechanical Quad8 UEL (quadratic u, bilinear T)

The mixed-order exemplar of the example pipeline. Displacement is
quadratic serendipity on all eight nodes while the temperature-like
scalar is bilinear on the four corners, so the generated element carries
node-dependent DOF maps: corners `(u1, u2, T)` (DOFs 1..12), midsides
`(u1, u2)` (DOFs 13..20), `NDOFEL = 20`. The scalar flux is pulled back
with `C^{-1}`, so the transport block is genuinely deformation
dependent.

## Theory and conventions

```text
Momentum:   P = G (F - F^{-T}) + K ln(J) F^{-T} - K alpha T F^{-T}
Transport:  storage = cT (T - T_old)/dt,   flux = -kappa C^{-1} grad_T
```

Properties in Abaqus *UEL PROPERTY order: `G, K, alpha, kappa, cT`.

Here `K` is the logarithmic volumetric (Lame-type) coefficient multiplying
`ln J`, and `alpha` is the implemented thermal-pressure coupling coefficient
in `P_th = -K alpha T F^{-T}`, not a linear-expansion coefficient.
No SVARS state. In a deck the corner temperature rides Abaqus DOF 11.

Deliberately omitted: boundary flux/convection terms (deck-side),
temperature-dependent properties, dissipation coupling.

## Independent oracles

Material point (`check_reference.py`): the pulled-back flux
`-kappa a / l^2` under uniaxial stretch (the `1/l^2` factor is the
point), storage forms, thermal stress `-K alpha T I` at `F = I`, and a
wrong-formulation control (a manually constructed incorrect
expression, not an injected code failure) rejecting a missing
`C^{-1}` pull-back.

Assembled element (`check_assembled.py`):

- zero reference residual and residual invariance under rigid
  translation of a heated state;
- a pulled-back-flux ELEMENT oracle: an exact homogeneous x-stretch `l`
  (representable by the quadratic basis) with a steady corner-linear
  field `T = aX` gives thermal rows `-(kappa a / l^2) g_a` from the
  bilinear shape-gradient integrals, and the `l`-run must differ from
  the undeformed run by exactly `1/l^2`;
- partition-of-unity heat balance on a distorted Quad8 (zero flux-term
  sum steady; exact `cT dT/dt x area` for a uniform step);
- assembled `AMATRX` vs `-dRHS/dU` over all 20 DOFs.

## Reproduce the pipeline

```bash
python build.py            # weak-form CS-vs-FD gates + generation
python check_reference.py  # material-point closed forms + broken control
python check_assembled.py  # mixed-order element oracles + FD tangent
python check_compiled.py   # regeneration parity, gfortran, f2py element calls
```

## Recorded evidence

| Facet | Command | Result |
|---|---|---|
| Weak-form tangent consistency | `python build.py` | all CS-vs-FD blocks pass; operator-sign gate OK |
| Material-point oracle | `python check_reference.py` | pull-back flux, storage, thermal stress closed forms exact |
| Broken control | `python check_reference.py` | missing `C^{-1}` pull-back rejected |
| Assembled invariants | `python check_assembled.py` | signed thermal rows `-(kappa a/l^2) g_a`; `1/l^2` ratio exact; heat balance on distorted element |
| Assembled tangent | `python check_assembled.py` | `AMATRX` vs `-dRHS/dU`, rel err 1.5e-07 |
| Deterministic generation | `python check_compiled.py` | byte-identical regeneration |
| Fortran compile | `python check_compiled.py` | gfortran fixed-form, no errors |
| Compiled parity | `python check_compiled.py` | f2py `RHS`/`AMATRX` vs reference assembly at 3 states, 7.6e-15 / 2.6e-16 |
| Compiled tangent | `python check_compiled.py` | compiled `AMATRX` vs `-dRHS/dU`, rel err 1.5e-07 |
| Abaqus execution | included one-element deck with current generated source | Abaqus/Standard 2022 completed 10 increments with zero cutbacks, errors, numerical-problem warnings, or negative-eigenvalue warnings |

The mixed-order DOF maps are exercised end to end: the same interleaved
layout must be produced by the Python reference assembly and consumed by
the generated Fortran for the machine-precision parity above to hold.

## Abaqus interface deck

`abaqus/job.inp` documents the intended `*USER ELEMENT` usage for the
mixed-order layout: data line `1, 2, 11` for corner positions 1-4 and a
second data line `5, 1, 2` so midside positions 5-8 carry displacement
only, with `Unsymm` on the element card and `unsymm=YES` on `*Step`.
The step is a `*Coupled Temperature-displacement` procedure: Abaqus
rejects DOF-11 boundary conditions under `*Static`, so the deck mirrors
the structure of an internal case that passed Abaqus 2022 with an
earlier generated UEL (near-zero-stiffness overlay plus a small remote
`CPE8T` element to activate DOF 11). Temperature is prescribed on the
two bottom corners only, leaving the top-corner scalar DOFs free for a
genuine transient solve. On 2026-07-30, the included deck and current
generated source completed in Abaqus/Standard 2022 with 10 increments,
zero cutbacks, and no numerical-problem or negative-eigenvalue warnings.

The generated UEL does not maintain the Abaqus `ENERGY` array, so
energy-balance output is not meaningful for this element. Property
values are illustrative in any consistent unit system; the checks in
this bundle define the tested parameter and state range. Mass,
damping, and initial-acceleration requests (`LFLAGS(3)=3,4,6`) return
zeroed arrays, and nonpositive `det F` states request a cutback.
