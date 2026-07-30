# Coupled Thermo-Mechanical Quad4 UEL (scalar diffusion + mechanics)

The first UEL exemplar of the example pipeline. Two nodal fields
(displacement and temperature), a cross-field coupling (thermal expansion
in the stress), and a backward-Euler transport balance enter through the
declared weak form; the generator emits a self-contained Abaqus UEL whose
seven coupled tangent blocks are evaluated by complex step.

## Theory and conventions

```text
Momentum:   P = G (F - F^{-T}) + K ln(J) F^{-T} - K alpha T F^{-T}
Transport:  storage = rho_cp (T - T_old)/dt,   flux = -k grad_T
```

Properties in Abaqus *UEL PROPERTY order: `G, K, alpha, k, rho_cp`.

Here `K` is the logarithmic volumetric (Lame-type) coefficient multiplying
`ln J`, and `alpha` is the implemented thermal-pressure coupling coefficient
in `P_th = -K alpha T F^{-T}`, not a linear-expansion coefficient.
Element Quad4, per-node DOFs `(u1, u2, T)`, `NDOFEL = 12`; temperature
rides Abaqus DOF 11 in a deck. No SVARS state. The transport equation is
declared through the recognized `solvent_storage`/`solvent_flux` material
methods (heat notation in the docstrings).

Deliberately omitted: radiation/convection boundary terms (deck-side),
temperature-dependent properties, and any dissipation coupling.

## Independent oracles

Material point (`check_reference.py`): Fourier flux `-k g`, storage
`rho_cp dT/dt` and its zero at `T = T_old`, thermal stress
`P = -K alpha T I` at `F = I`, all hand-evaluated; a wrong-formulation
control (a manually constructed incorrect expression, not an injected
code failure) rejects a sign-flipped Fourier law.

Assembled element (`check_assembled.py`), because `problem.verify()` is
not an element check:

- zero residual at the reference state;
- closed-form thermal-stress nodal forces for uniform T on the unit
  square (`RHS_mech = K alpha T [-1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2,
  1/2]`, from the shape-gradient integrals);
- rigid-translation invariance of the residual;
- heat balance by partition of unity, which kills the flux term exactly:
  zero thermal-row sum for a steady linear field on a DISTORTED element,
  and `rho_cp dT/dt x area` for a uniform temperature step;
- assembled `AMATRX` vs `-dRHS/dU` finite differences over all 12 DOFs.

## Reproduce the pipeline

```bash
python build.py            # weak-form CS-vs-FD gates + generation
python check_reference.py  # material-point closed forms + broken control
python check_assembled.py  # element invariants + assembled FD tangent
python check_compiled.py   # regeneration parity, gfortran, f2py element calls
```

## Recorded evidence

| Facet | Command | Result |
|---|---|---|
| Weak-form tangent consistency | `python build.py` | all CS-vs-FD blocks pass; operator-sign gate OK |
| Material-point oracle | `python check_reference.py` | flux/storage/thermal-stress closed forms exact |
| Broken control | `python check_reference.py` | sign-flipped Fourier law rejected |
| Assembled invariants | `python check_assembled.py` | thermal-stress nodal forces closed form; translation invariance; heat balance on distorted element |
| Assembled tangent | `python check_assembled.py` | `AMATRX` vs `-dRHS/dU`, rel err 1.0e-07 |
| Deterministic generation | `python check_compiled.py` | byte-identical regeneration |
| Fortran compile | `python check_compiled.py` | gfortran fixed-form, no errors |
| Compiled parity | `python check_compiled.py` | f2py `RHS`/`AMATRX` vs reference assembly at 3 states, 1.1e-16 / 2.9e-16 |
| Compiled tangent | `python check_compiled.py` | compiled `AMATRX` vs `-dRHS/dU`, rel err 1.0e-07 |
| Abaqus execution | included one-element deck with current generated source | Abaqus/Standard 2022 completed 10 increments with zero cutbacks, errors, numerical-problem warnings, or negative-eigenvalue warnings |

The compiled parity states include a generic deformed/heated state with a
nonzero increment on a distorted element, so the Voigt-free UEL path
(DOF maps, quadrature, coupled blocks, transport terms) is exercised end
to end in the exact Fortran that Abaqus would compile.

## Abaqus interface deck

`abaqus/job.inp` documents the intended `*USER ELEMENT` usage: DOF list
`1, 2, 11` on every node, `Unsymm` on the element card, `unsymm=YES` on
`*Step` (the complex-step tangent is the full nonsymmetric coupled
`dR/dU`), and `*UEL PROPERTY` in the declared props order. The step is a
`*Coupled Temperature-displacement` procedure: Abaqus rejects DOF-11
boundary conditions under `*Static`, so the deck mirrors the structure
of an internal case that passed Abaqus 2022 (near-zero-stiffness overlay
sharing the UEL nodes plus a small remote `CPE4T` element whose material
carries density/conductivity/specific heat to activate DOF 11).
Temperature is prescribed on the left edge only, leaving the right-edge
scalar DOFs free for a genuine transient solve. On 2026-07-30, the included
deck and current generated source completed in Abaqus/Standard 2022 with 10
increments, zero cutbacks, and no numerical-problem or negative-eigenvalue
warnings.

The generated UEL does not maintain the Abaqus `ENERGY` array, so
energy-balance output is not meaningful for this element. Property
values are illustrative in any consistent unit system; the checks in
this bundle define the tested parameter and state range. Mass,
damping, and initial-acceleration requests (`LFLAGS(3)=3,4,6`) return
zeroed arrays, and nonpositive `det F` states request a cutback.
