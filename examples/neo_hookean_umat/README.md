# Compressible Neo-Hookean UMAT (finite strain, stateless)

The finite-strain exemplar of the example pipeline. One hyperelastic law,
no state variables, generated as a self-contained Abaqus/Standard UMAT.

## Theory and conventions

First Piola-Kirchhoff stress, with shear modulus `G` and bulk-like modulus
`K` (Abaqus *User Material constants in this order):

```text
P = G (F - F^{-T}) + K ln(J) F^{-T},        J = det F
```

The generated UMAT receives `DFGRD1`, evaluates `P` and its complex-step
derivative `dP/dF`, and converts to Cauchy stress plus the Jaumann-rate
`DDSDDE` expected by Abaqus/Standard solid elements. Voigt order is the
Abaqus order (11, 22, 33, 12, 13, 23). There are no `STATEV` slots.

The generated UMAT is strictly three-dimensional: the wrapper requires
`NDI=3, NSHR=3, NTENS=6` and returns early with a `PNEWDT` cutback
otherwise, so plane-stress, plane-strain, and axisymmetric use are
unsupported.

Abaqus energy outputs (`SSE`, `SPD`, `SCD`) are not maintained by the
generated UMAT, so Abaqus energy-balance output is not meaningful for
this material. Property values are illustrative in any consistent
unit system; the checks in this bundle define the tested parameter
and state range.

Deliberately omitted: no viscous or rate effects, no anisotropy, and no
near-incompressible treatment (use a mixed or F-bar formulation when
`K/G` is large and locking matters).

## Independent oracle

`check_reference.py` evaluates two hand-derived closed forms without
calling the implemented model:

```text
uniaxial  F = diag(1.2, 1, 1):
  sigma_xx = G (lam - 1/lam) + K ln(lam)/lam
  sigma_yy = sigma_zz = K ln(lam)/lam

simple shear  F = I + 0.3 e_x e_y:
  sigma_xx = G gamma^2,  sigma_xy = G gamma   (normal-stress effect)
```

A broken control re-derives the uniaxial state with the volumetric term
deleted and asserts the oracle rejects it.

## Reproduce the pipeline

```bash
python build.py            # Python tangent consistency + generation
python check_reference.py  # closed-form oracle + broken control
python check_compiled.py   # regeneration parity, gfortran, f2py execution
```

## Recorded evidence

| Facet | Command | Result |
|---|---|---|
| Python tangent consistency | `python build.py` | CS vs FD err 4.1e-08; major symmetry 1.9e-16 |
| Independent physics oracle | `python check_reference.py` | uniaxial + simple shear closed forms, exact to 1e-12 |
| Broken control | `python check_reference.py` | volumetric-term deletion rejected |
| Deterministic generation | `python check_compiled.py` | byte-identical regeneration |
| Fortran compile | `python check_compiled.py` | gfortran fixed-form, no errors |
| Direct compiled execution | `python check_compiled.py` | f2py stress vs closed forms, max abs err 1.8e-15 |
| Compiled tangent | `python check_compiled.py` | `DDSDDE` vs Jaumann-corrected FD, rel err 1.3e-06 |
| Abaqus execution | not included in this bundle | none claimed |

No Abaqus run ships with this example; solver-level evidence is a
user-side step. The compiled f2py call exercises the exact generated
subroutine that Abaqus would compile.
