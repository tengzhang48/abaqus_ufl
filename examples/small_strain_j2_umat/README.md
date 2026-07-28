# Small-Strain J2 Plasticity UMAT (stateful, branching)

The stateful exemplar of the example pipeline. Radial-return J2
plasticity with linear isotropic hardening, one state variable, and a
genuine elastic/plastic branch.

## Theory and conventions

Compression-positive small-strain update with properties in Abaqus
*User Material order `G, lam, sigma_y, H`:

```text
sigma_trial = sigma_old + 2 G de + lam tr(de) I
f = q(sigma_trial) - (sigma_y + H ep_old),   q = sqrt(3/2 s:s)
f > 0:  dgamma = f / (3G + H),  n = 3 s / (2 q)
        sigma = sigma_trial - 2 G dgamma n,  ep = ep_old + dgamma
```

`STATEV(1)` holds the equivalent plastic strain `ep`. The declared zero
initializer fires only when the analysis starts (`TIME(2)==0`) with an
all-zero incoming `STATEV`; nonzero state supplied through `SDVINI` or
`*INITIAL CONDITIONS, TYPE=SOLUTION` is preserved. The generated UMAT converts between the
compression-positive tensor API and Abaqus' tension-positive Voigt
arrays, including the engineering-shear convention of `DSTRAN(4..6)`.

Deliberately omitted: kinematic hardening, rate dependence, thermal
effects, and finite rotation (small-strain theory only). The generated
UMAT is strictly three-dimensional (`NDI=3, NSHR=3, NTENS=6` required);
plane-stress, plane-strain, and axisymmetric use are unsupported.

Abaqus energy outputs (`SSE`, `SPD`, `SCD`) are not maintained by the
generated UMAT, so Abaqus energy-balance output is not meaningful for
this material. Property values are illustrative in any consistent
unit system; the checks in this bundle define the tested parameter
and state range.

## Independent oracle

Monotonic pure shear has an exact hand-derived solution because radial
return is exact under proportional deviatoric loading:

```text
e_y = sigma_y / (2 sqrt(3) G)
ep(e)  = (2 sqrt(3) G e - sigma_y) / (3G + H)     for e > e_y
tau(e) = (sigma_y + H ep) / sqrt(3)
```

`check_reference.py` drives the Python model through the elastic and
plastic branches and compares against these formulas; a broken control
confirms a perfect-plasticity mistake (H dropped from the return
denominator) is rejected.

## Reproduce the pipeline

```bash
python build.py            # Python tangent consistency + generation
python check_reference.py  # closed-form shear oracle + broken control
python check_compiled.py   # regeneration parity, gfortran, f2py path drive
```

## Recorded evidence

| Facet | Command | Result |
|---|---|---|
| Python tangent consistency | `python build.py` | CS vs FD err 5.9e-08 (field-argument blocks; state/dt evidence from the compiled history gates) |
| Independent physics oracle | `python check_reference.py` | elastic + plastic closed forms, err < 1e-10 |
| Regime coverage | `python check_reference.py` | both branches entered; plastic run yields strongly |
| Broken control | `python check_reference.py` | H-free return mapping rejected |
| Deterministic generation | `python check_compiled.py` | byte-identical regeneration |
| Fortran compile | `python check_compiled.py` | gfortran fixed-form, no errors |
| Direct compiled execution | `python check_compiled.py` | 40-increment shear path, stress vs closed form 5.6e-17 |
| `STATEV` round trip | `python check_compiled.py` | `STATEV(1)` = closed-form `ep` at every increment, err 1.4e-17 |
| Boundary conventions | `python check_compiled.py` | engineering shear `DSTRAN(4)`; elastic `DDSDDE(4,4)=G`; softened plastic tangent |
| Abaqus execution | not included in this bundle | none claimed |

The compiled path check doubles as a Voigt-boundary check. The Python
oracle uses tensor shear while the f2py drive passes engineering shear,
so agreement requires the generated conversion to be correct.
