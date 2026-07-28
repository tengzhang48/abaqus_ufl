# Small-Strain Viscoelastic SLS UMAT (history, tensor state)

The history exemplar of the example pipeline. A Standard Linear Solid
(one Maxwell branch in parallel with an equilibrium spring) whose
backward-Euler dashpot update is unconditionally stable in `dt/tau`. The
deviatoric viscous strain is a full 3x3 tensor state threaded through
`STATEV(1..9)`.

## Theory and conventions

Compression-positive small-strain update with properties in Abaqus
*User Material order `K, G_inf, G_v, tau`:

```text
eps_v_new = (eps_v_old + (dt/tau) dev(eps)) / (1 + dt/tau)
sigma     = K tr(eps) I + 2 G_inf dev(eps) + 2 G_v (dev(eps) - eps_v_new)
```

`STATEV(1..9)` stores `eps_v` column-major (slot `3(j-1)+i` for component
`(i,j)`), in the USER compression-positive convention. The generated UMAT
flips signs only at the `STRESS`/`STRAN`/`DSTRAN` boundary, never inside
the stored state, so a positive Abaqus shear strain produces a negative
stored `eps_v` component. The compiled gate asserts this documented
behavior explicitly.

The generated UMAT is strictly three-dimensional (`NDI=3, NSHR=3,
NTENS=6` required); plane-stress, plane-strain, and axisymmetric use
are unsupported.

Abaqus energy outputs (`SSE`, `SPD`, `SCD`) are not maintained by the
generated UMAT, so Abaqus energy-balance output is not meaningful for
this material. Property values are illustrative in any consistent
unit system; the checks in this bundle define the tested parameter
and state range.

Deliberately omitted: multiple Maxwell branches, volumetric relaxation,
temperature dependence (thermorheological simplicity), and finite strain.

## Independent oracle

Step-shear relaxation has an exact hand-derived closed form for the
DISCRETE backward-Euler recursion, with `r = dt/tau`:

```text
tau_n = 2 e [ G_inf + G_v / (1+r)^n ]
eps_v_xy(n) = e [ 1 - 1/(1+r)^n ]
```

so the response relaxes geometrically from the instantaneous limit
`2e(G_inf+G_v)` toward the equilibrium `2e G_inf`, and `(1+r)^{-n}`
converges to the continuous `exp(-t/tau)` as `dt` shrinks (checked at two
`dt` levels). A broken control confirms an explicit (forward-Euler-like)
dashpot update is rejected.

## Reproduce the pipeline

```bash
python build.py            # Python tangent consistency + generation
python check_reference.py  # discrete relaxation closed form + limits + control
python check_compiled.py   # regeneration parity, gfortran, f2py history drive
```

## Recorded evidence

| Facet | Command | Result |
|---|---|---|
| Python tangent consistency | `python build.py` | CS vs FD err 1.4e-07 (field-argument blocks; state/dt evidence from the compiled history gates) |
| Independent physics oracle | `python check_reference.py` | discrete relaxation closed form at every increment, err < 1e-13 |
| Physical limits | `python check_reference.py` | instantaneous/equilibrium bracket; >95% of relaxable stress shed |
| dt-convergence | `python check_reference.py` | backward-Euler factor converges to `exp(-t/tau)` |
| Broken control | `python check_reference.py` | explicit-dashpot update rejected |
| Deterministic generation | `python check_compiled.py` | byte-identical regeneration |
| Fortran compile | `python check_compiled.py` | gfortran fixed-form, no errors |
| Direct compiled execution | `python check_compiled.py` | 30-increment relaxation history, stress err 5.6e-17 |
| Tensor `STATEV` layout | `python check_compiled.py` | column-major slots, symmetric shear pair, empty diagonal, user-convention sign, err 5.2e-18 |
| Abaqus execution | not included in this bundle | none claimed |
