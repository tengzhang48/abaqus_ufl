# One-Term Compressible Ogden UMAT (explicit spectral path)

The spectral exemplar of the example pipeline. The strain energy lives in
principal stretches, so the material calls `eig` directly and
reconstructs stress with the eigenspace-invariant form
`V diag(f(lam)) inv(V)`. Repeated principal stretches occur already in
uniaxial stretch and pure dilation, so this example exercises the scale-
and rotation-safe `eig33z` fallbacks of the generated Fortran, on both
the value path and the complex-step tangent path.

## Theory and conventions

Energy and principal Kirchhoff stresses, with `lb_i = J^{-1/3} lambda_i`
and properties in Abaqus *User Material order `mu, alpha, K`:

```text
W = (2 mu / alpha^2) (lb1^a + lb2^a + lb3^a - 3) + (K/2) ln(J)^2
tau_i = (2 mu / alpha) (lb_i^a - (1/3) sum_j lb_j^a) + K ln(J)
S = V diag(tau_i / lam_i) inv(V)  on  eig(C),   P = F S
```

The generated UMAT converts to Cauchy stress and the Jaumann-rate
`DDSDDE`. No `STATEV` slots. Within the quasi-repeated eig band only
eigenspace-invariant reconstructions such as this one are supported;
derivatives of individual eigenvectors are outside the package contract.

The generated UMAT is strictly three-dimensional (`NDI=3, NSHR=3,
NTENS=6` required); plane-stress, plane-strain, and axisymmetric use
are unsupported.

Abaqus energy outputs (`SSE`, `SPD`, `SCD`) are not maintained by the
generated UMAT, so Abaqus energy-balance output is not meaningful for
this material. Property values are illustrative in any consistent
unit system; the checks in this bundle define the tested parameter
and state range.

Deliberately omitted: multiple Ogden terms (a straightforward extension
of the same reconstruction) and any inelasticity.

## Independent oracle

`check_reference.py` uses three hand-derived, eig-free checks:

1. pure dilation `F = cI`: `sigma = K ln(c^3)/c^3 I` (triple-repeated
   spectrum);
2. isochoric uniaxial: classical principal-stress formula
   `sigma_i = (2 mu/alpha)(l_i^a - mean)` (pair-repeated spectrum);
3. `alpha = 2` degeneracy: one-term Ogden equals isochoric neo-Hookean,
   checked at a rotated, sheared, non-isochoric `F` against the closed
   tensor formula `sigma = [mu dev(bbar) + K ln(J) I]/J`.

A broken control confirms that dropping the `J^{-1/3}` isochoric split is
rejected at `J != 1`.

## Reproduce the pipeline

```bash
python build.py            # Python tangent consistency + generation
python check_reference.py  # dilation/uniaxial/alpha=2 oracles + control
python check_compiled.py   # regeneration parity, gfortran, f2py spectral gates
```

## Recorded evidence

| Facet | Command | Result |
|---|---|---|
| Python tangent consistency | `python build.py` | CS vs FD err 2.0e-07; major symmetry 1.9e-15; `P(I)=0` exact |
| Independent physics oracle | `python check_reference.py` | dilation, uniaxial, and `alpha=2` eig-free closed forms |
| Repeated-spectrum coverage | `python check_reference.py` | pair- and triple-repeated value states pass |
| Broken control | `python check_reference.py` | missing isochoric split rejected |
| Deterministic generation | `python check_compiled.py` | byte-identical regeneration |
| Fortran compile | `python check_compiled.py` | gfortran fixed-form, no errors |
| Direct compiled execution | `python check_compiled.py` | stress parity at rotated/uniaxial/dilation states, 3.6e-15 |
| Compiled closed form | `python check_compiled.py` | runtime `alpha=2` props vs eig-free formula, 1.4e-14 |
| Compiled spectral tangent | `python check_compiled.py` | `DDSDDE` vs FD-Jaumann at distinct (1.4e-06) AND repeated (1.5e-06) spectra |
| Abaqus execution | not included in this bundle | none claimed |

The repeated-spectrum tangent row is the point of this example. The
complex-step perturbation of a repeated-eigenvalue state runs through the
rotation-invariant `eig33z` fallback in the generated Fortran; the
pre-correction near-diagonal guard returned `V = I` there and silently
corrupted exactly this tangent.
