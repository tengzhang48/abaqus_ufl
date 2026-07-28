# Examples

A released example is a self-contained verification bundle. It documents the
theory and scope, implements the Python model, generates committed Fortran,
checks the applicable regimes and assembled behavior, directly exercises the
compiled code, and records any solver/output-bridge evidence it actually has.

There is no single linear verification level. Each example reports separate
facets:

- independent physics oracle;
- Python tangent or assembled-element consistency;
- regime and deliberately broken-control coverage;
- reproducible generation and Fortran compile;
- direct f2py or checked FE execution;
- optional Abaqus execution;
- output-bridge identity, coverage, and parity; and
- quantitative benchmark or paper comparison.

The package does not prescribe general Abaqus model construction. Optional
`abaqus/` folders are example-owned integration checks. The shared tools are
appropriate for simple one-element extraction; difficult field mappings,
histories, projections, and derived observables remain with the example that
defines them.

To add an example, start with the working [`_template/`](_template/) and follow
the [`example verification pipeline`](../HOWTO_ADD_AN_EXAMPLE.md).

## Released examples

The table below is the public release allowlist. Add a row only after the
folder is present, license-clear, and its README records the evidence above.
The development repository's larger model inventory is intentionally not a
public v1 capability map.

| Example | Demonstrates | Generated target | Evidence |
|---|---|---|---|
| [`neo_hookean_umat/`](neo_hookean_umat/) | finite-strain stateless law; Cauchy + Jaumann `DDSDDE` conversion | UMAT | closed-form uniaxial/shear oracle + broken control; deterministic generation; compiled f2py stress 1.8e-15, FD-Jaumann tangent 1.3e-06; no Abaqus run |
| [`small_strain_j2_umat/`](small_strain_j2_umat/) | stateful radial return with elastic/plastic branch; `STATEV` threading; engineering-shear boundary | UMAT | exact shear closed form (radial return exact) + broken control; compiled 40-increment path, stress 5.6e-17, `STATEV` 1.4e-17; no Abaqus run |
| [`small_strain_viscoelastic_umat/`](small_strain_viscoelastic_umat/) | history/`dt` dependence; tensor state in column-major `STATEV(1..9)`; user-convention state sign | UMAT | exact discrete relaxation closed form + explicit-dashpot broken control; compiled 30-increment history, stress 5.6e-17, tensor `STATEV` 5.2e-18; no Abaqus run |
| [`ogden_umat/`](ogden_umat/) | explicit spectral (`eig`) constitutive path; eigenspace-invariant reconstruction; repeated-spectrum tangents | UMAT | eig-free dilation/uniaxial/`alpha=2` oracles + broken control; compiled parity 3.6e-15; FD-Jaumann tangent at distinct AND repeated spectra ~1.5e-06; no Abaqus run |
| [`scalar_diffusion_uel/`](scalar_diffusion_uel/) | coupled two-field UEL (u + T); cross-field tangent blocks; transport template; node-interleaved DOF maps | UEL | closed-form thermal-stress nodal forces + heat-balance invariants on a distorted element; assembled + compiled FD tangents 1.0e-07; f2py `RHS`/`AMATRX` parity vs reference assembly 2.9e-16; no Abaqus run |
| [`thermo_mechanics_quad8/`](thermo_mechanics_quad8/) | mixed-order fields (quadratic u, bilinear corner T); node-dependent DOF maps; deformation-dependent `C^{-1}` flux | UEL | pulled-back-flux element oracle with exact `1/l^2` ratio + broken control; heat balance on distorted Quad8; assembled/compiled FD tangents 1.5e-07; f2py parity 2.6e-16; no Abaqus run |

The `_template` is a working pipeline demonstration, not a released scientific
example.
