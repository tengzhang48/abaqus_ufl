# Examples

Each example is one self-contained folder: a Python model definition, the
generated Abaqus `.for`, and an `abaqus/` bundle (deck + independent reference)
that validates it. The generic run/extract/compare machinery is shared in
[`../tools/`](../tools/). To add your own, copy [`_template/`](_template/) and
see [`../HOWTO_ADD_AN_EXAMPLE.md`](../HOWTO_ADD_AN_EXAMPLE.md).

**Verification levels.** *Abaqus-validated* = job runs and matches an
independent reference within tolerance. *Abaqus-smoke* = compiles/links and
converges in Abaqus, no quantitative reference. *Python-verified* = tangent +
material-point reference, no Abaqus run. Every example states its own level.

## Released examples

| Example | Demonstrates | Type | Verification |
|---|---|---|---|
| `small_strain_elastic_umat` | minimal UMAT, stress/tangent | UMAT · C3D8 | Abaqus-validated |
| `neo_umat_c3d8` | finite-strain hyperelastic UMAT | UMAT · C3D8 | Abaqus-validated |
| `mooney_rivlin_umat` | two-term hyperelastic UMAT | UMAT · C3D8 | Abaqus-validated |
| `small_strain_j2_umat` | plasticity with history/state | UMAT · C3D8 | Abaqus-validated |
| `small_strain_viscoelastic_umat` | rate/time-dependent UMAT | UMAT · C3D8 | Abaqus-validated |
| `neo_fbar_quad4` | F-bar volumetric-locking control | UMAT · Quad4 | Abaqus-validated |
| `thermo_mechanics_quad8` | coupled weak-form UEL (ScalarField) | UEL · Quad8 | Abaqus-validated |
| `uel_scaffold_quad4` | minimal UEL skeleton | UEL · Quad4 | Abaqus-smoke |
| `uel_scaffold_mixed_quad8` | mixed-field UEL skeleton | UEL · Quad8 | Abaqus-smoke |

<!--
CURATION NOTE: this table is the positive release manifest. List ONLY examples
that are Abaqus-tested AND license-clear. Any internal or unpublished port, and
any model whose redistribution license is unresolved, must NOT be listed or
shipped. Grow the table only as new examples clear both the verification gate
and the license gate.
-->
