# abaqus_ufl Example Map

Use this map before opening old example folders. It tells an agent which
existing model is the closest reusable pattern.

Scope: this map lists the released `abaqus_ufl` examples. If a path is missing,
do not assume the repository is broken; run `find examples -maxdepth 2 -type d`
and choose the closest present pattern.

## Clean UMAT Models

| Example | Reuse for | Notes |
|---|---|---|
| `examples/small_strain_elastic_umat` | Minimal small-strain UMAT | Signature, compression-positive convention. |
| `examples/neo_hookean_umat` | Minimal finite-strain UMAT | `stress_PK1(F)` to Abaqus stress/DDSDDE. |
| `examples/mooney_rivlin_umat` | Level-1 finite-strain benchmark | Clean build smoke path. |
| `examples/small_strain_viscoelastic_umat` | Tensor `STATEV` | 9-slot tensor state and transient response. |

## Path-Dependent UMAT Models

| Example | Reuse for | Notes |
|---|---|---|
| `examples/small_strain_j2_umat` | Stateful plasticity | Return `(sigma, state_dict)`, yielding activation. |

## UEL Models

| Example | Reuse for | Notes |
|---|---|---|
| `examples/scalar_diffusion_uel` | Minimal scalar transport UEL | `transport_equation` returns `(storage, flux)`. |
| `examples/phasefield_fracture_uel` | Scalar phase-field UEL | Gradient damage requires UEL, not UMAT. |
| `examples/thermo_mechanics_quad8` | Coupled `u + T` UEL | Temperature-displacement scaffold and dummy thermal material requirements. |
| `examples/uel_scaffold_quad4` | Abaqus UEL input scaffold | Flat include structure, node-label extraction. |
| `examples/uel_scaffold_mixed_quad8` | Mixed UEL scaffold | Active DOF and interpolation sanity checks. |
| `examples/Fbar_uel` | Pure mechanics F-bar | Do not use as default for coupled gels. |

## Validation Workspaces

| Path | Purpose |
|---|---|
| `examples/<name>/abaqus/` | Per-example Abaqus validation decks and reference results. |
| `tools/` | Shared validation, extraction, and comparison scripts. |

## Authoritative Docs

Open these before deep changes:

- `docs/API_USAGE.md`
- `README.md`
- `HOWTO_ADD_AN_EXAMPLE.md`
- `examples/README.md`
