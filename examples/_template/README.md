# <Example name>

**Demonstrates:** <which package feature — e.g. finite-strain UMAT, ScalarField transport, mixed u-p element, F-bar>
**Element / formulation:** <e.g. C3D8 UMAT | Quad4 UEL, standard>
**Verification:** <Abaqus-validated | Abaqus-smoke | Python-verified | analytical>  ·  <reference, e.g. "vs closed-form uniaxial", "vs Abaqus C3D8 built-in">

<One paragraph: what physics this shows and why it is a useful demonstration.>

## Model

<Key constitutive equations, in the same notation as build.py.>

## Files

| File | Purpose |
|------|---------|
| `build.py` | defines the model, verifies the tangent, generates the `.for` |
| `template_umat.for` | generated UMAT/UEL (committed so it can be read without running) |
| `abaqus/job.inp` | single-element Abaqus deck |
| `abaqus/extract_config.json` | what to pull from the ODB |
| `abaqus/generate_reference.py` | independent Python reference → `reference.json` |
| `abaqus/reference.json` | expected values + tolerances |

## Reproduce

```bash
# 1. author: verify tangents + generate the Fortran
python build.py

# 2. reference: independent Python material-point values
python abaqus/generate_reference.py

# 3. validate in Abaqus (needs Abaqus/Standard on PATH)
cd abaqus && bash run.sh
python ../../../tools/compare_results.py --extracted job_extracted.json --reference reference.json
```
