# Working Example-Pipeline Template

This directory is a complete, deliberately simple small-strain elastic UMAT
example. The physics is not the template. The reusable template is the
evidence pipeline demonstrated by the directory:

```text
theory and scope
  -> Python implementation and tangent consistency
  -> model-specific independent check
  -> deterministic Fortran generation and compile
  -> direct compiled-subroutine execution
  -> optional solver execution
  -> output bridge and quantitative comparison
```

Real examples may use different filenames and different solvers. They should
preserve these responsibilities and show the exact commands that produced
their evidence. There is no universal Abaqus model builder or one-size-fits-all
verification script.

## Theory and conventions

The material uses compression-positive tensors in Python:

```text
sigma_(n+1) = sigma_n + 2 G Delta-epsilon
                        + lambda tr(Delta-epsilon) I
```

with `G = 10` and `lambda = 20`. The generated UMAT converts to and from
Abaqus' tension-positive Voigt convention at the subroutine boundary.

The independent benchmark is homogeneous constrained strain with final
stretch `lambda_x = 0.98`:

```text
epsilon_x^Abaqus = log(0.98)
epsilon_x^cp     = -log(0.98)
sigma_x^cp       = (2 G + lambda) epsilon_x^cp
sigma_y^cp       = sigma_z^cp = lambda epsilon_x^cp
sigma_mises      = 2 G epsilon_x^cp
```

The formula is evaluated independently in `check_reference.py`; it does not
call the implemented stress update to construct its expected values.

## Reproduce the local pipeline

From this directory:

```bash
# Python tangent consistency and deterministic source generation
python build.py

# Model-specific physics oracle
python check_reference.py

# Regeneration parity, gfortran compile, and an actual f2py UMAT call
python check_compiled.py
```

`check_compiled.py` uses the proven one-point f2py driver pattern. It regenerates
the source in a temporary directory, byte-compares it with the committed
`template_umat.for`, compiles it, calls the UMAT for ten increments, and checks
both `STRESS` and `DDSDDE` against the closed form.

## Adapting the pipeline

When copying this folder, replace the material and every model-specific oracle.
Do not keep the elastic check and merely change its expected numbers.

For a stateful UMAT, the direct compiled check should exercise every claimed
branch and compare `STRESS`, `DDSDDE`, and the documented `STATEV` layout.

For a UEL, `problem.verify()` alone is not an element check. The example must
also exercise the assembled residual and tangent, DOF/state layout, quadrature,
and an appropriate patch or invariant. Its compiled check should call the
generated UEL and compare `RHS`, `AMATRX`, `SVARS`, and `PNEWDT` with the Python
reference assembly. A model-specific f2py element driver replaces the UMAT
driver in this folder; the package API does not need to change.

An independent quantitative oracle remains required because tangent
consistency, code generation, and Python-versus-generated-Fortran agreement can
all reproduce the same wrong equation.

## Optional Abaqus and the output bridge

The `abaqus/` directory is an optional one-element integration check. It is not
a general Abaqus-analysis template. Users and examples own their mesh,
sections, procedures, contact, boundary conditions, solver controls, and
machine-specific launch setup.

The output bridge is in scope. This simple case uses one homogeneous C3D8R
element with one integration point, so the shared one-record extractor covers
the complete element-output set. Its section explicitly supplies the
UMAT hourglass stiffness `0.005 G = 0.05`, which Abaqus cannot infer from a
user-defined material:

```bash
cd abaqus
bash run.sh
```

The committed `reference.json` is a frozen closed-form expectation. Ordinary
verification consumes it; `generate_reference.py` is run only when the
documented benchmark changes.

For a multi-element, mixed-field, projected, or history-dependent case, use an
example-owned extractor. The example must document and check:

- the logical quantity-to-Abaqus field, component, slot, or active-DOF mapping;
- units, sign/offset convention, component order, and integration-point order;
- stable node or element/integration-point identity and exact expected
  coverage, including duplicates, missing records, and non-finite records;
- which field is the authoritative solved quantity and which is only a
  visualization bridge; and
- pointwise bridge parity when a `UVARM`, dummy-element, or projected field is
  used for plots.

Solver completion does not imply that every intended observable has a working
bridge. Record each observable as quantitatively checked, diagnostic-only, or
unavailable with a reason.

## Files

| File | Pipeline responsibility |
|---|---|
| `build.py` | Python model, tangent check, deterministic Fortran generation |
| `check_reference.py` | independent, model-specific physics oracle |
| `template_umat.for` | committed generated artifact |
| `f2py/drive_umat.f90` | direct compiled UMAT adapter |
| `check_compiled.py` | regeneration, compile, and compiled-runtime checks |
| `abaqus/job.inp` | optional example-owned one-element Abaqus model |
| `abaqus/extract_config.json` | simple-case ODB field selection |
| `abaqus/reference.json` | frozen expected quantities and tolerances |
| `abaqus/run.sh` | optional run, extraction, and comparison |

See [`../../HOWTO_ADD_AN_EXAMPLE.md`](../../HOWTO_ADD_AN_EXAMPLE.md) for the
release gates and the pipeline contract shared by all examples.
