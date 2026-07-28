# Verification Pipeline

Use the smallest pipeline gate that can expose the suspected bug. Do not use
Abaqus to debug a Python equation error, and do not call a model paper-validated
because a generated `.for` file compiles.

The layers below are an efficient execution order, not one scalar verification
level. Record the exact quantity, oracle, runtime, and output bridge that passed.

## Layer 0: Check the evidence record first

Before treating a discrepancy as a solver bug, establish the exact quantity and
configuration already checked, the oracle or comparison source, the runtime,
and the output path. Read the checked deck rather than re-deriving its
parameters so setup mismatches are separated from implementation defects.

- A constitutive or equilibrium check against an independent oracle is not the
  same as agreement with a new external run.
- Abaqus agreement is production-interface and convention evidence unless the
  comparison itself supplies an independent quantitative benchmark.
- A deliberate extreme may be useful as a stress test without defining the
  realistic operating envelope. Record both roles instead of conflating them.

## Layer 1: Python Reference

Typical commands:

```bash
python examples/<model>/build.py
python examples/<model>/check_reference.py
python examples/<model>/run_material_point.py
```

Good checks:

- zero stress/zero flux at reference state;
- elastic limit matches analytical response;
- state variables remain finite;
- plastic strain, damage, or hardening evolves monotonically when expected;
- the test actually enters the claimed regime, rather than staying elastic;
- storage is zero when old and current fields match;
- flux sign and magnitude match the strong form.
- for branch-dependent gradient damage, the gradient coefficient and reaction
  coefficient use the same branch parameters; measure `ell_eff` per branch.
- when an internal-variable evolution law has a closed-form regime solution,
  compare the accumulated state with that independent solution rather than
  checking only qualitative monotonicity.

## Layer 2: Framework Verification

Use:

```python
model.verify()
problem.verify()
```

For a material, this checks tangent consistency with the implemented
constitutive response. `WeakForm.verify()` currently delegates to material
verification; it does not assemble a UEL residual/tangent. It therefore does
not prove the weak form, DOF layout, quadrature, or generated element is
correct, and neither method proves the equations match the paper.

For branchy models, do not move straight from `verify()` to a solver run. Add a
material-point or element-level regime sweep that proves each branch is entered
and checks at least one branch-specific invariant. A branch mismatch between a
reaction term and its gradient coefficient can pass `verify()` because CS and
FD differentiate the same wrong residual.

For a UEL, add a Python reference-assembly check of `RHS` and `AMATRX`, an
element finite-difference tangent check, and an appropriate constant-field,
patch, rigid-body, or operator-sign test before generation is called verified.

## Layer 3: Generated Fortran Compile

Compile every generated source:

```bash
gfortran -c -ffixed-form -ffixed-line-length-none \
    examples/<model>/<generated>.for -o /tmp/<model>.o
```

Compile is necessary but not sufficient. Some invalid generated Fortran can
compile when helpers have no explicit interface and then fail at runtime. For
finite-strain UMAT/UEL tensor returns, grep for impossible component-indexed
matrix helper calls before moving to f2py or solver runs:

```bash
rg "det33z\\([^)]*\\(ii,jj\\)" examples/<model>/<generated>.for
```

The historical bad pattern was `det33z(F(ii,jj))`, caused by component-wise
tensorization of a direct return containing `log(det(F))`. If this appears,
fix/generate again before running Abaqus.

## Layer 4: f2py Point/Element Checks

Use f2py when generated Fortran compiles but solver behavior would add noise.
This layer calls a single subroutine:

| Target | Outputs |
|---|---|
| UMAT | `STRESS`, `DDSDDE`, `STATEV` |
| UEL | `RHS`, `AMATRX`, `SVARS`, `PNEWDT` |

With Python 3.12+ and NumPy 2.x, f2py needs Meson and Ninja.

The working `examples/_template/check_compiled.py` and
`examples/_template/f2py/drive_umat.f90` demonstrate regeneration parity,
compile, and a direct UMAT call. UELs need an example-owned element driver that
checks the generated `RHS`, `AMATRX`, `SVARS`, and `PNEWDT`.

## Layer 5: Optional Solver and Output Bridge

Use the per-example `abaqus/` bundle (see `examples/_template/abaqus/`) for
production-convention checks, with the shared harness in `tools/`.

Case layout:

```text
examples/<name>/abaqus/
  job.inp
  extract_config.json
  generate_reference.py
  run.sh
```

plus the generated `*.for` in the example folder.

Shared harness files:

- `tools/extract_odb.py`
- `tools/compare_results.py`
- `tools/reference_utils.py`
- `tools/run_case.sh`

UMATs usually extract `S`, `E`, and `SDV`. UELs often need nodal `U`/`RF`
or custom state extraction because standard stress output may not exist.

General Abaqus model construction is outside the package pipeline: mesh,
contact, sections, loading, procedures, solver controls, and launch setup
remain example/user-owned. The output bridge is in scope.

Treat solver completion and output verification as separate gates:

1. identify each intended logical quantity and its field/component/slot/active
   DOF, units, sign/offset, component order, and integration-point order;
2. identify the authoritative solved field separately from reconstructed or
   visualization-only output;
3. require complete, unique, finite node or element/integration-point coverage
   using stable identities;
4. audit `UVARM`, dummy-element, or projected visualization bridges pointwise
   against authoritative output where possible; and
5. only then compare the extracted quantity with a frozen independent
   reference.

The shared extractor is for simple one-element cases. Histories, named-set
reductions, projections, and derived observables belong in an example-owned
extractor. If a quantity has no checked bridge, label it unavailable or
diagnostic-only rather than inferring it from job completion.

## Layer 6: Paper Reproduction

Only attempt paper figures after smaller rungs are stable. Published models
often hide:

- missing parameters;
- implicit boundary schedules;
- frame choices;
- solver stabilization;
- geometry or scale effects;
- postprocessing definitions.

Use exact status labels:

- **code implementation complete**: scoped equations implemented, generated,
  compiled, and checked through an appropriate direct-runtime rung.
- **solver stabilization open**: residual/tangent exist but free solve is not
  robust.
- **paper reproduction complete**: geometry, boundary conditions, parameters,
  response metric, and figure comparison are all matched.

For softening/localization reproductions: the selected SHEAR-BAND PATTERN
(crossed-X vs single dominant band, and its orientation) is a DEGENERATE-
bifurcation outcome decided by perturbations (explicit round-off, seeded
imperfections, or a symmetric deterministic solve that cannot pick a winner)
— the reference paper may itself show different patterns across its own runs.
Gate the reproduction on pattern-insensitive observables (force-displacement
curve, peak, drop location, band-angle family), report the pattern and its
perturbation source honestly, and never pass/fail on band multiplicity alone.
