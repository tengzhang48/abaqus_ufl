# Validation Ladder

Use the smallest validation layer that can expose the suspected bug. Do not use
Abaqus to debug a Python equation error, and do not call a model paper-validated
because a generated `.for` file compiles.

## Layer 0: CHECK THE VALIDATION RECORD FIRST (before debugging a discrepancy)

Before treating a discrepancy as a solver bug, establish **what is already validated**:
the layer it was validated at, the **oracle** (Abaqus / analytical), and the
exact **quantity**. Then compare against *that* reference, don't re-litigate a passed
validation, and read the validated DECK rather than re-deriving its parameters (this is how
you catch setup mismatches like a wrong `theta` early). Two distinctions that save days:

- **"The MODEL is validated"** (constitutive / equilibrium, against an oracle) is NOT the
  same as **"this new comparison matches a specific external run."** The latter failing
  (e.g. a transient-rate gap vs one Abaqus run) is not a model-validation failure and is
  often lower-stakes.
- **Weight effort by how REALISTIC and how VALIDATED the regime is.** Benchmark decks often
  use deliberate extremes (e.g. a near-dry initial state that swells several-fold) as stress
  tests; those extremes drive most of the numerical difficulty. Don't pour runs into an
  extreme, already-validated, unrealistic corner -- test the realistic regime instead.

(Real example: a transient-rate gap chased across many runs as a "validation failure" was,
in fact, a new diffusion-transient comparison at an unrealistic, already-validated stress
case -- while the model itself was already validated against an independent oracle. Checking
the record first would have reframed it immediately.)

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
- when an internal-variable ODE has a closed-form regime solution (Voce
  hardening, exponential relaxation), integrate the increment WITH that exact
  form — then asserting the state against the paper's integrated equation is a
  free INDEPENDENT oracle. Prefer test states where the criterion under test
  collapses to a parameter (a yield locus in simple tension reducing to a single
  parameter), and build uniaxial-stress states exactly
  (Hencky: `lam_lat = lam_ax**(-nu)`) — high-exponent loci amplify a "roughly
  uniaxial" state into a double-digit threshold error.

## Layer 2: Framework Verification

Use:

```python
model.verify()
problem.verify()
```

This checks tangent consistency with the implemented residual. It does not
prove the residual matches the paper.

For branchy models, do not move straight from `verify()` to a solver run. Add a
material-point or element-level regime sweep that proves each branch is entered
and checks at least one branch-specific invariant. A gradient-damage UEL bug can
pass `verify()` because both CS and FD differentiate the same wrong residual:
`phase_storage` used the ductile `psi_star` in compression while `phase_flux`
used the brittle value. A per-branch regularization-length check caught it.

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
| UMATHT | `FLUX`, `DUDT`, `DUDG`, `DFDT`, `DFDG` |
| HETVAL | heat/source and tangent |
| UEL | `RHS`, `AMATRX`, `SVARS`, `PNEWDT` |

With Python 3.12+ and NumPy 2.x, f2py needs Meson and Ninja.

## Layer 5: Abaqus Validation

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
  compiled, and validated at a small rung.
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
