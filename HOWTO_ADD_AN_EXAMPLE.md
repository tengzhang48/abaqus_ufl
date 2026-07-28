# Example Verification Pipeline

An example is a verification bundle, not a single template file. Different
models can use different scripts, but a released example must make the same
evidence chain clear and reproducible:

```text
theory -> Python -> model checks -> generated Fortran -> compiled execution
       -> optional FE/Abaqus run -> output bridge -> quantitative evidence
```

Copy [`examples/_template/`](examples/_template/) for a working minimal UMAT
implementation of this pipeline. Reuse the responsibilities, not its elastic
physics or exact filenames.

## 1. Scope the theory

Before coding, document:

- equations implemented and deliberately omitted;
- fields, interpolation, and active DOFs for a UEL;
- properties, units, stress/flux signs, and tensor/Voigt conventions;
- complete `STATEV` or `SVARS` order and initialization; and
- the intended independent quantitative oracle.

Do not turn an Abaqus model choice into a package rule. Mesh conversion,
sections, contact, boundary conditions, procedures, solver controls, and launch
configuration remain example/user-owned.

## 2. Implement and check the Python model

Use `Material` or `SmallStrainMaterial` for local constitutive behavior and
`WeakForm` when additional nodal fields or gradient terms require a UEL.

Run `model.verify()` or `problem.verify()` as a tangent-consistency gate for the
implemented material methods. This is necessary but does not prove the
equations match the theory.

Add model-specific checks that:

- reach every claimed branch or regime;
- test physical invariants, signs, limits, and state evolution;
- compare at least one quantity with a closed form, hand calculation, or
  external benchmark that the implementation did not generate; and
- fail when a known defect is deliberately reintroduced.

For UELs, `problem.verify()` is not an assembled-element verification. Add a
reference-assembly residual/tangent check and an appropriate patch, rigid-body,
constant-field, or operator-sign test. Check the declared DOF and state layout
directly.

## 3. Generate reproducibly and compile

The example's build script should verify the applicable Python model and emit
the committed `.for` source. A release check should regenerate to a temporary
path and compare it with the committed artifact before compiling it:

```bash
gfortran -c -ffixed-form -ffixed-line-length-none \
  examples/<name>/<generated>.for -o /tmp/<name>.o
```

Compilation proves Fortran syntax and interfaces only; it is not runtime or
physics evidence.

## 4. Execute the generated subroutine directly

For a nontrivial model, call the generated code before adding solver noise:

- UMAT: compare `STRESS`, `DDSDDE`, and `STATEV`;
- UEL: compare `RHS`, `AMATRX`, `SVARS`, and `PNEWDT`;
- coupled subroutines: check their shared state layout and updates.

The working template ships a tested one-point f2py UMAT driver. A UEL or
coupled example needs its own checked element driver or an equivalent direct
compiled runtime. A feacheap run may serve as the compiled FE rung when that
case has a reliable, checked driver. Neither route removes the need for an
independent physics oracle.

## 5. Add solver evidence only when useful

A tiny feacheap or Abaqus run can check assembly, conventions, and the solver
interface. A paper-scale solve is not required for every example, and Abaqus
setup is not standardized by this package.

If an Abaqus bundle is included, keep its setup explicit and scoped. The shared
`tools/run_case.sh`, `tools/extract_odb.py`, and `tools/compare_results.py`
support simple one-element cases. Complex cases may and should use
example-owned extraction and comparison scripts.

## 6. Verify the output bridge

Output is a separate interface. A successful job does not prove that an ODB,
VTK, NPZ, CSV, or visualization field contains the intended quantity.

For every quantity used as evidence or shown in a plot, document:

- its logical meaning and source field/component/slot/active DOF;
- units, sign/offset, component order, integration-point order, and any
  interpolation or projection;
- its node or element/integration-point identity;
- expected, observed, unique, missing, duplicate, and non-finite coverage; and
- whether it is authoritative solution output, a reconstructed value, or a
  visualization-only bridge.

For `UVARM`, dummy-element, or projected visualization fields, compare the
bridge pointwise with the authoritative solved field when possible. Finite
values and plausible extrema are not coverage checks.

The generic extractor is intentionally limited to simple cases. Do not extend
it with model physics; put histories, named-set reductions, projections, and
derived observables in the example that defines them.

## 7. Record evidence honestly

The example README should list the exact command and result for each applicable
facet:

- theory/independent oracle;
- tangent or assembled-element consistency;
- regime and broken-control checks;
- deterministic generation and compile;
- direct compiled execution;
- solver execution;
- output-bridge coverage/parity; and
- quantitative benchmark or paper comparison.

These facets are not one linear "verification level." State precisely what
passed and what remains unavailable or diagnostic-only.

## Release gates

- [ ] Equations, non-scope, properties, conventions, fields/DOFs, and state
      layout are documented.
- [ ] Python checks enter every claimed regime and include an independent
      quantitative oracle.
- [ ] A UEL has assembled residual/tangent and patch/invariant evidence; it
      does not rely on `problem.verify()` alone.
- [ ] Generated Fortran is reproducible, committed, and compiled.
- [ ] Nontrivial generated code is directly exercised through f2py, feacheap,
      or an equivalent checked runtime.
- [ ] Every claimed solver output has a documented and checked output bridge.
- [ ] Frozen expected results are not regenerated during ordinary verification.
- [ ] No third-party source is included without its redistribution license.
- [ ] No PDFs, ODBs, large outputs, caches, private notes, unpublished models,
      or absolute HPC paths are shipped.
- [ ] The folder is present in the release and only then added to
      [`examples/README.md`](examples/README.md).

The public example manifest is an allowlist, not a map of every model ever
implemented in the development lab. Internal, unpublished, license-unclear, and
not-yet-ported capabilities stay out of the v1 gallery and AI example map.
