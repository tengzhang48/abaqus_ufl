# abaqus_ufl

`abaqus_ufl` helps researchers turn a supported constitutive model or
coupled-field element declaration written in Python into inspectable,
self-contained Fortran for Abaqus/Standard.

The practical aim is not simply to write Fortran faster. It is to separate the
parts of a user subroutine that are otherwise difficult to audit: field and
state definitions, constitutive responses, residual terms, tangent blocks,
local-variable treatment, DOF packing, and Abaqus interface conventions.

```text
Python declaration
  → consistency and model-specific checks
  → generated fixed-form Fortran
  → compiled element or material checks
  → user-owned Abaqus analysis
```

## What it is—and what it is not

The package currently generates:

- **UMAT** source from supported material-response declarations; and
- **UEL** source from supported field, interpolation, balance, material, and
  local-variable templates.

The name reflects an important design influence: the interface is inspired by
the declarative style of the Unified Form Language used in the FEniCS
ecosystem. `abaqus_ufl` is **not UFL-compatible**. It neither depends on nor
implements FEniCS UFL, and it is not a general compiler for arbitrary
variational expressions.

The package also does not generate an entire finite-element analysis. Meshes,
sections, contact, loads, boundary conditions, procedures, solver controls,
units, and physical validation remain the responsibility of the researcher
and the individual example.

## Who it is for

`abaqus_ufl` is intended for researchers who:

- need a custom Abaqus material or coupled-field user element;
- want the problem-specific declaration to be shorter and more reviewable than
  a monolithic hand-written UEL;
- need generated tangent blocks and repeatable source generation; and
- are willing to verify the material, assembled element, Abaqus setup, and
  scientific result at the level required by their claim.

It is a research tool, not a substitute for finite-element formulation
knowledge or an assurance that a declared model is physically valid.

## Try a complete local workflow

The Neo-Hookean example demonstrates the path from Python declaration through
an independent closed-form check to a compiled call of the generated UMAT.
From a clone of this repository:

```bash
pip install -e ".[dev]"
cd examples/neo_hookean_umat
python build.py
python check_reference.py
python check_compiled.py
```

The last command requires `gfortran`, f2py, Meson, and Ninja. The development
extra above installs the Python-side build tools. A conda environment
specification is also provided; create it from the repository root:

```bash
conda env create -f environment.yml
conda activate abaqus-ufl
```

The example's Python declaration is in
[`examples/neo_hookean_umat/build.py`](examples/neo_hookean_umat/build.py).
For a deliberately simple directory that can be copied and adapted, see
[`examples/_template/`](examples/_template/).

Requirements for Python-only generation are Python 3.8 or newer, NumPy, and
SymPy.

## What a passing check means

The package distinguishes several kinds of evidence because they catch
different failures:

- `Material.verify()` compares the implemented complex-step tangent with
  finite differences for the selected state.
- A model-specific oracle checks a limit, invariant, analytic solution, or
  independently implemented response.
- Deterministic regeneration checks that the committed Fortran matches the
  declaration and generator revision.
- A compiled call checks the actual generated subroutine boundary.
- An assembled UEL check tests residuals, tangent blocks, DOF/state layout,
  quadrature, and local-variable behavior.
- An Abaqus run checks the selected deck and solver path.
- An output-bridge audit checks that the field plotted or compared is the
  field the element actually computed.

These checks are complementary. `verify()` is a first consistency gate; it
does not establish that the governing equation is correct. Python and
generated Fortran can reproduce the same mistake. Likewise, solver completion
does not by itself establish quantitative reproduction or physical
validation.

The full example contract, including known-broken controls and clean-release
checks, is in [`HOWTO_ADD_AN_EXAMPLE.md`](HOWTO_ADD_AN_EXAMPLE.md).

## Examples

[`examples/README.md`](examples/README.md) describes the curated public
examples and the evidence each one actually carries. That allowlist is a
release subset, not a capability table for the larger research project.

[`paper_examples/`](paper_examples/), available in repository checkouts but
excluded from the Python source distribution, contains the declarations,
generated sources, decks, reduced data, and figure materials assembled around
the manuscript examples. Read each package README for its evidence level and
provenance; code-to-code reproduction, component verification, and execution
demonstration are not interchangeable labels.

**Release provenance note.** The accepted gel-bilayer deck is our deck for the
swell-induced bending problem of Chester, Di Leo, and Anand; its mesh
discretization follows their supplemental example, with attribution. Their
original supplemental files are not redistributed. Regenerating the deck with
`build.py` requires a separately obtained copy of their supplemental input as
the mesh seed. The corrosion comparison mesh retains third-party lineage for
which the precise BSD notice and redistribution status must still be recorded.
See [`CREDITS.md`](CREDITS.md).

## Documentation

- [API usage](docs/API_USAGE.md): entry points, declaration methods, tensor
  operations, UMAT/UEL generation, and example workflows.
- [Theory and conventions](docs/theory.md): residual, sign, field, and gel
  conventions.
- [Complex-step patterns and limits](docs/complex_step_patterns.md): analytic
  paths, branching, state, and verification.
- [Design documentation](docs/README.md): generator and tangent-engine
  internals.
- [Lessons learned](docs/lessons/): distilled Abaqus, Fortran, code-generation,
  validation, and release lessons.
- [Fresh-clone Abaqus validation](docs/ABAQUS_VALIDATION_2026-07-30.md):
  current-source smoke solves, paper-case datachecks, and remaining scope.
- [AI-assistant guide](ai_skills/README.md): operational project guidance for
  coding agents.

## Scientific lineage and credit

The framework grew from prior coupled-mechanics research and from careful
study of shared or published UEL/UMAT implementations. Important sources
include:

- Shawn A. Chester, Claudio V. Di Leo, and Lallit Anand's gel theory and
  supplemental UEL/decks;
- Bibekananda Datta and Thao D. Nguyen's modular hydrogel UEL;
- Professor Allan Bower's EN234_FEA teaching code at Brown University, the
  basis of the internal `feacheap` validation host;
- Chuanjie Cui, Rujin Ma, and Emilio Martínez-Pañeda's reference
  stress-corrosion formulation and UEL.

These sources contributed theory, conventions, execution infrastructure, or
independent benchmarks in different ways. They are not all incorporated into
this repository, and the MIT license of `abaqus_ufl` does not replace their
licenses. Exact roles, citations, redistribution boundaries, discussion
acknowledgements, and open provenance items are recorded in
[`CREDITS.md`](CREDITS.md).

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). Until a
release DOI is recorded there, cite the repository URL and the exact version
or commit used.

## License

Project-authored `abaqus_ufl` package source is MIT licensed; see
[`LICENSE`](LICENSE). Separately licensed, attribution-only, or
provenance-pending reference and derived artifacts are governed by their own
terms and are identified in [`CREDITS.md`](CREDITS.md).
