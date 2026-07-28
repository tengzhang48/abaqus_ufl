# abaqus_ufl

Generate **Abaqus** user subroutines (**UMAT** / **UEL**) from a compact symbolic
model definition written in Python, with automatic tangent generation and
material-method tangent verification.

You write the constitutive response (or an element weak form) once, in Python.
`abaqus_ufl` differentiates the implemented methods and emits fixed-form
Fortran. Material `verify()` checks complex-step tangents against finite
differences; a UEL still needs the assembled element checks defined below.

## Install

```bash
pip install -e .
```

or, for a full environment with a Fortran compiler (for `f2py`-compiled checks):

```bash
conda env create -f environment.yml
conda activate abaqus-ufl
```

Requires Python ≥ 3.8, NumPy, and SymPy.

## Quickstart — a UMAT in a dozen lines

```python
import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log


class NeoHookean(au.Material):
    props = dict(G=0.5, K=50.0)          # maps to *User Material constants, in order

    def stress_PK1(self, F):             # 1st Piola-Kirchhoff stress
        J = det(F)
        FinvT = inv(F).T
        return self.G * (F - FinvT) + self.K * log(J) * FinvT


model = NeoHookean()
assert model.verify()                    # complex-step vs finite-difference tangent
au.generate_umat(model, "neo_hookean_umat.for")
```

`verify()` is the first gate: if the tangent it derives doesn't match, it fails
before any Fortran is written. For an element-level weak form, subclass
`au.WeakForm` and call `au.generate_uel(...)` instead.

## The pipeline

The package sits inside a verification pipeline already exercised by the
development examples:

1. map the theory, assumptions, conventions, fields, and state layout;
2. implement the Python material or weak form;
3. check tangents plus model-specific regimes, invariants, and an independent
   quantitative oracle;
4. generate the Fortran reproducibly and compile it;
5. call nontrivial generated code directly through f2py or another checked
   compiled runtime;
6. add a solver run only where it provides useful evidence; and
7. verify the output bridge before comparing or plotting solver results.

`verify()` is one consistency gate in this chain. For a UEL,
`problem.verify()` does not replace an assembled residual/tangent check.
Abaqus mesh, contact, loading, procedure, solver-control, and launch choices
remain example/user-owned; the package does not attempt to automate the whole
FEM-analysis workflow.

The working [`examples/_template/`](examples/_template/) demonstrates the
local Python-to-compiled-UMAT loop. The complete shared contract is in
[`HOWTO_ADD_AN_EXAMPLE.md`](HOWTO_ADD_AN_EXAMPLE.md).

## Layout

```
abaqus_ufl/            the package: core/ (model API, tangents) + generators/ (UMAT/UEL codegen)
examples/              public allowlist + working pipeline demonstration
tools/                 shared Abaqus run / ODB-extract / compare machinery
docs/                  usage, theory, and design documentation
ai_skills/             an operational guide for AI coding assistants
```

Every released example reports separate evidence for theory/oracles, tangent
or element consistency, generated/compiled execution, solver runs, and output
bridges. See [`examples/README.md`](examples/README.md).

## Documentation

- **Usage:** [`docs/API_USAGE.md`](docs/API_USAGE.md) — entry points, the tensor DSL, and the example pipeline.
- **Theory:** [`docs/theory.md`](docs/theory.md), [`docs/complex_step_patterns.md`](docs/complex_step_patterns.md), [`docs/JAUMANN_RESOLUTION.md`](docs/JAUMANN_RESOLUTION.md).
- **Design:** the code-generator and tangent-engine internals — see the [`docs/` index](docs/README.md).
- **Lessons learned:** general Abaqus UMAT/UEL + codegen lessons — see [`docs/lessons/`](docs/lessons/).

## License

MIT — see [LICENSE](LICENSE).
