# abaqus_ufl

Generate **Abaqus** user subroutines (**UMAT** / **UEL**) from a compact symbolic
model definition written in Python, with automatic consistent-tangent
verification.

You write the constitutive response (or an element weak form) once, in Python.
`abaqus_ufl` differentiates it, checks the tangent (complex-step vs finite
difference), and emits ready-to-submit fixed-form Fortran.

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

## The workflow

1. **Author** — a `Material` (UMAT) or `WeakForm` (UEL) in Python.
2. **Verify** — `model.verify()` checks the generated tangent.
3. **Generate** — `generate_umat` / `generate_uel` emit the `.for`.
4. **Validate** — run the generated subroutine in Abaqus against an independent
   reference (see `examples/` and `tools/`).

Copy [`examples/_template/`](examples/_template/) to start your own; the full
recipe is in [`HOWTO_ADD_AN_EXAMPLE.md`](HOWTO_ADD_AN_EXAMPLE.md).

## Layout

```
abaqus_ufl/            the package: core/ (model API, tangents) + generators/ (UMAT/UEL codegen)
examples/              worked, verified examples + a copyable _template/
tools/                 shared Abaqus run / ODB-extract / compare machinery
```

Every example declares its **verification level** — Abaqus-validated,
Abaqus-smoke, Python-verified, or analytical — see
[`examples/README.md`](examples/README.md).

## License

MIT — see [LICENSE](LICENSE).
