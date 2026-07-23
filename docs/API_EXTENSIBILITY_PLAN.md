# abaqus_ufl API Extensibility

## Purpose

The package should stay easy to use for Abaqus UEL/UMAT generation while
remaining easy to extend to coupled multiphysics. The design goal is not a
fully generic finite-element framework. The goal is a small Abaqus-oriented DSL
where user code reads close to theory and the generator can still emit concrete,
debuggable Fortran.

The method-name style is central to this:

```python
def stress_PK1(self, F, p, mu): ...
def pressure_resid(self, F, p, mu): ...
def solvent_flux(self, F, p, mu, grad_mu): ...
def solvent_storage(self, F, F_old, p, p_old, dt): ...
```

This is clearer for constitutive modelers than a fully generic
`add_balance(...)` registry. Extensibility is added by removing hard-coded
assumptions first, then extracting repeated physics patterns from real examples.

## Design Assessment

### Strengths

- Material methods read like the equations in a paper.
- Complex-step tangent generation is a strong core design.
- Weak-form method names make common coupled models concise.
- Generated Fortran remains close enough to Abaqus UEL/UMAT conventions to
  debug.

### Design tensions

A few tensions shape the extensibility model:

- Element topology should live in one place rather than being assumed inside
  `WeakForm`.
- Local element variables, such as the condensed pressure in the local-pressure
  UEL, benefit from being first-class API objects rather than ad-hoc storage.
- Field argument names should not be hard-coded to `p`, `mu`, and `grad_mu`.
- Formulation selection should be explicit rather than split across flags and
  separate generator entry points.
- The coupled local-pressure path is still a narrow prototype rather than a
  general formulation path.

## Design Principle

Two stages:

1. **Remove hard-coded assumptions.** Fix the concrete design problems that
   block coupled multiphysics: element layout, local variables, field naming,
   and formulation selection.

2. **Extract physics templates.** After several concrete multiphysics examples
   exist, extract their repeated structure into first-class `Physics` templates.

A broad declarative equation registry is deliberately avoided until the patterns
are real. Generalization-by-extraction is safer than generalization-by-prediction.

## Extensibility Model

The first stage — removing hard-coded assumptions — is in place. It preserves
the method-name style and keeps `fbar` as a compatibility option while adding an
explicit formulation route.

### 1. Element layout in `ElementConfig`

`ElementConfig` owns:

- total node count,
- corner node count,
- field support by degree,
- node sets for high-order and low-order fields,
- default DOF ordering rules,
- shape-function template selection,
- quadrature and mapping routines.

`WeakForm` asks the selected `ElementConfig` for layout information instead of
assuming Quad8 topology internally:

```python
cfg = ELEMENT_CONFIGS["quad4"]
layout = cfg.build_layout(fields)
```

`WeakForm` may still store `_dof_map`, `_field_nodes`, and `_ndofel`, but these
are computed through the element config.

### 2. `LocalScalar`

Local variables are not nodal DOFs. They are element-level or integration-point
variables stored in `SVARS` and optionally condensed from the global system:

```python
self.p = au.LocalScalar("p", storage="SVARS", interpolation="constant")
```

Required metadata:

- name,
- shape/rank,
- storage location,
- interpolation (initially only `"constant"`),
- whether it is condensed,
- initial guess / state variable layout.

The first production use is local pressure:

```
global fields: u, mu
local field:   p
```

The local solve and condensation are intended to become a general formulation
feature rather than a one-off generator.

### 3. Generic field argument names

Fields can be declared with domain-specific names:

```python
self.T = au.ScalarField("T", degree=1)
self.c = au.ScalarField("c", degree=1)
self.eta = au.ScalarField("eta", degree=1)
```

The framework infers recognized arguments from the declared fields:

- `T` means scalar field value,
- `grad_T` means its gradient,
- `T_old` means previous-step value,
- similarly for `c`, `eta`, `mu`, etc.

The special mechanical arguments are kept:

- `F`,
- `F_old`,
- possibly `J` later only as a derived convenience, not an independent DOF.

This removes the hard-coded dependence on `mu` while preserving readable method
signatures.

### 4. Explicit `Formulation`

The formulation is selected explicitly:

```python
au.generate_uel(problem, path, element="quad4", formulation="standard")
au.generate_uel(problem, path, element="quad4", formulation="fbar_mechanics")
au.generate_uel(problem, path, element="quad4", formulation="local_pressure")
```

`fbar=True/False` remains as a compatibility alias.

Built-in formulations:

| Formulation | Purpose |
|-------------|---------|
| `standard` | Standard nodal UEL assembly |
| `fbar_mechanics` | Pure mechanical low-order F-bar |
| `local_pressure` | Global nodal fields plus condensed local pressure |

Important rule:

- `fbar_mechanics` is for pure mechanical low-order elements.
- Coupled swelling/diffusion should use `local_pressure` or a global mixed
  `u,p,mu` formulation.

## Future Direction: Physics Templates

The second stage — physics templates — is intended once several concrete coupled
examples exist beyond the initial coupled formulation.

Good extraction candidates:

- thermo-mechanics: `u + T`,
- mechanics plus scalar damage: `u + d`,
- two-species diffusion: `u + mu1 + mu2`,
- phase-field style scalar evolution if the weak form stays second order.

The goal is **physics-as-template**, not **equations-as-registry**.

### What a Physics Template Declares

A `Physics` template declares the mathematical shape the generator needs:

```python
class Transport(au.Physics):
    primary_field = "field"
    test_rank = 0
    residual_form = "storage_plus_flux"
    value_args = ["F", "field"]
    gradient_args = ["grad_field"]
    history_args = ["F_old", "field_old", "dt"]
```

Built-in templates might include:

- `Momentum`,
- `LocalConstraint`,
- `Transport`,
- `ThermalDiffusion`,
- `ReactionDiffusion`,
- `ScalarEvolution`.

These templates are meant to replace hard-coded equation-name branches over
time, while the user-facing method-name style stays. For example, a built-in
`transport_equation(...)` method can be backed by a `Transport` template
internally, so users still write code close to theory.

### Coupling Payoff

Once each physics template declares its value and gradient dependencies, the
generator can compose coupled tangents mechanically:

- if momentum depends on `T`, assemble `K_uT`,
- if transport depends on `F`, assemble `K_mu_u`,
- if damage depends on strain energy, assemble `K_d_u`,
- if stress depends on damage, assemble `K_ud`.

The user defines separate physics; the generator detects cross-blocks and uses
complex-step derivatives to fill them.

### Escape Hatch

Some models will not fit the built-in templates. Examples:

- Cahn-Hilliard with mixed splitting,
- higher-order gradient theories,
- phase-field fracture with nonstandard history handling,
- user-specific stabilization terms.

The API includes an explicit custom path:

```python
class MyPhysics(au.Physics):
    residual_form = "custom"
    assembly_hook = my_custom_assembly
```

This keeps the abstraction from pretending to cover every possible PDE.

## What Not To Do

Do not replace the theory-like method names with a generic `add_balance(...)`
registry.

Do not try to make the framework a general finite-element DSL before the next
target physics cases are known.

Do not make coupled F-bar the production route for swelling/diffusion. Keep it
as an experimental/comparison branch.

## Design Goals

The extensibility model satisfies:

- the coupled three-field formulation still generates correctly,
- local-pressure Quad4 is selected through `formulation="local_pressure"`,
- pure mechanical F-bar is selected through `formulation="fbar_mechanics"`,
  with validation that it is single-field mechanics,
- a scalar field named something other than `mu` can participate in value and
  gradient arguments (for example `T`, `T_old`, and `grad_T`),
- element layout is not assumed to be Quad8 inside `WeakForm`; it is built by
  `ElementConfig.build_layout(...)`.

The physics-template direction targets:

- at least one built-in equation method backed by a `Physics` template,
- a user-defined physics template that can add a new second-order scalar balance
  without editing the core generator,
- cross-coupling tangent blocks still generated automatically.
