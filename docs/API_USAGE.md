# abaqus_ufl API Usage Guide

This note is the compact user-facing API reference for writing models with
`abaqus_ufl`. It focuses on the stable authoring and generation API. Evidence
for any released scientific example belongs in the public example manifest,
not in this API guide.

## Import Pattern

Use the package namespace for framework objects and import tensor operations
from the tensor module:

```python
import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log, exp, sqrt
from abaqus_ufl.core.tensor import eye, trace, sym, dev, eig, logm, expm
```

The top-level namespace exposes the main generation API:

```python
au.Material
au.SmallStrainMaterial
au.WeakForm
au.VectorField
au.ScalarField
au.LocalScalar

au.generate_umat(...)
au.generate_small_strain_umat(...)
au.generate_uel(...)

au.generate_inp_scaffold(...)
au.write_job_inp(...)
```

## Choose The Target

Use a finite-strain UMAT when Abaqus built-in elements can handle the fields
and the user subroutine only supplies a local constitutive law:

```python
class NeoHookean(au.Material):
    props = dict(G=1.0, K=100.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * (F - inv(F).T) + self.K * log(J) * inv(F).T

model = NeoHookean()
model.verify()
au.generate_umat(model, "neo_hookean_umat.for")
```

Use a small-strain UMAT when the model is naturally written as an incremental
stress update:

```python
class SmallStrainElastic(au.SmallStrainMaterial):
    props = dict(G=10.0, lam=20.0)
    stress_convention = "compression_positive"

    def stress_update(self, sigma_old, strain_old, dstrain, dt):
        I = eye(3)
        return sigma_old + 2.0 * self.G * dstrain + self.lam * trace(dstrain) * I

model = SmallStrainElastic()
model.verify()
au.generate_small_strain_umat(model, "small_strain_elastic.for")
```

Use a UEL when the model needs extra nodal fields or gradient physics, such as
temperature, concentration, chemical potential, phase field, pore pressure, or
mixed displacement-pressure fields:

```python
class ThermoMaterial(au.Material):
    props = dict(G=1.0, K=10.0, alpha=1e-3, rho_cp=1.0, k=0.5)

    def stress_PK1(self, F, T):
        J = det(F)
        FinvT = inv(F).T
        return self.G * (F - FinvT) + self.K * log(J) * FinvT - self.K * self.alpha * T * FinvT

    # The storage/flux pair MUST use the recognized method names below;
    # the generator emits calls to exactly these names, and the WeakForm
    # constructor rejects a material that lacks them. (Heat notation is
    # fine in comments; 'heat_storage'/'heat_flux' are NOT recognized.)
    def solvent_storage(self, F, F_old, T, T_old, dt):
        return self.rho_cp * (T - T_old) / dt

    def solvent_flux(self, F, T, grad_T):
        return -self.k * grad_T

class ThermoProblem(au.WeakForm):
    material = ThermoMaterial
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)
        self.T = au.ScalarField("T", degree=1, test="theta")

    def momentum_equation(self, v, F, T):
        return self.material.stress_PK1(F, T)

    def transport_equation(self, theta, F, T, grad_T, F_old, T_old, dt):
        return self.material.solvent_storage(F, F_old, T, T_old, dt), self.material.solvent_flux(F, T, grad_T)

problem = ThermoProblem()
problem.verify()
au.generate_uel(problem, "thermo_quad4.for", element="quad4", formulation="standard")
```

Scalar tuple equations use the generator convention
`storage * test - flux . grad(test)`. For conserved transport equations this
`flux` is naturally the physical flux, for example `-k * grad_T` or
`-D * grad_c`. For phase-field damage or other non-conserved gradient-flow
models, do not copy a diffusion flux pattern into `phase_flux`. First write the
weak form with the positive gradient coefficient, then return the adapter value
needed by the tuple convention. For AT2 damage,
`storage = Gc/ell*d - 2*(1-d)*H` and the positive weak-form term
`+Gc*ell*grad(d).grad(eta)` require `phase_flux = -Gc*ell*grad_d`.

When a scalar gradient model has branch-dependent coefficients, derive the
storage/reaction and gradient terms from the same branch selection and thread
the required state into both methods. A CS-vs-FD check can differentiate the
same wrong branch consistently, so add a branch-specific invariant that checks
the declared reaction/gradient relationship. If a production run uses a
diagonal or staggered tangent, make that approximation explicit rather than
changing the physics to reduce coupling.

UEL material methods may also request `time`, which maps to Abaqus total time
`TIME(2)`. Treat it like `dt`: it is a real, non-differentiated parameter, not
a nodal field. This is useful for prescribed schedules, for example:

```python
def prescribed_temperature(self, time):
    T = self.T_start - self.cooling_rate * time
    if T.real < self.T_min:
        T = self.T_min
    return T
```

## Generator Entry Points

### `generate_umat`

```python
au.generate_umat(material, output_path, mat_prefix=None,
                 matrix_backend="iterative")
```

Use this for finite-strain `Material` classes that define `stress_PK1(self, F)`.
The generated UMAT computes Cauchy stress and Abaqus `DDSDDE` from the PK1
stress and the complex-step `dP/dF` tangent.

`matrix_backend` controls generated Fortran matrix functions:

| Backend | Use case |
|---------|----------|
| `"iterative"` | **Default.** Fixed-count iterative helpers (`sqrtm33z_iter`, `logm33z_iter`, `expm33z_iter`, `polar33z_iter`). CS-safe at diagonal-F states. |
| `"eig"` | Eig-based `sqrtm33z`, `logm33z`, `expm33z`, `polar33z`. The eig guards are scale- and rotation-safe (repeated spectra route through a rotation-invariant fallback); available for compatibility/debugging, with compiled spectral gates covering the states in scope. |

Example:

```python
au.generate_umat(model, "model_umat.for")  # default: iterative
au.generate_umat(model, "model_umat_eig.for",
                 matrix_backend="eig")
```

### `generate_small_strain_umat`

```python
au.generate_small_strain_umat(
    material,
    output_path,
    mat_prefix=None,
    update_method="stress_update",
    extra_fortran_files=None,
)
```

The update method signature is:

```python
def stress_update(self, sigma_old, strain_old, dstrain, state_old..., dt):
    return sigma_new, {"state_name": state_new, ...}
```

State variables are declared with `state_vars`:

```python
class MyPlasticity(au.SmallStrainMaterial):
    props = dict(E=100.0, nu=0.3, sigy0=0.2, H=1.0)
    state_vars = dict(ep=0.0)
    stress_convention = "compression_positive"

    def stress_update(self, sigma_old, strain_old, dstrain, ep_old, dt):
        ...
        return sigma_new, {"ep": ep_new}
```

Tensor state variables consume 9 `STATEV` slots. Scalar state variables consume
1 slot. The generator passes old state values as `<name>_old` and writes the
returned dictionary back to `STATEV`.

For complex research models, keep a readable Python reference and add a
generator-friendly method:

```python
au.generate_small_strain_umat(model, "model.for",
                              update_method="stress_update_codegen")
```

Use `@au.fortran_helper` plus `extra_fortran_files=[...]` only when a helper is
truly outside the Python-to-Fortran DSL. Prefer 3x3 tensor rewrites first. The
sidecar path is currently part of the small-strain UMAT generator workflow.

### `generate_uel`

```python
au.generate_uel(
    weakform,
    output_path,
    element="quad8",
    mat_prefix=None,
    fbar=None,
    element_config=None,
    formulation=None,
)
```

Supported element names are defined by `au.ELEMENT_CONFIGS`:

| Element | Typical use |
|---------|-------------|
| `quad4` | Low-order 2D mechanical or scalar UELs |
| `quad8` | 2D mixed/coupled UELs |
| `quad8r` | Reduced-integration 2D coupled UELs |
| `hex8` | Low-order 3D mechanical UELs |
| `hex20` | 3D mixed/coupled UELs |
| `tet4` | Linear tetrahedral UELs, 4-point quadrature (Abaqus C3D4 node ordering) |
| `tet4r` | Linear tetrahedral UELs, single-point quadrature (e.g. stabilized equal-order pairs) |

Use explicit `formulation=` in new examples:

| Formulation | Use case |
|-------------|----------|
| `"standard"` | General UEL path for mechanics, scalar transport, and coupled problems |
| `"fbar_mechanics"` | Pure mechanical low-order F-bar UELs |
| `"local_pressure"` | Prototype Quad4/Hex8 `u,mu` path with one condensed element-local pressure |

Do not use F-bar as the default route for coupled gel, diffusion, or thermal
problems. Use a mixed field or a standard coupled UEL unless the theory
explicitly justifies F-bar averaging.

For generated local-pressure UELs, `SVARS(1)` is the condensed pressure state.
When `Variables >= 1 + 2*NGP`, the generator also mirrors diagnostics as
`SVARS(1+i)=phi_i` and `SVARS(1+NGP+i)=p_i`. Abaqus `UVARM1/UVARM2` expose the
same `phi` and pressure diagnostics through the generated UVARM bridge.

### Abaqus Input Scaffolding

For converting an Abaqus deck so selected elements use a generated UEL:

```python
config = au.UELModelConfig(...)
report = au.generate_inp_scaffold(config, "source.inp", "output_dir")
au.write_job_inp("output_dir")
```

The signatures are `generate_inp_scaffold(config, mesh_path, out_dir, dry_run=False)`
(config first) and `write_job_inp(out_dir, heading=..., mesh_include=...)`, which
writes the top-level include driver into `out_dir`.

The scaffold path is intentionally conservative. It helps split mesh, material,
UEL, step, and output includes, but the user still owns Abaqus modeling choices
such as boundary conditions, element sets, output requests, hourglass controls,
and `UNSYMM`.

## Tensor DSL Rules

The translator is a 3x3 continuum-mechanics DSL, not a general NumPy compiler.

Prefer:

```python
F.T @ F
trace(A)
dev(A)
0.5*(A + A.T)
det(F)
inv(F)
eig(A)
logm(A)
expm(A)
```

Prefer explicit scalar/tensor temporaries in generated finite-strain UMAT/UEL
methods when a tensor expression contains matrix-to-scalar calls such as
`det(F)`.

Prefer:

```python
def stress_PK1(self, F):
    J = det(F)
    FinvT = inv(F).T
    return self.G * (F - FinvT) + self.K * log(J) * FinvT
```

Avoid writing the same formula as one direct tensor return:

```python
def stress_PK1(self, F):
    return self.G * (F - inv(F).T) + self.K * log(det(F)) * inv(F).T
```

The second form is mathematically fine. It is now covered by regression tests,
but it historically exposed a finite-strain UMAT/UEL code-generation bug while
converting tensor returns into component-wise Fortran. The confirmed symptom
was invalid generated code such as
`det33z(F(ii,jj))`: the generator has indexed the full tensor `F` inside the
scalar determinant call. The safe temporary style forces `J = det33z(F)` to be
computed once as a scalar before the component loop, which is easier to inspect.

Avoid in generated methods:

- dynamic list/dict construction inside the hot update,
- `while` loops,
- data-dependent array sizes,
- `np.einsum`, `np.tensordot`, general dynamic-size `np.linalg.solve`,
- arbitrary or ragged `np.array([[...]])` construction,
- unguarded `abs`, `max`, and sign switches when derivatives matter.

Supported small-array exceptions:

- `np.array([...])` literal vectors and square literal matrices with static
  sizes can be translated.
- `np.linalg.solve(A, b)` can be translated for statically inferred small
  systems with N in `{2, 3, 4, 5, 6, 9}`.

Use bounded `for k in range(N): ... break` for local iteration. Use helper
methods for readability, but keep signatures flat: pass tensors and scalars,
not dataclasses.

## Verification Ladder

Use this order for every new model:

1. Document the theory, non-scope, conventions, fields/DOFs, and state layout.
2. Add an independent quantitative oracle and model-specific regime checks.
3. Run material-method tangent consistency with `verify()`.
4. For a UEL, check assembled `RHS`/`AMATRX`, DOF/state layout, and an
   appropriate patch or invariant.
5. Regenerate deterministically and compile the Fortran with `gfortran`.
6. Directly call nontrivial generated code through f2py or another checked
   compiled runtime.
7. Add a small solver run only when useful, then validate its output bridge.
8. Attempt a published benchmark or collaborator comparison last.

`Material.verify()` checks the tangent of implemented material methods at its
chosen state. `WeakForm.verify()` currently delegates to material verification;
it does not assemble a UEL residual/tangent. Neither proves that the equations
match the paper, that all branches were entered, or that a gradient coefficient
matches the corresponding reaction coefficient. For branchy models, add tests
that assert the regime was reached and measure branch-specific invariants.

Do not interpret a failed solver run before the earlier rungs pass. Solver
nonconvergence can come from the material law, load step, tangent symmetry,
boundary conditions, hourglass control, or the solver path.

### Solver-output bridge

Abaqus setup is example/user-owned. When a solver result is used as evidence,
the example must still define and test its output bridge:

- map each logical quantity to its field/component/slot/active DOF;
- declare units, signs/offsets, component and integration-point order;
- require stable identity, exact coverage, uniqueness, and finiteness;
- distinguish authoritative solved fields from reconstructed or
  visualization-only bridges; and
- audit a `UVARM`, dummy-element, or projected bridge pointwise where possible.

The shared ODB extractor is intentionally for simple one-element cases. Keep
model-specific histories, projections, named-set reductions, and derived
observables in the example.

## Reference Examples

The public positive manifest is [`../examples/README.md`](../examples/README.md).
Use only folders present there and in this checkout. The working
[`../examples/_template/`](../examples/_template/) demonstrates the complete
minimal UMAT pipeline; it is not a scientific capability claim.
