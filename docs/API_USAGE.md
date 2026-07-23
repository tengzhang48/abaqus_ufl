# abaqus_ufl API Usage Guide

Status: 2026-06-27

This note is the compact user-facing API reference for writing models with
`abaqus_ufl`. It focuses on the stable workflows that have been exercised by
the current UMAT, UEL, f2py, and Abaqus-validation examples.

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
au.generate_uinter(...)

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

    def heat_storage(self, F, F_old, T, T_old, dt):
        return self.rho_cp * (T - T_old) / dt

    def heat_flux(self, F, T, grad_T):
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
        return self.material.heat_storage(F, F_old, T, T_old, dt), self.material.heat_flux(F, T, grad_T)

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

When a scalar gradient model has mechanism-dependent coefficients, thread the
same branch data through both the storage/reaction term and the gradient term.
For example, a damage law with tensile and compressive mechanisms should not use
`psi_star_damage` in `phase_storage` but a hard-coded brittle `psi_star` in
`phase_flux`. Pass the state needed to choose the mechanism into the flux method:

```python
def phase_storage(self, F, damage, damage_old, Fp_old, dt):
    H, psi_star_m, eta_m = self.damage_params(F, Fp_old)
    return eta_m * (damage - damage_old) / dt - 2.0 * (1.0 - damage) * H + 2.0 * psi_star_m * damage

def phase_flux(self, F, grad_damage, Fp_old):
    _, psi_star_m, _ = self.damage_params(F, Fp_old)
    return -2.0 * psi_star_m * self.ell * self.ell * grad_damage
```

This adds `dflux/dF` blocks to the full coupled UEL tangent. If a production
run needs a diagonal or staggered tangent, use `drop_tangent_coupling` or a
separate diagonal variant deliberately; do not make the physics branch-agnostic
just to keep the tangent smaller.

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
| `"eig"` | Eig-based `sqrtm33z`, `logm33z`, `expm33z`, `polar33z`. Not CS-safe at diagonal-F states due to `eig33z` near-diagonal guard; available for backward compatibility. |

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
report = au.generate_inp_scaffold("source.inp", "output_dir", config)
au.write_job_inp("output_dir/job.inp", report)
```

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
sym(A)
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

1. Python material or weak-form verification with `verify()`.
2. Material-point or element-level regime checks that prove the tested states
   enter every claimed branch.
3. Generated Fortran compile check with `gfortran`.
4. Optional f2py material-point or element-level check.
5. Abaqus one-element validation.
6. Published benchmark or collaborator comparison.

`verify()` checks the tangent of the implemented residual at its chosen state.
It does not prove that the residual matches the paper, that all branches were
entered, or that a gradient coefficient matches the corresponding reaction
coefficient. For branchy damage/plasticity models, add tests that assert the
regime was reached and measure branch-specific invariants such as an effective
regularization length from the ratio of gradient and reaction stiffnesses.

Do not interpret a failed solver run before the earlier rungs pass. Solver
nonconvergence can come from the material law, load step, tangent symmetry,
boundary conditions, hourglass control, or the solver path.

## Reference Examples

Use these examples as templates:

| Need | Start from |
|------|------------|
| finite-strain elastic UMAT | `examples/neo_hookean_umat`, `examples/mooney_rivlin_umat` |
| viscoelastic UMAT | `examples/small_strain_viscoelastic_umat` |
| small-strain plastic UMAT | `examples/small_strain_j2_umat` |
| scalar UEL | `examples/scalar_diffusion_uel` |
| phase-field UEL | `examples/phasefield_fracture_uel` |
| coupled displacement-temperature UEL | `examples/thermo_mechanics_quad8` |
| Abaqus deck validation | `examples/_template`, `tools/` |
