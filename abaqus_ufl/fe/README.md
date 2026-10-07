# Optional serial FE runtime

`abaqus_ufl.fe` assembles small verification meshes and solves their free DOFs
with SciPy. Its compiled adapter calls the generated Fortran `UEL` through
f2py. A Fortran compiler is required; Abaqus is not.

```bash
pip install -e ".[fe]"
python tools/check_fe_runtime.py
```

The supported compiled path is deliberately small: **2D standard Quad4**,
uniform degree-one nodal fields, interleaved node/component DOFs, and
stateless materials. Previous nodal values are passed through `DU`; transient
transport is supported. Material `SVARS`, condensed local fields, mixed
interpolation, F-bar, other topologies, contact, dynamics, external force
assembly, automatic time cutbacks, PETSc, and MPI are outside this API.
Unsupported layouts are rejected before compilation.

## API

```python
from abaqus_ufl.fe import (
    build_compiled_uel, CompiledAbaqusElement,
    structured_quad_mesh, newton_solve,
)

# problem is an instantiated, supported WeakForm.
module = build_compiled_uel(problem)
element = CompiledAbaqusElement(module, problem, dt=0.1)
nodes, elems = structured_quad_mesh(4, 2)

# conditions maps global DOF indices to prescribed values. Components at
# node i occupy i*dof_per_node + component. For the public heat example,
# the components are (ux, uy, T).
U = newton_solve(nodes, elems, element.dof_per_node, element.element_rk,
                 conditions, U0=initial, U_previous=previous)
```

`build_compiled_uel` generates and compiles once per module. The f2py wrapper
is package data, so installed wheels do not depend on a checkout's tests
directory. Supply `workdir=` to retain the generated source and build log.
Otherwise the returned module owns its temporary build directory. Unique
module names prevent a prior build from being imported accidentally.

`CompiledAbaqusElement` uses the material instance's property values,
including constructor overrides. It requests both residual and tangent and
converts Abaqus's `RHS` into `R=-RHS`, retaining `K=AMATRX=dR/dU`. A UEL
cutback, nonfinite output, singular global matrix, or unconverged Newton
solve raises an error. No rejected trial is returned as a solution.
Newton's `tol` is relative to the larger of 1 and the largest residual entry, reactions included.

For successive steps, retain the accepted `U` as `U_previous`, call
`element.begin_increment(kinc, time_begin=previous_time, dt=dt)`, and solve
with the next boundary values. `TIME` is the beginning of the increment,
matching the generated UEL's convention; the runtime uses one analysis step.
The previous solution remains fixed during Newton and line search.

Direct `assemble` calls must supply `DU=U-U_previous`. For an explicitly
from-zero call, use `du_is_total=True`. There is no silent history default.

`newton_solve(..., record_iteration=callback)` optionally reports the accepted
iterate's iteration number, free residual infinity norm, and scaled relative
residual, including the converged iterate. Rejected line-search trials are not
recorded as accepted steps.

## Verification

`tools/check_fe_runtime.py` compiles the public thermo-mechanical Quad4
declaration and checks:

- an affine mechanics patch with solved interior displacement DOFs;
- transient heat diffusion on 4x2, 8x2, and 16x2 meshes against the independent
  Q1 consistent-mass eigenmode formula;
- second-order spatial refinement against the continuum backward-Euler
  eigenmode, which is distinct from an exact continuous-time solution.

The livebench records these numerical errors alongside the six existing
material/element bundles. These mesh checks establish the stated host and
element contracts; they do not reproduce the paper's full Abaqus examples.

The livebench also runs actual multi-increment boundary-value problems:
`tools/check_heated_plate.py` solves spatially two-dimensional heating with
prescribed edge temperatures; `tools/check_thermal_bending.py` solves
thermally driven bending with a clamped end and traction-free remaining
mechanical edges. Each compiles once, advances accepted nodal history,
retains the solved fields, and checks independent continuum or reduced
references. See [the BVP definitions](../../docs/BOUNDARY_VALUE_SIMULATIONS.md).

This project-authored serial runtime and its wrapper use the repository's
MIT license. It does not include the separately attributed research host
described in `CREDITS.md`.
