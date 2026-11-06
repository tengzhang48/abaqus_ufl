# Optional FE runtimes

`abaqus_ufl.fe` assembles small verification meshes and solves their free DOFs
with SciPy. Its compiled adapter calls the generated Fortran `UEL` through
f2py. A Fortran compiler is required; Abaqus is not.

```bash
pip install -e ".[fe]"
python tools/check_fe_runtime.py
```

The supported compiled path is **2D standard Quad4**, uniform degree-one
nodal fields, and interleaved node/component DOFs. Previous nodal values are
passed through `DU`; transient transport and explicit material `SVARS`
commits are supported. Condensed local fields, mixed interpolation, F-bar,
other topologies, contact and dynamics are outside this API. `newton_solve`
solves one increment; the serial [stepping drivers](#stepping-drivers-and-bounded-solves)
add external loads, automatic time cutbacks, staggered and bounded solves.
PETSc and MPI are optional hosts for this same element contract.
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

## Exact batches and material state

The scalar path remains the default. To evaluate all elements in one f2py
call, pass `batch_fn=element.element_rk_batch` to `newton_solve`, or call
`assemble_batched` with the same explicit `DU` contract as `assemble`.

Construct `CompiledAbaqusElement(..., cache_batches=True)` to reuse the last
successful full batch. The default is `False`. This is safe only when the
kernel's complete result is determined by its arguments and no external
mutable state or required per-call side effect is involved. Exact snapshots
include geometry, U, DU, committed SVARS, properties, IDs, flags, parameters,
time and increment counters; module identity is also checked. Returned arrays
cannot modify the saved result. `batch_evaluations` counts actual attempts;
`batch_cache_hits` counts reuse. Scalar calls are never cached.

For a material declaring `state_vars`, supply `n_elem=len(elems)` to allocate
committed state. Every evaluation starts from that state and discards its
trial updates. Only after convergence, call
`element.commit_step(nodes, elems, dof_per_node, U, U_previous)`. An exact
cached trial may supply the commit without another kernel call. Commit,
`begin_increment`, `reject_step` and `clear_batch_cache` release the cache;
rejected/nonfinite kernel output is never cached or committed. A failed
Newton solve clears the bound adapter's cache. Custom solver rejection or
restart code should call `reject_step()` as well.

Generated UELs initialize SVARS from the material declaration when `TIME(2)=0`.
For a manually seeded state, use the appropriate nonzero beginning time so
that initialization does not replace it. Retain the accepted nodal solution
separately and set the beginning-of-increment clock before each step.

## Stepping drivers and bounded solves

`abaqus_ufl.fe.stepping` and `abaqus_ufl.fe.bounds` advance many increments
serially. They take the same element callables (`element.element_rk`, optional
`batch_fn=element.element_rk_batch`, `commit_fn=element.commit_step`) and
`dirichlet_fn(frac)`, which returns prescribed global DOF values at the load
fraction `frac`.

| Function | What it does |
|---|---|
| `solve_steps` | One monolithic Newton with line search per load step; optional `external_fn` load and `step_callback`; no cutback |
| `solve_steps_staggered` | Alternating Newton over disjoint per-node component groups (one group of all components is monolithic), with per-field scaled residual and update gates, adaptive cutbacks (`max_cutbacks`), restart dictionaries and a `StaggeredReport`; optional `field_bounds` clip or `BoxSpec`, `feasible_trials` and `newton_method="bounded_monolithic"` |
| `assemble_groups`, `solve_steps_multimaterial` | Element groups that touch different components of one uniform per-node layout |
| `solve_steps_bounded` | Projected Newton (primal-dual active set) inside a `BoxSpec` box, with optional `irreversible` or `monotone_decreasing` components, cutbacks and a `BoundedReport` |

**A step is committed only below `tol`.** The test is the free-DOF residual
2-norm (`solve_steps`, multimaterial), the per-field scaled residual and update
norms (staggered), or the weighted projected (KKT) residual (bounded). Weights
come from `residual_scale`; tolerances are not relative to the first residual.
The state produced by the last permitted Newton update is assessed. There is
no looser fallback, and a line search that cannot reduce the residual in 30
bisections fails the step instead of applying a vanishing update. A clip
or projection is followed by a fresh projected-residual certificate.

**Bounds.** `BoxSpec` and `field_bounds` must name existing components
(integer indices below `dof_per_node`), or the solve raises `ValueError`. On
every bounded route a prescribed value outside the bounds by more than 1e-12
raises `ValueError`; clipping and snapping never change a prescribed value.

**Time.** Pass `state=CompiledElementState([element])`. Before each attempt
the cutback drivers call `state.begin_increment(kinc, time_begin=t, dt=dt_sub)`:
`TIME` is the beginning of the attempted increment (accumulated accepted time,
starting at 0), `DTIME` is that attempt's own length, reduced on a cutback, and
`KINC` counts accepted increments from 1; a retry keeps its `KINC`.
`solve_steps` and `solve_steps_multimaterial` call
`state.begin_increment(kinc=step)`, so `TIME = (step - 1) * element.dt`. At
`TIME(2) = 0` a generated UEL initializes its SVARS from the declaration.

**Transactions.** SVARS and the element clock are saved before each attempt.
A rejected attempt restores both and discards any cached trial batch before
`attempt_callback` sees its record. An attempt is rejected by
`IncrementTooLarge` (`PNEWDT < 1`) at an iterate or during commit, by an
exhausted Newton or line search, or by a failed certificate. At a line-search
trial, `IncrementTooLarge` also rejects the attempt in unprojected staggered
blocks; projected solves backtrack past it unless
`cutback_on_element_rejection=True` (staggered driver). A nonfinite
line-search trial counts as an infinite residual. `NonfiniteElementOutput` at
an iterate rejects a staggered attempt and stops `solve_steps_bounded`. Any
other exception from an element or commit propagates, after SVARS and the
clock have been restored.

**Restart** (`solve_steps_staggered` only). Its `accepted_callback` receives a
restart dictionary (`U`, state snapshot, step cursor, `time_end`, `next_kinc`,
report) to pass back as `restart=`. `solve_steps_bounded` and the other
drivers have no restart parameters, and file checkpoints are not provided.

The stepping assembler accepts repeated node indices in an element, such as a
collapsed Quad4 at a polar centre, and scatter-adds them; the element must
still be admissible. Public `assemble` keeps rejecting such connectivity.
`newton_method="bounded_monolithic"` classifies bound activity exactly, so
`bound_atol` applies only to `feasible_trials` blocks and `solve_steps_bounded`.

```python
from abaqus_ufl.fe import CompiledElementState, solve_steps_staggered

U, report = solve_steps_staggered(
    nodes, elems, element.dof_per_node, element.element_rk, dirichlet_fn,
    [[0, 1], [2]], n_steps, batch_fn=element.element_rk_batch,
    commit_fn=element.commit_step, state=CompiledElementState([element]),
    dt=dt, tol=1e-8, max_cutbacks=4, return_report=True)
```

## PETSc and MPI

The package and MPI modules import without `petsc4py` or `mpi4py`. PETSc is
loaded only when requested; these hosts use PETSc collectives and do not
require mpi4py. Use a launcher that matches the MPI library linked to PETSc.
No MPI vendor or blanket ban on mpi4py is part of the framework contract.

Three interfaces serve different sizes of problem:

| Interface | Distribution | Input and return |
|---|---|---|
| `make_petsc_solver` / `make_fieldsplit_solver` | Single-rank communicator, default `COMM_SELF` | Serial CSR matrix and full RHS; full solution |
| `MPIBatchedElement` | Partitioned element evaluation; replicated mesh, state and assembly; rank-zero linear solve | Full batches on every rank; full R/K and solution |
| `solve_steps_mpi` / `solve_steps_mpi_local` | Distributed PETSc matrix and owned/ghost solution | Replicated mesh convenience interface, or local partitions; final full solution gathered on every rank |

For a serial CSR solve, pass `linear_solve=make_petsc_solver(pc="lu",
solver="mumps")` to `newton_solve`. Iterative PCs such as `gamg` default to
GMRES for nonsymmetric coupled tangents. `make_fieldsplit_solver` takes named
per-node component groups, for example `[("u", [0, 1]), ("T", [2])]`.
It preserves original component identities when `newton_solve` eliminates
individual constrained DOFs. Custom reduced-system callers should use
`solve.solve_reduced(K_free, b_free, free_dofs)`, where `free_dofs` gives the
original global DOF index of each reduced row. Calling `solve(K, b)` directly
assumes a full interleaved layout. Fully constrained fields are omitted; if
only one field remains, its `sub_pc` acts on that remaining system.
MUMPS, hypre and other external PETSc solver packages must exist in the local
PETSc build when selected; their failures are reported, never treated as a
converged solve.

The public parallel example compiles once on rank zero into a shared folder,
then checks two-increment mechanics and heat diffusion against serial and
independent affine/discrete-eigenmode oracles:

```bash
mpiexec -n 2 python examples/scalar_diffusion_uel/check_parallel.py
```

Run from the repository root on a shared filesystem. The default
`.parallel-build/` directory retains source and logs. Use `--build-dir` to
choose another shared build directory, and use a distinct directory for each
simultaneous job. The MPI batch adapter's
`root_linear_solver()` is a collective callable; every rank must enter it.
Batch input conversion and validation errors are voted collectively before
any rank raises, including preparation for `commit_group` and `commit_step`.
Do not use `MPIBatchedElement` as the local evaluator of the distributed host.

`solve_steps_mpi(..., begin_step=callback)` calls `callback(step, fraction)`
on every rank before assembly; use it to call `element.begin_increment(step)`.
It takes global element indices and an optional full-mesh `commit_fn`.
`solve_steps_mpi_ghost` is a compatibility spelling for this ghost-vector host.
`solve_steps_mpi_local` takes `(ndof, my_gm, my_coords, dof_per_node, ...)`,
local element indices, optional local `batch_fn` and `commit_element_fn`.
`partition_structured_quad` constructs local structured meshes directly;
`partition_unstructured` returns a node permutation and ownership count for
local assembly. Preserve original global one-based `elem_ids` on each compiled
adapter so JELEM does not collide between partitions.

The distributed host checks KSP convergence, line search, finite output and
Newton convergence before commit. Errors and cutbacks are voted collectively
before data gathers. Bound compiled adapters' SVARS are restored if any rank
rejects a commit; arbitrary external side effects in custom callbacks require
their own transaction. The distributed host has no adaptive cutback,
bounded/staggered solve or restart service; those are serial only. Mesh convenience inputs, final output, and full-mesh
commit callbacks are replicated; this is not an end-to-end scaling claim.

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
