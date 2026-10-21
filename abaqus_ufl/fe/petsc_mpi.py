"""Optional distributed PETSc assembly and Newton solves with ghost exchange.

Adapted from the generic MPI host: element partitions insert off-rank entries
into MPIAIJ matrices. Only owned and needed ghost values enter each Newton
iteration. The final result is gathered for verification. The convenience
mesh interface retains replicated mesh inputs; the local interface does not.
No MPI implementation is prescribed; the launcher and libraries must match.
"""

import time

import numpy as np

from ._collective import collective_status, gather
from .compiled_element import IncrementTooLarge
from .driver import _mesh_inputs
from .petsc_backend import _field_components, _petsc


def _positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(name + " must be a positive integer")


def partition_structured_quad(nx, ny, dof_per_node, rank, size, Lx=1.0, Ly=1.0):
    """Return ``(ndof, local_global_dofs, local_coords)`` without a full mesh.

    Numbering matches ``structured_quad_mesh``. Empty element partitions are
    supported, including runs with more ranks than elements.
    """
    for value, name in ((nx, 'nx'), (ny, 'ny'), (dof_per_node, 'dof_per_node'), (size, 'size')):
        _positive_int(value, name)
    if not 0 <= rank < size or not np.isfinite([Lx, Ly]).all() or min(Lx, Ly) <= 0:
        raise ValueError("invalid partition rank or mesh dimensions")
    ne = nx * ny
    first, last = ne * rank // size, ne * (rank + 1) // size
    e = np.arange(first, last)
    j, i = np.divmod(e, nx)
    n0 = j * (nx + 1) + i
    ids = np.column_stack((n0, n0 + 1, n0 + nx + 2, n0 + nx + 1))
    coords = np.stack(((ids % (nx + 1)) * (Lx / nx), (ids // (nx + 1)) * (Ly / ny)), axis=2)
    dofs = (ids[:, :, None] * dof_per_node + np.arange(dof_per_node)).reshape(len(e), 4*dof_per_node)
    return (nx + 1) * (ny + 1) * dof_per_node, dofs, coords


def partition_rcb(centroids, size):
    """Balanced recursive coordinate bisection, including empty partitions."""
    _positive_int(size, 'size')
    centroids = np.asarray(centroids, float)
    if centroids.ndim != 2 or centroids.shape[1] < 1 or not np.isfinite(centroids).all():
        raise ValueError("centroids must be a finite (N, dimension) array")
    owners = np.empty(len(centroids), dtype=np.int32)

    def split(indices, first, count):
        if not len(indices):
            return
        if count == 1:
            owners[indices] = first
            return
        points = centroids[indices]
        axis = int(np.argmax(np.ptp(points, axis=0)))
        order = indices[np.argsort(points[:, axis], kind='stable')]
        left = count // 2
        k = int(round(len(order) * left / count))
        split(order[:k], first, left)
        split(order[k:], first + left, count - left)

    split(np.arange(len(centroids)), 0, size)
    return owners


def partition_unstructured(nodes, elems, dof_per_node, size, rank, elem_rank=None):
    """Return a local Quad4 mesh with contiguous node ownership and renumbering.

    Returns ``(ndof, local_dofs, local_coords, local_ndof, old_to_new_nodes)``.
    Translate boundary conditions and unpermute the returned solution with the
    node permutation. Preserve original global element IDs separately for UELs.
    """
    nodes, elems, ndof = _mesh_inputs(nodes, elems, dof_per_node)
    _positive_int(size, 'size')
    if not 0 <= rank < size:
        raise ValueError("rank must be within the communicator")
    owners = partition_rcb(nodes[elems].mean(axis=1), size) if elem_rank is None else np.asarray(elem_rank)
    if (owners.shape != (len(elems),) or not np.issubdtype(owners.dtype, np.integer)
            or np.any(owners < 0) or np.any(owners >= size)):
        raise ValueError("elem_rank must assign every element to a valid rank")
    votes = np.zeros((len(nodes), size), dtype=np.int64)
    np.add.at(votes, (elems.ravel(), np.repeat(owners, 4)), 1)
    node_rank = votes.argmax(axis=1)
    order = np.argsort(node_rank, kind='stable')
    old_to_new = np.empty(len(nodes), dtype=np.int64)
    old_to_new[order] = np.arange(len(nodes))
    mine = elems[owners == rank]
    dofs = old_to_new[mine]
    dofs = (dofs[:, :, None] * dof_per_node + np.arange(dof_per_node)).reshape(len(mine), 4*dof_per_node)
    return ndof, dofs, nodes[mine], int((node_rank == rank).sum()) * dof_per_node, old_to_new


def _reject(functions):
    for fn in functions:
        owner = getattr(fn, '__self__', None)
        if hasattr(owner, 'reject_step'):
            owner.reject_step()


def _commit(callback, args, comm):
    """Restore a compiled adapter's old state if any partition rejects commit."""
    owner = getattr(callback, '__self__', None)
    old = getattr(owner, 'svars', None)
    old = old.copy() if old is not None else None
    error = None
    try:
        callback(*args)
    except Exception as exc:
        error = exc
    try:
        collective_status(error, getattr(error, 'pnewdt', 1.0), comm)
    except Exception:
        if old is not None:
            owner.svars = old
        _reject([callback])
        raise


def solve_steps_mpi(nodes, elems, dof_per_node, element_fn, dirichlet_fn, n_steps,
                    max_newton=40, tol=1e-7, line_search=True, commit_fn=None,
                    pc="gamg", ksp_type="gmres", solver="mumps", rtol=1e-10,
                    verbose=False, *, comm=None, begin_step=None, u0=None):
    """Distributed assembly/solve from a replicated Quad4 mesh.

    ``element_fn`` takes global element indices. ``commit_fn`` receives the
    full accepted nodal solution on every rank once per step. No state is
    committed after a failed solve. ``begin_step(step, fraction)`` may set the
    adapter's increment clock. This host does not perform automatic cutbacks.
    """
    PETSc = _petsc()
    comm = PETSc.COMM_WORLD if comm is None else comm
    error = None
    try:
        nodes, elems, ndof = _mesh_inputs(nodes, elems, dof_per_node)
    except Exception as exc:
        error = exc
    collective_status(error, comm=comm)
    rank, size = comm.getRank(), comm.getSize()
    first, last = len(elems) * rank // size, len(elems) * (rank + 1) // size
    mine = elems[first:last]
    gm = (mine[:, :, None] * dof_per_node + np.arange(dof_per_node)).reshape(len(mine), 4*dof_per_node)

    def local(ei, coords, U, DU):
        return element_fn(first + ei, coords, U, DU)

    def commit(U, previous):
        _commit(commit_fn, (nodes, elems, dof_per_node, U, previous), comm)

    try:
        return solve_steps_mpi_local(
            ndof, gm, nodes[mine], dof_per_node, local, dirichlet_fn, n_steps,
            max_newton=max_newton, tol=tol, line_search=line_search,
            pc=pc, ksp_type=ksp_type, solver=solver, rtol=rtol, verbose=verbose,
            comm=comm, begin_step=begin_step, u0=u0,
            _commit_full=commit if commit_fn is not None else None)
    except Exception:
        _reject([element_fn, commit_fn])
        raise


def solve_steps_mpi_ghost(*args, **kwargs):
    """Compatibility spelling: the mesh host also uses ghost vectors."""
    return solve_steps_mpi(*args, **kwargs)


def solve_steps_mpi_local(ndof, my_gm, my_coords, dof_per_node, element_fn,
                          dirichlet_fn, n_steps, max_newton=40, tol=1e-7,
                          line_search=True, commit_element_fn=None, pc="gamg",
                          ksp_type="gmres", solver="mumps", rtol=1e-10,
                          verbose=False, local_ndof=None, batch_fn=None,
                          field_components=None, fieldsplit_type="additive",
                          fieldsplit_sub_pc="gamg", u0=None, *, comm=None,
                          begin_step=None, _commit_full=None):
    """Solve using local elements and distributed owned/ghost solution values.

    Callbacks take local element indices. ``batch_fn`` evaluates this rank's
    entire local partition; use a regular CompiledAbaqusElement, with global
    ``elem_ids``, rather than MPIBatchedElement. The final solution is gathered
    on all ranks. Generic commit callbacks must supply their own transaction
    for external side effects; bound compiled adapters' SVARS are rolled back
    collectively when commit fails. Every rank must enter the same calls.
    """
    PETSc = _petsc()
    comm = PETSc.COMM_WORLD if comm is None else comm
    rank, size = comm.getRank(), comm.getSize()
    error, fields = None, None
    try:
        for value, name in ((ndof, 'ndof'), (dof_per_node, 'dof_per_node'),
                            (n_steps, 'n_steps'), (max_newton, 'max_newton')):
            _positive_int(value, name)
        if not np.isfinite([tol, rtol]).all() or min(tol, rtol) <= 0:
            raise ValueError("Newton and KSP tolerances must be positive and finite")
        my_gm, my_coords = np.asarray(my_gm), np.asarray(my_coords, float)
        if (my_gm.ndim != 2 or not np.issubdtype(my_gm.dtype, np.integer)
                or my_gm.shape[1] != 4*dof_per_node
                or np.any(my_gm < 0) or np.any(my_gm >= ndof)
                or my_coords.shape != (len(my_gm), 4, 2) or not np.isfinite(my_coords).all()):
            raise ValueError("local Quad4 DOF and coordinate arrays are invalid")
        if pc == 'fieldsplit':
            fields = _field_components([(str(i), c) for i, c in enumerate(field_components or [])], dof_per_node)
            if fieldsplit_type not in ('additive', 'multiplicative', 'schur'):
                raise ValueError("invalid FieldSplit type")
            if fieldsplit_type == 'schur' and len(fields) != 2:
                raise ValueError("Schur splitting requires two fields")
    except Exception as exc:
        error = exc
    collective_status(error, comm=comm)
    objects = []

    def keep(obj):
        objects.append(obj)
        return obj

    def full(vec):
        return gather(vec.getArray(readonly=True), comm)

    try:
        vsize = ndof if local_ndof is None else (local_ndof, ndof)
        U = keep(PETSc.Vec().createMPI(vsize, comm=comm))
        U.set(0.0)
        rs, re = U.getOwnershipRange()
        msize = [(re-rs, ndof), (re-rs, ndof)]
        error = None
        try:
            if u0 is not None:
                initial = np.asarray(u0, float)
                if initial.shape == (ndof,):
                    initial = initial[rs:re]
                if initial.shape != (re-rs,) or not np.isfinite(initial).all():
                    raise ValueError("u0 must contain finite owned or global DOF values")
                U.getArray()[:] = initial
        except Exception as exc:
            error = exc
        collective_status(error, comm=comm)
        needed = np.unique(my_gm).astype(PETSc.IntType)
        local_gm = np.searchsorted(needed, my_gm)  # O(local DOFs), no global lookup array
        needed_is = keep(PETSc.IS().createGeneral(needed, comm=comm))
        Uloc = keep(PETSc.Vec().createSeq(len(needed), comm=PETSc.COMM_SELF))
        Ploc = keep(Uloc.duplicate())
        scatter = keep(PETSc.Scatter().create(U, needed_is, Uloc, None))
        previous = keep(U.duplicate())
        rows = my_gm.astype(PETSc.IntType)
        owned = np.arange(rs, re, dtype=PETSc.IntType)

        def ghost(vec, into):
            scatter.scatter(vec, into, addv=PETSc.InsertMode.INSERT_VALUES,
                            mode=PETSc.ScatterMode.FORWARD)
            return into.getArray(readonly=True)

        def assemble(current):
            values, old = ghost(current, Uloc), ghost(previous, Ploc)
            ue, due = values[local_gm], values[local_gm] - old[local_gm]
            error = None
            try:
                if batch_fn is not None:
                    r, k = (np.asarray(a, float) for a in batch_fn(my_coords, ue, due))
                else:
                    r = np.empty_like(ue)
                    k = np.empty((len(ue), ue.shape[1], ue.shape[1]))
                    for ei in range(len(ue)):
                        r[ei], k[ei] = element_fn(ei, my_coords[ei], ue[ei], due[ei])
                if r.shape != ue.shape or k.shape != (len(ue), ue.shape[1], ue.shape[1]):
                    raise ValueError("element batch output shapes disagree with local DOFs")
                if not np.isfinite(r).all() or not np.isfinite(k).all():
                    raise RuntimeError("nonfinite element residual or tangent")
            except Exception as exc:
                error = exc
            collective_status(error, getattr(error, 'pnewdt', 1.0), comm)
            A = PETSc.Mat().createAIJ(msize, nnz=dof_per_node * 18, comm=comm)
            R = U.duplicate()
            try:
                A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
                R.set(0.0)
                for ei in range(len(ue)):
                    A.setValues(rows[ei], rows[ei], np.ascontiguousarray(k[ei]), addv=PETSc.InsertMode.ADD_VALUES)
                    R.setValues(rows[ei], r[ei], addv=PETSc.InsertMode.ADD_VALUES)
                A.assemble()
                R.assemble()
                return A, R
            except Exception:
                A.destroy()
                R.destroy()
                raise

        def constrain(R, bc_owned):
            scale = max(1.0, R.norm(PETSc.NormType.NORM_INFINITY))
            R.setValues(bc_owned, np.zeros(len(bc_owned)), addv=PETSc.InsertMode.INSERT_VALUES)
            R.assemble()
            norm = R.norm(PETSc.NormType.NORM_INFINITY)
            if not np.isfinite([norm, scale]).all():
                raise RuntimeError("nonfinite assembled residual")
            return norm, scale

        t_asm, ksp_its = 0.0, 0
        for step in range(1, n_steps + 1):
            U.copy(previous)
            error = None
            try:
                if begin_step is not None:
                    begin_step(step, step / n_steps)
                bc = dirichlet_fn(step / n_steps)
                if any(isinstance(i, bool) or not isinstance(i, (int, np.integer))
                       or i < 0 or i >= ndof or not np.isfinite(v) for i, v in bc.items()):
                    raise ValueError("invalid Dirichlet conditions")
                bc_owned = owned[np.isin(owned, list(bc))]
                U.getArray()[bc_owned - rs] = [bc[int(i)] for i in bc_owned]
            except Exception as exc:
                error = exc
            collective_status(error, comm=comm)
            converged = False
            for iteration in range(max_newton + 1):
                started = time.monotonic()
                A, R = assemble(U)
                t_asm += time.monotonic() - started
                work = [A, R]
                try:
                    norm, scale = constrain(R, bc_owned)
                    if norm <= tol * scale:
                        converged = True
                        break
                    if iteration == max_newton:
                        raise RuntimeError("MPI Newton solve did not converge")
                    A.zeroRows(bc_owned, diag=1.0)
                    rhs, delta = R.copy(), R.duplicate()
                    work.extend([rhs, delta])
                    rhs.scale(-1)
                    ksp = PETSc.KSP().create(comm)
                    work.append(ksp)
                    ksp.setOperators(A)
                    ksp.setType(ksp_type)
                    preconditioner = ksp.getPC()
                    if fields is None:
                        preconditioner.setType(pc)
                        if pc == 'lu' and solver:
                            preconditioner.setFactorSolverType(solver)
                    else:
                        preconditioner.setType('fieldsplit')
                        types = {'additive': PETSc.PC.CompositeType.ADDITIVE,
                                 'multiplicative': PETSc.PC.CompositeType.MULTIPLICATIVE,
                                 'schur': PETSc.PC.CompositeType.SCHUR}
                        preconditioner.setFieldSplitType(types[fieldsplit_type])
                        sets = [(name, PETSc.IS().createGeneral(
                            owned[np.isin(owned % dof_per_node, comps)], comm=comm)) for name, comps in fields]
                        work.extend(s for _, s in sets)
                        preconditioner.setFieldSplitIS(*sets)
                        preconditioner.setUp()
                        for sub in preconditioner.getFieldSplitSubKSP():
                            sub.getPC().setType(fieldsplit_sub_pc)
                    ksp.setTolerances(rtol=rtol, max_it=2000)
                    ksp.solve(rhs, delta)
                    ksp_its = ksp.getIterationNumber()
                    if ksp.getConvergedReason() <= 0 or not np.isfinite(delta.norm()):
                        raise RuntimeError("MPI PETSc KSP failed (reason {})".format(ksp.getConvergedReason()))
                    alpha = 1.0
                    if line_search:
                        trial = U.duplicate()
                        work.append(trial)
                        for search in range(12):
                            U.copy(trial)
                            trial.axpy(alpha, delta)
                            # Preserve prescribed values exactly even after an iterative solve.
                            trial.getArray()[bc_owned-rs] = U.getArray(readonly=True)[bc_owned-rs]
                            try:
                                At, Rt = assemble(trial)
                            except IncrementTooLarge:
                                alpha *= 0.5
                                continue
                            try:
                                trial_norm, _ = constrain(Rt, bc_owned)
                                if trial_norm < norm or trial_norm <= tol:
                                    trial.copy(U)
                                    break
                            finally:
                                At.destroy()
                                Rt.destroy()
                            alpha *= 0.5
                        else:
                            raise RuntimeError("MPI line search failed to reduce the free residual")
                    else:
                        delta.getArray()[bc_owned-rs] = 0.0
                        U.axpy(1.0, delta)
                finally:
                    for obj in reversed(work):
                        obj.destroy()
            if not converged:
                raise RuntimeError("MPI Newton solve did not converge")
            if commit_element_fn is not None:
                values, old = ghost(U, Uloc), ghost(previous, Ploc)
                owner = getattr(commit_element_fn, '__self__', None)
                old_state = getattr(owner, 'svars', None)
                old_state = old_state.copy() if old_state is not None else None
                error = None
                try:
                    for ei, indices in enumerate(local_gm):
                        commit_element_fn(ei, my_coords[ei], values[indices], values[indices] - old[indices])
                except Exception as exc:
                    error = exc
                try:
                    collective_status(error, getattr(error, 'pnewdt', 1.0), comm)
                except Exception:
                    if old_state is not None:
                        owner.svars = old_state
                    raise
            if _commit_full is not None:
                _commit_full(full(U), full(previous))
            if verbose and rank == 0:
                print("step {}/{} converged".format(step, n_steps), flush=True)
        n_ghost = int(np.count_nonzero((needed < rs) | (needed >= re)))
        return full(U), dict(rank=rank, size=size, t_asm=t_asm, my_ne=len(my_gm),
                            n_owned=re-rs, n_ghost=n_ghost,
                            n_ghost_total=int(gather([n_ghost], comm).sum()),
                            ksp_its=ksp_its, ksp_diverged=False)
    except Exception:
        _reject([element_fn, batch_fn, commit_element_fn])
        raise
    finally:
        for obj in reversed(objects):
            obj.destroy()
