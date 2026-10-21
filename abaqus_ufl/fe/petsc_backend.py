"""Optional PETSc linear solvers for the serial host's CSR solve seam.

Distributed element assembly and solves live in :mod:`petsc_mpi`. No PETSc
initialization or import occurs until a solver is requested.
"""

import numpy as np
from scipy.sparse import csr_matrix


def _petsc():
    try:
        from petsc4py import PETSc
    except ImportError as error:
        raise ImportError("this solver requires the optional petsc4py dependency") from error
    if np.issubdtype(PETSc.ScalarType, np.complexfloating):
        raise ValueError("the FE host requires a real-scalar PETSc build")
    return PETSc


def _field_components(fields, dof_per_node):
    if (isinstance(dof_per_node, bool) or not isinstance(dof_per_node, (int, np.integer))
            or dof_per_node < 1):
        raise ValueError("dof_per_node must be a positive integer")
    fields = [(str(name), tuple(comps)) for name, comps in fields]
    components = [c for _, comps in fields for c in comps]
    if (not fields or any(not comps for _, comps in fields)
            or len(set(name for name, _ in fields)) != len(fields)
            or any(isinstance(c, bool) or not isinstance(c, (int, np.integer)) for c in components)
            or sorted(components) != list(range(dof_per_node))):
        raise ValueError("fields must uniquely partition all per-node DOF components")
    return fields


def make_petsc_solver(pc="lu", solver="mumps", ksp=None, rtol=1e-12,
                      atol=1e-50, maxit=2000, comm=None):
    """Return ``solve(K_csr, b) -> x`` on a single-rank communicator.

    Defaults to COMM_SELF, including when called by rank zero of an MPI job.
    Use ``solve_steps_mpi`` for distributed assembly/linear algebra. Iterative
    PCs (``gamg``, ``hypre``, etc.) default to nonsymmetric-safe GMRES.
    """
    PETSc = _petsc()
    comm = PETSc.COMM_SELF if comm is None else comm
    if comm.getSize() != 1:
        raise ValueError("CSR solve requires a single-rank communicator; use solve_steps_mpi")
    if not np.isfinite([rtol, atol]).all() or rtol <= 0 or atol < 0 or maxit < 1:
        raise ValueError("invalid PETSc solver tolerances or iteration limit")
    ksp_type = ksp if ksp is not None else ("preonly" if pc == "lu" else "gmres")

    def setup(ksp_obj, matrix):
        ksp_obj.setType(ksp_type)
        pc_obj = ksp_obj.getPC()
        pc_obj.setType(pc)
        if pc == "lu" and solver:
            pc_obj.setFactorSolverType(solver)
        return []

    return _csr_solver(PETSc, comm, setup, rtol, atol, maxit)


def _csr_solver(PETSc, comm, setup, rtol, atol, maxit):
    def linear_solve(K, b):
        A, b = csr_matrix(K), np.asarray(b, dtype=float)
        if A.shape[0] != A.shape[1] or b.shape != (A.shape[0],):
            raise ValueError("PETSc requires a square matrix and matching RHS vector")
        if not np.isfinite(A.data).all() or not np.isfinite(b).all():
            raise ValueError("PETSc matrix and RHS must be finite")
        if not len(b):
            return b.copy()
        objects = []
        try:
            mat = PETSc.Mat().createAIJ(
                size=A.shape, csr=(A.indptr.astype(PETSc.IntType),
                                  A.indices.astype(PETSc.IntType),
                                  A.data.astype(PETSc.ScalarType)), comm=comm)
            objects.append(mat)
            mat.assemble()
            rhs = PETSc.Vec().createWithArray(np.ascontiguousarray(b, PETSc.ScalarType), comm=comm)
            objects.append(rhs)
            x = rhs.duplicate()
            objects.append(x)
            ksp_obj = PETSc.KSP().create(comm)
            objects.append(ksp_obj)
            ksp_obj.setOperators(mat)
            objects.extend(setup(ksp_obj, mat))
            ksp_obj.setTolerances(rtol=rtol, atol=atol, max_it=maxit)
            ksp_obj.solve(rhs, x)
            reason = ksp_obj.getConvergedReason()
            out = x.getArray().copy()
            if reason <= 0 or not np.isfinite(out).all():
                raise RuntimeError("PETSc KSP failed (reason {})".format(reason))
            return out
        finally:
            for obj in reversed(objects):
                obj.destroy()
    return linear_solve


def make_fieldsplit_solver(field_components, dof_per_node, ksp_type="gmres",
                           split_type="additive", sub_pc="gamg",
                           rtol=1e-10, atol=1e-50, maxit=2000, comm=None):
    """Block preconditioner for interleaved, coupled nodal fields on COMM_SELF.

    ``field_components`` is ``[(name, [components...]), ...]``. Each component
    must appear exactly once. Schur splitting requires exactly two fields.
    The returned callable solves full interleaved systems. Its
    ``solve_reduced(K, b, free_dofs)`` method preserves original global DOF
    identities after constraint elimination; ``newton_solve`` uses it
    automatically. A fully constrained field is omitted from the solve.
    """
    PETSc = _petsc()
    comm = PETSc.COMM_SELF if comm is None else comm
    if comm.getSize() != 1:
        raise ValueError("CSR FieldSplit requires a single-rank communicator")
    fields = _field_components(field_components, dof_per_node)
    types = {"additive": PETSc.PC.CompositeType.ADDITIVE,
             "multiplicative": PETSc.PC.CompositeType.MULTIPLICATIVE,
             "schur": PETSc.PC.CompositeType.SCHUR}
    if split_type not in types or (split_type == "schur" and len(fields) != 2):
        raise ValueError("invalid FieldSplit type or Schur field count")

    def solve_reduced(K, b, free_dofs):
        def setup(ksp_obj, matrix):
            size = matrix.getSize()[0]
            if free_dofs is None:
                if size % dof_per_node:
                    raise ValueError("matrix size must be divisible by dof_per_node")
                components = np.arange(size) % dof_per_node
            else:
                dofs = np.asarray(free_dofs)
                if (dofs.shape != (size,) or not np.issubdtype(dofs.dtype, np.integer)
                        or np.any(dofs < 0) or len(np.unique(dofs)) != size):
                    raise ValueError("free_dofs must contain one distinct nonnegative integer per row")
                components = dofs % dof_per_node
            active = [(name, np.flatnonzero(np.isin(components, comps)).astype(PETSc.IntType))
                      for name, comps in fields]
            active = [(name, indices) for name, indices in active if len(indices)]
            ksp_obj.setType(ksp_type)
            pc_obj = ksp_obj.getPC()
            if len(active) == 1:
                # No Schur/block operation remains when every other field is fixed.
                pc_obj.setType(sub_pc)
                return []
            pc_obj.setType("fieldsplit")
            pc_obj.setFieldSplitType(types[split_type])
            index_sets = [PETSc.IS().createGeneral(indices, comm=comm) for _, indices in active]
            try:
                pc_obj.setFieldSplitIS(*[(name, indices)
                                        for (name, _), indices in zip(active, index_sets)])
                pc_obj.setUp()
                for sub in pc_obj.getFieldSplitSubKSP():
                    sub.setType("preonly")
                    sub.getPC().setType(sub_pc)
            except Exception:
                for indices in index_sets:
                    indices.destroy()
                raise
            return index_sets

        return _csr_solver(PETSc, comm, setup, rtol, atol, maxit)(K, b)

    def linear_solve(K, b):
        return solve_reduced(K, b, None)

    linear_solve.solve_reduced = solve_reduced
    return linear_solve
