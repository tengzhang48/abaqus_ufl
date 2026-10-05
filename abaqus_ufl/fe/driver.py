"""Serial sparse assembly and Newton solve for small verification meshes."""

import warnings

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve


def _mesh_inputs(nodes, elems, dof_per_node):
    nodes = np.asarray(nodes, dtype=float)
    elems = np.asarray(elems)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or not len(nodes) or not np.isfinite(nodes).all():
        raise ValueError("nodes must be a finite, nonempty (N,2) array")
    if (elems.ndim != 2 or elems.shape[1] != 4 or not len(elems)
            or not np.issubdtype(elems.dtype, np.integer)
            or elems.min() < 0 or elems.max() >= len(nodes)
            or any(len(set(row)) != 4 for row in elems)):
        raise ValueError("elems must contain four distinct, valid node indices per Quad4")
    if isinstance(dof_per_node, bool) or not isinstance(dof_per_node, (int, np.integer)) or dof_per_node < 1:
        raise ValueError("dof_per_node must be a positive integer")
    return nodes, elems.astype(int, copy=False), len(nodes) * dof_per_node


def _vector(value, size, name):
    value = np.asarray(value, dtype=float)
    if value.shape != (size,) or not np.isfinite(value).all():
        raise ValueError(name + " must be a finite vector of the global DOF size")
    return value


def assemble(nodes, elems, dof_per_node, element_fn, U, DU=None, *, du_is_total=False):
    """Assemble (R,K) with K=dR/dU in interleaved node/component DOF order.

    ``element_fn(index, coordinates, U_element, DU_element)`` returns the
    residual and tangent. ``DU`` is measured from the previous accepted nodal
    solution. Pass ``du_is_total=True`` for an explicitly from-zero evaluation.
    Essential-boundary reactions remain in the returned residual.
    """
    nodes, elems, size = _mesh_inputs(nodes, elems, dof_per_node)
    U = _vector(U, size, "U")
    if DU is None and not du_is_total:
        raise ValueError("pass DU=U-U_previous, or du_is_total=True for a from-zero evaluation")
    if DU is not None and du_is_total:
        raise ValueError("specify DU or du_is_total, not both")
    DU = U if DU is None else _vector(DU, size, "DU")
    residual = np.zeros(size)
    rows, cols, values = [], [], []
    for index, connectivity in enumerate(elems):
        dofs = (connectivity[:, None] * dof_per_node
                + np.arange(dof_per_node)).ravel()
        local_r, local_k = element_fn(index, nodes[connectivity], U[dofs], DU[dofs])
        local_r, local_k = np.asarray(local_r), np.asarray(local_k)
        if local_r.shape != (len(dofs),) or local_k.shape != (len(dofs), len(dofs)):
            raise ValueError("element output shapes disagree with the DOF layout")
        if not np.isfinite(local_r).all() or not np.isfinite(local_k).all():
            raise RuntimeError("nonfinite element residual or tangent")
        np.add.at(residual, dofs, local_r)
        rows.extend(np.repeat(dofs, len(dofs)))
        cols.extend(np.tile(dofs, len(dofs)))
        values.extend(local_k.ravel())
    tangent = coo_matrix((values, (rows, cols)), shape=(size, size)).tocsr()
    return residual, tangent


def newton_solve(nodes, elems, dof_per_node, element_fn, dirichlet,
                 U0=None, U_previous=None, max_iter=25, tol=1e-9):
    """Solve free-DOF equilibrium with prescribed global DOF values.

    Returns U only on convergence. Previous nodal history is fixed throughout
    all Newton and line-search trials. Singular/nonfinite solves and failed
    convergence raise an error; this small driver has no automatic time
    cutback, distributed assembly, contact, or material-state transaction.
    ``tol`` is relative to the larger of 1 and the largest residual entry, reactions included.
    """
    nodes, elems, size = _mesh_inputs(nodes, elems, dof_per_node)
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be positive and finite")
    if isinstance(max_iter, bool) or not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
    U = np.zeros(size) if U0 is None else _vector(U0, size, "U0").copy()
    previous = np.zeros(size) if U_previous is None else _vector(U_previous, size, "U_previous").copy()
    free_mask = np.ones(size, dtype=bool)
    for dof, value in dirichlet.items():
        if (isinstance(dof, bool) or not isinstance(dof, (int, np.integer))
                or dof < 0 or dof >= size or not np.isfinite(value)):
            raise ValueError("Dirichlet conditions require valid DOF indices and finite values")
        U[dof] = value
        free_mask[dof] = False
    free = np.flatnonzero(free_mask)
    for iteration in range(max_iter + 1):
        residual, tangent = assemble(nodes, elems, dof_per_node, element_fn, U, U - previous)
        norm = np.linalg.norm(residual[free], ord=np.inf) if len(free) else 0.0
        # Reactions carry the force scale, so tol is relative to them (or 1).
        if norm <= tol * max(1.0, np.linalg.norm(residual, ord=np.inf)):
            return U
        if iteration == max_iter:
            break
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", MatrixRankWarning)
                delta = spsolve(tangent[free][:, free], -residual[free])
        except (MatrixRankWarning, RuntimeError) as error:
            raise RuntimeError("serial FE linear solve failed") from error
        if not np.isfinite(delta).all():
            raise RuntimeError("serial FE linear solve returned a nonfinite increment")
        for search in range(12):
            trial = U.copy()
            trial[free] += (0.5 ** search) * delta
            try:
                trial_r, _ = assemble(nodes, elems, dof_per_node, element_fn,
                                      trial, trial - previous)
            except RuntimeError:
                continue
            trial_norm = np.linalg.norm(trial_r[free], ord=np.inf)
            if trial_norm < norm or trial_norm <= tol:
                U = trial
                break
        else:
            raise RuntimeError("serial FE line search failed to reduce the free residual")
    raise RuntimeError("serial FE Newton solve did not converge in {} iterations".format(max_iter))
