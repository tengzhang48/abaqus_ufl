"""Rectangular meshes with the generated Quad4's counterclockwise ordering."""

import numpy as np


def structured_quad_mesh(nx, ny, Lx=1.0, Ly=1.0):
    """Return nodes (N,2) and Quad4 connectivity (M,4) on a rectangle."""
    for value in (nx, ny):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError("nx and ny must be positive integers")
    if not np.isfinite([Lx, Ly]).all() or Lx <= 0 or Ly <= 0:
        raise ValueError("rectangle lengths must be positive and finite")
    nodes = np.array([(x, y) for y in np.linspace(0, Ly, ny + 1)
                      for x in np.linspace(0, Lx, nx + 1)])
    elems = []
    for j in range(ny):
        for i in range(nx):
            first = j * (nx + 1) + i
            elems.append([first, first + 1, first + nx + 2, first + nx + 1])
    return nodes, np.asarray(elems, dtype=int)


def boundary_nodes(nodes, tol=1e-9):
    """Indices on a rectangular mesh's bounding box, not an arbitrary boundary."""
    nodes = np.asarray(nodes, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or not len(nodes) or not np.isfinite(nodes).all():
        raise ValueError("nodes must be a finite, nonempty (N,2) array")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be positive and finite")
    lower, upper = nodes.min(axis=0), nodes.max(axis=0)
    return np.flatnonzero(np.any((abs(nodes - lower) <= tol)
                                | (abs(nodes - upper) <= tol), axis=1))
