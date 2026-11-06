"""Small closed-form elements shared by the stepping and bounded-solver tests.

They constrain the solvers, not any physics: a Quad4 Poisson element, its L2
error, a strip mesh with prescribed ends, and a fake compiled module that lets
the public ``CompiledAbaqusElement`` run without a Fortran compiler while it
records every ``TIME``, ``DTIME`` and ``KINC`` it is called with.
"""

import numpy as np

from abaqus_ufl.fe import CompiledAbaqusElement, structured_quad_mesh

_GP = np.array([-1.0, 1.0]) / np.sqrt(3.0)
_SIGNS = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])


def _shape(xi, eta):
    N = 0.25 * (1 + _SIGNS[:, 0] * xi) * (1 + _SIGNS[:, 1] * eta)
    dN = 0.25 * np.column_stack([_SIGNS[:, 0] * (1 + _SIGNS[:, 1] * eta),
                                 _SIGNS[:, 1] * (1 + _SIGNS[:, 0] * xi)])
    return N, dN


def poisson_element(source_fn):
    """Quad4 element for ``-div grad u = f``: ``R_e = K_e U_e - f_e``, ``K_e = dR_e/dU_e``."""
    def element_fn(ei, coords, U_e, DU_e=None):
        coords = np.asarray(coords, dtype=float)
        K = np.zeros((4, 4))
        f = np.zeros(4)
        for xi in _GP:
            for eta in _GP:
                N, dN = _shape(xi, eta)
                J = dN.T @ coords
                detJ = np.linalg.det(J)
                if detJ <= 0.0:
                    raise ValueError("degenerate element in poisson_element")
                B = dN @ np.linalg.inv(J).T
                K += (B @ B.T) * detJ
                f += N * source_fn(N @ coords) * detJ
        return K @ np.asarray(U_e, dtype=float) - f, K
    return element_fn


def l2_error(nodes, elems, U, exact_fn):
    """L2 norm of a scalar nodal field minus ``exact_fn`` (2x2 Gauss)."""
    err2 = 0.0
    for el in elems:
        X = nodes[el]
        for xi in _GP:
            for eta in _GP:
                N, dN = _shape(xi, eta)
                diff = float(N @ U[el]) - exact_fn(N @ X)
                err2 += diff * diff * np.linalg.det(dN.T @ X)
    return float(np.sqrt(err2))


def convergence_rate(hs, errors):
    """Least-squares slope of log(error) against log(h)."""
    return float(np.polyfit(np.log(hs), np.log(errors), 1)[0])


def strip(nx, ny=2):
    """Unit-square mesh with u = 0 prescribed on x = 0 and x = 1 (one DOF per node)."""
    nodes, elems = structured_quad_mesh(nx, ny)
    ends = np.nonzero((nodes[:, 0] <= 1e-12) | (nodes[:, 0] >= 1.0 - 1e-12))[0]
    return nodes, elems, {int(n): 0.0 for n in ends}


class RecordingKernel:
    """Fake f2py module for a one-DOF-per-node Quad4 with one SVARS slot.

    The residual is ``K u - f + (mass / DTIME) * DU`` on each element (a
    Quad4 Poisson element plus a backward-Euler mass term); ``SVARS`` accumulate ``DTIME`` at every call (so a committed
    value counts accepted time).  ``calls`` records ``(TIME(1), TIME(2),
    DTIME, KSTEP, KINC)`` of every scalar or batch call.  ``reject`` is an
    optional ``(time, dtime, U, DU) -> bool`` returning PNEWDT = 0.5;
    ``nonfinite`` likewise returns a NaN residual.
    """

    def __init__(self, source=8.0, mass=0.0):
        self.source = source
        self.mass = mass
        self.calls = []
        self.reject = None
        self.nonfinite = None
        self._element = poisson_element(lambda xy: source)

    def _one(self, sv, coords, u, du, time, dtime, kstep, kinc):
        r, k = self._element(0, coords.T, u)
        if self.mass:
            r = r + self.mass / dtime * du
            k = k + self.mass / dtime * np.eye(4)
        pnewdt = 1.0
        if self.reject is not None and self.reject(time, dtime, u, du):
            pnewdt = 0.5
        if self.nonfinite is not None and self.nonfinite(time, dtime, u, du):
            r = r * np.nan
        sv = sv + dtime
        return -r, np.asfortranarray(k), sv, pnewdt

    def drive_uel(self, sv, coords, u, du, props, jprops, time, dtime, pnewdt,
                  lflags, params, jtype, kstep, kinc, jelem, period):
        self.calls.append((float(time[0]), float(time[1]), float(dtime), kstep, kinc))
        return self._one(np.array(sv, dtype=float), coords, u, du, time, dtime,
                         kstep, kinc)

    def drive_uel_batch(self, sv, coords, U, DU, props, jprops, time, dtime,
                        pnewdt, lflags, params, jtype, kstep, kinc, period, ids):
        self.calls.append((float(time[0]), float(time[1]), float(dtime), kstep, kinc))
        ne = U.shape[1]
        rhs, amatrx, state = np.zeros((4, ne)), np.zeros((4, 4, ne)), np.zeros_like(sv)
        pmin = 1.0
        for e in range(ne):
            r, k, s, p = self._one(sv[:, e], coords[:, :, e], U[:, e], DU[:, e],
                                   time, dtime, kstep, kinc)
            rhs[:, e], amatrx[:, :, e], state[:, e] = r, k, s
            pmin = min(pmin, p)
        return rhs, amatrx, state, pmin


def recording_element(n_elem, dt=1.0, time=None, source=8.0, mass=0.0,
                      cache_batches=False):
    """A public ``CompiledAbaqusElement`` on a :class:`RecordingKernel`.

    It bypasses the WeakForm signature check (``__new__`` plus the public
    attributes), as a host-local element may do.
    """
    kernel = RecordingKernel(source=source, mass=mass)
    element = CompiledAbaqusElement.__new__(CompiledAbaqusElement)
    element.module = kernel
    element.dof_per_node = 1
    element.props = np.zeros(1)
    element.ndofel = 4
    element.svars_size = 1
    element.svars = np.zeros((n_elem, 1))
    element.elem_ids = None
    element._jprops = np.array([0], dtype=np.int32)
    element._lflags = np.array([1, 0, 1, 0, 0, 0], dtype=np.int32)
    element._params = np.zeros(3)
    element.cache_batches = cache_batches
    element._batch_cache = None
    element.batch_evaluations = 0
    element.batch_cache_hits = 0
    element.dt = float(dt)
    element.begin_increment(1, time_begin=0.0 if time is None else time)
    return element, kernel
