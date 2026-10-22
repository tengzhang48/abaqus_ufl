"""Run under mpiexec: failure paths must terminate collectively and safely."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from petsc4py import PETSc
from scipy.sparse import eye as sparse_eye

from abaqus_ufl.fe import (IncrementTooLarge, make_petsc_solver,
                           partition_structured_quad, solve_steps_mpi_local,
                           structured_quad_mesh)
from abaqus_ufl.fe.mpi_batch import MPIBatchedElement, build_module, gather
from test_fe_batch_cache import AccumulatorProblem, Kernel

comm = PETSc.COMM_WORLD
rank, size = comm.getRank(), comm.getSize()
assert size == 2


def collective_failure(fn, exception=RuntimeError, match=None):
    caught = False
    try:
        fn()
    except exception as error:
        caught = match is None or match in str(error)
    assert np.all(gather([int(caught)], comm) == 1), 'not every rank saw the expected failure'


def check_local_host():
    ndof, gm, coords = partition_structured_quad(1, 1, 1, rank, size)
    assert len(gm) == rank  # rank zero owns no elements
    K = np.asfortranarray([[4., -1., 0., -0.5], [-0.2, 3., -1., 0.],
                          [0., -0.5, 3., -1.], [-1., 0., -0.3, 4.]])

    def element(ei, xy, u, du):
        exact = 1 + 2*xy[:, 0] + 3*xy[:, 1]
        return K @ (u-exact), K

    def batch(xy, u, du):
        exact = 1 + 2*xy[:, :, 0] + 3*xy[:, :, 1]
        return (u-exact) @ K.T, np.broadcast_to(K, (len(u), 4, 4))

    def run(fn=element, **kw):
        return solve_steps_mpi_local(ndof, gm, coords, 1, fn, lambda frac: {}, 1,
                                      pc='lu', solver='mumps', tol=1e-12, **kw)

    expected = np.array([1., 3., 4., 6.])
    for batched in (None, batch):
        U, info = run(batch_fn=batched)
        assert np.allclose(U, expected, rtol=0, atol=1e-12)
        assert info['my_ne'] == rank

    def broken(*args):
        raise RuntimeError('deliberate local error')

    collective_failure(lambda: run(broken), match='MPI operation failed')
    collective_failure(lambda: run(lambda *args: (np.full(4, np.nan), K)), match='MPI operation failed')

    def cutback(*args):
        error = IncrementTooLarge('deliberate cutback')
        error.pnewdt = 0.25
        raise error

    collective_failure(lambda: run(cutback), IncrementTooLarge, 'cutback')
    commits = []
    never_reduces = lambda *args: (np.ones(4), np.eye(4))
    collective_failure(lambda: run(never_reduces, commit_element_fn=lambda *args: commits.append(1)),
                       match='line search')
    collective_failure(lambda: run(never_reduces, max_newton=1, line_search=False,
                                   commit_element_fn=lambda *args: commits.append(1)),
                       match='did not converge')
    assert not commits

    class State:
        def __init__(self):
            self.svars = np.zeros((len(gm), 1))
            self.cache = None

        def commit(self, ei, xy, u, du):
            self.svars[ei] += 1
            raise IncrementTooLarge('commit refused')

        def reject_step(self):
            self.cache = None

    state = State()
    collective_failure(lambda: run(commit_element_fn=state.commit), IncrementTooLarge)
    assert np.count_nonzero(state.svars) == 0

    bad_coords = coords.copy()
    if rank == 1:
        bad_coords[0, 0, 0] = np.nan
    collective_failure(lambda: solve_steps_mpi_local(
        ndof, gm, bad_coords, 1, element, lambda frac: {}, 1), match='MPI operation failed')
    collective_failure(lambda: make_petsc_solver(comm=comm), ValueError, 'single-rank')


def check_batch_cache():
    problem = AccumulatorProblem()
    kernel = Kernel(problem)
    adapter = MPIBatchedElement(kernel, problem, n_elem=1, cache_batches=True)
    coords = np.array([[[0., 0.], [1., 0.], [1., 1.], [0., 1.]]])
    U, DU = np.zeros((1, 8)), np.zeros((1, 8))
    first = adapter.element_rk_batch(coords, U, DU)
    again = adapter.element_rk_batch(coords, U, DU)
    assert all(np.array_equal(a, b) for a, b in zip(first, again))
    assert adapter.batch_calls == 1 and adapter.batch_cache_hits == 1
    assert kernel.calls == rank
    if rank == 0:
        adapter.clear_batch_cache()
    adapter.element_rk_batch(coords, U, DU)
    assert adapter.batch_calls == 2
    adapter.cache_batches = rank != 0
    adapter.element_rk_batch(coords, U, DU)
    adapter.element_rk_batch(coords, U, DU)
    assert adapter.batch_calls == 4
    adapter.cache_batches = True
    adapter.element_rk_batch(coords, U, DU)
    adapter.element_rk_batch(coords, U, DU)
    assert adapter.batch_calls == 5
    assert adapter.batch_cache_hits == 2

    # Input conversion/preparation must vote before cache consensus, even on
    # the rank with no owned elements, with caching enabled or disabled.
    for cache in (False, True):
        adapter.cache_batches = cache
        for commit in (False, True):
            for argument in range(3):
                adapter.element_rk_batch(coords, U, DU)
                previous = adapter.svars.copy()
                invalid = [coords, U, DU]
                if rank == argument % size:
                    invalid[argument] = invalid[argument].astype(object)
                    invalid[argument].flat[0] = 'invalid-number'
                call = adapter.commit_group if commit else adapter.element_rk_batch
                collective_failure(lambda: call(*invalid), match='MPI operation failed')
                assert adapter._batch_cache is None
                assert np.array_equal(adapter.svars, previous)
            invalid_u = np.array(0.) if rank == 0 else U
            collective_failure(lambda: adapter.element_rk_batch(coords, invalid_u, DU),
                               match='MPI operation failed')
            invalid_u = U.copy()
            if rank == 1:
                invalid_u[0, 0] = np.nan
            collective_failure(lambda: adapter.element_rk_batch(coords, invalid_u, DU),
                               match='MPI operation failed')
        nodes, elems = structured_quad_mesh(1, 1)
        values = np.zeros(8)
        for failure in ('coords', 'U', 'layout'):
            adapter.element_rk_batch(coords, U, DU)
            previous = adapter.svars.copy()
            commit_nodes, commit_u, dpn = nodes, values, 2
            if rank == 0:
                if failure == 'coords':
                    commit_nodes = nodes.astype(object)
                    commit_nodes[0, 0] = 'invalid-number'
                elif failure == 'U':
                    commit_u = values.astype(object)
                    commit_u[0] = 'invalid-number'
                else:
                    dpn = 3
            collective_failure(lambda: adapter.commit_step(commit_nodes, elems, dpn, commit_u, values),
                               match='MPI operation failed')
            assert adapter._batch_cache is None
            assert np.array_equal(adapter.svars, previous)
    adapter.cache_batches = True
    adapter.element_rk_batch(coords, U, DU)  # recover after the rejected inputs

    U += 1
    if rank == 1:
        kernel.failure = 'cutback'
    collective_failure(lambda: adapter.commit_group(coords, U, DU), IncrementTooLarge)
    assert adapter._batch_cache is None
    assert np.count_nonzero(adapter.svars) == 0
    kernel.failure = None
    adapter.elem_ids = np.array([17]) if rank else np.array([1, 2])
    collective_failure(lambda: adapter.element_rk_batch(coords, U, DU), match='MPI operation failed')
    adapter.elem_ids = None
    root_solve = adapter.root_linear_solver(pc='invalid_pc')
    collective_failure(lambda: root_solve(sparse_eye(2), np.ones(2)), match='MPI operation failed')
    collective_failure(lambda: build_module(object(), '.unused-mpi-build'), match='MPI operation failed')


check_local_host()
check_batch_cache()
if rank == 0:
    print('MPI contracts passed', flush=True)
