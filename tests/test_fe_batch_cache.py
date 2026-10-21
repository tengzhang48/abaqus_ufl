"""Exact batch reuse must preserve residuals, tangents and state transactions."""

from types import SimpleNamespace

import numpy as np
import pytest

import abaqus_ufl as au
from abaqus_ufl.core.tensor import eye, trace
from abaqus_ufl.fe import (CompiledAbaqusElement, IncrementTooLarge, assemble,
                           assemble_batched, build_compiled_uel, newton_solve,
                           structured_quad_mesh)
from abaqus_ufl.fe.compiled_element import _signature


class Accumulator(au.Material):
    props = dict(G=2.0)
    state_vars = dict(a=0.0)

    def stress_PK1(self, F, a_old, dt):
        strain = F - eye(3)
        P = self.G * strain
        a_new = a_old + dt * trace(strain)
        return P, {'a': a_new}


class AccumulatorProblem(au.WeakForm):
    material = Accumulator
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField('u', degree=1)

    def momentum_equation(self, v, F, a_old, dt):
        return self.material.stress_PK1(F, a_old, dt)


class Kernel:
    """Emulate f2py's in-place trial-state update and a nonsymmetric tangent."""
    def __init__(self, problem):
        self._abaqus_ufl_signature = _signature(problem)
        self.calls = 0
        self.failure = None

    def drive_uel_batch(self, sv, coords, U, DU, props, jprops, time, dt,
                        pnewdt, flags, params, jtype, kstep, kinc, period, ids):
        self.calls += 1
        rhs = U + 2 * DU + props[0] + coords.sum(axis=(0, 1)) + sv.sum(axis=0)
        tangent = np.repeat(np.arange(64.).reshape(8, 8, 1), U.shape[1], axis=2)
        sv += 1 + U.sum(axis=0) + dt + ids
        if self.failure == 'exception':
            raise RuntimeError('kernel failed')
        if self.failure == 'nonfinite':
            tangent[0, 0, 0] = np.nan
        return rhs, tangent, sv, 0.5 if self.failure == 'cutback' else 1.0


def case(cache=True):
    problem = AccumulatorProblem()
    kernel = Kernel(problem)
    element = CompiledAbaqusElement(kernel, problem, n_elem=2, cache_batches=cache)
    nodes, elems = structured_quad_mesh(2, 1)
    coords = nodes[elems]
    U = np.arange(16.).reshape(2, 8) / 20
    DU = U / 2
    return element, kernel, (coords, U, DU)


def test_cache_opt_in_and_identical_results():
    cached, kc, inputs = case()
    plain, kp, _ = case(False)
    for _ in range(3):
        expected = plain.element_rk_batch(*inputs)
        actual = cached.element_rk_batch(*inputs)
        for a, b in zip(actual, expected):
            assert np.array_equal(a, b)
        assert np.array_equal(cached.svars, plain.svars)
    assert kc.calls == cached.batch_evaluations == 1
    assert cached.batch_cache_hits == 2
    assert kp.calls == plain.batch_evaluations == 3
    assert plain.batch_cache_hits == 0
    assert actual[1].flags.c_contiguous


def test_cache_outputs_and_inputs_cannot_alias_entry():
    element, kernel, inputs = case()
    expected = tuple(a.copy() for a in element.element_rk_batch(*inputs))
    r, k = element.element_rk_batch(*inputs)
    r[:] = -99
    k[:] = -99
    for a, b in zip(element.element_rk_batch(*inputs), expected):
        assert np.array_equal(a, b)
    inputs[1][0, 0] += 0.1
    element.element_rk_batch(*inputs)
    assert kernel.calls == 2
    assert np.array_equal(element.svars, np.zeros((2, 4)))


@pytest.mark.parametrize('changed', [
    'coords', 'U', 'DU', 'svars', 'props', 'dt', 'time', 'kinc',
    'elem_ids', '_jprops', '_lflags', '_params', 'module', 'signed_zero',
])
def test_every_variable_uel_input_invalidates_cache(changed):
    element, kernel, inputs = case()
    inputs[1][0, 0] = 0.0
    element.element_rk_batch(*inputs)
    if changed in ('coords', 'U', 'DU'):
        inputs[('coords', 'U', 'DU').index(changed)].flat[0] += 0.01
    elif changed == 'signed_zero':
        inputs[1][0, 0] = -0.0
    elif changed == 'module':
        element.module = Kernel(AccumulatorProblem())
    elif changed == 'elem_ids':
        element.elem_ids = np.array([7, 11], dtype=np.int32)
    elif changed in ('dt', 'time', 'kinc'):
        setattr(element, changed, getattr(element, changed) + 1)
    else:
        getattr(element, changed).flat[0] += 1
    element.element_rk_batch(*inputs)
    assert element.batch_evaluations == 2
    assert element.batch_cache_hits == 0


def test_commit_reuses_trial_then_clears_cache():
    cached, kernel, inputs = case()
    plain, _, _ = case(False)
    cached.element_rk_batch(*inputs)
    assert np.count_nonzero(cached.svars) == 0
    cached.commit_group(*inputs)
    plain.commit_group(*inputs)
    assert kernel.calls == 1
    assert cached.batch_cache_hits == 1
    assert cached._batch_cache is None
    assert np.array_equal(cached.svars, plain.svars)
    cached.element_rk_batch(*inputs)
    assert kernel.calls == 2  # the committed old state is now different


@pytest.mark.parametrize('failure', ['cutback', 'nonfinite', 'exception'])
def test_rejection_does_not_cache_or_commit_trial_state(failure):
    element, kernel, inputs = case()
    element.element_rk_batch(*inputs)
    inputs[1][0, 0] += 1
    kernel.failure = failure
    error = IncrementTooLarge if failure == 'cutback' else RuntimeError
    with pytest.raises(error):
        element.commit_group(*inputs)
    assert element._batch_cache is None
    assert np.array_equal(element.svars, np.zeros((2, 4)))
    kernel.failure = None
    element.element_rk_batch(*inputs)
    assert kernel.calls == 3


def test_manual_clear_rejection_and_increment_preserve_counters():
    element, kernel, inputs = case()
    for clear in (element.clear_batch_cache, element.reject_step,
                  lambda: element.begin_increment(2)):
        element.element_rk_batch(*inputs)
        clear()
        assert element._batch_cache is None
    assert kernel.calls == element.batch_evaluations == 3
    assert element.batch_cache_hits == 0


@pytest.mark.parametrize('failure', ['value', 'type', 'runtime', 'result', 'callback'])
def test_failed_newton_clears_trial_cache(failure):
    element, _, _ = case()
    nodes, elems = structured_quad_mesh(2, 1)
    previous = element.svars.copy()
    error_type = {'type': TypeError, 'runtime': RuntimeError}.get(failure, ValueError)

    def fail(*args):
        if failure == 'result':
            return ['invalid-number'] * len(args[1])
        raise error_type('deliberate failure')

    options = {'record_iteration': fail} if failure == 'callback' else {'linear_solve': fail}
    with pytest.raises(error_type):
        newton_solve(nodes, elems, 2, None, {}, batch_fn=element.element_rk_batch, **options)
    assert element._batch_cache is None
    assert np.array_equal(element.svars, previous)


@pytest.mark.parametrize('argument', [0, 1, 2])
@pytest.mark.parametrize('commit', [False, True])
def test_batch_conversion_failure_clears_cache_without_committing(argument, commit):
    element, _, inputs = case()
    element.element_rk_batch(*inputs)
    previous = element.svars.copy()
    invalid = list(inputs)
    invalid[argument] = invalid[argument].astype(object)
    invalid[argument].flat[0] = 'invalid-number'
    call = element.commit_group if commit else element.element_rk_batch
    with pytest.raises(ValueError, match='could not convert'):
        call(*invalid)
    assert element._batch_cache is None
    assert np.array_equal(element.svars, previous)


def test_commit_step_preparation_failure_clears_cache_without_committing():
    element, _, inputs = case()
    element.element_rk_batch(*inputs)
    nodes, elems = structured_quad_mesh(2, 1)
    previous = element.svars.copy()
    U = np.zeros(len(nodes) * 2)
    with pytest.raises(ValueError, match='commit DOF layout'):
        element.commit_step(nodes, elems, 3, U, U)
    assert element._batch_cache is None
    assert np.array_equal(element.svars, previous)


def test_scalar_calls_stay_uncached_and_scalar_commit_invalidates_batch():
    element, kernel, inputs = case()
    calls = []

    def scalar(sv, coords, u, du, *args):
        calls.append(1)
        sv += 2
        return u.copy(), np.asfortranarray(np.arange(64.).reshape(8, 8)), sv, 1.0

    kernel.drive_uel = scalar
    element.element_rk_batch(*inputs)
    for _ in range(2):
        r, k = element.element_rk(0, inputs[0][0], inputs[1][0], inputs[2][0])
        assert k.flags.c_contiguous
    assert len(calls) == 2 and element.batch_cache_hits == 0
    assert np.count_nonzero(element.svars) == 0
    element.commit_element(0, inputs[0][0], inputs[1][0], inputs[2][0])
    assert np.array_equal(element.svars[0], np.full(4, 2.0))
    assert element._batch_cache is None


def test_batch_assembly_matches_scalar_for_nonsymmetric_operator():
    nodes, elems = structured_quad_mesh(3, 2)
    U = np.arange(len(nodes), dtype=float)
    K = np.arange(16.).reshape(4, 4)
    scalar = lambda ei, coords, u, du: (K @ u + du, K)
    batch = lambda coords, u, du: (u @ K.T + du, np.broadcast_to(K, (len(u), 4, 4)))
    r, k = assemble(nodes, elems, 1, scalar, U, U / 3)
    rb, kb = assemble_batched(nodes, elems, 1, batch, U, U / 3)
    assert np.array_equal(r, rb)
    assert np.array_equal(k.toarray(), kb.toarray())


def test_compiled_scalar_batch_cache_and_state_commit(tmp_path):
    problem = AccumulatorProblem()
    module = build_compiled_uel(problem, workdir=tmp_path)
    nodes, elems = structured_quad_mesh(2, 1)
    gm = (elems[:, :, None] * 2 + np.arange(2)).reshape(2, 8)
    U = (nodes @ np.array([[0.05, 0.03], [0.02, -0.01]])).ravel()
    element = CompiledAbaqusElement(module, problem, n_elem=2, cache_batches=True)
    rhs, tangent = element.element_rk_batch(nodes[elems], U[gm], U[gm])
    for ei in range(2):
        r, k = element.element_rk(ei, nodes[elems[ei]], U[gm[ei]], U[gm[ei]])
        assert np.array_equal(rhs[ei], r)
        assert np.array_equal(tangent[ei], k)
    assert element.batch_evaluations == 1
    assert np.count_nonzero(element.svars) == 0
    element.commit_step(nodes, elems, 2, U, np.zeros_like(U))
    assert element.batch_cache_hits == 1
    assert element._batch_cache is None
    assert np.allclose(element.svars, 0.04, atol=1e-15)
    element.begin_increment(2)
    element.commit_step(nodes, elems, 2, U, U)
    assert np.allclose(element.svars, 0.08, atol=1e-15)
