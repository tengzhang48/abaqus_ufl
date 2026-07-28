"""State-order regression for the reference assembler (release blocker B).

The generated UEL packs a 3x3 tensor state variable into SVARS
column-major: slot ``offset + 3*(j-1) + i`` holds component ``(i, j)``
(1-based; see ``uel_gen._append_svars_read``). The reference assembler
claims Fortran parity, so its flatten/reshape pair must use the SAME
layout. These tests pin that layout with a NONSYMMETRIC tensor (a
symmetric or diagonal test state cannot distinguish C from Fortran
order), mixed with scalar and vector state, across multiple Gauss
points.
"""

from collections import OrderedDict

import numpy as np

from abaqus_ufl.core.reference_assembly import (
    _flatten_gp_state,
    _state_args_for_gp,
    _state_var_sizes,
)


class _MixedStateStub:
    """Minimal object exposing the state_vars interface the helpers read."""

    state_vars = OrderedDict([
        ('a', 0.0),
        ('v', np.zeros(3)),
        ('T', np.eye(3)),
    ])


def _gp_state(gp):
    """Distinct, deliberately nonsymmetric per-GP state values."""
    base = 100.0 * (gp + 1)
    T = base + np.arange(9, dtype=float).reshape(3, 3)  # T[i,j]=base+3i+j
    assert np.max(np.abs(T - T.T)) > 1.0  # the discriminator must exist
    return {
        'a': base + 50.0,
        'v': base + np.array([1.0, 2.0, 3.0]),
        'T': T,
    }


def test_tensor_slots_match_generated_svars_layout():
    """flat[o + 3*j + i] == T[i, j] (0-based), per the uel_gen packing."""
    mat = _MixedStateStub()
    n_gp = 3
    state_by_gp = {gp: _gp_state(gp) for gp in range(n_gp)}
    flat = _flatten_gp_state(mat, state_by_gp, n_gp, STATEV=None)

    per_gp = sum(s for _, s, _ in _state_var_sizes(mat))
    assert per_gp == 1 + 3 + 9
    for gp in range(n_gp):
        expected = _gp_state(gp)
        off = gp * per_gp
        assert flat[off] == expected['a']
        assert np.array_equal(flat[off + 1:off + 4], expected['v'])
        T = expected['T']
        for i in range(3):
            for j in range(3):
                assert flat[off + 4 + 3 * j + i] == T[i, j], (
                    "tensor slot {} at GP {} is not column-major".format(
                        3 * j + i, gp))

    # Discriminator: a C-order flatten of the nonsymmetric tensor must NOT
    # reproduce the stored slots.
    wrong = _gp_state(0)['T'].ravel(order='C')
    assert not np.array_equal(flat[4:13], wrong)


def test_mixed_state_round_trip_across_gps():
    """flatten -> read-back returns the exact per-GP values for every GP."""
    mat = _MixedStateStub()
    n_gp = 4
    state_by_gp = {gp: _gp_state(gp) for gp in range(n_gp)}
    flat = _flatten_gp_state(mat, state_by_gp, n_gp, STATEV=None)

    for gp in range(n_gp):
        back = _state_args_for_gp(mat, flat, gp)
        expected = _gp_state(gp)
        assert back['a_old'] == expected['a']
        assert np.array_equal(back['v_old'], expected['v'])
        assert np.array_equal(back['T_old'], expected['T']), (
            "GP {} tensor state did not round-trip (transposed read?)".format(
                gp))


def test_gp_without_update_falls_back_to_incoming_statev():
    mat = _MixedStateStub()
    n_gp = 2
    full = _flatten_gp_state(
        mat, {gp: _gp_state(gp) for gp in range(n_gp)}, n_gp, STATEV=None)

    # Only GP 1 updates; GP 0 must keep the incoming STATEV block.
    partial = _flatten_gp_state(mat, {1: _gp_state(1)}, n_gp, STATEV=full)
    per_gp = sum(s for _, s, _ in _state_var_sizes(mat))
    assert np.array_equal(partial[:per_gp], full[:per_gp])
    back = _state_args_for_gp(mat, partial, 0)
    assert np.array_equal(back['T_old'], _gp_state(0)['T'])
