"""Gates for the bounded/irreversible (KKT) solver and its transactional state.

The oracle is the classical **obstacle problem**, which has a closed-form
solution that is *not* the clipped unconstrained solution, so post-convergence
clipping passes a "field stays in [lo, hi]" check while failing these gates:

    minimise  1/2 |grad u|^2 - f u    over  u <= c,   u(0) = u(1) = 0

On a strip with natural top/bottom boundaries the solution is one-dimensional.
With a constant source ``f`` and a cap ``c <= f/8`` the contact set is
``[s, 1-s]`` with the smooth-fit free boundary

    s = sqrt(2c/f),   u(x) = f x (2s - x) / 2  for x <= s,   u = c on [s, 1-s]

(``u'(s) = 0`` is the free-boundary condition; the clipped parabola violates it).
"""

import numpy as np
import pytest
import scipy.sparse as sp

from abaqus_ufl.fe import (BoxSpec, CompiledElementState, IncrementTooLarge,
                           StaggeredSolveFailed, active_set,
                           activity_equivalent_vi_certificate, assemble,
                           boundary_nodes, projected_residual, solve_steps,
                           solve_steps_bounded, solve_steps_staggered)
from abaqus_ufl.fe.bounds import _bounded_newton
from abaqus_ufl.fe.stepping import _scipy_solve
from fe_stepping_toys import (convergence_rate, l2_error, poisson_element,
                              recording_element, strip)

F_SOURCE = 8.0        # constant body source
CAP = 0.5             # obstacle:  u <= CAP  (unconstrained peak is f/8 = 1.0)


def test_activity_equivalent_vi_certificate_preserves_sign_and_threshold():
    """A sub-activity-tolerance VI correction is certified, while a resolvable
    growth direction is not; classical projected-residual signs stay unchanged."""
    lo, hi, is_con, weights = np.array([0.0]), np.array([1.0]), np.array([False]), np.ones(1)

    def asm(U, _DU):
        return np.array([10.0 * U[0] - 0.05]), sp.csr_matrix([[10.0]])

    U = np.array([0.0])
    R, K = asm(U, U)
    classical, _, _ = projected_residual(R, U, lo, hi, at_tol=0.01)
    np.testing.assert_allclose(classical, [-0.05])
    assert activity_equivalent_vi_certificate(
        asm, U, U, is_con, lo, hi, weights, 0.01, 1e-12, _scipy_solve, R=R, K=K) is not None
    assert activity_equivalent_vi_certificate(
        asm, U, U, is_con, lo, hi, weights, 0.001, 1e-12, _scipy_solve, R=R, K=K) is None


def test_activity_vi_newton_hook_is_explicit_and_default_inert():
    # Exact root is 5e-13 above the floor: activity-equivalent only when the
    # opt-in certificate is selected.
    def asm(U, _DU):
        return np.array([U[0] - 5.0e-13]), sp.csr_matrix([[1.0]])

    common = dict(
        step_asm=asm, U=np.array([0.0]), U_prev=np.array([0.0]),
        con=np.array([], dtype=int), is_con=np.array([False]),
        lo=np.array([0.0]), hi=np.array([1.0]), w_res=np.ones(1),
        at_tol=1.0e-12, max_newton=1, tol=1.0e-14,
        linear_solve=_scipy_solve, line_search=False, context="unit")
    classical = _bounded_newton(**common)
    certified = _bounded_newton(**common, activity_vi_certificate=True)
    assert not classical["converged"]
    assert certified["converged"] and certified["activity_vi_certified"]
    assert certified["certified_update_norm"] == pytest.approx(5.0e-13)
    # Return the certified nearby VI point, not the floor representative.
    np.testing.assert_allclose(certified["U"], [5.0e-13], rtol=0.0, atol=0.0)


# --------------------------------------------------------------------------
# problem set-up (shared)
# --------------------------------------------------------------------------

def _exact_obstacle(x, f=F_SOURCE, c=CAP):
    s = np.sqrt(2.0 * c / f)
    x = np.minimum(x, 1.0 - x)            # symmetric about x = 1/2
    return np.where(x >= s, c, f * x * (2.0 * s - x) / 2.0)


def _load_vector(nodes, elems, f=F_SOURCE):
    """Consistent nodal load for a constant source (R = K U - F at U = 0)."""
    R, _ = assemble(nodes, elems, 1, poisson_element(lambda xy: f),
                    np.zeros(len(nodes)), du_is_total=True)
    return -R


def _pure_laplacian():
    return poisson_element(lambda xy: 0.0)


def _external(F_vec, profile=lambda frac: frac):
    return lambda step, frac, U, U_prev: (profile(frac) * F_vec, None)


def _solve_bounded(nx, n_steps=1, box=None, **kw):
    nodes, elems, dirichlet = strip(nx)
    box = box or BoxSpec(bounds={0: (-np.inf, CAP)})
    U, report = solve_steps_bounded(
        nodes, elems, 1, _pure_laplacian(), lambda frac: dirichlet, n_steps,
        box, external_fn=_external(_load_vector(nodes, elems)), tol=1e-10, **kw)
    return nodes, elems, U, report


def _kkt_norm_of(nodes, elems, U, dirichlet, lo, hi, frac=1.0):
    """Projected-residual norm of an arbitrary state (used on broken controls)."""
    R, _ = assemble(nodes, elems, 1, _pure_laplacian(), U, du_is_total=True)
    R = R - frac * _load_vector(nodes, elems)
    R[np.array(sorted(dirichlet), dtype=int)] = 0.0
    Rp, _, _ = projected_residual(R, U, lo, hi)
    return float(np.linalg.norm(Rp)), float(np.linalg.norm(R))


# --------------------------------------------------------------------------
# 1. the bounded solve reproduces the exact obstacle solution
# --------------------------------------------------------------------------

def test_bounded_solve_matches_exact_obstacle_solution():
    nodes, elems, U, report = _solve_bounded(48)
    # The nodal error is set by whether the last contact node falls inside the
    # exact contact interval: about one node's worth of the profile.
    assert np.max(np.abs(U - _exact_obstacle(nodes[:, 0]))) < 2.5e-3
    assert U.max() <= CAP + 1e-12, "primal feasibility violated"
    assert report.steps_accepted == 1 and report.cutbacks == 0


def test_bounded_solve_converges_under_refinement():
    """L2 convergence at about O(h^2), and a free boundary within one element."""
    exact_1d = lambda xy: float(_exact_obstacle(np.array([xy[0]]))[0])
    s = np.sqrt(2.0 * CAP / F_SOURCE)
    hs, errs = [], []
    for nx in (8, 16, 32, 64):
        nodes, elems, U, _ = _solve_bounded(nx)
        hs.append(1.0 / nx)
        errs.append(l2_error(nodes, elems, U, exact_1d))
        row = nodes[:, 1] == 0.0
        contact_x = nodes[row, 0][U[row] >= CAP - 1e-12]
        assert abs(contact_x.min() - s) <= 1.0 / nx + 1e-12
    assert all(b < a for a, b in zip(errs, errs[1:])), errs
    assert 1.6 < convergence_rate(hs, errs) < 2.6


def test_bounded_solution_is_a_kkt_point_not_merely_feasible():
    """The projected residual vanishes while the raw residual (the contact
    reaction) does not."""
    nodes, elems, U, _ = _solve_bounded(32)
    _, _, dirichlet = strip(32)
    lo, hi = np.full(len(nodes), -np.inf), np.full(len(nodes), CAP)
    kkt, raw = _kkt_norm_of(nodes, elems, U, dirichlet, lo, hi)
    assert kkt < 1e-8
    assert raw > 1e-2, "the constraint is inactive; this gate would not discriminate"


def test_active_set_matches_the_analytic_contact_interval():
    nx = 48
    nodes, elems, U, _ = _solve_bounded(nx)
    _, _, dirichlet = strip(nx)
    R, _ = assemble(nodes, elems, 1, _pure_laplacian(), U, du_is_total=True)
    R = R - _load_vector(nodes, elems)
    R[np.array(sorted(dirichlet), dtype=int)] = 0.0
    act = active_set(R, U, np.full(len(nodes), -np.inf), np.full(len(nodes), CAP),
                     at_tol=1e-9)
    s = np.sqrt(2.0 * CAP / F_SOURCE)
    x = nodes[:, 0]
    assert act[(x > s + 2.0 / nx) & (x < 1.0 - s - 2.0 / nx)].all()
    assert not act[(x < s - 2.0 / nx) | (x > 1.0 - s + 2.0 / nx)].any()


# --------------------------------------------------------------------------
# 2. broken controls: clipping after convergence is NOT a constrained solve
# --------------------------------------------------------------------------

def test_clipping_after_convergence_is_not_a_kkt_point():
    nx = 48
    nodes, elems, dirichlet = strip(nx)
    U_free = solve_steps(nodes, elems, 1, _pure_laplacian(), lambda frac: dirichlet, 1,
                         tol=1e-12, external_fn=_external(_load_vector(nodes, elems)))
    U_clipped = np.minimum(U_free, CAP)
    lo, hi = np.full(len(nodes), -np.inf), np.full(len(nodes), CAP)
    kkt_clipped, _ = _kkt_norm_of(nodes, elems, U_clipped, dirichlet, lo, hi)
    _, _, U_bounded, _ = _solve_bounded(nx)
    kkt_bounded, _ = _kkt_norm_of(nodes, elems, U_bounded, dirichlet, lo, hi)
    assert kkt_clipped > 1e3 * max(kkt_bounded, 1e-14)
    assert np.max(np.abs(U_clipped - _exact_obstacle(nodes[:, 0]))) > 0.05


def test_clipping_error_does_not_vanish_under_refinement():
    errs = []
    for nx in (12, 24, 48):
        nodes, elems, dirichlet = strip(nx)
        U_free = solve_steps(nodes, elems, 1, _pure_laplacian(), lambda frac: dirichlet, 1,
                             tol=1e-12, external_fn=_external(_load_vector(nodes, elems)))
        errs.append(np.max(np.abs(np.minimum(U_free, CAP) - _exact_obstacle(nodes[:, 0]))))
    assert min(errs) > 0.05, errs


def test_field_bounds_path_refuses_to_commit_a_clipped_state():
    """The staggered clip is followed by a fresh projected-residual
    certificate, so a clip that invalidates equilibrium is rejected."""
    nx = 24
    nodes, elems, dirichlet = strip(nx)
    source = poisson_element(lambda xy: F_SOURCE)
    U_free = solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet,
                                   [[0]], 1, tol=1e-10)
    U_clipped = np.minimum(U_free, CAP)
    R, _ = assemble(nodes, elems, 1, source, U_clipped, du_is_total=True)
    R[np.array(sorted(dirichlet), dtype=int)] = 0.0
    Rp, _, _ = projected_residual(R, U_clipped, np.full(len(nodes), -np.inf),
                                  np.full(len(nodes), CAP))
    assert np.linalg.norm(Rp) > 1e-2
    with pytest.raises(StaggeredSolveFailed,
                       match="post-clamp projected residual failed") as caught:
        solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1,
                              tol=1e-10, field_bounds={0: (-np.inf, CAP)})
    assert caught.value.reason == "post-clamp projected residual failed"


# --------------------------------------------------------------------------
# 3. irreversibility (dynamic step-start lower bound)
# --------------------------------------------------------------------------

def _load_unload_history(box, n_steps=8, nx=24):
    """Ramp the load up over the first half of the history, then back down."""
    nodes, elems, dirichlet = strip(nx)
    profile = lambda frac: (2.0 * frac if frac <= 0.5 else 2.0 * (1.0 - frac))
    history = []
    U, report = solve_steps_bounded(
        nodes, elems, 1, _pure_laplacian(), lambda frac: dirichlet, n_steps, box,
        external_fn=_external(_load_vector(nodes, elems), profile), tol=1e-10,
        step_callback=lambda s, f, Uc, rep: history.append(Uc))
    return nodes, np.array(history), U, report


def test_irreversible_field_never_heals_on_unloading():
    box = BoxSpec(bounds={0: (0.0, CAP)}, irreversible=(0,))
    _, hist, _, _ = _load_unload_history(box)
    assert np.min(np.diff(hist, axis=0)) >= -1e-10
    assert hist[-1].max() > 0.1, "nothing was ever driven up; gate is vacuous"


def test_static_box_alone_permits_healing_broken_control():
    _, hist, _, _ = _load_unload_history(BoxSpec(bounds={0: (0.0, CAP)}))
    assert np.min(np.diff(hist, axis=0)) < -1e-3


def test_frozen_dofs_carry_the_correct_multiplier_sign():
    """At a dynamic lower bound the KKT condition is R >= 0."""
    box = BoxSpec(bounds={0: (0.0, CAP)}, irreversible=(0,))
    nodes, hist, U, _ = _load_unload_history(box)
    _, elems, dirichlet = strip(24)
    R, _ = assemble(nodes, elems, 1, _pure_laplacian(), U, du_is_total=True)
    R[np.array(sorted(dirichlet), dtype=int)] = 0.0       # final load amplitude is 0
    lo = np.maximum(0.0, hist[-2])
    interior_frozen = (U <= lo + 1e-10) & (U > 1e-6)
    assert interior_frozen.any()
    assert R[interior_frozen].min() > -1e-9


def test_step_start_bound_is_not_ratcheted_within_the_solve():
    box = BoxSpec(bounds={0: (0.0, CAP)}, irreversible=(0,))
    nodes, elems, dirichlet = strip(24)
    U, _ = solve_steps_bounded(nodes, elems, 1, _pure_laplacian(), lambda frac: dirichlet,
                               1, box, external_fn=_external(_load_vector(nodes, elems)),
                               tol=1e-10)
    assert U.max() > 0.4


# --------------------------------------------------------------------------
# 4. BoxSpec contract
# --------------------------------------------------------------------------

def test_box_spec_requires_u_prev_for_a_dynamic_bound():
    with pytest.raises(ValueError, match="STEP-START"):
        BoxSpec(bounds={0: (0.0, 1.0)}, irreversible=(0,)).arrays(4, 1, None)


def test_box_spec_rejects_infeasible_and_contradictory_declarations():
    with pytest.raises(ValueError, match="empty box"):
        BoxSpec(bounds={0: (1.0, 0.0)})
    with pytest.raises(ValueError, match="irreversible and monotone"):
        BoxSpec(bounds={0: (0.0, 1.0)}, irreversible=(0,), monotone_decreasing=(0,))
    box = BoxSpec(bounds={0: (0.0, 1.0)}, irreversible=(0,))
    with pytest.raises(ValueError, match="empty box"):
        box.arrays(2, 1, np.array([2.0, 0.0]))     # start already above hi


def test_infeasible_dirichlet_value_is_rejected():
    nodes, elems, _ = strip(6)
    bad = {int(n): 2.0 * CAP for n in boundary_nodes(nodes)}
    with pytest.raises(ValueError, match="outside the box"):
        solve_steps_bounded(nodes, elems, 1, _pure_laplacian(), lambda frac: bad, 1,
                            BoxSpec(bounds={0: (-np.inf, CAP)}))


def test_monotone_decreasing_is_the_mirror_of_irreversible():
    lo, hi = BoxSpec(bounds={0: (-1.0, 1.0)}, monotone_decreasing=(0,)).arrays(
        3, 1, np.array([0.5, -0.2, 1.0]))
    np.testing.assert_allclose(hi, [0.5, -0.2, 1.0])
    np.testing.assert_allclose(lo, [-1.0, -1.0, -1.0])


# --------------------------------------------------------------------------
# 5. transactional state: a rejected attempt leaves nothing behind
# --------------------------------------------------------------------------

class _FlakyStatefulElement:
    """Scalar element that rejects a large increment and records its state.

    ``state`` accumulates the committed increment, like an SVARS array.  An
    increment larger than ``du_max`` raises :class:`IncrementTooLarge` (the
    ``PNEWDT < 1`` signal), so the solver must cut back rather than fail.
    The object is its own transactional state (``snapshot``/``restore``/
    ``begin_increment``) without a clock.
    """

    def __init__(self, du_max):
        self.du_max = du_max
        self.state = np.zeros(1)
        self.commits = 0
        self.rejections = 0
        self.increments = []
        self.inner = poisson_element(lambda xy: 0.0)

    def element_rk(self, ei, coords, U_e, DU_e):
        if np.max(np.abs(DU_e)) > self.du_max:
            self.rejections += 1
            raise IncrementTooLarge(f"element {ei}: |DU| > {self.du_max}")
        return self.inner(ei, coords, U_e, DU_e)

    def snapshot(self):
        return self.state.copy()

    def restore(self, snap):
        self.state = np.array(snap, dtype=float, copy=True)

    def begin_increment(self, kinc, time_begin=None, dt=None):
        self.increments.append((kinc, time_begin, dt))

    def commit(self, nodes, elems, dof_per_node, U, U_prev):
        self.state = self.state + np.abs(U - U_prev).max()
        self.commits += 1


def _flaky_solve(el, **kw):
    nodes, elems, dirichlet = strip(12)
    options = dict(external_fn=_external(_load_vector(nodes, elems)), tol=1e-10,
                   state=el, commit_fn=el.commit, max_cutbacks=4)
    options.update(kw)
    return solve_steps_bounded(nodes, elems, 1, el.element_rk, lambda frac: dirichlet, 1,
                               BoxSpec(bounds={0: (-np.inf, CAP)}), **options)


def test_failed_attempt_is_rolled_back_and_cut_back():
    el = _FlakyStatefulElement(du_max=0.30)
    U, report = _flaky_solve(el)
    rejected = [a for a in report.attempts if not a.converged]
    accepted = [a for a in report.attempts if a.converged]
    assert el.rejections > 0 and report.cutbacks >= 1 and rejected and accepted
    # Element PNEWDT requests and certified non-descent directions are cutback
    # signals; a grind to the iteration cap is not expected here.
    assert all(("smaller increment" in a.reason or "line search exhausted" in a.reason)
               for a in rejected)
    assert el.commits == len(accepted)
    dts = [inc[2] for inc in el.increments]
    assert min(dts) < 1.0 and max(dts) == 1.0
    assert np.isclose(sum(a.dt for a in accepted), 1.0)
    # TIME is the beginning of each attempt; KINC advances only on acceptance.
    accepted_increments, kinc, time = [], 1, 0.0
    for a, inc in zip(report.attempts, el.increments):
        assert inc == (kinc, pytest.approx(time), pytest.approx(a.dt))
        if a.converged:
            accepted_increments.append(inc)
            kinc, time = kinc + 1, time + a.dt
    assert [inc[0] for inc in accepted_increments] == list(range(1, len(accepted) + 1))


def test_commit_time_pnewdt_rejects_and_restores_the_bounded_attempt():
    """Assembly converges; a partially written commit alone requests a cutback."""
    class CommitCutback(_FlakyStatefulElement):
        def element_rk(self, ei, coords, U_e, DU_e):
            return self.inner(ei, coords, U_e, DU_e)

        def commit(self, nodes, elems, dof_per_node, U, U_prev):
            if np.max(np.abs(U - U_prev)) > self.du_max:
                self.state[:] = 999.0
                self.rejections += 1
                raise IncrementTooLarge("commit-side discriminator")
            super().commit(nodes, elems, dof_per_node, U, U_prev)

    element = CommitCutback(du_max=0.30)
    U, report = _flaky_solve(element)
    rejected = [a for a in report.attempts if not a.converged]
    accepted = [a for a in report.attempts if a.converged]
    assert element.rejections and any("during commit" in a.reason for a in rejected)
    assert element.commits == len(accepted)
    assert np.max(element.state) < 1.0          # the failed 999 write was rolled back
    assert U.max() <= CAP + 1e-12
    assert sum(a.dt for a in accepted) == pytest.approx(1.0)


def test_rejected_attempt_is_rolled_back_before_its_record_is_reported():
    """``attempt_callback`` sees the restored state, as in the staggered driver:
    a commit that writes 999 and then requests a cutback is already undone."""
    class CommitCutback(_FlakyStatefulElement):
        def element_rk(self, ei, coords, U_e, DU_e):
            return self.inner(ei, coords, U_e, DU_e)

        def commit(self, nodes, elems, dof_per_node, U, U_prev):
            if np.max(np.abs(U - U_prev)) > self.du_max:
                self.state[:] = 999.0
                raise IncrementTooLarge("commit-side discriminator")
            super().commit(nodes, elems, dof_per_node, U, U_prev)

    el = CommitCutback(du_max=0.30)
    seen = []
    _flaky_solve(el, attempt_callback=lambda rec: seen.append((rec.converged,
                                                               el.state.copy())))
    assert any(not ok for ok, _ in seen), "no attempt was rejected; gate is vacuous"
    committed = 0.0
    for ok, state in seen:
        if not ok:
            np.testing.assert_allclose(state, committed)
        committed = state[0]


class _TrialMutatingElement(_FlakyStatefulElement):
    """Writes its state during ASSEMBLY, the failure mode the transaction is
    for: ``state`` is the largest increment ever evaluated (a max-history)."""

    def element_rk(self, ei, coords, U_e, DU_e):
        self.state[0] = max(self.state[0], float(np.max(np.abs(DU_e))))
        return super().element_rk(ei, coords, U_e, DU_e)

    def commit(self, nodes, elems, dof_per_node, U, U_prev):
        self.commits += 1


def test_trial_state_mutation_is_undone_by_the_rollback():
    """Broken control: without the transaction the rejected trial's history stays."""
    clean = _TrialMutatingElement(du_max=0.30)
    _flaky_solve(clean)
    polluted = _TrialMutatingElement(du_max=0.30)
    _flaky_solve(polluted, state=None)
    assert polluted.rejections > 0
    assert clean.state.max() <= polluted.state.max() - 1e-6
    assert clean.state.max() <= 0.30 + 1e-12


def test_cutback_budget_is_enforced():
    el = _FlakyStatefulElement(du_max=1e-6)   # never satisfiable
    with pytest.raises(RuntimeError, match="cutbacks"):
        _flaky_solve(el, max_cutbacks=2)
    assert el.commits == 0


def test_compiled_element_state_snapshot_restore_clock_and_python_elements():
    element, kernel = recording_element(2, dt=0.5, cache_batches=True)
    python_element = poisson_element(lambda xy: 0.0)     # no svars, clock or methods
    tx = CompiledElementState([element, python_element])
    element.svars[:] = 1.0
    snap, clock = tx.snapshot(), tx.clock()
    assert snap[1] is None and clock[1] == {"time": None, "dt": None, "kinc": None}
    element.svars[0, 0] = 99.0
    tx.begin_increment(4, time_begin=2.5, dt=0.25)
    assert (element.time, element.dt, element.kinc) == (2.5, 0.25, 4)
    nodes = np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]])
    element.element_rk_batch(np.stack([nodes, nodes]), np.zeros((2, 4)), np.zeros((2, 4)))
    assert element._batch_cache is not None
    tx.restore(snap)
    assert element._batch_cache is None, "a restore must discard the cached trial"
    assert element.svars[0, 0] == 1.0, "snapshot aliased the live SVARS array"
    tx.restore_clock(clock)
    assert (element.time, element.dt, element.kinc) == (0.0, 0.5, 1)
    with pytest.raises(ValueError, match="element count"):
        tx.restore(snap[:1])


def test_bounded_rollback_restores_svars_and_clock_of_the_public_adapter():
    """A PNEWDT rejection at dt = 1 is retried at dt = 0.5 from the same clock."""
    nodes, elems, dirichlet = strip(4)
    element, kernel = recording_element(len(elems), dt=1.0, mass=1.0)
    kernel.reject = lambda time, dtime, u, du: dtime > 0.6
    element.begin_increment(5, time_begin=7.0, dt=0.3)   # a clock no attempt uses
    before = (element.time, element.dt, element.kinc, element.svars.copy())
    rolled_back = []

    def on_attempt(rec):
        if not rec.converged:
            rolled_back.append((element.time, element.dt, element.kinc,
                                element.svars.copy()))

    U, report = solve_steps_bounded(
        nodes, elems, 1, element.element_rk, lambda frac: dirichlet, 1,
        BoxSpec(bounds={0: (-np.inf, 10.0)}), batch_fn=element.element_rk_batch,
        commit_fn=element.commit_step, state=CompiledElementState([element]),
        dt=1.0, tol=1e-10, attempt_callback=on_attempt)
    assert [a.converged for a in report.attempts] == [False, True, True]
    assert len(rolled_back) == 1
    assert rolled_back[0][:3] == before[:3]
    np.testing.assert_array_equal(rolled_back[0][3], before[3])
    clocks = sorted({(c[0], c[2], c[4]) for c in kernel.calls})
    assert clocks == [(0.0, 0.5, 1), (0.0, 1.0, 1), (0.5, 0.5, 2)]
    np.testing.assert_allclose(element.svars, 1.0)      # committed DTIME 0.5 + 0.5


# --------------------------------------------------------------------------
# 6. scaling, Newton details and diagnostics
# --------------------------------------------------------------------------

def _two_field_element(scale):
    """Two decoupled Laplacians per node; component 1 is ``scale`` times weaker."""
    base = poisson_element(lambda xy: 0.0)

    def element_fn(ei, coords, U_e, DU_e=None):
        R, K = np.zeros(8), np.zeros((8, 8))
        for comp, fac in ((0, 1.0), (1, scale)):
            idx = np.arange(comp, 8, 2)
            Rc, Kc = base(ei, coords, U_e[idx])
            R[idx] = fac * Rc
            K[np.ix_(idx, idx)] = fac * Kc
        return R, K

    return element_fn


def test_residual_scale_makes_a_weak_field_enforceable():
    scale = 1e-10
    nodes, elems, _ = strip(8)
    F1 = _load_vector(nodes, elems)
    F_vec = np.zeros(2 * len(nodes))
    F_vec[0::2], F_vec[1::2] = F1, scale * F1
    edge = np.nonzero((nodes[:, 0] <= 1e-12) | (nodes[:, 0] >= 1.0 - 1e-12))[0]
    dirichlet = {int(2 * n + c): 0.0 for n in edge for c in (0, 1)}
    U, rep = solve_steps_bounded(
        nodes, elems, 2, _two_field_element(scale), lambda frac: dirichlet, 1,
        BoxSpec(bounds={0: (-np.inf, CAP), 1: (-np.inf, CAP)}),
        external_fn=_external(F_vec), tol=1e-8, residual_scale={1: 1.0 / scale})
    np.testing.assert_allclose(U[1::2], U[0::2], atol=1e-6)
    assert rep.attempts[-1].kkt_residual < 1e-8


def test_field_scale_controls_bound_activity_detection():
    lo, hi = np.zeros(3), np.full(3, 1.0e-9)
    U = np.array([0.0, 1.0e-9, 5.0e-10])
    R = np.array([-1.0, 1.0, 1.0])
    assert projected_residual(R, U, lo, hi, at_tol=0.0)[0][2] == 1.0
    assert projected_residual(R, U, lo, hi, at_tol=1.0e-9)[0][2] == 0.0


def test_initial_bound_adjacent_kkt_point_is_returned_exactly_on_bound():
    result = _bounded_newton(
        lambda U, DU: (np.array([1.0]), np.array([[1.0]])), np.array([5.0e-13]),
        np.array([0.0]), np.array([], dtype=int), np.array([False]),
        np.array([0.0]), np.array([1.0]), w_res=np.ones(1), at_tol=1.0e-12,
        max_newton=2, tol=1.0e-12, linear_solve=lambda K, R: np.linalg.solve(K, R),
        line_search=True, context="initial-snap regression")
    assert result["converged"] and result["U"][0] == 0.0


def test_exhausted_line_search_rejects_iteration_instead_of_crawling():
    """A certified non-descent direction fails the attempt at once, with U
    unchanged, instead of applying a 2^-30 update and crawling."""
    n_assemblies = []

    def asm(U, DU):
        n_assemblies.append(1)
        return np.array([U[0]]), sp.csr_matrix([[-1.0]])     # du = +U: ascent

    result = _bounded_newton(
        asm, np.array([1.0]), np.array([0.0]), np.array([], dtype=int),
        np.array([False]), np.array([-np.inf]), np.array([np.inf]), w_res=np.ones(1),
        at_tol=1.0e-12, max_newton=5, tol=1.0e-10, linear_solve=_scipy_solve,
        line_search=True, context="line-search-exhaustion regression")
    assert not result["converged"] and result["U"][0] == 1.0
    assert result["iterations"] == 1 and result["residual"] == 1.0
    assert "line search exhausted" in result["reason"]
    assert len(n_assemblies) == 31


def test_bounded_newton_assesses_the_state_of_its_final_update():
    """Defect regression: the last permitted update reached the exact root of
    a linear problem and was returned as a failure without being evaluated."""
    calls = []

    def asm(U, DU):
        calls.append(U[0])
        return np.array([2.0 * U[0] - 1.0]), sp.csr_matrix([[2.0]])

    result = _bounded_newton(
        asm, np.array([0.0]), np.array([0.0]), np.array([], dtype=int),
        np.array([False]), np.array([-np.inf]), np.array([np.inf]), w_res=np.ones(1),
        at_tol=0.0, max_newton=1, tol=1e-12, linear_solve=_scipy_solve,
        line_search=True, context="final update")
    assert result["converged"] and result["U"][0] == 0.5
    assert result["iterations"] == 2 and result["residual"] == 0.0
    np.testing.assert_array_equal(result["last_update"], [0.5])
    # Still no extra update: one iterate, one trial, one final assessment.
    assert calls == [0.0, 0.5, 0.5]

    def nonlinear(U, DU):
        return np.array([U[0] ** 3 + U[0] - 0.5]), sp.csr_matrix([[3 * U[0] ** 2 + 1]])

    failed = _bounded_newton(
        nonlinear, np.array([0.0]), np.array([0.0]), np.array([], dtype=int),
        np.array([False]), np.array([-np.inf]), np.array([np.inf]), w_res=np.ones(1),
        at_tol=0.0, max_newton=1, tol=1e-12, linear_solve=_scipy_solve,
        line_search=True, context="final update")
    assert not failed["converged"] and failed["iterations"] == 1
    assert failed["residual"] == pytest.approx(0.125, rel=1e-12)  # assessed at u = 0.5


def test_solve_steps_bounded_accepts_a_step_converged_by_its_last_update():
    """Through the driver: a linear solve with max_newton=1 needs no cutback."""
    nodes, elems, U, report = _solve_bounded(8, box=BoxSpec(bounds={0: (-np.inf, 10.0)}),
                                             max_newton=1)
    assert report.cutbacks == 0 and report.attempts[-1].iterations == 2
    assert report.attempts[-1].kkt_residual < 1e-10


@pytest.mark.parametrize("cutback", [None, False, True])
def test_element_rejection_during_line_search_is_an_opt_in_cutback(cutback):
    """Default backtracking recovers; the opt-in returns the material reason."""
    trials = []

    def asm(U, DU):
        x = float(U[0])
        trials.append(x)
        if x > 2.0:
            raise IncrementTooLarge("local material solve declined this trial")
        return np.array([x * x - 1.0]), sp.csr_matrix([[2.0 * x]])

    kw = {} if cutback is None else {"cutback_on_element_rejection": cutback}
    U0 = np.array([0.1])
    result = _bounded_newton(
        asm, U0, U0.copy(), np.array([], dtype=int), np.array([False]),
        np.array([0.0]), np.array([10.0]), w_res=np.ones(1), at_tol=0.0,
        max_newton=12, tol=1e-12, linear_solve=_scipy_solve,
        line_search=True, context="material line-search rejection", **kw)
    assert trials[1] > 2.0   # the first Newton trial (5.05) is rejected
    np.testing.assert_array_equal(U0, [0.1])
    if cutback:
        assert not result["converged"] and result["iterations"] == 1
        assert len(trials) == 2
        np.testing.assert_array_equal(result["U"], U0)
        assert np.isinf(result["residual"])
        assert "element requested a smaller increment during line search" in result["reason"]
        assert "local material solve declined this trial" in result["reason"]
    else:
        assert result["converged"] and len(trials) > 2
        assert result["U"][0] == pytest.approx(1.0, abs=1e-12)
        assert result["residual"] < 1e-12


def test_nonfinite_trial_is_an_infinite_kkt_residual_but_an_iterate_raises():
    from abaqus_ufl.fe import NonfiniteElementOutput

    def asm(U, DU):
        if U[0] > 2.0:
            raise NonfiniteElementOutput("trial outside the element's domain")
        return np.array([U[0] * U[0] - 1.0]), sp.csr_matrix([[2.0 * U[0]]])

    common = dict(con=np.array([], dtype=int), is_con=np.array([False]),
                  lo=np.array([0.0]), hi=np.array([10.0]), w_res=np.ones(1),
                  at_tol=0.0, max_newton=12, tol=1e-12, linear_solve=_scipy_solve,
                  line_search=True, context="nonfinite trial")
    result = _bounded_newton(asm, np.array([0.1]), np.array([0.1]), **common)
    assert result["converged"] and result["U"][0] == pytest.approx(1.0, abs=1e-12)
    with pytest.raises(NonfiniteElementOutput):
        _bounded_newton(asm, np.array([3.0]), np.array([3.0]), **common)

    def buggy(U, DU):
        if U[0] > 2.0:
            raise RuntimeError("element bug")
        return asm(U, DU)

    with pytest.raises(RuntimeError, match="element bug"):
        _bounded_newton(buggy, np.array([0.1]), np.array([0.1]), **common)


def test_report_records_active_counts_and_accepted_steps():
    _, _, U, report = _solve_bounded(24, n_steps=3)
    assert report.steps_accepted == 3
    assert report.attempts[-1].n_active_upper > 0
    summary = report.summary()
    assert summary["accepted"] == 3 and summary["cutbacks"] == 0
    assert summary["max_kkt_residual"] < 1e-10


def test_degenerate_touching_bound_still_terminates():
    """Cap exactly at the unconstrained peak: a single borderline active DOF."""
    nodes, elems, dirichlet = strip(16)
    box = BoxSpec(bounds={0: (-np.inf, F_SOURCE / 8.0)})
    U, report = solve_steps_bounded(
        nodes, elems, 1, _pure_laplacian(), lambda frac: dirichlet, 1, box,
        external_fn=_external(_load_vector(nodes, elems)), tol=1e-10, max_newton=25)
    assert report.attempts[-1].converged
    assert U.max() <= F_SOURCE / 8.0 + 1e-12


# --------------------------------------------------------------------------
# 7. prescribed values on bounded routes
# --------------------------------------------------------------------------

def _coupled_to_prescribed(ei, coords, U_e, DU_e):
    """Component 1 must equal 1e13 times the prescribed component 0, so any
    change of the prescribed value is magnified into an O(1) error."""
    return (np.array([0.0, U_e[1] - 1e13 * U_e[0]]),
            np.array([[1.0, 0.0], [-1e13, 1.0]]))


def _bounded_route(route, element_fn, dirichlet, dpn=2, **kw):
    nodes, elems = np.zeros((1, 2)), np.array([[0]])
    box = BoxSpec({0: (0.0, 1.0)})
    if route == "bounded":
        return solve_steps_bounded(nodes, elems, dpn, element_fn, lambda f: dirichlet, 1,
                                   box, tol=1e-10, **kw)[0]
    options = {"feasible": dict(box=box, feasible_trials=True),
               "bounded_monolithic": dict(box=box, newton_method="bounded_monolithic"),
               "static_clip": dict(field_bounds={0: (0.0, 1.0)}),
               "box_clip": dict(box=box)}[route]
    return solve_steps_staggered(nodes, elems, dpn, element_fn, lambda f: dirichlet,
                                 [list(range(dpn))], 1, tol=1e-10, **options, **kw)


@pytest.mark.parametrize("route", ["bounded", "feasible", "bounded_monolithic",
                                   "static_clip", "box_clip"])
@pytest.mark.parametrize("value", [5e-13, 1.0 + 5e-13])
def test_bounded_routes_keep_prescribed_values_exactly(route, value):
    """A prescribed DOF is never snapped or clipped onto a bound: a value
    within the activity tolerance of a bound, or inside the 1e-12 validation
    slack above it, is kept bit for bit."""
    U = _bounded_route(route, _coupled_to_prescribed, {0: value})
    assert U[0] == value
    assert U[1] == pytest.approx(1e13 * value, rel=1e-15)


@pytest.mark.parametrize("route", ["bounded", "feasible"])
def test_prescribed_value_inside_a_large_activity_tolerance_is_kept(route):
    linear = lambda ei, c, u, du: (u - 2.0, np.eye(len(u)))
    U = _bounded_route(route, linear, {0: 0.05}, dpn=1, bound_atol=0.1)
    assert U[0] == 0.05


@pytest.mark.parametrize("route", ["bounded", "feasible", "bounded_monolithic",
                                   "static_clip", "box_clip"])
def test_incompatible_prescribed_value_is_rejected_on_every_bounded_route(route):
    """Prescribing 2 inside the bounds [0, 1] must raise; the default static
    clip used to return an accepted 1."""
    linear = lambda ei, c, u, du: (u - 2.0, np.eye(len(u)))
    with pytest.raises(ValueError, match="outside the box"):
        _bounded_route(route, linear, {0: 2.0}, dpn=1)


# --------------------------------------------------------------------------
# 8. an exception that escapes an attempt is raised after the rollback
# --------------------------------------------------------------------------

def _clocked_element():
    """Public adapter with committed SVARS 0.25 and a clock no attempt uses."""
    nodes, elems, _ = strip(2)
    element, kernel = recording_element(len(elems), dt=0.3, time=7.0)
    element.begin_increment(5, time_begin=7.0, dt=0.3)
    element.svars[:] = 0.25
    return element, kernel


def _assert_untouched(element):
    np.testing.assert_array_equal(element.svars, 0.25)
    assert (element.time, element.dt, element.kinc) == (7.0, 0.3, 5)


def _drive(route, element, dirichlet=None, **kw):
    nodes, elems, bc = strip(2)
    bc = bc if dirichlet is None else dirichlet
    options = dict(state=CompiledElementState([element]),
                   batch_fn=kw.pop("batch_fn", element.element_rk_batch), dt=1.0, tol=1e-10)
    options.update(kw)
    if route == "bounded":
        return solve_steps_bounded(nodes, elems, 1, element.element_rk, lambda f: bc, 1,
                                   BoxSpec({0: (-np.inf, 10.0)}), **options)
    return solve_steps_staggered(nodes, elems, 1, element.element_rk, lambda f: bc, [[0]], 1,
                                 **options)


@pytest.mark.parametrize("route", ["bounded", "staggered"])
@pytest.mark.parametrize("error", ["NonfiniteElementOutput", "ValueError"])
def test_an_exception_escaping_a_commit_is_raised_after_rollback(route, error):
    from abaqus_ufl.fe import NonfiniteElementOutput
    error = {"NonfiniteElementOutput": NonfiniteElementOutput, "ValueError": ValueError}[error]
    element, _ = _clocked_element()

    def commit(*args):
        element.svars[:] = 999.0            # a partial write, then a failure
        raise error("commit failure")

    with pytest.raises(error, match="commit failure"):
        _drive(route, element, commit_fn=commit)
    _assert_untouched(element)


@pytest.mark.parametrize("route", ["bounded", "staggered"])
def test_an_element_error_escaping_a_trial_is_raised_after_rollback(route):
    element, _ = _clocked_element()

    def batch(coords, U, DU):
        if np.max(np.abs(DU)) > 0.0:
            raise RuntimeError("element bug")
        return element.element_rk_batch(coords, U, DU)

    with pytest.raises(RuntimeError, match="element bug"):
        _drive(route, element, batch_fn=batch)
    _assert_untouched(element)


def test_a_bounded_nonfinite_iterate_is_raised_after_rollback():
    from abaqus_ufl.fe import NonfiniteElementOutput
    element, kernel = _clocked_element()
    kernel.nonfinite = lambda *args: True
    with pytest.raises(NonfiniteElementOutput):
        _drive("bounded", element)
    _assert_untouched(element)


@pytest.mark.parametrize("route", ["bounded", "staggered"])
def test_an_invalid_dirichlet_value_is_raised_after_rollback(route):
    element, _ = _clocked_element()
    nodes, _, bc = strip(2)
    with pytest.raises(ValueError, match="outside the box|out of range"):
        _drive(route, element, dirichlet={**bc, 10 ** 6: 0.0} if route == "staggered"
               else {**bc, 0: 20.0})
    _assert_untouched(element)
