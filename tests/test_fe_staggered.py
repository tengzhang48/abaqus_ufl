"""Contract gates for the transactional, scaled staggered solver.

1. **Irreversibility.**  A static box lets a damage field HEAL whenever its
   unconstrained equilibrium falls, which happens at *frozen* history whenever
   the toughness rises.  The analytic element below reproduces that mechanism
   in closed form, so these gates are a solver statement, not a physics one.
   A monotone history variable does NOT imply a monotone damage field, and
   ``max_i d_i`` does not see the defect (it *rises* here while a third of the
   nodes heal), so every gate measures ``min_i (d_i[n+1] - d_i[n])``.
2. **Transactional adaptive stepping**: per-field scaled residual/update gates,
   the fresh post-clamp certificate, budgeted active-mask repair, PNEWDT
   cutbacks (from assembly AND from commit), rollback of ``U``, SVARS and the
   increment clock, and the TIME/DTIME/KINC passed to the public adapter.
3. **Restart** from an in-memory restart dictionary, including the dynamic
   bounds reconstructed from the restored state.
"""

import inspect
from dataclasses import asdict

import numpy as np
import pytest

from abaqus_ufl.fe import (BoxSpec, CompiledElementState, IncrementTooLarge,
                           StaggeredSolveFailed, assemble,
                           projected_residual, solve_steps_bounded,
                           solve_steps_staggered, structured_quad_mesh)
from abaqus_ufl.fe.stepping import _scipy_solve
from fe_stepping_toys import poisson_element, recording_element, strip

F_SOURCE = 8.0
CAP = 0.5

# analytic-element component layout
UC, TC, SC, DC = 0, 1, 2, 3
DPN = 4


class _Clock:
    """The public adapter's increment clock, for Python test elements."""

    def _init_clock(self):
        self.time, self.dt, self.kinc = 0.0, 1.0, 1

    def begin_increment(self, kinc, time_begin=None, dt=None):
        self.dt = self.dt if dt is None else float(dt)
        self.time = (kinc - 1) * self.dt if time_begin is None else float(time_begin)
        self.kinc = int(kinc)


# --------------------------------------------------------------------------
# analytic healing-damage element
# --------------------------------------------------------------------------

class HealingDamageElement(_Clock):
    """Lumped AT2-style damage element with PRESCRIBED drive and toughness.

    Per node ``psi(d) = (1 - d)^2 H + Gc/(2 ell) d^2`` (+ optional gradient),
    so ``R_d = Gc/ell * d - 2 (1 - d) H`` and the unconstrained nodal
    equilibrium is ``d* = 2 H / (Gc/ell + 2 H)``.  ``H`` is the prescribed
    nodal field ``S`` and ``Gc(T) = Gc0 (1 + gamma (T_ref - T))`` with the
    prescribed field ``T``: RAISING Gc AT FROZEN H LOWERS ``d*``.

    A companion field obeys ``R_u = k_u (u - beta d)``, so clamping ``d``
    leaves a residual only a re-sweep against the active mask can clear.
    ``grad`` adds ``grad * Lap(d)``, spreading the bound interaction.  The
    element carries SVARS and the public clock and can decline an increment
    (``max_dt``) or a commit (``reject_commit_above``) as ``PNEWDT < 1``.
    """

    def __init__(self, nodes, elems, *, Gc0=1.0, ell=1.0, gamma=0.5, T_ref=0.0,
                 k_u=1.0, beta=1.0, grad=0.0, max_dt=None, reject_commit_above=None):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elems = np.asarray(elems, dtype=int)
        self.Gc0, self.ell, self.gamma, self.T_ref = Gc0, ell, gamma, T_ref
        self.k_u, self.beta, self.grad = k_u, beta, grad
        self.max_dt = max_dt
        self.reject_commit_above = reject_commit_above
        self._lap = poisson_element(lambda xy: 0.0)
        self.svars = np.zeros((len(self.elems), 1))
        self._init_clock()
        self.commits = 0
        self.rejections = 0
        self.committed_moves = []

    def equilibrium(self, S, T):
        Gc = self.Gc0 * (1.0 + self.gamma * (self.T_ref - np.asarray(T, float)))
        S = np.asarray(S, dtype=float)
        return 2.0 * S / (Gc / self.ell + 2.0 * S)

    @staticmethod
    def _area(coords):
        x, y = coords[:, 0], coords[:, 1]
        return 0.5 * abs(x[0] * (y[1] - y[3]) + x[1] * (y[2] - y[0])
                         + x[2] * (y[3] - y[1]) + x[3] * (y[0] - y[2]))

    def element_rk(self, ei, coords, U_e, DU_e=None):
        if self.max_dt is not None and self.dt > self.max_dt + 1e-12:
            self.rejections += 1
            raise IncrementTooLarge(f"dt={self.dt:.4g} > {self.max_dt:.4g}")
        coords = np.asarray(coords, dtype=float)
        U_e = np.asarray(U_e, dtype=float)
        w = self._area(coords) / 4.0
        R = np.zeros(4 * DPN)
        K = np.zeros((4 * DPN, 4 * DPN))
        for a in range(4):
            iu, iT, iS, idd = 4 * a + UC, 4 * a + TC, 4 * a + SC, 4 * a + DC
            u, T, S, d = U_e[iu], U_e[iT], U_e[iS], U_e[idd]
            Gc = self.Gc0 * (1.0 + self.gamma * (self.T_ref - T))
            R[iu] = w * self.k_u * (u - self.beta * d)
            K[iu, iu] = w * self.k_u
            K[iu, idd] = -w * self.k_u * self.beta
            R[idd] = w * ((Gc / self.ell) * d - 2.0 * (1.0 - d) * S)
            K[idd, idd] = w * (Gc / self.ell + 2.0 * S)
            K[idd, iT] = -w * (d / self.ell) * self.Gc0 * self.gamma
            K[idd, iS] = -w * 2.0 * (1.0 - d)
        if self.grad:
            Rg, Kg = self._lap(ei, coords, U_e[DC::DPN])
            R[DC::DPN] += self.grad * np.asarray(Rg)
            for a in range(4):
                for b in range(4):
                    K[4 * a + DC, 4 * b + DC] += self.grad * Kg[a, b]
        return R, K

    def commit(self, nodes, elems, dof_per_node, U, U_prev):
        moved = np.abs(U[DC::DPN] - U_prev[DC::DPN])[np.asarray(elems)].max(axis=1)
        if self.reject_commit_above is not None and moved.max() > self.reject_commit_above:
            self.rejections += 1
            raise IncrementTooLarge(f"commit: |dd|={moved.max():.4g}")
        self.svars[:, 0] += moved
        self.committed_moves.append(float(moved.max()))
        self.commits += 1


def _heal_bc(nodes, S_left, S_mid, S_right, T_left, T_mid, T_right):
    """Prescribe (T, S) at EVERY node; u and d stay free."""
    x = nodes[:, 0]
    left, right = x < 0.25, x > 0.75
    S = np.where(left, S_left, np.where(right, S_right, S_mid))
    T = np.where(left, T_left, np.where(right, T_right, T_mid))
    bc = {}
    for i in range(len(nodes)):
        bc[i * DPN + TC] = float(T[i])
        bc[i * DPN + SC] = float(S[i])
    return bc


def _heal_dirichlet(nodes):
    def dirichlet_fn(frac):
        if frac <= 0.5 + 1e-14:
            return _heal_bc(nodes, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
        return _heal_bc(nodes, 3.0, 1.0, 1.0, 0.0, 0.0, -1.0)
    return dirichlet_fn


def _heal_history(**solver_kw):
    """Step 1 equilibrates uniformly; step 2 raises H on the left and Gc on the
    right, so the right column's unconstrained equilibrium FALLS."""
    nodes, elems = structured_quad_mesh(2, 2)
    el = solver_kw.pop("element", None) or HealingDamageElement(nodes, elems)
    history = []
    U = solve_steps_staggered(
        nodes, elems, DPN, el.element_rk, _heal_dirichlet(nodes),
        [[UC, TC, SC], [DC]], 2, tol=1e-11,
        step_callback=lambda s, f, Uc, rep: history.append(Uc[DC::DPN].copy()),
        **solver_kw)
    return nodes, el, np.array(history), U


def _column_masks(nodes):
    x = nodes[:, 0]
    return x < 0.25, (x >= 0.25) & (x <= 0.75), x > 0.75


# --------------------------------------------------------------------------
# 0. signature stability
# --------------------------------------------------------------------------

#: Positional parameters; every later parameter is keyword-only, so an older
#: positional call can never silently rebind (``update_tol`` landing on
#: ``line_search``, ``box`` on ``field_bounds``, ...).
_STAGGERED_POSITIONAL = [
    "nodes", "elems", "dof_per_node", "element_fn", "dirichlet_fn",
    "field_components", "n_steps", "max_stagger", "max_field_newton",
    "tol", "line_search", "commit_fn", "linear_solve", "batch_fn",
    "field_bounds",
]
_BOUNDED_POSITIONAL = ["nodes", "elems", "dof_per_node", "element_fn",
                       "dirichlet_fn", "n_steps", "box"]


def test_staggered_positional_signature_is_stable():
    params = inspect.signature(solve_steps_staggered).parameters
    assert [n for n, p in params.items() if p.kind is p.POSITIONAL_OR_KEYWORD] \
        == _STAGGERED_POSITIONAL
    assert {n for n, p in params.items() if p.kind is p.KEYWORD_ONLY} == {
        "update_tol", "box", "step_callback", "verbose", "U0",
        "residual_scale", "field_scale", "state", "dt", "max_cutbacks",
        "min_dt", "cutback_factor", "grow_after", "max_projection_rounds",
        "attempt_callback", "accepted_callback", "restart", "early_stop",
        "return_report", "feasible_trials", "bound_atol", "newton_method",
        "cutback_on_element_rejection"}


def test_bounded_positional_signature_and_keyword_order():
    params = list(inspect.signature(solve_steps_bounded).parameters.items())
    assert [n for n, p in params if p.kind is p.POSITIONAL_OR_KEYWORD] == _BOUNDED_POSITIONAL
    assert [n for n, p in params if p.kind is p.KEYWORD_ONLY] == [
        "batch_fn", "U0", "commit_fn", "state", "max_newton", "tol",
        "linear_solve", "residual_scale", "field_scale", "bound_atol",
        "line_search", "max_cutbacks", "dt", "external_fn", "step_callback",
        "attempt_callback", "verbose", "early_stop"]


def test_all_positional_staggered_call_binds_correctly():
    nodes, elems, dirichlet = strip(12)
    seen = []
    U = solve_steps_staggered(
        nodes, elems, 1, poisson_element(lambda xy: F_SOURCE), lambda frac: dirichlet,
        [[0]], 1, 60, 15, 1e-10, True,
        lambda nodes_, elems_, dpn, U, U_prev: seen.append(float(U.max())),
        _scipy_solve, None, {0: (-np.inf, 2.0)})      # a NON-active clip
    assert seen, "commit_fn did not bind to the 12th positional argument"
    assert 0.5 < U.max() <= 2.0 + 1e-12
    exact = 0.5 * F_SOURCE * nodes[:, 0] * (1.0 - nodes[:, 0])
    assert np.max(np.abs(U - exact)) < 2.0e-3


def test_positional_field_bounds_reaches_the_certificate():
    nodes, elems, dirichlet = strip(24)
    with pytest.raises(StaggeredSolveFailed, match="post-clamp projected residual failed"):
        solve_steps_staggered(
            nodes, elems, 1, poisson_element(lambda xy: F_SOURCE), lambda frac: dirichlet,
            [[0]], 1, 60, 15, 1e-10, True, None, _scipy_solve, None, {0: (-np.inf, CAP)})


def test_additive_parameters_are_keyword_only():
    nodes, elems, dirichlet = strip(6)
    with pytest.raises(TypeError, match="positional argument"):
        solve_steps_staggered(
            nodes, elems, 1, poisson_element(lambda xy: F_SOURCE), lambda frac: dirichlet,
            [[0]], 1, 60, 15, 1e-10, True, None, _scipy_solve, None, None, 1e-8)


# --------------------------------------------------------------------------
# 1. field groups (defect regression) and irreversibility
# --------------------------------------------------------------------------

def _two_component_linear(ei, coords, U_e, DU_e):
    """Per node R = U - 2 for both components (decoupled)."""
    return np.asarray(U_e, dtype=float) - 2.0, np.eye(len(U_e))


@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_a_free_dof_outside_every_field_group_is_rejected(method):
    """Defect regression: omitting a component used to 'converge' with that
    component's residual untouched (norm 2 per DOF here)."""
    nodes, elems = structured_quad_mesh(1, 1)
    with pytest.raises(ValueError, match="belongs to no field group"):
        solve_steps_staggered(nodes, elems, 2, _two_component_linear, lambda f: {},
                              [[0]], 1, tol=1e-10, newton_method=method)
    # A component that is prescribed at every node needs no group.
    U = solve_steps_staggered(nodes, elems, 2, _two_component_linear,
                              lambda f: {2 * n + 1: 2.0 for n in range(4)},
                              [[0]], 1, tol=1e-10, newton_method=method)
    np.testing.assert_allclose(U, 2.0)


@pytest.mark.parametrize("groups, message", [
    ([[0], [0, 1]], "more than one field group"),
    ([[0, 0], [1]], "more than one field group"),
    ([[0], []], "empty"),
    ([[0], [2]], "not a component index"),
    ([[0], [True]], "not a component index"),
    ([], "at least one"),
])
def test_invalid_or_overlapping_field_groups_are_rejected(groups, message):
    nodes, elems = structured_quad_mesh(1, 1)
    with pytest.raises(ValueError, match=message):
        solve_steps_staggered(nodes, elems, 2, _two_component_linear, lambda f: {},
                              groups, 1)


def test_static_box_heals_broken_control():
    """BROKEN CONTROL: the static ``field_bounds`` route lets a third of the
    nodes heal, exactly to the closed-form unconstrained equilibrium."""
    nodes, el, hist, _ = _heal_history(field_bounds={DC: (0.0, 1.0)})
    left, mid, right = _column_masks(nodes)
    delta = hist[1] - hist[0]
    assert delta.min() < -1e-3
    np.testing.assert_allclose(hist[0], el.equilibrium(1.0, 0.0), atol=1e-9)
    np.testing.assert_allclose(hist[1][right], el.equilibrium(1.0, -1.0), atol=1e-9)
    np.testing.assert_allclose(hist[1][left], el.equilibrium(3.0, 0.0), atol=1e-9)
    np.testing.assert_allclose(delta[mid], 0.0, atol=1e-9)
    assert hist[1].max() > hist[0].max() + 1e-3       # d_max rises: blind monitor


@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_irreversible_box_blocks_healing_with_nonuniform_nodal_changes(method):
    box = BoxSpec(bounds={DC: (0.0, 1.0)}, irreversible=(DC,))
    nodes, el, hist, U = _heal_history(box=box, max_projection_rounds=4,
                                       newton_method=method)
    left, mid, right = _column_masks(nodes)
    delta = hist[1] - hist[0]
    assert delta.min() >= -1e-12, f"a node healed by {delta.min():.3e}"
    np.testing.assert_array_equal(hist[1][right], hist[0][right])
    np.testing.assert_array_equal(hist[1][mid], hist[0][mid])
    assert delta[left].min() > 1e-2
    np.testing.assert_allclose(hist[1][left], el.equilibrium(3.0, 0.0), atol=1e-9)
    np.testing.assert_allclose(U[UC::DPN], el.beta * U[DC::DPN], atol=1e-9)


@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_irreversible_accepted_state_is_a_kkt_point_not_a_clip(method):
    box = BoxSpec(bounds={DC: (0.0, 1.0)}, irreversible=(DC,))
    nodes, el, hist, U = _heal_history(box=box, max_projection_rounds=4,
                                       newton_method=method)
    _, _, right = _column_masks(nodes)
    R, _ = assemble(nodes, el.elems, DPN, el.element_rk, U, U - U)
    lo, hi = np.full(len(U), -np.inf), np.full(len(U), np.inf)
    lo[DC::DPN], hi[DC::DPN] = hist[0], 1.0
    R[np.array(sorted(_heal_bc(nodes, 3.0, 1.0, 1.0, 0.0, 0.0, -1.0)), dtype=int)] = 0.0
    Rp, at_lo, _ = projected_residual(R, U, lo, hi)
    assert np.linalg.norm(Rp) < 1e-9
    held = np.zeros(len(U), dtype=bool)
    held[DC::DPN] = right
    assert at_lo[held].all() and R[held].min() > 1e-6
    assert np.linalg.norm(R) > 1e-3


def test_irreversible_without_repair_rounds_refuses_rather_than_heals():
    box = BoxSpec(bounds={DC: (0.0, 1.0)}, irreversible=(DC,))
    with pytest.raises(StaggeredSolveFailed, match="post-clamp projected residual failed"):
        _heal_history(box=box, max_projection_rounds=0)


@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_irreversible_box_with_gradient_coupling_is_spatially_selective(method):
    nodes, elems = structured_quad_mesh(2, 2)
    el = HealingDamageElement(nodes, elems, grad=0.05)
    box = BoxSpec(bounds={DC: (0.0, 1.0)}, irreversible=(DC,))
    nodes, el, hist, U = _heal_history(element=el, box=box, max_projection_rounds=6,
                                       newton_method=method)
    delta = hist[1] - hist[0]
    assert delta.min() >= -1e-12
    at_bound = np.isclose(delta, 0.0, atol=1e-14)
    assert at_bound.any() and not at_bound.all()
    np.testing.assert_allclose(U[UC::DPN], el.beta * U[DC::DPN], atol=1e-9)


def test_monotone_decreasing_is_supported_by_the_staggered_route():
    nodes, elems = structured_quad_mesh(2, 2)
    el = HealingDamageElement(nodes, elems)
    d0, d_target = float(el.equilibrium(1.0, -1.0)), float(el.equilibrium(3.0, 0.0))
    assert d_target > d0 + 1e-2
    U0 = np.zeros(len(nodes) * DPN)
    U0[DC::DPN], U0[UC::DPN], U0[SC::DPN], U0[TC::DPN] = d0, el.beta * d0, 1.0, -1.0

    def run(box):
        return solve_steps_staggered(
            nodes, elems, DPN, el.element_rk,
            lambda frac: _heal_bc(nodes, 3.0, 3.0, 3.0, 0.0, 0.0, 0.0),
            [[UC, TC, SC], [DC]], 1, tol=1e-11, box=box, max_projection_rounds=4,
            U0=U0.copy())

    held = run(BoxSpec(bounds={DC: (0.0, 1.0)}, monotone_decreasing=(DC,)))
    free = run(BoxSpec(bounds={DC: (0.0, 1.0)}))
    np.testing.assert_allclose(free[DC::DPN], d_target, atol=1e-9)
    np.testing.assert_array_equal(held[DC::DPN], np.full(len(nodes), d0))
    np.testing.assert_allclose(held[UC::DPN], el.beta * held[DC::DPN], atol=1e-9)


def test_box_and_field_bounds_are_mutually_exclusive():
    nodes, elems = structured_quad_mesh(2, 2)
    el = HealingDamageElement(nodes, elems)
    with pytest.raises(ValueError, match="mutually exclusive"):
        solve_steps_staggered(nodes, elems, DPN, el.element_rk, lambda frac: {},
                              [[UC, TC, SC], [DC]], 1, box=BoxSpec(bounds={DC: (0.0, 1.0)}),
                              field_bounds={DC: (0.0, 1.0)})
    with pytest.raises(TypeError, match="must be a BoxSpec"):
        solve_steps_staggered(nodes, elems, DPN, el.element_rk, lambda frac: {},
                              [[UC, TC, SC], [DC]], 1, box={DC: (0.0, 1.0)})


def test_infeasible_prescribed_value_under_a_box_is_rejected():
    nodes, elems = structured_quad_mesh(2, 2)
    el = HealingDamageElement(nodes, elems)

    def dirichlet_fn(frac):
        bc = _heal_bc(nodes, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
        bc[DC] = 2.0              # outside [0, 1]
        return bc

    with pytest.raises(ValueError, match="outside the box"):
        solve_steps_staggered(nodes, elems, DPN, el.element_rk, dirichlet_fn,
                              [[UC, TC, SC], [DC]], 1, box=BoxSpec(bounds={DC: (0.0, 1.0)}))


# --------------------------------------------------------------------------
# 2. transactional adaptive stepping
# --------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_irreversible_bounds_survive_a_cutback_and_are_not_ratcheted(method):
    """A rejected attempt restores U, SVARS and the clock, and the retried
    sub-increments re-derive the SAME step-start bound."""
    nodes, elems = structured_quad_mesh(2, 2)
    el = HealingDamageElement(nodes, elems, max_dt=0.6)
    state = CompiledElementState([el])
    el.time, el.kinc = 4.0, 7
    last = {"svars": el.svars.copy(), "clock": (el.time, el.dt, el.kinc)}
    seen = []

    def on_attempt(rec):
        # Rollback is complete before the callback: a rejected attempt leaves
        # the SVARS and clock of the last accepted state (or the initial one).
        seen.append(rec)
        if rec.converged:
            last.update(svars=el.svars.copy(), clock=(el.time, el.dt, el.kinc))
        else:
            np.testing.assert_array_equal(el.svars, last["svars"])
            assert (el.time, el.dt, el.kinc) == last["clock"]

    box = BoxSpec(bounds={DC: (0.0, 1.0)}, irreversible=(DC,))
    nodes, el, hist, U = _heal_history(
        element=el, box=box, max_projection_rounds=4, state=state, dt=1.0,
        max_cutbacks=3, commit_fn=el.commit, attempt_callback=on_attempt,
        newton_method=method)
    rejected = [a for a in seen if not a.converged]
    accepted = [a for a in seen if a.converged]
    assert rejected and all("smaller increment" in a.reason for a in rejected)
    assert el.commits == len(accepted)
    left, _, right = _column_masks(nodes)
    delta = hist[1] - hist[0]
    assert delta.min() >= -1e-12
    np.testing.assert_array_equal(hist[1][right], hist[0][right])
    assert delta[left].min() > 1e-2


@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_pnewdt_at_commit_rejects_the_attempt_transactionally(method):
    nodes, elems = structured_quad_mesh(2, 2)
    el = HealingDamageElement(nodes, elems, reject_commit_above=0.45)
    attempts = []
    U, report = solve_steps_staggered(
        nodes, elems, DPN, el.element_rk,
        lambda frac: _heal_bc(nodes, 3.0 * frac, 3.0 * frac, 3.0 * frac, 0.0, 0.0, 0.0),
        [[UC, TC, SC], [DC]], 1, tol=1e-11,
        box=BoxSpec(bounds={DC: (0.0, 1.0)}, irreversible=(DC,)),
        max_projection_rounds=4, state=CompiledElementState([el]), dt=1.0,
        max_cutbacks=4, commit_fn=el.commit, attempt_callback=attempts.append,
        return_report=True, newton_method=method)
    rejected = [a for a in attempts if not a.converged]
    assert rejected and all("smaller increment" in a.reason for a in rejected)
    assert report.cutbacks == len(rejected)
    assert el.committed_moves and max(el.committed_moves) <= 0.45 + 1e-12
    assert len(el.committed_moves) == len([a for a in attempts if a.converged])
    np.testing.assert_allclose(U[UC::DPN], el.beta * U[DC::DPN], atol=1e-9)


def test_final_inner_newton_update_is_assessed_not_discarded():
    nodes, elems, dirichlet = strip(8)
    source = poisson_element(lambda xy: F_SOURCE)
    U = solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1,
                              max_field_newton=1, tol=1e-10)
    exact = 0.5 * F_SOURCE * nodes[:, 0] * (1.0 - nodes[:, 0])
    assert np.max(np.abs(U - exact)) < 2.0e-3
    with pytest.raises(StaggeredSolveFailed, match="inner Newton exhausted"):
        solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1,
                              max_stagger=4, max_field_newton=0, tol=1e-12)


def test_staggered_failure_reports_per_component_residual():
    nodes, elems, dirichlet = strip(8)
    with pytest.raises(StaggeredSolveFailed, match=r"per-component \|R\|.*dominant comp 0"):
        solve_steps_staggered(nodes, elems, 1, poisson_element(lambda xy: F_SOURCE),
                              lambda frac: dirichlet, [[0]], 1, max_stagger=1,
                              max_field_newton=0, tol=1e-12)


def _two_field_poisson(scale1=1.0):
    """Uncoupled pair of Quad4 Poissons; field 1 residual/tangent times scale1."""
    e = poisson_element(lambda xy: F_SOURCE)

    def element_fn(ei, coords, U_e, DU_e=None):
        U_e = np.asarray(U_e, dtype=float)
        R0, K0 = e(ei, coords, U_e[0::2])
        R1, K1 = e(ei, coords, U_e[1::2])
        R, K = np.empty(8), np.zeros((8, 8))
        R[0::2], R[1::2] = R0, scale1 * R1
        K[0::2, 0::2], K[1::2, 1::2] = K0, scale1 * K1
        return R, K

    return element_fn


def _strip2(nx):
    nodes, elems, d1 = strip(nx)
    return nodes, elems, {2 * n + c: 0.0 for n in d1 for c in (0, 1)}


def test_scaled_residual_gate_enforces_a_weak_field():
    nodes, elems, dirichlet = _strip2(16)
    U = solve_steps_staggered(nodes, elems, 2, _two_field_poisson(1.0e-10),
                              lambda frac: dirichlet, [[0], [1]], 1, tol=1e-8,
                              residual_scale={1: 1.0e10})
    exact = 0.5 * F_SOURCE * nodes[:, 0] * (1.0 - nodes[:, 0])
    assert np.max(np.abs(U[0::2] - exact)) < 2.0e-3
    assert np.max(np.abs(U[1::2] - exact)) < 2.0e-3


def test_unscaled_weak_field_is_not_solved_broken_control():
    nodes, elems, dirichlet = _strip2(16)
    U = solve_steps_staggered(nodes, elems, 2, _two_field_poisson(1.0e-10),
                              lambda frac: dirichlet, [[0], [1]], 1, tol=1e-8)
    assert np.max(np.abs(U[1::2])) < 1.0e-8


def test_update_gate_requires_a_second_sweep():
    nodes, elems, dirichlet = strip(8)
    source = poisson_element(lambda xy: F_SOURCE)
    with pytest.raises(StaggeredSolveFailed, match=r"scaled \|dU\|"):
        solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1,
                              max_stagger=1, max_field_newton=8, tol=1e-10)
    U = solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1,
                              max_stagger=2, max_field_newton=8, tol=1e-10)
    assert U.max() > 0.5


def test_update_gate_passes_for_stationary_fields():
    """A field that does not move (a hold) passes the update gate at once."""
    nodes, elems, dirichlet = strip(8)
    source = poisson_element(lambda xy: F_SOURCE)
    U1 = solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1,
                               tol=1e-10)
    U, report = solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet,
                                      [[0]], 1, tol=1e-10, U0=U1, return_report=True)
    assert report.attempts[0].sweeps == 1 and report.attempts[0].scaled_update == 0.0


def test_zero_update_is_not_convergence():
    nodes, elems, dirichlet = strip(8)
    with pytest.raises(StaggeredSolveFailed, match="zero update with residual remaining"):
        solve_steps_staggered(nodes, elems, 1, poisson_element(lambda xy: F_SOURCE),
                              lambda frac: dirichlet, [[0]], 1, max_field_newton=4,
                              tol=1e-10, linear_solve=lambda K, b: np.zeros(K.shape[0]))


def test_line_search_exhaustion_is_a_hard_failure():
    from scipy.sparse.linalg import spsolve
    nodes, elems, dirichlet = strip(8)
    with pytest.raises(StaggeredSolveFailed, match="line-search exhausted"):
        solve_steps_staggered(nodes, elems, 1, poisson_element(lambda xy: F_SOURCE),
                              lambda frac: dirichlet, [[0]], 1, max_field_newton=4,
                              tol=1e-10, linear_solve=lambda K, b: -spsolve(K, b))


def test_callback_residual_matches_a_fresh_assembly():
    nodes, elems, dirichlet = strip(8)
    source = poisson_element(lambda xy: F_SOURCE)
    seen = []
    U = solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1,
                              tol=1e-10, step_callback=lambda s, f, u, r: seen.append(dict(r)))
    R, _ = assemble(nodes, elems, 1, source, U, du_is_total=True)
    R[np.array(sorted(dirichlet), dtype=int)] = 0.0
    assert abs(seen[-1]["residual_norm"] - np.linalg.norm(R)) < 1e-14
    assert seen[-1]["scaled_residual_norms"][0] < 1e-10


def test_bounded_monolithic_reports_fresh_residual_and_actual_last_update():
    """Bounded Newton accepts by the KKT residual; it reports the measured,
    nonzero last update rather than a fabricated zero."""
    nodes, elems, dirichlet = strip(8)
    source = poisson_element(lambda xy: F_SOURCE)
    seen = []
    U, report = solve_steps_staggered(
        nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1, tol=1e-10,
        update_tol=1e-30, newton_method="bounded_monolithic", field_scale={0: 2.0},
        return_report=True, step_callback=lambda s, f, u, r: seen.append(r))
    R, _ = assemble(nodes, elems, 1, source, U, du_is_total=True)
    R[np.array(sorted(dirichlet), dtype=int)] = 0.0
    rec = report.attempts[-1]
    assert rec.converged and rec.sweeps == 0 and rec.newton_iterations == 2
    assert rec.newton_method == "bounded_monolithic"
    assert rec.scaled_update == pytest.approx(np.linalg.norm(U / 2.0)) and rec.scaled_update > .1
    assert rec.raw_residual == pytest.approx(np.linalg.norm(R), abs=1e-15)
    assert seen[-1]["residual_norm"] == rec.raw_residual
    assert seen[-1]["scaled_update_norms"] == [rec.scaled_update]
    assert 0.0 <= rec.scaled_residual < 1e-10


def test_bounded_monolithic_converges_on_its_final_update():
    """Defect regression through the driver: one permitted coupled update of
    a linear problem is converged, not a cutback."""
    nodes, elems, dirichlet = strip(8)
    U, report = solve_steps_staggered(
        nodes, elems, 1, poisson_element(lambda xy: F_SOURCE), lambda frac: dirichlet,
        [[0]], 1, tol=1e-10, max_field_newton=1, newton_method="bounded_monolithic",
        return_report=True)
    assert report.cutbacks == 0 and report.attempts[-1].newton_iterations == 2


@pytest.mark.parametrize("bad_update", ["zero", "nonfinite"])
def test_bounded_monolithic_rejects_without_finite_kkt_certificate(bad_update):
    nodes, elems, dirichlet = strip(4)
    commits = []
    fill = 0.0 if bad_update == "zero" else np.nan
    with pytest.raises(StaggeredSolveFailed) as failed:
        solve_steps_staggered(
            nodes, elems, 1, poisson_element(lambda xy: F_SOURCE), lambda frac: dirichlet,
            [[0]], 1, tol=1e-10, newton_method="bounded_monolithic",
            linear_solve=lambda K, b: np.full(K.shape[0], fill),
            commit_fn=lambda *args: commits.append(True), return_report=True)
    assert not commits
    rec = failed.value.report.attempts[-1]
    assert not rec.converged and not np.isfinite(rec.scaled_residual)
    assert ("line search exhausted" if bad_update == "zero" else "non-finite") in rec.reason


def test_certified_nonactive_clip_reports_the_fresh_residual():
    nodes, elems, dirichlet = strip(12)
    source = poisson_element(lambda xy: F_SOURCE)
    seen = []
    U, report = solve_steps_staggered(
        nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 1, tol=1e-10,
        field_bounds={0: (-np.inf, 2.0)},
        step_callback=lambda s, f, u, r: seen.append(dict(r)), return_report=True)
    R, _ = assemble(nodes, elems, 1, source, U, du_is_total=True)
    R[np.array(sorted(dirichlet), dtype=int)] = 0.0
    assert report.cutbacks == 0
    assert seen[-1]["scaled_residual_norms"][0] < 1e-10
    assert abs(seen[-1]["residual_norm"] - np.linalg.norm(R)) < 1e-14


class _DtGatedPoisson(_Clock):
    """Poisson that mutates SVARS and rejects increments above ``max_dt``."""

    def __init__(self, max_dt, fail_if_moved=False, always_fail=False):
        self.inner = poisson_element(lambda xy: F_SOURCE)
        self.max_dt, self.fail_if_moved, self.always_fail = max_dt, fail_if_moved, always_fail
        self.svars = np.zeros(1)
        self._init_clock()
        self.seen_unorm = []

    def __call__(self, ei, coords, U_e, DU_e=None):
        self.svars[0] += 1.0
        self.seen_unorm.append(float(np.linalg.norm(U_e)))
        if self.always_fail:
            raise IncrementTooLarge("always")
        moved = np.linalg.norm(U_e) > 1e-12
        if self.dt > self.max_dt + 1e-12 and (moved if self.fail_if_moved else True):
            raise IncrementTooLarge(f"dt={self.dt}")
        return self.inner(ei, coords, U_e, DU_e)


def _history_poisson():
    """Path-dependent Poisson: R += DU, so subincrement history matters."""
    inner = poisson_element(lambda xy: F_SOURCE)

    def element_fn(ei, coords, U_e, DU_e=None):
        R, K = inner(ei, coords, U_e)
        du = np.zeros_like(U_e) if DU_e is None else np.asarray(DU_e, dtype=float)
        return R + du, K + np.eye(len(U_e))

    return element_fn


class _DtGatedHistory(_DtGatedPoisson):
    def __init__(self, max_dt):
        super().__init__(max_dt)
        self.inner = _history_poisson()


def test_bounded_line_search_cutback_restores_transaction_and_completes():
    """A rejected trial is rolled back before retrying the physical half steps."""
    nodes, elems, dirichlet = strip(4)

    class TrialGatedHistory(_Clock):
        def __init__(self, reject=True):
            self.inner = _history_poisson()
            self.reject = reject
            self.rejections = 0
            self.svars = np.array([0.25])
            self.time, self.dt, self.kinc = 2.0, 0.125, 11
            self.commits = []

        def __call__(self, ei, coords, U_e, DU_e):
            if self.reject and self.dt > 0.6 and np.linalg.norm(DU_e) > 1e-12:
                self.rejections += 1
                # Corrupt every transactional field; only a nonzero line-search
                # trial rejects, the initial assembly succeeds.
                self.svars[:] = -123.0
                self.time, self.dt, self.kinc = -99.0, -1.0, -2
                raise IncrementTooLarge("material trial requires dt <= 0.6")
            return self.inner(ei, coords, U_e, DU_e)

        def commit(self, nodes_, elems_, dpn_, U, U_prev):
            self.svars[0] += np.linalg.norm(U - U_prev)
            self.commits.append((self.dt, self.time, self.kinc, U.copy()))

    live = TrialGatedHistory()
    rejected, accepted = [], []

    def on_attempt(rec):
        if rec.converged:
            return
        rejected.append(rec)
        np.testing.assert_array_equal(live.svars, [0.25])
        assert (live.time, live.dt, live.kinc) == (2.0, 0.125, 11)
        assert not live.commits

    def run(el, n_steps, dt, **kw):
        return solve_steps_staggered(
            nodes, elems, 1, el, lambda frac: dirichlet, [[0]], n_steps,
            newton_method="bounded_monolithic", cutback_on_element_rejection=True,
            tol=1e-10, box=BoxSpec({0: (0.0, 2.0)}, irreversible=(0,)),
            state=CompiledElementState([el]), commit_fn=el.commit,
            dt=dt, max_cutbacks=1, return_report=True, **kw)

    U, report = run(live, 1, 1.0, attempt_callback=on_attempt,
                    accepted_callback=accepted.append)
    reference = TrialGatedHistory(reject=False)
    U_ref, ref_report = run(reference, 2, 0.5)
    assert live.rejections == len(rejected) == report.cutbacks == 1
    assert rejected[0].newton_iterations == 1
    assert "during line search" in rejected[0].reason
    assert "material trial requires dt <= 0.6" in rejected[0].reason
    assert [a["frac_done"] for a in accepted] == [0.5, 1.0]
    assert [a["time_end"] for a in accepted] == [0.5, 1.0]
    assert [a["next_kinc"] for a in accepted] == [2, 3]
    assert [c[:3] for c in live.commits] == [(0.5, 0.0, 1), (0.5, 0.5, 2)]
    assert all(a.converged and a.scaled_residual < 1e-10 for a in report.attempts[1:])
    assert ref_report.cutbacks == 0
    np.testing.assert_array_equal(U, U_ref)
    np.testing.assert_array_equal(live.svars, reference.svars)
    assert (live.time, live.dt, live.kinc) == (reference.time, reference.dt, reference.kinc)
    assert U.max() > 0.1


def test_rejected_attempt_restores_u_svars_and_clock():
    nodes, elems, dirichlet = strip(8)
    el = _DtGatedPoisson(max_dt=0.6, fail_if_moved=True)
    el.time, el.kinc = 3.0, 11
    state = CompiledElementState([el])
    svars0, clock0 = el.svars.copy(), (el.time, el.dt, el.kinc)
    rejected = []

    def on_attempt(rec):
        if not rec.converged:
            rejected.append(True)
            assert np.array_equal(el.svars, svars0)
            assert (el.time, el.dt, el.kinc) == clock0

    U, report = solve_steps_staggered(
        nodes, elems, 1, el, lambda frac: dirichlet, [[0]], 1, tol=1e-9, state=state,
        dt=1.0, max_cutbacks=2, attempt_callback=on_attempt, return_report=True)
    assert rejected and report.cutbacks >= 1 and U.max() > 0.5
    nz = next(i for i, v in enumerate(el.seen_unorm) if v > 1e-8)
    assert any(v < 1e-12 for v in el.seen_unorm[nz + 1:]), (
        "next trial did not restart from the accepted U")


def test_no_scientific_output_for_a_rejected_attempt():
    nodes, elems, dirichlet = strip(8)
    el = _DtGatedPoisson(max_dt=0.6)
    records, commits, attempts = [], [], []
    solve_steps_staggered(
        nodes, elems, 1, el, lambda frac: dirichlet, [[0]], 1, tol=1e-9,
        state=CompiledElementState([el]), dt=1.0, max_cutbacks=3,
        step_callback=lambda s, f, U, r: records.append(f),
        commit_fn=lambda *a: commits.append(True), attempt_callback=attempts.append)
    rejected = [a for a in attempts if not a.converged]
    accepted = [a for a in attempts if a.converged]
    assert rejected and (rejected[0].frac_start, rejected[0].frac_end) == (0.0, 1.0)
    assert [(a.frac_start, a.frac_end) for a in accepted] == [(0.0, 0.5), (0.5, 1.0)]
    assert len(commits) == len(accepted)
    assert records == [1.0]


def test_cutback_budget_and_nominal_output_times():
    nodes, elems, dirichlet = strip(6)
    el = _DtGatedPoisson(max_dt=0.6, always_fail=True)
    with pytest.raises(StaggeredSolveFailed, match="element requested a smaller increment"):
        solve_steps_staggered(nodes, elems, 1, el, lambda frac: dirichlet, [[0]], 1,
                              tol=1e-9, state=CompiledElementState([el]), dt=1.0,
                              max_cutbacks=2)
    fracs = []
    ok = _DtGatedPoisson(max_dt=0.6)
    solve_steps_staggered(nodes, elems, 1, ok, lambda frac: dirichlet, [[0]], 2, tol=1e-9,
                          state=CompiledElementState([ok]), dt=1.0, max_cutbacks=3,
                          step_callback=lambda s, f, U, r: fracs.append(f))
    assert fracs == [0.5, 1.0]


def test_forced_half_steps_match_the_two_step_reference():
    nodes, elems, dirichlet = strip(12)
    U_ref = solve_steps_staggered(nodes, elems, 1, _history_poisson(),
                                  lambda frac: dirichlet, [[0]], 2, tol=1e-10, dt=0.5)
    el = _DtGatedHistory(max_dt=0.5)
    U_ad = solve_steps_staggered(nodes, elems, 1, el, lambda frac: dirichlet, [[0]], 1,
                                 tol=1e-10, state=CompiledElementState([el]), dt=1.0,
                                 max_cutbacks=3)
    assert np.max(np.abs(U_ad - U_ref)) < 1e-8


def test_early_stop_ends_the_loop_and_is_recorded():
    nodes, elems, dirichlet = strip(8)
    source = poisson_element(lambda xy: F_SOURCE)
    fracs = []
    U, report = solve_steps_staggered(
        nodes, elems, 1, source, lambda frac: dirichlet, [[0]], 4, tol=1e-10,
        step_callback=lambda s, f, u, r: fracs.append(f),
        early_stop=lambda U, frac, state: frac >= 0.5, return_report=True)
    assert fracs == [0.25, 0.5]
    assert report.stopped_early and report.steps_accepted == 2
    _, full = solve_steps_staggered(nodes, elems, 1, source, lambda frac: dirichlet,
                                    [[0]], 4, tol=1e-10, return_report=True)
    assert not full.stopped_early and full.steps_accepted == 4


# --------------------------------------------------------------------------
# 3. feasible trials
# --------------------------------------------------------------------------

U_F, D_F, DPN_F = 0, 1, 2


class RecordingToy:
    """Per node: R_u = u - 1 + 0.1 d and R_d = d**2 - 4 + 0.1 u.

    The unconstrained root has d near 2, outside the box [0, 1]; with d held
    at its upper bound the solution is d = 1, u = 0.9.
    """

    def __init__(self):
        self.max_d, self.min_d, self.calls = -np.inf, np.inf, 0

    def __call__(self, ei, coords, U_e, DU_e):
        U_e = np.asarray(U_e, dtype=float)
        u, d = U_e[U_F::DPN_F], U_e[D_F::DPN_F]
        self.calls += 1
        self.max_d = max(self.max_d, float(d.max()))
        self.min_d = min(self.min_d, float(d.min()))
        R, K = np.zeros_like(U_e), np.zeros((U_e.size, U_e.size))
        for a in range(len(u)):
            iu, idd = DPN_F * a + U_F, DPN_F * a + D_F
            R[iu] = u[a] - 1.0 + 0.1 * d[a]
            R[idd] = d[a] ** 2 - 4.0 + 0.1 * u[a]
            K[iu, iu], K[iu, idd] = 1.0, 0.1
            K[idd, idd], K[idd, iu] = 2.0 * d[a], 0.1
        return R, K


def _feasible_solve(element, **kw):
    nodes, elems = structured_quad_mesh(2, 2)
    U0 = np.zeros(len(nodes) * DPN_F)
    U0[D_F::DPN_F] = 0.5
    return solve_steps_staggered(
        nodes, elems, DPN_F, element, lambda frac: {}, [[U_F], [D_F]], 1,
        box=BoxSpec(bounds={D_F: (0.0, 1.0)}), U0=U0, tol=1e-10, return_report=True, **kw)


def test_feasible_trials_never_evaluates_outside_the_box():
    element = RecordingToy()
    U, report = _feasible_solve(element, feasible_trials=True)
    assert element.calls > 0 and element.max_d <= 1.0 and element.min_d >= 0.0
    np.testing.assert_allclose(U[D_F::DPN_F], 1.0, atol=0.0)
    np.testing.assert_allclose(U[U_F::DPN_F], 0.9, atol=1e-10)
    assert report.attempts[-1].converged


def test_default_staggered_solve_still_evaluates_outside_the_box():
    """Broken control: without feasible trials the iterates leave the box."""
    element = RecordingToy()
    U, _ = _feasible_solve(element, max_projection_rounds=2)
    assert element.max_d > 1.5
    np.testing.assert_allclose(U[D_F::DPN_F], 1.0, atol=0.0)
    np.testing.assert_allclose(U[U_F::DPN_F], 0.9, atol=1e-10)


def test_feasible_trials_needs_bounds():
    nodes, elems = structured_quad_mesh(2, 2)
    with pytest.raises(ValueError, match="feasible_trials"):
        solve_steps_staggered(nodes, elems, DPN_F, RecordingToy(), lambda frac: {},
                              [[U_F], [D_F]], 1, feasible_trials=True)


def _linear(ei, coords, U_e, DU_e):
    U_e = np.asarray(U_e, dtype=float)
    return U_e - 2.0, np.eye(U_e.size)


def test_incompatible_prescribed_value_is_rejected_with_static_bounds():
    nodes, elems = structured_quad_mesh(1, 1)
    with pytest.raises(ValueError, match="outside the box"):
        solve_steps_staggered(nodes, elems, 1, _linear, lambda frac: {0: 2.0}, [[0]], 1,
                              field_bounds={0: (0.0, 1.0)}, feasible_trials=True)


def test_compatible_prescribed_value_is_preserved_with_static_bounds():
    nodes, elems = structured_quad_mesh(1, 1)
    U = solve_steps_staggered(nodes, elems, 1, _linear, lambda frac: {0: 0.5}, [[0]], 1,
                              field_bounds={0: (0.0, 1.0)}, feasible_trials=True, tol=1e-10)
    assert U[0] == 0.5
    np.testing.assert_allclose(U[1:], 1.0)


class _OvershootToy:
    """Per node (u, L, d): R_u = u + u^3 - 2 L, R_L = L, R_d = d - 0.5.

    ``L`` is prescribed (the load).  At the full load the first Newton trial
    of the u block from 0 is u = 2, which the element rejects (u > 1.5);
    half that step is the exact root u = 1.
    """

    def __init__(self):
        self.rejections = 0

    def __call__(self, ei, coords, U_e, DU_e):
        U_e = np.asarray(U_e, dtype=float)
        u, L, d = U_e[0::3], U_e[1::3], U_e[2::3]
        if u.max() > 1.5:
            self.rejections += 1
            raise IncrementTooLarge("u above the element's admissible range")
        R, K = np.zeros(12), np.zeros((12, 12))
        for a in range(4):
            R[3 * a] = u[a] + u[a] ** 3 - 2.0 * L[a]
            K[3 * a, 3 * a], K[3 * a, 3 * a + 1] = 1.0 + 3.0 * u[a] ** 2, -2.0
            R[3 * a + 1], K[3 * a + 1, 3 * a + 1] = L[a], 1.0
            R[3 * a + 2], K[3 * a + 2, 3 * a + 2] = d[a] - 0.5, 1.0
        return R, K


@pytest.mark.parametrize("flag", [False, True])
def test_feasible_trials_forward_cutback_on_element_rejection(flag):
    """Defect regression: with both options the feasible-trial blocks ignored
    the flag, backtracked past the rejected trial and recorded no cutback."""
    nodes, elems = structured_quad_mesh(1, 1)
    element = _OvershootToy()
    U, report = solve_steps_staggered(
        nodes, elems, 3, element, lambda frac: {3 * n + 1: frac for n in range(4)},
        [[0, 1], [2]], 1, tol=1e-10, box=BoxSpec(bounds={2: (0.0, 1.0)}),
        feasible_trials=True, cutback_on_element_rejection=flag, max_cutbacks=2,
        return_report=True)
    np.testing.assert_allclose(U[0::3], 1.0, atol=1e-10)
    assert element.rejections >= 1
    if flag:
        assert report.cutbacks == 1
        assert "during line search" in report.attempts[0].reason
    else:
        assert report.cutbacks == 0 and len(report.attempts) == 1


# --------------------------------------------------------------------------
# 4. TIME / DTIME / KINC on the public adapter, cutbacks and restart
# --------------------------------------------------------------------------

def _clock_run(method="staggered", restart=None, accepted=None, cache=False, n_steps=2):
    """Two nominal steps of dt = 1; the kernel rejects DTIME > 0.6 (PNEWDT)."""
    nodes, elems, dirichlet = strip(4)
    element, kernel = recording_element(len(elems), dt=1.0, mass=1.0, cache_batches=cache)
    kernel.reject = lambda time, dtime, u, du: dtime > 0.6
    kw = {} if method == "staggered" else {"newton_method": method}
    U, report = solve_steps_staggered(
        nodes, elems, 1, element.element_rk, lambda frac: dirichlet, [[0]], n_steps,
        batch_fn=element.element_rk_batch, commit_fn=element.commit_step,
        state=CompiledElementState([element]), dt=1.0, max_cutbacks=3, tol=1e-10,
        return_report=True, restart=restart, accepted_callback=accepted, **kw)
    return U, report, element, kernel


def _distinct(calls):
    out = []
    for c in calls:
        key = (c[0], c[2], c[4])
        if not out or out[-1] != key:
            out.append(key)
    return out


@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_public_adapter_sees_beginning_time_reduced_dtime_and_kinc(method):
    U, report, element, kernel = _clock_run(method)
    assert all(c[0] == c[1] and c[3] == 1 for c in kernel.calls)   # TIME(1)=TIME(2), KSTEP=1
    # (TIME at the start of the attempt, DTIME, KINC); the rejected full step
    # at TIME 0 is retried as two halves with KINC 1 then 2.
    assert _distinct(kernel.calls) == [(0.0, 1.0, 1), (0.0, 0.5, 1), (0.5, 0.5, 2),
                                       (1.0, 0.5, 3), (1.5, 0.5, 4)]
    assert [(a.frac_start, a.frac_end, a.dt, a.converged) for a in report.attempts] == [
        (0.0, 0.5, 1.0, False), (0.0, 0.25, 0.5, True), (0.25, 0.5, 0.5, True),
        (0.5, 0.75, 0.5, True), (0.75, 1.0, 0.5, True)]
    assert report.steps_accepted == 2 and report.cutbacks == 1
    np.testing.assert_allclose(element.svars, 2.0)        # four committed DTIMEs of 0.5
    assert (element.time, element.dt, element.kinc) == (1.5, 0.5, 4)


def test_public_adapter_batch_cache_gives_identical_results():
    U0, _, el0, k0 = _clock_run(cache=False)
    U1, _, el1, k1 = _clock_run(cache=True)
    np.testing.assert_array_equal(U0, U1)
    np.testing.assert_array_equal(el0.svars, el1.svars)
    assert el1.batch_cache_hits >= 4                      # one per accepted commit
    assert len(k1.calls) == len(k0.calls) - el1.batch_cache_hits
    assert el1._batch_cache is None


@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
@pytest.mark.parametrize("legacy", [False, True])
def test_restart_after_a_rejected_subincrement_reproduces_the_run(method, legacy):
    """Resume from the first accepted subincrement after the rejection: the
    resumed solve passes the same TIME/DTIME/KINC and ends identically."""
    saved = []
    U_full, report_full, el_full, k_full = _clock_run(method, accepted=saved.append)
    point = saved[0]
    assert (point["frac_done"], point["time_end"], point["next_kinc"],
            point["cutbacks_this_step"]) == (0.25, 0.5, 2, 1)
    if legacy:   # a dictionary saved without the newer fields
        point = {k: v for k, v in point.items() if k != "next_kinc"}
        point["report"] = dict(point["report"], attempts=[
            {k: v for k, v in a.items()
             if k not in ("newton_method", "newton_iterations", "projection_rounds")}
            for a in point["report"]["attempts"]])
    U, report, element, kernel = _clock_run(method, restart=point)
    n_before = next(i for i, c in enumerate(k_full.calls) if c[4] == 2)
    assert kernel.calls == k_full.calls[n_before:]
    np.testing.assert_array_equal(U, U_full)
    np.testing.assert_array_equal(element.svars, el_full.svars)
    assert (element.time, element.dt, element.kinc) == (el_full.time, el_full.dt,
                                                        el_full.kinc)
    if legacy:
        assert all(a.newton_method == "staggered" for a in report.attempts[:2])
        assert [asdict(a) for a in report.attempts[2:]] == [
            asdict(a) for a in report_full.attempts[2:]]
    else:
        assert report == report_full


def test_restart_reconstructs_the_irreversible_bound():
    """A resumed run with the dynamic box is indistinguishable from the
    uninterrupted one and still never heals."""
    nodes, elems = structured_quad_mesh(2, 2)
    box = BoxSpec(bounds={DC: (0.0, 1.0)}, irreversible=(DC,))

    def run(el, rows, **kw):
        return solve_steps_staggered(
            nodes, elems, DPN, el.element_rk, _heal_dirichlet(nodes),
            [[UC, TC, SC], [DC]], 2, tol=1e-11, box=box, max_projection_rounds=4,
            state=CompiledElementState([el]), dt=1.0, max_cutbacks=3, commit_fn=el.commit,
            step_callback=lambda s, f, U, r: rows.append(
                {"step": int(s), "frac": float(f), "d_min": float(U[DC::DPN].min())}),
            return_report=True, **kw)

    full_el, full_rows = HealingDamageElement(nodes, elems, max_dt=0.6), []
    U_full, report_full = run(full_el, full_rows)

    class StopHere(Exception):
        pass

    saved = {}

    def save(restart_state):
        saved["state"] = restart_state
        saved["rows"] = list(live_rows)
        if restart_state["report"]["steps_accepted"] == 1:
            raise StopHere

    live_rows = []
    with pytest.raises(StopHere):
        run(HealingDamageElement(nodes, elems, max_dt=0.6), live_rows,
            accepted_callback=save)
    rows = saved["rows"]
    assert [r["step"] for r in rows] == [1]
    resumed_el = HealingDamageElement(nodes, elems, max_dt=0.6)
    U_resumed, report_resumed = run(resumed_el, rows, restart=saved["state"])
    np.testing.assert_array_equal(U_resumed, U_full)
    np.testing.assert_array_equal(resumed_el.svars, full_el.svars)
    assert (resumed_el.time, resumed_el.dt, resumed_el.kinc) == (
        full_el.time, full_el.dt, full_el.kinc)
    assert rows == full_rows and report_resumed == report_full
    assert U_resumed[DC::DPN].min() >= rows[0]["d_min"] - 1e-12


# --------------------------------------------------------------------------
# 5. nonfinite output and unrelated exceptions
# --------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["staggered", "bounded_monolithic"])
def test_nonfinite_output_at_an_iterate_is_a_cutback(method):
    nodes, elems, dirichlet = strip(4)
    element, kernel = recording_element(len(elems), dt=1.0, mass=1.0)
    kernel.nonfinite = lambda time, dtime, u, du: dtime > 0.6
    kw = {} if method == "staggered" else {"newton_method": method}
    U, report = solve_steps_staggered(
        nodes, elems, 1, element.element_rk, lambda frac: dirichlet, [[0]], 1,
        batch_fn=element.element_rk_batch, commit_fn=element.commit_step,
        state=CompiledElementState([element]), dt=1.0, max_cutbacks=2, tol=1e-10,
        return_report=True, **kw)
    assert report.cutbacks == 1
    assert "non-finite" in report.attempts[0].reason
    assert all(a.converged for a in report.attempts[1:])


def test_nonfinite_line_search_trial_backtracks_in_a_staggered_block():
    nodes, elems, dirichlet = strip(4)
    element, kernel = recording_element(len(elems))
    reference = solve_steps_staggered(nodes, elems, 1, element.element_rk,
                                      lambda frac: dirichlet, [[0]], 1,
                                      batch_fn=element.element_rk_batch, tol=1e-10)
    kernel.calls.clear()
    kernel.nonfinite = lambda t, dt, u, du: len(kernel.calls) == 2   # first trial only
    U, report = solve_steps_staggered(nodes, elems, 1, element.element_rk,
                                      lambda frac: dirichlet, [[0]], 1,
                                      batch_fn=element.element_rk_batch, tol=1e-10,
                                      return_report=True)
    np.testing.assert_allclose(U, reference, atol=1e-10)
    assert report.cutbacks == 0


@pytest.mark.parametrize("method, feasible", [
    ("staggered", False), ("staggered", True), ("bounded_monolithic", False)])
def test_unrelated_element_exceptions_propagate(method, feasible):
    nodes, elems, dirichlet = strip(4)
    source = poisson_element(lambda xy: F_SOURCE)

    def element_fn(ei, coords, U_e, DU_e):
        if np.max(U_e) > 0.0:
            raise RuntimeError("element bug")
        return source(ei, coords, U_e)

    kw = {"field_bounds": {0: (-10.0, 10.0)}, "feasible_trials": True} if feasible else {}
    with pytest.raises(RuntimeError, match="element bug"):
        solve_steps_staggered(nodes, elems, 1, element_fn, lambda frac: dirichlet, [[0]], 1,
                              tol=1e-10, max_cutbacks=2, newton_method=method, **kw)


# --------------------------------------------------------------------------
# 6. compiled UEL through the cutback driver
# --------------------------------------------------------------------------

def test_compiled_uel_sees_beginning_time_and_reduced_dtime_on_cutback(tmp_path):
    """A generated UEL with history ``a_new = a_old + DTIME * tr(F - I)``.

    The full step is rejected (PNEWDT, from the host here) and retried as two
    halves; the result must equal a uniform two-step run bit for bit.  SVARS
    are seeded with 5: the first increment starts at TIME(2) = 0, so the
    generated initialization branch ignores the seed, as in Abaqus.
    """
    import abaqus_ufl as au
    from abaqus_ufl.core.tensor import eye, trace
    from abaqus_ufl.fe import CompiledAbaqusElement, build_compiled_uel

    class Accumulator(au.Material):
        props = dict(G=2.0)
        state_vars = dict(a=0.0)

        def stress_PK1(self, F, a_old, dt):
            strain = F - eye(3)
            return self.G * strain, {"a": a_old + dt * trace(strain)}

    class AccumulatorProblem(au.WeakForm):
        material = Accumulator
        ndim = 2

        def define_fields(self):
            self.u = au.VectorField("u", degree=1)

        def momentum_equation(self, v, F, a_old, dt):
            return self.material.stress_PK1(F, a_old, dt)

    problem = AccumulatorProblem()
    module = build_compiled_uel(problem, workdir=tmp_path)
    nodes, elems = structured_quad_mesh(2, 1)
    left, right = np.flatnonzero(nodes[:, 0] == 0.0), np.flatnonzero(nodes[:, 0] == 1.0)

    def dirichlet(frac):
        bc = {2 * int(n) + c: 0.0 for n in left for c in (0, 1)}
        bc.update({2 * int(n): 0.02 * frac for n in right})
        return bc

    def run(n_steps, dt, gate):
        element = CompiledAbaqusElement(module, problem, dt=dt, n_elem=len(elems))
        element.svars[:] = 5.0
        times = []

        def batch(coords, U, DU):
            times.append((element.time, element.dt, element.kinc))
            if gate and element.dt > 0.6:
                raise IncrementTooLarge("host limit on DTIME")
            return element.element_rk_batch(coords, U, DU)

        U, report = solve_steps_staggered(
            nodes, elems, 2, element.element_rk, dirichlet, [[0, 1]], n_steps,
            batch_fn=batch, commit_fn=element.commit_step,
            state=CompiledElementState([element]), dt=dt, max_cutbacks=2, tol=1e-10,
            return_report=True)
        return U, report, element, times

    U, report, element, times = run(1, 1.0, gate=True)
    U_ref, report_ref, element_ref, _ = run(2, 0.5, gate=False)
    assert report.cutbacks == 1 and report_ref.cutbacks == 0
    assert sorted(set(times)) == [(0.0, 0.5, 1), (0.0, 1.0, 1), (0.5, 0.5, 2)]
    np.testing.assert_array_equal(U, U_ref)
    np.testing.assert_array_equal(element.svars, element_ref.svars)
    assert np.all(element.svars < 1.0), "the seeded SVARS were read at TIME(2) = 0"
    assert np.all(element.svars > 0.0)
