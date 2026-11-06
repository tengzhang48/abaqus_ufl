"""Gates for the serial stepping drivers: solve_steps and the multi-material solve.

The elements are closed-form toys (``fe_stepping_toys``): one-node springs and
a Quad4 Poisson element.  A committed step must meet the requested ``tol``;
an exhausted line search must raise; the state produced by the final permitted
update must be assessed; the stepping assembler scatter-adds repeated
connectivity; nonfinite trials backtrack while other exceptions propagate.
"""

import numpy as np
import pytest

from abaqus_ufl.fe import (CompiledElementState, NonfiniteElementOutput,
                           assemble, assemble_groups, boundary_nodes,
                           solve_steps, solve_steps_multimaterial,
                           structured_quad_mesh)
from abaqus_ufl.fe.stepping import _make_assembler
from fe_stepping_toys import (convergence_rate, l2_error, poisson_element,
                              recording_element)

ONE_NODE = np.array([[0.0, 0.0]])
ONE_ELEM = np.array([[0]])


def _grounded_spring(k):
    """One-node element with internal residual ``k*u``."""
    def element_fn(ei, coords, U_e, DU_e):
        return np.array([k * U_e[0]]), np.array([[k]])
    return element_fn


def _cubic_spring(u):
    """One-node residual ``u^3 + u`` with its exact tangent."""
    return np.array([u ** 3 + u]), np.array([[3.0 * u ** 2 + 1.0]])


# --------------------------------------------------------------------------
# solve_steps: load, tangent, history, acceptance
# --------------------------------------------------------------------------

def test_solve_steps_poisson_converges_at_second_order():
    pi = np.pi
    exact = lambda xy: np.sin(pi * xy[0]) * np.sin(pi * xy[1])
    source = lambda xy: 2.0 * pi * pi * exact(xy)
    hs, errs = [], []
    for nx in (4, 8, 16):
        nodes, elems = structured_quad_mesh(nx, nx)
        bc = {int(n): 0.0 for n in boundary_nodes(nodes)}
        U = solve_steps(nodes, elems, 1, poisson_element(source),
                        lambda frac: bc, 1, tol=1e-10)
        hs.append(1.0 / nx)
        errs.append(l2_error(nodes, elems, U, exact))
    assert 1.8 < convergence_rate(hs, errs) < 2.2


def test_solve_steps_constant_dead_load_and_callback_report():
    """An external callback supplies the actual load; it is not ramped by frac."""
    seen_loads, reports = [], []

    def external(step, frac, U, U_prev):
        seen_loads.append(2.0)
        return np.array([2.0]), None

    U = solve_steps(ONE_NODE, ONE_ELEM, 1, _grounded_spring(4.0), lambda frac: {}, 4,
                    external_fn=external,
                    step_callback=lambda s, f, U, r: reports.append((s, f, U[0], r)))
    assert U[0] == pytest.approx(0.5)
    assert [item[2] for item in reports] == pytest.approx([0.5] * 4)
    assert set(seen_loads) == {2.0}
    assert [item[0] for item in reports] == [1, 2, 3, 4]
    assert all(item[3]["iterations"] >= 1 for item in reports)
    assert all(item[3]["residual_norm"] < 1e-12 for item in reports)


def test_solve_steps_uses_solution_dependent_external_tangent():
    """Solve ``4u - (3 + u) = 0`` in one Newton update using ``Kext=1``."""
    reports = []
    U = solve_steps(ONE_NODE, ONE_ELEM, 1, _grounded_spring(4.0), lambda frac: {}, 1,
                    max_newton=2,
                    external_fn=lambda s, f, U, Up: (np.array([3.0 + U[0]]), np.array([[1.0]])),
                    step_callback=lambda s, f, U, r: reports.append(r))
    assert U[0] == pytest.approx(1.0)
    assert reports == [{"iterations": 2, "residual_norm": 0.0}]


def test_external_load_broken_sign_control_is_detected():
    """Reversing the load gives the opposite response, so a sign error is visible."""
    U_bad = solve_steps(ONE_NODE, ONE_ELEM, 1, _grounded_spring(4.0), lambda frac: {}, 1,
                        external_fn=lambda s, f, U, Up: (np.array([-2.0]), None))
    assert U_bad[0] == pytest.approx(-0.5)


def test_history_is_measured_from_the_accepted_step_start():
    """``R = u - DU/2`` reaches different states with one step and with two."""
    def element_fn(ei, coords, U_e, DU_e):
        return np.array([U_e[0] - 1.0 + 0.5 * DU_e[0]]), np.array([[1.5]])

    one = solve_steps(ONE_NODE, ONE_ELEM, 1, element_fn, lambda f: {}, 1, tol=1e-14)
    two = solve_steps(ONE_NODE, ONE_ELEM, 1, element_fn, lambda f: {}, 2, tol=1e-14)
    assert one[0] == pytest.approx(2.0 / 3.0)
    assert two[0] == pytest.approx(8.0 / 9.0)   # 1.5 u = 1 + (2/3)/2


def test_solve_steps_rejects_a_step_that_misses_the_requested_tol():
    """No fallback acceptance: one Newton update of ``u^3 + u = 0.002`` leaves
    the residual 8e-9, above ``tol = 1e-12`` but below a looser absolute 1e-6."""
    element_fn = lambda ei, c, U_e, DU_e: _cubic_spring(U_e[0])
    external = lambda s, f, U, Up: (np.array([0.002]), None)
    with pytest.raises(RuntimeError, match=r"did not converge in 1 Newton"):
        solve_steps(ONE_NODE, ONE_ELEM, 1, element_fn, lambda frac: {}, 1,
                    max_newton=1, tol=1e-12, external_fn=external)
    reports = []
    U = solve_steps(ONE_NODE, ONE_ELEM, 1, element_fn, lambda frac: {}, 1, max_newton=5,
                    tol=1e-12, external_fn=external,
                    step_callback=lambda s, f, U, r: reports.append(r))
    assert abs(U[0] ** 3 + U[0] - 0.002) < 1e-12
    assert reports[0]["residual_norm"] < 1e-12


def test_solve_steps_rejects_an_exhausted_line_search():
    """A non-descent direction raises after one assembly and 30 trials."""
    calls = []

    def element_fn(ei, coords, U_e, DU_e):
        calls.append(U_e[0])
        return np.array([4.0 * U_e[0] - 1.0]), np.array([[-4.0]])

    with pytest.raises(RuntimeError, match="line search exhausted"):
        solve_steps(ONE_NODE, ONE_ELEM, 1, element_fn, lambda frac: {}, 1, max_newton=40)
    assert len(calls) == 31


def test_solve_steps_accepts_a_step_converged_by_its_last_update():
    reports = []
    U = solve_steps(ONE_NODE, ONE_ELEM, 1, _grounded_spring(4.0), lambda frac: {}, 1,
                    max_newton=1, tol=1e-10,
                    external_fn=lambda s, f, U, Up: (np.array([2.0]), None),
                    step_callback=lambda s, f, U, r: reports.append(r))
    assert U[0] == pytest.approx(0.5)
    assert reports == [{"iterations": 1, "residual_norm": 0.0}]


def test_solve_steps_state_hook_sets_kinc_and_beginning_time_per_step():
    calls = []

    class StubState:
        def begin_increment(self, kinc, time_begin=None, dt=None):
            calls.append((kinc, time_begin, dt))

    solve_steps(np.array([[0.0, 0.0], [1.0, 0.0]]), [[0, 1]], 1,
                lambda ei, c, U_e, DU_e: (U_e.copy(), np.eye(2)),
                lambda frac: {0: 0.0}, 3, state=StubState())
    assert calls == [(1, None, None), (2, None, None), (3, None, None)]

    # With the public adapter, TIME is the beginning of each step and KINC the step.
    nodes, elems = structured_quad_mesh(2, 1)
    element, kernel = recording_element(len(elems), dt=0.25)
    bc = {0: 0.0, 3: 0.0, 2: 0.0, 5: 0.0}
    solve_steps(nodes, elems, 1, element.element_rk, lambda frac: bc, 3,
                batch_fn=element.element_rk_batch, commit_fn=element.commit_step,
                state=CompiledElementState([element]), tol=1e-10)
    assert {(c[0], c[2], c[4]) for c in kernel.calls} == {
        (0.0, 0.25, 1), (0.25, 0.25, 2), (0.5, 0.25, 3)}
    np.testing.assert_allclose(element.svars, 0.75)   # committed DTIME, three steps


# --------------------------------------------------------------------------
# the stepping assembler
# --------------------------------------------------------------------------

def test_stepping_assembler_accumulates_repeated_connectivity():
    """A collapsed Quad4 (repeated node) scatter-adds both local slots; the
    public ``assemble`` stays strict and rejects that connectivity."""
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
    elems = np.array([[0, 1, 2, 0]])
    element_fn = lambda ei, c, U_e, DU_e: (np.array([1.0, 2.0, 3.0, 4.0]), np.eye(4))
    batch_fn = lambda c, U, DU: (np.array([[1.0, 2.0, 3.0, 4.0]]), np.eye(4)[None])
    for fn, batch in ((element_fn, None), (None, batch_fn)):
        R, K = _make_assembler(nodes, elems, 1, fn, batch)(np.zeros(3), np.zeros(3))
        np.testing.assert_allclose(R, [5.0, 2.0, 3.0])
        assert K[0, 0] == 2.0
    with pytest.raises(ValueError, match="distinct"):
        assemble(nodes, elems, 1, element_fn, np.zeros(3), du_is_total=True)


def test_collapsed_element_solves_through_solve_steps():
    """A polar-centre style collapsed Quad4 with a linear element converges."""
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
    elems = np.array([[0, 1, 2, 0]])
    K = np.array([[2.0, -1, 0, -1], [-1, 2, -1, 0], [0, -1, 2, -1], [-1, 0, -1, 2]])
    element_fn = lambda ei, c, U_e, DU_e: (K @ U_e, K)
    U = solve_steps(nodes, elems, 1, element_fn, lambda f: {1: 1.0, 2: 1.0}, 1, tol=1e-12)
    # Node 0 occupies slots 0 and 3: rows 0 + 3 give 2 u0 - u1 - u2 = 0.
    assert U[0] == pytest.approx(1.0)


@pytest.mark.parametrize("bad", ["residual_shape", "tangent_shape"])
def test_stepping_assembler_validates_output_shapes(bad):
    nodes, elems = structured_quad_mesh(1, 1)
    R = np.zeros(3) if bad == "residual_shape" else np.zeros(4)
    K = np.eye(4) if bad == "residual_shape" else np.eye(3)
    with pytest.raises(ValueError, match="output shapes"):
        solve_steps(nodes, elems, 1, lambda *a: (R, K), lambda f: {}, 1)


def test_stepping_drivers_validate_mesh_and_initial_state():
    nodes, elems = structured_quad_mesh(1, 1)
    fn = lambda ei, c, U_e, DU_e: (U_e.copy(), np.eye(4))
    with pytest.raises(ValueError, match="valid node"):
        solve_steps(nodes, np.array([[0, 1, 2, 4]]), 1, fn, lambda f: {}, 1)
    with pytest.raises(ValueError, match="U0 has shape"):
        solve_steps(nodes, elems, 1, fn, lambda f: {}, 1, U0=np.zeros(3))
    with pytest.raises(ValueError, match="n_steps"):
        solve_steps(nodes, elems, 1, fn, lambda f: {}, 0)


# --------------------------------------------------------------------------
# nonfinite output and unrelated exceptions
# --------------------------------------------------------------------------

def _trial_nan_spring(threshold):
    """Linear spring ``u - 1`` whose residual is NaN for ``u > threshold``."""
    def element_fn(ei, coords, U_e, DU_e):
        r = np.array([U_e[0] - 1.0])
        if U_e[0] > threshold:
            r[:] = np.nan
        return r, np.array([[0.4]])          # too soft: the full step overshoots
    return element_fn


def test_nonfinite_trial_backtracks_but_nonfinite_iterate_raises():
    # The first trial u = 2.5 is nonfinite; backtracking finds finite descent.
    U = solve_steps(ONE_NODE, ONE_ELEM, 1, _trial_nan_spring(2.0), lambda f: {}, 1,
                    tol=1e-10, max_newton=60)
    assert U[0] == pytest.approx(1.0)
    # A nonfinite iterate (here the initial state) is a failed step.
    with pytest.raises(NonfiniteElementOutput):
        solve_steps(ONE_NODE, ONE_ELEM, 1, _trial_nan_spring(-1.0), lambda f: {}, 1)


def test_public_adapter_nonfinite_trial_is_backtracked():
    """The compiled adapter raises on NaN output; a trial backtracks past it."""
    nodes, elems = structured_quad_mesh(2, 1)
    element, kernel = recording_element(len(elems))
    bc = {0: 0.0, 3: 0.0, 2: 0.0, 5: 0.0}
    reference = solve_steps(nodes, elems, 1, element.element_rk, lambda f: bc, 1,
                            batch_fn=element.element_rk_batch, tol=1e-10)
    kernel.calls.clear()
    # Only the first line-search trial (the second batch call) is nonfinite:
    # the half step is accepted, and Newton then converges to the same state.
    kernel.nonfinite = lambda t, dt, u, du: len(kernel.calls) == 2
    U = solve_steps(nodes, elems, 1, element.element_rk, lambda f: bc, 1,
                    batch_fn=element.element_rk_batch, tol=1e-10)
    np.testing.assert_allclose(U, reference, atol=1e-12)
    assert len(kernel.calls) == 6   # iterate, NaN trial, half step, iterate, trial, iterate
    kernel.nonfinite = lambda t, dt, u, du: True
    with pytest.raises(NonfiniteElementOutput):
        solve_steps(nodes, elems, 1, element.element_rk, lambda f: bc, 1,
                    batch_fn=element.element_rk_batch)


@pytest.mark.parametrize("error", [ValueError, RuntimeError, ZeroDivisionError])
def test_unrelated_element_exceptions_propagate_from_trials(error):
    def element_fn(ei, coords, U_e, DU_e):
        if U_e[0] > 0.0:
            raise error("element bug")
        return np.array([U_e[0] - 1.0]), np.array([[1.0]])

    with pytest.raises(error, match="element bug"):
        solve_steps(ONE_NODE, ONE_ELEM, 1, element_fn, lambda f: {}, 1)


# --------------------------------------------------------------------------
# multi-material assembly and solve
# --------------------------------------------------------------------------

def _two_material_mesh():
    """Two Quad4s sharing an edge; left = coupled (u, c), right = u only."""
    nodes, elems = structured_quad_mesh(2, 1)
    return nodes, elems[:1], elems[1:]


def test_assemble_groups_maps_components_and_accumulates_shared_nodes():
    nodes, left, right = _two_material_mesh()
    U = np.arange(12, dtype=float) / 10

    def coupled(ei, coords, U_e, DU_e):     # comps [0, 1], 8 local DOFs
        return U_e.copy(), np.diag(np.arange(1.0, 9.0))

    def u_only(ei, coords, U_e, DU_e):      # comps [0], 4 local DOFs
        return 2.0 * U_e, 2.0 * np.eye(4)

    groups = [{"elems": left, "comps": [0, 1], "element_fn": coupled},
              {"elems": right, "comps": [0], "element_fn": u_only}]
    R, K = assemble_groups(nodes, len(nodes), 2, groups, U, du_is_total=True)
    shared = [1, 4]                          # nodes on the common edge
    for n in range(len(nodes)):
        expected_u = (U[2 * n] if n in (0, 1, 3, 4) else 0.0) + (
            2.0 * U[2 * n] if n in (1, 2, 4, 5) else 0.0)
        assert R[2 * n] == pytest.approx(expected_u)
        assert R[2 * n + 1] == pytest.approx(U[2 * n + 1] if n in (0, 1, 3, 4) else 0.0)
    for n in shared:                         # both materials add to the u row
        local = [0, 1, 4, 3].index(n)
        assert K[2 * n, 2 * n] == pytest.approx(1.0 + 2 * local + 2.0)
    assert K[5, 5] == 0.0                    # the c row of a right-only node is empty


def test_assemble_groups_batch_matches_scalar_and_refuses_implicit_du():
    nodes, left, right = _two_material_mesh()
    Kl = np.arange(64.0).reshape(8, 8)
    scalar = lambda ei, c, U_e, DU_e: (Kl @ U_e + DU_e, Kl)
    batch = lambda c, U, DU: (U @ Kl.T + DU, np.broadcast_to(Kl, (len(U), 8, 8)))
    U = np.linspace(0.0, 1.0, 12)
    out = []
    for key, fn in (("element_fn", scalar), ("batch_fn", batch)):
        groups = [{"elems": np.vstack([left, right]), "comps": [0, 1], key: fn}]
        out.append(assemble_groups(nodes, len(nodes), 2, groups, U, U / 3))
    np.testing.assert_allclose(out[0][0], out[1][0], rtol=1e-14)
    np.testing.assert_array_equal(out[0][1].toarray(), out[1][1].toarray())
    with pytest.raises(ValueError, match="du_is_total"):
        assemble_groups(nodes, len(nodes), 2, groups, U)
    with pytest.raises(ValueError, match="comps"):
        assemble_groups(nodes, len(nodes), 2, [{**groups[0], "comps": [0, 0]}], U, U)


def test_multimaterial_commits_only_below_tol():
    """Defect regression: the solve used to accept a residual up to 1e-6
    whatever ``tol`` was, and to apply a 2^-30 update after an exhausted line
    search.  Both must now raise, as in ``solve_steps``."""
    def cubic(ei, coords, U_e, DU_e):
        R, K = _cubic_spring(U_e[0])
        return R - 0.002, K

    groups = [{"elems": ONE_ELEM, "comps": [0], "element_fn": cubic}]
    with pytest.raises(RuntimeError, match=r"did not converge in 1 Newton"):
        solve_steps_multimaterial(ONE_NODE, 1, 1, groups, lambda f: {}, 1,
                                  max_newton=1, tol=1e-12)
    U = solve_steps_multimaterial(ONE_NODE, 1, 1, groups, lambda f: {}, 1,
                                  max_newton=5, tol=1e-12)
    assert abs(U[0] ** 3 + U[0] - 0.002) < 1e-12

    calls = []

    def ascent(ei, coords, U_e, DU_e):
        calls.append(U_e[0])
        return np.array([4.0 * U_e[0] - 1.0]), np.array([[-4.0]])

    with pytest.raises(RuntimeError, match="line search exhausted"):
        solve_steps_multimaterial(ONE_NODE, 1, 1, [{"elems": ONE_ELEM, "comps": [0],
                                                     "element_fn": ascent}],
                                  lambda f: {}, 1, max_newton=40)
    assert len(calls) == 31


def test_multimaterial_final_update_commit_hook_and_state():
    commits, kincs = [], []

    class StubState:
        def begin_increment(self, kinc, time_begin=None, dt=None):
            kincs.append(kinc)

    groups = [{"elems": ONE_ELEM, "comps": [0],
               "element_fn": lambda ei, c, U_e, DU_e: (np.array([4.0 * U_e[0] - 2.0]),
                                                        np.array([[4.0]])),
               "commit": lambda U, U_prev: commits.append((U[0], U_prev[0]))}]
    U = solve_steps_multimaterial(ONE_NODE, 1, 1, groups, lambda f: {}, 2,
                                  max_newton=1, tol=1e-12, state=StubState())
    assert U[0] == pytest.approx(0.5)
    assert commits == [(0.5, 0.0), (0.5, 0.5)]
    assert kincs == [1, 2]
