"""Broken controls attack compiled cyclic J2 host conventions and transactions."""

import copy
from types import SimpleNamespace

import numpy as np
import pytest

from abaqus_ufl.fe import newton_solve, structured_quad_mesh
from tools import plasticity_bvp_common as host_module
from tools.plasticity_bvp_common import (
    J2Quad4Host, assert_trial_purity, build_compiled_j2, conditions,
    endpoint_work_balance, oracle_assembly, point_closed_form_gate, project_plastic_strain,
    solve_history, validate_history,
)


@pytest.fixture(scope="module")
def compiled_j2():
    return build_compiled_j2()


@pytest.fixture(scope="module")
def cycle(compiled_j2):
    result = solve_history(compiled_j2, 4)
    validate_history(result)
    return result


def test_compiled_shear_and_normal_abi_have_independent_closed_forms(compiled_j2):
    assert point_closed_form_gate(compiled_j2) < 1e-10


def assert_noncoaxial_plastic_tangent(module):
    """Mixed-state element tangent against FD of independent fixed-old forces."""
    nodes, elements = structured_quad_mesh(1, 1)
    host = J2Quad4Host(module, nodes, elements)
    x, y = nodes.T
    old = np.column_stack((.015 * y, np.zeros(len(nodes)))).ravel()
    zero = np.zeros_like(old)
    old_reference = oracle_assembly(nodes, elements, old, zero,
                                    np.zeros((1, 4, 6)), np.zeros((1, 4)))
    # A local material-state tangent check needs a valid old constitutive
    # state, not an equilibrated global BVP. The old pure-shear point state is
    # set entirely by the independent oracle, with all six stress components.
    host.stress, host.ep, host.strain = old_reference[1:]
    host.accepted_displacement = old.copy()
    host._protect_accepted_state()
    increment = np.column_stack((.01 * x - .004 * y + .0015 * x * y,
                                  -.002 * x - .004 * y - .001 * x * y)).ravel()
    current = old + increment
    connectivity = elements[0]
    dofs = (2 * connectivity[:, None] + np.arange(2)).ravel()
    candidate = host.evaluate_element(0, current[dofs], increment[dofs])
    assert np.all(candidate[3] > host.ep[0] + 1e-4)
    # Stress eigendirections actually change; this gate exercises coupling
    # components absent from the existing proportional pure-shear check.
    def tensor(s):
        return np.array([[s[0], s[3], s[4]], [s[3], s[1], s[5]], [s[4], s[5], s[2]]])

    before, after = tensor(host.stress[0, 0]), tensor(candidate[2][0])
    assert np.linalg.norm(before @ after - after @ before) > 1e-4
    derivative = np.zeros((8, 8))
    delta = 1e-7
    for column, dof in enumerate(dofs):
        plus, minus = current.copy(), current.copy()
        plus[dof] += delta
        minus[dof] -= delta
        positive = oracle_assembly(nodes, elements, plus, old, host.stress, host.ep)[0]
        negative = oracle_assembly(nodes, elements, minus, old, host.stress, host.ep)[0]
        derivative[:, column] = (positive[dofs] - negative[dofs]) / (2 * delta)
    if np.max(abs(candidate[1] - derivative)) > 2e-7:
        raise AssertionError("mixed plastic tangent disagrees with independent force derivative")


def test_noncoaxial_multiaxial_assembled_tangent_matches_independent_forces(compiled_j2):
    assert_noncoaxial_plastic_tangent(compiled_j2)


def test_dropped_tangent_cross_components_fail_mixed_state_oracle(compiled_j2, monkeypatch):
    original = host_module.call_umat

    def drop_cross_coupling(*args, **kwargs):
        stress, ep, tangent = original(*args, **kwargs)
        return stress, ep, np.diag(np.diag(tangent))

    monkeypatch.setattr(host_module, "call_umat", drop_cross_coupling)
    with pytest.raises(AssertionError, match="mixed plastic tangent"):
        assert_noncoaxial_plastic_tangent(compiled_j2)


def test_cycle_has_free_equilibrium_yield_unload_reverse_and_work_oracles(cycle):
    metrics = validate_history(cycle)
    assert metrics["yielded_gauss_points_at_peak"] == 64
    assert metrics["elastic_unloading_gauss_points"] == 64
    assert metrics["reverse_plastic_gauss_points"] >= 16
    assert cycle["endpoint_work_balance_errors"].max() < 1e-10
    assert cycle["state_oracle_error"] < 1e-10


def test_half_engineering_shear_fails_independent_global_oracle(compiled_j2, monkeypatch):
    original = host_module.quad4_kinematics

    def tensor_shear_in_engineering_slot(coordinates):
        B, weights, points = original(coordinates)
        B[:, 3, :] *= .5
        return B, weights, points

    monkeypatch.setattr(host_module, "quad4_kinematics", tensor_shear_in_engineering_slot)
    with pytest.raises(AssertionError, match="independent radial return"):
        solve_history(compiled_j2, 4)


def test_lost_accepted_state_cannot_pass_the_cyclic_oracle(compiled_j2, monkeypatch):
    original = host_module.call_umat

    def reset_state(module, stress, ep, strain, increment, time, dt, *args):
        return original(module, stress, 0., strain, increment, time, dt, *args)

    monkeypatch.setattr(host_module, "call_umat", reset_state)
    with pytest.raises(AssertionError, match="independent radial return|plastic strain decreased"):
        solve_history(compiled_j2, 4)


def test_rejected_trial_has_exact_rollback_and_detects_early_commit(compiled_j2, monkeypatch):
    nodes, elements = structured_quad_mesh(4, 4)
    host = J2Quad4Host(compiled_j2, nodes, elements)
    zero = np.zeros(2 * len(nodes))
    assert_trial_purity(host, zero, zero)
    assert host.version == 0
    original = host.evaluate_element

    def commit_during_trial(ei, current, increment):
        output = original(ei, current, increment)
        new_state = host.ep.copy()
        new_state[ei] = output[3]
        host.ep = new_state
        host.version += 1
        return output

    monkeypatch.setattr(host, "evaluate_element", commit_during_trial)
    with pytest.raises(AssertionError, match="rejected plasticity trial changed accepted state"):
        assert_trial_purity(host, zero, zero)


def test_failed_newton_and_cutback_never_commit_material_state(compiled_j2):
    nodes, elements = structured_quad_mesh(4, 4)
    host = J2Quad4Host(compiled_j2, nodes, elements)
    zero = np.zeros(2 * len(nodes))
    before = [value.copy() for value in (host.stress, host.strain, host.ep)]
    prescribed = conditions(nodes, 1.)
    with pytest.raises(RuntimeError, match="did not converge"):
        newton_solve(nodes, elements, 2, host.element_rk, prescribed,
                     U0=zero, U_previous=zero, max_iter=1, tol=1e-10)
    with pytest.raises(AssertionError, match="unconverged plasticity trial"):
        candidate = zero.copy()
        for dof, value in prescribed.items():
            candidate[dof] = value
        host.commit_equilibrium(candidate, prescribed)
    assert host.version == 0
    assert all(np.array_equal(a, b) for a, b in zip(before, (host.stress, host.strain, host.ep)))

    def cutback(*args):
        output = list(compiled_j2.drive_umat(*args))
        output[-1] = .5
        return tuple(output)

    host.module = SimpleNamespace(drive_umat=cutback)
    with pytest.raises(RuntimeError, match="cutback"):
        host.element_rk(0, nodes[elements[0]], np.ones(8) * .01, np.ones(8) * .01)
    assert host.version == 0
    assert all(np.array_equal(a, b) for a, b in zip(before, (host.stress, host.strain, host.ep)))


def test_fixing_every_displacement_cannot_claim_a_mesh_solve(compiled_j2):
    def all_prescribed(nodes, time):
        return {i: 0. for i in range(2 * len(nodes))}

    with pytest.raises(AssertionError, match="must solve free displacement DOFs"):
        solve_history(compiled_j2, 4, boundary_conditions=all_prescribed)


def test_dropping_plane_strain_normal_stress_fails_global_work_identity(cycle):
    step = 8
    stress, old_stress = cycle["stress"][step].copy(), cycle["stress"][step - 1].copy()
    stress[:, :, 2] = old_stress[:, :, 2] = 0.
    with pytest.raises(AssertionError, match="endpoint-work balance"):
        endpoint_work_balance(cycle["solutions"][step], cycle["solutions"][step - 1],
                              cycle["residuals"][step - 1], conditions(cycle["nodes"], cycle["times"][step]),
                              cycle["strain"][step], stress, cycle["ep"][step],
                              cycle["strain"][step - 1], old_stress, cycle["ep"][step - 1], cycle["weights"])


def test_visual_projection_averages_incident_reference_element_means():
    nodes, elements = structured_quad_mesh(2, 1, 2., 1.)
    result = {"nodes": nodes, "elements": elements, "weights": np.ones((2, 4)) / 4,
              "ep": np.array([[[0., 1., 2., 3.], [4., 5., 6., 7.]]])}
    assert np.allclose(project_plastic_strain(result), [[1.5, 3.5, 5.5, 1.5, 3.5, 5.5]])


def test_nonfinite_gauss_point_history_cannot_be_exported(cycle):
    broken = copy.deepcopy(cycle)
    broken["ep"][3, 0, 0] = np.nan
    with pytest.raises(AssertionError, match="invalid plasticity history: ep"):
        validate_history(broken)
