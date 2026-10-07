"""Compiled Ogden mesh gates and deliberately broken physics controls."""

import copy

import numpy as np
import pytest

from abaqus_ufl.fe import CompiledAbaqusElement, build_compiled_uel
from tools.check_ogden_bvp import (
    OgdenShearProblem, PROPERTIES, SHEAR, audit_alpha_two_limit,
    audit_element_tangent, audit_history, conditions, shear_patch, solve_history,
)
from tools.ogden_bvp_reference import alpha_two_piola, piola_and_energy


@pytest.fixture(scope="module")
def compiled_problem():
    problem = OgdenShearProblem(**PROPERTIES)
    return problem, build_compiled_uel(problem)


@pytest.fixture(scope="module")
def baseline(compiled_problem):
    problem, compiled = compiled_problem
    result = solve_history(problem, compiled, 4, np.linspace(0, 1, 7))
    audit_history(result)
    return result


def test_compiled_noncoaxial_element_energy_tangent_and_closed_form_shear(compiled_problem):
    problem, compiled = compiled_problem
    audit_element_tangent(compiled, problem)
    shear_patch(compiled, problem)


def test_alpha_two_limit_at_nonisochoric_noncoaxial_states(compiled_problem, baseline):
    _, compiled = compiled_problem
    audit_alpha_two_limit(compiled, baseline)
    F = np.array([[1.18, 0.35, 0.0], [-0.12, 0.93, 0.0], [0.0, 0.0, 1.0]])
    assert abs(np.linalg.det(F) - 1) > 0.1
    svd, _, _ = piola_and_energy(F, mu=1.0, alpha=2.0, K=10.0)
    eig_free = alpha_two_piola(F, mu=1.0, K=10.0)
    assert np.allclose(svd, eig_free, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("angle", [0.1, 0.51, 1.0])
def test_rigid_rotation_has_zero_force_and_objective_compiled_tangent(compiled_problem, angle):
    problem, compiled = compiled_problem
    element = CompiledAbaqusElement(compiled, problem)
    nodes = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    _, initial_tangent = element.element_rk(0, nodes, np.zeros(8), np.zeros(8))
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    u = (nodes @ (rotation - np.eye(2)).T).ravel()
    residual, tangent = element.element_rk(0, nodes, u, u)
    assert np.max(abs(residual)) < 1e-12
    nodal_rotation = np.kron(np.eye(4), rotation)
    expected = nodal_rotation @ initial_tangent @ nodal_rotation.T
    assert np.allclose(tangent, expected, rtol=2e-10, atol=2e-11)


def test_wrong_compiled_ogden_exponent_fails_independent_mesh_oracle(compiled_problem, monkeypatch):
    problem, compiled = compiled_problem
    original = CompiledAbaqusElement.__init__

    def wrong_exponent(self, module, weakform, *args, **kwargs):
        original(self, module, weakform, *args, **kwargs)
        self.props[weakform._mat.props_names.index("alpha")] = 2.0

    monkeypatch.setattr(CompiledAbaqusElement, "__init__", wrong_exponent)
    result = solve_history(problem, compiled, 4, np.linspace(0, 1, 7))
    with pytest.raises(AssertionError, match="principal-stretch equilibrium"):
        audit_history(result)


def test_prescribed_affine_illustration_cannot_claim_a_solved_block(compiled_problem):
    problem, compiled = compiled_problem

    def affine_everywhere(nodes, parameter):
        result = conditions(nodes, parameter)
        result.update({2 * i: SHEAR * parameter * nodes[i, 1] for i in range(len(nodes))})
        result.update({2 * i + 1: 0.0 for i in range(len(nodes))})
        return result

    result = solve_history(problem, compiled, 4, np.linspace(0, 1, 7), boundary_conditions=affine_everywhere)
    with pytest.raises(AssertionError, match="free side and interior"):
        audit_history(result)


def test_wrong_grip_schedule_fails_boundary_value_gate(compiled_problem):
    problem, compiled = compiled_problem
    result = solve_history(problem, compiled, 4, np.linspace(0, 1, 7),
                           boundary_conditions=lambda nodes, s: conditions(nodes, s / 2))
    with pytest.raises(AssertionError, match="grip boundary conditions"):
        audit_history(result)


def test_nonfinite_history_cannot_be_published(baseline):
    result = copy.deepcopy(baseline)
    result["solutions"][2, 0] = np.nan
    with pytest.raises(AssertionError, match="invalid Ogden"):
        audit_history(result)
