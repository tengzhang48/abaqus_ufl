"""Deliberately broken compiled BVPs must fail the physical acceptance gates."""

import copy

import numpy as np
import pytest

from abaqus_ufl.fe import CompiledAbaqusElement, build_compiled_uel
from tools.boundary_value_common import heat_problem, solve_history, thermal_norms
from tools.check_heated_plate import assert_plate_accuracy, conditions as plate_conditions
from tools.check_thermal_bending import assert_bending_accuracy, conditions as strip_conditions


@pytest.fixture(scope="module")
def compiled_problem():
    problem = heat_problem(K=2.0)
    return problem, build_compiled_uel(problem)


@pytest.fixture(scope="module")
def plate_baseline(compiled_problem):
    problem, compiled = compiled_problem
    result = solve_history(problem, compiled, 6, 6, np.linspace(0, 1, 5), plate_conditions)
    assert_plate_accuracy(result)
    return result


def test_independent_q1_norm_quadrature_matches_polynomial_integrals():
    nodes = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    elements = np.array([[0, 1, 2, 3]])
    assert np.allclose(thermal_norms(nodes, elements, nodes[:, 0] + 2 * nodes[:, 1]), [8 / 3, 5])
    assert np.allclose(thermal_norms(nodes, elements, nodes[:, 0] * nodes[:, 1]), [1 / 9, 2 / 3])


def test_lost_old_temperature_fails_independent_heat_balance(compiled_problem, plate_baseline, monkeypatch):
    original = CompiledAbaqusElement.element_rk

    def lose_history(self, index, coordinates, current, increment):
        return original(self, index, coordinates, current, current)

    monkeypatch.setattr(CompiledAbaqusElement, "element_rk", lose_history)
    problem, compiled = compiled_problem
    with pytest.raises(AssertionError, match="heat balance"):
        solve_history(problem, compiled, 6, 6, np.linspace(0, 1, 5), plate_conditions)


def test_delayed_boundary_temperature_fails_fourier_gate(compiled_problem, plate_baseline):
    problem, compiled = compiled_problem
    result = solve_history(problem, compiled, 6, 6, np.linspace(0, 1, 5),
                           lambda nodes, time: plate_conditions(nodes, max(0, time - 0.25)))
    with pytest.raises(AssertionError, match="Fourier"):
        assert_plate_accuracy(result)


def test_prescribed_strip_displacements_cannot_claim_a_solved_bending_problem(compiled_problem):
    problem, compiled = compiled_problem

    def fixed_displacements(nodes, time):
        values = strip_conditions(nodes, time)
        values.update({3 * i + j: 0.0 for i in range(len(nodes)) for j in (0, 1)})
        return values

    result = solve_history(problem, compiled, 24, 4, np.linspace(0, 1, 5), fixed_displacements, length=8)
    with pytest.raises(AssertionError, match="free temperature and displacement"):
        assert_bending_accuracy(result)


def test_disabled_thermal_expansion_fails_independent_beam_gate(compiled_problem, monkeypatch):
    problem, compiled = compiled_problem
    baseline = solve_history(problem, compiled, 24, 4, np.linspace(0, 1, 5), strip_conditions, length=8)
    assert_bending_accuracy(baseline)
    original = CompiledAbaqusElement.__init__

    def disable_expansion(self, module, weakform, *args, **kwargs):
        original(self, module, weakform, *args, **kwargs)
        self.props[weakform._mat.props_names.index("alpha")] = 0

    monkeypatch.setattr(CompiledAbaqusElement, "__init__", disable_expansion)
    result = solve_history(problem, compiled, 24, 4, np.linspace(0, 1, 5), strip_conditions, length=8)
    with pytest.raises(AssertionError, match="reference"):
        assert_bending_accuracy(result)


def test_nonfinite_history_is_never_accepted(plate_baseline):
    result = copy.deepcopy(plate_baseline)
    result["solutions"][1, 2] = np.nan
    with pytest.raises(AssertionError, match="invalid"):
        assert_plate_accuracy(result)
