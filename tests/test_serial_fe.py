"""Serial FE failures must not become apparent equilibrium or stale history."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from abaqus_ufl.fe import (CompiledAbaqusElement, assemble, boundary_nodes,
                           newton_solve, structured_quad_mesh)
from abaqus_ufl.fe.compiled_element import _signature

ROOT = Path(__file__).resolve().parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_serial_fe_mesh_oracles():
    module = _load("_serial_mesh_oracles", ROOT / "tools/check_fe_runtime.py")
    metrics = module.check()
    assert metrics["affine_patch_max_abs_error"] < 1e-9
    assert metrics["transient_mode_max_abs_error"] < 1e-10
    assert 1.9 < metrics["spatial_convergence_rate"] < 2.1


@pytest.mark.parametrize("stiffness", [1.0, 2e11])
def test_convergence_tolerance_scales_with_reactions(stiffness):
    # Each Quad4 couples its four edges; a linear field is the exact solution.
    ring = np.array([[2, -1, 0, -1], [-1, 2, -1, 0], [0, -1, 2, -1], [-1, 0, -1, 2]], float)
    nodes, elems = structured_quad_mesh(3, 3)
    exact = nodes[:, 0] + 2 * nodes[:, 1]
    conditions = {int(i): float(exact[i]) for i in boundary_nodes(nodes)}
    element = lambda index, coordinates, U, DU: (stiffness * ring @ U, stiffness * ring)
    history = []
    U = newton_solve(nodes, elems, 1, element, conditions, record_iteration=history.append)
    assert np.allclose(U, exact, rtol=0, atol=1e-12)
    assert history[0]["iteration"] == 0
    assert history[0]["relative_residual"] > 1e-9
    assert history[-1]["relative_residual"] <= 1e-9
    assert all(np.isfinite(h["free_residual_inf"]) for h in history)


def _fake_element(response):
    model = _load("_serial_fake_model", ROOT / "examples/scalar_diffusion_uel/build.py")
    problem = model.HeatDiffusionProblem(G=3.5)
    module = SimpleNamespace(_abaqus_ufl_signature=_signature(problem), drive_uel=response)
    return CompiledAbaqusElement(module, problem, dt=0.2)


def test_compiled_adapter_sign_flags_instance_properties_and_increment_time():
    calls = []

    def response(*args):
        calls.append(args)
        return np.arange(12, dtype=float), np.eye(12), np.zeros(1), 1.0

    element = _fake_element(response)
    element.begin_increment(3)
    nodes, elems = structured_quad_mesh(1, 1)
    r, k = element.element_rk(0, nodes[elems[0]], np.zeros(12), np.zeros(12))
    assert np.array_equal(r, -np.arange(12))
    assert np.array_equal(k, np.eye(12))
    assert calls[0][4][0] == 3.5  # constructor override, not the class defaults
    assert np.array_equal(calls[0][9], [1, 0, 1, 0, 0, 0])
    assert np.allclose(calls[0][6], [0.4, 0.4])  # start, not end, of increment 3
    assert calls[0][13] == 3


@pytest.mark.parametrize("pnewdt", [0.5, 0.0, float("nan")])
def test_cutback_zero_work_is_not_equilibrium(pnewdt):
    element = _fake_element(lambda *args: (np.zeros(12), np.zeros((12, 12)), np.zeros(1), pnewdt))
    nodes, elems = structured_quad_mesh(1, 1)
    with pytest.raises(RuntimeError, match="cutback"):
        newton_solve(nodes, elems, 3, element.element_rk, {})


def test_nonfinite_output_is_not_equilibrium():
    element = _fake_element(lambda *args: (np.full(12, np.nan), np.eye(12), np.zeros(1), 1.0))
    nodes, elems = structured_quad_mesh(1, 1)
    with pytest.raises(RuntimeError, match="nonfinite"):
        newton_solve(nodes, elems, 3, element.element_rk, {})


def test_assembly_requires_explicit_history():
    nodes, elems = structured_quad_mesh(1, 1)
    fn = lambda index, coordinates, U, DU: (DU.copy(), np.eye(4))
    U = np.ones(4)
    with pytest.raises(ValueError, match="U_previous"):
        assemble(nodes, elems, 1, fn, U)
    residual, _ = assemble(nodes, elems, 1, fn, U, np.zeros(4))
    assert np.array_equal(residual, np.zeros(4))
    from_zero, _ = assemble(nodes, elems, 1, fn, U, du_is_total=True)
    assert np.array_equal(from_zero, U)


def test_nonconvergence_and_singular_matrix_raise():
    nodes, elems = structured_quad_mesh(1, 1)
    constant_residual = lambda *args: (np.ones(4), np.eye(4))
    with pytest.raises(RuntimeError, match="line search"):
        newton_solve(nodes, elems, 1, constant_residual, {})
    singular = lambda *args: (np.ones(4), np.zeros((4, 4)))
    with pytest.raises(RuntimeError, match="linear solve"):
        newton_solve(nodes, elems, 1, singular, {})


def test_unsupported_layouts_are_rejected_before_build():
    from abaqus_ufl.fe import build_compiled_uel
    model = _load("_serial_unsupported_model", ROOT / "examples/thermo_mechanics_quad8/build.py")
    with pytest.raises(ValueError, match="degree-one"):
        build_compiled_uel(model.ThermoMechanicalQuad8())
