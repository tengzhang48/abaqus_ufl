"""Transient thermal bending of a clamped, traction-free plane-strain strip."""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.boundary_value_common import (
    beam_tip_reference, boundary_region, heat_problem, history_mesh, ramp,
    simulation_setup, solve_history, strip_reference, validate_history,
)
from tools.livebench_data import comparison_curve

LENGTH, HEIGHT, RATE, TEMPERATURE, DT = 8.0, 1.0, 4.0, 1.0, 0.05


def conditions(nodes, time):
    left = np.flatnonzero(np.isclose(nodes[:, 0], 0))
    result = {3 * int(i) + j: 0.0 for i in left for j in (0, 1)}
    for i in np.flatnonzero(np.isclose(nodes[:, 1], 0) | np.isclose(nodes[:, 1], HEIGHT)):
        result[3 * int(i) + 2] = TEMPERATURE * (2 * nodes[i, 1] / HEIGHT - 1) * ramp(time, RATE)
    return result


def tip_values(result):
    tip = np.flatnonzero(np.all(np.isclose(result["nodes"], [LENGTH, HEIGHT / 2]), axis=1))
    return result["solutions"][:, 3 * int(tip[0]) + 1]


def assert_bending_accuracy(result, *, finest=False):
    validate_history(result)
    temperature = result["solutions"][:, 2::3]
    exact = strip_reference(result["nodes"], result["times"], backward_euler=True)
    nx = len(np.unique(result["nodes"][:, 0])) - 1
    ny = len(np.unique(result["nodes"][:, 1])) - 1
    expected_dofs = {"temperature": (nx + 1) * (ny - 1), "displacement": 2 * nx * (ny + 1)}
    if len(result["free_dof_counts"]) != len(result["times"]) - 1 or any(c != expected_dofs for c in result["free_dof_counts"]):
        raise AssertionError("strip must solve free temperature and displacement DOFs")
    temperature_error = float(np.max(abs(temperature - exact)))
    guide = beam_tip_reference(result["times"])
    relative_tip_error = float(abs(tip_values(result)[-1] / guide[-1] - 1))
    if (temperature_error > 0.032 * (4 / ny) ** 2 or relative_tip_error > 0.08
            or (finest and (temperature_error > 2e-3 or relative_tip_error > 0.03))):
        raise AssertionError("strip disagrees with the independent heat or slender-beam reference")
    if tip_values(result)[-1] >= -0.01 or max(result["dissipation"]) < 0.1:
        raise AssertionError("thermal expansion did not drive the prescribed bending regime")
    return temperature_error, relative_tip_error


def check(*, record=None):
    from abaqus_ufl.fe import build_compiled_uel
    problem = heat_problem(K=2.0)
    compiled = build_compiled_uel(problem)
    times = np.linspace(0, 1, 21)
    coarse = solve_history(problem, compiled, 24, 4, times, conditions, length=LENGTH, height=HEIGHT)
    main = solve_history(problem, compiled, 48, 8, times, conditions, length=LENGTH, height=HEIGHT)
    coarse_error, coarse_tip_error = assert_bending_accuracy(coarse)
    temperature_error, tip_error = assert_bending_accuracy(main, finest=True)
    if not (temperature_error < coarse_error / 3 and tip_error < coarse_tip_error):
        raise AssertionError("thermal bending did not improve with mesh refinement")
    nodes, elems = main["nodes"], main["elements"]
    left = np.flatnonzero(np.isclose(nodes[:, 0], 0))
    top = np.flatnonzero(np.isclose(nodes[:, 1], HEIGHT))
    bottom = np.flatnonzero(np.isclose(nodes[:, 1], 0))
    force_errors, moment_errors = [], []
    for state, residual in zip(main["solutions"][1:], main["residuals"]):
        reactions = residual.reshape(-1, 3)[left, :2]
        positions = nodes[left] + state.reshape(-1, 3)[left, :2]
        force_errors.append(float(np.max(abs(np.sum(reactions, axis=0)))))
        moment_errors.append(float(abs(np.sum(positions[:, 0] * reactions[:, 1]
                                              - positions[:, 1] * reactions[:, 0]))))
        if np.sum(residual[3 * top + 2]) <= 0 or np.sum(residual[3 * bottom + 2]) >= 0:
            raise AssertionError("hot and cold baths have invalid supplied-heat signs")
    if max(force_errors + moment_errors) > 1e-8:
        raise AssertionError("clamp reactions violate force or moment balance")
    print("[PASS] coupled free temperature/displacement DOFs solved for 20 increments")
    print("[PASS] temperature Fourier series, thermal dissipation and heat balance")
    print("[PASS] clamp reactions balance force and moment; refined tip beam error {:.2%}".format(tip_error))

    boundary = [boundary_region(nodes, 0, 0, "Clamped: u₁ = u₂ = 0", "displacement"),
                boundary_region(nodes, 1, HEIGHT, "Hot bath: T = 1 − exp(−4t)", "temperature"),
                boundary_region(nodes, 1, 0, "Cold bath: T = −(1 − exp(−4t))", "temperature")]
    guide = beam_tip_reference(times)
    row = np.flatnonzero(np.isclose(nodes[:, 1], HEIGHT / 2))
    column = np.flatnonzero(np.isclose(nodes[:, 0], LENGTH / 2))
    temperature_reference = strip_reference(nodes, times, backward_euler=True)
    iterations = [h[-1]["iteration"] for h in main["iterations"]]
    relative_residual = [h[-1]["relative_residual"] for h in main["iterations"]]
    plots = [
        history_mesh("deformation", "Thermally driven bending", main,
            "One end is clamped; all other mechanical edges are traction-free. "
            "Temperature and displacement are solved together across 384 Quad4 elements. "
            "Colors show authoritative nodal displacement magnitude; deformation is displayed at ×10.",
            field="displacement", boundaries=boundary),
        history_mesh("temperature", "Heat diffusion through the bending strip", main,
            "Initially T = 0. Opposite baths ramp to ±1; both end faces are insulated. "
            "Heat transfer drives quasi-static finite-strain deformation. This material has no mechanical-heating feedback.",
            boundaries=boundary[1:]),
        comparison_curve("tip-history", "Deflection of the free tip", "Time (model units)", "Vertical tip displacement",
            times, tip_values(main), guide,
            "The independent small-strain slender-beam approximation uses the transient Fourier temperature moment. "
            "All recorded samples are within 8%; the final tip is required within 3%. "
            "Full end clamping and finite strain prevent exact beam agreement.",
            reference_label="Slender beam approximation", rtol=0.08, atol=1e-7),
        comparison_curve("temperature-profile", "Temperature through the thickness", "Reference Y", "Temperature rise",
            nodes[column, 1], main["solutions"][-1, 2::3][column], temperature_reference[-1, column],
            "Independent continuum backward-Euler Fourier series at the strip midpoint, t = 1.",
            reference_label="Continuum backward Euler", rtol=0, atol=2e-3),
        {"kind": "curve", "id": "centerline", "title": "Deformed centerline", "description":
            "Solved centerline displacement compared with the approximate curvature from the independent thermal moment. "
            "The beam guide omits clamp end effects and finite strain.",
            "x_label": "Reference X", "y_label": "Vertical displacement", "x": nodes[row, 0].tolist(),
            "x_scale": "linear", "y_scale": "linear", "series": [
                {"label": "Compiled FE", "role": "computed", "values": main["solutions"][-1, 3 * row + 1].tolist()},
                {"label": "Slender beam guide", "role": "guide", "values": (guide[-1] * (nodes[row, 0] / LENGTH) ** 2).tolist()}]},
        {"kind": "curve", "id": "equilibrium", "title": "Equilibrium at each accepted step",
            "description": "Maximum free residual divided by the force/heat scale at each accepted step. "
                           "Newton solves use tolerance 10⁻¹⁰; no unconverged trial is retained.",
            "x_label": "Time (model units)", "y_label": "Relative free residual", "x": times[1:].tolist(),
            "x_scale": "linear", "y_scale": "linear", "series": [
                {"label": "Compiled FE", "role": "computed", "values": relative_residual},
                {"label": "Acceptance tolerance", "role": "guide", "values": [1e-10] * len(relative_residual)}]},
    ]
    setup = simulation_setup(problem, main, "Plane-strain strip: length 8, thickness 1. Consistent model units.", [
        "Left end: u₁ = u₂ = 0 (clamped).",
        "All other mechanical edges: zero traction.",
        "Top/bottom: T = ±(1 − exp(−4t)). Left/right: zero normal heat flux."],
        ["Div P = 0; quasi-static finite-strain mechanics.",
         "ρcp ∂T/∂t = k ΔT. Heat transfer drives thermal expansion; no heating feedback."])
    for plot in plots:
        plot["setup"] = setup
        if record is not None:
            record(plot)
    return {"temperature_max_abs_error": temperature_error, "final_tip_relative_beam_error": tip_error,
            "coarse_tip_relative_beam_error": coarse_tip_error, "final_tip_displacement": float(tip_values(main)[-1]),
            "force_balance_max_abs_error": max(force_errors), "moment_balance_max_abs_error": max(moment_errors),
            "heat_balance_max_abs_error": max(main["heat_balance_errors"]),
            "thermal_norm_balance_max_abs_error": max(main["thermal_norm_balance_errors"]),
            "boundary_condition_max_abs_error": main["boundary_condition_error"],
            "max_newton_iterations": max(iterations), "increments": len(times) - 1}


if __name__ == "__main__":
    for name, value in check().items():
        print("{} = {:.8e}".format(name, value))
