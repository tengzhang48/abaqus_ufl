"""A transient, spatially two-dimensional heated-plate boundary-value problem."""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.boundary_value_common import (
    boundary_region, heat_content, heat_problem, history_mesh, plate_heat_reference,
    plate_reference, ramp, simulation_setup, solve_history, validate_history,
)
from tools.livebench_data import comparison_curve

RATE, END_TIME, DT = 2.0, 1.0, 0.05


def conditions(nodes, time):
    from abaqus_ufl.fe import boundary_nodes
    result = {3 * i + j: 0.0 for i in range(len(nodes)) for j in (0, 1)}
    for i in boundary_nodes(nodes):
        result[3 * int(i) + 2] = (np.sin(np.pi * nodes[i, 0]) * ramp(time, RATE)
                                  if np.isclose(nodes[i, 1], 1) else 0.0)
    return result


def assert_plate_accuracy(result, *, finest=False):
    validate_history(result)
    times, nodes = result["times"], result["nodes"]
    computed = result["solutions"][:, 2::3]
    discrete = plate_reference(nodes, times, rate=RATE, backward_euler=True)
    continuous = plate_reference(nodes, times, rate=RATE)
    nx, ny = len(np.unique(nodes[:, 0])) - 1, len(np.unique(nodes[:, 1])) - 1
    expected_dofs = {"temperature": (nx - 1) * (ny - 1), "displacement": 0}
    if len(result["free_dof_counts"]) != len(times) - 1 or any(c != expected_dofs for c in result["free_dof_counts"]):
        raise AssertionError("plate must solve every interior temperature DOF")
    spatial = float(np.max(abs(computed[-1] - discrete[-1])))
    continuous_error = float(np.max(abs(computed[-1] - continuous[-1])))
    if spatial > 0.012 * (6 / nx) ** 2 or (finest and (spatial > 6e-4 or continuous_error > 1.5e-3)):
        raise AssertionError("heated plate disagrees with the independent Fourier solution")
    if np.max(computed[-1]) < 0.5 or np.min(computed) < -1e-9:
        raise AssertionError("the plate did not enter the prescribed heating regime")
    center = np.flatnonzero(np.all(np.isclose(nodes, [0.5, 0.5]), axis=1))
    row = np.flatnonzero(np.isclose(nodes[:, 1], 0.5))
    if not len(center) or computed[-1, center[0]] < 0.05 or np.ptp(computed[-1, row]) < 0.05:
        raise AssertionError("heat did not diffuse into the two-dimensional interior")
    return spatial, continuous_error


def _convergence(identifier, title, x, error, order, description):
    return {"kind": "curve", "id": identifier, "title": title, "description": description,
            "x_label": "Mesh spacing h" if identifier == "spatial-convergence" else "Time step Δt",
            "y_label": "Maximum nodal temperature difference", "x": list(x),
            "x_scale": "log", "y_scale": "log", "series": [
                {"label": "Compiled FE", "role": "computed", "values": list(error)},
                {"label": "Order {} guide".format(order), "role": "guide",
                 "values": [error[0] * (v / x[0]) ** order for v in x]}]}


def check(*, record=None):
    from abaqus_ufl.fe import build_compiled_uel
    problem = heat_problem()
    compiled = build_compiled_uel(problem)
    times = np.linspace(0, END_TIME, round(END_TIME / DT) + 1)
    runs, spatial_errors = {}, []
    for n in (6, 12, 24):
        run = solve_history(problem, compiled, n, n, times, conditions)
        error, continuous_error = assert_plate_accuracy(run, finest=n == 24)
        runs[n] = run
        spatial_errors.append(error)
        print("[PASS] {}x{} plate, 20 accepted increments, spatial error {:.6e}".format(n, n, error))
    spatial_rates = np.log(np.array(spatial_errors[:-1]) / spatial_errors[1:]) / np.log(2)
    if not np.all((spatial_rates > 1.8) & (spatial_rates < 2.2)):
        raise AssertionError("plate spatial refinement is not second order")

    # Same-mesh step differences remove the spatial error floor. The separate
    # Fourier comparisons above still bound actual continuum accuracy.
    temporal = {DT: runs[12]}
    for step in (0.2, 0.1, 0.025):
        temporal[step] = solve_history(problem, compiled, 12, 12,
            np.linspace(0, END_TIME, round(END_TIME / step) + 1), conditions)
    time_steps = [0.05, 0.1, 0.2]
    time_errors = [float(np.max(abs(temporal[step]["solutions"][-1, 2::3]
                          - temporal[step / 2]["solutions"][-1, 2::3]))) for step in time_steps]
    time_rates = np.log(np.array(time_errors[1:]) / time_errors[:-1]) / np.log(2)
    if not np.all((time_rates > 0.8) & (time_rates < 1.2)):
        raise AssertionError("compiled plate time refinement is not first order")
    print("[PASS] independent Fourier solution, second-order space and first-order time refinement")

    main = runs[24]
    nodes, elems = main["nodes"], main["elements"]
    reference = plate_reference(nodes, times, rate=RATE, backward_euler=True)
    center = int(np.flatnonzero(np.all(np.isclose(nodes, [0.5, 0.5]), axis=1))[0])
    row = np.flatnonzero(np.isclose(nodes[:, 0], 0.5))
    content = [heat_content(nodes, elems, state[2::3]) for state in main["solutions"]]
    content_reference = plate_heat_reference(times, rate=RATE, backward_euler=True)
    if content[-1] < 0.1 or max(main["dissipation"]) < 0.1:
        raise AssertionError("heated plate has no heat input or diffusion dissipation")
    plots = [
        history_mesh("temperature", "Heating of a 2D plate", main,
            "A spatially varying top-edge bath heats an initially cold plate. "
            "All 529 interior temperature DOFs are solved at each of 20 increments; the other three edges remain at zero.",
            boundaries=[boundary_region(nodes, 1, 1, "T = sin(πX)(1 − exp(−2t))", "temperature"),
                        boundary_region(nodes, 1, 0, "T = 0", "temperature"),
                        boundary_region(nodes, 0, 0, "T = 0", "temperature"),
                        boundary_region(nodes, 0, 1, "T = 0", "temperature")]),
        comparison_curve("center-history", "Temperature at the plate center", "Time (model units)",
            "Temperature rise (model units)", times, main["solutions"][:, 3 * center + 2], reference[:, center],
            "The independent sine-series solution of the continuum backward-Euler equations isolates spatial error. "
            "A continuous-time series is shown separately.", reference_label="Continuum backward Euler", atol=6e-4, rtol=0,
            extra_series=[{"label": "Continuous time", "role": "guide",
                           "values": plate_reference(nodes[[center]], times, rate=RATE)[:, 0].tolist()}]),
        comparison_curve("temperature-profile", "Vertical profile through the center", "Reference Y", "Temperature rise",
            nodes[row, 1], main["solutions"][-1, 2::3][row], reference[-1, row],
            "Final profile at X = 0.5, t = 1, Δt = 0.05, 24 × 24 mesh.",
            reference_label="Continuum backward Euler", atol=6e-4, rtol=0),
        comparison_curve("heat-content", "Heat stored in the plate", "Time (model units)", "Integrated heat content",
            times, content, content_reference,
            "Reference-domain heat content ∫ρcp T dA is independently integrated from the Fourier series. "
            "Its increment also balances the supplied heat from all prescribed thermal boundaries.",
            reference_label="Integrated Fourier solution", atol=2e-4, rtol=0),
        _convergence("spatial-convergence", "Spatial convergence of the plate",
            [1 / 24, 1 / 12, 1 / 6], spatial_errors[::-1], 2,
            "Final nodal errors against the continuum backward-Euler Fourier series on 6², 12² and 24² meshes. "
            "Measured orders: {:.3f}, {:.3f}. Time step stays at 0.05.".format(*spatial_rates)),
        _convergence("temporal-convergence", "Time convergence of the plate", time_steps, time_errors, 1,
            "Final differences between Δt and Δt/2 on the same 12 × 12 mesh. "
            "Measured orders: {:.3f}, {:.3f}. Spatial bias cancels in this comparison.".format(*time_rates)),
    ]
    setup = simulation_setup(problem, main, "Unit square: 0 ≤ X ≤ 1, 0 ≤ Y ≤ 1. Consistent model units.", [
        "Top edge: T = sin(πX)(1 − exp(−2t)).",
        "Bottom, left and right edges: T = 0.",
        "Displacement is fixed throughout this heat-transfer problem."],
        ["ρcp ∂T/∂t = k ΔT; no volumetric heat source."])
    for plot in plots:
        plot["setup"] = setup
        if record is not None:
            record(plot)
    return {"spatial_error_n24": spatial_errors[-1], "spatial_order": float(min(spatial_rates)),
            "temporal_order": float(min(time_rates)), "continuous_time_error_n24": continuous_error,
            "heat_balance_max_abs_error": max(main["heat_balance_errors"]),
            "thermal_norm_balance_max_abs_error": max(main["thermal_norm_balance_errors"]),
            "boundary_condition_max_abs_error": main["boundary_condition_error"],
            "final_heat_content": content[-1], "increments": len(times) - 1}


if __name__ == "__main__":
    for name, value in check().items():
        print("{} = {:.8e}".format(name, value))
