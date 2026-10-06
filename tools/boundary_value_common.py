"""Compiled boundary-value solves and independently derived heat oracles.

The existing public thermo-mechanical Quad4 supplies all element residuals and
tangents. The Fourier references below solve the continuum PDE directly; they
do not call the declaration, generated element, or FE assembly.
"""

import importlib.util
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def heat_problem(**properties):
    spec = importlib.util.spec_from_file_location(
        "_boundary_value_heat_model", str(ROOT / "examples/scalar_diffusion_uel/build.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.HeatDiffusionProblem(**properties)


def ramp(time, rate=4.0):
    return -math.expm1(-rate * time)


def _modes(times, coefficients, decay, rate, backward_euler):
    times = np.asarray(times, dtype=float)
    if times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("reference times must start at zero and increase")
    result = np.zeros((len(times), len(decay)))
    if backward_euler:
        for j in range(1, len(times)):
            increment = ramp(times[j], rate) - ramp(times[j - 1], rate)
            result[j] = (result[j - 1] - coefficients * increment) / (
                1 + (times[j] - times[j - 1]) * decay)
    else:
        difference = decay - rate
        resonant = np.isclose(difference, 0.0, rtol=0, atol=1e-12)
        t = times[:, None]
        ordinary = ~resonant
        result[:, ordinary] = (-coefficients[ordinary] * rate
            * (np.exp(-rate * t) - np.exp(-t * decay[ordinary])) / difference[ordinary])
        result[:, resonant] = -coefficients[resonant] * rate * t * np.exp(-rate * t)
    return result


def plate_reference(nodes, times, diffusivity=0.5, rate=4.0, *, backward_euler=False):
    """Unit square: top sin(pi X) g(t), other sides zero, initially zero.

    The harmonic lifting is sin(pi X) sinh(pi Y)/sinh(pi). Its sine-series
    coefficients in Y are 2*(-1)^(n+1)*n/[pi*(1+n^2)]. A continuum BE
    recurrence separates spatial error from time-discretization error.
    """
    nodes = np.asarray(nodes)
    n = np.arange(1, 257)
    coefficients = 2 * (-1.0) ** (n + 1) * n / (math.pi * (1 + n * n))
    modes = _modes(times, coefficients, diffusivity * math.pi ** 2 * (1 + n * n),
                   rate, backward_euler)
    lifting = np.sinh(math.pi * nodes[:, 1]) / math.sinh(math.pi)
    return np.sin(math.pi * nodes[:, 0])[None, :] * (
        np.array([ramp(t, rate) for t in times])[:, None] * lifting
        + modes @ np.sin(math.pi * n[:, None] * nodes[:, 1]))


def plate_heat_reference(times, capacity=1.0, diffusivity=0.5, rate=4.0, *, backward_euler=False):
    """Exact spatial integral of the independent plate sine series."""
    n = np.arange(1, 257)
    coefficients = 2 * (-1.0) ** (n + 1) * n / (math.pi * (1 + n * n))
    modes = _modes(times, coefficients, diffusivity * math.pi ** 2 * (1 + n * n),
                   rate, backward_euler)
    lifting_integral = (math.cosh(math.pi) - 1) / (math.pi * math.sinh(math.pi))
    return capacity * 2 / math.pi * (
        np.array([ramp(t, rate) for t in times]) * lifting_integral
        + modes @ ((1 - (-1.0) ** n) / (math.pi * n)))


def strip_reference(nodes, times, diffusivity=0.5, rate=4.0, height=1.0,
                    temperature=1.0, *, backward_euler=False):
    """Insulated ends, bottom -T0*g(t), top +T0*g(t), initially zero."""
    n = np.arange(1, 257)
    coefficients = -2 * (1 + (-1.0) ** n) / (n * math.pi)
    modes = _modes(times, coefficients, diffusivity * (n * math.pi / height) ** 2,
                   rate, backward_euler)
    y = np.asarray(nodes)[:, 1] / height
    return temperature * (np.array([ramp(t, rate) for t in times])[:, None] * (2 * y - 1)
                          + modes @ np.sin(math.pi * n[:, None] * y))


def beam_tip_reference(times, *, shear=1.0, volumetric=2.0, alpha=1e-3,
                       diffusivity=0.5, rate=4.0, height=1.0, length=8.0,
                       temperature=1.0, backward_euler=True):
    """Small-strain slender-beam guide, including the transient thermal moment.

    The actual material linearizes to 2G eps+K tr(eps) I-K alpha T I.
    Eliminating transverse stress in plane strain gives beta=K alpha/[2(G+K)].
    Zero bending moment gives curvature=-12 beta/H^3 * integral(y_center*T).
    Full end clamping introduces an end effect; this is an approximation,
    not an exact continuum solution for the finite-strain clamped strip.
    """
    n = np.arange(1, 257)
    coefficients = -2 * (1 + (-1.0) ** n) / (n * math.pi)
    modes = _modes(times, coefficients, diffusivity * (n * math.pi / height) ** 2,
                   rate, backward_euler)
    moment = temperature * height ** 2 * (
        np.array([ramp(t, rate) for t in times]) / 6
        - modes @ ((1 + (-1.0) ** n) / (2 * n * math.pi)))
    beta = volumetric * alpha / (2 * (shear + volumetric))
    return -6 * beta * length ** 2 / height ** 3 * moment


def heat_content(nodes, elems, temperature, capacity=1.0):
    """Independent integral of the rectangular Q1 temperature field."""
    corners = nodes[elems]
    following = np.roll(corners, -1, axis=1)
    areas = 0.5 * np.sum(corners[:, :, 0] * following[:, :, 1]
                         - following[:, :, 0] * corners[:, :, 1], axis=1)
    if np.any(areas <= 0):
        raise AssertionError("nonpositive reference element area")
    return float(capacity * np.sum(areas * np.mean(temperature[elems], axis=1)))


def thermal_norms(nodes, elems, temperature):
    """Independent rectangular Q1 quadrature of T^2 and |Grad T|^2."""
    squared, gradient_squared = 0.0, 0.0
    for coordinates, values in zip(nodes[elems], temperature[elems]):
        dx = coordinates[1, 0] - coordinates[0, 0]
        dy = coordinates[3, 1] - coordinates[0, 1]
        for xi in (-1 / math.sqrt(3), 1 / math.sqrt(3)):
            for eta in (-1 / math.sqrt(3), 1 / math.sqrt(3)):
                weights = np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                                    (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)]) / 4
                gradient = np.array([
                    [-(1 - eta), 1 - eta, 1 + eta, -(1 + eta)],
                    [-(1 - xi), -(1 + xi), 1 + xi, 1 - xi],
                ]) / np.array([[2 * dx], [2 * dy]])
                squared += float((weights @ values) ** 2 * dx * dy / 4)
                gradient_squared += float(np.sum((gradient @ values) ** 2) * dx * dy / 4)
    return squared, gradient_squared


def solve_history(problem, compiled, nx, ny, times, boundary_conditions,
                  *, length=1.0, height=1.0):
    """Solve every free nodal DOF; preserve accepted history during Newton."""
    from abaqus_ufl.fe import CompiledAbaqusElement, assemble, newton_solve, structured_quad_mesh

    nodes, elems = structured_quad_mesh(nx, ny, length, height)
    element = CompiledAbaqusElement(compiled, problem)
    previous = np.zeros(3 * len(nodes))
    solutions, residuals, iterations, balances = [previous.copy()], [], [], []
    dissipation, norm_balance_errors = [], []
    free_dof_counts = []
    condition_error = 0.0
    for j, time in enumerate(times[1:], 1):
        dt = time - times[j - 1]
        conditions = boundary_conditions(nodes, time)
        element.begin_increment(j, time_begin=float(times[j - 1]), dt=float(dt))
        history = []
        solved = newton_solve(nodes, elems, 3, element.element_rk, conditions,
                             U0=previous, U_previous=previous, tol=1e-10,
                             record_iteration=history.append)
        residual, _ = assemble(nodes, elems, 3, element.element_rk, solved, solved - previous)
        free = np.array(sorted(set(range(len(solved))) - set(conditions)), dtype=int)
        free_dof_counts.append({"temperature": int(np.sum(free % 3 == 2)),
                                "displacement": int(np.sum(free % 3 != 2))})
        if not len(free) or np.max(abs(residual[free])) > 1e-9:
            raise AssertionError("boundary-value free residual did not converge")
        condition_error = max(condition_error, max(abs(solved[i] - v) for i, v in conditions.items()))
        thermal_dofs = [i for i in conditions if i % 3 == 2]
        rate = (heat_content(nodes, elems, solved[2::3], problem._mat.rho_cp)
                - heat_content(nodes, elems, previous[2::3], problem._mat.rho_cp)) / dt
        balances.append(abs(rate - np.sum(residual[thermal_dofs])))
        squared, gradient_squared = thermal_norms(nodes, elems, solved[2::3])
        old_squared, _ = thermal_norms(nodes, elems, previous[2::3])
        increment_squared, _ = thermal_norms(nodes, elems, (solved - previous)[2::3])
        diffusion = problem._mat.k * gradient_squared
        norm_rate = problem._mat.rho_cp * (squared - old_squared + increment_squared) / (2 * dt)
        boundary_work = sum(solved[i] * residual[i] for i in thermal_dofs)
        norm_balance_errors.append(abs(norm_rate + diffusion - boundary_work))
        dissipation.append(diffusion)
        solutions.append(solved.copy())
        residuals.append(residual)
        iterations.append(history)
        previous = solved
    if condition_error > 1e-12 or max(balances + norm_balance_errors) > 1e-9:
        raise AssertionError("boundary conditions or supplied-heat balance failed")
    return {"nodes": nodes, "elements": elems, "times": np.asarray(times),
            "solutions": np.asarray(solutions), "residuals": residuals,
            "iterations": iterations, "heat_balance_errors": balances,
            "boundary_condition_error": condition_error,
            "dissipation": dissipation, "thermal_norm_balance_errors": norm_balance_errors,
            "free_dof_counts": free_dof_counts}


def validate_history(result):
    nodes, times, solutions = result["nodes"], result["times"], result["solutions"]
    if (nodes.ndim != 2 or nodes.shape[1] != 2 or len(times) < 2
            or solutions.shape != (len(times), 3 * len(nodes))
            or not np.isfinite(nodes).all() or not np.isfinite(times).all()
            or not np.isfinite(solutions).all() or times[0] != 0
            or np.any(np.diff(times) <= 0)):
        raise AssertionError("invalid boundary-value solution history")
    dissipation = np.asarray(result["dissipation"])
    if dissipation.shape != (len(times) - 1,) or not np.isfinite(dissipation).all() or np.any(dissipation < 0):
        raise AssertionError("invalid thermal dissipation history")


def history_mesh(identifier, title, result, description, *, field="temperature", boundaries=()):
    """Export authoritative nodal history with a final-frame static fallback."""
    from tools.livebench_data import mesh_field
    solutions = result["solutions"].reshape(len(result["times"]), -1, 3)
    displacement = solutions[:, :, :2]
    values = solutions[:, :, 2] if field == "temperature" else np.linalg.norm(displacement, axis=2)
    plot = mesh_field(identifier, title, result["nodes"], result["elements"], values[-1],
                      description, displacements=displacement[-1] if field != "temperature" else None)
    plot["value_label"] = "Temperature rise (model units)" if field == "temperature" else "Displacement magnitude (model units)"
    plot["frames"] = [{"time": float(time), "values": v.tolist(),
                       **({"displacements": u.tolist()} if field != "temperature" else {})}
                      for time, v, u in zip(result["times"], values, displacement)]
    plot["boundaries"] = list(boundaries)
    return plot


def boundary_region(nodes, axis, value, label, field, kind="prescribed"):
    indices = np.flatnonzero(np.isclose(nodes[:, axis], value))
    indices = indices[np.argsort(nodes[indices, 1 - axis])]
    return {"label": label, "field": field, "kind": kind, "nodes": indices.tolist()}


def simulation_setup(problem, result, domain, boundary_conditions, equations):
    return {"domain": domain, "equations": equations, "initial_condition": "T = 0 and u = 0 throughout the mesh.",
            "boundary_conditions": boundary_conditions,
            "properties": [{"name": name, "value": float(getattr(problem._mat, name))}
                           for name in problem._mat.props_names],
            "time_step": float(result["times"][1] - result["times"][0]),
            "steps": len(result["times"]) - 1,
            "nodes": len(result["nodes"]), "elements": len(result["elements"]),
            "free_temperature_dofs": result["free_dof_counts"][0]["temperature"],
            "free_displacement_dofs": result["free_dof_counts"][0]["displacement"]}
