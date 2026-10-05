"""Independent mesh-level oracles for the optional serial compiled UEL runtime."""

import importlib.util
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _problem_type():
    path = ROOT / "examples/scalar_diffusion_uel/build.py"
    spec = importlib.util.spec_from_file_location("_livebench_heat_model", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.HeatDiffusionProblem


def check(*, record=None):
    from abaqus_ufl.fe import (CompiledAbaqusElement, assemble, boundary_nodes,
                               build_compiled_uel, newton_solve, structured_quad_mesh)

    problem = _problem_type()()
    compiled = build_compiled_uel(problem)
    element = CompiledAbaqusElement(compiled, problem, dt=0.1)

    # A constant-stress affine deformation has zero interior momentum balance.
    # Only boundary displacements and temperatures are prescribed; interior
    # displacement DOFs start away from the exact solution and are solved.
    nodes, elems = structured_quad_mesh(3, 3)
    exact = np.column_stack([0.02 * nodes[:, 0], -0.01 * nodes[:, 1],
                             np.full(len(nodes), 5.0)]).ravel()
    conditions = {int(3 * i + 2): 5.0 for i in range(len(nodes))}
    for i in boundary_nodes(nodes):
        for component in (0, 1):
            conditions[int(3 * i + component)] = exact[3 * i + component]
    initial = np.zeros_like(exact)
    initial[2::3] = 5.0
    solved = newton_solve(nodes, elems, 3, element.element_rk, conditions,
                          U0=initial, U_previous=exact, tol=1e-11)
    affine_error = float(np.max(abs(solved - exact)))
    if affine_error > 1e-9:
        raise AssertionError("compiled affine patch error: {}".format(affine_error))
    print("[PASS] compiled 3x3-mesh affine mechanics patch (interior DOFs solved)")
    if record is not None:
        from tools.livebench_data import comparison_curve, mesh_field
        patch = mesh_field(
            "patch", "Affine displacement patch", nodes, elems,
            np.linalg.norm(solved.reshape(-1, 3)[:, :2], axis=1),
            "Boundary displacements u₁ = 0.02X, u₂ = −0.01Y prescribe an affine patch. "
            "Interior DOFs are solved by the compiled UEL. Colors show nodal displacement magnitude; "
            "the display magnifies deformation by 10.", displacements=solved.reshape(-1, 3)[:, :2])
        profiles = []

    # For T_old=sin(pi*x), all y directions insulated, Q1 consistent mass gives
    # lambda_h = 6(1-cos(pi*h)) / (h^2 (2+cos(pi*h))). Backward Euler therefore
    # gives amplitude 1/(1 + dt*k/rho_cp*lambda_h). This scalar closed form is
    # independent of the generated UEL and the global FE assembly.
    dt = 0.1
    conductivity, capacity = problem._mat.k, problem._mat.rho_cp
    discrete_errors, continuum_errors, free_residuals = [], [], []
    for n in (4, 8, 16):
        nodes, elems = structured_quad_mesh(n, 2)
        previous = np.zeros(len(nodes) * 3)
        previous[2::3] = np.sin(math.pi * nodes[:, 0])
        conditions = {int(3 * i + c): 0.0 for i in range(len(nodes)) for c in (0, 1)}
        for i in np.flatnonzero((nodes[:, 0] == 0.0) | (nodes[:, 0] == 1.0)):
            conditions[int(3 * i + 2)] = 0.0
        element.begin_increment(1, time_begin=0.0, dt=dt)
        solved = newton_solve(nodes, elems, 3, element.element_rk, conditions,
                              U0=previous, U_previous=previous, tol=1e-12)
        h = 1.0 / n
        cosine = math.cos(math.pi * h)
        eigenvalue = 6.0 * (1.0 - cosine) / (h * h * (2.0 + cosine))
        discrete = previous[2::3] / (1.0 + dt * conductivity / capacity * eigenvalue)
        continuum = previous[2::3] / (1.0 + dt * conductivity / capacity * math.pi ** 2)
        discrete_errors.append(float(np.max(abs(solved[2::3] - discrete))))
        continuum_errors.append(float(np.max(abs(solved[2::3] - continuum))))
        residual, _ = assemble(nodes, elems, 3, element.element_rk, solved, solved - previous)
        free = np.array([i for i in range(len(solved)) if i not in conditions])
        free_residuals.append(float(np.max(abs(residual[free]))))
        if not np.allclose(solved[2::3], discrete, rtol=1e-10, atol=1e-11):
            raise AssertionError("compiled transient mode disagrees with discrete closed form")
        # Known-broken control: zero-work output or lost old-field storage cannot
        # satisfy the quantitative amplitude test.
        if np.allclose(np.zeros_like(discrete), discrete, atol=1e-3):
            raise AssertionError("transient mode oracle accepted zero output")
        print("[PASS] compiled {}x2 diffusion mesh vs discrete closed form".format(n))
        if record is not None:
            row = np.flatnonzero(np.isclose(nodes[:, 1], 0.0))
            row = row[np.argsort(nodes[row, 0])]
            profiles.append(comparison_curve(
                "temperature-{}".format(n), "Temperature profile · {} × 2 mesh".format(n),
                "Position X / L", "Temperature / initial peak", nodes[row, 0], solved[2::3][row], discrete[row],
                "Initial temperature sin(πX); T = 0 at x = 0 and 1, insulated y edges. "
                "Displacements are fixed. One backward-Euler step uses Δt = 0.1, k = 0.5, ρcp = 1. "
                "The discrete Q1 consistent-mass eigenmode is exact for this mesh. "
                "The continuum backward-Euler mode measures spatial error, not time-integration error.",
                reference_label="Discrete closed form", rtol=1e-10, atol=1e-11,
                extra_series=[
                    {"label": "Initial temperature", "role": "guide", "values": previous[2::3][row].tolist()},
                    {"label": "Continuum backward Euler", "role": "guide", "values": continuum[row].tolist()},
                ]))
    rate = math.log(continuum_errors[-2] / continuum_errors[-1], 2.0)
    if not 1.9 < rate < 2.1:
        raise AssertionError("unexpected spatial refinement rate: {}".format(rate))
    print("[PASS] spatial convergence to continuum backward-Euler mode (rate {:.3f})".format(rate))
    if record is not None:
        for profile in reversed(profiles):
            record(profile)
        record(mesh_field(
            "temperature-field", "Temperature on the 16 × 2 mesh", nodes, elems, solved[2::3],
            "Nodal temperatures from the same compiled diffusion solve. "
            "Cell colors interpolate the nodal values; mesh lines show the Q1 element boundaries."))
        record({
            "kind": "curve", "id": "convergence", "title": "Spatial mesh convergence",
            "description": "Maximum nodal error against the continuum backward-Euler mode, "
                           "on 4 × 2, 8 × 2, and 16 × 2 meshes. The measured order is {:.4f}. "
                           "The h² guide is anchored at the coarsest mesh; it is not an exact-error reference.".format(rate),
            "x_label": "Mesh spacing h / L", "y_label": "Maximum temperature error",
            "x_scale": "log", "y_scale": "log", "x": [1 / 16, 1 / 8, 1 / 4],
            "series": [
                {"label": "Compiled FE error", "role": "computed", "values": list(reversed(continuum_errors))},
                {"label": "h² guide", "role": "guide", "values": [continuum_errors[0] / 16, continuum_errors[0] / 4, continuum_errors[0]]},
            ],
        })
        record(patch)
    return {
        "affine_patch_max_abs_error": affine_error,
        "transient_mode_max_abs_error": max(discrete_errors),
        "free_residual_max_abs": max(free_residuals),
        "continuum_mode_error_n16": continuum_errors[-1],
        "spatial_convergence_rate": rate,
    }


if __name__ == "__main__":
    for name, value in check().items():
        print("  {} = {:.6e}".format(name, value))
