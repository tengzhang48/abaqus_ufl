"""Nonhomogeneous, finite-strain Ogden shear-block boundary-value problem.

The shipped one-term Ogden material supplies a displacement-only standard
Quad4 verification UEL. Every solved element calls generated Fortran through
the packaged f2py interface; this is a license-free serial verification, not
an Abaqus solver run or a claim that the original UMAT is a UEL.

Omega0=[0,1]², plane strain, zero body force and unit out-of-plane thickness.
Bottom: u=(0,0). Top: u=(0.6*s,0), 0<=s<=1. Side faces: P N=0.
Weak form: integral_Omega0 P:Grad_X(v) dV=0 for v=0 on top/bottom.
The unknowns are Q1 u1/u2; there are no STATEV/SVARS or inertia terms.
Loading parameter s is quasi-static, not physical time. Parameters are
mu=1, alpha=3.5, K=10 in consistent model units. K explicitly differs from
the shipped material's K=100 to avoid a near-incompressible locking example.

Independent real-SVD stress/energy and independent rectangular-Q1 integration
audit equilibrium and reactions. A closed-form isochoric shear traction tests
the nonlinear spectral law. Refinement concerns this mesh BVP and is not a
claim of an exact continuum solution for the clamped/free-corner problem.
"""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import abaqus_ufl as au
from examples.ogden_umat.build import OgdenOneTerm
from tools.boundary_value_common import boundary_region
from tools.livebench_data import comparison_curve, mesh_field
from tools.ogden_bvp_reference import alpha_two_element, element_reference, mesh_reference, shear_traction

SHEAR = 0.6
PROPERTIES = dict(mu=1.0, alpha=3.5, K=10.0)


class OgdenShearProblem(au.WeakForm):
    """Only displacement is solved; P comes from the released Ogden model."""

    material = OgdenOneTerm
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)

    def momentum_equation(self, v, F):
        return self.material.stress_PK1(F)


def conditions(nodes, parameter):
    """End grips constrain two DOFs; side and interior DOFs remain free."""
    bottom = np.flatnonzero(np.isclose(nodes[:, 1], 0))
    top = np.flatnonzero(np.isclose(nodes[:, 1], 1))
    result = {2 * int(i) + j: 0.0 for i in bottom for j in (0, 1)}
    for i in top:
        result[2 * int(i)] = SHEAR * parameter
        result[2 * int(i) + 1] = 0.0
    return result


def solve_history(problem, compiled, subdivisions, parameters, *, boundary_conditions=conditions):
    """Keep only converged nodal solutions and fix accepted history in Newton."""
    from abaqus_ufl.fe import CompiledAbaqusElement, assemble, newton_solve, structured_quad_mesh

    nodes, elements = structured_quad_mesh(subdivisions, subdivisions)
    element = CompiledAbaqusElement(compiled, problem)
    previous = np.zeros(2 * len(nodes))
    solutions, residuals, iterations, free_counts = [previous.copy()], [], [], []
    for index, parameter in enumerate(parameters[1:], 1):
        element.begin_increment(index, time_begin=float(parameters[index - 1]),
                                dt=float(parameter - parameters[index - 1]))
        boundary = boundary_conditions(nodes, parameter)
        history = []
        solved = newton_solve(nodes, elements, 2, element.element_rk, boundary,
                             U0=previous, U_previous=previous, tol=1e-10,
                             record_iteration=history.append)
        residual, _ = assemble(nodes, elements, 2, element.element_rk, solved, solved - previous)
        solutions.append(solved.copy())
        residuals.append(residual)
        iterations.append(history)
        free_counts.append(len(solved) - len(boundary))
        previous = solved
    return {"nodes": nodes, "elements": elements, "times": np.asarray(parameters),
            "solutions": np.asarray(solutions), "residuals": np.asarray(residuals),
            "iterations": iterations, "free_dof_counts": free_counts,
            "properties": {name: float(getattr(problem._mat, name)) for name in PROPERTIES}}


def audit_history(result):
    """Reject inaccurate equations, wrong BCs, frozen DOFs, and invalid output."""
    nodes, elements, parameters = result["nodes"], result["elements"], result["times"]
    solutions, residuals = result["solutions"], result["residuals"]
    if (parameters.ndim != 1 or len(parameters) < 2 or parameters[0] != 0
            or np.any(np.diff(parameters) <= 0) or solutions.shape != (len(parameters), 2 * len(nodes))
            or residuals.shape != (len(parameters) - 1, 2 * len(nodes))
            or not all(np.isfinite(v).all() for v in (nodes, parameters, solutions, residuals))):
        raise AssertionError("invalid Ogden boundary-value solution history")
    subdivisions = len(np.unique(nodes[:, 0])) - 1
    expected_free = 2 * (subdivisions + 1) * (subdivisions - 1)
    if not expected_free or result["free_dof_counts"] != [expected_free] * (len(parameters) - 1):
        raise AssertionError("Ogden BVP must solve free side and interior displacement DOFs")
    if np.max(abs(solutions[0])) > 1e-12:
        raise AssertionError("Ogden initial state must be undeformed")
    top = np.flatnonzero(np.isclose(nodes[:, 1], 1))
    grip_nodes = np.flatnonzero(np.isclose(nodes[:, 1], 0) | np.isclose(nodes[:, 1], 1))
    stress_errors, free_errors, force_errors, moment_errors, work_errors = [], [], [], [], []
    energy, top_force, minimum_J, maximum_stretch, all_samples = [0.0], [0.0], [1.0], [1.0], []
    boundary_error = 0.0
    for parameter, solution, residual in zip(parameters[1:], solutions[1:], residuals):
        boundary = conditions(nodes, parameter)
        free = np.array(sorted(set(range(len(solution))) - set(boundary)), dtype=int)
        boundary_error = max(boundary_error, max(abs(solution[dof] - value) for dof, value in boundary.items()))
        independent, W, samples = mesh_reference(nodes, elements, solution, **PROPERTIES)
        if samples.shape != (len(elements), 4, 3) or not np.isfinite(samples).all():
            raise AssertionError("invalid independent Gauss-point coverage")
        stress_errors.append(float(np.max(abs(residual - independent))))
        free_errors.append(float(np.max(abs(independent[free]))))
        reactions = independent.reshape(-1, 2)[grip_nodes]
        positions = nodes[grip_nodes] + solution.reshape(-1, 2)[grip_nodes]
        force_errors.append(float(np.max(abs(np.sum(reactions, axis=0)))))
        moment_errors.append(float(abs(np.sum(positions[:, 0] * reactions[:, 1]
                                                - positions[:, 1] * reactions[:, 0]))))
        top_force.append(float(np.sum(independent.reshape(-1, 2)[top, 0])))
        # The virtual field v=Y e1 is zero on the bottom and unit on the top.
        # At a converged state its free-node work vanishes, so dE/d(a) is the
        # summed top reaction. Compute this derivative from the independent
        # energy, not from the stress/residual expression being checked.
        virtual = np.column_stack([nodes[:, 1], np.zeros(len(nodes))]).ravel()
        h = 2e-6
        plus_w = mesh_reference(nodes, elements, solution + h * virtual, **PROPERTIES)[1]
        minus_w = mesh_reference(nodes, elements, solution - h * virtual, **PROPERTIES)[1]
        work_errors.append(abs((plus_w - minus_w) / (2 * h) - top_force[-1]))
        energy.append(W)
        minimum_J.append(float(np.min(samples[:, :, 0])))
        maximum_stretch.append(float(np.max(samples[:, :, 1])))
        all_samples.append(samples)
    if boundary_error > 1e-12:
        raise AssertionError("Ogden grip boundary conditions disagree with prescribed loading")
    if max(stress_errors + free_errors) > 2e-9:
        raise AssertionError("compiled Ogden BVP disagrees with independent principal-stretch equilibrium")
    if max(force_errors + moment_errors) > 2e-8:
        raise AssertionError("Ogden grip reactions violate force or current-position moment balance")
    if max(work_errors) > 3e-8:
        raise AssertionError("Ogden grip reactions violate independent global energy derivatives")
    affine_u1 = SHEAR * parameters[-1] * nodes[:, 1]
    nonaffinity = float(np.max(abs(solutions[-1, 0::2] - affine_u1)))
    if (min(minimum_J) <= 0.5 or max(maximum_stretch) < 1.2 or nonaffinity < 0.01
            or top_force[-1] <= 0.1 or np.any(np.diff(energy) <= 0)):
        raise AssertionError("Ogden BVP did not reach the nonhomogeneous finite-strain regime")
    return {"oracle_residual_max_abs_error": max(stress_errors),
            "independent_free_residual_max_abs_error": max(free_errors),
            "force_balance_max_abs_error": max(force_errors),
            "moment_balance_max_abs_error": max(moment_errors),
            "energy_work_max_abs_error": max(work_errors),
            "boundary_condition_max_abs_error": boundary_error,
            "minimum_j": min(minimum_J), "maximum_principal_stretch": max(maximum_stretch),
            "nonaffine_displacement_max_abs": nonaffinity,
            "stored_energy": np.asarray(energy), "top_force": np.asarray(top_force),
            "free_residual_history": np.asarray(free_errors),
            "gauss_samples": np.asarray(all_samples)}


def shear_patch(compiled, problem):
    """Compiled one-element shear response against an independent closed form."""
    from abaqus_ufl.fe import CompiledAbaqusElement, structured_quad_mesh

    nodes, _ = structured_quad_mesh(1, 1)
    element = CompiledAbaqusElement(compiled, problem)
    gamma = np.linspace(0, SHEAR, 13)
    observed = []
    for value in gamma:
        u = np.column_stack([value * nodes[:, 1], np.zeros(len(nodes))]).ravel()
        residual, _ = element.element_rk(0, nodes[[0, 1, 3, 2]], u.reshape(-1, 2)[[0, 1, 3, 2]].ravel(),
                                         u.reshape(-1, 2)[[0, 1, 3, 2]].ravel())
        observed.append(float(np.sum(residual.reshape(4, 2)[[2, 3], 0])))
    exact = shear_traction(gamma, mu=PROPERTIES["mu"], alpha=PROPERTIES["alpha"])
    if not np.allclose(observed, exact, rtol=2e-10, atol=2e-11):
        raise AssertionError("compiled Ogden shear patch disagrees with closed-form principal stretches")
    if float(exact[-1] / (PROPERTIES["mu"] * SHEAR) - 1) < 0.05:
        raise AssertionError("Ogden shear patch did not exercise nonlinear alpha != 2 response")
    return gamma, np.asarray(observed), exact


def audit_element_tangent(compiled, problem):
    """Check a noncoaxial element against independent energy/residual derivatives."""
    from abaqus_ufl.fe import CompiledAbaqusElement

    nodes = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    u = np.array([0., 0., 0.08, -0.04, 0.42, 0.03, 0.27, -0.01])
    element = CompiledAbaqusElement(compiled, problem)
    residual, tangent = element.element_rk(0, nodes, u, u)
    independent, _, _ = element_reference(nodes, u, **PROPERTIES)
    difference = np.empty_like(tangent)
    energy_gradient = np.zeros(8)
    h = 2e-6
    for column in range(8):
        du = np.zeros(8)
        du[column] = h
        plus, plus_w, _ = element_reference(nodes, u + du, **PROPERTIES)
        minus, minus_w, _ = element_reference(nodes, u - du, **PROPERTIES)
        difference[:, column] = (plus - minus) / (2 * h)
        energy_gradient[column] = (plus_w - minus_w) / (2 * h)
    error = float(np.max(abs(tangent - difference)))
    if (np.max(abs(residual - independent)) > 2e-10 or error > 3e-7
            or np.max(abs(residual - energy_gradient)) > 2e-9
            or np.max(abs(tangent - tangent.T)) > 2e-9):
        raise AssertionError("compiled Ogden residual/tangent violates independent energy derivatives")
    return error


def audit_alpha_two_limit(compiled, result):
    """Execute α=2 in the compiled spectral law at all final deformed elements."""
    from abaqus_ufl.fe import CompiledAbaqusElement

    problem = OgdenShearProblem(mu=PROPERTIES["mu"], alpha=2.0, K=PROPERTIES["K"])
    element = CompiledAbaqusElement(compiled, problem)
    error = 0.0
    for index, connectivity in enumerate(result["elements"]):
        coordinates = result["nodes"][connectivity]
        u = result["solutions"][-1].reshape(-1, 2)[connectivity].ravel()
        generated, _ = element.element_rk(index, coordinates, u, u)
        exact = alpha_two_element(coordinates, u, mu=PROPERTIES["mu"], K=PROPERTIES["K"])
        error = max(error, float(np.max(abs(generated - exact))))
    if error > 2e-10:
        raise AssertionError("compiled Ogden alpha=2 limit disagrees with eig-free stress")
    return error


def check(*, record=None):
    from abaqus_ufl.fe import build_compiled_uel

    problem = OgdenShearProblem(**PROPERTIES)
    if not problem.verify():
        raise AssertionError("Ogden declaration tangent verification failed")
    compiled = build_compiled_uel(problem)
    tangent_error = audit_element_tangent(compiled, problem)
    gamma, observed_patch, exact_patch = shear_patch(compiled, problem)
    parameters = np.linspace(0, 1, 13)
    results = [solve_history(problem, compiled, n, parameters) for n in (4, 8, 16)]
    audits = [audit_history(result) for result in results]
    forces = np.array([audit["top_force"][-1] for audit in audits])
    energies = np.array([audit["stored_energy"][-1] for audit in audits])
    changes = abs(np.diff(forces))
    if (changes[1] >= changes[0] / 1.5 or changes[1] / abs(forces[-1]) > 0.03
            or np.any(np.diff(energies) >= 0)):
        raise AssertionError("Ogden BVP did not improve with 4/8/16 mesh refinement")
    main, audit = results[-1], audits[-1]
    alpha_two_error = audit_alpha_two_limit(compiled, main)
    nodes, elements = main["nodes"], main["elements"]
    displacement = main["solutions"].reshape(len(parameters), len(nodes), 2)
    values = np.linalg.norm(displacement, axis=2)
    boundary = [boundary_region(nodes, 1, 0, "Clamped: u₁ = u₂ = 0", "displacement"),
                boundary_region(nodes, 1, 1, "Moving grip: u₁ = 0.6s, u₂ = 0", "displacement"),
                boundary_region(nodes, 0, 0, "Traction-free side: P N = 0", "displacement", "natural"),
                boundary_region(nodes, 0, 1, "Traction-free side: P N = 0", "displacement", "natural")]
    deformation = mesh_field("deformation", "Ogden block under finite shear", nodes, elements, values[-1],
        "Bottom grip is clamped and the top grip moves by 0.6. Side faces are traction-free. "
        "All 510 free side/interior displacement DOFs are solved across 256 Quad4 elements. "
        "Colors are authoritative nodal displacement magnitudes; deformation is displayed at ×1.",
        displacements=displacement[-1])
    deformation["value_label"] = "Displacement magnitude (model units)"
    deformation["deformation_scale"] = 1.0
    deformation["frames"] = [{"time": float(s), "values": v.tolist(), "displacements": u.tolist()}
                              for s, v, u in zip(parameters, values, displacement)]
    deformation["boundaries"] = boundary
    plots = [deformation,
        comparison_curve("stretch", "Independent finite-shear material check", "Homogeneous shear γ", "Nominal shear P₁₂",
            gamma, observed_patch, exact_patch,
            "A separate one-element affine shear patch checks the generated spectral law. "
            "The principal stretches (sqrt(γ²+4) ± γ)/2 give an independent closed-form traction. "
            "The block simulation above has free sides and a nonhomogeneous displacement field.",
            reference_label="Closed-form principal stretches", rtol=2e-10, atol=2e-11),
        {"kind": "curve", "id": "load-history", "title": "Reaction of the moving grip",
            "description": "Total horizontal grip reaction per unit out-of-plane thickness. "
                "The second curve independently assembles SVD principal-stretch stresses at every Gauss point; "
                "it audits the discretized equations, rather than providing an exact continuum solution.",
            "x_label": "Loading parameter s", "y_label": "Horizontal reaction force", "x": parameters.tolist(),
            "x_scale": "linear", "y_scale": "linear", "series": [
                {"label": "Compiled FE", "role": "computed", "values": [0.0] + [
                    float(np.sum(r.reshape(-1, 2)[np.isclose(nodes[:, 1], 1), 0])) for r in main["residuals"]]},
                {"label": "Independent SVD stress assembly", "role": "reference", "values": audit["top_force"].tolist()}]},
        {"kind": "curve", "id": "refinement", "title": "Grip reaction under mesh refinement",
            "description": "4×4, 8×8 and 16×16 meshes solve the same nonhomogeneous BVP. "
                "The final refinement must reduce the reaction change by at least 1.5 and change it by less than 3%. "
                "This is a refinement check, not an exact continuum error estimate.",
            "x_label": "Elements per edge", "y_label": "Final horizontal grip reaction", "x": [4, 8, 16],
            "x_scale": "linear", "y_scale": "linear", "series": [
                {"label": "Compiled FE", "role": "computed", "values": forces.tolist()}]},
        {"kind": "curve", "id": "equilibrium", "title": "Independent free-DOF equilibrium",
            "description": "The maximum free residual is recomputed with independent real-SVD stress and "
                "rectangular-Q1 integration at each accepted step. Grip force and current-position moment "
                "balances are also checked. No unconverged trial is recorded.",
            "x_label": "Loading parameter s", "y_label": "Independent maximum free residual", "x": parameters[1:].tolist(),
            "x_scale": "linear", "y_scale": "linear", "series": [
                {"label": "Independent audit of compiled FE", "role": "computed", "values": [
                    float(value) for value in audit["free_residual_history"]]},
                {"label": "Acceptance tolerance", "role": "guide", "values": [2e-9] * 12}]},
    ]
    setup = {"domain": "Unit square, plane strain, unit out-of-plane thickness. Consistent model units.",
             "equations": ["Div_X P = 0; finite-strain one-term compressible Ogden elasticity.",
                           "Loading parameter s is quasi-static; no inertia or time-dependent material state."],
             "initial_condition": "u = 0; F = I throughout the mesh.",
             "boundary_conditions": ["Bottom: u₁ = u₂ = 0.", "Top: u₁ = 0.6s and u₂ = 0.",
                                     "Left/right: P N = 0 (traction-free)."],
             "properties": [{"name": name, "value": value} for name, value in PROPERTIES.items()],
             "time_step": float(parameters[1]), "steps": 12, "nodes": len(nodes), "elements": len(elements),
             "free_temperature_dofs": 0, "free_displacement_dofs": main["free_dof_counts"][0]}
    for plot in plots:
        plot["setup"] = setup
        if record is not None:
            record(plot)
    print("[PASS] Ogden verification UEL: 12 nonhomogeneous finite-strain mesh increments")
    print("[PASS] independent principal-stretch stress, energy, shear law, force and moment balance")
    print("[PASS] 4/8/16 refinement; final reaction change {:.3%}".format(changes[-1] / abs(forces[-1])))
    metrics = {name: value for name, value in audit.items() if np.isscalar(value)}
    metrics.update({"independent_tangent_max_abs_error": tangent_error,
                    "alpha_two_limit_max_abs_error": alpha_two_error,
                    "shear_patch_max_abs_error": float(np.max(abs(observed_patch - exact_patch))),
                    "refinement_final_reaction_relative_change": float(changes[-1] / abs(forces[-1])),
                    "final_grip_reaction": float(forces[-1]), "increments": 12,
                    "max_newton_iterations": max(h[-1]["iteration"] for h in main["iterations"])})
    return metrics


if __name__ == "__main__":
    for name, value in check().items():
        print("{} = {:.8e}".format(name, value))
