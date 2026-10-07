"""Compiled J2 plasticity on a non-homogeneous load/unload boundary-value mesh.

The shipped SmallStrainJ2 declaration generates a local UMAT. An explicitly
example-owned full-integration Quad4 host supplies its weak equilibrium,
displacement boundary conditions and accepted Gauss-point history. This is a
license-free plane-strain verification simulation, not an Abaqus execution or
a general stateful abaqus_ufl.fe API. All quantities use consistent model units.
"""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.livebench_data import comparison_curve
from tools.plasticity_bvp_common import (
    PROPERTIES, build_compiled_j2, mesh_record, point_closed_form_gate,
    solve_history, top_displacement, top_reactions, validate_history,
)


def check(*, record=None):
    compiled = build_compiled_j2()
    point_error = point_closed_form_gate(compiled)
    meshes = (4, 8, 12)
    results = [solve_history(compiled, size) for size in meshes]
    diagnostics = [validate_history(result) for result in results]
    main = results[-1]
    peak = int(np.argmax([top_displacement(time) for time in main["times"]]))
    peak_force = np.array([top_reactions(result)[peak] for result in results])
    changes = abs(np.diff(peak_force))
    last_relative_change = float(changes[-1] / abs(peak_force[-1]))
    if not changes[-1] < .5 * changes[0] or last_relative_change > .025:
        raise AssertionError("J2 peak shear reaction did not settle with mesh refinement")
    print("[PASS] regenerated UMAT byte parity and compiled normal/shear closed forms")
    print("[PASS] non-homogeneous J2 mesh: elastic, plastic, elastic unloading and reverse yield")
    print("[PASS] independent radial return at every Gauss point, free equilibrium and force/moment balance")
    print("[PASS] rejected trials leave state unchanged; accepted endpoint work and dissipation balance")
    print("[PASS] 4/8/12 mesh refinement: last peak-reaction change {:.2%}".format(last_relative_change))
    weights = main["weights"]
    mean_ep = np.sum(main["ep"] * weights[None, :, :], axis=(1, 2)) / np.sum(weights)
    reference_ep = np.sum(main["reference_ep"] * weights[None, :, :], axis=(1, 2)) / np.sum(weights)
    times = main["times"]
    free_residual, oracle_residual = [], []
    for history, computed, reference in zip(main["iterations"], main["residuals"], main["reference_residuals"]):
        free_residual.append(history[-1]["relative_residual"])
        oracle_residual.append(float(np.max(abs(computed - reference)) / max(1., np.max(abs(computed)))))
    plots = [
        mesh_record("plastic-strain", "Irreversible plastic zones through the cycle", main, plastic=True),
        mesh_record("deformation", "Displacement during shear loading and unloading", main),
        comparison_curve("load-history", "Reaction-force hysteresis", "Prescribed top displacement u₁", "Top shear reaction (model units)",
            [top_displacement(time) for time in times], top_reactions(main), top_reactions(main, reference=True),
            "The displacement-controlled cycle loads from 0 to 0.02 and returns to zero. "
            "Each top reaction is independently reassembled using the scalar/Voigt radial-return consistency solution "
            "and an independent rectangular-Q1 strain/force calculation, retaining accepted history. "
            "Late unloading enters reverse plasticity; this is rate-independent pseudo-time, with no inertia.",
            reference_label="Independent radial-return assembly", rtol=0, atol=1e-10),
        comparison_curve("state-history", "Accumulated plastic strain", "Load-cycle pseudo-time", "Equivalent plastic strain",
            times, mean_ep, reference_ep,
            "The reference-area-weighted mean uses raw Gauss-point STATEV(1), not the nodal visualization projection. "
            "The independent oracle evolves its own accepted state. Plastic strain remains fixed during elastic unloading "
            "and increases again when the reverse-shear yield surface is reached.",
            reference_label="Independent accepted-state oracle", rtol=0, atol=1e-10),
        {"kind": "curve", "id": "equilibrium", "title": "Accepted equilibrium and independent force agreement",
            "description": "Free equilibrium is checked at every accepted increment. The second curve measures the maximum "
                           "difference between compiled and independently assembled nodal force vectors, including boundary reactions. "
                           "Both use the larger of one and the maximum nodal force as their scale. "
                           "The exact endpoint-work identity is separately gated at 10⁻¹⁰ and includes σzz, stored hardening, "
                           "physical plastic dissipation, and nonnegative endpoint numerical dissipation.",
            "x_label": "Load-cycle pseudo-time", "y_label": "Normalized force error", "x": times[1:].tolist(),
            "x_scale": "linear", "y_scale": "linear", "series": [
                {"label": "Free equilibrium residual", "role": "computed", "values": free_residual},
                {"label": "Independent nodal-force difference", "role": "computed", "values": oracle_residual},
                {"label": "Acceptance tolerance", "role": "guide", "values": [1e-10] * (len(times) - 1)}]},
        {"kind": "curve", "id": "refinement", "title": "Peak shear reaction with mesh refinement",
            "description": "Uniform full-integration Quad4 meshes contain 16, 64 and 144 elements. "
                           "The same twenty accepted load increments are solved on each mesh. The final mesh change in "
                           "peak shear reaction is below 2.5% and is less than half the preceding change. "
                           "This is a measured refinement study, not a claimed exact continuum solution or a convergence order.",
            "x_label": "Elements per edge", "y_label": "Top reaction at peak displacement", "x": list(meshes),
            "x_scale": "linear", "y_scale": "linear", "series": [
                {"label": "Compiled UMAT / Quad4 host", "role": "computed", "values": peak_force.tolist()}]},
    ]
    # This response is a chronological hysteresis path, so its abscissa is
    # deliberately revisited while loading reverses. Other curves retain the
    # report schema's default strictly increasing abscissa.
    plots[2]["x_order"] = "recorded"
    setup = {
        "domain": "Unit square, plane strain (εzz = γxz = γyz = 0); consistent model units.",
        "equations": ["Div σ = 0, ε = sym Grad u; quasi-static small-strain equilibrium.",
                      "q ≤ σy + H εp, associated J2 radial return, and irreversible STATEV(1) = εp.",
                      "∫Ω σ : ε(v) dΩ = 0 for tests vanishing on the prescribed top and bottom."],
        "initial_condition": "u = 0, all six stress components = 0, and STATEV(1) = 0 at every Gauss point.",
        "boundary_conditions": ["Bottom: u₁ = u₂ = 0.",
                                "Top: u₁ = 0.02t for 0 ≤ t ≤ 1, then 0.02(2−t) for 1 < t ≤ 2; u₂ = 0.",
                                "Left/right sides: zero traction. No body force."],
        "properties": [{"name": key, "value": value} for key, value in PROPERTIES.items()],
        "time_step": float(times[1] - times[0]), "steps": len(times) - 1,
        "nodes": len(main["nodes"]), "elements": len(main["elements"]),
        "free_temperature_dofs": 0, "free_displacement_dofs": main["free_dof_counts"][0],
    }
    for plot in plots:
        plot["setup"] = setup
        if record is not None:
            record(plot)
    return {
        **diagnostics[-1], "point_closed_form_max_abs_error": point_error,
        "stress_oracle_max_abs_error": main["stress_oracle_error"],
        "state_oracle_max_abs_error": main["state_oracle_error"],
        "force_oracle_max_abs_error": main["force_oracle_error"],
        "endpoint_work_balance_max_abs_error": float(np.max(main["endpoint_work_balance_errors"])),
        "yield_consistency_max_abs_error": float(np.max(main["yield_consistency_errors"])),
        "boundary_condition_max_abs_error": main["boundary_condition_error"],
        "peak_shear_reaction": float(peak_force[-1]), "last_mesh_peak_reaction_relative_change": last_relative_change,
        "final_mean_equivalent_plastic_strain": float(mean_ep[-1]),
        "peak_max_equivalent_plastic_strain": float(np.max(main["ep"][peak])),
        "physical_plastic_dissipation": float(np.sum(main["dissipation"])),
        "max_newton_iterations": max(history[-1]["iteration"] for history in main["iterations"]),
        "increments": len(times) - 1, "nodes": len(main["nodes"]), "elements": len(main["elements"]),
        "free_displacement_dofs": main["free_dof_counts"][0], "gauss_points": 4 * len(main["elements"]),
    }


if __name__ == "__main__":
    for name, value in check().items():
        print("{} = {:.8e}".format(name, value))
