"""Example-owned plane-strain Quad4 host for the public generated J2 UMAT.

This is a verification host, not a new material-state facility in abaqus_ufl.fe
and not an Abaqus job. Two nodal displacement DOFs use 2x2 Gauss integration.
The UMAT receives the full three-dimensional tension-positive Abaqus vectors
[xx, yy, zz, xy, xz, yz]; strain shear slots are engineering shear. Plane strain
sets eps_zz=gamma_xz=gamma_yz=0 while retaining the computed sigma_zz.

Each trial starts from copies of ACCEPTED stress, strain and STATEV(1), the
equivalent plastic strain. Only commit_equilibrium may accept a new state.
The independently written radial-return oracle uses no material/generator
helpers; closed forms additionally check the elastic limit and pure shear.
"""

import importlib.util
from importlib.machinery import EXTENSION_SUFFIXES
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

PROPERTIES = dict(G=10.0, lam=20.0, sigma_y=0.1, H=5.0)
GAUSS = tuple((xi, eta) for xi in (-1 / math.sqrt(3), 1 / math.sqrt(3))
              for eta in (-1 / math.sqrt(3), 1 / math.sqrt(3)))


def build_compiled_j2(workdir=None):
    """Regenerate the shipped law, require byte parity, and compile its UMAT."""
    import abaqus_ufl as au

    if shutil.which("gfortran") is None:
        raise RuntimeError("gfortran is required for the J2 mesh check")
    source_dir = ROOT / "examples/small_strain_j2_umat"
    spec = importlib.util.spec_from_file_location("_plasticity_bvp_model", source_dir / "build.py")
    declaration = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(declaration)
    model = declaration.SmallStrainJ2()
    if not model.verify(verbose=False):
        raise AssertionError("J2 Python tangent verification failed")
    temporary = tempfile.TemporaryDirectory(prefix="abaqus_ufl_j2_bvp_") if workdir is None else None
    build = Path(temporary.name if temporary else workdir).resolve()
    build.mkdir(parents=True, exist_ok=True)
    generated = build / "small_strain_j2.for"
    au.generate_small_strain_umat(model, str(generated))
    if generated.read_bytes() != (source_dir / generated.name).read_bytes():
        raise AssertionError("committed J2 Fortran is stale")
    module_name = "_abaqus_ufl_j2_bvp_" + uuid.uuid4().hex
    command = [sys.executable, "-m", "numpy.f2py", "-c",
               str(source_dir / "f2py/drive_umat.f90"), str(generated),
               "-m", module_name, "only:", "drive_umat", ":"]
    if np.lib.NumpyVersion(np.__version__) >= "1.26.0":
        command.extend(["--backend", "meson"])
    process = subprocess.run(command, cwd=str(build), capture_output=True, text=True,
                             encoding="utf-8", errors="replace",
                             env={**os.environ, "PYTHONUTF8": "1"}, check=False)
    log = process.stdout + "\n" + process.stderr
    (build / "compile.log").write_text(log, encoding="utf-8")
    if process.returncode:
        raise RuntimeError("f2py J2 build failed:\n" + log)
    extension = next((build / (module_name + suffix) for suffix in EXTENSION_SUFFIXES
                      if (build / (module_name + suffix)).is_file()), None)
    if extension is None:
        raise RuntimeError("f2py did not produce the J2 extension")
    spec = importlib.util.spec_from_file_location(module_name, extension)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module._plasticity_bvp_build_directory = temporary
    return module


def call_umat(module, stress_old, ep_old, strain_old, strain_increment,
              time_begin, dt, properties=PROPERTIES):
    """Copy old state before an in/out Fortran call; reject cutbacks or NaNs."""
    stress, state, tangent, pnewdt = module.drive_umat(
        np.asarray(stress_old, dtype=float).copy(), np.array([ep_old], dtype=float),
        np.eye(3), np.eye(3), np.asarray(strain_increment, dtype=float),
        np.asarray(strain_old, dtype=float), np.array(list(properties.values()), dtype=float),
        np.array([time_begin, time_begin]), float(dt), 0., 0., np.array([0.]),
        np.array([0.]), 1., np.zeros(3), np.eye(3), 3, 3, b"J2_MESH" + b" " * 73)
    if (not all(np.isfinite(x).all() for x in (stress, state, tangent))
            or not np.isfinite(pnewdt) or pnewdt < 1):
        raise RuntimeError("compiled J2 UMAT returned nonfinite output or requested cutback")
    return np.asarray(stress), float(state[0]), np.asarray(tangent)


def radial_return_reference(stress_old, ep_old, strain_increment, properties=PROPERTIES):
    """Independent scalar/Voigt J2 return in the host's tensile convention.

    q=sqrt(3/2*(sxx²+syy²+szz²+2*(txy²+txz²+tyz²))).
    The consistency solution is dep=max(q-sigy-H*ep_old,0)/(3G+H).
    No shipped material method, generated tangent or tensor helper is used.
    """
    G, lam, sigy, hardening = (properties[k] for k in ("G", "lam", "sigma_y", "H"))
    increment = np.asarray(strain_increment)
    trial = np.asarray(stress_old).copy()
    trial[:3] += 2 * G * increment[:3] + lam * np.sum(increment[:3])
    trial[3:] += G * increment[3:]
    mean = np.sum(trial[:3]) / 3
    deviator = trial.copy()
    deviator[:3] -= mean
    q = math.sqrt(1.5 * (np.sum(deviator[:3] ** 2) + 2 * np.sum(deviator[3:] ** 2)))
    excess = q - sigy - hardening * ep_old
    dep = max(excess, 0.) / (3 * G + hardening)
    if dep > 0:
        trial -= (3 * G * dep / q) * deviator
    return trial, float(ep_old + dep)


def elastic_strain(stress, properties=PROPERTIES):
    """Independent compliance; engineering shear is stress/G."""
    stress = np.asarray(stress)
    G, lam = properties["G"], properties["lam"]
    value = stress / (2 * G)
    value[..., :3] -= (lam * np.sum(stress[..., :3], axis=-1)
                      / (2 * G * (2 * G + 3 * lam)))[..., None]
    value[..., 3:] = stress[..., 3:] / G
    return value


def point_closed_form_gate(module):
    """Hand-derived elastic normal and loading/unloading pure-shear checks."""
    G, lam, sigy, hardening = (PROPERTIES[k] for k in ("G", "lam", "sigma_y", "H"))
    stress, ep, tangent = call_umat(module, np.zeros(6), 0., np.zeros(6),
                                   np.array([1e-4, 0., 0., 0., 0., 0.]), 0., .1)
    normal_expected = np.array([(2 * G + lam) * 1e-4, lam * 1e-4,
                                lam * 1e-4, 0., 0., 0.])
    errors = [float(np.max(abs(stress - normal_expected))), abs(ep), abs(tangent[3, 3] - G)]
    stress, ep, strain = np.zeros(6), 0., np.zeros(6)
    yield_gamma = sigy / (math.sqrt(3) * G)
    # The first increment remains elastic; the last loading state is plastic.
    gamma_path = [yield_gamma / 2, yield_gamma, .01, .02, .018]
    unloading_ep = None
    for i, gamma in enumerate(gamma_path):
        new_strain = np.array([0., 0., 0., gamma, 0., 0.])
        stress, ep, tangent = call_umat(module, stress, ep, strain, new_strain - strain, i * .1, .1)
        if i < 4:
            ep_expected = max(math.sqrt(3) * G * gamma - sigy, 0.) / (3 * G + hardening)
            tau_expected = G * gamma if ep_expected == 0 else (sigy + hardening * ep_expected) / math.sqrt(3)
            errors.extend([abs(stress[3] - tau_expected), abs(ep - ep_expected)])
            if i == 3:
                errors.append(abs(tangent[3, 3] - G * hardening / (3 * G + hardening)))
                unloading_ep = ep_expected
                unloading_stress = tau_expected - .002 * G
        else:
            errors.extend([abs(stress[3] - unloading_stress), abs(ep - unloading_ep),
                           abs(tangent[3, 3] - G)])
        strain = new_strain
    if max(errors) > 1e-10:
        raise AssertionError("compiled UMAT disagrees with elastic/pure-shear closed forms")
    return max(errors)


def quad4_kinematics(coordinates):
    """B maps ux,uy to six-component engineering strain, 2x2 full integration."""
    coordinates = np.asarray(coordinates)
    sx, sy = np.array([-1, 1, 1, -1]), np.array([-1, -1, 1, 1])
    matrices, weights, points = [], [], []
    for xi, eta in GAUSS:
        shape = (1 + sx * xi) * (1 + sy * eta) / 4
        natural_gradient = np.column_stack((sx * (1 + sy * eta), sy * (1 + sx * xi))) / 4
        jacobian = coordinates.T @ natural_gradient
        weight = float(np.linalg.det(jacobian))
        if weight <= 0:
            raise ValueError("J2 host requires positively oriented Quad4 elements")
        gradient = natural_gradient @ np.linalg.inv(jacobian)
        B = np.zeros((6, 8))
        B[0, 0::2], B[1, 1::2] = gradient[:, 0], gradient[:, 1]
        B[3, 0::2], B[3, 1::2] = gradient[:, 1], gradient[:, 0]
        matrices.append(B)
        weights.append(weight)
        points.append(shape @ coordinates)
    return np.asarray(matrices), np.asarray(weights), np.asarray(points)


class J2Quad4Host:
    """Case-local accepted-state ownership around the stateless FE driver."""

    def __init__(self, module, nodes, elements):
        self.module = module
        self.nodes, self.elements = np.asarray(nodes), np.asarray(elements)
        kinematics = [quad4_kinematics(self.nodes[c]) for c in self.elements]
        self.B = np.asarray([entry[0] for entry in kinematics])
        self.weights = np.asarray([entry[1] for entry in kinematics])
        self.points = np.asarray([entry[2] for entry in kinematics])
        self.stress = np.zeros((len(elements), 4, 6))
        self.strain = np.zeros_like(self.stress)
        self.ep = np.zeros((len(elements), 4))
        self.accepted_displacement = np.zeros(2 * len(nodes))
        self.time, self.dt, self.version = 0., 1., 0
        self._protect_accepted_state()

    def _protect_accepted_state(self):
        for value in (self.stress, self.strain, self.ep, self.accepted_displacement):
            value.flags.writeable = False

    def begin_increment(self, time_begin, dt):
        if not np.isfinite([time_begin, dt]).all() or time_begin < 0 or dt <= 0:
            raise ValueError("finite nonnegative time and positive increment required")
        self.time, self.dt = float(time_begin), float(dt)

    def evaluate_element(self, ei, current, increment):
        """Return a disposable candidate; accepted arrays never receive output."""
        old_strain = self.B[ei] @ (current - increment)
        if not np.allclose(old_strain, self.strain[ei], rtol=0, atol=1e-12):
            raise AssertionError("nodal and accepted Gauss-point strain history disagree")
        residual, tangent = np.zeros(8), np.zeros((8, 8))
        stresses, states = [], []
        for point, (B, weight) in enumerate(zip(self.B[ei], self.weights[ei])):
            stress, ep, C = call_umat(self.module, self.stress[ei, point], self.ep[ei, point],
                                     self.strain[ei, point], B @ increment, self.time, self.dt)
            residual += weight * B.T @ stress
            tangent += weight * B.T @ C @ B
            stresses.append(stress)
            states.append(ep)
        return residual, tangent, np.asarray(stresses), np.asarray(states), self.B[ei] @ current

    def element_rk(self, ei, coordinates, current, increment):
        if not np.allclose(coordinates, self.nodes[self.elements[ei]], rtol=0, atol=1e-14):
            raise ValueError("element identity changed during the J2 history")
        return self.evaluate_element(ei, current, increment)[:2]

    def commit_equilibrium(self, displacement, conditions, tol=1e-10):
        """Accept all GP states atomically only after free equilibrium and BCs."""
        displacement = np.asarray(displacement)
        free = np.array(sorted(set(range(len(displacement))) - set(conditions)), dtype=int)
        if not len(free):
            raise AssertionError("plasticity BVP must solve free displacement DOFs")
        residual = np.zeros_like(displacement)
        stress, ep, strain = [], [], []
        for ei, connectivity in enumerate(self.elements):
            dofs = (2 * connectivity[:, None] + np.arange(2)).ravel()
            local = self.evaluate_element(ei, displacement[dofs],
                                          displacement[dofs] - self.accepted_displacement[dofs])
            np.add.at(residual, dofs, local[0])
            stress.append(local[2])
            ep.append(local[3])
            strain.append(local[4])
        if not np.isfinite(residual).all() or np.max(abs(residual[free])) > tol:
            raise AssertionError("unconverged plasticity trial cannot commit state")
        if max(abs(displacement[i] - value) for i, value in conditions.items()) > 1e-12:
            raise AssertionError("plasticity boundary values are not satisfied")
        if np.any(np.asarray(ep) < self.ep - 1e-12):
            raise AssertionError("equivalent plastic strain decreased")
        self.stress, self.ep, self.strain = np.asarray(stress), np.asarray(ep), np.asarray(strain)
        self.accepted_displacement = displacement.copy()
        self.version += 1
        self._protect_accepted_state()
        return residual


def oracle_assembly(nodes, elements, displacement, old_displacement, stress_old, ep_old):
    """Independent rectangular-Q1 force assembly and radial-return evolution.

    Spatial derivatives are written from the four corner values and rectangle
    lengths rather than reusing the host's B matrices or its assembly routine.
    """
    force = np.zeros(2 * len(nodes))
    stresses, states, strains = [], [], []
    for ei, connectivity in enumerate(elements):
        coordinates = nodes[connectivity]
        dx, dy = coordinates[1, 0] - coordinates[0, 0], coordinates[3, 1] - coordinates[0, 1]
        current = displacement.reshape(-1, 2)[connectivity]
        increment = (displacement - old_displacement).reshape(-1, 2)[connectivity]
        e_stress, e_state, e_strain = [], [], []
        for p, (xi, eta) in enumerate(GAUSS):
            gx = np.array([-(1 - eta), 1 - eta, 1 + eta, -(1 + eta)]) / (2 * dx)
            gy = np.array([-(1 - xi), -(1 + xi), 1 + xi, 1 - xi]) / (2 * dy)
            deps = np.array([gx @ increment[:, 0], gy @ increment[:, 1], 0.,
                             gy @ increment[:, 0] + gx @ increment[:, 1], 0., 0.])
            total = np.array([gx @ current[:, 0], gy @ current[:, 1], 0.,
                              gy @ current[:, 0] + gx @ current[:, 1], 0., 0.])
            stress, ep = radial_return_reference(stress_old[ei, p], ep_old[ei, p], deps)
            point_force = np.column_stack((stress[0] * gx + stress[3] * gy,
                                           stress[3] * gx + stress[1] * gy)) * (dx * dy / 4)
            np.add.at(force.reshape(-1, 2), connectivity, point_force)
            e_stress.append(stress)
            e_state.append(ep)
            e_strain.append(total)
        stresses.append(e_stress)
        states.append(e_state)
        strains.append(e_strain)
    return force, np.asarray(stresses), np.asarray(states), np.asarray(strains)


def plastic_dissipation(strain, stress, ep, old_strain, old_stress, old_ep, weights):
    """Gate endpoint plastic work, return physical dissipation sigy*Delta ep.

    Endpoint plastic work less the change of stored hardening energy includes
    H*(Delta ep)^2/2 numerical dissipation as well as sigy*Delta ep. The latter
    is the rate-independent physical plastic dissipation recorded in histories.
    """
    old_plastic = old_strain - elastic_strain(old_stress)
    new_plastic = strain - elastic_strain(stress)
    work = np.sum(stress * (new_plastic - old_plastic), axis=-1)
    hardening_energy = .5 * PROPERTIES["H"] * (ep ** 2 - old_ep ** 2)
    dissipated = work - hardening_energy
    dep = ep - old_ep
    expected = PROPERTIES["sigma_y"] * dep + .5 * PROPERTIES["H"] * dep ** 2
    if (not np.isfinite(dissipated).all() or np.min(dissipated) < -1e-11
            or np.max(abs(dissipated - expected)) > 1e-10):
        raise AssertionError("plastic dissipation violates the independent work identity")
    return float(np.sum(weights * PROPERTIES["sigma_y"] * dep))


def endpoint_work_balance(displacement, old_displacement, residual, prescribed,
                          strain, stress, ep, old_strain, old_stress, old_ep, weights):
    """Independent exact BE/radial-return work identity, including sigma_zz.

    R_new*Delta u = Delta(We + H ep^2/2) + sigy Delta ep
                    + Delta sigma:C^-1:Delta sigma/2 + H(Delta ep)^2/2.
    The last two terms are numerical endpoint-work dissipation; they are not
    included in the reported physical plastic dissipation.
    """
    elastic_energy = .5 * np.sum(stress * elastic_strain(stress), axis=-1)
    old_elastic_energy = .5 * np.sum(old_stress * elastic_strain(old_stress), axis=-1)
    dep, dsigma = ep - old_ep, stress - old_stress
    storage = elastic_energy - old_elastic_energy + .5 * PROPERTIES["H"] * (ep ** 2 - old_ep ** 2)
    plastic = PROPERTIES["sigma_y"] * dep
    numerical = .5 * np.sum(dsigma * elastic_strain(dsigma), axis=-1) + .5 * PROPERTIES["H"] * dep ** 2
    internal_work = float(np.sum(weights * np.sum(stress * (strain - old_strain), axis=-1)))
    expected_work = float(np.sum(weights * (storage + plastic + numerical)))
    delta_u = displacement - old_displacement
    boundary_work = float(sum(residual[i] * delta_u[i] for i in prescribed))
    error = max(abs(boundary_work - expected_work), abs(internal_work - expected_work))
    if (not np.isfinite([boundary_work, internal_work, expected_work]).all()
            or np.min(plastic) < -1e-12 or np.min(numerical) < -1e-12 or error > 1e-10):
        raise AssertionError("independent discrete plasticity endpoint-work balance failed")
    return error


def yield_gate(stress, ep, old_ep):
    """Admissibility everywhere, consistency on each accepted plastic branch."""
    deviator = stress.copy()
    deviator[..., :3] -= np.mean(stress[..., :3], axis=-1)[..., None]
    q = np.sqrt(1.5 * (np.sum(deviator[..., :3] ** 2, axis=-1)
                       + 2 * np.sum(deviator[..., 3:] ** 2, axis=-1)))
    yield_function = q - PROPERTIES["sigma_y"] - PROPERTIES["H"] * ep
    plastic = ep - old_ep > 1e-12
    error = max(float(np.max(yield_function)),
                float(np.max(abs(yield_function[plastic]))) if np.any(plastic) else 0.)
    if not np.isfinite(yield_function).all() or error > 1e-10:
        raise AssertionError("J2 admissibility or plastic consistency failed")
    return max(error, 0.)


def top_displacement(time, peak=.02):
    """A slow load/unload pseudo-time cycle; no inertia or viscosity."""
    return peak * (time if time <= 1 else 2 - time)


def conditions(nodes, time):
    bottom = np.flatnonzero(np.isclose(nodes[:, 1], 0))
    top = np.flatnonzero(np.isclose(nodes[:, 1], 1))
    result = {2 * int(i) + component: 0. for i in bottom for component in (0, 1)}
    result.update({2 * int(i): top_displacement(time) for i in top})
    result.update({2 * int(i) + 1: 0. for i in top})
    return result


def assert_trial_purity(host, displacement, old_displacement):
    """Probe and reject a larger arbitrary trial, then verify exact rollback."""
    from abaqus_ufl.fe import assemble

    snapshot = tuple(value.copy() for value in (host.stress, host.strain, host.ep,
                                               host.accepted_displacement))
    version = host.version
    # This unconverged trial deliberately traverses plastic branches. Its state
    # is discarded, as every failed Newton/line-search candidate must be.
    trial = displacement.copy()
    trial[0::2] += .01 * np.sin(math.pi * host.nodes[:, 0]) * np.sin(math.pi * host.nodes[:, 1])
    assemble(host.nodes, host.elements, 2, host.element_rk, trial, trial - old_displacement)
    if (host.version != version or any(not np.array_equal(before, after)
            for before, after in zip(snapshot, (host.stress, host.strain, host.ep,
                                                host.accepted_displacement)))):
        raise AssertionError("rejected plasticity trial changed accepted state")


def solve_history(module, nx, ny=None, times=None, boundary_conditions=conditions):
    """Solve a non-homogeneous shear cycle and accept only certified states."""
    from abaqus_ufl.fe import newton_solve, structured_quad_mesh

    ny = nx if ny is None else ny
    times = np.linspace(0., 2., 21) if times is None else np.asarray(times)
    if len(times) < 3 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("cycle times must start at zero and strictly increase")
    nodes, elements = structured_quad_mesh(nx, ny)
    host = J2Quad4Host(module, nodes, elements)
    previous = np.zeros(2 * len(nodes))
    reference_stress, reference_ep = np.zeros_like(host.stress), np.zeros_like(host.ep)
    solutions, stresses, strains, states = [previous.copy()], [host.stress.copy()], [host.strain.copy()], [host.ep.copy()]
    reference_states = [reference_ep.copy()]
    residuals, reference_residuals, iterations, free_counts, dissipations = [], [], [], [], [0.]
    work_errors, yield_errors = [], []
    boundary_error, stress_error, ep_error, residual_error = 0., 0., 0., 0.
    for step, time in enumerate(times[1:], 1):
        host.begin_increment(float(times[step - 1]), float(time - times[step - 1]))
        prescribed = boundary_conditions(nodes, time)
        history = []
        solved = newton_solve(nodes, elements, 2, host.element_rk, prescribed,
                             U0=previous, U_previous=previous, max_iter=30, tol=1e-10,
                             record_iteration=history.append)
        assert_trial_purity(host, solved, previous)
        independent = oracle_assembly(nodes, elements, solved, previous, reference_stress, reference_ep)
        old_stress, old_strain, old_ep = host.stress, host.strain, host.ep
        residual = host.commit_equilibrium(solved, prescribed)
        free = np.array(sorted(set(range(len(solved))) - set(prescribed)), dtype=int)
        stress_error = max(stress_error, float(np.max(abs(host.stress - independent[1]))))
        ep_error = max(ep_error, float(np.max(abs(host.ep - independent[2]))))
        residual_error = max(residual_error, float(np.max(abs(residual - independent[0]))))
        if stress_error > 1e-9 or ep_error > 1e-10 or residual_error > 1e-10:
            raise AssertionError("compiled mesh history disagrees with independent radial return")
        if np.max(abs(independent[0][free])) > 1e-9:
            raise AssertionError("independent plasticity mesh residual is not in equilibrium")
        dissipations.append(plastic_dissipation(host.strain, host.stress, host.ep,
                                               old_strain, old_stress, old_ep, host.weights))
        work_errors.append(endpoint_work_balance(solved, previous, residual, prescribed,
                           host.strain, host.stress, host.ep, old_strain, old_stress, old_ep, host.weights))
        yield_errors.append(yield_gate(host.stress, host.ep, old_ep))
        boundary_error = max(boundary_error, max(abs(solved[i] - v) for i, v in prescribed.items()))
        solutions.append(solved.copy())
        stresses.append(host.stress.copy())
        strains.append(host.strain.copy())
        states.append(host.ep.copy())
        reference_states.append(independent[2].copy())
        residuals.append(residual)
        reference_residuals.append(independent[0])
        iterations.append(history)
        free_counts.append(len(free))
        reference_stress, reference_ep = independent[1], independent[2]
        previous = solved
    return {"nodes": nodes, "elements": elements, "times": np.asarray(times),
            "solutions": np.asarray(solutions), "stress": np.asarray(stresses),
            "strain": np.asarray(strains), "ep": np.asarray(states),
            "reference_ep": np.asarray(reference_states), "weights": host.weights,
            "gauss_points": host.points, "residuals": np.asarray(residuals),
            "reference_residuals": np.asarray(reference_residuals), "iterations": iterations,
            "free_dof_counts": free_counts, "dissipation": np.asarray(dissipations),
            "boundary_condition_error": boundary_error, "stress_oracle_error": stress_error,
            "state_oracle_error": ep_error, "force_oracle_error": residual_error,
            "endpoint_work_balance_errors": np.asarray(work_errors),
            "yield_consistency_errors": np.asarray(yield_errors)}


def top_reactions(result, reference=False):
    top = np.flatnonzero(np.isclose(result["nodes"][:, 1], 1))
    residual = result["reference_residuals"] if reference else result["residuals"]
    return np.r_[0., np.sum(residual[:, 2 * top], axis=1)]


def validate_history(result):
    """Finite bridge, state monotonicity, free DOFs, force/moment and regimes."""
    nodes, elements, times = result["nodes"], result["elements"], result["times"]
    nt, ne = len(times), len(elements)
    expected = {"solutions": (nt, 2 * len(nodes)), "stress": (nt, ne, 4, 6),
                "strain": (nt, ne, 4, 6), "ep": (nt, ne, 4), "weights": (ne, 4),
                "reference_ep": (nt, ne, 4),
                "gauss_points": (ne, 4, 2), "residuals": (nt - 1, 2 * len(nodes)),
                "reference_residuals": (nt - 1, 2 * len(nodes)), "dissipation": (nt,),
                "endpoint_work_balance_errors": (nt - 1,), "yield_consistency_errors": (nt - 1,)}
    if (times[0] != 0 or np.any(np.diff(times) <= 0) or not np.isfinite(times).all()
            or nodes.shape != (len(nodes), 2) or not np.isfinite(nodes).all()):
        raise AssertionError("invalid plasticity time/node bridge")
    for name, shape in expected.items():
        value = np.asarray(result[name])
        if value.shape != shape or not np.isfinite(value).all():
            raise AssertionError("invalid plasticity history: " + name)
    nx, ny = len(np.unique(nodes[:, 0])) - 1, len(np.unique(nodes[:, 1])) - 1
    if result["free_dof_counts"] != [2 * (nx + 1) * (ny - 1)] * (nt - 1):
        raise AssertionError("plasticity BVP must solve the expected free displacement DOFs")
    increments = np.diff(result["ep"], axis=0)
    if np.min(increments) < -1e-12 or np.min(result["dissipation"]) < -1e-12:
        raise AssertionError("plastic strain or cumulative dissipation decreased")
    force_error, moment_error, free_error = [], [], []
    for time, residual in zip(times[1:], result["residuals"]):
        prescribed = conditions(nodes, time)
        free = sorted(set(range(len(residual))) - set(prescribed))
        reactions = np.zeros((len(nodes), 2))
        for dof in prescribed:
            reactions.ravel()[dof] = residual[dof]
        # Symmetric small-strain stress balances moments about REFERENCE nodes.
        force_error.append(float(np.max(abs(np.sum(reactions, axis=0)))))
        moment_error.append(float(abs(np.sum(nodes[:, 0] * reactions[:, 1]
                                             - nodes[:, 1] * reactions[:, 0]))))
        free_error.append(float(np.max(abs(residual[free]))))
    if max(force_error + moment_error + free_error) > 1e-9:
        raise AssertionError("plasticity force, moment or free equilibrium failed")
    peak = int(np.argmax([top_displacement(time) for time in times]))
    if (peak == 0 or peak == nt - 1 or np.max(result["ep"][peak]) < .002
            or np.sum(result["ep"][peak] > 1e-8) < ne
            or np.ptp(result["ep"][peak]) < .001):
        raise AssertionError("non-homogeneous yielded Gauss-point regime was not reached")
    # At least one accepted initial state must be wholly elastic; the first
    # unloading interval must genuinely sample elastic unloading after yield.
    elastic = np.any(np.max(result["ep"][1:peak], axis=(1, 2)) < 1e-12)
    unloading = np.sum((result["ep"][peak] > 1e-8)
                       & (abs(increments[peak]) < 1e-12))
    reverse_plastic = np.any((increments[peak:] > 1e-10)
                            & (result["stress"][peak + 1:, :, :, 3] < -1e-8), axis=0)
    if not elastic or unloading < ne or top_reactions(result)[peak + 1] >= top_reactions(result)[peak]:
        raise AssertionError("elastic loading/plastic loading/elastic unloading not all sampled")
    if np.sum(reverse_plastic) < ne:
        raise AssertionError("reverse-shear plastic loading was not sampled")
    if result["dissipation"].sum() <= 0 or result["boundary_condition_error"] > 1e-12:
        raise AssertionError("plastic dissipation or boundary-condition gate failed")
    if np.max(result["endpoint_work_balance_errors"]) > 1e-10 or np.max(result["yield_consistency_errors"]) > 1e-10:
        raise AssertionError("independent work or yield history gate failed")
    return {"force_balance_max_abs_error": max(force_error),
            "moment_balance_max_abs_error": max(moment_error),
            "free_residual_max_abs": max(free_error),
            "yielded_gauss_points_at_peak": int(np.sum(result["ep"][peak] > 1e-8)),
            "elastic_unloading_gauss_points": int(unloading),
            "reverse_plastic_gauss_points": int(np.sum(reverse_plastic))}


def project_plastic_strain(result):
    """Visualization-only volume-weighted corner projection of raw GP ep.

    Each incident element contributes its Gauss-weighted mean to each corner,
    weighted by reference element area. This bounded averaging is NOT an
    authoritative nodal internal variable and does not extrapolate Gauss data.
    """
    elements, weights, ep = result["elements"], result["weights"], result["ep"]
    area = np.sum(weights, axis=1)
    numerator = np.zeros((len(ep), len(result["nodes"])))
    denominator = np.zeros(len(result["nodes"]))
    element_mean = np.sum(ep * weights[None, :, :], axis=2) / area[None, :]
    for element, connectivity in enumerate(elements):
        numerator[:, connectivity] += (element_mean[:, element] * area[element])[:, None]
        denominator[connectivity] += area[element]
    if np.any(denominator <= 0):
        raise AssertionError("plastic-strain projection has missing node coverage")
    projected = numerator / denominator
    if (not np.isfinite(projected).all() or np.min(projected) < -1e-12
            or np.any(projected.max(axis=1) > ep.max(axis=(1, 2)) + 1e-12)):
        raise AssertionError("invalid plastic-strain projection")
    return projected


def mesh_record(identifier, title, result, *, plastic=False):
    """Authoritative nodal U; projected ep includes full raw GP provenance."""
    from tools.livebench_data import mesh_field

    displacement = result["solutions"].reshape(len(result["times"]), -1, 2)
    values = project_plastic_strain(result) if plastic else np.linalg.norm(displacement, axis=2)
    description = ("Colors show a visualization-only volume-weighted nodal average of STATEV(1). "
                   "The authoritative equivalent plastic strain is stored at four Gauss points per element; "
                   "raw point identities and every accepted value accompany this plot. "
                   "The projection performs no extrapolation. " if plastic else
                   "Colors show authoritative nodal displacement magnitude. ")
    description += ("The bottom is clamped and the top loads to u₁ = 0.02 then returns to zero, with u₂ = 0. "
                    "Both sides are traction-free. Displacements are displayed at ×1.")
    plot = mesh_field(identifier, title, result["nodes"], result["elements"], values[-1], description,
                      displacements=displacement[-1])
    plot["deformation_scale"] = 1
    plot["value_label"] = "Projected equivalent plastic strain (visualization)" if plastic else "Displacement magnitude (model units)"
    plot["frames"] = [{"time": float(time), "values": value.tolist(), "displacements": u.tolist()}
                      for time, value, u in zip(result["times"], values, displacement)]
    if plastic:
        plot["gauss_point_provenance"] = {
            "source": "Generated UMAT STATEV(1), after accepted equilibrium only",
            "quantity": "Equivalent plastic strain", "units": "dimensionless",
            "element_index_base": 0, "gauss_point_index_base": 0,
            "identity_order": "frames.values[element_index][gauss_point_index], in plotted element order",
            "natural_coordinate_order": [list(point) for point in GAUSS],
            "reference_coordinates": result["gauss_points"].tolist(),
            "weights": result["weights"].tolist(),
            "frames": [{"time": float(time), "values": value.tolist()}
                       for time, value in zip(result["times"], result["ep"])],
            "projection": "Reference-area-weighted average of incident element Gauss-weighted means",
            "expected_points_per_frame": 4 * len(result["elements"]),
            "observed_points_per_frame": int(result["ep"].shape[1] * result["ep"].shape[2]),
            "unique_points_per_frame": 4 * len(result["elements"]),
            "missing_points_per_frame": 0, "duplicate_points_per_frame": 0,
            "nonfinite_points": int(np.sum(~np.isfinite(result["ep"]))),
        }
    for y, label in ((0., "Clamp: u₁ = u₂ = 0"), (1., "Top: prescribed load/unload, u₂ = 0")):
        indices = np.flatnonzero(np.isclose(result["nodes"][:, 1], y))
        plot.setdefault("boundaries", []).append({"label": label, "field": "displacement",
                                                "kind": "prescribed", "nodes": indices.tolist()})
    return plot
