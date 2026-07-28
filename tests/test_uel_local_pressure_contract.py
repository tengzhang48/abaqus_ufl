"""Contract and condensed-Jacobian tests for the local-pressure generator.

The element-local pressure route condenses one internal pressure per
element behind a local Newton solve, so it must obey the same UEL
contract as the other generators: unsupported procedure requests
(``LFLAGS(3)=3,4,6``) return zeroed arrays, and nonpositive ``det F``
states request a cutback instead of feeding an invalid deformation into
the material evaluation.

The condensed tangent is verified by a black-box central-difference
comparison on BOTH supported elements (Quad4, 12 DOFs, and Hex8,
32 DOFs): every global element variable is perturbed, the wrapper
re-solves the element pressure from the same committed state after
every perturbation, and the derivative of the returned condensed
residual is compared with the returned tangent, so the comparison
covers the implicit pressure sensitivity carried by the Schur term. A
frozen-pressure negative control (the same source with only the Schur
subtraction removed) must FAIL the identical comparison, proving the
test discriminates, and a companion check confirms the local solve is
idempotent at the converged pressure.

Uses a minimal pressure-coupled gel declaration (the public workflow
example). Build dependencies: gfortran + a working ``numpy.f2py``
backend; environments without the toolchain skip.
"""
import os
import subprocess
import sys
from collections import OrderedDict

import numpy as np
import pytest

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, exp, inv, log

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DRIVER = os.path.join(ROOT, 'examples', 'scalar_diffusion_uel', 'f2py',
                      'drive_uel.f90')

PROPS = OrderedDict([
    ('G', 1.0), ('K', 100.0), ('chi', 0.1), ('D', 5.0e-9), ('mu0', 0.0),
    ('Omega', 1.0e-4), ('Rgas', 8.314), ('theta', 298.0), ('phi0', 0.5),
])
RT = PROPS['Rgas'] * PROPS['theta']
MU_EQ = PROPS['mu0'] + RT * (np.log(1.0 - PROPS['phi0']) + PROPS['phi0']
                             + PROPS['chi'] * PROPS['phi0'] ** 2)

TOOLCHAIN_SIGNATURES = (
    "Compiler.__init__()",
    "No module named 'distutils'",
    "No module named 'numpy.distutils'",
)


class PressureGelMaterial(au.Material):
    props = dict(PROPS)

    def stress_PK1(self, F, p, mu):
        return self.G * (F - inv(F).T) + p * inv(F).T

    def pressure_resid(self, F, p, mu):
        J = det(F)
        Je = exp(p / self.K)
        phi = self.phi0 * Je / J
        return (mu - self.mu0
                - self.Rgas * self.theta
                * (log(1.0 - phi) + phi + self.chi * phi ** 2)
                + (phi / self.phi0) * self.Omega * p)

    def solvent_flux(self, F, p, mu, grad_mu):
        J = det(F)
        Je = exp(p / self.K)
        phi = self.phi0 * Je / J
        Cinv = inv(F.T @ F)
        cR0 = (1.0 - self.phi0) / self.Omega
        cR = cR0 + (self.phi0 - phi) / (self.Omega * phi)
        M = self.D * cR / (self.Rgas * self.theta)
        return -M * (Cinv @ grad_mu)

    def solvent_storage(self, F, F_old, p, p_old, dt):
        J = det(F)
        J_old = det(F_old)
        Je = exp(p / self.K)
        Je_old = exp(p_old / self.K)
        return (J / Je - J_old / Je_old) / (self.Omega * dt)


class PressureGelLocalPressureQuad4(au.WeakForm):
    material = PressureGelMaterial
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField('u', degree=1)
        self.mu = au.ScalarField('mu', degree=1)
        self.p = au.LocalScalar('p', storage='SVARS',
                                condensed=True, initial=0.0)

    def momentum_equation(self, v, F, p, mu):
        return self.material.stress_PK1(F, p, mu)

    def pressure_equation(self, q, F, p, mu):
        return self.material.pressure_resid(F, p, mu)

    def transport_equation(self, w, F, p, mu, grad_mu, F_old, p_old, dt):
        c_dot = self.material.solvent_storage(F, F_old, p, p_old, dt)
        j_R = self.material.solvent_flux(F, p, mu, grad_mu)
        return c_dot, j_R


class PressureGelLocalPressureHex8(PressureGelLocalPressureQuad4):
    ndim = 3


CASES = {
    'quad4': dict(
        problem=PressureGelLocalPressureQuad4,
        element='Quad4',
        ndim=2,
        coords=np.array([[0.0, 1.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0, 1.0]]),
    ),
    'hex8': dict(
        problem=PressureGelLocalPressureHex8,
        element='Hex8',
        ndim=3,
        coords=np.array([[0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0],
                         [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]]),
    ),
}


@pytest.fixture(scope='module', params=sorted(CASES))
def generated_source(request, tmp_path_factory):
    name = request.param
    case = CASES[name]
    tmp = tmp_path_factory.mktemp('localp_contract_' + name)
    path = tmp / ('localp_gel_%s.for' % name)
    problem = case['problem']()
    au.generate_uel(problem, str(path), element=case['element'],
                    formulation='local_pressure')
    return name, case, path


def test_generated_source_carries_contract_guards(generated_source):
    _, _, path = generated_source
    src = path.read_text()
    assert 'LFLAGS(3) .NE. 1 .AND. LFLAGS(3) .NE. 2' in src
    assert 'LFLAGS(1) .NE. 1 .AND. LFLAGS(1) .NE. 2' in src
    assert 'LFLAGS(4) .NE. 0' in src
    assert 'DABS(K_pp) .GT. 1.0d30' in src
    assert src.count('det33d(F) .LE. 0.0d0') >= 2
    assert 'RHS(i,2) = 0.0d0' in src
    assert 'conv_p .EQ. 0' in src


def _build(source_path, module_name, cwd):
    cmd = [
        sys.executable, '-m', 'numpy.f2py', '-c',
        DRIVER, str(source_path),
        '-m', module_name, 'only:', 'drive_uel', ':',
    ]
    result = subprocess.run(
        cmd, cwd=str(cwd), capture_output=True, text=True, check=False)
    if result.returncode != 0:
        combined = result.stdout + '\n' + result.stderr
        if any(sig in combined for sig in TOOLCHAIN_SIGNATURES):
            pytest.skip('f2py backend unavailable')
        raise RuntimeError('f2py build failed:\n' + combined[-2000:])
    sys.path.insert(0, str(cwd))
    try:
        sys.modules.pop(module_name, None)
        module = __import__(module_name)
    finally:
        sys.path.pop(0)
    return module


@pytest.fixture(scope='module')
def compiled(generated_source):
    name, case, path = generated_source
    module = _build(path, 'localp_contract_' + name, path.parent)
    return name, case, module


def _state(case):
    """Generic deformed/chemical state near the free-swell datum."""
    ndim = case['ndim']
    nnode = case['coords'].shape[1]
    ndofel = (ndim + 1) * nnode
    rng = np.linspace(-0.02, 0.02, ndim * nnode)
    U = np.zeros(ndofel)
    k = 0
    for a in range(nnode):
        for i in range(ndim):
            U[a * (ndim + 1) + i] = rng[k]
            k += 1
        U[a * (ndim + 1) + ndim] = MU_EQ * (1.0 + 0.01 * (a - 1))
    return U, 0.3 * U


def _call(module, case, U, DU, svars_in=0.0, lflags3=1, lflags1=72,
          lflags4=0):
    svars = np.array([svars_in], dtype=float)
    props = np.array(list(PROPS.values()))
    rhs, amatrx, svars_out, pnewdt = module.drive_uel(
        svars, case['coords'], np.asarray(U, float), np.asarray(DU, float),
        props, np.zeros(1, dtype=np.int32), np.array([0.1, 0.1]), 0.1, 1.0,
        np.array([lflags1, 1, lflags3, lflags4, 0, 0], dtype=np.int32),
        np.zeros(3), 1, 1, 1, 1, 0.0,
    )
    return rhs, amatrx, float(svars_out[0]), float(pnewdt)


def test_unsupported_requests_return_zeroed_arrays(compiled):
    _, case, module = compiled
    U, DU = _state(case)
    for req in (3, 4, 6, 100):
        rhs, amatrx, _, _ = _call(module, case, U, DU, lflags3=req)
        assert np.max(np.abs(rhs)) == 0.0, req
        assert np.max(np.abs(amatrx)) == 0.0, req
    rhs, amatrx, _, _ = _call(module, case, U, DU, lflags1=99)
    assert np.max(np.abs(rhs)) == 0.0
    assert np.max(np.abs(amatrx)) == 0.0
    rhs, amatrx, _, _ = _call(module, case, U, DU, lflags4=1)
    assert np.max(np.abs(rhs)) == 0.0, "perturbation step not rejected"
    assert np.max(np.abs(amatrx)) == 0.0


def test_nonfinite_local_solve_requests_cutback(compiled):
    """An overflowing incoming pressure (exp(p/K) -> Inf) must cut back."""
    _, case, module = compiled
    U, DU = _state(case)
    _, _, _, pnewdt = _call(module, case, U, DU, svars_in=1.0e5)
    assert pnewdt < 1.0, pnewdt


def test_stiffness_only_call_does_not_commit_state(compiled):
    """LFLAGS(3)=2 must not mutate SVARS (state commits are type-1 only)."""
    _, case, module = compiled
    U, DU = _state(case)
    _, _, p_normal, _ = _call(module, case, U, DU, lflags3=1)
    assert abs(p_normal) > 1.0e-12          # normal call does commit
    _, _, p_stiff, _ = _call(module, case, U, DU, lflags3=2)
    assert p_stiff == 0.0, p_stiff          # incoming state preserved


def test_nonpositive_detF_requests_cutback(compiled):
    _, case, module = compiled
    ndim = case['ndim']
    U, DU = _state(case)
    U_bad = U.copy()
    for node in range(case['coords'].shape[1]):
        U_bad[node * (ndim + 1)] = -2.0 * case['coords'][0][node]
    _, _, _, pnewdt = _call(module, case, U_bad, DU)
    assert pnewdt < 1.0, pnewdt


def test_condensed_jacobian_black_box(compiled):
    name, case, module = compiled
    U, DU = _state(case)
    ndofel = U.size
    rhs0, amatrx0, p_star, pnewdt = _call(module, case, U, DU)
    assert pnewdt >= 1.0 and np.isfinite(rhs0).all()

    fd = np.zeros((ndofel, ndofel))
    dp = np.zeros(ndofel)
    for j in range(ndofel):
        eps = 1.0e-6 * max(1.0, abs(U[j]))
        Up, Um = U.copy(), U.copy()
        Up[j] += eps
        Um[j] -= eps
        DUp, DUm = DU.copy(), DU.copy()
        DUp[j] += eps
        DUm[j] -= eps
        rhs_p, _, pp, _ = _call(module, case, Up, DUp)
        rhs_m, _, pm, _ = _call(module, case, Um, DUm)
        fd[:, j] = (rhs_p - rhs_m) / (2.0 * eps)
        dp[j] = (pp - pm) / (2.0 * eps)
    rel = (float(np.max(np.abs(amatrx0 - (-fd))))
           / float(np.max(np.abs(amatrx0))))
    assert np.isfinite(rel) and rel < 5.0e-5, (name, rel)

    # Discriminator: the re-solved pressure must respond to the global
    # perturbations, otherwise the Schur term is dormant and the test
    # could not distinguish the condensed tangent from the frozen-p one.
    assert float(np.max(np.abs(dp))) > 1.0e-6, (name, float(np.max(np.abs(dp))))


def test_pressure_solve_is_idempotent(compiled):
    """Re-entering with the converged pressure returns the same pressure."""
    name, case, module = compiled
    U, DU = _state(case)
    _, _, p1, _ = _call(module, case, U, DU, svars_in=0.0)
    _, _, p2, _ = _call(module, case, U, DU, svars_in=p1)
    assert abs(p2 - p1) < 1.0e-8 * (1.0 + abs(p1)), (name, p1, p2)


FROZEN_SCHUR_LINE = 'AMATRX(i,j) = Kxx(i,j) - Kxp(i)*Kpx(j)/K_pp'
FROZEN_SCHUR_REPL = 'AMATRX(i,j) = Kxx(i,j)'


def test_frozen_pressure_negative_control(tmp_path):
    """The same comparison must FAIL when the Schur term is removed.

    Injected-failure control: generate the Quad4 source, delete only the
    Schur subtraction from the returned tangent (leaving the residual
    and the local solve untouched), compile, and run the identical
    black-box comparison. A large relative error here proves the main
    test discriminates the condensed tangent from the frozen-pressure
    one, so its pass is not vacuous.
    """
    case = CASES['quad4']
    path = tmp_path / 'localp_gel_quad4_control.for'
    au.generate_uel(case['problem'](), str(path), element=case['element'],
                    formulation='local_pressure')
    src = path.read_text()
    assert src.count(FROZEN_SCHUR_LINE) == 1, 'Schur line not found uniquely'
    broken_path = tmp_path / 'frozen_schur_quad4.for'
    broken_path.write_text(src.replace(FROZEN_SCHUR_LINE, FROZEN_SCHUR_REPL))
    module = _build(broken_path, 'localp_frozen_schur_control', tmp_path)

    U, DU = _state(case)
    ndofel = U.size
    rhs0, amatrx0, _, _ = _call(module, case, U, DU)
    fd = np.zeros((ndofel, ndofel))
    for j in range(ndofel):
        eps = 1.0e-6 * max(1.0, abs(U[j]))
        Up, Um = U.copy(), U.copy()
        Up[j] += eps
        Um[j] -= eps
        DUp, DUm = DU.copy(), DU.copy()
        DUp[j] += eps
        DUm[j] -= eps
        rhs_p, _, _, _ = _call(module, case, Up, DUp)
        rhs_m, _, _, _ = _call(module, case, Um, DUm)
        fd[:, j] = (rhs_p - rhs_m) / (2.0 * eps)
    rel = (float(np.max(np.abs(amatrx0 - (-fd))))
           / float(np.max(np.abs(amatrx0))))
    assert rel > 1.0e-3, (
        'frozen-pressure control PASSED the comparison (rel err %.3e); '
        'the main test would not discriminate' % rel)
