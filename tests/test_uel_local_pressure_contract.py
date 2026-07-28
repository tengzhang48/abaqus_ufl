"""Contract and condensed-Jacobian tests for the local-pressure generator.

The element-local pressure route condenses one internal pressure per
element behind a local Newton solve, so it must obey the same UEL
contract as the other generators: unsupported procedure requests
(``LFLAGS(3)=3,4,6``) return zeroed arrays, and nonpositive ``det F``
states request a cutback instead of feeding an invalid deformation into
the material evaluation. A black-box central-difference comparison of
the returned condensed tangent (with the pressure re-solved from the
same committed state at every perturbation) closes the loop.

Uses a minimal pressure-coupled gel declaration (the public workflow
example) on a Quad4. Build dependencies: gfortran + a working
``numpy.f2py`` backend; environments without the toolchain skip.
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


@pytest.fixture(scope='module')
def generated_source(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('localp_contract')
    path = tmp / 'localp_gel_quad4.for'
    problem = PressureGelLocalPressureQuad4()
    au.generate_uel(problem, str(path), element='Quad4',
                    formulation='local_pressure')
    return path


def test_generated_source_carries_contract_guards(generated_source):
    src = generated_source.read_text()
    assert 'LFLAGS(3) .NE. 1 .AND. LFLAGS(3) .NE. 2' in src
    assert 'LFLAGS(1) .NE. 1 .AND. LFLAGS(1) .NE. 2' in src
    assert 'LFLAGS(4) .NE. 0' in src
    assert 'DABS(K_pp) .GT. 1.0d30' in src
    assert src.count('det33d(F) .LE. 0.0d0') >= 2
    assert 'RHS(i,2) = 0.0d0' in src
    assert 'conv_p .EQ. 0' in src


@pytest.fixture(scope='module')
def compiled(generated_source):
    tmp = generated_source.parent
    module_name = 'localp_contract_quad4'
    cmd = [
        sys.executable, '-m', 'numpy.f2py', '-c',
        DRIVER, str(generated_source),
        '-m', module_name, 'only:', 'drive_uel', ':',
    ]
    result = subprocess.run(
        cmd, cwd=str(tmp), capture_output=True, text=True, check=False)
    if result.returncode != 0:
        combined = result.stdout + '\n' + result.stderr
        if any(sig in combined for sig in TOOLCHAIN_SIGNATURES):
            pytest.skip('f2py backend unavailable')
        raise RuntimeError('f2py build failed:\n' + combined[-2000:])
    sys.path.insert(0, str(tmp))
    try:
        sys.modules.pop(module_name, None)
        module = __import__(module_name)
    finally:
        sys.path.pop(0)
    return module


COORDS = np.array([[0.0, 1.0, 1.0, 0.0],
                   [0.0, 0.0, 1.0, 1.0]])


def _state():
    rng = np.linspace(-0.02, 0.02, 8)
    U = np.zeros(12)
    k = 0
    for a in range(4):
        for i in range(2):
            U[a * 3 + i] = rng[k]
            k += 1
        U[a * 3 + 2] = MU_EQ * (1.0 + 0.01 * (a - 1))
    return U, 0.3 * U


def _call(module, U, DU, svars_in=0.0, lflags3=1, lflags1=72, lflags4=0):
    svars = np.array([svars_in], dtype=float)
    props = np.array(list(PROPS.values()))
    rhs, amatrx, svars_out, pnewdt = module.drive_uel(
        svars, COORDS, np.asarray(U, float), np.asarray(DU, float), props,
        np.zeros(1, dtype=np.int32), np.array([0.1, 0.1]), 0.1, 1.0,
        np.array([lflags1, 1, lflags3, lflags4, 0, 0], dtype=np.int32),
        np.zeros(3), 1, 1, 1, 1, 0.0,
    )
    return rhs, amatrx, float(svars_out[0]), float(pnewdt)


def test_unsupported_requests_return_zeroed_arrays(compiled):
    U, DU = _state()
    for req in (3, 4, 6, 100):
        rhs, amatrx, _, _ = _call(compiled, U, DU, lflags3=req)
        assert np.max(np.abs(rhs)) == 0.0, req
        assert np.max(np.abs(amatrx)) == 0.0, req
    rhs, amatrx, _, _ = _call(compiled, U, DU, lflags1=99)
    assert np.max(np.abs(rhs)) == 0.0
    assert np.max(np.abs(amatrx)) == 0.0
    rhs, amatrx, _, _ = _call(compiled, U, DU, lflags4=1)
    assert np.max(np.abs(rhs)) == 0.0, "perturbation step not rejected"
    assert np.max(np.abs(amatrx)) == 0.0


def test_nonfinite_local_solve_requests_cutback(compiled):
    """An overflowing incoming pressure (exp(p/K) -> Inf) must cut back."""
    U, DU = _state()
    _, _, _, pnewdt = _call(compiled, U, DU, svars_in=1.0e5)
    assert pnewdt < 1.0, pnewdt


def test_stiffness_only_call_does_not_commit_state(compiled):
    """LFLAGS(3)=2 must not mutate SVARS (state commits are type-1 only)."""
    U, DU = _state()
    _, _, p_normal, _ = _call(compiled, U, DU, lflags3=1)
    assert abs(p_normal) > 1.0e-12          # normal call does commit
    _, _, p_stiff, _ = _call(compiled, U, DU, lflags3=2)
    assert p_stiff == 0.0, p_stiff          # incoming state preserved


def test_nonpositive_detF_requests_cutback(compiled):
    U, DU = _state()
    U_bad = U.copy()
    for node in range(4):
        U_bad[node * 3] = -2.0 * COORDS[0][node]
    _, _, _, pnewdt = _call(compiled, U_bad, DU)
    assert pnewdt < 1.0, pnewdt


def test_condensed_jacobian_black_box(compiled):
    U, DU = _state()
    rhs0, amatrx0, p_star, pnewdt = _call(compiled, U, DU)
    assert pnewdt >= 1.0 and np.isfinite(rhs0).all()

    fd = np.zeros((12, 12))
    for j in range(12):
        eps = 1.0e-6 * max(1.0, abs(U[j]))
        Up, Um = U.copy(), U.copy()
        Up[j] += eps
        Um[j] -= eps
        DUp, DUm = DU.copy(), DU.copy()
        DUp[j] += eps
        DUm[j] -= eps
        rhs_p, _, _, _ = _call(compiled, Up, DUp)
        rhs_m, _, _, _ = _call(compiled, Um, DUm)
        fd[:, j] = (rhs_p - rhs_m) / (2.0 * eps)
    rel = float(np.max(np.abs(amatrx0 - (-fd)))) / float(np.max(np.abs(amatrx0)))
    assert np.isfinite(rel) and rel < 5.0e-5, rel
