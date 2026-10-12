"""Runtime guards in generated Fortran, checked by compiling and calling it.

- A finite-strain UMAT asked to evaluate det F <= 0 requests a cutback, as
  the UEL path already does (audit 2026-10-10: G6).
- UMATs and UELs request a cutback instead of returning NaN or Inf
  (audit 2026-10-10: G6, G7).
- det(F) inside a one-line return compiles and evaluates correctly (G1).

The singular materials divide by (det F - 1), which is 0/0 in the
undeformed state while det F stays positive. They are generated with
verify=False because they are deliberately singular at F = I.
"""

import numpy as np
import pytest

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, eye, inv, log, trace
from abaqus_ufl.generators.uel_gen import generate_uel
from abaqus_ufl.generators.umat_gen import generate_small_strain_umat, generate_umat

from compiled_helpers import build_uel_module, build_umat_module, call_uel, call_umat

G, K = 1.0, 10.0
SHEAR = np.array([[1.10, 0.30, 0.00],
                  [0.00, 1.00, 0.20],
                  [0.00, 0.00, 0.95]])
INVERTED = np.diag([1.0, 1.0, -0.9])
SQUARE = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])


class DirectReturnNeoHookean(au.Material):
    props = dict(G=G, K=K)

    def stress_PK1(self, F):
        return self.G * (F - inv(F).T) + self.K * log(det(F)) * inv(F).T


class SingularAtIdentity(au.Material):
    props = dict(G=G)

    def stress_PK1(self, F):
        return self.G * (F - inv(F).T) / (det(F) - 1.0)


class SingularSmallStrain(au.SmallStrainMaterial):
    props = dict(G=10.0)
    state_vars = dict(a=0.0)

    def stress_update(self, sigma_old, strain_old, dstrain, a_old, dt):
        scale = trace(dstrain) / trace(dstrain)
        return sigma_old + 2.0 * self.G * scale * dstrain, {'a': a_old}


def neo_hookean_cauchy(F):
    J = np.linalg.det(F)
    P = G * (F - np.linalg.inv(F).T) + K * np.log(J) * np.linalg.inv(F).T
    sigma = P @ F.T / J
    return np.array([sigma[0, 0], sigma[1, 1], sigma[2, 2], sigma[0, 1], sigma[0, 2], sigma[1, 2]])


@pytest.fixture(scope='module')
def direct_umat(tmp_path_factory):
    work = tmp_path_factory.mktemp('direct_umat')
    generate_umat(DirectReturnNeoHookean(), str(work / 'umat.for'))
    return build_umat_module(work / 'umat.for', work, 'guard_direct_umat')


@pytest.fixture(scope='module')
def singular_umat(tmp_path_factory):
    work = tmp_path_factory.mktemp('singular_umat')
    generate_umat(SingularAtIdentity(), str(work / 'umat.for'), verify=False)
    return build_umat_module(work / 'umat.for', work, 'guard_singular_umat')


def test_direct_return_det_evaluates_correctly(direct_umat):
    stress, _, ddsdde, pnewdt = call_umat(direct_umat, SHEAR, [G, K])
    assert pnewdt == 1.0
    np.testing.assert_allclose(stress, neo_hookean_cauchy(SHEAR), rtol=1e-12, atol=1e-12)
    assert np.isfinite(ddsdde).all()


def test_umat_requests_a_cutback_for_an_inverted_element(direct_umat):
    _, _, _, pnewdt = call_umat(direct_umat, INVERTED, [G, K])
    assert pnewdt < 1.0


def test_umat_requests_a_cutback_for_nonfinite_output(singular_umat):
    _, _, _, pnewdt = call_umat(singular_umat, np.eye(3), [G])
    assert pnewdt < 1.0
    stress, _, ddsdde, pnewdt = call_umat(singular_umat, SHEAR, [G])
    assert pnewdt == 1.0 and np.isfinite(stress).all() and np.isfinite(ddsdde).all()


def test_small_strain_umat_requests_a_cutback_for_nonfinite_output(tmp_path):
    generate_small_strain_umat(SingularSmallStrain(), str(tmp_path / 'umat.for'), verify=False)
    module = build_umat_module(tmp_path / 'umat.for', tmp_path, 'guard_small_strain_umat')
    _, _, _, pnewdt = call_umat(module, np.eye(3), [10.0], dstran=np.zeros(6))
    assert pnewdt < 1.0
    stress, _, ddsdde, pnewdt = call_umat(module, np.eye(3), [10.0],
                                          dstran=[1e-3, 0, 0, 0, 0, 0])
    assert pnewdt == 1.0 and np.isfinite(stress).all() and np.isfinite(ddsdde).all()


class SingularProblem(au.WeakForm):
    material = SingularAtIdentity
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField('u', degree=1)

    def momentum_equation(self, v, F):
        return self.material.stress_PK1(F)


@pytest.mark.parametrize('formulation', ['standard', 'fbar_mechanics'])
def test_uel_requests_a_cutback_for_nonfinite_output(tmp_path, formulation):
    generate_uel(SingularProblem(), str(tmp_path / 'uel.for'), element='Quad4',
                 formulation=formulation, verify=False)
    module = build_uel_module(tmp_path / 'uel.for', tmp_path, f'guard_{formulation}_uel')
    _, _, _, pnewdt = call_uel(module, SQUARE, np.zeros(8), [G])
    assert pnewdt < 1.0
    stretched = np.array([0.0, 0.0, 0.1, 0.0, 0.1, 0.05, 0.0, 0.05])
    rhs, amatrx, _, pnewdt = call_uel(module, SQUARE, stretched, [G], du=stretched)
    assert pnewdt == 1.0 and np.isfinite(rhs).all() and np.isfinite(amatrx).all()
