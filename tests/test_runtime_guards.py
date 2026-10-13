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


class LoopIndexArithmetic(au.Material):
    """INTEGER loop indices in k**(-1) and i/j: s = 11/6 + 9/2 as in Python."""
    props = dict(G=G)

    def stress_PK1(self, F):
        s = 0.0
        for k in range(1, 4):
            s = s + k ** (-1)
        for i in range(1, 3):
            for j in range(1, 3):
                s = s + i / j
        return self.G * s * (F - inv(F).T)


class SmallStrainRealArguments(au.SmallStrainMaterial):
    """abs() and a comparison on dt (real) and abs() on a_old (unperturbed)."""
    props = dict(G=10.0)
    state_vars = dict(a=0.0)

    def stress_update(self, sigma_old, strain_old, dstrain, a_old, dt):
        scale = 1.0 + 0.1 * abs(dt - 0.5) + 0.1 * abs(a_old)
        if dt > 0.25:
            scale = scale + 0.2
        return sigma_old + 2.0 * self.G * scale * dstrain, {'a': a_old}


class HeatWithRealArguments(au.Material):
    """abs() of X, time and the history value T_old in a UEL material."""
    props = dict(k=0.5, rho_cp=1.0, x0=0.5)

    def stress_PK1(self, F, T):
        return F - inv(F).T

    def solvent_storage(self, F, F_old, T, T_old, dt, time, X):
        ramp = 1.0 + 0.1 * abs(X[0] - self.x0) + 0.01 * abs(time) + 0.01 * abs(T_old)
        return self.rho_cp * ramp * (T - T_old) / dt

    def solvent_flux(self, F, T, grad_T):
        return -self.k * grad_T


class HeatWithRealArgumentsProblem(au.WeakForm):
    material = HeatWithRealArguments
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField('u', degree=1)
        self.T = au.ScalarField('T', degree=1, test='theta')

    def momentum_equation(self, v, F, T):
        return self.material.stress_PK1(F, T)

    def transport_equation(self, theta, F, T, grad_T, F_old, T_old, dt, time, X):
        return (self.material.solvent_storage(F, F_old, T, T_old, dt, time, X),
                self.material.solvent_flux(F, T, grad_T))


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


def test_integer_loop_arithmetic_evaluates_like_python(tmp_path):
    generate_umat(LoopIndexArithmetic(), str(tmp_path / 'umat.for'))
    module = build_umat_module(tmp_path / 'umat.for', tmp_path, 'guard_loop_arithmetic_umat')
    stress, _, _, pnewdt = call_umat(module, SHEAR, [G])
    J = np.linalg.det(SHEAR)
    s = (1.0 + 1.0 / 2 + 1.0 / 3) + (1.0 / 1 + 1.0 / 2 + 2.0 / 1 + 2.0 / 2)
    sigma = s * G * (SHEAR - np.linalg.inv(SHEAR).T) @ SHEAR.T / J
    expected = [sigma[0, 0], sigma[1, 1], sigma[2, 2], sigma[0, 1], sigma[0, 2], sigma[1, 2]]
    assert pnewdt == 1.0
    np.testing.assert_allclose(stress, expected, rtol=1e-12, atol=1e-12)


def test_abs_and_comparison_of_real_arguments_compile_and_run(tmp_path):
    """dt is DOUBLE PRECISION and a_old is never perturbed, so abs() of both
    and a comparison on dt are allowed and must compile."""
    generate_small_strain_umat(SmallStrainRealArguments(), str(tmp_path / 'umat.for'))
    module = build_umat_module(tmp_path / 'umat.for', tmp_path, 'guard_real_args_umat')
    stress, _, ddsdde, pnewdt = call_umat(module, np.eye(3), [10.0], dstran=[1e-3, 0, 0, 0, 0, 0])
    scale = 1.0 + 0.1 * abs(1.0 - 0.5) + 0.2           # dtime = 1 in call_umat
    assert pnewdt == 1.0
    np.testing.assert_allclose(stress[0], 2.0 * 10.0 * scale * 1e-3, rtol=1e-12)
    assert np.isfinite(ddsdde).all()


def test_abs_of_position_time_and_history_compiles_in_a_uel(tmp_path):
    generate_uel(HeatWithRealArgumentsProblem(), str(tmp_path / 'uel.for'),
                 element='Quad4', formulation='standard')
    module = build_uel_module(tmp_path / 'uel.for', tmp_path, 'guard_real_args_uel')
    u = np.array([0.0, 0.0, 0.1, 0.0, 0.0, 0.2, 0.1, 0.05, 0.3, 0.0, 0.0, 0.1])
    rhs, amatrx, _, pnewdt = call_uel(module, SQUARE, u, [0.5, 1.0, 0.5], du=0.5 * u)
    assert pnewdt == 1.0 and np.isfinite(rhs).all() and np.isfinite(amatrx).all()


def test_nonfinite_check_survives_fast_math(tmp_path):
    """isnonfinite_d tests exponent bits, so it still detects NaN and Inf when
    the compiler may assume no NaNs (where x .NE. x is folded to .FALSE.)."""
    import shutil
    import subprocess
    from pathlib import Path

    import abaqus_ufl.generators as generators

    gfortran = shutil.which('gfortran')
    if gfortran is None:
        pytest.skip('gfortran not available')
    template = Path(generators.__file__).with_name('templates') / 'tensor_ops.for'
    program = tmp_path / 'probe.f90'
    program.write_text(
        "program probe\n"
        "  implicit none\n"
        "  logical :: isnonfinite_d\n"
        "  double precision :: x(5)\n"
        "  read (*, *) x\n"
        "  print '(5l2)', isnonfinite_d(x(1)), isnonfinite_d(x(2)), isnonfinite_d(x(3)),"
        " isnonfinite_d(x(4)), isnonfinite_d(x(5))\n"
        "end program probe\n")
    exe = tmp_path / 'probe'
    subprocess.run([gfortran, '-O2', '-ffast-math', '-ffixed-line-length-none',
                    '-o', str(exe), str(program), str(template)], check=True,
                   capture_output=True, text=True)
    # Values arrive at run time, so the compiler cannot fold them.
    out = subprocess.run([str(exe)], input='NaN Inf -Inf 1.5 1.7976931348623157d308\n',
                         capture_output=True, text=True, check=True)
    assert out.stdout.split() == ['T', 'T', 'T', 'F', 'F']
