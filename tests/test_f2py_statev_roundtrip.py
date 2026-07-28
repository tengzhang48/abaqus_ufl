"""Compiled nonsymmetric-tensor STATEV round-trip sentinel.

Generates a finite-strain UMAT whose state variable is a genuinely
NONSYMMETRIC 3x3 tensor, compiles it with f2py + ``drive_umat.f90``, and
checks two things a symmetric or diagonal test state cannot see:

1. **Layout** — the STATEV slots hold the tensor in column-major
   (Fortran) order: ``STATEV[3*j + i] == A[i, j]``, i.e. the layout the
   generator emits and the layout ``A.flatten(order='F')`` produces on
   the Python side.
2. **Round trip** — a second UMAT call with ``TIME(2) > 0`` reads the
   state back and continues the recursion ``A_new = A_old @ F`` exactly.
   A row/column-major mismatch would transpose ``A`` between steps; the
   test asserts the correct product AND that the transposed-read result
   is far away, so the gate genuinely discriminates.

Build dependencies: gfortran + a working ``numpy.f2py`` backend
(meson/ninja on NumPy 2.x). Environments with a broken f2py toolchain
skip rather than fail.
"""
import os
import subprocess
import sys
from collections import OrderedDict

import numpy as np
import pytest

from abaqus_ufl import Material
from abaqus_ufl.core.tensor import det, inv, log
from abaqus_ufl.generators.umat_gen import generate_umat

HERE = os.path.dirname(os.path.abspath(__file__))
WRAPPER_DIR = os.path.join(HERE, 'f2py_wrappers')


class FabricNeoHookean(Material):
    """Compressible neo-Hookean stress + a nonsymmetric tensor state.

    The state update ``A_new = A_old @ F`` is deliberately order
    sensitive: for a shear ``F``, ``A_old.T @ F != A_old @ F``, so a
    transposed STATEV read cannot masquerade as correct.
    """

    props = OrderedDict([
        ('mu', 10.0),
        ('kappa', 100.0),
    ])

    state_vars = OrderedDict([
        ('A', np.eye(3)),
    ])

    def stress_PK1(self, F, A_old, dt):
        J = det(F)
        Finv = inv(F)
        FinvT = Finv.T
        P = self.mu * (F - FinvT) + self.kappa * log(J) * FinvT
        A_new = A_old @ F
        return P, {'A': A_new}


@pytest.fixture(scope='module')
def umat_module(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('statev_roundtrip')
    umat_path = os.path.join(str(tmp), 'fabric_umat.for')
    generate_umat(FabricNeoHookean(), umat_path)

    module_name = 'statev_roundtrip_f2py'
    cmd = [
        sys.executable, '-m', 'numpy.f2py',
        '-c',
        os.path.join(WRAPPER_DIR, 'drive_umat.f90'),
        umat_path,
        '-m', module_name,
        'only:', 'drive_umat', ':',
    ]
    result = subprocess.run(
        cmd, cwd=str(tmp), capture_output=True, text=True, check=False)
    if result.returncode != 0:
        combined = f'{result.stdout}\n{result.stderr}'
        if any(sig in combined for sig in (
                'Compiler.__init__()',
                "No module named 'distutils'",
                "No module named 'numpy.distutils'")):
            pytest.skip(
                'f2py build backend unavailable in this environment; '
                f'toolchain signature:\n{combined[-600:]}')
        raise RuntimeError(f'f2py statev sentinel build failed:\n{combined}')

    if str(tmp) not in sys.path:
        sys.path.insert(0, str(tmp))
    if module_name in sys.modules:
        del sys.modules[module_name]
    return __import__(module_name)


def _call_umat(module, dfgrd1, statev, time):
    dstran = np.zeros(6, dtype=float)
    stran = np.zeros(6, dtype=float)
    props = np.array([10.0, 100.0], dtype=float)
    predef = np.array([0.0], dtype=float)
    dpred = np.array([0.0], dtype=float)
    coords = np.zeros(3, dtype=float)
    drot = np.eye(3, dtype=float)
    cmname = b'fabric'.ljust(80)
    stress, statev_out, ddsdde, pnewdt = module.drive_umat(
        statev, dfgrd1, np.eye(3, dtype=float), dstran, stran, props,
        np.asarray(time, dtype=float), 1.0, 0.0, 0.0, predef, dpred, 1.0,
        coords, drot, 3, 3, cmname)
    return stress, statev_out, ddsdde, pnewdt


# Shear + unequal stretches: F1 and A stay nonsymmetric, J != 1.
F1 = np.array([
    [1.10, 0.30, 0.00],
    [0.00, 1.00, 0.20],
    [0.00, 0.00, 0.95],
])
F2 = np.array([
    [1.00, 0.00, 0.10],
    [0.25, 1.05, 0.00],
    [0.00, 0.00, 1.00],
])


def test_statev_layout_is_column_major(umat_module):
    """After the TIME(2)==0 init step, STATEV[3*j+i] == A_new[i, j]."""
    statev = np.zeros(9, dtype=float)
    _, statev_out, _, _ = _call_umat(
        umat_module, F1, statev, time=(0.0, 0.0))

    A1 = np.eye(3) @ F1
    assert np.max(np.abs(A1 - A1.T)) > 0.1, (
        'test state degenerated to symmetric; the layout gate would be blind')

    expected = A1.flatten(order='F')
    np.testing.assert_allclose(statev_out, expected, rtol=0, atol=1e-12)
    # The same data in row-major order must NOT match.
    wrong = A1.flatten(order='C')
    assert not np.allclose(statev_out, wrong, atol=1e-6)


def test_statev_round_trip_two_steps(umat_module):
    """Step 2 must read back A1 (not A1.T) and produce A1 @ F2."""
    statev = np.zeros(9, dtype=float)
    _, statev_1, _, _ = _call_umat(
        umat_module, F1, statev, time=(0.0, 0.0))

    _, statev_2, _, _ = _call_umat(
        umat_module, F2, statev_1.copy(), time=(1.0, 1.0))

    A1 = np.eye(3) @ F1
    A2 = A1 @ F2
    A2_transposed_read = A1.T @ F2

    assert np.max(np.abs(A2 - A2_transposed_read)) > 1e-3, (
        'discriminator collapsed: transposed read would be indistinguishable')

    np.testing.assert_allclose(
        statev_2, A2.flatten(order='F'), rtol=0, atol=1e-12)
    assert not np.allclose(
        statev_2, A2_transposed_read.flatten(order='F'), atol=1e-6)


def test_sdvini_state_is_honored_at_time_zero(umat_module):
    """Nonzero incoming STATEV at TIME(2)==0 must be READ, not reset.

    Abaqus users can initialize state via SDVINI (or *INITIAL CONDITIONS,
    TYPE=SOLUTION), which delivers a nonzero STATEV on the very first
    increment. The generated guard may fall back to the declared
    initializers only when TIME(2)==0 AND the incoming STATEV is all
    zeros; a nonzero incoming state must win.
    """
    A_in = np.array([
        [1.00, 0.40, 0.00],
        [0.10, 1.00, 0.00],
        [0.00, 0.05, 1.00],
    ])
    statev = A_in.flatten(order='F').copy()
    _, statev_out, _, _ = _call_umat(
        umat_module, F1, statev, time=(0.0, 0.0))

    honored = A_in @ F1          # guard read the SDVINI-provided state
    reset = np.eye(3) @ F1       # guard clobbered it with the declared init
    assert np.max(np.abs(honored - reset)) > 1e-2, (
        'discriminator collapsed: reset would be indistinguishable')

    np.testing.assert_allclose(
        statev_out, honored.flatten(order='F'), rtol=0, atol=1e-12)
    assert not np.allclose(statev_out, reset.flatten(order='F'), atol=1e-6)


def test_all_zero_statev_at_time_zero_uses_declared_init(umat_module):
    """All-zero STATEV at TIME(2)==0 falls back to the declared A=I."""
    statev = np.zeros(9, dtype=float)
    _, statev_out, _, _ = _call_umat(
        umat_module, F1, statev, time=(0.0, 0.0))
    np.testing.assert_allclose(
        statev_out, (np.eye(3) @ F1).flatten(order='F'), rtol=0, atol=1e-12)


def test_stress_finite_and_tangent_shape(umat_module):
    statev = np.zeros(9, dtype=float)
    stress, _, ddsdde, pnewdt = _call_umat(
        umat_module, F1, statev, time=(0.0, 0.0))
    assert np.all(np.isfinite(stress))
    assert ddsdde.shape == (6, 6)
    assert np.all(np.isfinite(ddsdde))
    assert pnewdt == pytest.approx(1.0)
