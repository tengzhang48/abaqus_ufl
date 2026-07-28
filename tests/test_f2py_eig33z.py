"""Compiled gates for the near-diagonal ``eig33z`` CS fallback."""

import os
import subprocess
import sys

import numpy as np
import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
WRAPPER_DIR = os.path.join(HERE, 'f2py_wrappers')
TENSOR_OPS = os.path.join(
    REPO, 'abaqus_ufl', 'generators', 'templates', 'tensor_ops.for')


@pytest.fixture(scope='session')
def eig33z_module():
    module_name = 'eig33z_f2py'
    cmd = [
        sys.executable,
        '-m',
        'numpy.f2py',
        '-c',
        os.path.join(WRAPPER_DIR, 'drive_eig33z.f90'),
        TENSOR_OPS,
        '-m',
        module_name,
        'only:',
        'drive_sqrtm33z_cs',
        'drive_sqrtm33z_direction_cs',
        'drive_sqrtm33z_value',
        ':',
    ]
    result = subprocess.run(
        cmd,
        cwd=WRAPPER_DIR,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        combined = f'{result.stdout}\n{result.stderr}'
        if any(sig in combined for sig in (
                'Compiler.__init__()',
                "No module named 'distutils'",
                "No module named 'numpy.distutils'")):
            pytest.skip(
                'f2py build backend unavailable in this environment; '
                f'toolchain signature:\n{combined[-600:]}')
        raise RuntimeError(
            f'f2py eig33z build failed:\n{combined}')

    if WRAPPER_DIR not in sys.path:
        sys.path.insert(0, WRAPPER_DIR)
    if module_name in sys.modules:
        del sys.modules[module_name]
    return __import__(module_name)


@pytest.mark.parametrize(
    'diag, ii, jj',
    [
        ([1.0, 1.0, 1.0], 0, 1),
        ([4.0, 4.0, 4.0], 0, 2),
        ([1.0, 2.0, 3.0], 0, 1),
        ([1.0, 2.0, 3.0], 1, 2),
        ([2.0, 2.0, 3.0], 0, 1),
        ([2.0, 2.0, 3.0], 0, 2),
        ([2.0, 2.0 + 1.0e-8, 3.0], 0, 1),
        ([2.0, 2.0, 2.0 + 1.0e-8], 0, 1),
    ],
)
def test_compiled_eig_sqrtm_divided_difference(
        eig33z_module, diag, ii, jj):
    """Compiled spectral derivatives survive identity/diagonal degeneracy."""
    diag = np.asarray(diag, dtype=float)
    ds = eig33z_module.drive_sqrtm33z_cs(
        np.diag(diag), ii + 1, jj + 1, 1.0e-20)
    exact = 1.0 / (np.sqrt(diag[ii]) + np.sqrt(diag[jj]))

    assert ds[ii, jj] == pytest.approx(exact, rel=2e-12, abs=2e-14)
    assert ds[jj, ii] == pytest.approx(exact, rel=2e-12, abs=2e-14)
    assert np.all(np.isfinite(ds))


def test_compiled_eig_sqrtm_dense_direction_at_triple_eigenvalue(
        eig33z_module):
    """The three-vector Jacobi branch retains a dense CS direction."""
    value = 4.0
    direction = np.array([
        [0.30, -0.20, 0.15],
        [-0.20, -0.10, 0.25],
        [0.15, 0.25, 0.40],
    ])
    ds = eig33z_module.drive_sqrtm33z_direction_cs(
        value * np.eye(3), direction, 1.0e-20)

    exact = direction / (2.0 * np.sqrt(value))
    np.testing.assert_allclose(ds, exact, rtol=3e-12, atol=3e-14)


def test_compiled_eig_sqrtm_dense_direction_at_near_triple_spectrum(
        eig33z_module):
    """A small real split must not hide a dense derivative direction."""
    diag = np.array([4.0, 4.0 + 1.0e-8, 4.0 + 3.0e-8])
    direction = np.array([
        [0.30, -0.20, 0.15],
        [-0.20, -0.10, 0.25],
        [0.15, 0.25, 0.40],
    ])
    ds = eig33z_module.drive_sqrtm33z_direction_cs(
        np.diag(diag), direction, 1.0e-20)
    exact = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            exact[i, j] = (
                direction[i, j]
                / (np.sqrt(diag[i]) + np.sqrt(diag[j]))
            )

    # Inside the declared 1e-6 quasi-degenerate band, the fallback uses the
    # limiting eigenspace derivative. Its error scales with the tiny real
    # split; this tolerance still detects the former zero-direction guard.
    np.testing.assert_allclose(ds, exact, rtol=3e-8, atol=3e-12)


@pytest.mark.parametrize('gap', [0.0, 1.0e-8])
def test_compiled_eig_sqrtm_dense_direction_at_rotated_repeated_spectrum(
        eig33z_module, gap):
    """The compiled CS eig path must not depend on a diagonal real basis."""
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.71
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    Q = (
        np.eye(3)
        + np.sin(angle) * K
        + (1.0 - np.cos(angle)) * (K @ K)
    )
    diag = np.array([2.0, 2.0 + gap, 3.0])
    A0 = Q @ np.diag(diag) @ Q.T
    direction = np.array([
        [0.30, -0.20, 0.15],
        [-0.20, -0.10, 0.25],
        [0.15, 0.25, 0.40],
    ])
    got = eig33z_module.drive_sqrtm33z_direction_cs(
        A0, direction, 1.0e-20)

    h_principal = Q.T @ direction @ Q
    exact_principal = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            exact_principal[i, j] = (
                h_principal[i, j]
                / (np.sqrt(diag[i]) + np.sqrt(diag[j]))
            )
    exact = Q @ exact_principal @ Q.T
    np.testing.assert_allclose(got, exact, rtol=3e-8, atol=3e-10)


@pytest.mark.parametrize('gap', [0.0, 1.0e-12, 1.0e-8])
def test_compiled_eig_sqrtm_value_at_rotated_repeated_spectrum(
        eig33z_module, gap):
    """The real eig value path must remain nonsingular after rotation."""
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.71
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    Q = (
        np.eye(3)
        + np.sin(angle) * K
        + (1.0 - np.cos(angle)) * (K @ K)
    )
    diag = np.array([2.0, 2.0 + gap, 3.0])
    A0 = Q @ np.diag(diag) @ Q.T
    got = eig33z_module.drive_sqrtm33z_value(A0)
    exact = Q @ np.diag(np.sqrt(diag)) @ Q.T
    np.testing.assert_allclose(got, exact, rtol=3e-12, atol=3e-12)


@pytest.mark.parametrize(
    'scale',
    [1.0, 1.0e-9, 1.0e-12, 1.0e-13, 1.0e-20, 1.0e-30, 1.0e-31],
)
def test_compiled_eig_sqrtm_value_at_scaled_rotated_distinct_spectrum(
        eig33z_module, scale):
    """Compiled value guards must remain correct at very small scales."""
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.71
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    Q = (
        np.eye(3)
        + np.sin(angle) * K
        + (1.0 - np.cos(angle)) * (K @ K)
    )
    diag = scale * np.array([0.1, 1.0, 10.0])
    A0 = Q @ np.diag(diag) @ Q.T
    got = eig33z_module.drive_sqrtm33z_value(A0)
    exact = Q @ np.diag(np.sqrt(diag)) @ Q.T
    rel = np.linalg.norm(got - exact) / np.linalg.norm(exact)
    assert rel < 2.0e-12
