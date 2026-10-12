"""Shared helpers for tests that compile generated Fortran through f2py."""

import os
import subprocess
import sys

import numpy as np
import pytest

WRAPPER_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'f2py_wrappers')
TOOLCHAIN_SIGNATURES = (
    'Compiler.__init__()',
    "No module named 'distutils'",
    "No module named 'numpy.distutils'",
)


def build_umat_module(umat_path, workdir, module_name):
    """Compile a generated UMAT with ``drive_umat.f90`` and import it."""
    cmd = [sys.executable, '-m', 'numpy.f2py', '-c',
           os.path.join(WRAPPER_DIR, 'drive_umat.f90'), str(umat_path),
           '-m', module_name, 'only:', 'drive_umat', ':']
    result = subprocess.run(cmd, cwd=str(workdir), capture_output=True,
                            text=True, check=False)
    if result.returncode != 0:
        combined = f'{result.stdout}\n{result.stderr}'
        if any(sig in combined for sig in TOOLCHAIN_SIGNATURES):
            pytest.skip('f2py build backend unavailable; toolchain signature:\n'
                        + combined[-600:])
        raise RuntimeError(f'f2py UMAT build failed:\n{combined}')
    if str(workdir) not in sys.path:
        sys.path.insert(0, str(workdir))
    sys.modules.pop(module_name, None)
    return __import__(module_name)


def call_umat(module, dfgrd1, props, nstatv=1, time=(0.0, 0.0), ntens=6, dstran=None):
    """One UMAT call through ``drive_umat``; returns (stress, statev, ddsdde, pnewdt)."""
    dstran = np.zeros(ntens) if dstran is None else np.asarray(dstran, dtype=float)
    return module.drive_umat(
        np.zeros(nstatv), np.asarray(dfgrd1, dtype=float), np.eye(3),
        dstran, np.zeros(ntens), np.asarray(props, dtype=float),
        np.asarray(time, dtype=float), 1.0, 0.0, 0.0, np.zeros(1), np.zeros(1),
        1.0, np.zeros(3), np.eye(3), 3, ntens - 3, b'test'.ljust(80))


def build_uel_module(uel_path, workdir, module_name):
    """Compile a generated Quad4 UEL with the packaged ``drive_uel.f90``."""
    from abaqus_ufl.fe import compiled_element
    wrapper = os.path.join(os.path.dirname(compiled_element.__file__), 'drive_uel.f90')
    cmd = [sys.executable, '-m', 'numpy.f2py', '-c', wrapper, str(uel_path),
           '-m', module_name, 'only:', 'drive_uel', ':']
    if np.lib.NumpyVersion(np.__version__) >= '1.26.0':
        cmd.extend(['--backend', 'meson'])
    result = subprocess.run(cmd, cwd=str(workdir), capture_output=True,
                            text=True, check=False)
    if result.returncode != 0:
        combined = f'{result.stdout}\n{result.stderr}'
        if any(sig in combined for sig in TOOLCHAIN_SIGNATURES):
            pytest.skip('f2py build backend unavailable; toolchain signature:\n'
                        + combined[-600:])
        raise RuntimeError(f'f2py UEL build failed:\n{combined}')
    if str(workdir) not in sys.path:
        sys.path.insert(0, str(workdir))
    sys.modules.pop(module_name, None)
    return __import__(module_name)


def call_uel(module, coords, u, props, du=None, lflags=(1, 0, 1, 0, 0, 0)):
    """One normal UEL call through ``drive_uel``; returns (rhs, amatrx, svars, pnewdt)."""
    u = np.asarray(u, dtype=float)
    du = np.zeros_like(u) if du is None else np.asarray(du, dtype=float)
    return module.drive_uel(
        np.zeros(1), np.asfortranarray(np.asarray(coords, dtype=float).T), u, du,
        np.asarray(props, dtype=float), np.array([0], dtype=np.int32), np.array([0.0, 0.0]),
        1.0, 1.0, np.array(lflags, dtype=np.int32), np.zeros(3), 1, 1, 1, 1, 1.0)
