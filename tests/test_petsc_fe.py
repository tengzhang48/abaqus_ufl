"""PETSc is optional; when present, solves must match independent serial checks."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest
from scipy.sparse import csr_matrix, diags

from abaqus_ufl.fe import (make_fieldsplit_solver, make_petsc_solver,
                           partition_rcb, partition_structured_quad,
                           partition_unstructured, structured_quad_mesh)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def petsc():
    pytest.importorskip('petsc4py', reason='optional runtime: petsc4py is not installed')
    from petsc4py import PETSc
    return PETSc


def test_serial_package_and_mpi_modules_import_without_optional_libraries():
    code = '''
import importlib.abc
import sys
class Missing(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('petsc4py', 'mpi4py'):
            raise ModuleNotFoundError(fullname)
sys.meta_path.insert(0, Missing())
import abaqus_ufl
import abaqus_ufl.fe
import abaqus_ufl.fe.mpi_batch
import abaqus_ufl.fe.petsc_mpi
assert 'petsc4py' not in sys.modules and 'mpi4py' not in sys.modules
try:
    abaqus_ufl.fe.make_petsc_solver()
except ImportError as error:
    assert 'petsc4py' in str(error)
else:
    raise AssertionError('missing dependency was not reported')
'''
    result = subprocess.run([sys.executable, '-c', code], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('pc,solver', [('lu', 'mumps'), ('gamg', None)])
def test_petsc_matches_scipy_and_rejects_bad_inputs(petsc, pc, solver):
    K = diags([-np.ones(7), 3*np.ones(8), -0.5*np.ones(7)], [-1, 0, 1]).tocsr()
    exact = np.arange(8.) / 7
    solve = make_petsc_solver(pc=pc, solver=solver)
    assert np.allclose(solve(K, K @ exact), exact, atol=1e-10, rtol=1e-10)
    with pytest.raises(ValueError, match='finite'):
        solve(K, np.full(8, np.nan))
    with pytest.raises(ValueError, match='matching'):
        solve(K, np.ones(7))


@pytest.mark.parametrize('split', ['additive', 'multiplicative', 'schur'])
def test_fieldsplit_on_nonsymmetric_coupled_system(petsc, split):
    K = np.diag([3., 5., 4., 6., 5., 7.])
    K += 0.1 * np.triu(np.ones((6, 6)), 1)
    exact = np.linspace(-0.2, 0.7, 6)
    solve = make_fieldsplit_solver([('u', [0]), ('q', [1])], 2,
                                    split_type=split, sub_pc='lu')
    assert np.allclose(solve(csr_matrix(K), K @ exact), exact, atol=1e-10, rtol=1e-10)


def test_ksp_failure_is_not_a_solution(petsc):
    K = diags(np.arange(1., 21.)).tocsr()
    solve = make_petsc_solver(pc='none', ksp='gmres', maxit=1, rtol=1e-15)
    with pytest.raises(RuntimeError, match='KSP failed'):
        solve(K, np.ones(20))


def test_partitions_reconstruct_mesh_and_handle_empty_ranks():
    nodes, elems = structured_quad_mesh(3, 2)
    all_gm = (elems[:, :, None] * 2 + np.arange(2)).reshape(6, 8)
    partitions = [partition_structured_quad(3, 2, 2, rank, 8) for rank in range(8)]
    assert all(ndof == len(nodes)*2 for ndof, _, _ in partitions)
    assert np.array_equal(np.concatenate([p[1] for p in partitions]), all_gm)
    assert np.allclose(np.concatenate([p[2] for p in partitions]), nodes[elems])
    owners = partition_rcb(nodes[elems].mean(axis=1), 8)
    assert np.all((0 <= owners) & (owners < 8))
    local = [partition_unstructured(nodes, elems, 2, 8, rank, owners) for rank in range(8)]
    assert sum(p[3] for p in local) == len(nodes)*2
    assert sum(len(p[1]) for p in local) == len(elems)
    for ndof, gm, coords, owned, permutation in local:
        assert np.array_equal(np.sort(permutation), np.arange(len(nodes)))
        assert gm.shape == (len(coords), 8)
        assert 0 <= owned <= ndof


def _mpiexec(petsc, script, args, timeout=90):
    launcher = shutil.which('mpiexec')
    if launcher is None:
        pytest.skip('optional runtime: mpiexec is not installed')
    # This subprocess is also the MPI ABI test: failures/timeouts must fail.
    result = subprocess.run([launcher, '-n', '2', sys.executable, str(script), *map(str, args)],
                            cwd=ROOT, capture_output=True, text=True, timeout=timeout,
                            env={**os.environ, 'PYTHONUTF8': '1'})
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def test_parallel_compiled_example_matches_serial_and_oracles(petsc, tmp_path):
    if shutil.which('gfortran') is None:
        pytest.skip('gfortran is required for compiled UEL tests')
    output = _mpiexec(petsc, ROOT / 'examples/scalar_diffusion_uel/check_parallel.py',
                      ['--build-dir', tmp_path])
    report = json.loads(next(line for line in output.splitlines() if line.startswith('{')))
    assert report['ranks'] == 2
    assert max(report['errors'].values()) < 1e-10
    assert report['cache_hits'] > 0
    assert report['ghost_dofs'] > 0


def test_parallel_rejections_empty_partitions_and_cache_consensus(petsc):
    output = _mpiexec(petsc, ROOT / 'tests/mpi_fe_contract.py', [])
    assert 'MPI contracts passed' in output
