"""The CI skip policy in tests/conftest.py (audit 2026-10-10: G9).

Runs pytest on small probe files, with and without ABAQUS_UFL_FAIL_ON_SKIP=1,
using the real conftest.py.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

CONFTEST = Path(__file__).with_name('conftest.py')

TEST_SKIPS = """
import pytest

def test_skip_in_body():
    pytest.skip('toolchain probe')

@pytest.mark.skipif(True, reason='marker probe')
def test_skip_marker():
    pass

@pytest.mark.xfail(reason='expected failure')
def test_expected_failure():
    assert False
"""

MODULE_SKIP = """
import pytest
pytest.skip('module probe', allow_module_level=True)

def test_never_collected():
    pass
"""


def run_probe(tmp_path, source, fail_on_skip):
    shutil.copy(CONFTEST, tmp_path / 'conftest.py')
    (tmp_path / 'test_probe.py').write_text(source)
    env = {k: v for k, v in os.environ.items() if k != 'ABAQUS_UFL_FAIL_ON_SKIP'}
    if fail_on_skip:
        env['ABAQUS_UFL_FAIL_ON_SKIP'] = '1'
    result = subprocess.run(
        [sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider', str(tmp_path)],
        cwd=str(tmp_path), capture_output=True, text=True, env=env, check=False)
    return result.returncode, result.stdout.strip().splitlines()[-1]


def test_skips_stay_skips_without_the_variable(tmp_path):
    code, summary = run_probe(tmp_path, TEST_SKIPS, fail_on_skip=False)
    assert code == 0, summary
    assert '2 skipped' in summary and '1 xfailed' in summary


def test_skips_fail_with_the_variable_but_xfail_does_not(tmp_path):
    code, summary = run_probe(tmp_path, TEST_SKIPS, fail_on_skip=True)
    assert code != 0, summary
    assert '1 failed' in summary and '1 error' in summary and '1 xfailed' in summary
    assert 'skipped' not in summary


def test_module_level_skip_fails_with_the_variable(tmp_path):
    code, summary = run_probe(tmp_path, MODULE_SKIP, fail_on_skip=False)
    assert code == 5 and '1 skipped' in summary, summary    # 5: no tests ran
    code, summary = run_probe(tmp_path, MODULE_SKIP, fail_on_skip=True)
    assert code != 0 and 'error' in summary, summary


def test_absent_optional_runtime_can_skip_in_the_mpi_free_ci_job(tmp_path):
    source = "import pytest\ndef test_optional():\n    pytest.skip('optional runtime: petsc4py is not installed')\n"
    code, summary = run_probe(tmp_path, source, fail_on_skip=True)
    assert code == 0 and '1 skipped' in summary, summary
