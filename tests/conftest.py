"""Suite-wide policy hooks.

Compiled tests skip when the local f2py toolchain is unusable, so a missing
compiler does not block development. In CI that would silently remove a gate:
on 2026-07-30 f2py was skipped after a setuptools conflict
(docs/ABAQUS_VALIDATION_2026-07-30.md). CI therefore sets
ABAQUS_UFL_FAIL_ON_SKIP=1, which turns every skip into a failure.
"""

import os

import pytest


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if report.skipped and not hasattr(report, 'wasxfail') \
            and os.environ.get('ABAQUS_UFL_FAIL_ON_SKIP') == '1':
        reason = report.longrepr[2] if isinstance(report.longrepr, tuple) else report.longrepr
        report.outcome = 'failed'
        report.longrepr = f'skipped while ABAQUS_UFL_FAIL_ON_SKIP=1: {reason}'
