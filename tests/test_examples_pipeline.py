"""Collect every example verification bundle into the pytest run.

Each example ships ``check_reference.py`` (independent closed-form
oracle), ``check_compiled.py`` (deterministic regeneration + gfortran +
f2py execution), and, for the UELs, ``check_assembled.py`` (assembled
residual/tangent oracles). They are self-contained scripts so users can
run them directly from the example directory; this module runs the same
scripts as subprocesses so that plain ``pytest tests`` (and CI) covers
the full license-free pipeline of every shipped example.

The compiled gates need gfortran plus a working ``numpy.f2py`` backend
(meson/ninja on NumPy 2.x); environments without that toolchain skip
those gates rather than fail.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"

# Explicit release manifest: every listed example MUST ship the listed
# checks; a new tracked example directory must be added here or the suite
# fails. This keeps the pipeline fail-closed (no silent skips for
# missing evidence scripts).
MANIFEST = {
    "neo_hookean_umat": {"assembled": False},
    "small_strain_j2_umat": {"assembled": False},
    "small_strain_viscoelastic_umat": {"assembled": False},
    "ogden_umat": {"assembled": False},
    "scalar_diffusion_uel": {"assembled": True},
    "thermo_mechanics_quad8": {"assembled": True},
}
BUNDLES = sorted(MANIFEST)


def test_manifest_matches_tracked_examples():
    on_disk = sorted(
        d.name for d in EXAMPLES.iterdir()
        if d.is_dir() and not d.name.startswith("_")
        and not d.name.startswith(".")
    )
    assert on_disk == BUNDLES, (
        "example directories and release manifest disagree; update "
        "MANIFEST in this file: on disk {} vs manifest {}".format(
            on_disk, BUNDLES))


@pytest.mark.parametrize("example", BUNDLES)
def test_required_evidence_scripts_exist(example):
    required = ["build.py", "check_reference.py", "check_compiled.py",
                "README.md"]
    if MANIFEST[example]["assembled"]:
        required.append("check_assembled.py")
    missing = [f for f in required if not (EXAMPLES / example / f).exists()]
    assert not missing, "{} is missing {}".format(example, missing)

TOOLCHAIN_SKIP_SIGNATURES = (
    "Compiler.__init__()",
    "No module named 'distutils'",
    "No module named 'numpy.distutils'",
)


def _run_check(example, script):
    path = EXAMPLES / example / script
    assert path.exists(), "{} is missing {}".format(example, script)
    result = subprocess.run(
        [sys.executable, script],
        cwd=str(EXAMPLES / example),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        combined = result.stdout + "\n" + result.stderr
        if any(sig in combined for sig in TOOLCHAIN_SKIP_SIGNATURES):
            pytest.skip(
                "f2py build backend unavailable; toolchain signature:\n"
                + combined[-600:])
        raise AssertionError(
            "{} {} failed:\n{}".format(example, script, combined))


@pytest.mark.parametrize("example", BUNDLES)
def test_reference_oracle(example):
    _run_check(example, "check_reference.py")


@pytest.mark.parametrize(
    "example", [e for e in BUNDLES if MANIFEST[e]["assembled"]])
def test_assembled_oracle(example):
    _run_check(example, "check_assembled.py")


@pytest.mark.parametrize("example", BUNDLES)
def test_compiled_pipeline(example):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    _run_check(example, "check_compiled.py")
