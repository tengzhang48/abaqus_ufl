"""Generic release-boundary checks on the public tree.

This public test enforces only NON-identifying rules: no binary or
private-format artifacts may be tracked, and no absolute internal mount
paths may appear in tracked text. The full release-boundary policy,
including the identifier denylist, is maintained and enforced from the
private development repository, which scans this repository's checkout;
keeping that list here would disclose the very names it protects.
"""

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Tracked artifact types that must never ship.
DENY_SUFFIXES = {".pdf", ".odb", ".bundle", ".so", ".o", ".pyc"}

# Absolute internal mount prefixes that must not leak into tracked text.
DENY_PATH_MARKERS = ["/media/volume", "/mnt/project"]


def _tracked_files():
    result = subprocess.run(
        ["git", "ls-files"], cwd=str(ROOT),
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0 or not result.stdout.strip():
        # Installed from an sdist/wheel: no git checkout to audit. The
        # boundary gate is a repository-release check, not a runtime one.
        pytest.skip("release-boundary audit requires a git checkout")
    return [line for line in result.stdout.splitlines() if line]


def _is_texty(path):
    try:
        (ROOT / path).read_text(encoding="utf-8")
        return True
    except (UnicodeDecodeError, FileNotFoundError):
        return False


def test_no_denied_artifact_types_tracked():
    bad = [f for f in _tracked_files()
           if Path(f).suffix.lower() in DENY_SUFFIXES]
    assert not bad, "denied artifact types tracked: {}".format(bad)


SELF = "tests/test_release_boundary.py"


def test_no_internal_mount_paths_in_content():
    violations = []
    for f in _tracked_files():
        if f == SELF:      # this file names the markers it scans for
            continue
        if not _is_texty(f):
            continue
        content = (ROOT / f).read_text(encoding="utf-8")
        for marker in DENY_PATH_MARKERS:
            if marker in content:
                violations.append((f, marker))
    assert not violations, (
        "internal mount paths in tracked content: {}".format(violations))
