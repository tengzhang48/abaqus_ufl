"""The agent guidance files keep working links, anchors and named tests (tools/check_guidance.py)."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("check_guidance", ROOT / "tools" / "check_guidance.py")
cg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cg)


def test_guidance_links_anchors_and_named_tests_resolve():
    nlinks, problems = cg.check()
    assert nlinks > 0
    assert not problems, "\n".join(problems)


def test_slug_matches_github_anchors():
    # Expected values are GitHub's anchors for these headings.
    assert cg.slug("Comments use ** not bare *") == "comments-use--not-bare-"
    assert cg.slug("Perturbation size (`CS_H`)") == "perturbation-size-cs_h"
    assert cg.slug("`_scratch` arrays") == "_scratch-arrays"
    assert cg.slug("F-bar: contract over in-plane / out-of-plane") == "f-bar-contract-over-in-plane--out-of-plane"
    assert cg.slug("`constants=N` must match the UMAT contract") == "constantsn-must-match-the-umat-contract"


def test_checker_rejects_broken_links_and_handles_duplicates_and_fences(tmp_path):
    lessons = tmp_path / "docs" / "lessons"
    lessons.mkdir(parents=True)
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_x.py").write_text("def test_real():\n    pass\n")
    (lessons / "a.md").write_text(
        "## Same\n\n## Same\n\n~~~~\n## Inside a fence\n~~~\n~~~~\n\n<a id=\"old-name\"></a>\n")
    (lessons / "b.md").write_text(
        "[ok](a.md#same) [ok](a.md#same-1) [ok](a.md#old-name)\n"
        "[bad](a.md#same-2) [bad](a.md#inside-a-fence) [bad](missing.md)\n"
        "tests/test_x.py::test_real tests/test_x.py::test_missing\n")
    nlinks, problems = cg.check(str(tmp_path))
    assert nlinks == 6
    assert sorted(problems) == sorted([
        "docs/lessons/b.md: missing anchor a.md#same-2",
        "docs/lessons/b.md: missing anchor a.md#inside-a-fence",
        "docs/lessons/b.md: missing file missing.md",
        "docs/lessons/b.md: missing test tests/test_x.py::test_missing",
    ])
