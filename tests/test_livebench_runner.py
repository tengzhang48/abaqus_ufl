"""Check failure reporting with small adversarial scripts, without compilation."""

import json
import runpy

import pytest

from tools import run_livebench as runner


def _case(tmp_path, monkeypatch, body):
    directory = tmp_path / "examples" / "broken"
    directory.mkdir(parents=True)
    (directory / "check_reference.py").write_text(body)
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    monkeypatch.setattr(runner, "CASES", {"broken": {
        "title": "Broken control", "target": "UMAT", "directory": "examples/broken",
        "scripts": ["check_reference.py", "check_compiled.py"],
    }})
    output = tmp_path / "results"
    (output / "logs").mkdir(parents=True)
    return output


def test_failed_and_missing_checks_produce_failed_report(tmp_path, monkeypatch):
    output = _case(tmp_path, monkeypatch, "def check():\n    raise AssertionError('known defect')\n")
    report = runner.run_benchmarks(["broken"], output, timeout=5)
    assert report["status"] == "failed"
    assert all(c["status"] == "failed" for c in report["cases"][0]["checks"])
    logs = [output / c["log"] for c in report["cases"][0]["checks"]]
    assert "known defect" in logs[0].read_text()
    assert "Missing required check" in logs[1].read_text()
    assert json.loads((output / "report.json").read_text())["status"] == "failed"


def test_timeout_is_recorded_as_failure(tmp_path, monkeypatch):
    output = _case(tmp_path, monkeypatch, "import time\ndef check():\n    time.sleep(60)\n")
    result = runner.run_check("broken", "check_reference.py", output, 0.1)
    assert result["status"] == "timed_out"
    assert result["exit_code"] != 0
    assert "TIMEOUT" in (output / result["log"]).read_text()


def test_metric_precision_is_retained_and_nan_is_rejected(tmp_path, monkeypatch):
    output = _case(tmp_path, monkeypatch, "def check():\n    return {'error': 1.234567890123456e-14}\n")
    result = runner.run_check("broken", "check_reference.py", output, 5)
    assert result["status"] == "passed"
    assert result["metrics"]["error"] == 1.234567890123456e-14
    path = tmp_path / "examples/broken/check_reference.py"
    path.write_text("def check():\n    return {'error': float('nan')}\n")
    result = runner.run_check("broken", "check_reference.py", output, 5)
    assert result["status"] == "failed"


def test_response_data_is_recorded_during_the_same_check(tmp_path, monkeypatch):
    output = _case(tmp_path, monkeypatch,
                   "calls = 0\ndef check(*, record=None):\n"
                   "    global calls\n    calls += 1\n"
                   "    record({'kind': 'curve', 'id': 'curve', 'title': 'Response', 'description': 'Test samples', "
                   "'x': [0, 0.1], 'x_label': 'strain', 'y_label': 'stress', 'x_scale': 'linear', 'y_scale': 'linear', "
                   "'series': [{'label': 'Compiled', 'role': 'computed', 'values': [0, 1.234567890123456]}]})\n"
                   "    return {'calls': calls}\n")
    result = runner.run_check("broken", "check_reference.py", output, 5)
    assert result["status"] == "passed"
    assert result["metrics"]["calls"] == 1
    assert result["plots"][0]["series"][0]["values"][1] == 1.234567890123456
    assert (output / result["plots"][0]["figure"]).is_file()


@pytest.mark.parametrize("timeout", ["0", "-1", "nan", "inf"])
def test_invalid_timeout_is_rejected(timeout):
    with pytest.raises(SystemExit) as error:
        runner.main(["--timeout", timeout])
    assert error.value.code == 2


def test_manifest_tracks_the_existing_bundles_and_one_mesh_gate():
    MANIFEST = runpy.run_path(str(runner.ROOT / "tests/test_examples_pipeline.py"))["MANIFEST"]
    assert set(runner.CASES) == set(MANIFEST) | {"serial_fe"}
    for case_id, specification in runner.CASES.items():
        for script in specification["scripts"]:
            assert (runner.ROOT / specification["directory"] / script).is_file()
