"""Run the public verification bundles through Python and compiled f2py calls.

No Abaqus installation is used. Each check runs in a fresh process so the
example-local ``build`` and ``check_reference`` modules cannot cross-contaminate
other examples. Missing tools, failed assertions, and timeouts fail the run.
"""

import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
REPOSITORY = "https://github.com/tengzhang48/abaqus_ufl"
CASES = json.loads((ROOT / "tools/livebench_cases.json").read_text(encoding="utf-8"))
CHECK_LABELS = {
    "check_reference.py": "Independent reference",
    "check_assembled.py": "Assembled element",
    "check_compiled.py": "Generated Fortran / f2py",
    "check_fe_runtime.py": "Mesh solves / f2py",
}

# Call the existing check function rather than scraping rounded numbers from
# stdout. Numerical metrics come directly from its return value; assertions and
# the example's numerical tolerances remain the acceptance gate.
CHECK_COMMAND = """
import inspect
import json
from pathlib import Path
import runpy
import sys
sys.path.insert(0, str(Path(sys.argv[1]).parent))
namespace = runpy.run_path(sys.argv[1])
plots = []
check = namespace['check']
result = check(record=plots.append) if 'record' in inspect.signature(check).parameters else check()
metrics = result if isinstance(result, dict) else {}
Path(sys.argv[2]).write_text(json.dumps({'metrics': metrics, 'plots': plots}, allow_nan=False), encoding='utf-8')
"""


def _utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _command_line(command):
    try:
        result = subprocess.run(command, cwd=str(ROOT), capture_output=True,
                                text=True, timeout=10, check=False)
        return result.stdout.strip() if result.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        return None


def _environment():
    versions = {}
    for name in ("numpy", "sympy", "scipy", "meson", "ninja", "matplotlib"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    compiler = _command_line(["gfortran", "--version"])
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "gfortran": compiler.splitlines()[0] if compiler else None,
        "packages": versions,
    }


def _source():
    commit = _command_line(["git", "rev-parse", "HEAD"])
    changes = _command_line(["git", "status", "--porcelain"])
    return {
        "repository": REPOSITORY,
        "commit": commit,
        "dirty": bool(changes) if changes is not None else None,
        "ref": os.environ.get("GITHUB_REF"),
    }


def _run_url():
    repository = os.environ.get("GITHUB_REPOSITORY")
    run_id = os.environ.get("GITHUB_RUN_ID")
    if repository and run_id:
        return "{}/{}/actions/runs/{}".format(
            os.environ.get("GITHUB_SERVER_URL", "https://github.com"),
            repository, run_id)
    return None


def _validate_metrics(metrics):
    if not isinstance(metrics, dict):
        raise ValueError("check metrics must be an object")
    for name, value in metrics.items():
        if (not isinstance(name, str) or isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)):
            raise ValueError("invalid numerical metric: {}".format(name))
    return metrics


def run_check(case_id, script, output, timeout, source=None):
    started = time.perf_counter()
    log_path = Path("logs") / (case_id + "-" + Path(script).stem + ".log")
    check = {
        "script": script,
        "label": CHECK_LABELS[script],
        "status": "failed",
        "exit_code": None,
        "duration_seconds": 0.0,
        "metrics": {},
        "plots": [],
        "log": log_path.as_posix(),
    }
    path = ROOT / CASES[case_id]["directory"] / script
    check["source_path"] = path.relative_to(ROOT).as_posix()
    log = ""
    try:
        if not path.is_file():
            raise FileNotFoundError("Missing required check: " + str(path))
        with tempfile.TemporaryDirectory(prefix="abaqus_ufl_livebench_") as temp:
            metrics_path = Path(temp) / "metrics.json"
            process = subprocess.Popen(
                [sys.executable, "-c", CHECK_COMMAND, str(path), str(metrics_path)],
                cwd=str(path.parent), stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                errors="replace", start_new_session=os.name == "posix",
                env={**os.environ, "PYTHONOPTIMIZE": "0", "PYTHONUTF8": "1"},
            )
            try:
                log, _ = process.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                # Also stop any compiler/build process started by the check.
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
                log, _ = process.communicate()
                check["status"] = "timed_out"
                log += "\n[TIMEOUT] check exceeded {} seconds\n".format(timeout)
            check["exit_code"] = process.returncode
            if process.returncode == 0 and check["status"] != "timed_out":
                result = json.loads(metrics_path.read_text(encoding="utf-8"))
                check["metrics"] = _validate_metrics(result["metrics"])
                check["plots"] = result["plots"]
                if check["plots"]:
                    from tools.livebench_figures import write_figures
                    write_figures(check["plots"], output, case_id, CASES[case_id]["title"], source or _source())
                check["status"] = "passed"
    except (OSError, ValueError, RuntimeError, ImportError, KeyError, TypeError) as error:
        log += "\n[ERROR] {}: {}\n".format(type(error).__name__, error)
    check["duration_seconds"] = round(time.perf_counter() - started, 3)
    (output / log_path).write_text(log, encoding="utf-8")
    return check


def run_benchmarks(case_ids, output, timeout=300.0):
    output = Path(output).resolve()
    (output / "logs").mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    report = {
        "schema_version": 1,
        "started_at": _utc_now(),
        "source": _source(),
        "environment": _environment(),
        "run_url": _run_url(),
        "scope": "Python reference checks and generated Fortran via f2py; no Abaqus solves or paper-scale reproduction.",
        "cases": [],
    }
    for case_id in case_ids:
        case = CASES[case_id]
        scripts = case["scripts"]
        checks = []
        for script in scripts:
            check = run_check(case_id, script, output, timeout, source=report["source"])
            checks.append(check)
            print("[{}] {} / {} ({:.2f}s)".format(
                check["status"].upper(), case_id, script,
                check["duration_seconds"]), flush=True)
        report["cases"].append({
            "id": case_id,
            "title": case["title"],
            "target": case["target"],
            "status": "passed" if all(c["status"] == "passed" for c in checks) else "failed",
            "checks": checks,
        })
    report["finished_at"] = _utc_now()
    report["duration_seconds"] = round(time.perf_counter() - started, 3)
    report["status"] = "passed" if all(c["status"] == "passed" for c in report["cases"]) else "failed"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                     dir=str(output), delete=False) as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
        temporary_path = Path(stream.name)
    temporary_path.replace(output / "report.json")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(CASES), action="append",
                        help="Run one case (repeatable); the default runs all seven.")
    parser.add_argument("--output", type=Path, default=ROOT / "benchmark-results")
    parser.add_argument("--timeout", type=float, default=300.0,
                        help="Maximum seconds per check (default: 300).")
    args = parser.parse_args(argv)
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be positive and finite")
    case_ids = list(dict.fromkeys(args.case or CASES))
    report = run_benchmarks(case_ids, args.output, args.timeout)
    print("Report: " + str(args.output.resolve() / "report.json"))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
