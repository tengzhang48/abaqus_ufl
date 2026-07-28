"""Compare extracted Abaqus ODB results against reference values.

Usage:
    python compare_results.py --extracted results/elastic_umat.json \
                              --reference references/elastic_umat.json \
                              --output report.json

Supports:
    - Scalar comparisons with relative tolerance (default 1e-4)
    - Array comparisons element-wise
    - Nested key paths like "fields.S.data.0" or "fields.S.mises"
    - Per-key tolerance overrides

Exit code 0 if all checks pass, 1 if any fail.
"""
import argparse
import json
import sys
from pathlib import Path


def get_nested(obj, path):
    """Get a nested value by dot-separated path. Supports integer indices."""
    parts = path.split(".")
    for part in parts:
        if isinstance(obj, dict):
            obj = obj.get(part)
        elif isinstance(obj, list):
            obj = obj[int(part)]
        else:
            return None
        if obj is None:
            return None
    return obj


def compare_scalar(observed, expected, rtol=1e-4, atol=1e-10):
    """Return (passed, message)."""
    if observed is None or expected is None:
        return False, f"missing value (obs={observed}, exp={expected})"
    if isinstance(observed, (str, bool)) or isinstance(expected, (str, bool)):
        if observed == expected:
            return True, "exact match"
        return False, f"value mismatch (obs={observed!r}, exp={expected!r})"
    denom = max(abs(expected), 1e-15)
    rel_err = abs(observed - expected) / denom
    abs_err = abs(observed - expected)
    if rel_err <= rtol or abs_err <= atol:
        return True, f"rel_err={rel_err:.3e} <= {rtol}"
    return False, f"rel_err={rel_err:.3e} > {rtol} (obs={observed:.12e}, exp={expected:.12e})"


def compare_array(observed, expected, rtol=1e-4, atol=1e-10):
    """Compare arrays element-wise. Return (passed, messages_list)."""
    if observed is None or expected is None:
        return False, ["missing array value"]
    if len(observed) != len(expected):
        return False, [f"length mismatch: {len(observed)} vs {len(expected)}"]
    messages = []
    all_pass = True
    for i, (o, e) in enumerate(zip(observed, expected)):
        p, msg = compare_scalar(o, e, rtol=rtol, atol=atol)
        if not p:
            all_pass = False
            messages.append(f"[{i}]: {msg}")
    if all_pass:
        messages.append("all elements within tolerance")
    return all_pass, messages


def run_comparison(extracted, reference, default_rtol=1e-4):
    """Run all comparisons defined in reference["checks"].

    Reference JSON format:
    {
      "checks": {
        "p": {"path": "fields.S.data.0", "expected": 1.23, "rtol": 1e-4},
        "q": {"path": "fields.S.mises", "expected": 0.45, "rtol": 1e-4}
      }
    }
    """
    checks = reference.get("checks", {})
    if not isinstance(checks, dict) or not checks:
        raise ValueError(
            "reference.json must contain at least one named quantitative check"
        )
    results = {}
    all_pass = True

    for check_name, check_spec in checks.items():
        path = check_spec["path"]
        expected = check_spec["expected"]
        rtol = check_spec.get("rtol", default_rtol)
        atol = check_spec.get("atol", 1e-10)

        observed = get_nested(extracted, path)

        if isinstance(expected, list):
            passed, messages = compare_array(observed, expected, rtol=rtol, atol=atol)
        else:
            passed, msg = compare_scalar(observed, expected, rtol=rtol, atol=atol)
            messages = [msg]

        if not passed:
            all_pass = False

        results[check_name] = {
            "passed": passed,
            "path": path,
            "expected": expected,
            "observed": observed,
            "rtol": rtol,
            "atol": atol,
            "messages": messages
        }

    return all_pass, results


def print_report(all_pass, results):
    print("=" * 70)
    print("Solver Output Comparison Report")
    print("=" * 70)
    for name, res in results.items():
        status = "PASS" if res["passed"] else "FAIL"
        print(f"\n[{status}] {name}")
        print(f"  path     : {res['path']}")
        print(f"  expected : {res['expected']}")
        print(f"  observed : {res['observed']}")
        for msg in res["messages"]:
            print(f"  >> {msg}")
    print("\n" + "=" * 70)
    print(f"Overall: {'PASS' if all_pass else 'FAIL'} ({sum(1 for r in results.values() if r['passed'])}/{len(results)} checks passed)")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description="Compare extracted Abaqus results to reference values.")
    parser.add_argument("--extracted", required=True, help="Path to extracted JSON from extract_odb.py")
    parser.add_argument("--reference", required=True, help="Path to reference JSON with checks")
    parser.add_argument("--output", default=None, help="Optional JSON output for the report")
    parser.add_argument("--rtol", type=float, default=1e-4, help="Default relative tolerance")
    args = parser.parse_args()

    with open(args.extracted) as f:
        extracted = json.load(f)
    with open(args.reference) as f:
        reference = json.load(f)

    all_pass, results = run_comparison(extracted, reference, default_rtol=args.rtol)
    print_report(all_pass, results)

    if args.output:
        report = {
            "passed": all_pass,
            "extracted_file": args.extracted,
            "reference_file": args.reference,
            "results": results
        }
        with open(args.output, "w") as f:
            json.dump(report, f, indent=2)
        print(f"\nWrote report to {args.output}")

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
