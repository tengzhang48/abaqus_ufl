"""Write the frozen closed-form reference for the optional Abaqus check.

This script intentionally does not import the Python material implementation.
It evaluates the constrained-strain isotropic-elastic solution independently.
Run it only when the documented benchmark itself changes; ordinary pipeline
checks consume the committed ``reference.json``.
"""

from pathlib import Path
import math
import sys


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.reference_utils import (
    strain_checks,
    stress_checks,
    write_reference,
)


G = 10.0
LAM = 20.0
FINAL_STRETCH_X = 0.98


def main():
    eps_abq = math.log(FINAL_STRETCH_X)
    eps_cp = -eps_abq
    ref = {
        "stress_abq": [
            -(2.0 * G + LAM) * eps_cp,
            -LAM * eps_cp,
            -LAM * eps_cp,
            0.0,
            0.0,
            0.0,
        ],
        "mises": 2.0 * G * eps_cp,
        "strain_abq": [eps_abq, 0.0, 0.0, 0.0, 0.0, 0.0],
    }
    checks = {}
    checks.update(stress_checks(ref, rtol=1.0e-5))
    checks.update(strain_checks(ref, rtol=1.0e-5))
    checks["step_name"] = {
        "path": "step",
        "expected": "LOAD",
    }
    checks["frame_value"] = {
        "path": "frame_value",
        "expected": 1.0,
        "rtol": 0.0,
        "atol": 1.0e-12,
    }
    for field_name in ("S", "LE"):
        key = field_name.lower()
        checks[key + "_record_count"] = {
            "path": "fields.{}.record_count".format(field_name),
            "expected": 1,
            "rtol": 0.0,
            "atol": 0.0,
        }
        checks[key + "_element_label"] = {
            "path": "fields.{}.element_label".format(field_name),
            "expected": 1,
            "rtol": 0.0,
            "atol": 0.0,
        }
        checks[key + "_integration_point"] = {
            "path": "fields.{}.integration_point".format(field_name),
            "expected": 1,
            "rtol": 0.0,
            "atol": 0.0,
        }
    write_reference(
        Path(__file__).parent,
        (
            "Closed-form small-strain elastic response for one constrained "
            "C3D8R under 2% x-compression."
        ),
        checks,
    )


if __name__ == "__main__":
    main()
