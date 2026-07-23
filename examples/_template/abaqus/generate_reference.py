"""Produce reference.json: the values the Abaqus run must match.

The reference is an INDEPENDENT Python material-point integration of the same
model along the same deformation path as job.inp — not a copy of a prior
Abaqus result. That independence is what makes the comparison a real check.

Run once to (re)generate the reference:

    python generate_reference.py

Then, after the Abaqus job:

    python ../../../tools/compare_results.py job_extracted.json reference.json
"""
from pathlib import Path
import sys

# repo root = .../examples/_template/abaqus -> up 3
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from tools.reference_utils import (  # shared helpers, see tools/reference_utils.py
    load_symbol,
    small_strain_path,
    stress_checks,
    strain_checks,
    write_reference,
)


def main():
    model_cls = load_symbol("examples/_template/build.py", "TemplateMaterial")

    # Drive the same path as job.inp (here: 2% x-compression over 10 increments).
    ref = small_strain_path(model_cls, final_stretch_x=0.98, ninc=10)

    checks = {}
    checks.update(stress_checks(ref))
    checks.update(strain_checks(ref))

    write_reference(
        Path(__file__).parent,
        "TEMPLATE: <one-line description of model + loading path>.",
        checks,
    )


if __name__ == "__main__":
    main()
