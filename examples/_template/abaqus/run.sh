#!/bin/bash
# Run this optional case in Abaqus, extract the ODB, and compare the frozen
# closed-form reference. Abaqus model construction remains example-owned.
# Requires an Abaqus/Standard install on PATH (set ABQ_CMD to override).
#
#   bash run.sh
#
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="$(cd "$HERE/../../../tools" && pwd)"

# The generated UMAT/UEL lives one level up (examples/<name>/template_umat.for).
USER_FOR="$HERE/../template_umat.for"

cd "$HERE"
bash "$TOOLS/run_case.sh" job "$USER_FOR" "$HERE/extract_config.json"
python "$TOOLS/compare_results.py" \
  --extracted "$HERE/job_extracted.json" \
  --reference "$HERE/reference.json" \
  --output "$HERE/job_report.json"
