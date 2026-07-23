#!/bin/bash
# Run this case in Abaqus and extract the ODB, using the shared runner.
# Requires an Abaqus/Standard install on PATH (set ABQ_CMD to override).
#
#   bash run.sh
#
# Then compare against the reference:
#   python ../../../tools/compare_results.py job_extracted.json reference.json
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="$(cd "$HERE/../../../tools" && pwd)"

# The generated UMAT/UEL lives one level up (examples/<name>/template_umat.for).
USER_FOR="$HERE/../template_umat.for"

bash "$TOOLS/run_case.sh" job "$USER_FOR" "$HERE/extract_config.json"
