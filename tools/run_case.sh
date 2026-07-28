#!/bin/bash
# Per-case Abaqus runner + ODB extractor for abaqus_ufl integration checks.
#
# Usage (from inside a case directory):
#   bash /path/to/run_case.sh <job_name> <user_for> [extract_config.json]
#
# Environment variables:
#   ABQ_CMD       - Abaqus command (default: abaqus)
#   ABQ_CPUS      - CPUs for Abaqus (default: 1)
#   ABQ_MP_MODE   - mp_mode for Abaqus (default: threads)
#   VALIDATION_SCRIPTS - Directory containing extract_odb.py (auto-detected if unset)
#
# Example:
#   cd examples/<name>/abaqus
#   bash ../../../tools/run_case.sh job ../<name>_umat.for extract_config.json

set -e

ABQ_CMD="${ABQ_CMD:-abaqus}"
ABQ_CPUS="${ABQ_CPUS:-1}"
ABQ_MP_MODE="${ABQ_MP_MODE:-threads}"

# Auto-detect scripts directory if not set
if [ -z "$VALIDATION_SCRIPTS" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    VALIDATION_SCRIPTS="$SCRIPT_DIR"
fi

JOB_NAME="${1:-job}"
USER_FOR="${2:-}"
EXTRACT_CONFIG="${3:-extract_config.json}"

if [ -z "$USER_FOR" ]; then
    echo "Usage: run_case.sh <job_name> <user_for> [extract_config.json]"
    exit 1
fi

echo "=========================================="
echo "Running Abaqus integration check"
echo "  job_name = $JOB_NAME"
echo "  user_for = $USER_FOR"
echo "  extract  = $EXTRACT_CONFIG"
echo "=========================================="

# Run Abaqus
if [ "$ABQ_CPUS" -gt 1 ]; then
    $ABQ_CMD job="$JOB_NAME" user="$USER_FOR" interactive \
        ask_delete=OFF mp_mode="$ABQ_MP_MODE" cpus="$ABQ_CPUS" \
        > "${JOB_NAME}.log" 2>&1
else
    $ABQ_CMD job="$JOB_NAME" user="$USER_FOR" interactive \
        ask_delete=OFF \
        > "${JOB_NAME}.log" 2>&1
fi

# Check convergence
if grep -Eq "THE ANALYSIS HAS (BEEN )?COMPLETED|THE ANALYSIS HAS COMPLETED SUCCESSFULLY" "${JOB_NAME}.sta" 2>/dev/null; then
    echo "[PASS] Abaqus analysis completed successfully"
else
    echo "[FAIL] Abaqus analysis did not complete. Check ${JOB_NAME}.sta and ${JOB_NAME}.msg"
    exit 1
fi

# Extract ODB if config exists
if [ -f "$EXTRACT_CONFIG" ]; then
    echo "[INFO] Extracting ODB with config: $EXTRACT_CONFIG"
    $ABQ_CMD python "$VALIDATION_SCRIPTS/extract_odb.py" \
        --config "$EXTRACT_CONFIG" \
        --output "${JOB_NAME}_extracted.json"
    echo "[PASS] ODB extracted to ${JOB_NAME}_extracted.json"
else
    echo "[WARN] No extract config found ($EXTRACT_CONFIG); skipping ODB extraction"
fi

echo "=========================================="
echo "Case complete: $JOB_NAME"
echo "=========================================="
