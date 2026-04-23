#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
VENV_PATH="${REPO_ROOT}/.venv"
PYTHON_BIN="${VENV_PATH}/bin/python"
TARGET_SCRIPT="${SCRIPT_DIR}/ganak_odd_degree.py"
LOG_HELPER="${REPO_ROOT}/reproduce/utils/run_logging.sh"
USE_UV=false
PYTHON_CMD=()

print_usage() {
  cat <<'EOF'
Usage: run_ganak_odd_degree.sh [args...]

Runs reproduce/OEISsequence/ganak_odd_degree.py from the repository root.

All arguments are passed through to the Python script.
Examples:
  ./run_ganak_odd_degree.sh
  ./run_ganak_odd_degree.sh --ganak-path /path/to/ganak --approxmc-path /path/to/approxmc

Options:
  --disable-python-logging
  --enable-python-logging
  -h, --help

Optional environment variables:
  GANAK_BIN
  APPROXMC_BIN
  REPRO_DISABLE_PYTHON_LOGGING
EOF
}

PYTHON_LOGGING_DISABLED="${REPRO_DISABLE_PYTHON_LOGGING:-0}"
FORWARDED_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      print_usage
      exit 0
      ;;
    --disable-python-logging)
      PYTHON_LOGGING_DISABLED="1"
      shift
      ;;
    --enable-python-logging)
      PYTHON_LOGGING_DISABLED="0"
      shift
      ;;
    *)
      FORWARDED_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ -x "${PYTHON_BIN}" ]]; then
  PYTHON_CMD=("${PYTHON_BIN}")
elif command -v uv >/dev/null 2>&1; then
  USE_UV=true
  PYTHON_CMD=(uv run python)
  echo "Info: ${VENV_PATH} not found; falling back to 'uv run python'."
else
  echo "Error: Python runtime not found." >&2
  echo "Expected ${PYTHON_BIN}, or install uv to enable fallback execution." >&2
  exit 1
fi

if [[ ! -f "${TARGET_SCRIPT}" ]]; then
  echo "Error: target script not found at ${TARGET_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${LOG_HELPER}" ]]; then
  echo "Error: log helper script not found at ${LOG_HELPER}" >&2
  exit 1
fi

# shellcheck source=/dev/null
source "${LOG_HELPER}"

RUN_START_TIME="$(repro_now_iso)"
RUN_LOG_DIR="$(repro_prepare_log_dir "${REPO_ROOT}" "single")"
RUN_LOG_FILE="${RUN_LOG_DIR}/oeissequence_table2.log"
RUN_META_FILE="${RUN_LOG_DIR}/oeissequence_table2.meta.json"

trap 'exit_code=$?; repro_write_meta "${RUN_META_FILE}" "${REPO_ROOT}" "reproduce/OEISsequence/run_ganak_odd_degree.sh" "${RUN_START_TIME}" "${RUN_LOG_DIR}" "${RUN_LOG_FILE}" "${exit_code}"' EXIT

if [[ "${USE_UV}" == "false" ]]; then
  source "${VENV_PATH}/bin/activate"
fi

cd "${REPO_ROOT}"
echo "Run logs directory: ${RUN_LOG_DIR}"
echo "Step log: ${RUN_LOG_FILE}"
echo "REPRO_DISABLE_PYTHON_LOGGING=${PYTHON_LOGGING_DISABLED}"
if ! repro_run_and_tee "reproduce/OEISsequence/ganak_odd_degree.py" "${RUN_LOG_FILE}" env PYTHONPATH=src:. REPRO_DISABLE_PYTHON_LOGGING="${PYTHON_LOGGING_DISABLED}" "${PYTHON_CMD[@]}" "${TARGET_SCRIPT}" "${FORWARDED_ARGS[@]}"; then
  exit 1
fi
