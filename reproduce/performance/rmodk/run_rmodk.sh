#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
VENV_PATH="${REPO_ROOT}/.venv"
PYTHON_BIN="${VENV_PATH}/bin/python"
LOG_HELPER="${REPO_ROOT}/reproduce/utils/run_logging.sh"
USE_UV=false
PYTHON_CMD=()

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

if [[ ! -f "${LOG_HELPER}" ]]; then
  echo "Error: log helper script not found at ${LOG_HELPER}" >&2
  exit 1
fi

# shellcheck source=/dev/null
source "${LOG_HELPER}"

RUN_START_TIME="$(repro_now_iso)"
RUN_LOG_DIR="$(repro_prepare_log_dir "${REPO_ROOT}" "single")"
RUN_LOG_FILE="${RUN_LOG_DIR}/performance_rmodk.log"
RUN_META_FILE="${RUN_LOG_DIR}/performance_rmodk.meta.json"

trap 'exit_code=$?; repro_write_meta "${RUN_META_FILE}" "${REPO_ROOT}" "reproduce/performance/rmodk/run_rmodk.sh" "${RUN_START_TIME}" "${RUN_LOG_DIR}" "${RUN_LOG_FILE}" "${exit_code}"' EXIT

if [[ "${USE_UV}" == "false" ]]; then
  source "${VENV_PATH}/bin/activate"
fi

cd "${REPO_ROOT}"
echo "Run logs directory: ${RUN_LOG_DIR}"
echo "Step log: ${RUN_LOG_FILE}"
if ! repro_run_and_tee "reproduce.performance.rmodk.main" "${RUN_LOG_FILE}" env PYTHONPATH=src:. "${PYTHON_CMD[@]}" -m reproduce.performance.rmodk.main "$@"; then
  exit 1
fi
