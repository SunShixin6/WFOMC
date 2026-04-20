#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
VENV_PATH="${REPO_ROOT}/.venv"
PYTHON_BIN="${VENV_PATH}/bin/python"
LOG_HELPER="${REPO_ROOT}/reproduce/utils/run_logging.sh"

if [[ ! -f "${VENV_PATH}/bin/activate" ]]; then
  echo "Error: virtual environment not found at ${VENV_PATH}" >&2
  echo "Please create it first, e.g. python -m venv .venv" >&2
  exit 1
fi

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Error: Python executable not found at ${PYTHON_BIN}" >&2
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
RUN_LOG_FILE="${RUN_LOG_DIR}/performance_odd_degree.log"
RUN_META_FILE="${RUN_LOG_DIR}/performance_odd_degree.meta.json"

trap 'exit_code=$?; repro_write_meta "${RUN_META_FILE}" "${REPO_ROOT}" "reproduce/performance/odd_degree/run_odd_degree.sh" "${RUN_START_TIME}" "${RUN_LOG_DIR}" "${RUN_LOG_FILE}" "${exit_code}"' EXIT

# Activate project virtual environment
source "${VENV_PATH}/bin/activate"

cd "${REPO_ROOT}"
echo "Run logs directory: ${RUN_LOG_DIR}"
echo "Step log: ${RUN_LOG_FILE}"
if ! repro_run_and_tee "reproduce.performance.odd_degree.main" "${RUN_LOG_FILE}" env PYTHONPATH=src:. "${PYTHON_BIN}" -m reproduce.performance.odd_degree.main "$@"; then
  exit 1
fi
