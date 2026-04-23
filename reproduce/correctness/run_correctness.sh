#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
VENV_PATH="${REPO_ROOT}/.venv"
PYTHON_BIN="${VENV_PATH}/bin/python"
LOG_HELPER="${REPO_ROOT}/reproduce/utils/run_logging.sh"
USE_UV=false
PYTHON_CMD=()

print_usage() {
  cat <<'EOF'
Usage: run_correctness.sh [args...]

Run reproduce.correctness.main only.

For odd-degree correctness, use:
  bash reproduce/correctness/odd_degree/run_odd_degree.sh [args...]

All remaining arguments are passed through to reproduce.correctness.main.
EOF
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    --odd-degree)
      echo "Error: --odd-degree is not supported by run_correctness.sh." >&2
      echo "Use: bash reproduce/correctness/odd_degree/run_odd_degree.sh" >&2
      exit 1
      ;;
    --main)
      # Backward-compatible no-op.
      shift
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
  esac
fi

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
RUN_LOG_FILE="${RUN_LOG_DIR}/correctness_main.log"
RUN_META_FILE="${RUN_LOG_DIR}/correctness_main.meta.json"

trap 'exit_code=$?; repro_write_meta "${RUN_META_FILE}" "${REPO_ROOT}" "reproduce/correctness/run_correctness.sh" "${RUN_START_TIME}" "${RUN_LOG_DIR}" "${RUN_LOG_FILE}" "${exit_code}"' EXIT

if [[ "${USE_UV}" == "false" ]]; then
  source "${VENV_PATH}/bin/activate"
fi

cd "${REPO_ROOT}"

echo "[1/1] Running reproduce.correctness.main"
echo "Run logs directory: ${RUN_LOG_DIR}"
echo "Step log: ${RUN_LOG_FILE}"
if ! repro_run_and_tee "reproduce.correctness.main" "${RUN_LOG_FILE}" env PYTHONPATH=src:. "${PYTHON_CMD[@]}" -m reproduce.correctness.main "$@"; then
  exit 1
fi

echo "Correctness main experiment finished."
