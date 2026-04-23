#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
VENV_PATH="${REPO_ROOT}/.venv"
PYTHON_BIN="${VENV_PATH}/bin/python"
CLEAN_MODULE="reproduce.utils.clean_results"
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

if [[ ! -f "${SCRIPT_DIR}/clean_results.py" ]]; then
  echo "Error: cleanup script not found at ${SCRIPT_DIR}/clean_results.py" >&2
  exit 1
fi

if [[ "${USE_UV}" == "false" ]]; then
  source "${VENV_PATH}/bin/activate"
fi

cd "${REPO_ROOT}"
PYTHONPATH=src:. "${PYTHON_CMD[@]}" -m "${CLEAN_MODULE}" "$@"