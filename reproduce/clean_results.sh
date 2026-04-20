#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
VENV_PATH="${REPO_ROOT}/.venv"
PYTHON_BIN="${VENV_PATH}/bin/python"
CLEAN_SCRIPT="${REPO_ROOT}/reproduce/clean_results.py"

if [[ ! -f "${VENV_PATH}/bin/activate" ]]; then
  echo "Error: virtual environment not found at ${VENV_PATH}" >&2
  echo "Please create it first, e.g. python -m venv .venv" >&2
  exit 1
fi

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Error: Python executable not found at ${PYTHON_BIN}" >&2
  exit 1
fi

if [[ ! -f "${CLEAN_SCRIPT}" ]]; then
  echo "Error: cleanup script not found at ${CLEAN_SCRIPT}" >&2
  exit 1
fi

source "${VENV_PATH}/bin/activate"

cd "${REPO_ROOT}"
PYTHONPATH=src:. "${PYTHON_BIN}" "${CLEAN_SCRIPT}" "$@"
