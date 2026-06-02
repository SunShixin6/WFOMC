#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
VENV_PATH="${REPO_ROOT}/.venv"
PYTHON_BIN="${VENV_PATH}/bin/python"
OEIS_RUNNER="${REPO_ROOT}/reproduce/OEISsequence/run_ganak_odd_degree.sh"
LOG_HELPER="${REPO_ROOT}/reproduce/utils/run_logging.sh"

USE_UV=false
PYTHON_CMD=()

print_usage() {
	cat <<'EOF'
Usage: bash reproduce/run_all.sh [--smoke] [--dry-run] [args...]

Run all reproduction workflows from repository root.

Options:
	--smoke        Enable smoke profile (small benchmark ranges).
	--dry-run      Print planned commands and exit without executing.
	--verbose
	                       Enable Python logging in reproduce runners.
	--disable-python-logging
	                       Keep Python logging disabled (default).
	--enable-python-logging
	                       Alias for --verbose.
	-h, --help     Show this help message.

Argument forwarding:
	Forwarded to:
		- reproduce.correctness.main
		- reproduce.correctness.odd_degree.main
		- reproduce.performance.odd_degree.main
		- reproduce.performance.rmodk.main
		- (compatible generic args only)

	OEIS-only arguments:
		- --ganak-path
		- --approxmc-path
		- --max-n
	Not forwarded to:
		- reproduce.performance.main

Examples:
	bash reproduce/run_all.sh
	bash reproduce/run_all.sh --smoke
	bash reproduce/run_all.sh --dry-run
	bash reproduce/run_all.sh --smoke --ganak-path /path/to/ganak --approxmc-path /path/to/approxmc --max-n 4
EOF
}

DRY_RUN=false
SMOKE_MODE=false
DISABLE_PYTHON_LOGGING=true
MODULE_ARGS=()
OEIS_ARGS=()

if [[ -n "${REPRO_DISABLE_PYTHON_LOGGING:-}" ]]; then
        if [[ "${REPRO_DISABLE_PYTHON_LOGGING}" == "0" ]]; then
                DISABLE_PYTHON_LOGGING=false
        else
                DISABLE_PYTHON_LOGGING=true
        fi
fi


while [[ $# -gt 0 ]]; do
	case "$1" in
		-h|--help)
			print_usage
			exit 0
			;;
                --dry-run)
                        DRY_RUN=true
                        shift
                        ;;
                --smoke)
                        SMOKE_MODE=true
                        shift
                        ;;
                --disable-python-logging)
                        DISABLE_PYTHON_LOGGING=true
                        shift
                        ;;
                --verbose|--enable-python-logging)
                        DISABLE_PYTHON_LOGGING=false
                        shift
                        ;;
		--ganak-path|--approxmc-path|--max-n)
			if [[ $# -lt 2 ]]; then
				echo "Error: option '$1' requires a value" >&2
				exit 1
			fi
			OEIS_ARGS+=("$1" "$2")
			shift 2
			;;
		--ganak-path=*|--approxmc-path=*|--max-n=*)
			OEIS_ARGS+=("$1")
			shift
			;;
		*)
			MODULE_ARGS+=("$1")
			shift
			;;
	esac
done

if [[ "${SMOKE_MODE}" == "true" ]]; then
	export REPRO_SMOKE=1
	echo "Smoke mode enabled (REPRO_SMOKE=1)"
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

if [[ ! -f "${OEIS_RUNNER}" ]]; then
	echo "Error: OEIS runner script not found at ${OEIS_RUNNER}" >&2
	exit 1
fi

if [[ ! -f "${LOG_HELPER}" ]]; then
	echo "Error: log helper script not found at ${LOG_HELPER}" >&2
	exit 1
fi

# shellcheck source=/dev/null
source "${LOG_HELPER}"

RUN_START_TIME="$(repro_now_iso)"
RUN_LOG_DIR="$(repro_prepare_log_dir "${REPO_ROOT}" "run")"
RUN_LOG_FILE="${RUN_LOG_DIR}/run_all.log"
RUN_META_FILE="${RUN_LOG_DIR}/run_meta.json"

trap 'exit_code=$?; repro_write_meta "${RUN_META_FILE}" "${REPO_ROOT}" "reproduce/run_all.sh" "${RUN_START_TIME}" "${RUN_LOG_DIR}" "${RUN_LOG_FILE}" "${exit_code}"' EXIT

if [[ "${USE_UV}" == "false" ]]; then
	source "${VENV_PATH}/bin/activate"
fi

cd "${REPO_ROOT}"

mkdir -p "${RUN_LOG_DIR}"
exec > >(tee -a "${RUN_LOG_FILE}") 2>&1

echo "Run logs directory: ${RUN_LOG_DIR}"
echo "Main run log: ${RUN_LOG_FILE}"

if [[ "${DISABLE_PYTHON_LOGGING}" == "true" ]]; then
	export REPRO_DISABLE_PYTHON_LOGGING=1
else
	export REPRO_DISABLE_PYTHON_LOGGING=0
fi
echo "REPRO_DISABLE_PYTHON_LOGGING=${REPRO_DISABLE_PYTHON_LOGGING}"

export REPRO_LOG_DIR="${RUN_LOG_DIR}"

MODULES=(
	"reproduce.correctness.main"
	"reproduce.correctness.odd_degree.main"
	"reproduce.performance.main"
	"reproduce.performance.odd_degree.main"
	"reproduce.performance.rmodk.main"
)
TOTAL_STEPS="$(( ${#MODULES[@]} + 1 ))"

for i in "${!MODULES[@]}"; do
	STEP="$((i + 1))"
	MODULE="${MODULES[$i]}"
	MODULE_LOG_NAME="${MODULE//./_}"
	STEP_LOG_FILE="${RUN_LOG_DIR}/$(printf '%02d' "${STEP}")_${MODULE_LOG_NAME}.log"
	echo "[${STEP}/${TOTAL_STEPS}] Running ${MODULE}"
	echo "[${STEP}/${TOTAL_STEPS}] Step log: ${STEP_LOG_FILE}"

	if [[ "${DRY_RUN}" == "true" ]]; then
		if [[ "${MODULE}" == "reproduce.performance.main" ]]; then
			echo "[dry-run] PYTHONPATH=src:. ${PYTHON_CMD[*]} -m ${MODULE}" | tee -a "${STEP_LOG_FILE}"
		else
			echo "[dry-run] PYTHONPATH=src:. ${PYTHON_CMD[*]} -m ${MODULE} ${MODULE_ARGS[*]}" | tee -a "${STEP_LOG_FILE}"
		fi
		continue
	fi

	if [[ "${MODULE}" == "reproduce.performance.main" ]]; then
		if ! repro_run_and_tee "${MODULE}" "${STEP_LOG_FILE}" env PYTHONPATH=src:. "${PYTHON_CMD[@]}" -m "${MODULE}"; then
			exit 1
		fi
	else
		if ! repro_run_and_tee "${MODULE}" "${STEP_LOG_FILE}" env PYTHONPATH=src:. "${PYTHON_CMD[@]}" -m "${MODULE}" "${MODULE_ARGS[@]}"; then
			exit 1
		fi
	fi
done

echo "[${TOTAL_STEPS}/${TOTAL_STEPS}] Running reproduce/OEISsequence/run_ganak_odd_degree.sh"
OEIS_STEP_LOG="${RUN_LOG_DIR}/$(printf '%02d' "${TOTAL_STEPS}")_oeissequence_ganak_odd_degree.log"
echo "[${TOTAL_STEPS}/${TOTAL_STEPS}] Step log: ${OEIS_STEP_LOG}"

if [[ "${DRY_RUN}" == "true" ]]; then
	echo "[dry-run] bash ${OEIS_RUNNER} ${OEIS_ARGS[*]}" | tee -a "${OEIS_STEP_LOG}"
	echo "Dry run completed. No commands were executed."
	exit 0
fi

if ! repro_run_and_tee "reproduce/OEISsequence/run_ganak_odd_degree.sh" "${OEIS_STEP_LOG}" bash "${OEIS_RUNNER}" "${OEIS_ARGS[@]}"; then
	exit 1
fi

echo "All reproduce main experiments finished."
