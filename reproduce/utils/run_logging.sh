# Shared helpers for reproduce shell scripts.

repro_now_stamp() {
    date +"%Y%m%d_%H%M%S"
}


repro_now_iso() {
    date -u +"%Y-%m-%dT%H:%M:%SZ"
}


repro_prepare_log_dir() {
    local repo_root="$1"
    local mode="$2"

    if [[ -n "${REPRO_LOG_DIR:-}" ]]; then
        mkdir -p "${REPRO_LOG_DIR}"
        echo "${REPRO_LOG_DIR}"
        return 0
    fi

    local stamp
    stamp="$(repro_now_stamp)"
    local log_dir="${repo_root}/reproduce/logs/${mode}_${stamp}"
    mkdir -p "${log_dir}"
    echo "${log_dir}"
}


repro_write_meta() {
    local meta_path="$1"
    local repo_root="$2"
    local script_name="$3"
    local start_time="$4"
    local log_dir="$5"
    local log_file="$6"
    local exit_code="$7"

    local end_time
    end_time="$(repro_now_iso)"

    local status="failed"
    if [[ "${exit_code}" -eq 0 ]]; then
        status="success"
    fi

    local git_commit="unknown"
    if git -C "${repo_root}" rev-parse --short HEAD >/dev/null 2>&1; then
        git_commit="$(git -C "${repo_root}" rev-parse --short HEAD)"
    fi

    cat > "${meta_path}" <<EOF
{
  "script": "${script_name}",
  "start_time_utc": "${start_time}",
  "end_time_utc": "${end_time}",
  "status": "${status}",
  "exit_code": ${exit_code},
  "repo_root": "${repo_root}",
  "git_commit": "${git_commit}",
  "log_dir": "${log_dir}",
  "log_file": "${log_file}"
}
EOF
}


repro_run_and_tee() {
    local step_name="$1"
    local step_log="$2"
    shift 2

    echo "[$(repro_now_iso)] START ${step_name}" | tee -a "${step_log}"
    if "$@" 2>&1 | tee -a "${step_log}"; then
        echo "[$(repro_now_iso)] END ${step_name} status=0" | tee -a "${step_log}"
        return 0
    fi

    local status=$?
    echo "[$(repro_now_iso)] END ${step_name} status=${status}" | tee -a "${step_log}"
    return "${status}"
}
