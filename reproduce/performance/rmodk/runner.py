from __future__ import annotations

import csv
import multiprocessing as mp
import threading
import time
from datetime import datetime
from multiprocessing import Process, Queue
from pathlib import Path
from zoneinfo import ZoneInfo

from tqdm import tqdm

from reproduce.performance.runner import monitor_memory as shared_monitor_memory
from reproduce.utils.logging_control import (
    apply_python_logging_policy,
    env_disables_python_logging,
    silence_process_stdio,
)
from reproduce.utils.model_resolver import resolve_model_file

from .benchmarks import (
    CSV_FIELDNAMES,
    MEMORY_POLL_INTERVAL,
    TIMEZONE,
    RuntimeConfig,
)
from .fo2cnf import CNFContext
from .render import plot_metric


def _normalize_status(result: object) -> str:
    if result in ("timeout", "hard_timeout", "error"):
        return str(result)

    if isinstance(result, str) and result.startswith("python_error:"):
        return "error"

    numeric_types: tuple[type, ...] = (int, float)
    try:
        from symengine import Number as SymNumber

        numeric_types = numeric_types + (SymNumber,)
    except Exception:
        pass

    try:
        from flint.types.fmpq import fmpq

        numeric_types = numeric_types + (fmpq,)
    except Exception:
        pass

    return "completed" if isinstance(result, numeric_types) else "error"


def _is_failure_status(status: str) -> bool:
    normalized = status.lower()
    return "timeout" in normalized or "error" in normalized


def print_result(result_data: dict[str, object], timeout_seconds: int) -> None:
    result_str = f"""
{'-' * 20} START {'-' * 20}
Timestamp: {datetime.now(ZoneInfo(TIMEZONE)).strftime('%Y-%m-%d %H:%M:%S')}
Timeout setting: {timeout_seconds} seconds
Formula: {result_data['formula']}
Domain size: {result_data['domain_size']}
Algorithm: {result_data['algorithm']}
Result: {result_data['result']}
Time taken: {result_data['time_sec']} seconds
Memory usage: {result_data['memory_bytes']} bytes ({result_data['memory_kb']} KB, {result_data['memory_mb']} MB)
Status: {result_data['status']}
{'-' * 20} END {'-' * 20}
"""
    print(result_str)


def _resolve_model_path(models_path: Path, model_filename: str) -> Path:
    return resolve_model_file(models_path, model_filename)


def _run_fo2cnf_counter(
    config: RuntimeConfig,
    model_name: str,
    model_filename: str,
    domain_size: int,
    counter: str,
) -> int:
    model_path = _resolve_model_path(config.models_path, model_filename)
    model_cnf_dir = config.cnf_results_root / model_name

    context = CNFContext(model_path, domain_size=domain_size, output_dir=model_cnf_dir)
    if not context.exists_cnf_file():
        context.convert()
        context.dump()

    normalized_counter = counter.lower()
    if normalized_counter == "ganak":
        return CNFContext.model_count_ganak(context.cnf_path)
    if normalized_counter == "approxmc":
        return CNFContext.model_count_approxmc(
            context.cnf_path,
            epsilon=config.epsilon,
            delta=config.delta,
        )

    raise ValueError(f"Unknown counter: {counter}")


def python_worker(
    q: Queue,
    config: RuntimeConfig,
    model_name: str,
    model_filename: str,
    domain_size: int,
    algorithm_name: str,
) -> None:
    """Execute one rmodk benchmark point in an isolated subprocess."""
    try:
        quiet = env_disables_python_logging()
        silence_process_stdio(quiet)
        apply_python_logging_policy(quiet)

        normalized_algorithm = algorithm_name.lower()

        if normalized_algorithm in {"ganak", "approxmc"}:
            result = _run_fo2cnf_counter(
                config=config,
                model_name=model_name,
                model_filename=model_filename,
                domain_size=domain_size,
                counter=normalized_algorithm,
            )
            q.put(result)
            return

        from wfomc import Algo as WfomcAlgo
        from wfomc import Const, parse_input, wfomc

        apply_python_logging_policy(quiet)

        model_path = _resolve_model_path(config.models_path, model_filename)
        problem = parse_input(str(model_path))
        problem.domain = {Const(f"d{i}") for i in range(domain_size)}

        wfomc_algo = WfomcAlgo(normalized_algorithm)
        result = wfomc(problem, wfomc_algo)
        q.put(result)
    except Exception as exc:
        q.put(f"python_error: {exc}")


def run_python_single(
    config: RuntimeConfig,
    model_name: str,
    model_filename: str,
    domain_size: int,
    algorithm_name: str,
) -> dict[str, object]:
    """Run one measurement point with timeout + memory tracking."""
    start_time = time.time()
    q: Queue = Queue()

    process = Process(
        target=python_worker,
        args=(
            q,
            config,
            model_name,
            model_filename,
            domain_size,
            algorithm_name,
        ),
    )
    process.start()

    peak_mem = mp.Value("Q", 0)
    stop_event = threading.Event()
    monitor_thread = threading.Thread(
        target=shared_monitor_memory,
        args=(
            process.pid,
            MEMORY_POLL_INTERVAL,
            peak_mem,
            stop_event,
        ),
        daemon=True,
    )
    monitor_thread.start()

    process.join(config.timeout_seconds)
    if process.is_alive():
        process.terminate()
        process.join()
        result: object = "timeout"
    else:
        result = q.get() if not q.empty() else "error"

    stop_event.set()
    monitor_thread.join(timeout=1)

    end_time = time.time()
    memory_used = peak_mem.value
    status = _normalize_status(result)

    return {
        "timestamp": datetime.now(ZoneInfo(TIMEZONE)).strftime("%Y%m%d_%H%M%S"),
        "formula": model_name,
        "domain_size": domain_size,
        "algorithm": algorithm_name,
        "result": result,
        "time_sec": round(end_time - start_time, 4),
        "memory_bytes": memory_used,
        "memory_kb": memory_used / 1024,
        "memory_mb": memory_used / (1024 * 1024),
        "status": status,
    }


def _count_total_iterations(config: RuntimeConfig) -> int:
    total = 0
    for group in config.groups:
        total += len(group["models"]) * len(group["domain_sizes"]) * len(group["algorithms"])
    return total


def run_experiment(
    config: RuntimeConfig,
    results_path: Path,
    flush_every: int = 20,
    disable_python_logging: bool = False,
) -> None:
    """Run rmodk performance experiments and write CSV + Figure 9 plots."""
    apply_python_logging_policy(disable_python_logging, update_env=True)

    total_iterations = _count_total_iterations(config)

    with tqdm(total=total_iterations, desc="Rmodk Progress") as progress_bar:
        for group in config.groups:
            group_name = group["name"]
            domain_sizes = list(group["domain_sizes"])
            algorithms = [str(algo).lower() for algo in group["algorithms"]]
            models = group["models"]

            group_results_path = results_path / group_name
            group_results_path.mkdir(parents=True, exist_ok=True)

            for model_name, model_filename in models.items():
                model_results_path = group_results_path / model_name
                model_results_path.mkdir(parents=True, exist_ok=True)

                csv_path = model_results_path / f"{model_name}.csv"

                with csv_path.open("w", newline="", encoding="utf-8") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDNAMES)
                    writer.writeheader()

                    skip_remaining_by_algo = {algo: False for algo in algorithms}

                    for algorithm_name in algorithms:
                        for domain_size in domain_sizes:
                            if skip_remaining_by_algo[algorithm_name]:
                                progress_bar.update(1)
                                continue

                            result_data = run_python_single(
                                config=config,
                                model_name=model_name,
                                model_filename=model_filename,
                                domain_size=domain_size,
                                algorithm_name=algorithm_name,
                            )

                            writer.writerow(result_data)
                            progress_bar.update(1)

                            if progress_bar.n % flush_every == 0:
                                csvfile.flush()

                            print_result(result_data, config.timeout_seconds)

                            if _is_failure_status(str(result_data["status"])):
                                skip_remaining_by_algo[algorithm_name] = True
                                print(
                                    f"Algorithm {algorithm_name} failed at domain size "
                                    f"{domain_size}; skipping larger sizes for this algorithm."
                                )

                plot_metric(str(csv_path), config=config, metric="time_sec")
                print(
                    "\nThe rmodk performance experiment is completed! "
                    f"The results are saved in: {csv_path}"
                )
