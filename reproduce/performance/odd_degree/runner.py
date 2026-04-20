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

from .benchmarks import (
    CSV_FIELDNAMES,
    MEMORY_POLL_INTERVAL,
    TIMEZONE,
    RuntimeConfig,
)
from .ganak_odd_degree import fo2_count_odd_degree
from .render import plot_metric
from .transforms import build_problem_for_wfomc, resolve_model_path
from reproduce.utils.logging_control import (
    apply_python_logging_policy,
    env_disables_python_logging,
    silence_process_stdio,
)


def _normalize_status(result: object) -> str:
    if isinstance(result, str):
        return result

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
K value: {result_data['k_value']}
M value: {result_data['m_value']}
Algorithm: {result_data['algorithm']}
Result: {result_data['result']}
Time taken: {result_data['time_sec']} seconds
Memory usage: {result_data['memory_bytes']} bytes ({result_data['memory_kb']} KB, {result_data['memory_mb']} MB)
Status: {result_data['status']}
{'-' * 20} END {'-' * 20}
"""
    print(result_str)


def python_worker(
    q: Queue,
    config: RuntimeConfig,
    model_name: str,
    model_filename: str,
    domain_size: int,
    m_value: int,
    directed_k_value: int,
    algorithm_name: str,
) -> None:
    """Execute a single odd-degree benchmark point in an isolated subprocess."""
    try:
        quiet = env_disables_python_logging()
        silence_process_stdio(quiet)
        apply_python_logging_policy(quiet)

        normalized_algorithm = algorithm_name.lower()
        model_path = resolve_model_path(config.models_path, model_filename)

        if normalized_algorithm in {"ganak", "approxmc"}:
            undirected_k_value = directed_k_value // 2
            result = fo2_count_odd_degree(
                file_path=model_path,
                n=domain_size,
                m=m_value,
                k=undirected_k_value,
                counter=normalized_algorithm,
                epsilon=config.epsilon,
                delta=config.delta,
                output_dir=config.cnf_results_root / model_name,
            )
            q.put(result)
            return

        from wfomc import Algo, wfomc

        # wfomc may initialize its own log handlers during import.
        # Re-apply policy after import to keep quiet mode effective.
        apply_python_logging_policy(quiet)

        wfomc_algo = next((algo for algo in Algo if str(algo) == normalized_algorithm), None)
        if wfomc_algo is None:
            raise ValueError(f"Unsupported algorithm: {algorithm_name}")

        problem = build_problem_for_wfomc(
            model_path,
            domain_size=domain_size,
            m_value=m_value,
            directed_k_value=directed_k_value,
        )
        result = wfomc(problem, wfomc_algo)
        q.put(result)
    except Exception as exc:
        q.put(f"python_error: {exc}")


def run_python_single(
    config: RuntimeConfig,
    model_name: str,
    model_filename: str,
    domain_size: int,
    m_value: int,
    directed_k_value: int,
    algorithm_name: str,
) -> dict[str, object]:
    """Run one odd-degree benchmark point and collect runtime/memory stats."""
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
            m_value,
            directed_k_value,
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
        "k_value": directed_k_value,
        "m_value": m_value,
        "algorithm": algorithm_name,
        "result": result,
        "time_sec": round(end_time - start_time, 4),
        "memory_bytes": memory_used,
        "memory_kb": memory_used / 1024,
        "memory_mb": memory_used / (1024 * 1024),
        "status": status,
    }


def _count_total_iterations(config: RuntimeConfig) -> int:
    return len(config.domain_sizes) * len(config.m_values) * len(config.algorithms)


def run_experiment(
    config: RuntimeConfig,
    results_path: Path,
    flush_every: int = 20,
    disable_python_logging: bool = False,
) -> None:
    """Run fixed odd-degree benchmark experiments and write CSV + plots."""
    apply_python_logging_policy(disable_python_logging, update_env=True)

    total_iterations = _count_total_iterations(config)
    model_results_path = results_path / config.group_name / config.model_name
    model_results_path.mkdir(parents=True, exist_ok=True)

    with tqdm(total=total_iterations, desc="Odd-degree Progress") as progress_bar:
        for m_value in config.m_values:
            csv_path = model_results_path / f"{config.model_name}_m{m_value}.csv"

            with csv_path.open("w", newline="", encoding="utf-8") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDNAMES)
                writer.writeheader()

                for algorithm_name in config.algorithms:
                    normalized_algorithm = str(algorithm_name).lower()
                    skip_remaining_domain = False

                    for domain_size in config.domain_sizes:
                        if skip_remaining_domain:
                            progress_bar.update(1)
                            continue

                        directed_k_value = config.k_multiplier * domain_size
                        result_data = run_python_single(
                            config=config,
                            model_name=config.model_name,
                            model_filename=config.model_filename,
                            domain_size=domain_size,
                            m_value=m_value,
                            directed_k_value=directed_k_value,
                            algorithm_name=normalized_algorithm,
                        )

                        writer.writerow(result_data)
                        progress_bar.update(1)
                        if progress_bar.n % flush_every == 0:
                            csvfile.flush()

                        print_result(result_data, config.timeout_seconds)

                        if _is_failure_status(str(result_data["status"])):
                            skip_remaining_domain = True
                            print(
                                "Algorithm failed at domain size "
                                f"{domain_size}; skipping larger sizes for this setup."
                            )

            plot_metric(str(csv_path), config=config, metric="time_sec")
            print(f"\nOdd-degree m={m_value} experiment completed. Results saved to: {csv_path}")
