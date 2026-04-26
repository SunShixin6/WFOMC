from __future__ import annotations

import csv
import multiprocessing as mp
import threading
import time
from datetime import datetime
from multiprocessing import Process, Queue
from pathlib import Path
from zoneinfo import ZoneInfo

import psutil
from tqdm import tqdm

from .benchmarks import (
    CSV_FIELDNAMES,
    MEMORY_POLL_INTERVAL,
    TIMEZONE,
    RuntimeConfig,
)
from .render import plot_metric
from .transforms import apply_benchmark_transforms
from reproduce.utils.logging_control import (
    apply_python_logging_policy,
    env_disables_python_logging,
    silence_process_stdio,
)
from reproduce.utils.model_resolver import resolve_model_file


def monitor_memory(
    pid: int, interval: float, peak_mem: mp.Value, stop_evt: threading.Event
) -> None:
    """
    Poll the memory usage of the specified process in a separate thread and record the peak value.
    """
    try:
        proc = psutil.Process(pid)
        while proc.is_running() and not stop_evt.is_set():
            try:
                rss = proc.memory_full_info().uss
            except psutil.AccessDenied:
                rss = proc.memory_info().rss

            for child in proc.children(recursive=True):
                try:
                    rss += child.memory_full_info().uss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue

            with peak_mem.get_lock():
                if rss > peak_mem.value:
                    peak_mem.value = rss

            stop_evt.wait(interval)
    except psutil.NoSuchProcess:
        pass


def print_result(result_data: dict, timeout_seconds: int) -> None:
    result_str = f"""
    {'-' * 20} START {'-' * 20}
    Timestamp of code execution: {datetime.now(ZoneInfo(TIMEZONE)).strftime('%Y-%m-%d %H:%M:%S')}
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


def normalize_algo_name(algo: object) -> str:
    """Normalize algorithm identifiers across process boundaries.

    We may receive enum members from different module namespaces
    (e.g., `src.wfomc.Algo` vs `wfomc.Algo`), but workers always
    resolve using `wfomc.Algo`.
    """
    value = getattr(algo, "value", None)
    if isinstance(value, str):
        return value
    return str(algo).strip()


def python_worker(
    q: Queue,
    model_name: str,
    model_csv_name: str,
    domain_size: int,
    algo_name: str,
    models_path: Path,
) -> None:
    """Call wfomc in the independent subprocess"""
    try:
        quiet = env_disables_python_logging()
        silence_process_stdio(quiet)
        apply_python_logging_policy(quiet)

        from wfomc import Algo as WfomcAlgo
        from wfomc import Const, parse_input, wfomc

        # wfomc may initialize its own log handlers during import.
        # Re-apply policy after import to keep quiet mode effective.
        apply_python_logging_policy(quiet)

        file_path = resolve_model_file(models_path, model_csv_name)
        problem = parse_input(str(file_path))
        problem.domain = {Const(f"d{i}") for i in range(domain_size)}

        apply_benchmark_transforms(problem, model_name, domain_size)

        wfomc_algo = WfomcAlgo(algo_name)

        res = wfomc(problem, wfomc_algo)
        q.put(res)
    except Exception as e:
        q.put(f"python_error: {e}")


def run_python_single(
    config: RuntimeConfig,
    model_name: str,
    model_csv_name: str,
    domain_size: int,
    algo,
    r=None,
    k=None,
) -> dict:
    """Run a single experiment"""
    start_time = time.time()
    q: Queue = Queue()
    algo_name = normalize_algo_name(algo)

    p = Process(
        target=python_worker,
        args=(
            q,
            model_name,
            model_csv_name,
            domain_size,
            algo_name,
            config.models_path,
        ),
    )
    p.start()

    peak_mem = mp.Value("Q", 0)
    stop_evt = threading.Event()
    monitor_thread = threading.Thread(
        target=monitor_memory,
        args=(
            p.pid,
            MEMORY_POLL_INTERVAL,
            peak_mem,
            stop_evt,
        ),
        daemon=True,
    )
    monitor_thread.start()

    p.join(config.timeout_seconds)
    if p.is_alive():
        p.terminate()
        p.join()
        result = "timeout"
    else:
        result = q.get() if not q.empty() else "error"

    monitor_thread.join()
    end_time = time.time()

    stop_evt.set()
    monitor_thread.join(timeout=1)
    memory_used = peak_mem.value

    status = "completed"
    if result in ("timeout", "hard_timeout", "error"):
        status = str(result)
    elif isinstance(result, str) and result.startswith("python_error:"):
        status = "error"

    return {
        "timestamp": datetime.now(ZoneInfo(TIMEZONE)).strftime("%Y%m%d_%H%M%S"),
        "formula": model_name,
        "domain_size": domain_size,
        "algorithm": algo_name,
        "result": result,
        "time_sec": round(end_time - start_time, 4),
        "memory_bytes": memory_used,
        "memory_kb": memory_used / 1024,
        "memory_mb": memory_used / (1024 * 1024),
        "status": status,
    }


def run_experiment(
    config: RuntimeConfig,
    results_path: Path,
    flush_every: int = 20,
    disable_python_logging: bool = False,
) -> None:
    """Run the entire experimental batch and generate CSV files and graphs"""
    apply_python_logging_policy(disable_python_logging, update_env=True)

    total_iterations = 0
    for group in config.groups:
        domain_size = group["domain_size"]
        algorithms = group["algorithms"]
        group_iterations = len(domain_size) * len(algorithms) * len(group["models"])
        total_iterations += group_iterations

    with tqdm(total=total_iterations, desc="Overall Progress") as pbar:
        for group in config.groups:
            group_name = group["name"]
            domain_size = group["domain_size"]
            algorithms = group["algorithms"]
            models = group["models"]

            group_results_path = results_path / group_name
            group_results_path.mkdir(parents=True, exist_ok=True)
            

            for model_name, model_csv_name in models.items():
                single_graph_results_path = group_results_path / model_name
                single_graph_results_path.mkdir(parents=True, exist_ok=True)

                csv_name = f"{model_name}.csv"
                csv_path = single_graph_results_path / csv_name

                with open(csv_path, "w", newline="") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDNAMES)
                    writer.writeheader()

                    for algo_name in [str(a).lower() for a in algorithms]:
                        skip_domain_size = False
                        for n in domain_size:
                            if skip_domain_size:
                                pbar.update(1)
                                continue

                            single_result = run_python_single(
                                config=config,
                                model_name=model_name,
                                model_csv_name=model_csv_name,
                                domain_size=n,
                                algo=algo_name,
                            )

                            if (
                                single_result["status"] == "timeout"
                                or single_result["status"] == "hard_timeout"
                                or single_result["status"] == "error"
                            ):
                                skip_domain_size = True
                                print(
                                    f"Algorithm timed out at domain size {n}; skipping remaining domain sizes"
                                )

                            result_data = single_result
                            if result_data:
                                writer.writerow(result_data)
                                pbar.update(1)
                                if pbar.n % flush_every == 0:
                                    csvfile.flush()
                                print_result(result_data, config.timeout_seconds)

                for metric in ("time_sec", "memory_mb"):
                    plot_metric(str(csv_path), config=config, metric=metric)
                print(f"\nThe experiment is completed! The results are saved in...: {csv_path}")
