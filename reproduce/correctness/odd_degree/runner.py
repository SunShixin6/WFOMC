from __future__ import annotations

import csv
from datetime import datetime
from multiprocessing import Process, Queue
from pathlib import Path
from zoneinfo import ZoneInfo

from tqdm import tqdm

from reproduce.performance.odd_degree.ganak_odd_degree import fo2_count_odd_degree
from reproduce.performance.odd_degree.transforms import build_problem_for_wfomc, resolve_model_path
from reproduce.utils.logging_control import (
    apply_python_logging_policy,
    env_disables_python_logging,
    silence_process_stdio,
)

from .benchmarks import CSV_FIELDNAMES, TIMEZONE, RuntimeConfig, SweepSpec
from .render import plot_model_count


def _normalize_status(result: object) -> str:
    if result in ("timeout", "hard_timeout", "error"):
        return str(result)

    if isinstance(result, str):
        if result.startswith("python_error:"):
            return "error"
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


def _undirected_k_value(raw_k_value: int) -> int:
    return raw_k_value


def _to_directed_k_value(undirected_k_value: int) -> int:
    return 2 * undirected_k_value


def print_result(result_data: dict[str, object], timeout_seconds: int) -> None:
    result_str = f"""
{'-' * 20} START {'-' * 20}
Timestamp: {datetime.now(ZoneInfo(TIMEZONE)).strftime('%Y-%m-%d %H:%M:%S')}
Timeout setting: {timeout_seconds} seconds
Sweep: {result_data['sweep']}
Variable: {result_data['variable_name']} = {result_data['variable_value']}
Formula: {result_data['formula']}
Domain size: {result_data['domain_size']}
K value (undirected): {result_data['k_value']}
M value: {result_data['m_value']}
Algorithm: {result_data['algorithm']}
Model count: {result_data['model_count']}
Status: {result_data['status']}
{'-' * 20} END {'-' * 20}
"""
    print(result_str)


def python_worker(
    q: Queue,
    config: RuntimeConfig,
    domain_size: int,
    m_value: int,
    undirected_k_value: int,
    algorithm_name: str,
) -> None:
    """Execute one Figure13 point in an isolated subprocess."""
    try:
        quiet = env_disables_python_logging()
        silence_process_stdio(quiet)
        apply_python_logging_policy(quiet)

        normalized_algorithm = algorithm_name.lower()
        model_path = resolve_model_path(config.models_path, config.model_filename)
        directed_k_value = _to_directed_k_value(undirected_k_value)

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
                output_dir=config.cnf_results_root / config.model_name,
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
    sweep_name: str,
    variable_name: str,
    variable_value: int,
    domain_size: int,
    m_value: int,
    undirected_k_value: int,
    algorithm_name: str,
) -> dict[str, object]:
    """Run one Figure13 point with timeout handling."""
    q: Queue = Queue()

    process = Process(
        target=python_worker,
        args=(
            q,
            config,
            domain_size,
            m_value,
            undirected_k_value,
            algorithm_name,
        ),
    )
    process.start()

    process.join(config.timeout_seconds)
    if process.is_alive():
        process.terminate()
        process.join()
        result: object = "timeout"
    else:
        result = q.get() if not q.empty() else "error"

    status = _normalize_status(result)

    return {
        "timestamp": datetime.now(ZoneInfo(TIMEZONE)).strftime("%Y%m%d_%H%M%S"),
        "sweep": sweep_name,
        "variable_name": variable_name,
        "variable_value": variable_value,
        "formula": config.model_name,
        "domain_size": domain_size,
        "k_value": undirected_k_value,
        "m_value": m_value,
        "algorithm": algorithm_name,
        "model_count": result,
        "status": status,
    }


def _count_total_iterations(config: RuntimeConfig) -> int:
    total = 0
    for sweep in config.sweeps:
        if sweep.variable_name == "domain_size":
            total += (
                len(sweep.k_values)
                * len(sweep.m_values)
                * len(sweep.domain_sizes)
                * len(config.algorithms)
            )
        elif sweep.variable_name == "m_value":
            total += (
                len(sweep.domain_sizes)
                * len(sweep.k_values)
                * len(sweep.m_values)
                * len(config.algorithms)
            )
        elif sweep.variable_name == "k_value":
            total += (
                len(sweep.domain_sizes)
                * len(sweep.m_values)
                * len(sweep.k_values)
                * len(config.algorithms)
            )
        else:
            raise ValueError(f"Unsupported sweep variable: {sweep.variable_name}")
    return total


def _run_one_csv(
    config: RuntimeConfig,
    sweep: SweepSpec,
    csv_path: Path,
    variable_values: list[int],
    fixed_domain_size: int,
    fixed_m_value: int,
    fixed_k_value: int,
    progress_bar: tqdm,
    flush_every: int,
) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()

        for algorithm_name in config.algorithms:
            normalized_algorithm = str(algorithm_name).lower()
            skip_remaining = False

            for variable_value in variable_values:
                if skip_remaining:
                    progress_bar.update(1)
                    continue

                if sweep.variable_name == "domain_size":
                    domain_size = variable_value
                    m_value = fixed_m_value
                    undirected_k_value = _undirected_k_value(fixed_k_value)
                elif sweep.variable_name == "m_value":
                    domain_size = fixed_domain_size
                    m_value = variable_value
                    undirected_k_value = _undirected_k_value(fixed_k_value)
                elif sweep.variable_name == "k_value":
                    domain_size = fixed_domain_size
                    m_value = fixed_m_value
                    undirected_k_value = _undirected_k_value(variable_value)
                else:
                    raise ValueError(f"Unsupported sweep variable: {sweep.variable_name}")

                result_data = run_python_single(
                    config=config,
                    sweep_name=sweep.name,
                    variable_name=sweep.variable_name,
                    variable_value=variable_value,
                    domain_size=domain_size,
                    m_value=m_value,
                    undirected_k_value=undirected_k_value,
                    algorithm_name=normalized_algorithm,
                )

                writer.writerow(result_data)
                progress_bar.update(1)

                if progress_bar.n % flush_every == 0:
                    csvfile.flush()

                print_result(result_data, config.timeout_seconds)

                if _is_failure_status(str(result_data["status"])):
                    skip_remaining = True
                    print(
                        "Algorithm failed at "
                        f"{sweep.variable_name}={variable_value}; "
                        "skipping larger values for this series."
                    )


def run_experiment(
    config: RuntimeConfig,
    results_path: Path,
    flush_every: int = 20,
    disable_python_logging: bool = False,
) -> None:
    """Run all Figure13 sweeps and write CSV files + plots."""
    apply_python_logging_policy(disable_python_logging, update_env=True)

    total_iterations = _count_total_iterations(config)

    with tqdm(total=total_iterations, desc="Figure13 Model-Count Progress") as progress_bar:
        for sweep in config.sweeps:
            sweep_results_path = results_path / sweep.name / config.model_name
            sweep_results_path.mkdir(parents=True, exist_ok=True)

            if sweep.variable_name == "domain_size":
                for m_value in sweep.m_values:
                    for k_value in sweep.k_values:
                        csv_path = (
                            sweep_results_path
                            / f"{config.model_name}_{sweep.name}_m{m_value}_k{k_value}.csv"
                        )
                        _run_one_csv(
                            config=config,
                            sweep=sweep,
                            csv_path=csv_path,
                            variable_values=sweep.domain_sizes,
                            fixed_domain_size=0,
                            fixed_m_value=m_value,
                            fixed_k_value=k_value,
                            progress_bar=progress_bar,
                            flush_every=flush_every,
                        )
                        plot_model_count(str(csv_path), config=config)
                        print(
                            "\nSweep "
                            f"{sweep.name} (m={m_value}, k={k_value}) completed: {csv_path}"
                        )

            elif sweep.variable_name == "m_value":
                for domain_size in sweep.domain_sizes:
                    for k_value in sweep.k_values:
                        csv_path = (
                            sweep_results_path
                            / f"{config.model_name}_{sweep.name}_n{domain_size}_k{k_value}.csv"
                        )
                        _run_one_csv(
                            config=config,
                            sweep=sweep,
                            csv_path=csv_path,
                            variable_values=sweep.m_values,
                            fixed_domain_size=domain_size,
                            fixed_m_value=0,
                            fixed_k_value=k_value,
                            progress_bar=progress_bar,
                            flush_every=flush_every,
                        )
                        plot_model_count(str(csv_path), config=config)
                        print(
                            "\nSweep "
                            f"{sweep.name} (n={domain_size}, k={k_value}) completed: {csv_path}"
                        )

            elif sweep.variable_name == "k_value":
                for domain_size in sweep.domain_sizes:
                    for m_value in sweep.m_values:
                        csv_path = (
                            sweep_results_path
                            / f"{config.model_name}_{sweep.name}_n{domain_size}_m{m_value}.csv"
                        )
                        _run_one_csv(
                            config=config,
                            sweep=sweep,
                            csv_path=csv_path,
                            variable_values=sweep.k_values,
                            fixed_domain_size=domain_size,
                            fixed_m_value=m_value,
                            fixed_k_value=0,
                            progress_bar=progress_bar,
                            flush_every=flush_every,
                        )
                        plot_model_count(str(csv_path), config=config)
                        print(
                            "\nSweep "
                            f"{sweep.name} (n={domain_size}, m={m_value}) completed: {csv_path}"
                        )

            else:
                raise ValueError(f"Unsupported sweep variable: {sweep.variable_name}")
