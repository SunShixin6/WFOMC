from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence
from zoneinfo import ZoneInfo

# Allow running this file directly via
# `python reproduce/correctness/odd_degree/main.py`.
if __package__ is None or __package__ == "":
    repo_root = Path(__file__).resolve().parents[3]
    for candidate in (repo_root, repo_root / "src"):
        candidate_str = str(candidate)
        if candidate_str not in sys.path:
            sys.path.insert(0, candidate_str)

from reproduce.correctness.odd_degree.benchmarks import (
    DEFAULT_DELTA,
    DEFAULT_EPSILON,
    DEFAULT_FLUSH_EVERY,
    DEFAULT_TIMEOUT_SECONDS,
    RuntimeConfig,
    TIMEZONE,
    build_runtime_config,
)
from reproduce.correctness.odd_degree.runner import run_experiment
from reproduce.utils.smoke_profile import (
    apply_smoke_correctness_odd_degree_config,
    is_smoke_mode,
)


def build_arg_parser() -> argparse.ArgumentParser:
    disable_logging_default = os.getenv("REPRO_DISABLE_PYTHON_LOGGING", "0") == "1"

    parser = argparse.ArgumentParser(
        description="Run Figure13 odd-degree model-count experiments."
    )
    parser.add_argument("--models-path", type=str, default=None)
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--delta", type=float, default=DEFAULT_DELTA)
    parser.add_argument("--flush-every", type=int, default=DEFAULT_FLUSH_EVERY)
    parser.add_argument(
        "--disable-python-logging",
        action="store_true",
        default=disable_logging_default,
        help="Disable Python logging output emitted by benchmark workers.",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Enable smoke profile (also enabled by REPRO_SMOKE=1).",
    )
    return parser


def _format_values(values: list[int]) -> str:
    if not values:
        return "[]"
    if len(values) == 1:
        return str(values[0])

    is_consecutive = all(
        values[idx + 1] == values[idx] + 1 for idx in range(len(values) - 1)
    )
    if is_consecutive:
        return f"[{values[0]}, {values[-1]}]"
    return "{" + ", ".join(str(value) for value in values) + "}"


def _print_setup(config: RuntimeConfig) -> None:
    print("Figure13 setup (from benchmarks.py):")
    print("  (paper uses undirected k; internal computation uses directed |E| = 2k)")
    for sweep in config.sweeps:
        if sweep.variable_name == "domain_size":
            print(
                f"  - {sweep.name}: n in {_format_values(sweep.domain_sizes)}, "
                f"m in {_format_values(sweep.m_values)}, "
                f"k (undirected) = {_format_values(sweep.k_values)}"
            )
        elif sweep.variable_name == "m_value":
            print(
                f"  - {sweep.name}: n = {_format_values(sweep.domain_sizes)}, "
                f"k (undirected) = {_format_values(sweep.k_values)}, "
                f"m in {_format_values(sweep.m_values)}"
            )
        elif sweep.variable_name == "k_value":
            print(
                f"  - {sweep.name}: n = {_format_values(sweep.domain_sizes)}, "
                f"m = {_format_values(sweep.m_values)}, "
                f"k (undirected) in {_format_values(sweep.k_values)}"
            )
        else:
            print(f"  - {sweep.name}: unsupported sweep variable {sweep.variable_name}")


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    models_path = Path(args.models_path).resolve() if args.models_path else None

    config = build_runtime_config(
        models_path=models_path,
        timeout_seconds=args.timeout_seconds,
        epsilon=args.epsilon,
        delta=args.delta,
    )

    if is_smoke_mode(args.smoke):
        config = apply_smoke_correctness_odd_degree_config(config)
        print("Smoke mode enabled for correctness.odd_degree.main")

    timestamp = datetime.now(ZoneInfo(TIMEZONE)).strftime("%Y%m%d_%H%M%S")
    results_path = config.csv_results_root / f"results_{timestamp}"

    results_path.mkdir(parents=True, exist_ok=True)
    config.full_fig_results_root.mkdir(parents=True, exist_ok=True)
    config.publish_fig_results_root.mkdir(parents=True, exist_ok=True)
    config.cnf_results_root.mkdir(parents=True, exist_ok=True)

    _print_setup(config)
    print(f"Figure13 model-count results will be saved in: {results_path}")

    run_experiment(
        config=config,
        results_path=results_path,
        flush_every=args.flush_every,
        disable_python_logging=args.disable_python_logging,
    )


if __name__ == "__main__":
    main()
