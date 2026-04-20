from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence
from zoneinfo import ZoneInfo

# Allow running this file directly via
# `python reproduce/performance/odd_degree/main.py`.
if __package__ is None or __package__ == "":
    repo_root = Path(__file__).resolve().parents[3]
    for candidate in (repo_root, repo_root / "src"):
        candidate_str = str(candidate)
        if candidate_str not in sys.path:
            sys.path.insert(0, candidate_str)

from reproduce.performance.odd_degree.benchmarks import DEFAULT_FLUSH_EVERY, TIMEZONE, build_runtime_config
from reproduce.performance.odd_degree.runner import run_experiment
from reproduce.utils.smoke_profile import (
    apply_smoke_performance_odd_degree_config,
    is_smoke_mode,
)


def build_arg_parser() -> argparse.ArgumentParser:
    disable_logging_default = os.getenv("REPRO_DISABLE_PYTHON_LOGGING", "0") == "1"

    parser = argparse.ArgumentParser(
        description="Run fixed odd-degree benchmark experiments (m=2/4/6, directed k=2n)."
    )
    parser.add_argument("--models-path", type=str, default=None)
    parser.add_argument("--timeout-seconds", type=int, default=100)
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


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    models_path = Path(args.models_path).resolve() if args.models_path else None

    config = build_runtime_config(
        models_path=models_path,
        timeout_seconds=args.timeout_seconds,
    )

    if is_smoke_mode(args.smoke):
        config = apply_smoke_performance_odd_degree_config(config)
        print("Smoke mode enabled for performance.odd_degree.main")

    timestamp = datetime.now(ZoneInfo(TIMEZONE)).strftime("%Y%m%d_%H%M%S")
    results_path = config.csv_results_root / f"results_{timestamp}"

    results_path.mkdir(parents=True, exist_ok=True)
    config.full_fig_results_root.mkdir(parents=True, exist_ok=True)
    config.publish_fig_results_root.mkdir(parents=True, exist_ok=True)
    config.cnf_results_root.mkdir(parents=True, exist_ok=True)

    print(
        "Fixed odd-degree setup: "
        f"n in {config.domain_sizes}, directed k={config.k_multiplier}*n, "
        f"m in {config.m_values}, algorithms={','.join(map(str, config.algorithms))}"
    )
    print(f"Odd-degree experiment results will be saved in: {results_path}")
    run_experiment(
        config=config,
        results_path=results_path,
        flush_every=args.flush_every,
        disable_python_logging=args.disable_python_logging,
    )


if __name__ == "__main__":
    main()
