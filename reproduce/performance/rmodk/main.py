from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence
from zoneinfo import ZoneInfo

# Allow direct execution:
# python reproduce/performance/rmodk/main.py
if __package__ is None or __package__ == "":
    repo_root = Path(__file__).resolve().parents[3]
    for candidate in (repo_root, repo_root / "src"):
        candidate_str = str(candidate)
        if candidate_str not in sys.path:
            sys.path.insert(0, candidate_str)

from reproduce.performance.rmodk.benchmarks import (
    DEFAULT_DELTA,
    DEFAULT_EPSILON,
    DEFAULT_FLUSH_EVERY,
    DEFAULT_TIMEOUT_SECONDS,
    TIMEZONE,
    build_runtime_config,
)
from reproduce.performance.rmodk.runner import run_experiment
from reproduce.utils.smoke_profile import apply_smoke_performance_rmodk_config, is_smoke_mode


def build_arg_parser() -> argparse.ArgumentParser:
    disable_logging_default = os.getenv("REPRO_DISABLE_PYTHON_LOGGING", "1") != "0"
    verbose_default = not disable_logging_default

    parser = argparse.ArgumentParser(
        description="Run Figure 9 modulo-counting regular graph performance experiments."
    )
    parser.add_argument("--models-path", type=str, default=None)
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--delta", type=float, default=DEFAULT_DELTA)
    parser.add_argument("--flush-every", type=int, default=DEFAULT_FLUSH_EVERY)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--verbose",
        action="store_true",
        default=verbose_default,
        help="Enable Python logging output emitted by benchmark workers.",
    )
    group.add_argument(
        "--disable-python-logging",
        action="store_true",
        default=disable_logging_default,
        help="Disable Python logging output emitted by benchmark workers (default).",
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
        epsilon=args.epsilon,
        delta=args.delta,
    )

    if is_smoke_mode(args.smoke):
        config = apply_smoke_performance_rmodk_config(config)
        print("Smoke mode enabled for performance.rmodk.main")

    timestamp = datetime.now(ZoneInfo(TIMEZONE)).strftime("%Y%m%d_%H%M%S")
    results_path = config.csv_results_root / f"results_{timestamp}"

    results_path.mkdir(parents=True, exist_ok=True)
    config.full_fig_results_root.mkdir(parents=True, exist_ok=True)
    config.publish_fig_results_root.mkdir(parents=True, exist_ok=True)
    config.cnf_results_root.mkdir(parents=True, exist_ok=True)

    print(f"Rmodk performance results will be saved in: {results_path}")
    run_experiment(
        config=config,
        results_path=results_path,
        flush_every=args.flush_every,
        disable_python_logging=(args.disable_python_logging or not args.verbose),
    )


if __name__ == "__main__":
    main()
