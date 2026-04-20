from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

TIMEZONE = "Asia/Shanghai"
MEMORY_POLL_INTERVAL = 0.05
DEFAULT_FLUSH_EVERY = 1
DEFAULT_TIMEOUT_SECONDS = 10000
DEFAULT_EPSILON = 0.05
DEFAULT_DELTA = 0.1

CSV_FIELDNAMES = [
    "timestamp",
    "formula",
    "domain_size",
    "algorithm",
    "result",
    "time_sec",
    "memory_bytes",
    "memory_kb",
    "memory_mb",
    "status",
    "incremental3_within_approx_confidence",
]


@dataclass(frozen=True)
class RuntimeConfig:
    """Runtime configuration for correctness checking experiments."""

    dir_path: Path
    repo_root: Path
    models_path: Path
    raw_data_results_path: Path
    appendix_b1_results_path: Path
    csv_results_path: Path
    cnf_results_root: Path
    timeout_seconds: int
    epsilon: float
    delta: float
    groups: list[dict[str, Any]]


def default_groups() -> list[dict[str, Any]]:
    """Default correctness groups migrated from experiment/check/experiment_script.py."""
    return [
        {
            "name": "0mod2-regular-graph",
            "domain_sizes": list(range(2, 9, 1)),
            "algorithms": ["incremental3", "ganak", "approxmc"],
            "models": {
                "0mod2-regular-graph": "0mod2-regular-graph.wfomcs",
            },
        },
        {
            "name": "1mod2-regular-graph",
            # All vertices are constrained to odd degree; by handshaking lemma,
            # odd domain sizes are unsat (model count 0). Keep step=1 when you
            # want to show these zero points, or use step=2 to sample even n only.
            "domain_sizes": list(range(2, 9, 1)),
            "algorithms": ["incremental3", "ganak", "approxmc"],
            "models": {
                "1mod2-regular-graph": "1mod2-regular-graph.wfomcs",
            },
        },
        {
            "name": "2mod4-regular-graph",
            "domain_sizes": list(range(3, 9, 1)),
            "algorithms": ["incremental3", "ganak", "approxmc"],
            "models": {
                "2mod4-regular-graph": "2mod4-regular-graph.wfomcs",
            },
        },
        {
            "name": "2-regular-graph",
            "domain_sizes": list(range(3, 9, 1)),
            "algorithms": ["incremental3", "ganak", "approxmc"],
            "models": {
                "2-regular-graph": "2-regular-graph-sc2.wfomcs",
            },
        },
    {
            "name": "3-regular-graph",
            "domain_sizes": list(range(4, 9, 1)),
            "algorithms": ["incremental3", "ganak", "approxmc"],
            "models": {
                "3-regular-graph": "3-regular-graph.wfomcs",
            },
        },
        {
            "name": "4-regular-graph",
            "domain_sizes": list(range(5, 9, 1)),
            "algorithms": ["incremental3", "ganak", "approxmc"],
            "models": {
                "4-regular-graph": "4-regular-graph.wfomcs",
            },
        },
    ]


def build_runtime_config(
    models_path: Path | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    epsilon: float = DEFAULT_EPSILON,
    delta: float = DEFAULT_DELTA,
) -> RuntimeConfig:
    """Build runtime config with repository-relative defaults."""
    dir_path = Path(__file__).resolve().parent
    repo_root = dir_path.parents[1]

    resolved_models_path = models_path or (repo_root / "models")
    correctness_results_root = dir_path / "correctness_results"
    raw_data_results_path = correctness_results_root / "raw_data"
    appendix_b1_results_path = correctness_results_root / "AppendixB.1"
    csv_results_path = raw_data_results_path
    cnf_results_root = raw_data_results_path / "cnf"

    return RuntimeConfig(
        dir_path=dir_path,
        repo_root=repo_root,
        models_path=resolved_models_path,
        raw_data_results_path=raw_data_results_path,
        appendix_b1_results_path=appendix_b1_results_path,
        csv_results_path=csv_results_path,
        cnf_results_root=cnf_results_root,
        timeout_seconds=timeout_seconds,
        epsilon=epsilon,
        delta=delta,
        groups=default_groups(),
    )
