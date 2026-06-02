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
]


@dataclass(frozen=True)
class RuntimeConfig:
    """Runtime configuration for Figure 9 rmodk performance experiments."""

    dir_path: Path
    repo_root: Path
    models_path: Path
    csv_results_root: Path
    full_fig_results_root: Path
    publish_fig_results_root: Path
    cnf_results_root: Path
    timeout_seconds: int
    epsilon: float
    delta: float
    groups: list[dict[str, Any]]


def default_groups() -> list[dict[str, Any]]:
    """Figure 9 modulo-counting regular graph benchmarks."""
    return [
        {
            "name": "rmodk-regular-graphs",
            "domain_sizes": list(range(2, 20, 1)),
            "algorithms": ["incremental3", "ganak", "approxmc"],
            "models": {
                "0mod2-regular-graph": "0mod2-regular-graph.wfomcs",
                "1mod2-regular-graph": "1mod2-regular-graph.wfomcs",
                "2mod4-regular-graph": "2mod4-regular-graph-sc2.wfomcs",
            },
        },
    ]


def build_runtime_config(
    models_path: Path | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    epsilon: float = DEFAULT_EPSILON,
    delta: float = DEFAULT_DELTA,
) -> RuntimeConfig:
    """Build runtime config using repository-relative paths."""
    dir_path = Path(__file__).resolve().parent
    repo_root = dir_path.parents[2]

    resolved_models_path = models_path or (repo_root / "models")
    full_results_root = repo_root / "reproduce" / "performance" / "performance_results"
    csv_results_root = full_results_root / "raw_data" / "rmodk"
    full_fig_results_root = full_results_root / "Section6"
    publish_fig_results_root = repo_root / "reproduce" / "results" / "Section6"
    cnf_results_root = full_results_root / "raw_data" / "cnf" / "rmodk"

    return RuntimeConfig(
        dir_path=dir_path,
        repo_root=repo_root,
        models_path=resolved_models_path,
        csv_results_root=csv_results_root,
        full_fig_results_root=full_fig_results_root,
        publish_fig_results_root=publish_fig_results_root,
        cnf_results_root=cnf_results_root,
        timeout_seconds=timeout_seconds,
        epsilon=epsilon,
        delta=delta,
        groups=default_groups(),
    )
