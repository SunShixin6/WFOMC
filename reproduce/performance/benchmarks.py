from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

TIMEZONE = "Asia/Shanghai"
MEMORY_POLL_INTERVAL = 0.05
DEFAULT_FLUSH_EVERY = 1
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
    """Runtime experiment configuration copied from the original performance script."""

    dir_path: Path
    csv_results_path: Path
    models_path: Path
    timeout_seconds: int
    groups: list[dict[str, Any]]
    full_fig_results_path: Path
    publish_fig_results_path: Path


def default_groups() -> list[dict[str, Any]]:
    """Benchmark groups preserved from experiment/performance/ouralgo_performance_script.py."""
    return [
        {
            "name": "regular-graphs",
            "domain_size": list(range(2, 100, 3)),
            # "domain_size": list(range(2, 20, 3)),  # test
            "algorithms": ["recursive", "fast", "incremental3"],
            "models": {
                "3-regular-graph": "3-regular-graph.wfomcs",
                "4-regular-graph": "4-regular-graph.wfomcs",
                "5-regular-graph": "5-regular-graph.wfomcs",
            },
        },
        {
            "name": "colored-graphs",
            "domain_size": list(range(2, 100, 3)),
            # "domain_size": list(range(2, 20, 3)),  # test
            "algorithms": ["recursive", "fast", "incremental3"],
            "models": {
                "3-regular-2-colored-graph": "3-regular-2-colored-graph.wfomcs",
                "4-regular-2-colored-graph": "4-regular-2-colored-graph.wfomcs",
                "5-regular-2-colored-graph": "5-regular-2-colored-graph.wfomcs",
            },
        },
        # Since the following domains are more detailed, they are tested separately
        {
            "name": "3-regular-3-colored-graph",
            "domain_size": list(range(2, 40, 1)),
            # "domain_size": list(range(2, 10, 1)),  # test
            "algorithms": ["recursive", "fast", "incremental3"],
            "models": {
                "3-regular-3-colored-graph": "3-regular-3-colored-graph.wfomcs",
            },
        },
        {
            "name": "3-regular-4-colored-graph",
            "domain_size": list(range(2, 30, 1)),
            # "domain_size": list(range(2, 10, 1)),
            "algorithms": ["recursive", "fast", "incremental3"],
            "models": {
                "3-regular-4-colored-graph": "3-regular-4-colored-graph.wfomcs",
            },
        },
        {
            "name": "4-regular-3-colored-graph",
            "domain_size": list(range(2, 40, 1)),
            # "domain_size": list(range(2, 10, 1)),
            "algorithms": ["recursive", "fast", "incremental3"],
            "models": {
                "4-regular-3-colored-graph": "4-regular-3-colored-graph.wfomcs",
            },
        },
        {
            "name": "directed-graphs",
            "domain_size": list(range(2, 50, 3)),
            # "domain_size": list(range(2, 20, 3)),  # test
            "algorithms": ["recursive", "fast", "incremental3"],
            "models": {
                "2-regular-directed-graph": "2-regular-directed-graph.wfomcs",
                "3-regular-directed-graph": "3-regular-directed-graph.wfomcs",
            },
        },
        {
            "name": "BA",  
            "domain_size": list(range(2,100, 3)), 
            # "domain_size": list(range(2, 20, 3)), 
            "algorithms": ["recursive", "incremental3"],  
            "models": {
                "BA_CC": "BA_CC.wfomcs",
                "BA": "BA.wfomcs",
            },
        },
        {
            "name": "modk-regular-graphs",
            # Even domains avoid trivial UNSAT for odd-degree parity families.
            "domain_size": list(range(2, 40, 2)),
            "algorithms": ["recursive", "fast", "incremental3"],
            "models": {
                "0mod2-regular-graph": "0mod2-regular-graph.wfomcs",
                "1mod2-regular-graph": "1mod2-regular-graph.wfomcs",
                "2mod4-regular-graph": "2mod4-regular-graph-sc2.wfomcs",
            },
        },
    ]


def build_runtime_config(
    models_path: Path | None = None,
    timeout_seconds: int = 10000,
) -> RuntimeConfig:
    """
    Build runtime config. Compared with the original script's absolute models path,
    this uses a repository-relative path for safer portability.
    """
    dir_path = Path(__file__).resolve().parent
    repo_root = dir_path.parents[1]
    resolved_models_path = models_path or (repo_root / "models")
    full_results_root = dir_path / "performance_results"
    csv_results_path = full_results_root / "raw_data"
    full_fig_results_path = full_results_root
    publish_fig_results_path = repo_root / "reproduce" / "results"

    return RuntimeConfig(
        dir_path=dir_path,
        models_path=resolved_models_path,
        timeout_seconds=timeout_seconds,
        groups=default_groups(),
        csv_results_path=csv_results_path,
        full_fig_results_path=full_fig_results_path,
        publish_fig_results_path=publish_fig_results_path,
    )
