from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

TIMEZONE = "Asia/Shanghai"
MEMORY_POLL_INTERVAL = 0.05
DEFAULT_FLUSH_EVERY = 1
DEFAULT_EPSILON = 0.05
DEFAULT_DELTA = 0.1
# DEFAULT_TIMEOUT_SECONDS = 10000
DEFAULT_TIMEOUT_SECONDS = 100

# Fixed paper-style odd-degree setup.
FIXED_DOMAIN_SIZES = list(range(5, 10))
FIXED_M_VALUES = [2, 4, 6]
FIXED_ALGORITHMS = ["incremental3", "ganak", "approxmc"]
FIXED_K_MULTIPLIER = 2
FIXED_GROUP_NAME = "odd-degree-n"
FIXED_MODEL_NAME = "m-odd-degree-graph-sc2"
FIXED_MODEL_FILENAME = "m-odd-degree-graph-sc2.wfomcs"

CSV_FIELDNAMES = [
    "timestamp",
    "formula",
    "domain_size",
    "k_value",
    "m_value",
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
    """Runtime configuration for fixed odd-degree benchmark experiments."""

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
    group_name: str
    model_name: str
    model_filename: str
    domain_sizes: list[int]
    m_values: list[int]
    algorithms: list[str]
    k_multiplier: int


def build_runtime_config(
    models_path: Path | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    epsilon: float = DEFAULT_EPSILON,
    delta: float = DEFAULT_DELTA,
) -> RuntimeConfig:
    """Build fixed odd-degree runtime config using repository-relative paths."""
    dir_path = Path(__file__).resolve().parent
    repo_root = dir_path.parents[2]

    resolved_models_path = models_path or (repo_root / "models")
    full_results_root = repo_root / "reproduce" / "performance" / "performance_results"
    csv_results_root = full_results_root / "raw_data" / "odd_degree"
    full_fig_results_root = full_results_root / "Section6"
    publish_fig_results_root = repo_root / "reproduce" / "results" / "Section6"
    cnf_results_root = full_results_root / "raw_data" / "cnf" / "odd_degree"

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
        group_name=FIXED_GROUP_NAME,
        model_name=FIXED_MODEL_NAME,
        model_filename=FIXED_MODEL_FILENAME,
        domain_sizes=FIXED_DOMAIN_SIZES,
        m_values=FIXED_M_VALUES,
        algorithms=FIXED_ALGORITHMS,
        k_multiplier=FIXED_K_MULTIPLIER,
    )
