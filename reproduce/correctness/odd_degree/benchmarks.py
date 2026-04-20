from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

TIMEZONE = "Asia/Shanghai"
DEFAULT_FLUSH_EVERY = 1
DEFAULT_TIMEOUT_SECONDS = 10000
DEFAULT_EPSILON = 0.05
DEFAULT_DELTA = 0.1

MODEL_NAME = "m-odd-degree-graph-sc2"
MODEL_FILENAME = "m-odd-degree-graph-sc2.wfomcs"
DEFAULT_ALGORITHMS = ["incremental3", "ganak", "approxmc"]

CSV_FIELDNAMES = [
    "timestamp",
    "sweep",
    "variable_name",
    "variable_value",
    "formula",
    "domain_size",
    "k_value",
    "m_value",
    "algorithm",
    "model_count",
    "status",
]


@dataclass(frozen=True)
class SweepSpec:
    """Define one parameter sweep used in Figure13 reproduction."""

    name: str
    variable_name: str
    domain_sizes: list[int]
    m_values: list[int]
    k_values: list[int]


@dataclass(frozen=True)
class RuntimeConfig:
    """Runtime configuration for Figure13 model-count experiments."""

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
    model_name: str
    model_filename: str
    algorithms: list[str]
    sweeps: list[SweepSpec]


def default_sweeps() -> list[SweepSpec]:
    """Three Figure13 sweeps: vary n, vary m, and vary k."""
    return [
        SweepSpec(
            name="vary_n",
            variable_name="domain_size",
            domain_sizes=list(range(2, 9)),
            m_values=[2],
            # Paper-facing k value (undirected |E|); used exactly as configured.
            k_values=[1],
        ),
        SweepSpec(
            name="vary_k",
            variable_name="k_value",
            domain_sizes=[5],
            m_values=[2],
            # Paper-facing k values (undirected |E|).
            k_values=list(range(1, 10)),
        ),
        SweepSpec(
            name="vary_m",
            variable_name="m_value",
            domain_sizes=[8],
            m_values=list(range(0, 7, 2)),
            # Paper-facing k value (undirected |E|).
            k_values=[3],
        ),
    ]


def build_runtime_config(
    models_path: Path | None = None,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    epsilon: float = DEFAULT_EPSILON,
    delta: float = DEFAULT_DELTA,
) -> RuntimeConfig:
    """Build repository-relative runtime config for Figure13 experiments."""
    dir_path = Path(__file__).resolve().parent
    repo_root = dir_path.parents[2]

    resolved_models_path = models_path or (repo_root / "models")
    full_results_root = repo_root / "reproduce" / "correctness" / "correctness_results"
    raw_data_results_root = full_results_root / "raw_data"
    csv_results_root = raw_data_results_root / "odd_degree"
    full_fig_results_root = full_results_root / "AppendixB.1"
    publish_fig_results_root = repo_root / "reproduce" / "results" / "AppendixB.1"
    cnf_results_root = raw_data_results_root / "cnf" / "odd_degree"

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
        model_name=MODEL_NAME,
        model_filename=MODEL_FILENAME,
        algorithms=DEFAULT_ALGORITHMS,
        sweeps=default_sweeps(),
    )
