"""Shared utilities for paper reproduction workflows."""

from .figure_naming import (
    resolve_publish_figure_filename,
    resolve_publish_figure_stem,
)
from .model_resolver import get_reproduce_models_dir, resolve_model_file
from .smoke_profile import (
    apply_smoke_correctness_config,
    apply_smoke_correctness_odd_degree_config,
    apply_smoke_performance_config,
    apply_smoke_performance_odd_degree_config,
    is_smoke_mode,
    resolve_oeis_max_n,
)

__all__ = [
    "resolve_publish_figure_filename",
    "resolve_publish_figure_stem",
    "get_reproduce_models_dir",
    "resolve_model_file",
    "apply_smoke_correctness_config",
    "apply_smoke_correctness_odd_degree_config",
    "apply_smoke_performance_config",
    "apply_smoke_performance_odd_degree_config",
    "is_smoke_mode",
    "resolve_oeis_max_n",
]
