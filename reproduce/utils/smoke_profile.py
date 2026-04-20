from __future__ import annotations

import os
from dataclasses import replace
from typing import Any, Sequence

_TRUTHY = {"1", "true", "yes", "y", "on"}


def is_smoke_mode(cli_flag: bool = False) -> bool:
    """Return True when smoke mode is enabled via CLI or environment."""
    if cli_flag:
        return True
    value = os.environ.get("REPRO_SMOKE", "")
    return value.strip().lower() in _TRUTHY


def _take(values: Sequence[int], count: int) -> list[int]:
    sliced = list(values)[:count]
    return sliced if sliced else list(values)


def _prefer_incremental3(algorithms: Sequence[str]) -> list[str]:
    normalized = [str(algo) for algo in algorithms]
    for preferred in ("incremental3", "dr"):
        for algo in normalized:
            if algo.lower() == preferred:
                return [algo]
    return normalized[:1]


def apply_smoke_correctness_config(config: Any) -> Any:
    """Apply lightweight overrides for reproduce.correctness.main."""
    smoke_groups = []
    for group in config.groups:
        new_group = dict(group)
        new_group["domain_sizes"] = _take(group["domain_sizes"], 2)
        new_group["algorithms"] = _prefer_incremental3(group["algorithms"])
        smoke_groups.append(new_group)

    return replace(
        config,
        timeout_seconds=min(int(config.timeout_seconds), 120),
        groups=smoke_groups,
    )


def apply_smoke_correctness_odd_degree_config(config: Any) -> Any:
    """Apply lightweight overrides for reproduce.correctness.odd_degree.main."""
    smoke_sweeps = []
    for sweep in config.sweeps:
        if sweep.variable_name == "domain_size":
            smoke_sweeps.append(
                replace(
                    sweep,
                    domain_sizes=_take(sweep.domain_sizes, 3),
                    m_values=_take(sweep.m_values, 1),
                    k_values=_take(sweep.k_values, 1),
                )
            )
        elif sweep.variable_name == "k_value":
            smoke_sweeps.append(
                replace(
                    sweep,
                    domain_sizes=_take(sweep.domain_sizes, 1),
                    m_values=_take(sweep.m_values, 1),
                    k_values=_take(sweep.k_values, 3),
                )
            )
        elif sweep.variable_name == "m_value":
            smoke_sweeps.append(
                replace(
                    sweep,
                    domain_sizes=_take(sweep.domain_sizes, 1),
                    m_values=_take(sweep.m_values, 2),
                    k_values=_take(sweep.k_values, 1),
                )
            )
        else:
            smoke_sweeps.append(sweep)

    return replace(
        config,
        timeout_seconds=min(int(config.timeout_seconds), 120),
        algorithms=_prefer_incremental3(config.algorithms),
        sweeps=smoke_sweeps,
    )


def apply_smoke_performance_config(config: Any) -> Any:
    """Apply lightweight overrides for reproduce.performance.main."""
    smoke_groups = []
    for group in config.groups:
        new_group = dict(group)
        new_group["domain_size"] = _take(group["domain_size"], 2)
        new_group["algorithms"] = _prefer_incremental3(group["algorithms"])
        # Keep all models so smoke still covers every benchmark type.
        new_group["models"] = dict(group["models"])
        smoke_groups.append(new_group)

    return replace(
        config,
        timeout_seconds=min(int(config.timeout_seconds), 120),
        groups=smoke_groups,
    )


def apply_smoke_performance_odd_degree_config(config: Any) -> Any:
    """Apply lightweight overrides for reproduce.performance.odd_degree.main."""
    return replace(
        config,
        timeout_seconds=min(int(config.timeout_seconds), 120),
        domain_sizes=_take(config.domain_sizes, 2),
        # Keep all m groups so Figure_10 a/b/c are all produced in smoke mode.
        m_values=list(config.m_values),
        algorithms=_prefer_incremental3(config.algorithms),
    )


def resolve_oeis_max_n(
    cli_max_n: int | None,
    cli_smoke: bool,
    default_max_n: int = 10,
    smoke_max_n: int = 4,
) -> int:
    """Resolve max_n for OEIS generation, honoring smoke mode and explicit override."""
    if cli_max_n is not None:
        return cli_max_n
    return smoke_max_n if is_smoke_mode(cli_smoke) else default_max_n
