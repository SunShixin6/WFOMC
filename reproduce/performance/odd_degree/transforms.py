from __future__ import annotations

import re
from pathlib import Path

from reproduce.utils.model_resolver import resolve_model_file

_ODD_CARDINALITY_PATTERN = re.compile(
    r"\\exists_\{=\d+\}\s*X:\s*\(Odd\(X\)\)"
)
_EDGE_CARDINALITY_PATTERN = re.compile(r"\|E\|\s*=\s*\d+")


def resolve_model_path(models_path: Path, model_filename: str) -> Path:
    """Resolve a model path and fail early with a clear error if missing."""
    return resolve_model_file(models_path, model_filename)


def load_model_template(model_path: Path) -> str:
    """Read the odd-degree model template from disk."""
    return model_path.read_text(encoding="utf-8")


def materialize_odd_degree_model(
    template: str,
    *,
    m_value: int,
    directed_k_value: int,
) -> str:
    """Inject m and directed |E| values into the odd-degree model template."""
    transformed, odd_sub_count = _ODD_CARDINALITY_PATTERN.subn(
        f"\\\\exists_{{={m_value}}} X: (Odd(X))",
        template,
        count=1,
    )
    transformed, edge_sub_count = _EDGE_CARDINALITY_PATTERN.subn(
        f"|E| = {directed_k_value}",
        transformed,
        count=1,
    )

    if odd_sub_count != 1:
        raise ValueError("Failed to locate odd-cardinality constraint in model template.")
    if edge_sub_count != 1:
        raise ValueError("Failed to locate edge-cardinality constraint in model template.")

    return transformed


def build_problem_for_wfomc(
    model_path: Path,
    *,
    domain_size: int,
    m_value: int,
    directed_k_value: int,
):
    """Build a wfomc problem object for a concrete odd-degree benchmark point."""
    from wfomc import Const
    from wfomc.parser.wfomcs_parser import parse as wfomcs_parse

    template = load_model_template(model_path)
    model_text = materialize_odd_degree_model(
        template,
        m_value=m_value,
        directed_k_value=directed_k_value,
    )
    problem = wfomcs_parse(model_text)
    problem.domain = {Const(f"d{i}") for i in range(domain_size)}
    return problem
