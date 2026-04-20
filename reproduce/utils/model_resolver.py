from __future__ import annotations

from functools import lru_cache
from pathlib import Path


def get_reproduce_models_dir() -> Path:
    """Return the repository-local fallback model directory: reproduce/models."""
    repo_root = Path(__file__).resolve().parents[2]
    return repo_root / "reproduce" / "models"


@lru_cache(maxsize=None)
def _find_nested_model(base_dir: Path, model_filename: str) -> Path | None:
    """Find the first matching model file under nested subdirectories."""
    if not base_dir.exists() or not base_dir.is_dir():
        return None

    matches = sorted(path for path in base_dir.rglob(model_filename) if path.is_file())
    if not matches:
        return None

    return matches[0]


def _resolve_from_base_dir(base_dir: Path, model_filename: str) -> Path | None:
    """Resolve model file from a base dir via direct path, then nested search."""
    direct = base_dir / model_filename
    if direct.exists() and direct.is_file():
        return direct

    return _find_nested_model(base_dir, model_filename)


def resolve_model_file(models_path: Path, model_filename: str) -> Path:
    """Resolve model file from primary models path, then reproduce/models fallback."""
    primary_resolved = _resolve_from_base_dir(models_path, model_filename)
    if primary_resolved is not None:
        return primary_resolved

    fallback_dir = get_reproduce_models_dir()
    fallback_resolved = _resolve_from_base_dir(fallback_dir, model_filename)
    if fallback_resolved is not None:
        return fallback_resolved

    raise FileNotFoundError(
        "Model file not found. "
        f"Searched under: {models_path} (direct and nested) and {fallback_dir} "
        "(direct and nested). "
        "Please place the file under repo-level models/ (flat or subfolders) "
        "or reproduce/models/ (flat or subfolders)."
    )
