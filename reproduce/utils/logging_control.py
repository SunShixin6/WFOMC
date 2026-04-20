from __future__ import annotations

import logging
import os
import sys
from importlib import import_module

_DEVNULL_STREAM = None

_TRUTHY = {"1", "true", "yes", "y", "on"}


def env_disables_python_logging() -> bool:
    """Return True when REPRO_DISABLE_PYTHON_LOGGING enables quiet mode."""
    value = os.environ.get("REPRO_DISABLE_PYTHON_LOGGING", "")
    return value.strip().lower() in _TRUTHY


def apply_python_logging_policy(
    disable_python_logging: bool,
    *,
    update_env: bool = False,
) -> None:
    """Apply reproducibility logging policy for std logging and loguru."""
    if update_env:
        os.environ["REPRO_DISABLE_PYTHON_LOGGING"] = "1" if disable_python_logging else "0"

    logging.disable(logging.CRITICAL if disable_python_logging else logging.NOTSET)

    # wfomc emits verbose logs through loguru. Silence its sinks when disabled.
    if not disable_python_logging:
        return

    try:
        loguru_module = import_module("loguru")
        loguru_logger = getattr(loguru_module, "logger", None)
    except Exception:
        return

    if loguru_logger is None:
        return

    try:
        loguru_logger.remove()
    except Exception:
        # No configured sink (or loguru internals changed) should not break runs.
        pass


def silence_process_stdio(disable_python_logging: bool) -> None:
    """Redirect child-process stdout/stderr to /dev/null in quiet mode."""
    if not disable_python_logging:
        return

    global _DEVNULL_STREAM
    if _DEVNULL_STREAM is None:
        _DEVNULL_STREAM = open(os.devnull, "w", encoding="utf-8")

    sys.stdout = _DEVNULL_STREAM
    sys.stderr = _DEVNULL_STREAM
