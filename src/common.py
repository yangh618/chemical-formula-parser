"""Shared logging and I/O helpers used across scripts and modules."""

from __future__ import annotations

import logging
import time
from pathlib import Path


class RateLimiter:
    """Minimal inter-request pause to be polite to public APIs.

    Also exposes a token/reset mechanism so callers can throttle aggressive
    pagination loops.
    """

    def __init__(self, min_interval: float = 0.0) -> None:
        self.min_interval = float(min_interval)
        self._last = 0.0

    def wait(self) -> None:
        if self.min_interval <= 0.0:
            return
        now = time.monotonic()
        delay = self._last + self.min_interval - now
        if delay > 0.0:
            time.sleep(delay)
        self._last = time.monotonic()


def configure_logging(
    name: str = "chemical_formula_parser",
    log_dir: str | Path = "logs",
    log_filename: str | None = None,
    level: int = logging.INFO,
) -> logging.Logger:
    """Return a logger writing to both the console and a per-script log file.

    Logs are stored under ``logs/`` (section 22). Idempotent: repeated calls
    for the same logger name do not duplicate handlers.
    """
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(level)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(name)s | %(message)s")

    console = logging.StreamHandler()
    console.setLevel(level)
    console.setFormatter(fmt)
    logger.addHandler(console)

    if log_filename:
        path = Path(log_dir)
        path.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(path / log_filename, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    logger.propagate = False
    return logger


def get_logger(name: str) -> logging.Logger:
    """Fetch a logger configured via :func:`configure_logging` or plain."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        # Add a default console handler so module-level logging never vanishes.
        configure_logging(name=name)
    return logger
