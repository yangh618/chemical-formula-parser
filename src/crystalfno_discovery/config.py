"""Configuration loading and project-path resolution.

The project keeps a single human-editable :file:`config.yaml` at the project
root. This module locates that root, loads the YAML (merging it over a
built-in default dict), and exposes the result as a plain nested mapping so it
can be read from any module without hard-coding values throughout the source.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml

# Name of the top-level configuration file.
CONFIG_FILENAME = "config.yaml"

# Environment variable that may point at an alternative configuration file.
CONFIG_ENV = "CRYSTALFNO_DISCOVERY_CONFIG"

# ---------------------------------------------------------------------------
# Built-in defaults. They mirror config.yaml so that the package is usable
# even when the config file is missing.
# ---------------------------------------------------------------------------
DEFAULTS: dict[str, Any] = {
    "elements": ["Al", "Cr", "Co", "Mo", "Nb", "Ta", "Ti", "V", "W", "Zr"],
    "smact": {
        "num_elements": [2, 3, 4],
        "max_stoich": 8,
        "oxidation_states": "icsd24",
        "charge_neutrality": True,
        "electronegativity": False,
        "include_alloys": True,
    },
    "filters": {
        "min_elements": 2,
        "max_elements": 6,
        "max_atoms": 200,
        "remove_partial_occupancy": True,
        "remove_disordered": True,
    },
    "databases": {"cod": True, "oqmd": True, "jarvis": True, "aflow": False},
    "materials_project": {"metadata_file": None},
    "network": {
        "timeout_seconds": 60,
        "retries": 5,
        "backoff_factor": 1.5,
        "rate_limit_seconds": 0.35,
        "page_limit": 100,
        "endpoints": {
            "cod_cif": "https://www.crystallography.net/cod/{cod_id}.cif",
            "cod_optimade": "https://www.crystallography.net/cod/optimade/v1/structures",
            "oqmd_optimade": "https://oqmd.org/optimade/v1/structures",
            "oqmd_rest": "https://oqmd.org/oqmdapi/formationenergy",
            "jarvis_optimade": "https://jarvis.nist.gov/optimade/jarvisdft/v1/structures",
            "aflow_aflux": "https://aflow.org/API/aflux/",
        },
    },
    "ranking": {"enabled": True, "criteria": ["n_elements", "stoich_complexity", "tm_count", "smact_validity_strength"]},
}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge ``override`` over ``base`` (both plain dicts)."""
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def get_project_root() -> Path:
    """Return the project root directory (the one holding ``config.yaml``)."""
    # This module lives at <root>/src/crystalfno_discovery/config.py
    here = Path(__file__).resolve()
    for candidate in (here.parents[2], here.parents[3], here.parents[4]):
        if (candidate / CONFIG_FILENAME).exists():
            return candidate
    # Fall back to the repository layout root.
    return here.parents[3]


def resolve_config_path() -> Path:
    """Return the path to the active configuration file."""
    env_path = os.environ.get(CONFIG_ENV)
    if env_path:
        return Path(env_path).expanduser().resolve()
    return get_project_root() / CONFIG_FILENAME


def load_config() -> dict[str, Any]:
    """Load configuration, merging ``config.yaml`` over the built-in defaults."""
    config: dict[str, Any] = _deep_merge(DEFAULTS, {})
    cfg_path = resolve_config_path()
    if cfg_path.exists():
        with cfg_path.open("r", encoding="utf-8") as handle:
            user_cfg = yaml.safe_load(handle) or {}
        config = _deep_merge(config, user_cfg)
    return config


# Module-level singleton, populated lazily so callers can simply
# ``from crystalfno_discovery.config import CFG``.
CFG: dict[str, Any] = load_config()


def get(key: str, default: Any = None) -> Any:
    """Dot-path access into the singleton config, e.g. ``get("smact.max_stoich")``."""
    node: Any = CFG
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node
