"""OQMD (Open Quantum Materials Database) client via its OPTIMADE endpoint.

OQMD exposes full structures (lattice + atomic positions) through OPTIMADE at
``https://oqmd.org/optimade/v1/structures`` (verified live, schema v1.2.0,
pagination with ``page_limit`` / ``page_offset``). A native REST fallback at
``https://oqmd.org/oqmdapi/formationenergy`` is documented but the OPTIMADE
endpoint returns the structure coordinates we need for CIF export.
"""

from __future__ import annotations

from typing import Any

from ..config import get
from . import OQMD, collect_optimade

SOURCE_NAME = "Open Quantum Materials Database (OQMD)"
SOURCE_URL = "https://oqmd.org"


def optimade_endpoint() -> str:
    return get("network.endpoints.oqmd_optimade", "https://oqmd.org/optimade/v1/structures")


def collect(
    client: Any,
    *,
    max_pages: int | None = None,
    page_limit: int | None = None,
) -> dict[str, Any]:
    """Enumerate OQMD structures through OPTIMADE and save them as CIFs.

    ``max_pages`` is used for cheap testing (e.g. ``--max-pages 10``); pass
    ``max_pages=None`` with ``--all`` for a complete extraction.
    """
    summary = collect_optimade(
        client, OQMD, optimade_endpoint(), max_pages=max_pages, page_limit=page_limit, save_raw=True
    )
    return summary
