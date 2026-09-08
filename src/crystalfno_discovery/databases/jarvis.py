"""JARVIS-DFT client via its OPTIMADE endpoint.

The JARVIS-DFT OPTIMADE server
(``https://jarvis.nist.gov/optimade/jarvisdft/v1/structures``, provider
``NIST-JARVIS``) returns full structures but, unlike OQMD/COD, **does not
enumerate the database on an unfiltered query** (an unfiltered request returns
zero rows). Enumeration must therefore be driven by ``elements`` filters.

This client enumerates by element combination for a user-supplied element
universe (default: the ``elements`` list in ``config.yaml``). A complete
``--all`` sweep therefore means *every* element combination up to a chosen
size over that universe, which is the documented, robust way to poll JARVIS.
"""

from __future__ import annotations

import itertools
import logging
from typing import Any

from .. import database_utils as du
from ..config import get
from . import JARVIS, save_optimade_structure, structures_dir, raw_dir

logger = logging.getLogger(__name__)

SOURCE_NAME = "NIST JARVIS-DFT (Joint Automated Repository for Various Integrated Simulations)"
SOURCE_URL = "https://jarvis.nist.gov/"


def optimade_endpoint() -> str:
    return get("network.endpoints.jarvis_optimade", "https://jarvis.nist.gov/optimade/jarvisdft/v1/structures")


def elements_filter(combo: tuple[str, ...]) -> str:
    """OPTIMADE ``elements HAS ALL ...`` filter for an element combination."""
    quoted = ",".join(f'"{el}"' for el in combo)
    return f"elements HAS ALL {quoted}"


def iter_filter_exprs(elements: list[str], sizes: tuple[int, ...]) -> Any:
    """Yield element-combination filter expressions for all sizes requested."""
    for size in sizes:
        if size < 1:
            continue
        for combo in itertools.combinations(sorted(elements), size):
            yield elements_filter(combo)


def collect(
    client: du.NetClient,
    elements: list[str],
    *,
    sizes: tuple[int, ...] = (1,),
    max_pages: int | None = None,
    page_limit: int | None = None,
    max_filters: int | None = None,
) -> dict[str, Any]:
    """Enumerate JARVIS-DFT by element filters and save structures as CIFs.

    A shared ``max_pages`` budget spans *all* filters so that ``--max-pages
    10`` really means ten HTTP pages in total (cheap test). Records whose CIF
    already exists are skipped, keeping interrupted runs resumable.
    """
    page_limit = int(page_limit or get("network.page_limit", 100))
    summary = {"downloaded": 0, "skipped": 0, "failed": 0, "rows": [], "status_counts": {}, "filters_used": 0}
    pages_used = 0
    filter_count = 0

    for filter_expr in iter_filter_exprs(elements, tuple(sizes)):
        if max_filters is not None and filter_count >= max_filters:
            break
        filter_count += 1
        summary["filters_used"] = filter_count
        offset = 0
        while True:
            if max_pages is not None and pages_used >= max_pages:
                return summary
            params: dict[str, Any] = {"page_limit": page_limit, "page_offset": offset, "filter": filter_expr}
            payload = client.get_json(optimade_endpoint(), params=params)
            meta = payload.get("meta", {}) or {}
            records = payload.get("data", []) or []
            pages_used += 1
            for record in records:
                attributes = record.get("attributes", {}) or {}
                db_id = str(record.get("id", ""))
                if not db_id:
                    summary["failed"] += 1
                    continue
                is_new, metadata, status = save_optimade_structure(
                    attributes, JARVIS, db_id,
                    struct_dir=structures_dir(JARVIS), raw_dir_path=raw_dir(JARVIS), save_raw=True,
                )
                summary["status_counts"][status] = summary["status_counts"].get(status, 0) + 1
                if is_new:
                    summary["downloaded"] += 1
                elif status == du.STATUS_DUPLICATE:
                    summary["skipped"] += 1
                else:
                    summary["failed"] += 1
                if metadata is not None:
                    metadata["database"] = JARVIS
                    metadata["database_id"] = db_id
                    summary["rows"].append(metadata)
            if not meta.get("more_data_available", False):
                break
            if not records:
                break
            offset += page_limit
    return summary
