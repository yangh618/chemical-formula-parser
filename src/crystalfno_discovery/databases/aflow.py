"""AFLOW client.

AFLOW's public, credential-free API is AFLUX (``https://aflow.org/API/aflux/``).
After inspecting it (help/paging/schema), this package uses AFLUX **only for
composition / metadata discovery**: ``compound`` strings are enumerated and
canonicalised into known formulas for the novelty screen.

Limitation (documented, not hidden): resolving full atomic structures
(positions + lattice) from AFLUX in a stable way requires pulling `catalog`
records and aflow entries that are large and frequently reorganised. Rather
than risk emitting broken structures, AFLOW is treated as a *formula-source
plugin only*. Structure downloads from AFLOW are therefore not implemented;
``resolve_structure`` returns ``None`` with an explanatory note. The plugin
interface is kept so a full AFLOW structure source can be added later without
touching the rest of the pipeline.

AFLOW data are free for scientific, academic and non-commercial purposes;
commercial use is prohibited (see AFLUX terms). Redistribution follows the
provider licence - store derived metadata only.
"""

from __future__ import annotations

import logging
import urllib.parse
from typing import Any

from .. import database_utils as du
from ..config import get
from . import AFLOW

logger = logging.getLogger(__name__)

SOURCE_NAME = "AFLOW (Automatic-FLOW)"
SOURCE_URL = "https://aflow.org"


def aflux_endpoint() -> str:
    return get("network.endpoints.aflow_aflux", "https://aflow.org/API/aflux/")


def build_summon(fields: list[str], paging: tuple[int, int] | None = None) -> str:
    """Build an AFLUX summons string.

    ``fields`` are AFLUX properties (e.g. ``["compound", "auid"]``); optional
    ``paging`` is ``(page_index, page_size)``.
    """
    parts = [*fields]
    if paging is not None:
        parts.append(f"paging({paging[0]},{paging[1]})")
    return ",".join(parts)


def query(
    client: du.NetClient,
    fields: list[str] | None = None,
    paging: tuple[int, int] | None = None,
) -> Any:
    """Run one AFLUX query and return the parsed JSON payload."""
    fields = fields or ["compound", "auid"]
    summon = build_summon(fields, paging)
    url = f"{aflux_endpoint()}?{urllib.parse.quote(summon, safe=',()')}"
    return client.get_json(url)


def extract_records(payload: Any) -> list[dict[str, Any]]:
    """Parse an AFLUX JSON response into a list of record dicts.

    AFLUX returns either a list of strings (each tab/newline separated when
    multiple fields are requested) or a nested list.
    """
    records: list[dict[str, Any]] = []
    if payload is None:
        return records
    rows = payload if isinstance(payload, list) else [payload]
    for row in rows:
        if isinstance(row, str):
            # A single-field query may return one string per record.
            records.append({"value": row})
        elif isinstance(row, list):
            # Multiple requested fields -> one list element per record.
            records.append({"fields": row})
        elif isinstance(row, dict):
            records.append(row)
    return records


def compound_to_metadata(compound: str, aflow_id: str = "") -> dict[str, Any] | None:
    """Canonicalise an AFLOW ``compound`` string (e.g. ``Al2O3``) to metadata."""
    from ..formula_utils import FormulaError, get_anonymous_formula, get_elements, normalize_formula

    compound = (compound or "").strip()
    if not compound:
        return None
    try:
        reduced = normalize_formula(compound)
    except FormulaError:
        return None
    return {
        "database": AFLOW,
        "database_id": aflow_id or compound,
        "reduced_formula": reduced,
        "anonymous_formula": get_anonymous_formula(reduced),
        "elements": " ".join(get_elements(reduced)),
        "nelements": len(get_elements(reduced)),
        "status": du.STATUS_VALID,
    }


def collect_compounds(
    client: du.NetClient,
    *,
    max_pages: int = 1,
    page_size: int = 64,
    database_version: str | None = None,
) -> list[dict[str, Any]]:
    """Fetch composition-only metadata from AFLUX (no structures).

    Returns a list of metadata rows with ``database == 'AFLOW'``. Set
    ``max_pages`` to zero/None to skip (AFLOW is disabled by default).
    """
    rows: list[dict[str, Any]] = []
    if not max_pages:
        return rows
    for page_index in range(int(max_pages)):
        payload = query(client, fields=["compound", "auid"], paging=(page_index, int(page_size)))
        for rec in extract_records(payload):
            compound = None
            aflow_id = ""
            if isinstance(rec, dict) and "fields" in rec:
                fields = rec["fields"]
                compound = fields[0] if len(fields) > 0 else ""
                aflow_id = fields[1] if len(fields) > 1 else ""
            elif isinstance(rec, dict):
                compound = rec.get("compound", rec.get("value", ""))
                aflow_id = rec.get("auid", rec.get("value", ""))
            else:
                compound = rec.get("compound", "")
                aflow_id = rec.get("auid", "")
            meta = compound_to_metadata(compound, aflow_id)
            if meta is not None:
                rows.append(meta)
    logger.info("Collected %d AFLOW compound records (metadata only).", len(rows))
    return rows


def resolve_structure(aflow_id: str) -> Any:
    """Resolve an AFLOW structure - intentionally unimplemented (see module doc).

    Returns ``None``. Kept as a stable plugin point so a full structure source
    can be added later without changing the pipeline.
    """
    logger.warning(
        "AFLOW structure resolution for '%s' is not implemented. AFLOW is used as a "
        "formula/metadata source only.", aflow_id)
    return None
