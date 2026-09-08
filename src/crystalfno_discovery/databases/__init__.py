"""Database clients (COD, OQMD, JARVIS, AFLOW).

Each client implements the current, publicly accessible interface:

* **OQMD**   - OPTIMADE v1 ``https://oqmd.org/optimade/v1/structures``
* **JARVIS** - OPTIMADE v1 ``https://jarvis.nist.gov/optimade/jarvisdft/v1/structures``
* **COD**    - per-id CIF ``https://www.crystallography.net/cod/<id>.cif`` plus
               OPTIMADE metadata ``.../cod/optimade/v1/structures``
* **AFLOW**  - AFLUX ``https://aflow.org/API/aflux/`` (formula/metadata only;
               structure coordinates are not resolved and remain documented)

Database records are converted to a unified pymatgen structure representation
and written as CIF + CSV metadata (see :mod:`crystalfno_discovery.database_utils`).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .. import database_utils as du
from ..config import get_project_root
from ..structure_utils import structure_to_metadata

# Short human names used as provenance keys / column values.
COD = "COD"
OQMD = "OQMD"
JARVIS = "JARVIS"
AFLOW = "AFLOW"


def data_root() -> Path:
    """Absolute path to the project ``data/`` directory."""
    return get_project_root() / "data"


def metadata_csv(db: str) -> Path:
    return data_root() / "metadata" / f"{db.lower()}.csv"


def raw_dir(db: str) -> Path:
    return data_root() / "raw" / db.lower()


def structures_dir(db: str) -> Path:
    return data_root() / "structures" / db.lower()


def cif_path(db: str, db_id: str) -> Path:
    return structures_dir(db) / f"{db_id}.cif"


def read_processed_ids(path: Path) -> set[str]:
    """Read previously processed structure IDs (one per line) for resumability."""
    if not path.exists():
        return set()
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def append_processed_id(path: Path, db_id: str) -> None:
    """Append a processed ID; creating the file as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{db_id}\n")


def record_to_metadata(
    attributes: dict[str, Any],
    database: str,
    database_id: str,
    structure_file: str = "",
) -> tuple[Any, dict[str, Any], str]:
    """Build ``(structure, metadata, status)`` from OPTIMADE structure attributes.

    Uses ``structure_from_optimade``; on success a pymatgen structure is
    returned together with unified metadata. Records that lack coordinates
    (e.g. COD OPTIMADE metadata-only rows) yield ``(None, None, status)``.
    """
    try:
        structure = du.structure_from_optimade(attributes)
    except ValueError:
        return None, None, du.STATUS_MISSING_COORDINATES
    status = du.STATUS_VALID
    metadata = structure_to_metadata(structure, database=database, database_id=database_id,
                                     structure_file=structure_file)
    partial, disordered = du.site_issues(structure)
    metadata["partial_occupancy"] = partial
    metadata["disordered"] = disordered
    if disordered:
        status = du.STATUS_DISORDERED
    elif partial:
        status = du.STATUS_PARTIAL_OCCUPANCY
    metadata["status"] = status
    return structure, metadata, status


def save_optimade_structure(
    attributes: dict[str, Any],
    database: str,
    database_id: str,
    *,
    struct_dir: str | Path | None = None,
    raw_dir_path: str | Path | None = None,
    save_raw: bool = True,
) -> tuple[bool, dict[str, Any] | None, str]:
    """Persist one OPTIMADE structure as CIF (+ optional raw JSON).

    Returns ``(is_new, metadata, status)``. Skips (is_new=False) when the CIF
    already exists, which is what makes interrupted downloads resumable.
    """
    import json

    from ..structure_utils import save_cif

    struct_dir = Path(struct_dir) if struct_dir else structures_dir(database)
    struct_dir.mkdir(parents=True, exist_ok=True)
    cif = struct_dir / f"{database_id}.cif"
    if cif.exists() and cif.stat().st_size > 0:
        _s, meta, _st, _r = du.parse_structure_file(cif, database=database, database_id=database_id)
        return False, meta, du.STATUS_DUPLICATE

    structure, metadata, status = record_to_metadata(attributes, database, database_id)
    if structure is None:
        return False, metadata, status  # no coordinates -> nothing to save

    try:
        save_cif(structure, str(cif))
    except Exception as exc:  # noqa: BLE001
        return False, metadata, du.STATUS_PARSE_ERROR
    metadata["structure_file"] = str(cif)

    if save_raw and raw_dir_path is not None:
        raw = Path(raw_dir_path)
        raw.mkdir(parents=True, exist_ok=True)
        try:
            (raw / f"{database_id}.json").write_text(json.dumps(attributes, indent=2), encoding="utf-8")
        except OSError:
            logger = __import__("logging").getLogger(__name__)
            logger.warning("Could not write raw JSON for %s/%s", database, database_id)
    return True, metadata, status


def collect_optimade(
    client: du.NetClient,
    database: str,
    base_url: str,
    *,
    max_pages: int | None = None,
    page_limit: int | None = None,
    filter_expr: str | None = None,
    save_raw: bool = True,
) -> dict[str, Any]:
    """Enumerate an OPTIMADE ``/structures`` endpoint and save every structure.

    Resumable: records whose CIF already exists are skipped. Returns a summary
    dict including a ``rows`` list of unified metadata records.
    """
    import logging

    from ..config import get

    logger = logging.getLogger(__name__)
    page_limit = int(page_limit or get("network.page_limit", 100))
    summary = {
        "downloaded": 0,
        "skipped": 0,
        "failed": 0,
        "rows": [],
        "status_counts": {},
    }
    for rows, meta in du.optimade_pages(
        client, base_url, page_limit=page_limit, max_pages=max_pages, filter_expr=filter_expr
    ):
        for record in rows:
            attributes = record.get("attributes", {}) or {}
            db_id = str(record.get("id", ""))
            if not db_id:
                summary["failed"] += 1
                continue
            is_new, metadata, status = save_optimade_structure(
                attributes, database, db_id,
                struct_dir=structures_dir(database),
                raw_dir_path=raw_dir(database),
                save_raw=save_raw,
            )
            summary["status_counts"][status] = summary["status_counts"].get(status, 0) + 1
            if is_new:
                summary["downloaded"] += 1
            elif status == du.STATUS_DUPLICATE:
                summary["skipped"] += 1
            else:
                summary["failed"] += 1
            if metadata is not None:
                metadata["database"] = database
                metadata["database_id"] = db_id
                summary["rows"].append(metadata)
    return summary

