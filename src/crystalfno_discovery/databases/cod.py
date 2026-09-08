"""COD (Crystallography Open Database) client.

Official access is used (no scraping of individual HTML pages):

* structure files: ``https://www.crystallography.net/cod/<id>.cif`` (the flat
  per-record CIF distribution that COD maintains and mirrors);
* record/metadata discovery: the COD OPTIMADE endpoint
  ``https://www.crystallography.net/cod/optimade/v1/structures``
  (formula + lattice + database ID; no atomic coordinates).

Full-archive bulk dumps (e.g. ``cod.zip`` / ``cod-cifs-mysql.txz`` on
``https://www.crystallography.net/cod/archives/``) are tens of GB and are best
downloaded out-of-band; point at a local archive with ``--archive`` so records
are read from disk instead of the network.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from .. import database_utils as du
from ..config import get
from . import COD, metadata_csv, raw_dir, structures_dir

logger = logging.getLogger(__name__)

SOURCE_NAME = "Crystallography Open Database (COD)"
SOURCE_URL = "https://www.crystallography.net/cod/"


def optimade_endpoint() -> str:
    return get("network.endpoints.cod_optimade", "https://www.crystallography.net/cod/optimade/v1/structures")


def cif_url(cod_id: str | int) -> str:
    """Per-record CIF URL for a COD id."""
    template = get("network.endpoints.cod_cif", "https://www.crystallography.net/cod/{cod_id}.cif")
    return template.format(cod_id=int(cod_id))


def iter_records(
    client: du.NetClient,
    max_pages: int | None = None,
    page_limit: int | None = None,
    filter_expr: str | None = None,
):
    """Yield ``(cod_id, attributes)`` from the COD OPTIMADE metadata endpoint."""
    page_limit = int(page_limit or get("network.page_limit", 100))
    for rows, _meta in du.optimade_pages(
        client, optimade_endpoint(), page_limit=page_limit, max_pages=max_pages, filter_expr=filter_expr
    ):
        for record in rows:
            attribs = record.get("attributes", {}) or {}
            yield str(record.get("id", "")), attribs


def metadata_from_attributes(attributes: dict[str, Any], cod_id: str) -> dict[str, Any]:
    """Build a formula-level metadata row from COD OPTIMADE attributes.

    COD OPTIMADE rows carry formula + lattice metadata but **no atomic
    coordinates**, so this is a metadata-only record (structure_file blank).
    """
    from ..formula_utils import FormulaError, get_anonymous_formula, get_elements, normalize_formula

    reduced_raw = attributes.get("chemical_formula_reduced") or attributes.get("_cod_formula") or ""
    try:
        reduced = normalize_formula(str(reduced_raw))
    except FormulaError:
        reduced = str(reduced_raw).strip()
    elements = get_elements(reduced) if reduced else []
    meta: dict[str, Any] = {
        "database": COD,
        "database_id": cod_id,
        "reduced_formula": reduced,
        "anonymous_formula": get_anonymous_formula(reduced) if reduced else "",
        "elements": " ".join(elements),
        "nelements": len(elements),
    }
    lv = attributes.get("lattice_vectors")
    if lv:
        import numpy as np

        from pymatgen.core import Lattice

        try:
            lat = Lattice(np.asarray(lv, dtype=float))
            a, b, c, al, be, ga = lat.parameters
            meta.update(lattice_a=a, lattice_b=b, lattice_c=c,
                        alpha=al, beta=be, gamma=ga, volume=float(lat.volume))
        except Exception:  # noqa: BLE001
            pass
    elif attributes.get("_cod_a"):
        meta.update(
            lattice_a=attributes.get("_cod_a"), lattice_b=attributes.get("_cod_b"),
            lattice_c=attributes.get("_cod_c"), alpha=attributes.get("_cod_alpha"),
            beta=attributes.get("_cod_beta"), gamma=attributes.get("_cod_gamma"),
            volume=attributes.get("_cod_volume"),
        )
    sg = attributes.get("_cod_sg") or attributes.get("space_group_symbol_hermann_mauguin")
    if sg:
        meta["space_group"] = sg
    meta["status"] = du.STATUS_VALID
    return meta


def download_one(
    client: du.NetClient,
    cod_id: str,
    raw_out: str | Path,
    struct_out: str | Path,
) -> tuple[str, bool, dict[str, Any] | None, str]:
    """Download + normalize one COD record (resumable).

    Returns ``(status, downloaded_new, metadata, reason)``. If the normalized
    structure file already exists the download is skipped so a restarted run
    continues rather than restarting (sections 5, 23).
    """
    from ..structure_utils import save_cif

    raw_out = Path(raw_out)
    struct_out = Path(struct_out)
    raw_out.parent.mkdir(parents=True, exist_ok=True)
    struct_out.parent.mkdir(parents=True, exist_ok=True)

    if struct_out.exists() and struct_out.stat().st_size > 0:
        _s, meta, _st, _r = du.parse_structure_file(struct_out, database=COD, database_id=cod_id)
        return du.STATUS_DUPLICATE, False, meta, "already processed (skipped)"

    downloaded_new = client.download_to_file(cif_url(cod_id), raw_out)
    if not downloaded_new and not raw_out.exists():
        return du.STATUS_PARSE_ERROR, False, None, "download produced no file"

    try:
        structure, meta, status, reason = du.parse_structure_file(raw_out, database=COD, database_id=cod_id)
    except Exception as exc:  # noqa: BLE001
        return du.STATUS_PARSE_ERROR, downloaded_new, None, f"parse failure: {exc}"

    if structure is not None:
        try:
            save_cif(structure, str(struct_out))
        except Exception as exc:  # noqa: BLE001
            return du.STATUS_PARSE_ERROR, downloaded_new, None, f"could not write normalised CIF: {exc}"
        if meta is not None:
            meta["structure_file"] = str(struct_out)
    return status, downloaded_new, meta, reason


def read_ids_from_archive(archive_dir: str | Path) -> list[str]:
    """List COD ids available in a local archive of downloaded CIF files."""
    archive_dir = Path(archive_dir)
    if not archive_dir.exists():
        return []
    return sorted({p.stem for p in archive_dir.glob("*.cif")})


def ensure_dirs() -> dict[str, Path]:
    """Return the ``raw``, ``structures`` and ``metadata`` paths for COD."""
    return {"raw": raw_dir(COD), "structures": structures_dir(COD), "metadata": metadata_csv(COD)}
