"""Shared database utilities: HTTP client, OPTIMADE/Structure conversion,
CIF parsing with status classification, known-formula database assembly and
data provenance tracking.

All network access flows through :class:`NetClient`, which applies retries,
timeouts, back-off and (optionally) inter-request rate limiting so that a
multi-million-entry download does not die on a transient network failure and
does not hammer a public server (sections 6, 7, 23).
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any

import numpy as np
import requests
from pymatgen.core import Lattice, Structure

from .common import RateLimiter
from .formula_utils import FormulaError, get_anonymous_formula, get_elements, normalize_formula
from .structure_utils import (
    METADATA_COLUMNS,
    site_issues,
    structure_to_metadata,
    validate_structure,
)

logger = logging.getLogger(__name__)

# Structural statuses used to tag parsed files (never automatically deleted).
STATUS_VALID = "valid"
STATUS_PARTIAL_OCCUPANCY = "partial_occupancy"
STATUS_DISORDERED = "disordered"
STATUS_MISSING_COORDINATES = "missing_coordinates"
STATUS_INVALID_LATTICE = "invalid_lattice"
STATUS_PARSE_ERROR = "parse_error"
STATUS_DUPLICATE = "duplicate"

ALL_STATUSES = {
    STATUS_VALID,
    STATUS_PARTIAL_OCCUPANCY,
    STATUS_DISORDERED,
    STATUS_MISSING_COORDINATES,
    STATUS_INVALID_LATTICE,
    STATUS_PARSE_ERROR,
    STATUS_DUPLICATE,
}


class NetClient:
    """HTTP GET client with retries, back-off, timeouts and rate limiting."""

    def __init__(
        self,
        timeout: float = 60,
        retries: int = 5,
        backoff_factor: float = 1.5,
        rate_limit_seconds: float = 0.0,
        headers: dict[str, str] | None = None,
    ) -> None:
        self.timeout = float(timeout)
        self.retries = int(retries)
        self.backoff_factor = float(backoff_factor)
        self.rate_limiter = RateLimiter(rate_limit_seconds)
        self.session = requests.Session()
        self.session.headers.update(
            headers or {"User-Agent": "crystalfno-discovery/0.1 (+research; contact: local)"}
        )

    def _request(self, method: str, url: str, **kwargs: Any) -> requests.Response:
        last_exc: Exception | None = None
        for attempt in range(1, self.retries + 1):
            self.rate_limiter.wait()
            try:
                resp = self.session.request(method, url, timeout=self.timeout, **kwargs)
                if resp.status_code in (429, 500, 502, 503, 504):
                    # Transient server errors -> retry with back-off.
                    raise requests.HTTPError(f"HTTP {resp.status_code} for {url}", response=resp)
                resp.raise_for_status()
                return resp
            except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as exc:
                last_exc = exc
                wait = self.backoff_factor ** (attempt - 1)
                logger.warning("Request failed (%s) attempt %d/%d for %s; retrying in %.1fs",
                               type(exc).__name__, attempt, self.retries, url, wait)
                if attempt < self.retries:
                    time.sleep(wait)
        raise requests.RequestException(f"Giving up after {self.retries} attempts: {url}") from last_exc

    def get_json(self, url: str, params: dict[str, Any] | None = None, **kwargs: Any) -> Any:
        resp = self._request("GET", url, params=params, **kwargs)
        return resp.json()

    def get_text(self, url: str, params: dict[str, Any] | None = None, **kwargs: Any) -> str:
        resp = self._request("GET", url, params=params, **kwargs)
        return resp.text

    def download_to_file(self, url: str, out_path: str | Path, params: dict[str, Any] | None = None) -> bool:
        """Download *url* to *out_path*; returns False if the file already exists."""
        out_path = Path(out_path)
        if out_path.exists() and out_path.stat().st_size > 0:
            return False
        out_path.parent.mkdir(parents=True, exist_ok=True)
        resp = self._request("GET", url, params=params, stream=True)
        with out_path.open("wb") as handle:
            for chunk in resp.iter_content(chunk_size=1 << 16):
                handle.write(chunk)
        return True

    def close(self) -> None:
        self.session.close()


def structure_from_optimade(attributes: dict[str, Any]) -> Structure:
    """Build a pymatgen :class:`Structure` from OPTIMADE structure attributes.

    Requires ``lattice_vectors``, ``cartesian_site_positions`` and
    ``species_at_sites`` (or ``species``). Raises :class:`ValueError` when the
    required structural fields are absent (e.g. COD's metadata-only records).
    """
    lattice_vectors = attributes.get("lattice_vectors")
    positions = attributes.get("cartesian_site_positions")
    species_at_sites = attributes.get("species_at_sites")
    if lattice_vectors is None or positions is None or species_at_sites is None:
        raise ValueError("OPTIMADE record lacks structural coordinates")
    matrix = np.asarray(lattice_vectors, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError(f"lattice_vectors must be 3x3, got {matrix.shape}")
    coords = np.asarray(positions, dtype=float)
    species = [spec.split(" ")[0] if isinstance(spec, str) else spec for spec in species_at_sites]
    if len(species) != len(coords):
        raise ValueError("species_at_sites/cartesian_site_positions length mismatch")
    return Structure(Lattice(matrix), species, coords, coords_are_cartesian=True)


_ATOM_SITE_MARKERS = (
    "_atom_site_fract_x",
    "_atom_site_fract_y",
    "_atom_site_fract_z",
    "_atom_site_Cartn_x",
    "_atom_site_cartn_x",
    "_atom_site_frac_x",
)


def parse_structure_file(
    path: str | Path,
    database: str = "",
    database_id: str = "",
) -> tuple[Structure | None, dict[str, Any] | None, str, str]:
    """Parse one CIF into ``(structure, metadata, status, reason)``.

    Problematic structures are *never deleted*; instead a ``status`` is
    reported so downstream tools can filter them (section 5). Statuses belong
    to :data:`ALL_STATUSES`.
    """
    path = Path(path)
    reason = ""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        return None, None, STATUS_PARSE_ERROR, f"unreadable file: {exc}"

    try:
        structure = Structure.from_file(str(path))
    except Exception as exc:  # noqa: BLE001
        if not any(marker in text.lower() for marker in _ATOM_SITE_MARKERS):
            return None, None, STATUS_MISSING_COORDINATES, f"no atom coordinates found ({exc})"
        return None, None, STATUS_PARSE_ERROR, f"pymatgen parse error: {exc}"

    partial, disordered = site_issues(structure)
    val = validate_structure(structure)
    status = STATUS_VALID
    if disordered:
        status = STATUS_DISORDERED
    elif partial:
        status = STATUS_PARTIAL_OCCUPANCY
    if not val["valid"] and val["reason"] and str(val["reason"]).startswith("lattice"):
        status = STATUS_INVALID_LATTICE
        reason = val["reason"]
    elif status == STATUS_VALID and not val["valid"]:
        reason = val["reason"] or ""

    metadata = structure_to_metadata(structure, database=database, database_id=database_id,
                                     structure_file=str(path))
    metadata["status"] = status
    metadata["reason"] = reason
    metadata["partial_occupancy"] = partial
    metadata["disordered"] = disordered
    return structure, metadata, status, reason


def scan_structure_directory(
    structure_dir: str | Path,
    database: str,
    pattern: str = "*.cif",
    limit: int | None = None,
) -> list[dict[str, Any]]:
    """Parse every CIF under *structure_dir*, returning metadata rows.

    The ``status`` column reflects :data:`ALL_STATUSES`; nothing is deleted.
    """
    structure_dir = Path(structure_dir)
    rows: list[dict[str, Any]] = []
    files = sorted(structure_dir.glob(pattern)) if structure_dir.exists() else []
    for idx, cif in enumerate(files):
        if limit is not None and idx >= limit:
            break
        database_id = cif.stem
        _struct, metadata, _status, _reason = parse_structure_file(
            cif, database=database, database_id=database_id)
        if metadata is not None:
            rows.append(metadata)
        else:
            rows.append({
                "database": database, "database_id": database_id,
                "status": _status, "reason": _reason, "structure_file": str(cif),
            })
    return rows


def write_metadata_csv(rows: list[dict[str, Any]], out_path: str | Path,
                       columns: list[str] | None = None) -> Path:
    """Write *rows* to *out_path* as CSV, preserving database-specific columns."""
    import pandas as pd

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        pd.DataFrame(columns=columns or METADATA_COLUMNS).to_csv(out_path, index=False)
        return out_path
    df = pd.DataFrame(rows)
    if columns is not None:
        for col in columns:
            if col not in df.columns:
                df[col] = ""
        order = columns + [c for c in df.columns if c not in columns]
        df = df[order]
    df.to_csv(out_path, index=False)
    return out_path


def load_metadata_csv(path: str | Path) -> "pd.DataFrame":
    """Load a metadata CSV if present, otherwise return an empty DataFrame."""
    import pandas as pd

    path = Path(path)
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _resolve_formula_column(frame: "pd.DataFrame") -> str:
    """Pick the most reliable formula column present in *frame*."""
    for col in ("reduced_formula", "reduced composition", "formula", "composition", "pretty_formula", "name"):
        if col in frame.columns:
            return col
    raise ValueError(f"No formula column in dataframe with columns {list(frame.columns)}")


def _normalise_db_frame(frame: "pd.DataFrame", database: str) -> "pd.DataFrame":
    """Normalise a per-database metadata frame into known-formula rows."""
    import pandas as pd

    cols = ["formula", "database", "database_id", "anonymous_formula", "elements"]
    if frame.empty:
        return pd.DataFrame(columns=cols)
    formula_col = _resolve_formula_column(frame)
    if "database_id" not in frame.columns:
        raise ValueError(f"'{database}' metadata lacks a 'database_id' column")

    id_col = "material_id" if ("material_id" in frame.columns and database.lower() == "mp") else "database_id"
    rows: list[dict[str, Any]] = []
    for _i, record in frame.iterrows():
        raw = record[formula_col]
        if raw is None or (isinstance(raw, float) and np.isnan(raw)):
            continue
        raw = str(raw).strip()
        try:
            canonical = normalize_formula(raw)
        except FormulaError:
            continue
        db_id = record[id_col]
        rows.append({
            "formula": canonical,
            "database": database,
            "database_id": db_id,
            "anonymous_formula": get_anonymous_formula(canonical),
            "elements": " ".join(get_elements(canonical)),
        })
    return pd.DataFrame(rows, columns=cols)


def build_known_formula_database(
    metadata_dir: str | Path,
    sources: dict[str, str] | None = None,
    mp_metadata_file: str | Path | None = None,
    out_path: str | Path | None = None,
) -> "pd.DataFrame":
    """Assemble ``known_formulas.csv`` from per-database metadata CSVs.

    ``metadata_dir`` is ``data/metadata``; ``sources`` maps a database name to
    its CSV filename (defaults to ``<db>.csv``). ``mp_metadata_file`` is an
    optional Materials Project export (its formula column is auto-detected).

    Rows are *not* deduplicated across databases: a formula present in both
    COD and OQMD yields two provenance-preserving rows (section 9).
    """
    import pandas as pd

    metadata_dir = Path(metadata_dir)
    sources = sources or {"cod": "cod.csv", "oqmd": "oqmd.csv", "jarvis": "jarvis.csv", "aflow": "aflow.csv"}

    frames: list["pd.DataFrame"] = []
    for db, filename in sources.items():
        path = metadata_dir / filename
        if path.exists():
            frame = load_metadata_csv(path)
            norm = _normalise_db_frame(frame, database=db)
            if not norm.empty:
                frames.append(norm)
            logger.info("Known formulas from %s: %d", db.upper(), len(norm))

    if mp_metadata_file is not None and Path(mp_metadata_file).exists():
        mp_frame = load_metadata_csv(mp_metadata_file)
        if not mp_frame.empty:
            norm = _normalise_db_frame(mp_frame, database="MP")
            frames.append(norm)
            logger.info("Known formulas from Materials Project: %d", len(norm))
    else:
        logger.info("Materials Project metadata absent -> no MP overlap removal applied.")

    if frames:
        combined = pd.concat(frames, ignore_index=True)
    else:
        combined = pd.DataFrame(columns=["formula", "database", "database_id", "anonymous_formula", "elements"])
    if out_path is not None:
        write_metadata_csv(combined.to_dict("records"), out_path,
                           columns=["formula", "database", "database_id", "anonymous_formula", "elements"])
    return combined


# ---------------------------------------------------------------------------
# Provenance (section 28)
# ---------------------------------------------------------------------------

def default_provenance_path() -> Path:
    """Return ``data/metadata/provenance.json`` under the project root."""
    from .config import get_project_root

    return get_project_root() / "data" / "metadata" / "provenance.json"


def read_provenance(path: str | Path | None = None) -> dict[str, Any]:
    """Read the provenance JSON (empty dict if missing)."""
    path = Path(path) if path else default_provenance_path()
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    return {}


def add_provenance(
    database: str,
    *,
    source_name: str,
    source_url: str,
    retrieval_date: str,
    database_version: str | None = None,
    record_count: int | None = None,
    notes: str = "",
    path: str | Path | None = None,
) -> None:
    """Record data-provenance information for *database*."""
    data = read_provenance(path)
    entry = {
        "database": database,
        "source_name": source_name,
        "source_url": source_url,
        "retrieval_date": retrieval_date,
        "database_version": database_version,
        "record_count": record_count,
        "notes": notes,
    }
    data.setdefault("databases", {})[database] = entry
    out_path = Path(path) if path else default_provenance_path()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)


def get_retrieval_date() -> str:
    """ISO-8601 UTC timestamp used to mark a retrieval."""
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def optimade_pages(
    client: NetClient,
    base_url: str,
    *,
    page_limit: int = 100,
    max_pages: int | None = None,
    page_offset_start: int = 0,
    filter_expr: str | None = None,
) -> Any:
    """Yield ``(rows, meta)`` pages from an OPTIMADE ``/structures`` endpoint.

    Paginates with the standard ``page_limit``/``page_offset`` parameters and
    stops when the server reports ``more_data_available == false`` or after
    ``max_pages`` (both are honoured). Yields lazily so a ``--max-pages`` CLI
    flag can be used to test a download cheaply (sections 6, 7).
    """
    offset = int(page_offset_start)
    pages_done = 0
    while max_pages is None or pages_done < max_pages:
        params: dict[str, Any] = {"page_limit": int(page_limit), "page_offset": int(offset)}
        if filter_expr:
            params["filter"] = filter_expr
        payload = client.get_json(base_url, params=params)
        meta = payload.get("meta", {}) or {}
        rows = payload.get("data", []) or []
        yield rows, meta
        pages_done += 1
        if not meta.get("more_data_available", False):
            break
        if not rows:
            logger.warning("Server reports more data but returned an empty page at offset=%d; stopping.", offset)
            break
        offset += int(page_limit)

