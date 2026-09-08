#!/usr/bin/env python3
"""Download COD CIF structures (resumable) and refresh the COD metadata CSV.

Official COD access: per-record CIFs at ``crystallography.net/cod/<id>.cif``
plus the COD OPTIMADE endpoint for id/metadata enumeration.

Examples
--------
# parse a handful already available locally under a COD archive dir
python scripts/download_cod.py --archive data/raw/cod

# download 20 CIFs enumerated from the COD OPTIMADE endpoint
python scripts/download_cod.py --limit 20

# quick connectivity/parser test
python scripts/download_cod.py --test

# full extraction (very large; resumes where it left off)
python scripts/download_cod.py --all
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# --- make the package importable without installation ----------------------
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crystalfno_discovery import database_utils as du  # noqa: E402
from crystalfno_discovery.common import configure_logging  # noqa: E402
from crystalfno_discovery.config import get, get_project_root  # noqa: E402
from crystalfno_discovery.databases import COD, cod  # noqa: E402

ROOT = get_project_root()


def _client() -> du.NetClient:
    return du.NetClient(
        timeout=get("network.timeout_seconds", 60),
        retries=get("network.retries", 5),
        backoff_factor=get("network.backoff_factor", 1.5),
        rate_limit_seconds=get("network.rate_limit_seconds", 0.35),
    )


def _ids_from_args(client: du.NetClient, args: argparse.Namespace):
    """Return an iterator of COD ids to process."""
    from crystalfno_discovery.config import get

    if args.archive:
        ids = cod.read_ids_from_archive(ROOT / args.archive)
        logger.info("Using %d COD ids from local archive %s", len(ids), args.archive)
        return iter(ids)
    if args.cod_ids:
        return iter([x.strip() for x in args.cod_ids.split(",") if x.strip()])
    if args.ids_file:
        path = ROOT / args.ids_file
        ids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        return iter(ids)

    # Otherwise enumerate ids from the COD OPTIMADE endpoint.
    max_pages = None if args.all else (args.max_pages or (1 if args.test else 1))
    limit = args.limit or (5 if args.test else (get("network.page_limit", 100) if not args.all else None))
    enumerated = []
    for idx, (cod_id, _attrib) in enumerate(cod.iter_records(client, max_pages=max_pages)):
        if limit is not None and idx >= limit:
            break
        enumerated.append(cod_id)
    if not enumerated and not args.all:
        logger.warning("No COD ids resolved. Use --cod-ids, --ids-file, --archive, --limit or --all.")
    return iter(enumerated)


def main(argv=None) -> int:
    global logger
    parser = argparse.ArgumentParser(description="Download COD CIF structures (resumable).")
    parser.add_argument("--test", action="store_true", help="quick connectivity/parser test (5 records)")
    parser.add_argument("--all", action="store_true", help="download the entire enumerated COD collection")
    parser.add_argument("--max-pages", type=int, default=None, help="max OPTIMADE metadata pages to enumerate")
    parser.add_argument("--limit", type=int, default=None, help="max number of records to download")
    parser.add_argument("--cod-ids", default=None, help="comma-separated COD ids")
    parser.add_argument("--ids-file", default=None, help="text file of COD ids, one per line")
    parser.add_argument("--archive", default=None, help="dir of local COD CIFs to use instead of the network")
    args = parser.parse_args(argv)

    logger = configure_logging("download_cod", log_dir=ROOT / "logs", log_filename="download_cod.log")

    client = _client()
    dirs = cod.ensure_dirs()
    dirs["raw"].mkdir(parents=True, exist_ok=True)
    dirs["structures"].mkdir(parents=True, exist_ok=True)

    counts = {}
    downloaded = 0
    for cod_id in _ids_from_args(client, args):
        status, is_new, _meta, reason = cod.download_one(
            client, cod_id, dirs["raw"] / f"{cod_id}.cif", dirs["structures"] / f"{cod_id}.cif")
        counts[status] = counts.get(status, 0) + 1
        if is_new:
            downloaded += 1
        if reason and status in (du.STATUS_PARSE_ERROR, du.STATUS_INVALID_LATTICE):
            logger.warning("COD %s -> %s: %s", cod_id, status, reason)

    logger.info("COD download summary: %s (newly downloaded=%d)", counts, downloaded)

    # Rebuild the COD metadata CSV from every structure we now hold (idempotent).
    rows = du.scan_structure_directory(dirs["structures"], database=COD)
    du.write_metadata_csv(rows, dirs["metadata"])
    logger.info("Wrote COD metadata with %d rows -> %s", len(rows), dirs["metadata"])

    from crystalfno_discovery.database_utils import add_provenance, get_retrieval_date
    add_provenance(COD, source_name=cod.SOURCE_NAME, source_url=cod.SOURCE_URL,
                   retrieval_date=get_retrieval_date(), record_count=len(rows),
                   notes="Structures under data/structures/cod; raw CIF under data/raw/cod.")
    client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
