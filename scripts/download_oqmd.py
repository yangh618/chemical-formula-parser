#!/usr/bin/env python3
"""Download OQMD structures through its OPTIMADE endpoint (resumable).

Examples
--------
python scripts/download_oqmd.py --max-pages 10     # cheap sample
python scripts/download_oqmd.py --all              # complete extraction
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crystalfno_discovery import database_utils as du  # noqa: E402
from crystalfno_discovery.common import configure_logging  # noqa: E402
from crystalfno_discovery.config import get, get_project_root  # noqa: E402
from crystalfno_discovery.databases import OQMD, oqmd  # noqa: E402

ROOT = get_project_root()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Download OQMD structures via OPTIMADE.")
    parser.add_argument("--max-pages", type=int, default=None, help="max pages to fetch (testing)")
    parser.add_argument("--all", action="store_true", help="download the full OQMD collection")
    parser.add_argument("--page-limit", type=int, default=None)
    args = parser.parse_args(argv)

    logger = configure_logging("download_oqmd", log_dir=ROOT / "logs", log_filename="download_oqmd.log")

    client = du.NetClient(
        timeout=get("network.timeout_seconds", 60),
        retries=get("network.retries", 5),
        backoff_factor=get("network.backoff_factor", 1.5),
        rate_limit_seconds=get("network.rate_limit_seconds", 0.35),
    )

    max_pages = None if args.all else (args.max_pages or 5)
    summary = oqmd.collect(client, max_pages=max_pages, page_limit=args.page_limit)
    logger.info("OQMD summary: downloaded=%d skipped=%d failed=%d statuses=%s",
                summary["downloaded"], summary["skipped"], summary["failed"], summary["status_counts"])

    # Rebuild the OQMD metadata CSV from held structures (idempotent + resumable).
    struct_dir = ROOT / "data" / "structures" / "oqmd"
    rows = du.scan_structure_directory(struct_dir, database=OQMD)
    du.write_metadata_csv(rows, ROOT / "data" / "metadata" / "oqmd.csv")
    logger.info("Wrote OQMD metadata with %d rows.", len(rows))

    from crystalfno_discovery.database_utils import add_provenance, get_retrieval_date
    add_provenance(OQMD, source_name=oqmd.SOURCE_NAME, source_url=oqmd.SOURCE_URL,
                   retrieval_date=get_retrieval_date(), record_count=len(rows))
    client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
