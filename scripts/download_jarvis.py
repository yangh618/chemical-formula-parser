#!/usr/bin/env python3
"""Download JARVIS-DFT structures via its OPTIMADE endpoint (resumable).

The JARVIS OPTIMADE server only returns records for an explicit ``elements``
filter, so enumeration is driven by element combinations over the configured
element universe (default: the ``elements`` list in config.yaml).

Examples
--------
python scripts/download_jarvis.py --max-pages 10
python scripts/download_jarvis.py --all
python scripts/download_jarvis.py --all --elements Li O Fe
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crystalfno_discovery import database_utils as du  # noqa: E402
from crystalfno_discovery.common import configure_logging  # noqa: E402
from crystalfno_discovery.config import get, get_project_root  # noqa: E402
from crystalfno_discovery.databases import JARVIS, jarvis  # noqa: E402

ROOT = get_project_root()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Download JARVIS-DFT structures via OPTIMADE.")
    parser.add_argument("--max-pages", type=int, default=None)
    parser.add_argument("--all", action="store_true", help="full element-combination sweep")
    parser.add_argument("--page-limit", type=int, default=None)
    parser.add_argument("--elements", default=None, help="comma-separated element universe")
    parser.add_argument("--num-elements", default="1",
                        help="element-combination sizes to enumerate, e.g. '1' or '1 2'")
    parser.add_argument("--max-filters", type=int, default=None, help="cap on filter queries (testing)")
    args = parser.parse_args(argv)

    logger = configure_logging("download_jarvis", log_dir=ROOT / "logs", log_filename="download_jarvis.log")

    elements = [e.strip() for e in (args.elements or get("elements", [])) if e.strip()]
    sizes = tuple(int(n) for n in args.num_elements.split())
    if not elements:
        logger.error("No element universe given; set config 'elements' or pass --elements.")
        return 2

    # Cheap, safe default sample if neither --all nor --max-pages is given.
    max_pages = None if args.all else (args.max_pages or 5)

    client = du.NetClient(
        timeout=get("network.timeout_seconds", 60),
        retries=get("network.retries", 5),
        backoff_factor=get("network.backoff_factor", 1.5),
        rate_limit_seconds=get("network.rate_limit_seconds", 0.35),
    )
    summary = jarvis.collect(
        client, elements, sizes=sizes, max_pages=max_pages,
        page_limit=args.page_limit, max_filters=args.max_filters)
    logger.info("JARVIS summary: downloaded=%d skipped=%d failed=%d filters=%d statuses=%s",
                summary["downloaded"], summary["skipped"], summary["failed"],
                summary["filters_used"], summary["status_counts"])

    struct_dir = ROOT / "data" / "structures" / "jarvis"
    rows = du.scan_structure_directory(struct_dir, database=JARVIS)
    du.write_metadata_csv(rows, ROOT / "data" / "metadata" / "jarvis.csv")
    logger.info("Wrote JARVIS metadata with %d rows.", len(rows))

    from crystalfno_discovery.database_utils import add_provenance, get_retrieval_date
    add_provenance(JARVIS, source_name=jarvis.SOURCE_NAME, source_url=jarvis.SOURCE_URL,
                   retrieval_date=get_retrieval_date(), record_count=len(rows))
    client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
