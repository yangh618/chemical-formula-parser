#!/usr/bin/env python3
"""Download AFLOW composition metadata (formula source only - no structures).

AFLOW is disabled by default (config ``databases.aflow: false``). This script
only fetches composition strings via AFLUX and records them for the known-formula /
novelty screen. Structure resolution from AFLOW is intentionally not
implemented (see ``src/.../databases/aflow.py``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from crystalfno_discovery import database_utils as du  # noqa: E402
from crystalfno_discovery.common import configure_logging  # noqa: E402
from crystalfno_discovery.config import get, get_project_root  # noqa: E402
from crystalfno_discovery.databases import AFLOW, aflow  # noqa: E402

ROOT = get_project_root()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Fetch AFLOW composition metadata via AFLUX.")
    parser.add_argument("--max-pages", type=int, default=None,
                        help="number of AFLUX pages to fetch (default 1)")
    parser.add_argument("--page-size", type=int, default=64)
    args = parser.parse_args(argv)

    logger = configure_logging("download_aflow", log_dir=ROOT / "logs", log_filename="download_aflow.log")
    max_pages = args.max_pages if args.max_pages is not None else 1
    if max_pages <= 0:
        logger.info("AFLOW disabled (no pages requested); nothing to do.")
        return 0

    client = du.NetClient(
        timeout=get("network.timeout_seconds", 60),
        retries=get("network.retries", 5),
        backoff_factor=get("network.backoff_factor", 1.5),
        rate_limit_seconds=get("network.rate_limit_seconds", 0.35),
    )
    rows = aflow.collect_compounds(client, max_pages=max_pages, page_size=args.page_size)
    if rows:
        df = pd.DataFrame(rows)
        out = ROOT / "data" / "metadata" / "aflow.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        logger.info("Wrote %d AFLOW compound records -> %s", len(rows), out)
    else:
        logger.warning("AFLUX returned no parseable compounds.")

    from crystalfno_discovery.database_utils import add_provenance, get_retrieval_date
    add_provenance(AFLOW, source_name=aflow.SOURCE_NAME, source_url=aflow.SOURCE_URL,
                   retrieval_date=get_retrieval_date(), record_count=len(rows),
                   notes="Composition/metadata only; structure resolution not implemented.")
    client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
