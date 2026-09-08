#!/usr/bin/env python3
"""Rebuild a database metadata CSV by parsing every held structure CIF.

Structures are stored as ``data/structures/<db>/<id>.cif`` (see section 3).
This script scans a chosen directory and (re)writes the corresponding
``data/metadata/<db>.csv`` with per-structure status flags (valid /
partial_occupancy / disordered / parse_error / invalid_lattice / ...).
Nothing is deleted; problematic structures are tagged, not removed.

Examples
--------
python scripts/parse_structures.py --db cod
python scripts/parse_structures.py --db all --limit 1000
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crystalfno_discovery import database_utils as du  # noqa: E402
from crystalfno_discovery.common import configure_logging  # noqa: E402
from crystalfno_discovery.config import get_project_root  # noqa: E402
from crystalfno_discovery.databases import AFLOW, COD, JARVIS, OQMD  # noqa: E402

ROOT = get_project_root()

_DBS = {"cod": COD, "oqmd": OQMD, "jarvis": JARVIS, "aflow": AFLOW}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Parse structure CIFs into metadata CSVs.")
    parser.add_argument("--db", default="all", choices=["all", "cod", "oqmd", "jarvis", "aflow"])
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args(argv)

    logger = configure_logging("parse_structures", log_dir=ROOT / "logs", log_filename="download_cod.log")

    dbs = list(_DBS) if args.db == "all" else [args.db]
    for db in dbs:
        if db == "aflow":
            logger.info("AFLOW stores formula metadata (no structure CIFs); nothing to parse.")
            continue
        struct_dir = ROOT / "data" / "structures" / db
        rows = du.scan_structure_directory(struct_dir, database=_DBS[db], limit=args.limit)
        out = ROOT / "data" / "metadata" / f"{db}.csv"
        du.write_metadata_csv(rows, out)
        from collections import Counter

        counts = Counter(r.get("status", "") for r in rows)
        logger.info("%s: parsed %d structures -> %s  statuses=%s", db.upper(), len(rows), out, dict(counts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
