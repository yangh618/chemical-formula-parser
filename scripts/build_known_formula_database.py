#!/usr/bin/env python3
"""Build the union known-formula database from per-database metadata CSVs.

Reads ``data/metadata/{cod,oqmd,jarvis,aflow}.csv`` (whichever exist) plus an
optional Materials Project export (see ``materials_project.metadata_file`` in
config.yaml) and writes ``data/formulas/known_formulas.csv``.

Provenance is preserved: a formula present in multiple databases yields one
row per database (section 9). Rows are canonicalised reduced formulas.

Examples
--------
python scripts/build_known_formula_database.py
python scripts/build_known_formula_database.py --mp data/metadata/mp.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crystalfno_discovery import database_utils as du  # noqa: E402
from crystalfno_discovery.common import configure_logging  # noqa: E402
from crystalfno_discovery.config import get, get_project_root  # noqa: E402

ROOT = get_project_root()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build the union known-formula database.")
    parser.add_argument("--mp", default=None, help="optional Materials Project metadata CSV")
    parser.add_argument("--out", default=None, help="output path for known_formulas.csv")
    args = parser.parse_args(argv)

    logger = configure_logging("build_known_formula_database", log_dir=ROOT / "logs",
                               log_filename="novelty_filter.log")

    mp_file = args.mp or get("materials_project.metadata_file", None)
    out = ROOT / (args.out or "data/formulas/known_formulas.csv")
    metadata_dir = ROOT / "data" / "metadata"
    combined = du.build_known_formula_database(
        metadata_dir=metadata_dir,
        mp_metadata_file=mp_file,
        out_path=out,
    )
    logger.info("Known-formula database written: %d rows -> %s", len(combined), out)
    if not combined.empty:
        logger.info("Unique reduced formulas: %d", combined["formula"].nunique())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
