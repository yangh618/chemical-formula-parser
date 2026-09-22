#!/usr/bin/env python3
"""Merge the per-database formula CSVs into one de-duplicated formula CSV.

Reads the per-database files written by ``scripts/extract_database_formulas.py``
(``datasets/{cod,mp,oqmd}_formulas.csv``, one formula per row) and writes
``datasets/known_formulas.csv`` containing **every unique formula exactly once**::

    formula,databases,anonymous_formula,elements,num_elements
    Ac,cod;mp;oqmd,A,Ac,1
    AcAg2As2,cod;mp,ABC3,Ac Ag As,3
    NbAl2Cr,cod;mp;oqmd,ABC2,Al Cr Nb,3

The ``databases`` column records which sources contain the formula, so merging
loses no provenance. The derived ``anonymous_formula`` / ``elements`` columns make
the element-combination and anonymous-stoichiometry novelty levels meaningful
(disable them with ``--no-metadata`` for a faster formula-only build).
Formulas are stored in a canonical, notation-independent form (polyatomic groups
such as ``Ac(AgAs)2`` are flattened to ``AcAg2As2``), which is exactly the form
the SMACT candidate generator emits, so the ``--filter-known`` screen compares
formulas exactly.

Usage:
    python scripts/merge_known_formulas.py                 # cod + mp + oqmd
    python scripts/merge_known_formulas.py --db cod mp     # subset
    python scripts/merge_known_formulas.py -o datasets/known_formulas.csv
    python scripts/merge_known_formulas.py --datasets-dir /path/to/datasets
    python scripts/merge_known_formulas.py --no-metadata   # formula list only

Requires the package to be installed:
    pip install -e .
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from common import configure_logging  # noqa: E402
from config import get_project_root  # noqa: E402
from novelty_filter import merge_database_csv_files, merged_csv_path  # noqa: E402

ROOT = get_project_root()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Merge COD/MP/OQMD formula CSVs into one unique-formula CSV."
    )
    parser.add_argument(
        "--datasets-dir", default=None,
        help="Directory holding {cod,mp,oqmd}_formulas.csv "
             "(default: datasets/ at the project root)",
    )
    parser.add_argument(
        "--db", nargs="+", choices=["cod", "mp", "oqmd"], default=None,
        help="Which databases to merge (default: all available)",
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output CSV path (default: <datasets-dir>/known_formulas.csv)",
    )
    parser.add_argument(
        "--no-metadata", action="store_true",
        help="Skip the derived anonymous_formula/elements/num_elements columns "
             "(faster, but element-combination and anonymous-stoichiometry "
             "novelty can then no longer be assessed from the file)",
    )
    args = parser.parse_args(argv)

    logger = configure_logging("merge_known_formulas", log_dir=ROOT / "logs")

    datasets_dir = Path(args.datasets_dir) if args.datasets_dir else ROOT / "datasets"
    if not datasets_dir.is_dir():
        logger.error("Datasets directory not found: %s", datasets_dir)
        print()
        print("Run scripts/extract_database_formulas.py first, or pass "
              "--datasets-dir pointing at the folder holding "
              "{cod,mp,oqmd}_formulas.csv.")
        return 1

    print(f"Merging formulas from {datasets_dir} ...")
    if not args.no_metadata:
        print("  (deriving anonymous_formula/elements per unique formula; "
              "use --no-metadata to skip)")
    try:
        result = merge_database_csv_files(
            datasets_dir=datasets_dir,
            databases=args.db,
            output_path=args.output,
            with_metadata=not args.no_metadata,
        )
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1

    if not result.per_database:
        logger.error("No per-database formula CSVs found in %s", datasets_dir)
        return 1

    print()
    for db_label, count in result.per_database.items():
        print(f"  {db_label.upper():4s}: {count:>9,d} unique formulas")
    print(f"  Overlap removed : {result.overlap_removed:>9,d} duplicates "
          f"across databases")
    print(f"  Unique total    : {result.total_unique:>9,d} formulas "
          f"(written once each)")
    print()
    print(f"Saved -> {result.output_path}")
    print()
    print("Use it in any generation script with:")
    print("  --filter-known                       # picks up the merged CSV")
    print(f"  --filter-known --known-csv {result.output_path}")
    default_path = merged_csv_path(datasets_dir)
    if result.output_path != default_path:
        print()
        print(f"NOTE: --filter-known looks for {default_path} by default; "
              "pass --known-csv to use a different file.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
