#!/usr/bin/env python3
"""Screen SMACT candidates against the known-formula database for novelty.

Inputs: a candidates CSV (from generate_formulas.py) + a known-formulas CSV.
The known file may be either

* the merged, de-duplicated ``datasets/known_formulas.csv`` (columns
  ``formula, databases``), as built by ``scripts/merge_known_formulas.py``, or
* the per-database schema (``formula, database [, anonymous_formula, elements]``).

Output: ``data/formulas/novel_formulas.csv``.

Each output row carries the SMACT filter flags from the candidates file plus
membership/novelty flags. The most important field is ``novel_exact_formula``
which is true only when the canonical reduced formula is absent from *every*
searched database. The three novelty levels (sections 15-16) are reported
separately and never conflated.

Terminology is deliberately conservative: these are *database-unseen
candidate compositions*, not new materials or stable phases (section 29).

Examples
--------
python scripts/find_novel_formulas.py
python scripts/find_novel_formulas.py --known datasets/known_formulas.csv
python scripts/find_novel_formulas.py --novel-only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import database_utils as du  # noqa: E402
from common import configure_logging  # noqa: E402
from config import get_project_root  # noqa: E402
from novelty import (  # noqa: E402
    KnownIndex,
    annotate,
    canonicalize_formula,
    nearest_known_distance,
)
from novelty_filter import load_known_index_from_csvs, merged_csv_path  # noqa: E402

ROOT = get_project_root()

_OUTPUT_COLUMNS = [
    "formula", "anonymous_formula", "elements", "num_elements", "stoichiometry",
    "smact_valid", "charge_neutral", "electronegativity_valid",
    "found_in_mp", "found_in_cod", "found_in_oqmd", "found_in_jarvis", "found_in_aflow",
    "novel_exact_formula", "novel_element_combination", "novel_anonymous_formula",
    "nearest_known_formula", "min_composition_distance", "min_element_set_distance",
]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Novelty screening of SMACT candidates.")
    parser.add_argument("--candidates", default="data/formulas/candidates.csv")
    parser.add_argument("--known", default="data/formulas/known_formulas.csv",
                        help="known-formulas CSV: either the merged "
                             "datasets/known_formulas.csv (one unique formula per "
                             "row) or a per-database file "
                             "(formula, database [, anonymous_formula, elements])")
    parser.add_argument("--out", default="data/formulas/novel_formulas.csv")
    parser.add_argument("--novel-only", action="store_true",
                        help="only write rows where novel_exact_formula is True")
    parser.add_argument("--with-distance", action="store_true", help="add nearest-known distance columns")
    parser.add_argument("--limit", type=int, default=None, help="cap candidates (testing)")
    args = parser.parse_args(argv)

    logger = configure_logging("find_novel_formulas", log_dir=ROOT / "logs", log_filename="novelty_filter.log")

    candidates_path = ROOT / args.candidates
    known_path = ROOT / args.known
    if not candidates_path.exists():
        logger.error("Candidates file missing: %s (run generate_formulas.py first)", candidates_path)
        return 1

    if not known_path.exists() and args.known == parser.get_default("known"):
        # Fall back to the merged unique-formula CSV built by merge_known_formulas.py.
        merged = merged_csv_path()
        if merged.exists():
            logger.info("Default known file missing; using merged CSV %s", merged)
            known_path = merged

    candidates_df = du.load_metadata_csv(candidates_path)
    index: KnownIndex
    if known_path.exists():
        index = load_known_index_from_csvs(csv_path=known_path)
        logger.info("Loaded %d unique known formulas from %s.",
                    len(index.all_formulas), known_path)
    else:
        logger.warning("Known-formula file missing (%s) -> all candidates treated as unseen.",
                       known_path)
        index = KnownIndex()
        index.finalize()

    rows = []
    for _i, rec in candidates_df.iterrows():
        raw = rec.get("formula", "")
        try:
            c = canonicalize_formula(str(raw))
        except Exception:  # noqa: BLE001
            continue
        out = annotate(c, index)
        for flag in ("smact_valid", "charge_neutral", "electronegativity_valid"):
            out[flag] = bool(rec[flag]) if flag in candidates_df.columns and not _isna(rec[flag]) else None
        if args.with_distance:
            out.update(nearest_known_distance(c, index))
        rows.append(out)
        if args.limit is not None and len(rows) >= args.limit:
            break

    if args.novel_only:
        rows = [r for r in rows if r.get("novel_exact_formula")]
    logger.info("Rows after selection: %d (novel_exact=%d)",
                len(rows), sum(1 for r in rows if r.get("novel_exact_formula")))

    out_path = ROOT / args.out
    distance_cols = ["nearest_known_formula", "min_composition_distance", "min_element_set_distance"]
    cols = [c for c in _OUTPUT_COLUMNS if c not in distance_cols]
    if args.with_distance:
        cols = cols + distance_cols
    du.write_metadata_csv(rows, out_path, columns=cols)
    logger.info("Wrote novelty-screened formulas -> %s", out_path)
    return 0


def _isna(value) -> bool:
    import math

    if value is None:
        return True
    try:
        return bool(math.isnan(float(value)))
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    raise SystemExit(main())
