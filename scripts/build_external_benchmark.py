#!/usr/bin/env python3
"""Build formula-unseen external test sets for a frozen CrystalFNO model.

Reads the unified per-database metadata CSVs (cod/oqmd/jarvis, plus aflow if
present) and produces, per database, a CSV of **external-test structures**
whose canonical reduced formula is absent from the Materials Project training
set (when ``--exclude-mp``) and that pass structural-quality filters.

This does NOT train or retrain CrystalFNO - it only assembles validation input
for an already-trained, frozen model (sections 12, 29).

Examples
--------
python scripts/build_external_benchmark.py --exclude-mp --exclude-disordered \\
    --min-elements 2 --max-elements 6
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from crystalfno_discovery import database_utils as du  # noqa: E402
from crystalfno_discovery.common import configure_logging  # noqa: E402
from crystalfno_discovery.config import get_project_root  # noqa: E402
from crystalfno_discovery.formula_utils import FormulaError, get_num_elements, normalize_formula  # noqa: E402

ROOT = get_project_root()

_SKIP_STATUSES = {"parse_error", "missing_coordinates", "invalid_cif", "invalid_lattice"}


def _canonical(frame, i: int) -> str | None:
    rec = frame.iloc[i]
    for col in ("reduced_formula", "formula"):
        if col in frame.columns and rec[col] is not None:
            try:
                return normalize_formula(str(rec[col]))
            except FormulaError:
                continue
    return None


def _load_mp_formulas(known_path: Path, mp_path: Path | None) -> set[str]:
    """Return canonical MP formulas from the known DB (or a raw MP file)."""
    known = du.load_metadata_csv(known_path)
    mp: set[str] = set()
    if not known.empty and "database" in known.columns and "formula" in known.columns:
        sub = known[known["database"].astype(str).str.lower() == "mp"]
        for f in sub["formula"]:
            try:
                mp.add(normalize_formula(str(f)))
            except FormulaError:
                continue
    if mp_path is not None and mp_path.exists():
        frame = du.load_metadata_csv(mp_path)
        if not frame.empty:
            formula_col = next((c for c in ("reduced_formula", "formula", "pretty_formula", "composition")
                                if c in frame.columns), None)
            if formula_col:
                for f in frame[formula_col]:
                    try:
                        mp.add(normalize_formula(str(f)))
                    except FormulaError:
                        continue
    return mp


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build formula-unseen external benchmark sets.")
    parser.add_argument("--metadata-dir", default="data/metadata")
    parser.add_argument("--known", default="data/formulas/known_formulas.csv")
    parser.add_argument("--mp", default=None, help="optional raw Materials Project CSV")
    parser.add_argument("--dbs", default="cod oqmd jarvis", help="databases to benchmark")
    parser.add_argument("--out", default="data/benchmark")
    parser.add_argument("--exclude-mp", action="store_true", help="keep only formulas unseen in MP")
    parser.add_argument("--exclude-disordered", action="store_true")
    parser.add_argument("--min-elements", type=int, default=2)
    parser.add_argument("--max-elements", type=int, default=6)
    parser.add_argument("--max-atoms", type=int, default=200)
    args = parser.parse_args(argv)

    logger = configure_logging("build_external_benchmark", log_dir=ROOT / "logs",
                               log_filename="download_oqmd.log")

    metadata_dir = ROOT / args.metadata_dir
    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    mp_path = Path(args.mp) if args.mp else None
    mp_formulas = _load_mp_formulas(ROOT / args.known, mp_path) if args.exclude_mp else set()
    logger.info("MP formula set size (exclude_mp=%s): %d", args.exclude_mp, len(mp_formulas))

    summary_rows = []
    for db in args.dbs.split():
        path = metadata_dir / f"{db.lower()}.csv"
        if not path.exists():
            logger.warning("No metadata for %s at %s; skipping.", db.upper(), path)
            continue
        frame = du.load_metadata_csv(path)
        if frame.empty:
            logger.warning("%s metadata empty; skipping.", db.upper())
            continue

        has_status = "status" in frame.columns
        kept = []
        dropped = {"disordered": 0, "parse": 0, "size": 0, "mp_overlap": 0}
        for i in range(len(frame)):
            status = str(frame.iloc[i].get("status", "")) if has_status else ""
            if status in _SKIP_STATUSES or (args.exclude_disordered and status in
                                            ("disordered", "partial_occupancy")):
                dropped["disordered" if status in ("disordered", "partial_occupancy") else "parse"] += 1
                continue
            if args.exclude_disordered and not has_status:
                rec = frame.iloc[i]
                if bool(rec.get("disordered", False)) or bool(rec.get("partial_occupancy", False)):
                    dropped["disordered"] += 1
                    continue
            formula = _canonical(frame, i)
            if formula is None:
                dropped["parse"] += 1
                continue
            if not (args.min_elements <= get_num_elements(formula) <= args.max_elements):
                dropped["size"] += 1
                continue
            if "nsites" in frame.columns and _isna(frame.iloc[i]["nsites"]) is False:
                if float(frame.iloc[i]["nsites"]) > args.max_atoms:
                    dropped["size"] += 1
                    continue
            if args.exclude_mp and formula in mp_formulas:
                dropped["mp_overlap"] += 1
                continue
            row = frame.iloc[i].to_dict()
            row["reduced_formula"] = formula
            kept.append(row)

        out = du.write_metadata_csv(kept, out_dir / f"{db.lower()}_external_test.csv", columns=None)
        n_formulas = len({r["reduced_formula"] for r in kept})
        summary_rows.append({
            "database": db.upper(), "structures": len(kept), "unique_formulas": n_formulas,
            "dropped_disordered": dropped["disordered"], "dropped_parse": dropped["parse"],
            "dropped_size": dropped["size"], "dropped_mp_overlap": dropped["mp_overlap"],
        })
        logger.info("%s external test: %d structures, %d unique formulas -> %s",
                    db.upper(), len(kept), n_formulas, out)

    du.write_metadata_csv(summary_rows, out_dir / "summary.csv")
    logger.info("Benchmark summary written to %s", out_dir / "summary.csv")
    return 0


def _isna(value):
    if value is None:
        return True
    try:
        return bool(float(value) != float(value))
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    raise SystemExit(main())
