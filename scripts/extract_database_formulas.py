#!/usr/bin/env python3
"""Extract unique reduced formulas from COD, MP, and OQMD datasets.

Reads the locally downloaded datasets and produces one CSV per database
containing unique, normalized (reduced) formula strings. These CSVs can
then be used by the generation scripts to filter out known formulas.

Usage:
    python scripts/extract_database_formulas.py              # extract all three
    python scripts/extract_database_formulas.py --db mp      # extract MP only
    python scripts/extract_database_formulas.py --db cod mp  # extract COD and MP

Requires the package to be installed:
    pip install -e .
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

# Dataset root directories (adjust if your datasets live elsewhere)
DATASET_ROOTS = {
    "cod": Path.home() / "Git" / "COD-Dataset",
    "mp": Path.home() / "Git" / "MP-Dataset",
    "oqmd": Path.home() / "Git" / "OQMD-Dataset",
}

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "datasets"


def normalize_formula_safe(formula: str) -> str | None:
    """Normalize a formula string, returning None on failure."""
    try:
        from formula_utils import normalize_formula
        return normalize_formula(formula)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Materials Project
# ---------------------------------------------------------------------------

def extract_mp_formulas() -> set[str]:
    """Extract unique reduced formulas from the MP summary CSV."""
    mp_csv = DATASET_ROOTS["mp"] / "data" / "mp_full" / "mp_summary.csv"
    if not mp_csv.exists():
        print(f"  ERROR: MP summary not found at {mp_csv}")
        return set()

    print(f"  Reading {mp_csv} ...")
    formulas = set()
    import json

    with open(mp_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            comp_reduced = row.get("composition_reduced", "")
            if not comp_reduced:
                continue
            try:
                comp = json.loads(comp_reduced)
                # Build formula string from composition dict
                parts = []
                for el in sorted(comp):
                    count = int(comp[el])
                    if count == 1:
                        parts.append(el)
                    else:
                        parts.append(f"{el}{count}")
                formula_str = "".join(parts)
                normalized = normalize_formula_safe(formula_str)
                if normalized:
                    formulas.add(normalized)
            except (json.JSONDecodeError, ValueError):
                continue

    print(f"  MP: {len(formulas)} unique reduced formulas")
    return formulas


# ---------------------------------------------------------------------------
# COD (Crystallography Open Database)
# ---------------------------------------------------------------------------

def extract_cod_formulas() -> set[str]:
    """Extract unique reduced formulas from COD CIF files."""
    cif_dir = DATASET_ROOTS["cod"] / "cif"
    if not cif_dir.exists():
        print(f"  ERROR: COD CIF directory not found at {cif_dir}")
        return set()

    print(f"  Scanning CIF files in {cif_dir} ...")

    # Use grep to extract all _chemical_formula_sum lines at once (much faster)
    import subprocess
    result = subprocess.run(
        ["grep", "-rh", "_chemical_formula_sum", str(cif_dir)],
        capture_output=True, text=True, timeout=600,
    )

    formulas = set()
    # Pattern: _chemical_formula_sum  'C5 H17 Al N2 O8 P2'
    pattern = re.compile(r"_chemical_formula_sum\s+['\"]?([^'\"]+)['\"]?")

    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        m = pattern.match(line)
        if not m:
            continue
        raw = m.group(1).strip()
        # COD format: space-separated element+count, e.g. "C5 H17 Al N2 O8 P2"
        # Convert to standard formula: remove spaces
        formula_str = raw.replace(" ", "")
        normalized = normalize_formula_safe(formula_str)
        if normalized:
            formulas.add(normalized)

    print(f"  COD: {len(formulas)} unique reduced formulas")
    return formulas


# ---------------------------------------------------------------------------
# OQMD (Open Quantum Materials Database)
# ---------------------------------------------------------------------------

def extract_oqmd_formulas() -> set[str]:
    """Extract unique reduced formulas from the OQMD SQL dump."""
    sql_dir = DATASET_ROOTS["oqmd"] / "qmdb__v1_8__022026.sql"
    # Find the actual SQL file
    sql_file = None
    if sql_dir.is_dir():
        for f in sql_dir.iterdir():
            if f.suffix == ".sql":
                sql_file = f
                break
    elif sql_dir.is_file():
        sql_file = sql_dir

    if not sql_file or not sql_file.exists():
        print(f"  ERROR: OQMD SQL dump not found")
        return set()

    print(f"  Reading {sql_file} ...")

    formulas = set()
    # OQMD packs many rows per INSERT line:
    #   INSERT INTO `compositions` VALUES ('A1 B1','AB',2,...),('C1 D2','CD2',3,...);
    # Match every ('formula', pattern on compositions lines
    entry_pattern = re.compile(r"\('([^']*)',")

    with open(sql_file, "r", errors="ignore") as f:
        for line in f:
            if "INSERT INTO `compositions`" not in line:
                continue
            for m in entry_pattern.finditer(line):
                raw_formula = m.group(1)
                # OQMD format: "A1 Ag1 B1 Nd1 Se3" (space-separated)
                # Remove spaces to get standard formula
                formula_str = raw_formula.replace(" ", "")
                if not formula_str or formula_str in ("", "_"):
                    continue
                normalized = normalize_formula_safe(formula_str)
                if normalized:
                    formulas.add(normalized)

    print(f"  OQMD: {len(formulas)} unique reduced formulas")
    return formulas


# ---------------------------------------------------------------------------
# Save to CSV
# ---------------------------------------------------------------------------

def save_formulas_csv(formulas: set[str], output_path: Path) -> None:
    """Save a set of formulas to a CSV file with a single 'formula' column."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sorted_formulas = sorted(formulas)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["formula"])
        for formula in sorted_formulas:
            writer.writerow([formula])
    print(f"  Saved {len(sorted_formulas)} formulas to {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

EXTRACTORS = {
    "mp": ("mp_formulas.csv", extract_mp_formulas),
    "cod": ("cod_formulas.csv", extract_cod_formulas),
    "oqmd": ("oqmd_formulas.csv", extract_oqmd_formulas),
}


def main():
    parser = argparse.ArgumentParser(
        description="Extract unique formulas from COD, MP, and OQMD datasets."
    )
    parser.add_argument(
        "--db", nargs="+", choices=["cod", "mp", "oqmd"], default=["cod", "mp", "oqmd"],
        help="Which databases to extract (default: all three)",
    )
    parser.add_argument(
        "--output-dir", type=str, default=str(OUTPUT_DIR),
        help=f"Output directory for CSV files (default: {OUTPUT_DIR})",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Extracting formulas from: {', '.join(args.db)}")
    print(f"Output directory: {out_dir}")
    print()

    for db_name in args.db:
        csv_name, extractor = EXTRACTORS[db_name]
        print(f"[{db_name.upper()}]")
        formulas = extractor()
        if formulas:
            save_formulas_csv(formulas, out_dir / csv_name)
        print()

    print("Done.")


if __name__ == "__main__":
    main()
