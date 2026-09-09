#!/usr/bin/env python3
"""Generate chemically plausible candidate compositions with SMACT.

Enumerates SMACT-valid candidates over the configured element universe and
writes ``data/formulas/candidates.csv`` with per-candidate SMACT filter flags.

Examples
--------
python scripts/generate_formulas.py --num-elements 3 --max-stoich 8
python scripts/generate_formulas.py --num-elements 2 3 4 --max-stoich 8
python scripts/generate_formulas.py --elements Al Cr Co Mo Nb Ta Ti V W Zr \\
    --num-elements 3 4 --max-stoich 8
python scripts/generate_formulas.py --no-electronegativity-filter
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import database_utils as du  # noqa: E402
from common import configure_logging  # noqa: E402
from config import get, get_project_root  # noqa: E402
from novelty.canonicalize import canonicalize_formula  # noqa: E402
from smactlib import filters as smact_filters  # noqa: E402
from smactlib.generator import generate_formulas  # noqa: E402

ROOT = get_project_root()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="SMACT candidate-composition generation.")
    parser.add_argument("--elements", nargs="+", default=None,
                        help="element list (space or comma separated)")
    parser.add_argument("--num-elements", type=int, nargs="+", default=None, help="combination sizes")
    parser.add_argument("--max-stoich", type=int, default=None)
    parser.add_argument("--oxidation-states", default=None,
                        help="SMACT oxidation-state set: smact14/icsd16/icsd24/pymatgen_sp/wiki")
    parser.add_argument("--no-electronegativity-filter", action="store_true",
                        help="disable the Pauling electronegativity test")
    parser.add_argument("--no-alloys", action="store_true", help="disable the all-metal alloy pathway")
    parser.add_argument("--out", default=None, help="output candidates CSV path")
    args = parser.parse_args(argv)

    logger = configure_logging("generate_formulas", log_dir=ROOT / "logs", log_filename="generate_formulas.log")

    if args.elements:
        tokens = [tok for chunk in args.elements for tok in chunk.split(",")]
        elements = [e.strip() for e in tokens if e.strip()]
    else:
        elements = [e for e in (get("elements", []) or []) if e and str(e).strip()]
    num_elements = args.num_elements or get("smact.num_elements", [2, 3, 4])
    max_stoich = args.max_stoich or get("smact.max_stoich", 8)
    ox_set = args.oxidation_states or get("smact.oxidation_states", "icsd24")
    include_alloys = get("smact.include_alloys", True) and not args.no_alloys
    use_en = get("smact.electronegativity", False) and not args.no_electronegativity_filter

    logger.info("Generating candidates: elements=%d n_elements=%s max_stoich=%d ox=%s en=%s alloys=%s",
                len(elements), num_elements, max_stoich, ox_set, use_en, include_alloys)

    formulas = generate_formulas(
        elements,
        num_elements=num_elements,
        max_stoich=max_stoich,
        oxidation_states_set=ox_set,
        electronegativity_filter=use_en,
        include_alloys=include_alloys,
    )

    rows = []
    for formula in formulas:
        try:
            c = canonicalize_formula(formula)
        except Exception:  # noqa: BLE001
            continue
        row = c.to_row()
        row["smact_valid"] = smact_filters.smact_valid(
            formula, oxidation_set=ox_set, use_en=use_en, include_alloys=include_alloys)
        row["charge_neutral"] = smact_filters.charge_neutral(
            formula, oxidation_set=ox_set, include_alloys=include_alloys)
        row["electronegativity_valid"] = smact_filters.electronegativity_valid(
            formula, oxidation_set=ox_set, include_alloys=include_alloys)
        row["electronegativity_filter_applied"] = use_en
        rows.append(row)

    out = ROOT / (args.out or "data/formulas/candidates.csv")
    du.write_metadata_csv(rows, out)
    logger.info("Wrote %d SMACT candidates -> %s", len(rows), out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
