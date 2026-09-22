"""Shared command-line driver for the candidate-generation example scripts.

The three material-system example scripts differ only in which composition space
they load and what they call it. This module holds the common argument parsing,
generation, SMACT screening, novelty filtering, diversity selection and output
logic so the behaviour stays identical across systems.

Pipeline order: enumerate charge-balanced compositions from config.yaml ->
SMACT chemical-plausibility screen (``smact.screening.smact_validity``) ->
drop formulas already in COD/MP/OQMD -> farthest-point diversity selection ->
flatten and emit.

Output formulas are always written in *flat* form (no parenthesised groups):
``La4(SO)3`` is emitted as ``La4S3O3``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from composition import CompositionSpace
from diversity import diversity_report, select_diverse
from formula_utils import FormulaError, flatten_formula


def build_arg_parser(description: str, epilog: str | None = None) -> argparse.ArgumentParser:
    """Create the standard argument parser for a generation script."""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-n", "--num", type=int, default=None,
        help="Number of candidates to return (default: max_results from config)",
    )
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output file path (default: print to stdout)",
    )
    parser.add_argument(
        "--filter-known", action="store_true",
        help="Exclude formulas already present in the COD/MP/OQMD databases",
    )
    parser.add_argument(
        "--db", nargs="+", choices=["cod", "mp", "oqmd"], default=None,
        help="Databases to filter against (default: the merged "
             "datasets/known_formulas.csv; with this flag the individual "
             "per-database CSVs are used instead). Only used with --filter-known",
    )
    parser.add_argument(
        "--known-csv", type=str, default=None,
        help="Known-formulas CSV to filter against, one formula per row "
             "(e.g. datasets/known_formulas.csv built by "
             "scripts/merge_known_formulas.py). Only used with --filter-known",
    )
    parser.add_argument(
        "--no-diverse", action="store_true",
        help="Return candidates in enumeration order instead of selecting a "
             "maximally spread-out subset",
    )
    parser.add_argument(
        "--pool", type=int, default=None,
        help="Candidate pool size for diversity selection "
             "(default: 20x the requested number, at least 5000)",
    )
    parser.add_argument(
        "--per-combo", type=int, default=None,
        help="Override the config per_combo_limit: max formulas kept per "
             "element combination. Stratification is what lets the pool span "
             "many species rather than exhausting the first few.",
    )
    parser.add_argument(
        "--max-cations", type=int, default=None,
        help="Override the config cap on distinct cation species per formula",
    )
    parser.add_argument(
        "--max-anions", type=int, default=None,
        help="Override the config cap on distinct anion species per formula",
    )
    parser.add_argument(
        "--w-composition", type=float, default=1.0,
        help="Weight on stoichiometric distance in diversity selection (default: 1.0)",
    )
    parser.add_argument(
        "--w-element-set", type=float, default=1.0,
        help="Weight on element-identity (1-Jaccard) distance in diversity "
             "selection (default: 1.0)",
    )
    parser.add_argument(
        "--no-smact", action="store_true",
        help="Skip the SMACT chemical-plausibility screen (charge neutrality "
             "under an oxidation-state set plus the Pauling electronegativity "
             "test). Screening is on by default.",
    )
    parser.add_argument(
        "--smact-ox-set", default="icsd24",
        choices=["icsd24", "icsd16", "smact14", "pymatgen_sp", "wiki"],
        help="SMACT oxidation-state set used by the screen (default: icsd24)",
    )
    parser.add_argument(
        "--smact-no-en", action="store_true",
        help="Screen on charge neutrality only, without the Pauling "
             "electronegativity test",
    )
    return parser


def print_space_header(space: CompositionSpace, title: str) -> None:
    """Print the composition-space banner and its constraints."""
    constraints = space.constraints
    print("=" * 60)
    print(f"  {title}")
    print("=" * 60)
    print()
    print(f"  Cations ({len(space.cations)}): {', '.join(space.cations)}")
    print(f"  Anions  ({len(space.anions)}): {', '.join(space.anions)}")
    print()
    print("  Constraints:")
    print(f"    Required anions  : {constraints.get('required_anions') or 'none'}")
    print(f"    Max stoich       : {constraints.get('max_stoich', 'unlimited')}")
    print(f"    Max total atoms  : {constraints.get('max_total_atoms', 'unlimited')}")
    print(f"    Max cation species: {constraints.get('max_cation_species') or 'unlimited'}")
    print(f"    Per-combo limit  : {constraints.get('per_combo_limit') or 'unlimited'}")
    print()


def _generation_kwargs(args: argparse.Namespace) -> dict[str, Any]:
    """Collect user overrides for ``CompositionSpace.generate``."""
    kwargs: dict[str, Any] = {}
    if args.per_combo is not None:
        kwargs["per_combo_limit"] = args.per_combo
    if args.max_cations is not None:
        kwargs["max_cations"] = args.max_cations
    if args.max_anions is not None:
        kwargs["max_anions"] = args.max_anions
    return kwargs


def _flatten_all(formulas: list[str]) -> list[str]:
    """Flatten every formula and drop duplicates, preserving order."""
    flat: list[str] = []
    seen: set[str] = set()
    for f in formulas:
        try:
            text = flatten_formula(f)
        except FormulaError:
            continue
        if text not in seen:
            seen.add(text)
            flat.append(text)
    return flat


def _smact_validator(args: argparse.Namespace):
    """Return the SMACT screening predicate for *args*, or ``None`` if disabled."""
    if getattr(args, "no_smact", False):
        return None

    from functools import partial

    from smactlib.filters import smact_valid

    return partial(
        smact_valid,
        oxidation_set=getattr(args, "smact_ox_set", "icsd24"),
        use_en=not getattr(args, "smact_no_en", False),
        include_alloys=True,
    )


def _screen(pool: list[str], validator) -> list[str]:
    """Apply the SMACT screen to *pool*, printing pass/reject counts."""
    if validator is None or not pool:
        return pool
    kept = [f for f in pool if validator(f)]
    rejected = len(pool) - len(kept)
    pct = rejected / len(pool)
    print(f"  SMACT screen: kept {len(kept)}, rejected {rejected} ({pct:.1%})")
    return kept


def _load_known(args: argparse.Namespace) -> set[str]:
    """Load the known-formula set to filter against, reporting the source used."""
    from novelty_filter import load_known_formulas_set, merged_csv_path

    csv_path = Path(args.known_csv) if args.known_csv else None
    merged = merged_csv_path()

    if csv_path is not None:
        if args.db:
            print(f"  NOTE: --known-csv takes precedence over --db "
                  f"({' '.join(args.db)})")
        source = str(csv_path)
    elif args.db:
        source = "datasets/{" + ",".join(sorted(args.db)) + "}_formulas.csv"
    elif merged.exists():
        source = f"{merged} (unique formulas from COD+MP+OQMD)"
    else:
        source = "datasets/{cod,mp,oqmd}_formulas.csv (merged CSV not found)"

    print("Loading known formulas ...")
    try:
        known = load_known_formulas_set(csv_path=csv_path, databases=args.db)
    except FileNotFoundError as exc:
        raise SystemExit(
            f"ERROR: {exc}\n"
            "Build the merged CSV with: python scripts/merge_known_formulas.py"
        ) from exc

    print(f"  Loaded {len(known)} unique known formulas from {source}")
    print()
    return known


def run(space: CompositionSpace, args: argparse.Namespace, title: str) -> list[str]:
    """Execute the generation pipeline and return the final flat formulas."""
    print_space_header(space, title)

    num = args.num or space.constraints.get("max_results", 5000)
    gen_kwargs = _generation_kwargs(args)
    validator = _smact_validator(args)

    known: set[str] | None = None
    if args.filter_known:
        known = _load_known(args)

    if args.no_diverse:
        candidates = _generate_enumeration_order(space, num, known, gen_kwargs, validator)
    else:
        candidates = _generate_diverse(space, num, known, gen_kwargs, args, validator)

    flat = _flatten_all(candidates)

    _print_diversity_report(flat, args)
    _emit(flat, args)
    return flat


def _generate_enumeration_order(
    space: CompositionSpace,
    num: int,
    known: set[str] | None,
    gen_kwargs: dict[str, Any],
    validator=None,
) -> list[str]:
    """Candidates in raw enumeration order (clustered by element combination)."""
    if known is None:
        print(f"Generating up to {num} charge-balanced candidates ...")
        candidates = space.generate(max_results=num, **gen_kwargs)
        print(f"  Generated {len(candidates)} candidates")
        return _screen(candidates, validator)

    from novelty_filter import generate_novel_candidates

    print(f"Searching for {num} novel candidates (enumeration order) ...")
    result = generate_novel_candidates(
        space, num, known, verbose=False, validator=validator, **gen_kwargs)
    print(f"  Generated {result.total_generated} charge-balanced candidates")
    if result.rejected_by_validator:
        print(f"  SMACT-invalid     : {result.rejected_by_validator}")
    print(f"  Known in databases : {result.known_filtered}")
    print(f"  Novel (unseen)     : {result.novel_found} ({result.novel_fraction:.1%})")
    if not result.target_met:
        reason = ("composition space exhausted" if result.space_exhausted
                  else "generation limit reached")
        print(f"  NOTE: found {result.novel_found} of {result.target} requested ({reason})")
    return result.novel


def _generate_diverse(
    space: CompositionSpace,
    num: int,
    known: set[str] | None,
    gen_kwargs: dict[str, Any],
    args: argparse.Namespace,
    validator=None,
) -> list[str]:
    """Build a broad pool, then select a maximally spread-out subset."""
    pool_size = args.pool or max(num * 20, 5000)
    print(f"Building candidate pool (up to {pool_size}) ...")
    pool = space.generate(max_results=pool_size, **gen_kwargs)
    print(f"  Pool: {len(pool)} charge-balanced candidates")

    pool = _screen(pool, validator)

    if known is not None:
        before = len(pool)
        pool = [f for f in pool if f not in known]
        removed = before - len(pool)
        pct = (removed / before) if before else 0.0
        print(f"  Novel: {len(pool)}  (filtered {removed} known, {pct:.1%} of pool)")

    if not pool:
        print("  No novel candidates available.")
        return []

    print(f"Selecting up to {num} maximally diverse candidates ...")
    selected = select_diverse(
        pool,
        num,
        w_composition=args.w_composition,
        w_element_set=args.w_element_set,
    )
    print(f"  Selected {len(selected)}")
    return selected


def _print_diversity_report(flat: list[str], args: argparse.Namespace) -> None:
    """Print spread statistics for the final selection."""
    if not flat:
        return
    report = diversity_report(
        flat,
        w_composition=args.w_composition,
        w_element_set=args.w_element_set,
    )
    d = report.to_dict()
    print()
    print("  Diversity of final selection:")
    print(f"    Formulas              : {d['n_formulas']}")
    print(f"    Distinct element sets : {d['distinct_element_sets']}")
    print(f"    Elements covered      : {d['elements_covered']}")
    print(f"    Mean pairwise distance: {d['mean_pairwise_distance']}")
    print(f"    Min pairwise distance : {d['min_pairwise_distance']}")
    print(f"    Largest element-set group: "
          f"{d['most_common_element_set']} x{d['most_common_element_set_count']}")


def _emit(flat: list[str], args: argparse.Namespace) -> None:
    """Write to file or print to stdout."""
    if args.output:
        outpath = Path(args.output)
        outpath.write_text("\n".join(flat) + ("\n" if flat else ""))
        print()
        print(f"Saved {len(flat)} candidates to {outpath}")
    else:
        print()
        print(f"{len(flat)} candidates:")
        print("-" * 40)
        for i, f in enumerate(flat, 1):
            print(f"  {i:4d}. {f}")
