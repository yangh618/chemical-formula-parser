#!/usr/bin/env python3
"""Check a list of proposed formulas against the live databases (no bulk download).

Reads candidate formulas from a text file (one per line) or a CSV with a
formula column, then asks each enabled database (COD / OQMD / JARVIS, plus
AFLOW if enabled in config.yaml) whether that reduced formula already exists.
Network failures on a database are reported as ``error`` for that row rather
than aborting the run, so a single flaky source never blocks the rest.

This answers "is this proposed formula already known somewhere?" It does not
need the large offline archives, and it never claims a material is novel in an
absolute sense - only "not found in the databases that were reachable".

Examples
--------
python scripts/search_known_formulas.py candidates.txt
python scripts/search_known_formulas.py candidates.csv --output results.csv
python scripts/search_known_formulas.py --formula-column reduced_formula in.csv
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# Keep the underlying HTTP client quiet: a CLI lookup should not spam retry logs.
logging.getLogger("database_utils").setLevel(logging.ERROR)

from database_utils import NetClient  # noqa: E402
from formula_utils import safe_normalize_formula  # noqa: E402
from config import get  # noqa: E402

# provider -> endpoint config key
_STRUCTURE_PROVIDERS = {
    "cod": "network.endpoints.cod_optimade",
    "oqmd": "network.endpoints.oqmd_optimade",
    "jarvis": "network.endpoints.jarvis_optimade",
}

_FORMULA_COLUMNS = {"formula", "reduced_formula", "candidate", "composition", "candidate_formula"}


def enabled_databases() -> list[str]:
    """Return the list of live databases currently switched on in config.yaml."""
    db_config = get("databases", {}) or {}
    return [db for db in _STRUCTURE_PROVIDERS if db_config.get(db, False)]


def read_formulas(path: Path, formula_column: str | None = None) -> list[str]:
    """Read candidate formulas from a plain-text file or a CSV."""
    if path.suffix.lower() == ".csv":
        return _read_csv_formulas(path, formula_column)
    out = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            out.append(line)
    return out


def _read_csv_formulas(path: Path, formula_column: str | None) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        column = formula_column
        if column is None:
            column = next((c for c in reader.fieldnames if c.strip().lower() in _FORMULA_COLUMNS), None)
        if column is None:
            column = reader.fieldnames[0]
        return [(row.get(column) or "").strip() for row in reader if (row.get(column) or "").strip()]


def _provider_endpoint(provider: str) -> str:
    return get(_STRUCTURE_PROVIDERS[provider], "")


def formula_in_provider(client: NetClient, provider: str, reduced: str) -> bool:
    """True if *any* structure with this reduced formula exists in the provider."""
    endpoint = _provider_endpoint(provider)
    params = {"filter": f'chemical_formula_reduced = "{reduced}"', "page_limit": 5}
    payload = client.get_json(endpoint, params=params)
    rows = payload.get("data", []) or []
    return len(rows) > 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Check proposed formulas against live databases for existing entries.")
    parser.add_argument("input", help="path to a text file (one formula per line) or CSV with a formula column")
    parser.add_argument("--output", default=None, help="optional CSV path to write per-formula results")
    parser.add_argument("--formula-column", default=None, help="CSV column holding the formula (auto-detected if omitted)")
    parser.add_argument("--databases", default=None,
                        help="comma-separated list to query (default: those enabled in config.yaml)")
    args = parser.parse_args(argv)

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"error: no such file: {input_path}", file=sys.stderr)
        return 2

    formulas = read_formulas(input_path, args.formula_column)
    if not formulas:
        print("no formulas found in input", file=sys.stderr)
        return 2

    providers = [p for p in (args.databases or "").split(",") if p.strip()] if args.databases \
        else enabled_databases()
    if not providers:
        print("no databases enabled; enable cod/oqmd/jarvis in config.yaml or pass --databases",
              file=sys.stderr)
        return 2

    print(f"checking {len(formulas)} formulas against: {', '.join(providers)}")

    client = NetClient()
    rows = []
    known_in_any = 0
    errors = 0
    try:
        for raw in formulas:
            reduced = safe_normalize_formula(raw)
            row = {"input_formula": raw, "canonical_formula": reduced or raw}
            row["normalised"] = "yes" if reduced else "no"
            if reduced is None:
                for provider in providers:
                    row[provider] = "invalid"
                rows.append(row)
                print(f"{raw:30s} invalid formula", flush=True)
                continue
            found = []
            for provider in providers:
                try:
                    is_known = formula_in_provider(client, provider, reduced)
                except Exception as exc:  # noqa: BLE001 - never let one DB abort the batch
                    row[provider] = "error"
                    errors += 1
                    print(f"  {provider}: query error -> {type(exc).__name__}: {str(exc)[:70]}", flush=True)
                    continue
                row[provider] = "yes" if is_known else "no"
                if is_known:
                    found.append(provider)
            if found:
                known_in_any += 1
            status = "KNOWN in " + ",".join(found) if found else "not found (in reachable DBs)"
            rows.append(row)
            print(f"{raw:30s} {status}", flush=True)
    finally:
        client.close()

    if args.output:
        fieldnames = ["input_formula", "canonical_formula", "normalised", *providers]
        with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote results to {args.output}")

    print(f"\nsummary: {len(formulas)} checked | {known_in_any} found in at least one DB | "
          f"{errors} database query errors")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
