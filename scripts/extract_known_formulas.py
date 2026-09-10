#!/usr/bin/env python3
"""Custom extractors: build a known-formula database from local datasets.

No database server required. Two subcommands stream the raw data directly:

* ``cod``  - walk a COD CIF mirror (the rsync ``cif/<d>/<dd>/<dd>/<id>.cif``
  layout) and read each entry's ``_chemical_formula_sum`` header field.
* ``oqmd`` - stream the OQMD ``qmdb`` MySQL dump (plain ``.sql`` or ``.gz``) and
  pull the ``formula`` column from the ``compositions`` table's INSERT rows,
  without importing the ~100 GB database into a server.

Both canonicalise to reduced formulas and write an intermediate CSV with the
columns the novelty screen expects (``formula,database,database_id,
anonymous_formula,elements``). Use ``merge`` to combine them into
``data/formulas/known_formulas.csv``.

Examples
--------
python scripts/extract_known_formulas.py cod
python scripts/extract_known_formulas.py cod --cod-dir datasets/COD-Dataset/cif --limit 20000
python scripts/extract_known_formulas.py oqmd
python scripts/extract_known_formulas.py merge
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from formula_utils import get_anonymous_formula, get_elements, safe_normalize_formula  # noqa: E402
from config import get_project_root  # noqa: E402

ROOT = get_project_root()
COLUMNS = ["formula", "database", "database_id", "anonymous_formula", "elements"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
# The shared HTTP client logger is irrelevant here.
logging.getLogger("database_utils").setLevel(logging.ERROR)
logger = logging.getLogger("extract_known_formulas")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def canon_formula(raw: str) -> str | None:
    """Canonicalise a dataset formula string to a reduced formula (or None)."""
    if not raw:
        return None
    token = raw.strip().replace(" ", "").replace("_", "").replace("\t", "")
    if not token:
        return None
    try:
        return safe_normalize_formula(token)
    except Exception:  # noqa: BLE001 - exotic/partial records are skipped
        return None


def row_for(canon: str, database: str, database_id: str) -> dict[str, str]:
    """Build one known-formula row for a canonical formula."""
    return {
        "formula": canon,
        "database": database,
        "database_id": database_id,
        "anonymous_formula": get_anonymous_formula(canon),
        "elements": " ".join(get_elements(canon)),
    }


def write_rows(rows: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows -> %s", len(rows), out_path)


# ---------------------------------------------------------------------------
# COD: walk a CIF mirror, read the _chemical_formula_sum header field
# ---------------------------------------------------------------------------

def read_cif_formula(path: Path, header_bytes: int = 16384) -> str | None:
    """Return the raw ``_chemical_formula_sum`` (fallback ``_moiety``) value."""
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            head = handle.read(header_bytes)
    except OSError:
        return None
    moiety = None
    for line in head.splitlines():
        if line.startswith("_chemical_formula_sum"):
            value = line[len("_chemical_formula_sum"):].strip()
            return _strip_quotes(value)
        if line.startswith("_chemical_formula_moiety") and moiety is None:
            moiety = _strip_quotes(line[len("_chemical_formula_moiety"):].strip())
    return moiety


def _strip_quotes(value: str) -> str | None:
    value = value.strip()
    if not value:
        return None
    if value[0] in "'\"":
        quote = value[0]
        value = value[1:].split(quote, 1)[0]
    return value.strip() or None


def cmd_cod(args: argparse.Namespace) -> int:
    cif_dir = Path(args.cod_dir)
    if not cif_dir.exists():
        logger.error("COD CIF directory not found: %s", cif_dir)
        return 2

    raw_to_id: dict[str, str] = {}
    scanned = 0
    with_formula = 0
    for path in cif_dir.rglob("*.cif"):
        scanned += 1
        if args.limit and scanned > args.limit:
            break
        raw = read_cif_formula(path, header_bytes=args.header_bytes)
        if raw:
            with_formula += 1
            raw_to_id.setdefault(raw, path.stem)
        if args.log_every and scanned % args.log_every == 0:
            logger.info("COD: scanned %d files (%d with a formula, %d unique)",
                        scanned, with_formula, len(raw_to_id))

    logger.info("COD: scanned %d files; %d had a formula; %d unique raw strings",
                scanned, with_formula, len(raw_to_id))

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    invalid = 0
    for raw, db_id in raw_to_id.items():
        canon = canon_formula(raw)
        if canon is None:
            invalid += 1
            continue
        if canon in seen:
            continue
        seen.add(canon)
        rows.append(row_for(canon, "COD", db_id))
    logger.info("COD: %d unique canonical formulas (%d unparseable raw strings)", len(rows), invalid)
    write_rows(rows, Path(args.output))
    return 0


# ---------------------------------------------------------------------------
# OQMD: stream the qmdb SQL dump and pull compositions.formula values
# ---------------------------------------------------------------------------

def iter_first_fields(statement: str) -> list[str]:
    """Return the first column of every ``(...)`` tuple in an INSERT statement.

    The ``compositions`` table's first column is ``formula`` (e.g. ``Fe2 O3``),
    so this extracts the primary key directly without a SQL server.
    """
    values: list[str] = []
    idx = statement.find("VALUES")
    i = idx + 6 if idx != -1 else 0
    n = len(statement)
    while i < n:
        while i < n and statement[i] != "(":
            if statement[i] == ";":
                return values
            i += 1
        if i >= n:
            break
        i += 1  # skip '('
        while i < n and statement[i] == " ":
            i += 1
        if i < n and statement[i] == "'":
            i += 1
            chars: list[str] = []
            while i < n:
                ch = statement[i]
                if ch == "\\" and i + 1 < n:
                    chars.append(statement[i + 1])
                    i += 2
                    continue
                if ch == "'":
                    if i + 1 < n and statement[i + 1] == "'":
                        chars.append("'")
                        i += 2
                        continue
                    i += 1
                    break
                chars.append(ch)
                i += 1
            values.append("".join(chars))
        # skip the rest of this tuple (ignoring quoted commas/parens)
        depth = 1
        while i < n:
            ch = statement[i]
            if ch == "'":
                i += 1
                while i < n:
                    if statement[i] == "\\":
                        i += 2
                        continue
                    if statement[i] == "'":
                        if i + 1 < n and statement[i + 1] == "'":
                            i += 2
                            continue
                        i += 1
                        break
                    i += 1
                continue
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    i += 1
                    break
            i += 1
    return values


def iter_oqmd_formulas(dump_path: Path):
    """Yield raw ``compositions.formula`` strings from an OQMD SQL dump."""
    opener = gzip.open if str(dump_path).endswith(".gz") else open
    in_compositions = False
    with opener(dump_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not in_compositions:
                if line.startswith("-- Dumping data for table `compositions`"):
                    in_compositions = True
                continue
            if line.startswith("INSERT INTO `compositions`"):
                yield from iter_first_fields(line)
            elif line.startswith("-- ") or line.startswith("UNLOCK TABLES"):
                # Reached the next table's section -> compositions data is done.
                return


def cmd_oqmd(args: argparse.Namespace) -> int:
    dump = Path(args.dump)
    if not dump.exists():
        logger.error("OQMD dump not found: %s", dump)
        return 2

    raw_set: set[str] = set()
    count = 0
    for raw in iter_oqmd_formulas(dump):
        count += 1
        raw_set.add(raw)
        if args.limit and count >= args.limit:
            break
        if args.log_every and count % args.log_every == 0:
            logger.info("OQMD: read %d rows (%d unique)", count, len(raw_set))

    logger.info("OQMD: read %d composition rows; %d unique raw strings", count, len(raw_set))

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    invalid = 0
    for raw in raw_set:
        canon = canon_formula(raw)
        if canon is None:
            invalid += 1
            continue
        if canon in seen:
            continue
        seen.add(canon)
        rows.append(row_for(canon, "OQMD", raw.strip()))
    logger.info("OQMD: %d unique canonical formulas (%d unparseable)", len(rows), invalid)
    write_rows(rows, Path(args.output))
    return 0


# ---------------------------------------------------------------------------
# mp: query the Materials Project API (requires an API key + mp-api)
# ---------------------------------------------------------------------------

def cmd_mp(args: argparse.Namespace) -> int:
    api_key = args.api_key or os.environ.get("MP_API_KEY") or os.environ.get("PMG_MAPI_KEY")
    if not api_key:
        logger.error(
            "No Materials Project API key found. Set MP_API_KEY (or pass --api-key). "
            "Get a free key at https://materialsproject.org/api")
        return 2
    try:
        from mp_api.client import MPRester
    except ImportError:
        logger.error("The Materials Project client is not installed: pip install mp-api")
        return 2

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    seen_ids: set[str] = set()
    count = 0
    invalid = 0
    logger.info("Querying Materials Project summary docs (this can take a while)...")
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(fields=["material_id", "formula_pretty"])
        for doc in docs:
            count += 1
            if args.limit and count > args.limit:
                break
            mid = str(getattr(doc, "material_id", "") or "")
            formula = getattr(doc, "formula_pretty", None)
            if mid in seen_ids:
                continue
            seen_ids.add(mid)
            canon = canon_formula(formula) if formula else None
            if canon is None:
                invalid += 1
                continue
            if canon in seen:
                continue
            seen.add(canon)
            rows.append(row_for(canon, "MP", mid))
            if args.log_every and count % args.log_every == 0:
                logger.info("MP: %d docs (%d unique formulas)", count, len(rows))

    logger.info("MP: %d docs -> %d unique formulas (%d unparseable)", count, len(rows), invalid)
    write_rows(rows, Path(args.output))
    return 0


# ---------------------------------------------------------------------------
# merge: combine per-database extractor outputs into known_formulas.csv
# ---------------------------------------------------------------------------

def cmd_merge(args: argparse.Namespace) -> int:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for path_str in args.inputs:
        path = Path(path_str)
        if not path.exists():
            logger.warning("Skipping missing input: %s", path)
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            for rec in csv.DictReader(handle):
                formula = (rec.get("formula") or "").strip()
                database = (rec.get("database") or "").strip()
                if not formula or not database:
                    continue
                key = (formula, database)
                if key in seen:
                    continue
                seen.add(key)
                rows.append({col: (rec.get(col) or "") for col in COLUMNS})
    logger.info("Merged %d unique (formula, database) rows from %d inputs", len(rows), len(args.inputs))
    write_rows(rows, Path(args.output))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract known formulas from local COD/OQMD datasets (no DB server).")
    sub = parser.add_subparsers(dest="command", required=True)

    cod = sub.add_parser("cod", help="extract formulas from a local COD CIF mirror")
    cod.add_argument("--cod-dir", default=str(ROOT / "datasets" / "COD-Dataset" / "cif"))
    cod.add_argument("--output", default="data/formulas/known_cod.csv")
    cod.add_argument("--limit", type=int, default=None, help="max CIF files to scan (testing)")
    cod.add_argument("--header-bytes", type=int, default=16384, help="bytes of each CIF header to read")
    cod.add_argument("--log-every", type=int, default=100000)
    cod.set_defaults(func=cmd_cod)

    oqmd = sub.add_parser("oqmd", help="extract formulas from the OQMD SQL dump")
    oqmd.add_argument("--dump", default=_default_oqmd_dump())
    oqmd.add_argument("--output", default="data/formulas/known_oqmd.csv")
    oqmd.add_argument("--limit", type=int, default=None, help="max composition rows to read (testing)")
    oqmd.add_argument("--log-every", type=int, default=500000)
    oqmd.set_defaults(func=cmd_oqmd)

    mp = sub.add_parser("mp", help="query the Materials Project API for formulas (needs an API key)")
    mp.add_argument("--api-key", default=None, help="MP API key (or set MP_API_KEY)")
    mp.add_argument("--output", default="data/formulas/known_mp.csv")
    mp.add_argument("--limit", type=int, default=None, help="max docs to read (testing)")
    mp.add_argument("--log-every", type=int, default=20000)
    mp.set_defaults(func=cmd_mp)

    merge = sub.add_parser("merge", help="merge per-database CSVs into known_formulas.csv")
    merge.add_argument("--inputs", nargs="+",
                       default=["data/formulas/known_cod.csv", "data/formulas/known_oqmd.csv",
                                "data/formulas/known_mp.csv"])
    merge.add_argument("--output", default="data/formulas/known_formulas.csv")
    merge.set_defaults(func=cmd_merge)
    return parser


def _default_oqmd_dump() -> str:
    base = ROOT / "datasets" / "OQMD-Dataset"
    extracted = base / "qmdb__v1_8__022026.sql" / "qmdb_20260129_release.sql"
    if extracted.exists():
        return str(extracted)
    gz = base / "qmdb__v1_8__022026.sql.gz"
    return str(gz)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())


