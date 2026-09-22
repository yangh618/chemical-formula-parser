"""Novelty assessment for generated composition-space formulas.

Classifies each formula by comparison against a reference database, but
*never* rejects a formula solely because it is database-unknown.

Novelty categories (increasing novelty):

``known_formula``
    Exact reduced formula present in at least one reference database.
``known_element_set``
    Element combination has been reported, but this specific stoichiometry
    has not.
``known_stoichiometry_type``
    Same anonymous (reduced) stoichiometry type exists, but with different
    elements.
``novel_element_combination``
    This exact set of elements has never been reported together.
``fully_novel``
    Both element combination and stoichiometry are new.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator

from formula_utils import (
    FormulaError,
    get_anonymous_formula,
    get_elements,
    normalize_formula,
    safe_flatten_formula,
)
from novelty.canonicalize import CanonicalFormula, canonicalize_formula
from novelty.filter import KnownIndex

logger = logging.getLogger(__name__)

# Default location of extracted database CSVs
_DEFAULT_DATASETS_DIR = Path(__file__).resolve().parent.parent / "datasets"


@dataclass
class NoveltyAssessment:
    """Result of novelty assessment for a single formula."""

    formula: str
    canonical_formula: str
    known_formula: bool
    known_element_set: bool
    known_stoichiometry: bool
    novelty_category: str
    databases_found_in: list[str]
    nearest_known_formula: str | None = None
    composition_distance: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "formula": self.formula,
            "canonical_formula": self.canonical_formula,
            "known_formula": self.known_formula,
            "known_element_set": self.known_element_set,
            "known_stoichiometry": self.known_stoichiometry,
            "novelty_category": self.novelty_category,
            "databases_found_in": self.databases_found_in,
            "nearest_known_formula": self.nearest_known_formula,
            "composition_distance": self.composition_distance,
        }


def assess_novelty(
    formula: str,
    known_index: KnownIndex,
) -> NoveltyAssessment:
    """Assess the novelty of *formula* against a :class:`KnownIndex`."""
    canonical = normalize_formula(formula)
    elements = frozenset(get_elements(canonical))
    anonymous = get_anonymous_formula(canonical)

    dbs_found = known_index.databases_containing(canonical)
    known_formula = len(dbs_found) > 0
    known_element_set = elements in known_index.element_combos
    known_stoich = anonymous in known_index.anonymous

    category = _classify_novelty(known_formula, known_element_set, known_stoich)

    nearest = None
    dist = None
    if not known_formula and known_index.vector_matrix.shape[0] > 0:
        try:
            cf = canonicalize_formula(canonical)
            from novelty.filter import nearest_known_distance
            result = nearest_known_distance(cf, known_index)
            nearest = result.get("nearest_known_formula")
            dist = result.get("min_composition_distance")
        except Exception:
            pass

    return NoveltyAssessment(
        formula=formula,
        canonical_formula=canonical,
        known_formula=known_formula,
        known_element_set=known_element_set,
        known_stoichiometry=known_stoich,
        novelty_category=category,
        databases_found_in=dbs_found,
        nearest_known_formula=nearest,
        composition_distance=dist,
    )


def assess_novelty_batch(
    formulas: list[str],
    known_index: KnownIndex,
) -> list[NoveltyAssessment]:
    """Assess novelty for a list of formulas. Skips invalid formulas."""
    results = []
    for f in formulas:
        try:
            results.append(assess_novelty(f, known_index))
        except (FormulaError, Exception):
            continue
    return results


def _classify_novelty(
    known_formula: bool,
    known_element_set: bool,
    known_stoichiometry: bool,
) -> str:
    if known_formula:
        return "known_formula"
    if known_element_set and known_stoichiometry:
        return "known_stoichiometry_type"
    if known_element_set:
        return "known_element_set"
    if known_stoichiometry:
        return "novel_element_combination"
    return "fully_novel"


def filter_novel_only(
    assessments: list[NoveltyAssessment],
) -> list[NoveltyAssessment]:
    """Return only assessments where the formula is not in any database."""
    return [a for a in assessments if not a.known_formula]


def filter_by_category(
    assessments: list[NoveltyAssessment],
    categories: set[str],
) -> list[NoveltyAssessment]:
    """Return assessments matching any of the given novelty categories."""
    return [a for a in assessments if a.novelty_category in categories]


def novelty_summary(assessments: list[NoveltyAssessment]) -> dict[str, int]:
    """Count formulas in each novelty category."""
    counts: dict[str, int] = {}
    for a in assessments:
        counts[a.novelty_category] = counts.get(a.novelty_category, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Loading known-formula CSVs into a KnownIndex
# ---------------------------------------------------------------------------

# Mapping from CSV filename prefix to database label
_CSV_DB_MAP = {
    "cod_formulas": "cod",
    "mp_formulas": "mp",
    "oqmd_formulas": "oqmd",
}

# Name of the merged, de-duplicated known-formula CSV written by
# ``scripts/merge_known_formulas.py`` into ``datasets/``: one row per unique
# formula across all selected databases.
MERGED_CSV_NAME = "known_formulas.csv"

# Database label used when a row carries no usable source information.
_UNLABELLED_DB = "known"


def comparison_key(formula: str) -> str | None:
    """Return the canonical, notation-independent key used to compare formulas.

    Extracted datasets store pymatgen reduced formulas that keep polyatomic
    groups (``Ac(AgAs)2``), whereas generated candidates are flattened
    (``AcAg2As2``). Flattening both sides makes the comparison exact; strings
    that are already parenthesis-free are returned unchanged because they are
    already canonical, which keeps multi-million-row loads fast.

    Returns ``None`` for empty/unparseable input.
    """
    text = (formula or "").strip()
    if not text:
        return None
    if "(" not in text:
        return text
    return safe_flatten_formula(text)


def merged_csv_path(datasets_dir: str | Path | None = None) -> Path:
    """Path of the merged unique-formula CSV (whether or not it exists yet)."""
    base = Path(datasets_dir) if datasets_dir is not None else _DEFAULT_DATASETS_DIR
    return base / MERGED_CSV_NAME


@dataclass
class KnownCsvMergeResult:
    """Outcome of merging the per-database formula CSVs."""

    output_path: Path
    per_database: dict[str, int]
    total_unique: int

    @property
    def overlap_removed(self) -> int:
        """How many per-database entries collapsed away as duplicates."""
        return sum(self.per_database.values()) - self.total_unique

    def to_dict(self) -> dict[str, Any]:
        return {
            "output_path": str(self.output_path),
            "per_database": dict(self.per_database),
            "total_unique": self.total_unique,
            "overlap_removed": self.overlap_removed,
        }


def formula_metadata(formula: str) -> tuple[str, str]:
    """Return ``(anonymous_formula, "Element Element ...")`` for *formula*.

    Derived with pymatgen; returns ``("", "")`` for formulas that cannot be
    parsed. Used to enrich the merged CSV so element-combination and
    anonymous-stoichiometry novelty can be assessed without re-parsing every
    known formula at load time.
    """
    try:
        anonymous = get_anonymous_formula(formula)
        elements = " ".join(get_elements(formula))
    except FormulaError:
        return "", ""
    except Exception:  # noqa: BLE001 - pymatgen raises several types
        return "", ""
    return anonymous, elements


def merge_database_csv_files(
    datasets_dir: str | Path | None = None,
    databases: list[str] | None = None,
    output_path: str | Path | None = None,
    with_metadata: bool = True,
) -> "KnownCsvMergeResult":
    """Merge ``{cod,mp,oqmd}_formulas.csv`` into one unique-formula CSV.

    Every formula appears exactly once; the ``databases`` column records which
    source databases contain it. Formulas are stored as canonical
    :func:`comparison_key` values so they match flattened candidate formulas.

    Parameters
    ----------
    datasets_dir
        Directory holding the per-database CSVs. Defaults to ``datasets/`` at
        the project root.
    databases
        Which databases to merge (e.g. ``["cod", "mp"]``). Defaults to every
        available CSV.
    output_path
        Where to write the merged CSV. Defaults to
        ``<datasets_dir>/known_formulas.csv``.
    with_metadata
        Also derive the ``anonymous_formula``, ``elements`` and
        ``num_elements`` columns (one-off pymatgen pass over the unique
        formulas). These make element-combination and anonymous-stoichiometry
        novelty assessments meaningful; set to False for a fast formula-only
        file.

    Returns
    -------
    KnownCsvMergeResult
        Output path plus per-database and deduplication statistics.
    """
    base = Path(datasets_dir) if datasets_dir is not None else _DEFAULT_DATASETS_DIR
    if not base.is_dir():
        raise FileNotFoundError(f"Datasets directory not found: {base}")

    labels = [db_label for _prefix, db_label in _CSV_DB_MAP.items()]
    masks: dict[str, int] = {}
    per_database: dict[str, int] = {}

    for bit, (csv_prefix, db_label) in enumerate(_CSV_DB_MAP.items()):
        if databases and db_label not in databases:
            continue
        csv_path = base / f"{csv_prefix}.csv"
        if not csv_path.exists():
            logger.warning("CSV not found, skipping: %s", csv_path)
            continue

        db_formulas = _read_formulas_from_csv(csv_path)
        per_database[db_label] = len(db_formulas)
        for formula in db_formulas:
            masks[formula] = masks.get(formula, 0) | (1 << bit)
        logger.info("Merged %d unique formulas from %s (%s)",
                    len(db_formulas), csv_path.name, db_label)

    out = Path(output_path) if output_path is not None else merged_csv_path(base)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["formula", "databases", "anonymous_formula",
                         "elements", "num_elements"])
        for formula in sorted(masks):
            mask = masks[formula]
            sources = [labels[i] for i in range(len(labels)) if mask & (1 << i)]
            if with_metadata:
                anonymous, elements = formula_metadata(formula)
                num_elements = len(elements.split()) if elements else ""
            else:
                anonymous, elements, num_elements = "", "", ""
            writer.writerow([formula, ";".join(sources), anonymous, elements,
                             num_elements])

    return KnownCsvMergeResult(
        output_path=out,
        per_database=per_database,
        total_unique=len(masks),
    )


def _row_database_labels(row: dict[str, str], default_label: str) -> list[str]:
    """Database labels for a CSV row (``databases`` or ``database`` column)."""
    raw = (row.get("databases") or "").strip()
    if raw:
        return [part.strip().lower() for part in raw.replace(",", ";").split(";")
                if part.strip()]
    label = (row.get("database") or "").strip().lower()
    return [label or default_label]


def iter_known_csv(
    csv_path: str | Path,
    default_label: str = _UNLABELLED_DB,
) -> Iterator[tuple[str, list[str], str, str]]:
    """Stream a known-formulas CSV as ``(key, databases, anonymous, elements)``.

    Both supported schemas are accepted:

    * merged schema — ``formula, databases`` (``databases`` is ``;``-joined);
    * per-database schema — ``formula, database [, anonymous_formula, elements]``.

    Formulas are normalised through :func:`comparison_key`; rows without a
    usable formula are skipped.
    """
    with open(csv_path, "r", newline="", encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f):
            key = comparison_key(row.get("formula") or "")
            if key is None:
                continue
            anonymous = (row.get("anonymous_formula") or "").strip()
            elements = (row.get("elements") or "").strip()
            yield key, _row_database_labels(row, default_label), anonymous, elements


def _read_formulas_from_csv(csv_path: str | Path) -> set[str]:
    """Return the set of canonical formulas stored in *csv_path*."""
    return {key for key, _labels, _anon, _els in iter_known_csv(csv_path)}


def _read_known_index(csv_path: str | Path, default_label: str = _UNLABELLED_DB) -> KnownIndex:
    """Build a :class:`KnownIndex` from a single known-formulas CSV."""
    index = KnownIndex()
    for key, labels, anonymous, elements in iter_known_csv(csv_path, default_label):
        for label in labels:
            index.add(label, key, anonymous=anonymous, elements=elements.split())
    index.finalize()
    return index


def _resolve_known_source(
    datasets_dir: str | Path | None,
    databases: list[str] | None,
    csv_path: str | Path | None,
) -> Path | None:
    """Pick the single CSV to read, or ``None`` to fall back to per-database files.

    An explicit *csv_path* always wins; otherwise the merged CSV is used when it
    exists and no database subset was requested.
    """
    if csv_path is not None:
        explicit = Path(csv_path)
        if not explicit.exists():
            raise FileNotFoundError(f"Known-formulas CSV not found: {explicit}")
        return explicit

    merged = merged_csv_path(datasets_dir)
    if databases is None and merged.exists():
        return merged
    return None


def load_known_index_from_csvs(
    datasets_dir: str | Path | None = None,
    databases: list[str] | None = None,
    csv_path: str | Path | None = None,
) -> KnownIndex:
    """Build a :class:`KnownIndex` from known-formula CSV files.

    Parameters
    ----------
    datasets_dir
        Directory containing the extracted per-database CSVs
        (``{cod,mp,oqmd}_formulas.csv``) and/or the merged
        ``known_formulas.csv``. Defaults to ``datasets/`` at the project root.
    databases
        Which databases to load (e.g. ``["cod", "mp"]``). Defaults to all
        available CSVs. Mutually exclusive with *csv_path*.
    csv_path
        Read this single CSV instead of ``datasets_dir``. Both schemas are
        accepted: the merged ``formula, databases`` file produced by
        ``scripts/merge_known_formulas.py`` and the per-database
        ``formula, database [, anonymous_formula, elements]`` schema.
        Takes precedence over *databases*.

    Returns
    -------
    KnownIndex
        Ready-to-use index with all formulas loaded. Call
        :meth:`~KnownIndex.formula_in` to check membership.
    """
    source = _resolve_known_source(datasets_dir, databases, csv_path)
    if source is not None:
        logger.info("Loading known formulas from %s", source)
        return _read_known_index(source)

    base = Path(datasets_dir) if datasets_dir is not None else _DEFAULT_DATASETS_DIR
    if not base.is_dir():
        raise FileNotFoundError(f"Datasets directory not found: {base}")

    index = KnownIndex()

    for csv_prefix, db_label in _CSV_DB_MAP.items():
        if databases and db_label not in databases:
            continue
        csv_file = base / f"{csv_prefix}.csv"
        if not csv_file.exists():
            logger.warning("CSV not found: %s", csv_file)
            continue

        count = 0
        for key, _labels, anonymous, elements in iter_known_csv(csv_file, default_label=db_label):
            index.add(db_label, key, anonymous=anonymous, elements=elements.split())
            count += 1
        logger.info("Loaded %d formulas from %s (%s)", count, csv_file.name, db_label)

    index.finalize()
    total = len(index.all_formulas)
    logger.info("KnownIndex ready: %d unique formulas across %d databases",
                total, len(index.formulas))
    return index


def load_known_formulas_set(
    datasets_dir: str | Path | None = None,
    databases: list[str] | None = None,
    csv_path: str | Path | None = None,
) -> set[str]:
    """Load all known formulas into a plain set for fast membership testing.

    This is a lightweight alternative to :func:`load_known_index_from_csvs`
    when you only need exact-formula lookups (no element-combo or anonymous
    formula checks).

    If ``datasets/known_formulas.csv`` (the merged, de-duplicated file built by
    ``scripts/merge_known_formulas.py``) exists it is read directly and the
    individual ``{cod,mp,oqmd}_formulas.csv`` files are ignored; pass
    *databases* to force a subset of the per-database files, or *csv_path* to
    read one specific CSV.

    Formulas are returned as canonical :func:`comparison_key` strings, i.e.
    flattened, so they compare exactly against the flattened candidate formulas
    produced by the generation scripts.
    """
    source = _resolve_known_source(datasets_dir, databases, csv_path)
    if source is not None:
        logger.info("Loading known formulas from %s", source)
        return _read_formulas_from_csv(source)

    base = Path(datasets_dir) if datasets_dir is not None else _DEFAULT_DATASETS_DIR
    if not base.is_dir():
        raise FileNotFoundError(f"Datasets directory not found: {base}")

    formulas: set[str] = set()

    for csv_prefix, db_label in _CSV_DB_MAP.items():
        if databases and db_label not in databases:
            continue
        csv_file = base / f"{csv_prefix}.csv"
        if not csv_file.exists():
            continue
        formulas.update(_read_formulas_from_csv(csv_file))

    return formulas


@dataclass
class NovelGenerationResult:
    """Outcome of a novelty-filtered candidate generation run."""

    novel: list[str]
    total_generated: int
    known_filtered: int
    novel_found: int
    target: int
    space_exhausted: bool
    limit_reached: int
    target_met: bool
    rejected_by_validator: int = 0

    @property
    def novel_fraction(self) -> float:
        """Fraction of generated candidates that were database-unseen."""
        if not self.total_generated:
            return 0.0
        return self.novel_found / self.total_generated


def generate_novel_candidates(
    space,
    num_novel: int,
    known_formulas: set[str],
    base_limit: int | None = None,
    max_limit: int = 2_000_000,
    growth: float = 4.0,
    verbose: bool = True,
    validator: Callable[[str], bool] | None = None,
    **generate_kwargs,
) -> "NovelGenerationResult":
    """Generate *num_novel* candidates from *space* that are absent from *known_formulas*.

    Because generation is capped by ``max_results`` and the enumeration order is
    deterministic, this progressively scales the generation limit until enough
    novel formulas have been collected or *max_limit* is reached.

    Parameters
    ----------
    space
        A :class:`~composition.CompositionSpace` instance.
    num_novel
        Number of novel (database-unseen) formulas wanted.
    known_formulas
        Set of canonical reduced formulas already reported.
    base_limit
        Starting generation limit. Defaults to ``num_novel * 3``.
    max_limit
        Hard cap on the generation limit, to bound runtime.
    growth
        Multiplicative factor applied to the limit each round.
    verbose
        Log progress per round.
    validator
        Optional chemical-plausibility predicate applied to every generated
        candidate *before* the novelty test (e.g. ``smactlib.smact_valid``).
        Candidates for which it returns False are dropped and counted in
        ``rejected_by_validator``.
    generate_kwargs
        Extra keyword arguments forwarded to ``space.generate`` (e.g.
        ``per_combo_limit``, ``max_cations``, ``max_anions``).

    Returns
    -------
    NovelGenerationResult
        Holds the novel formulas (at most *num_novel*) plus statistics: how many
        candidates were generated, how many were dropped as known, and whether
        the composition space was fully exhausted before reaching the target.
    """
    limit = base_limit or max(num_novel * 3, 100)

    while True:
        candidates = space.generate(max_results=limit, **generate_kwargs)
        total_generated = len(candidates)

        if validator is not None:
            survivors = [f for f in candidates if validator(f)]
            rejected = total_generated - len(survivors)
        else:
            survivors, rejected = candidates, 0

        novel = [f for f in survivors if f not in known_formulas]
        known_filtered = len(survivors) - len(novel)

        if verbose:
            logger.info("limit=%d generated=%d valid=%d novel=%d (target=%d)",
                        limit, total_generated, len(survivors), len(novel), num_novel)

        exhausted = total_generated < limit
        reached_cap = limit >= max_limit
        enough = len(novel) >= num_novel

        if enough or exhausted or reached_cap:
            return NovelGenerationResult(
                novel=novel[:num_novel],
                total_generated=total_generated,
                known_filtered=known_filtered,
                novel_found=len(novel),
                target=num_novel,
                space_exhausted=exhausted,
                limit_reached=limit,
                target_met=enough,
                rejected_by_validator=rejected,
            )

        limit = min(int(limit * growth), max_limit)
