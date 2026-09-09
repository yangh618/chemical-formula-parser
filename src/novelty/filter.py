"""Membership-based novelty classification and heuristic candidate ranking."""

from __future__ import annotations

import logging
from typing import Any, Iterable

import numpy as np

from formula_utils import composition_distance, element_set_distance
from .canonicalize import CanonicalFormula

logger = logging.getLogger(__name__)

_DATABASE_COLUMNS = ("cod", "oqmd", "jarvis", "aflow", "mp")


class KnownIndex:
    """Precomputed lookups over the known-formula database for fast screening."""

    def __init__(self) -> None:
        self.formulas: dict[str, set[str]] = {}      # db -> canonical formulas
        self.all_formulas: set[str] = set()
        self.element_combos: set[frozenset[str]] = set()
        self.anonymous: set[str] = set()
        self._vectors: dict[str, np.ndarray] = {}
        self._known_formulas_order: list[str] = []

    def add(self, database: str, formula: str, anonymous: str = "", elements: Iterable[str] = ()) -> None:
        db = database.lower()
        self.formulas.setdefault(db, set()).add(formula)
        self.all_formulas.add(formula)
        if anonymous:
            self.anonymous.add(anonymous)
        els = frozenset(elements)
        if els:
            self.element_combos.add(els)

    def formula_in(self, formula: str, database: str | None = None) -> bool:
        if database is None:
            return formula in self.all_formulas
        return formula in self.formulas.get(database.lower(), set())

    def databases_containing(self, formula: str) -> list[str]:
        return [db for db, vals in self.formulas.items() if formula in vals]

    def finalize(self) -> None:
        """Cache fractional vectors for distance queries (sampled if huge)."""
        from formula_utils import get_elemental_fraction_vector

        order = sorted(self.all_formulas)
        if len(order) > 50000:
            logger.info("Known set large (%d); distance queries use a 50k sample.", len(order))
            order = order[:50000]
        self._known_formulas_order = order
        for formula in order:
            try:
                self._vectors[formula] = get_elemental_fraction_vector(formula)
            except Exception:  # noqa: BLE001
                continue

    @property
    def vector_matrix(self) -> np.ndarray:
        if not self._known_formulas_order:
            return np.zeros((0, 118))
        return np.vstack([self._vectors[f] for f in self._known_formulas_order if f in self._vectors])


def build_known_index(frame: "Any") -> KnownIndex:
    """Build a :class:`KnownIndex` from a known-formulas DataFrame."""
    index = KnownIndex()
    if frame is None or (hasattr(frame, "empty") and frame.empty):
        return index
    required = {"formula", "database"}
    if not required.issubset(frame.columns):
        raise ValueError(f"known-formula frame must contain {sorted(required)} columns")

    has_anon = "anonymous_formula" in frame.columns
    has_elements = "elements" in frame.columns
    for _i, record in frame.iterrows():
        formula = record["formula"]
        db = record["database"]
        anonymous = record["anonymous_formula"] if has_anon else ""
        elements = record["elements"].split() if has_elements and isinstance(record["elements"], str) else []
        index.add(str(db), str(formula), anonymous=str(anonymous or ""), elements=elements)
    index.finalize()
    return index


def annotate(candidate: CanonicalFormula, index: KnownIndex) -> dict[str, Any]:
    """Compute novelty flags for a canonical candidate (sections 15-16)."""
    row: dict[str, Any] = candidate.to_row()
    for db in _DATABASE_COLUMNS:
        row[f"found_in_{db}"] = index.formula_in(candidate.formula, database=db)
    row["found_in_any"] = index.formula_in(candidate.formula)
    row["novel_exact_formula"] = not row["found_in_any"]
    row["novel_element_combination"] = frozenset(candidate.elements) not in index.element_combos
    row["novel_anonymous_formula"] = candidate.anonymous_formula not in index.anonymous
    return row


def nearest_known_distance(
    candidate: CanonicalFormula,
    index: KnownIndex,
    max_sample: int = 50000,
) -> dict[str, Any]:
    """Distance to the nearest known composition (simple baselines, sec. 18).

    Returns ``nearest_known_formula``, ``min_composition_distance`` and
    ``min_element_set_distance``.
    """
    mat = index.vector_matrix
    if mat.shape[0] == 0:
        return {"nearest_known_formula": None, "min_composition_distance": None,
                "min_element_set_distance": None}
    if mat.shape[0] > max_sample:
        mat = mat[:max_sample]
    delta = mat - candidate.fraction_vector[None, :]
    dists = np.linalg.norm(delta, axis=1)
    j = int(np.argmin(dists))
    nearest = index._known_formulas_order[j]
    el_dist = element_set_distance(candidate.formula, nearest)
    return {
        "nearest_known_formula": nearest,
        "min_composition_distance": float(dists[j]),
        "min_element_set_distance": float(el_dist),
    }


# ---------------------------------------------------------------------------
# Heuristic, deterministic ranking (no machine learning, section 17).
# ---------------------------------------------------------------------------

def _criteria(candidate: CanonicalFormula) -> dict[str, float]:
    """Modular sub-scores; extend here and register in ``rank_candidates``."""
    els = candidate.elements
    n_elements = candidate.num_elements
    n_el_score = 1.0 if n_elements in (3, 4) else (0.4 if n_elements == 2 else 0.6)
    f = candidate.fraction_vector
    f = f[f > 0]
    p = f / f.sum()
    entropy = float(-(p * np.log(p)).sum())
    stoich_score = min(entropy / np.log(max(len(p), 2)), 1.0)
    tm = {"Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zr", "Nb", "Mo", "Hf", "Ta", "W", "Re"}
    tm_elements = {el for el in els if el in tm}
    tm_score = min(len(tm_elements) / max(n_elements, 1), 1.0)
    return {"n_elements": n_el_score, "stoich_complexity": stoich_score, "tm_count": tm_score}


def rank_candidates(
    candidates: list[CanonicalFormula],
    weights: dict[str, float] | None = None,
) -> list[dict[str, Any]]:
    """Return candidate dicts with a deterministic ``rank_score`` (0..1).

    The score is a weighted mean of modular criteria; highest first. No
    machine-learning model is used (section 17).
    """
    weights = weights or {"n_elements": 1.0, "stoich_complexity": 1.0, "tm_count": 1.0}
    denom = sum(weights.values()) or 1.0
    rows: list[dict[str, Any]] = []
    for cand in candidates:
        sub = _criteria(cand)
        total = sum(weights.get(k, 0.0) * v for k, v in sub.items())
        row = cand.to_row()
        row["rank_score"] = round(total / denom, 4)
        for k, v in sub.items():
            row[f"score_{k}"] = round(v, 4)
        rows.append(row)
    return sorted(rows, key=lambda r: r["rank_score"], reverse=True)
