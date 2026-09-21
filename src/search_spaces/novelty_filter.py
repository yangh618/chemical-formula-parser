"""Novelty-oriented filtering for search-space generated formulas.

This module sits between the search-space parser and the CrystalFNO structure
generator. It classifies each formula into a novelty category based on
comparison against a reference database, but *never* rejects a formula solely
because it is database-unknown.

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
``novel_stoichiometry``
    New stoichiometric ratio for a known element set.
``fully_novel``
    Both element combination and stoichiometry are new.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from formula_utils import (
    FormulaError,
    get_anonymous_formula,
    get_elements,
    normalize_formula,
)
from novelty.canonicalize import CanonicalFormula, canonicalize_formula
from novelty.filter import KnownIndex


# ---------------------------------------------------------------------------
# Novelty record
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Core assessment
# ---------------------------------------------------------------------------

def assess_novelty(
    formula: str,
    known_index: KnownIndex,
) -> NoveltyAssessment:
    """Assess the novelty of *formula* against a :class:`KnownIndex`.

    This function never raises for chemically valid formulas. If the formula
    is syntactically invalid, a :class:`FormulaError` propagates.
    """
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


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------

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
