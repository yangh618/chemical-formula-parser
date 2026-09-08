"""Canonical formula object used by novelty and benchmark code."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ..formula_utils import (
    FormulaError,
    get_anonymous_formula,
    get_elemental_fraction_vector,
    get_elements,
    get_num_elements,
    get_stoichiometric_vector,
    normalize_formula,
)


@dataclass(frozen=True)
class CanonicalFormula:
    """A fully canonicalised composition."""

    formula: str                    # canonical reduced formula, e.g. "NbAl2Cr"
    anonymous_formula: str          # e.g. "ABC2"
    elements: tuple[str, ...]       # element symbols sorted by atomic number
    num_elements: int
    stoichiometry: str              # Z-ordered counts joined by ":", e.g. "2:1:1"
    counts: tuple[int, ...]         # Z-ordered stoichiometric counts
    fraction_vector: np.ndarray     # 118-D elemental fraction vector

    def to_row(self) -> dict:
        return {
            "formula": self.formula,
            "anonymous_formula": self.anonymous_formula,
            "elements": " ".join(self.elements),
            "num_elements": self.num_elements,
            "stoichiometry": self.stoichiometry,
        }


def canonicalize_formula(raw: str) -> CanonicalFormula:
    """Canonicalise *raw* into a :class:`CanonicalFormula`; raises on error."""
    reduced = normalize_formula(raw)
    anonymous = get_anonymous_formula(reduced)
    elements = tuple(get_elements(reduced))
    vec = get_stoichiometric_vector(reduced)
    counts = tuple(int(x) for x in vec if x)
    stoich = ":".join(str(c) for c in counts)
    return CanonicalFormula(
        formula=reduced,
        anonymous_formula=anonymous,
        elements=elements,
        num_elements=len(elements),
        stoichiometry=stoich,
        counts=counts,
        fraction_vector=get_elemental_fraction_vector(reduced),
    )


def safe_canonicalize_formula(raw: str) -> Optional[CanonicalFormula]:
    """``canonicalize_formula`` returning ``None`` instead of raising."""
    try:
        return canonicalize_formula(raw)
    except FormulaError:
        return None
