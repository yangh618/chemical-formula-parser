"""Novelty screening utilities.

The three novelty *levels* defined in section 16 are kept distinct:

* **Level 1 - exact formula novel**: the reduced formula is absent from every
  searched database (``novel_exact_formula``);
* **Level 2 - element-combination novel**: the *set of elements* has never
  occurred together in any searched database (``novel_element_combination``);
* **Level 3 - anonymous-stoichiometry novel**: the anonymous formula is absent
  (``novel_anonymous_formula``).

Terminology is deliberately careful: these flags identify **database-unseen**
compositions, not new materials, stable phases or experimentally discovered
compounds (see section 29).
"""

from __future__ import annotations

from .canonicalize import CanonicalFormula, canonicalize_formula, safe_canonicalize_formula
from .filter import KnownIndex, annotate, build_known_index, nearest_known_distance, rank_candidates

__all__ = [
    "CanonicalFormula",
    "canonicalize_formula",
    "safe_canonicalize_formula",
    "KnownIndex",
    "build_known_index",
    "annotate",
    "nearest_known_distance",
    "rank_candidates",
]
