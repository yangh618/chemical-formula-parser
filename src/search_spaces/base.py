"""Base classes and shared utilities for chemical search spaces.

A *search space* defines a chemically meaningful region of composition space
(e.g. oxychalcogenides, phosphorus chalcogenides, double perovskites). Each
search space can:

* **parse** a formula string and return a structured :class:`FormulaRecord`,
* **validate** whether a formula belongs to the space,
* **generate** candidate compositions within the space,
* **canonicalize** formulas so that chemically equivalent inputs map to the
  same representation.

The primary design principle is that *database absence does not imply chemical
invalidity*. A formula being unknown to Materials Project / ICSD / COD / OQMD
is recorded as novelty information, never used as a rejection criterion.
"""

from __future__ import annotations

import itertools
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from pymatgen.core import Composition
from pymatgen.core.periodic_table import Element

from formula_utils import (
    FormulaError,
    _reduced_amounts,
    get_anonymous_formula,
    get_elemental_fraction_vector,
    get_elements,
    normalize_formula,
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class FormulaRecord:
    """Structured representation of a parsed formula within a search space."""

    formula: str
    canonical_formula: str
    system: str
    elements: list[str]
    element_count: int
    cation_species: list[str]
    anion_species: list[str]
    stoichiometry: dict[str, int]
    atomic_fractions: dict[str, float]
    cation_fractions: dict[str, float]
    anion_fractions: dict[str, float]
    prototype: str
    mixed_occupancy: bool
    application_domains: list[str]
    oxidation_states: dict[str, list[int]] = field(default_factory=dict)
    charge_balanced: bool | None = None
    warnings: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "formula": self.formula,
            "canonical_formula": self.canonical_formula,
            "system": self.system,
            "elements": self.elements,
            "element_count": self.element_count,
            "cation_species": self.cation_species,
            "anion_species": self.anion_species,
            "stoichiometry": self.stoichiometry,
            "atomic_fractions": self.atomic_fractions,
            "cation_fractions": self.cation_fractions,
            "anion_fractions": self.anion_fractions,
            "prototype": self.prototype,
            "mixed_occupancy": self.mixed_occupancy,
            "application_domains": self.application_domains,
            "oxidation_states": self.oxidation_states,
            "charge_balanced": self.charge_balanced,
            "warnings": self.warnings,
            **self.extra,
        }


# ---------------------------------------------------------------------------
# Abstract base class
# ---------------------------------------------------------------------------

class SearchSpace(ABC):
    """Abstract base for a chemical search space."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Short identifier, e.g. ``'oxychalcogenides'``."""

    @property
    @abstractmethod
    def application_domains(self) -> list[str]:
        """Labels describing target application areas."""

    @abstractmethod
    def parse(self, formula: str) -> FormulaRecord:
        """Parse *formula* and return a :class:`FormulaRecord`.

        Raises :class:`FormulaError` if the formula is syntactically invalid or
        contains elements outside the search-space pool.
        """

    @abstractmethod
    def is_valid(self, formula: str) -> bool:
        """Return ``True`` if *formula* belongs to this search space."""

    @abstractmethod
    def generate(self, **kwargs: Any) -> list[str]:
        """Generate candidate formula strings within this search space."""

    def validate_pool(self, formula: str) -> set[str]:
        """Return elements in *formula* that are outside the allowed pool."""
        try:
            elements = set(get_elements(formula))
        except FormulaError:
            return set()
        return elements - self.allowed_elements

    @property
    @abstractmethod
    def allowed_elements(self) -> set[str]:
        """Full set of element symbols permitted in this search space."""


# ---------------------------------------------------------------------------
# Oxidation-state data (common states for the elements used across all three
# search spaces). Used for charge-balance checking and stoichiometry generation.
# ---------------------------------------------------------------------------

COMMON_OXIDATION_STATES: dict[str, list[int]] = {
    # Rare earths
    "La": [3], "Ce": [3, 4], "Pr": [3, 4], "Nd": [3],
    "Sm": [2, 3], "Eu": [2, 3], "Gd": [3],
    # Alkaline earth
    "Ca": [2], "Sr": [2], "Ba": [2],
    # Transition metals
    "Ti": [3, 4], "V": [3, 4, 5], "Cr": [3, 6],
    "Mn": [2, 3, 4], "Fe": [2, 3], "Co": [2, 3],
    "Ni": [2, 3], "Cu": [1, 2],
    # Group 5/6
    "Nb": [5], "Mo": [4, 6], "W": [6],
    # Anions
    "O": [-2], "S": [-2], "Se": [-2], "Te": [-2],
    # Phosphorus
    "P": [5],
}


# ---------------------------------------------------------------------------
# Charge-balance utilities
# ---------------------------------------------------------------------------

def check_charge_balance(
    stoichiometry: dict[str, int],
    oxidation_states: dict[str, list[int]] | None = None,
) -> list[dict[str, int]]:
    """Find all oxidation-state assignments that achieve charge balance.

    Returns a list of ``{element: oxidation_state}`` dicts. An empty list means
    no charge-neutral assignment was found with the given states.
    """
    ox_states = oxidation_states or COMMON_OXIDATION_STATES
    elements = list(stoichiometry.keys())
    if not elements:
        return []

    possible_states = []
    for el in elements:
        states = ox_states.get(el, [])
        if not states:
            return []
        possible_states.append(states)

    balanced = []
    for combo in itertools.product(*possible_states):
        total = sum(stoichiometry[el] * ox for el, ox in zip(elements, combo))
        if abs(total) < 1e-6:
            balanced.append(dict(zip(elements, combo)))
    return balanced


def is_charge_balanced(
    stoichiometry: dict[str, int],
    oxidation_states: dict[str, list[int]] | None = None,
) -> bool:
    """Return ``True`` if at least one charge-neutral assignment exists."""
    return len(check_charge_balance(stoichiometry, oxidation_states)) > 0


# ---------------------------------------------------------------------------
# Stoichiometry generation
# ---------------------------------------------------------------------------

def generate_charge_balanced_stoichiometries(
    cations: list[str],
    anions: list[str],
    max_coeff: int = 8,
    max_total_atoms: int = 20,
    oxidation_states: dict[str, list[int]] | None = None,
) -> list[dict[str, int]]:
    """Enumerate small-integer, charge-neutral stoichiometries.

    For each combination of cation and anion coefficients (1 .. *max_coeff*),
    checks whether any oxidation-state assignment yields zero net charge.
    Results are reduced to simplest integer ratio and deduplicated.
    """
    ox_states = oxidation_states or COMMON_OXIDATION_STATES
    all_elements = cations + anions
    n = len(all_elements)
    if n == 0:
        return []

    seen: set[tuple[tuple[str, int], ...]] = set()
    results: list[dict[str, int]] = []

    for el in all_elements:
        if not ox_states.get(el):
            return []

    element_states = [ox_states[el] for el in all_elements]

    ranges = [range(1, max_coeff + 1)] * n
    for coeffs in itertools.product(*ranges):
        if sum(coeffs) > max_total_atoms:
            continue
        stoich = dict(zip(all_elements, coeffs))

        for ox_combo in itertools.product(*element_states):
            total = sum(c * o for c, o in zip(coeffs, ox_combo))
            if abs(total) < 1e-6:
                reduced = _reduce_stoichiometry(stoich)
                key = tuple(sorted(reduced.items()))
                if key not in seen:
                    seen.add(key)
                    results.append(reduced)
                break

    return results


def _reduce_stoichiometry(stoich: dict[str, int]) -> dict[str, int]:
    """Reduce a stoichiometry dict to simplest integer ratio."""
    from math import gcd
    values = [v for v in stoich.values() if v > 0]
    if not values:
        return dict(stoich)
    g = values[0]
    for v in values[1:]:
        g = gcd(g, v)
    if g <= 1:
        return dict(stoich)
    return {el: c // g for el, c in stoich.items()}


# ---------------------------------------------------------------------------
# Fraction helpers
# ---------------------------------------------------------------------------

def compute_fractions(stoichiometry: dict[str, int]) -> dict[str, float]:
    """Compute normalized atomic fractions from a stoichiometry dict."""
    total = sum(stoichiometry.values())
    if total <= 0:
        return {el: 0.0 for el in stoichiometry}
    return {el: count / total for el, count in stoichiometry.items()}


def compute_role_fractions(
    stoichiometry: dict[str, int],
    role_elements: set[str],
) -> dict[str, float]:
    """Compute fractions restricted to a subset of elements (e.g. cations only).

    Returns fractions that sum to 1 within the role subset.
    """
    role_counts = {el: c for el, c in stoichiometry.items() if el in role_elements}
    total = sum(role_counts.values())
    if total <= 0:
        return {el: 0.0 for el in role_counts}
    return {el: c / total for el, c in role_counts.items()}


# ---------------------------------------------------------------------------
# Mixed-occupancy expansion
# ---------------------------------------------------------------------------

def expand_mixed_occupancy(formula: str) -> list[str]:
    """Expand parenthesized mixed-occupancy notation into individual formulas.

    ``(La0.5Sr0.5)2FeMoO6`` -> ``LaFeMoO6``, ``SrFeMoO6`` (approximate).

    For formulas without parentheses, returns ``[formula]`` unchanged.
    """
    import re
    if "(" not in formula:
        return [formula]

    pattern = re.compile(r"\(([^)]+)\)(\d*\.?\d*)")
    expanded = [formula]

    for match in pattern.finditer(formula):
        inner = match.group(1)
        multiplier_str = match.group(2) or "1"
        multiplier = float(multiplier_str)

        parts = re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", inner)
        if not parts:
            continue

        new_expanded = []
        for base_formula in expanded:
            for elem, frac_str in parts:
                frac = float(frac_str) if frac_str else 1.0
                amount = frac * multiplier
                amount_int = int(round(amount))
                if amount_int <= 0:
                    amount_int = 1
                replacement = f"{elem}{amount_int}" if amount_int > 1 else elem
                new_formula = base_formula.replace(match.group(0), replacement, 1)
                new_expanded.append(new_formula)
        expanded = new_expanded

    return expanded


def has_mixed_occupancy(formula: str) -> bool:
    """Return ``True`` if *formula* contains parenthesized site notation."""
    return "(" in formula


# ---------------------------------------------------------------------------
# Canonicalization
# ---------------------------------------------------------------------------

def canonicalize_search_space_formula(formula: str) -> str:
    """Canonicalize a formula string (reduced, element-sorted).

    Delegates to :func:`formula_utils.normalize_formula` for the core work.
    """
    return normalize_formula(formula)


def safe_canonicalize(formula: str) -> str | None:
    """``canonicalize_search_space_formula`` returning ``None`` on failure."""
    try:
        return canonicalize_search_space_formula(formula)
    except (FormulaError, Exception):
        return None


# ---------------------------------------------------------------------------
# Oxidation-state plausibility check
# ---------------------------------------------------------------------------

def flag_unusual_oxidation_states(
    stoichiometry: dict[str, int],
    oxidation_states: dict[str, list[int]] | None = None,
    known_common: dict[str, list[int]] | None = None,
) -> list[str]:
    """Flag elements whose required oxidation state is unusual but not impossible.

    Returns a list of warning strings. Does NOT reject the formula.
    """
    ox = oxidation_states or COMMON_OXIDATION_STATES
    warnings = []
    balanced = check_charge_balance(stoichiometry, ox)
    if not balanced:
        warnings.append(
            "No common oxidation-state combination achieves charge balance"
        )
    return warnings
