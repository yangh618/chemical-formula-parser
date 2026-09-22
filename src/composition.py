"""Generic composition-space definition driven entirely by configuration.

A *composition space* is defined by a config dict (typically loaded from
``config.yaml``) that specifies:

* **cations** — element symbols acting as positive-charge species
* **anions** — element symbols acting as negative-charge species
* **oxidation_states** — mapping of element → list of allowed oxidation states
* **constraints** — optional rules like required elements, max stoich, etc.

The :class:`CompositionSpace` class provides parse / validate / generate
operations over any user-defined chemistry. No chemistry is hardcoded.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field
from math import gcd
from typing import Any

import numpy as np
from pymatgen.core import Composition
from pymatgen.core.periodic_table import Element

from config import get as _cfg
from formula_utils import (
    FormulaError,
    _reduced_amounts,
    get_anonymous_formula,
    get_elements,
    normalize_formula,
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class FormulaRecord:
    """Structured representation of a parsed formula within a composition space."""

    formula: str
    canonical_formula: str
    space_name: str
    elements: list[str]
    element_count: int
    cation_species: list[str]
    anion_species: list[str]
    stoichiometry: dict[str, int]
    atomic_fractions: dict[str, float]
    charge_balanced: bool
    oxidation_states_used: dict[str, list[int]]
    warnings: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "formula": self.formula,
            "canonical_formula": self.canonical_formula,
            "space_name": self.space_name,
            "elements": self.elements,
            "element_count": self.element_count,
            "cation_species": self.cation_species,
            "anion_species": self.anion_species,
            "stoichiometry": self.stoichiometry,
            "atomic_fractions": self.atomic_fractions,
            "charge_balanced": self.charge_balanced,
            "oxidation_states_used": self.oxidation_states_used,
            "warnings": self.warnings,
            **self.extra,
        }


# ---------------------------------------------------------------------------
# Generic composition space
# ---------------------------------------------------------------------------

class CompositionSpace:
    """A user-defined composition space loaded from configuration.

    Parameters
    ----------
    name
        Identifier for this space (e.g. ``'oxychalcogenides'``).
    config
        Dict with keys ``cations``, ``anions``, ``oxidation_states``, and
        optional ``constraints``. Typically loaded from ``config.yaml``.

    Example config::

        oxychalcogenides:
          cations: [La, Ce, Ti, V, Cr, Mn, Fe, Co, Ni]
          anions: [O, S, Se, Te]
          oxidation_states:
            La: [3]
            Ti: [3, 4]
            O: [-2]
            S: [-2]
          constraints:
            required_anions: [O]
            max_stoich: 8
            max_total_atoms: 20
    """

    def __init__(self, name: str, config: dict[str, Any]) -> None:
        self._name = name
        self._config = dict(config)

        self._cations: list[str] = list(config.get("cations", []))
        self._anions: list[str] = list(config.get("anions", []))
        self._oxidation_states: dict[str, list[int]] = {
            el: list(states)
            for el, states in config.get("oxidation_states", {}).items()
        }

        constraints = config.get("constraints", {}) or {}
        self._required_anions: list[str] = list(constraints.get("required_anions", []))
        self._required_cations: list[str] = list(constraints.get("required_cations", []))
        self._max_stoich: int = int(constraints.get("max_stoich", 8))
        self._max_total_atoms: int = int(constraints.get("max_total_atoms", 20))
        self._max_results: int = int(constraints.get("max_results", 5000))
        # Optional species-count caps (0 / absent = no cap). These encode the
        # chemistry of the system — e.g. a double perovskite A2BB'O6 has at most
        # four distinct cation species — and keep enumeration tractable.
        self._max_cation_species: int = int(constraints.get("max_cation_species", 0) or 0)
        self._max_anion_species: int = int(constraints.get("max_anion_species", 0) or 0)
        self._per_combo_limit: int = int(constraints.get("per_combo_limit", 0) or 0)

        self._all_elements = set(self._cations) | set(self._anions)
        self._cation_set = set(self._cations)
        self._anion_set = set(self._anions)

    @property
    def name(self) -> str:
        return self._name

    @property
    def cations(self) -> list[str]:
        return list(self._cations)

    @property
    def anions(self) -> list[str]:
        return list(self._anions)

    @property
    def allowed_elements(self) -> set[str]:
        return set(self._all_elements)

    @property
    def oxidation_states(self) -> dict[str, list[int]]:
        return dict(self._oxidation_states)

    @property
    def constraints(self) -> dict[str, Any]:
        """Return the constraints dict as loaded from config."""
        return {
            "required_anions": list(self._required_anions),
            "required_cations": list(self._required_cations),
            "max_stoich": self._max_stoich,
            "max_total_atoms": self._max_total_atoms,
            "max_results": self._max_results,
            "max_cation_species": self._max_cation_species,
            "max_anion_species": self._max_anion_species,
            "per_combo_limit": self._per_combo_limit,
        }

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    def parse(self, formula: str) -> FormulaRecord:
        """Parse *formula* and return a :class:`FormulaRecord`.

        Raises :class:`FormulaError` if the formula is invalid or contains
        elements outside this composition space.
        """
        if not formula or not formula.strip():
            raise FormulaError("Empty formula")

        text = formula.strip()
        elements = get_elements(text)
        element_set = set(elements)

        outside = element_set - self._all_elements
        if outside:
            raise FormulaError(
                f"Elements outside composition space {self._name!r}: {sorted(outside)}"
            )

        cations = sorted(element_set & self._cation_set, key=_by_z)
        anions = sorted(element_set & self._anion_set, key=_by_z)

        if not cations:
            raise FormulaError(f"No cation species found (pool: {self._cations})")
        if not anions:
            raise FormulaError(f"No anion species found (pool: {self._anions})")

        for req in self._required_anions:
            if req not in element_set:
                raise FormulaError(f"Required anion {req!r} not found")
        for req in self._required_cations:
            if req not in element_set:
                raise FormulaError(f"Required cation {req!r} not found")

        stoich = _reduced_amounts(text)
        canonical = normalize_formula(text)

        warnings: list[str] = []
        balanced = check_charge_balance(stoich, self._oxidation_states)
        charge_ok = len(balanced) > 0

        if not charge_ok:
            warnings.append("No common oxidation-state combination achieves charge balance")

        ox_used: dict[str, list[int]] = {}
        if balanced:
            for el in elements:
                states = {b[el] for b in balanced if el in b}
                ox_used[el] = sorted(states)

        total = sum(stoich.values())
        fractions = {el: c / total for el, c in stoich.items()} if total > 0 else {}

        return FormulaRecord(
            formula=text,
            canonical_formula=canonical,
            space_name=self._name,
            elements=elements,
            element_count=len(elements),
            cation_species=cations,
            anion_species=anions,
            stoichiometry=stoich,
            atomic_fractions=fractions,
            charge_balanced=charge_ok,
            oxidation_states_used=ox_used,
            warnings=warnings,
        )

    def is_valid(self, formula: str) -> bool:
        """Return ``True`` if *formula* belongs to this composition space."""
        try:
            self.parse(formula)
            return True
        except (FormulaError, Exception):
            return False

    # ------------------------------------------------------------------
    # Generation
    # ------------------------------------------------------------------

    def generate(
        self,
        min_cations: int = 1,
        max_cations: int | None = None,
        min_anions: int = 1,
        max_anions: int | None = None,
        max_coeff: int | None = None,
        max_total_atoms: int | None = None,
        max_results: int | None = None,
        per_combo_limit: int | None = None,
        max_scan: int = 20000,
    ) -> list[str]:
        """Generate candidate formulas within this composition space.

        Parameters
        ----------
        min_cations, max_cations
            Range for number of cation species. Default: 1 up to the config
            ``max_cation_species`` if set, otherwise ``len(cations)``.
        min_anions, max_anions
            Range for number of anion species. Default: 1 up to the config
            ``max_anion_species`` if set, otherwise ``len(anions)``.
        max_coeff
            Maximum stoichiometric coefficient. Default: from config constraints.
        max_total_atoms
            Maximum total atoms in formula. Default: from config constraints.
        max_results
            Cap on returned formulas. Default: from config constraints.
        per_combo_limit
            If given, keep at most this many formulas from each individual
            (cation-set, anion-set) combination; defaults to the config
            ``per_combo_limit`` when set. Enumeration otherwise exhausts
            every stoichiometry of the first few element combinations before
            moving on, so the head of the result list is dominated by a handful
            of species. Capping per combination spreads the output across the
            whole element-combination space, which is what diversity-oriented
            selection needs, and is also much faster because each combination
            stops early.
        max_scan
            Work budget: maximum number of coefficient tuples examined per
            element combination. Guards against combinatorial blow-up for
            combinations with many species.
        """
        max_cations = max_cations or self._max_cation_species or len(self._cations)
        max_anions = max_anions or self._max_anion_species or len(self._anions)
        max_coeff = max_coeff or self._max_stoich
        max_total_atoms = max_total_atoms or self._max_total_atoms
        max_results = max_results or self._max_results
        if per_combo_limit is None:
            per_combo_limit = self._per_combo_limit or None

        results: list[str] = []
        seen: set[str] = set()

        for n_cat in range(min_cations, max_cations + 1):
            for cat_combo in itertools.combinations(self._cations, n_cat):
                cations = list(cat_combo)

                for n_an in range(min_anions, max_anions + 1):
                    for an_combo in itertools.combinations(self._anions, n_an):
                        anions = list(an_combo)

                        # Check required elements
                        anion_set = set(anions)
                        if not set(self._required_anions).issubset(anion_set):
                            continue

                        stoichs = generate_charge_balanced_stoichiometries(
                            cations=cations,
                            anions=anions,
                            max_coeff=max_coeff,
                            max_total_atoms=max_total_atoms,
                            oxidation_states=self._oxidation_states,
                            max_results=per_combo_limit,
                            max_scan=max_scan,
                        )

                        for stoich in stoichs:
                            formula = _stoich_to_formula(stoich)
                            try:
                                canonical = normalize_formula(formula)
                            except FormulaError:
                                continue
                            if canonical not in seen:
                                seen.add(canonical)
                                results.append(canonical)
                                if len(results) >= max_results:
                                    return results

        return results


# ---------------------------------------------------------------------------
# Charge-balance utilities
# ---------------------------------------------------------------------------

def check_charge_balance(
    stoichiometry: dict[str, int],
    oxidation_states: dict[str, list[int]],
) -> list[dict[str, int]]:
    """Find all oxidation-state assignments that achieve charge balance.

    Returns a list of ``{element: oxidation_state}`` dicts. Empty list means
    no charge-neutral assignment exists.
    """
    elements = list(stoichiometry.keys())
    if not elements:
        return []

    possible_states = []
    for el in elements:
        states = oxidation_states.get(el, [])
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
    oxidation_states: dict[str, list[int]],
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
    max_results: int | None = None,
    max_scan: int | None = None,
) -> list[dict[str, int]]:
    """Enumerate small-integer, charge-neutral stoichiometries.

    For each combination of coefficients (1 .. *max_coeff*), checks whether
    any oxidation-state assignment yields zero net charge.

    Parameters
    ----------
    max_results
        Stop once this many distinct reduced stoichiometries are found.
    max_scan
        Stop after examining this many coefficient tuples, bounding the work
        spent on any single element combination.
    """
    if oxidation_states is None:
        return []

    all_elements = cations + anions
    n = len(all_elements)
    if n == 0:
        return []

    for el in all_elements:
        if not oxidation_states.get(el):
            return []

    element_states = [oxidation_states[el] for el in all_elements]

    seen: set[tuple[tuple[str, int], ...]] = set()
    results: list[dict[str, int]] = []

    ranges = [range(1, max_coeff + 1)] * n
    scanned = 0
    for coeffs in itertools.product(*ranges):
        if max_scan is not None and scanned >= max_scan:
            break
        scanned += 1

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
                    if max_results is not None and len(results) >= max_results:
                        return results
                break

    return results


def _reduce_stoichiometry(stoich: dict[str, int]) -> dict[str, int]:
    """Reduce a stoichiometry dict to simplest integer ratio."""
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
# Formula helpers
# ---------------------------------------------------------------------------

def _by_z(symbol: str) -> int:
    return Element(symbol).Z


def _stoich_to_formula(stoich: dict[str, int]) -> str:
    parts = []
    for el in sorted(stoich, key=_by_z):
        count = stoich[el]
        if count == 1:
            parts.append(el)
        else:
            parts.append(f"{el}{count}")
    return "".join(parts)


def has_mixed_occupancy(formula: str) -> bool:
    """Return ``True`` if *formula* contains parenthesized site notation."""
    return "(" in formula


def expand_mixed_occupancy(formula: str) -> list[str]:
    """Expand parenthesized mixed-occupancy notation into individual formulas."""
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
