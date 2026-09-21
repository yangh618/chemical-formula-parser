"""Double perovskite oxide search space: A2BB'O6.

Chemical space
--------------
* A-site: rare-earth + alkaline-earth cations
* B/B'-site: transition metals (ordered or mixed)
* Anion: oxygen only (fixed at 6 in the prototype)

Supports mixed-occupancy notation ``(La,Sr)2FeMoO6`` and distinguishes
A-site vs B-site species.
"""

from __future__ import annotations

import itertools
import re
from typing import Any

from formula_utils import FormulaError, _reduced_amounts, get_elements

from .base import (
    FormulaRecord,
    SearchSpace,
    canonicalize_search_space_formula,
    check_charge_balance,
    compute_fractions,
    compute_role_fractions,
    flag_unusual_oxidation_states,
    has_mixed_occupancy,
)

# ---------------------------------------------------------------------------
# Element pools
# ---------------------------------------------------------------------------

A_SITE_POOL: list[str] = [
    "La", "Nd", "Sm", "Eu", "Gd",
    "Ca", "Sr", "Ba",
]
B_SITE_POOL: list[str] = [
    "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Nb", "Mo", "W",
]

ALL_ALLOWED = set(A_SITE_POOL + B_SITE_POOL + ["O"])

_APPLICATION_DOMAINS = [
    "photovoltaics",
    "photocatalysis",
    "oxygen evolution catalysis",
    "hydrogen evolution catalysis",
    "magnetic materials",
    "spintronics",
    "thermoelectrics",
    "dielectric materials",
    "energy storage",
]

_OXIDATION_STATES = {
    "La": [3], "Nd": [3], "Sm": [2, 3], "Eu": [2, 3], "Gd": [3],
    "Ca": [2], "Sr": [2], "Ba": [2],
    "Ti": [3, 4], "V": [3, 4, 5], "Cr": [3, 6],
    "Mn": [2, 3, 4], "Fe": [2, 3], "Co": [2, 3],
    "Ni": [2, 3], "Cu": [1, 2],
    "Nb": [5], "Mo": [4, 6], "W": [6],
    "O": [-2],
}


class DoublePerovskiteSpace(SearchSpace):
    """Double perovskite oxide A2BB'O6 search space."""

    @property
    def name(self) -> str:
        return "double_perovskites"

    @property
    def application_domains(self) -> list[str]:
        return list(_APPLICATION_DOMAINS)

    @property
    def allowed_elements(self) -> set[str]:
        return ALL_ALLOWED

    @property
    def a_site_pool(self) -> list[str]:
        return list(A_SITE_POOL)

    @property
    def b_site_pool(self) -> list[str]:
        return list(B_SITE_POOL)

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    def parse(self, formula: str) -> FormulaRecord:
        if not formula or not formula.strip():
            raise FormulaError("Empty formula")

        text = formula.strip()
        mixed = has_mixed_occupancy(text)

        parseable = _expand_for_pymatgen(text) if mixed else text

        elements = get_elements(parseable)
        element_set = set(elements)

        outside = element_set - ALL_ALLOWED
        if outside:
            raise FormulaError(
                f"Elements outside double-perovskite pool: {sorted(outside)}"
            )

        if "O" not in element_set:
            raise FormulaError("Double perovskites require oxygen")

        a_species = sorted(
            element_set & set(A_SITE_POOL), key=lambda s: _by_z(s)
        )
        b_species = sorted(
            element_set & set(B_SITE_POOL), key=lambda s: _by_z(s)
        )

        if not a_species:
            raise FormulaError("No A-site species found")
        if not b_species:
            raise FormulaError("No B-site species found")
        if len(b_species) < 1:
            raise FormulaError("Double perovskite requires at least one B-site species")

        stoich = _reduced_amounts(parseable)
        canonical = canonicalize_search_space_formula(parseable)

        warnings: list[str] = []

        o_count = stoich.get("O", 0)
        total_a = sum(stoich.get(a, 0) for a in a_species)
        total_b = sum(stoich.get(b, 0) for b in b_species)

        if mixed:
            warnings.append("Mixed-occupancy notation detected")
            site_info = _analyze_mixed_sites(text)
            if site_info:
                warnings.append(f"Mixed sites: {site_info}")

        prototype = _classify_prototype(stoich, a_species, b_species, total_a, total_b, o_count)

        balanced = check_charge_balance(stoich, _OXIDATION_STATES)
        charge_ok = len(balanced) > 0

        if not charge_ok:
            warnings.append("Unusual oxidation-state combination flagged for review")

        ox_warnings = flag_unusual_oxidation_states(stoich, _OXIDATION_STATES)
        warnings.extend(ox_warnings)

        ox_states_used = {}
        if balanced:
            for el in elements:
                states = {b[el] for b in balanced if el in b}
                ox_states_used[el] = sorted(states)

        fractions = compute_fractions(stoich)
        cation_species = a_species + b_species
        cation_set = set(A_SITE_POOL) | set(B_SITE_POOL)
        cation_fracs = compute_role_fractions(stoich, cation_set)
        anion_fracs = {"O": 1.0}

        return FormulaRecord(
            formula=text,
            canonical_formula=canonical,
            system=self.name,
            elements=elements,
            element_count=len(elements),
            cation_species=cation_species,
            anion_species=["O"],
            stoichiometry=stoich,
            atomic_fractions=fractions,
            cation_fractions=cation_fracs,
            anion_fractions=anion_fracs,
            prototype=prototype,
            mixed_occupancy=mixed,
            application_domains=self.application_domains,
            oxidation_states=ox_states_used,
            charge_balanced=charge_ok,
            warnings=warnings,
            extra={
                "a_site_species": a_species,
                "b_site_species": b_species,
                "total_a": total_a,
                "total_b": total_b,
                "oxygen_count": o_count,
                "ordered": not mixed and len(b_species) >= 2,
            },
        )

    def is_valid(self, formula: str) -> bool:
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
        min_a: int = 1,
        max_a: int = 3,
        min_b: int = 2,
        max_b: int = 4,
        max_results: int = 5000,
    ) -> list[str]:
        """Generate candidate double-perovskite formulas A2BB'O6.

        Generates ordered compositions in the A2BB'O6 framework. Mixed-
        occupancy variants can be derived from these ordered endpoints.

        Parameters
        ----------
        min_a, max_a
            Range for the number of A-site species.
        min_b, max_b
            Range for the number of B-site species.
        max_results
            Cap on the number of returned formulas.
        """
        results: list[str] = []
        seen: set[str] = set()

        for n_a in range(min_a, max_a + 1):
            for a_combo in itertools.combinations(A_SITE_POOL, n_a):
                a_species = list(a_combo)

                for n_b in range(min_b, max_b + 1):
                    for b_combo in itertools.combinations(B_SITE_POOL, n_b):
                        b_species = list(b_combo)

                        stoichs = _generate_perovskite_stoichiometries(
                            a_species, b_species
                        )

                        for stoich in stoichs:
                            formula = _stoich_to_formula(stoich)
                            try:
                                canonical = canonicalize_search_space_formula(formula)
                            except FormulaError:
                                continue
                            if canonical not in seen:
                                seen.add(canonical)
                                results.append(canonical)
                                if len(results) >= max_results:
                                    return results

        return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _by_z(symbol: str) -> int:
    from pymatgen.core.periodic_table import Element
    return Element(symbol).Z


def _expand_for_pymatgen(formula: str) -> str:
    """Expand mixed-occupancy notation for pymatgen compatibility.

    ``(La,Sr)2FeMoO6`` -> ``LaSrFeMoO6`` (equal sharing assumed).
    ``(La0.5Sr0.5)2FeMoO6`` -> ``LaFeMoO6`` (fractional amounts).
    """
    result = formula
    for match in re.finditer(r"\(([^)]+)\)(\d*\.?\d*)", result):
        inner = match.group(1)
        multiplier_str = match.group(2) or "1"
        multiplier = float(multiplier_str)

        parts = re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", inner)
        replacement_parts = []
        for elem, frac_str in parts:
            frac = float(frac_str) if frac_str else 1.0
            amount = frac * multiplier
            amount_int = int(round(amount))
            if amount_int <= 0:
                amount_int = 1
            if amount_int == 1:
                replacement_parts.append(elem)
            else:
                replacement_parts.append(f"{elem}{amount_int}")
        result = result.replace(match.group(0), "".join(replacement_parts), 1)
    return result


def _stoich_to_formula(stoich: dict[str, int]) -> str:
    parts = []
    for el in sorted(stoich, key=_by_z):
        count = stoich[el]
        if count == 1:
            parts.append(el)
        else:
            parts.append(f"{el}{count}")
    return "".join(parts)


def _generate_perovskite_stoichiometries(
    a_species: list[str],
    b_species: list[str],
) -> list[dict[str, int]]:
    """Generate charge-balanced A2BB'O6 stoichiometries.

    The double perovskite framework has A_total=2, B_total=2, O=6.
    For multiple A-site species, the A count is distributed equally.
    For multiple B-site species, the B count is distributed equally.
    Only even distributions are used (e.g. 2 A species -> 1 each).
    """
    results: list[dict[str, int]] = []
    n_a = len(a_species)
    n_b = len(b_species)

    if n_b < 2:
        return []

    a_total = 2
    b_total = 2

    if a_total % n_a != 0:
        return []
    if b_total % n_b != 0:
        return []

    a_each = a_total // n_a
    b_each = b_total // n_b

    stoich: dict[str, int] = {}
    for a_el in a_species:
        stoich[a_el] = a_each
    for b_el in b_species:
        stoich[b_el] = b_each
    stoich["O"] = 6

    balanced = check_charge_balance(stoich, _OXIDATION_STATES)
    if balanced:
        from .base import _reduce_stoichiometry
        reduced = _reduce_stoichiometry(stoich)
        results.append(reduced)

    return results


def _classify_prototype(
    stoich: dict[str, int],
    a_species: list[str],
    b_species: list[str],
    total_a: int,
    total_b: int,
    o_count: int,
) -> str:
    """Classify the double perovskite prototype variant."""
    n_a = len(a_species)
    n_b = len(b_species)

    if n_a == 1 and n_b == 2 and total_a == 2 and total_b == 2 and o_count == 6:
        return "A2BB'O6"
    if n_a > 1 and n_b == 2 and total_b == 2 and o_count == 6:
        return "(A,A')2BB'O6"
    if n_a == 1 and n_b > 2 and total_a == 2 and o_count == 6:
        return "A2(B,B')(B',B'')O6"
    if n_a > 1 and n_b > 2 and o_count == 6:
        return "(A,A')2(B,B')(B',B'')O6"
    if o_count == 3 and total_a == 1 and total_b == 1:
        return "ABO3"
    return "substituted"


def _analyze_mixed_sites(formula: str) -> str:
    """Extract a human-readable description of mixed-occupancy sites."""
    parts = []
    pattern = re.compile(r"\(([^)]+)\)")
    for match in pattern.finditer(formula):
        parts.append(match.group(1))
    return ", ".join(parts) if parts else ""
