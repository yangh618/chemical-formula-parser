"""Multicomponent oxychalcogenide search space: A-B-O-X.

Chemical space
--------------
* Cations (A/B): rare-earth + 3d transition metals
* Anions: O + chalcogen (S, Se, Te)

The parser treats this as a *search space*, not a fixed stoichiometry. Generated
compositions may contain only a subset of the allowed species. Database absence
is recorded as novelty, never used as a rejection criterion.
"""

from __future__ import annotations

import itertools
from typing import Any

from formula_utils import FormulaError, get_elements

from .base import (
    COMMON_OXIDATION_STATES,
    FormulaRecord,
    SearchSpace,
    canonicalize_search_space_formula,
    check_charge_balance,
    compute_fractions,
    compute_role_fractions,
    flag_unusual_oxidation_states,
    generate_charge_balanced_stoichiometries,
    has_mixed_occupancy,
)

# ---------------------------------------------------------------------------
# Element pools
# ---------------------------------------------------------------------------

RARE_EARTHS: list[str] = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd"]
TRANSITION_METALS: list[str] = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
CATION_POOL: list[str] = RARE_EARTHS + TRANSITION_METALS
ANION_POOL: list[str] = ["O", "S", "Se", "Te"]
CHALCOGEN_POOL: list[str] = ["S", "Se", "Te"]

ALL_ALLOWED = set(CATION_POOL + ANION_POOL)

_APPLICATION_DOMAINS = [
    "optoelectronics",
    "photocatalysis",
    "thermoelectrics",
    "electronics",
    "catalysis",
    "magnetic materials",
    "energy materials",
]

_OXIDATION_STATES = {
    "La": [3], "Ce": [3, 4], "Pr": [3, 4], "Nd": [3],
    "Sm": [2, 3], "Eu": [2, 3], "Gd": [3],
    "Ti": [3, 4], "V": [3, 4, 5], "Cr": [3, 6],
    "Mn": [2, 3, 4], "Fe": [2, 3], "Co": [2, 3], "Ni": [2, 3],
    "O": [-2], "S": [-2], "Se": [-2], "Te": [-2],
}


class OxychalcogenideSpace(SearchSpace):
    """Multicomponent A-B-O-X oxychalcogenide search space."""

    @property
    def name(self) -> str:
        return "oxychalcogenides"

    @property
    def application_domains(self) -> list[str]:
        return list(_APPLICATION_DOMAINS)

    @property
    def allowed_elements(self) -> set[str]:
        return ALL_ALLOWED

    @property
    def cation_pool(self) -> list[str]:
        return list(CATION_POOL)

    @property
    def anion_pool(self) -> list[str]:
        return list(ANION_POOL)

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    def parse(self, formula: str) -> FormulaRecord:
        if not formula or not formula.strip():
            raise FormulaError("Empty formula")

        text = formula.strip()
        mixed = has_mixed_occupancy(text)

        elements = get_elements(text)
        element_set = set(elements)

        outside = element_set - ALL_ALLOWED
        if outside:
            raise FormulaError(
                f"Elements outside oxychalcogenide pool: {sorted(outside)}"
            )

        cations = sorted(element_set & set(CATION_POOL), key=lambda s: _by_z(s))
        anions = sorted(element_set & set(ANION_POOL), key=lambda s: _by_z(s))

        if not cations:
            raise FormulaError("No cation species (rare-earth or transition metal) found")
        if not anions:
            raise FormulaError("No anion species (O, S, Se, Te) found")

        has_oxygen = "O" in element_set
        has_chalcogen = bool(element_set & set(CHALCOGEN_POOL))

        if not has_oxygen:
            raise FormulaError(
                "Oxychalcogenide system requires at least oxygen; "
                "chalcogen-only compositions belong to a different system"
            )

        from formula_utils import _reduced_amounts
        stoich = _reduced_amounts(text)
        canonical = canonicalize_search_space_formula(text)

        warnings: list[str] = []
        if has_chalcogen:
            warnings.append("Mixed anion: oxychalcogenide (O + chalcogen)")
        else:
            warnings.append("Oxide-only: no chalcogen present")

        ox_warnings = flag_unusual_oxidation_states(stoich, _OXIDATION_STATES)
        warnings.extend(ox_warnings)

        balanced = check_charge_balance(stoich, _OXIDATION_STATES)
        charge_ok = len(balanced) > 0

        if not charge_ok:
            warnings.append("Unusual oxidation-state combination flagged for review")

        ox_states_used = {}
        if balanced:
            for el in elements:
                states = {b[el] for b in balanced if el in b}
                ox_states_used[el] = sorted(states)

        fractions = compute_fractions(stoich)
        cation_fracs = compute_role_fractions(stoich, set(CATION_POOL))
        anion_fracs = compute_role_fractions(stoich, set(ANION_POOL))

        prototype = _infer_prototype(stoich, cations, anions)

        return FormulaRecord(
            formula=text,
            canonical_formula=canonical,
            system=self.name,
            elements=elements,
            element_count=len(elements),
            cation_species=cations,
            anion_species=anions,
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
        min_cations: int = 1,
        max_cations: int = 4,
        include_chalcogen: bool = True,
        chalcogen_only: bool = False,
        max_coeff: int = 8,
        max_total_atoms: int = 20,
        max_results: int = 5000,
    ) -> list[str]:
        """Generate candidate oxychalcogenide formulas.

        Parameters
        ----------
        min_cations, max_cations
            Range for the number of cation species.
        include_chalcogen
            If ``True``, generate compositions with at least one chalcogen.
        chalcogen_only
            If ``True``, require at least one chalcogen (no oxide-only).
        max_coeff
            Maximum stoichiometric coefficient for any element.
        max_total_atoms
            Maximum total atoms in the reduced formula.
        max_results
            Cap on the number of returned formulas.
        """
        results: list[str] = []
        seen: set[str] = set()

        for n_cat in range(min_cations, max_cations + 1):
            for cat_combo in itertools.combinations(CATION_POOL, n_cat):
                cations = list(cat_combo)

                anion_combos = self._anion_combos(include_chalcogen, chalcogen_only)

                for anions in anion_combos:
                    stoichs = generate_charge_balanced_stoichiometries(
                        cations=cations,
                        anions=anions,
                        max_coeff=max_coeff,
                        max_total_atoms=max_total_atoms,
                        oxidation_states=_OXIDATION_STATES,
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

    def _anion_combos(
        self, include_chalcogen: bool, chalcogen_only: bool
    ) -> list[list[str]]:
        combos: list[list[str]] = []
        if chalcogen_only or include_chalcogen:
            for n_chalc in range(1, len(CHALCOGEN_POOL) + 1):
                for chalc in itertools.combinations(CHALCOGEN_POOL, n_chalc):
                    combos.append(["O"] + list(chalc))
        if not chalcogen_only:
            combos.append(["O"])
        return combos


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _by_z(symbol: str) -> int:
    from pymatgen.core.periodic_table import Element
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


def _infer_prototype(
    stoich: dict[str, int],
    cations: list[str],
    anions: list[str],
) -> str:
    """Assign a loose prototype label based on stoichiometry pattern."""
    total_cat = sum(stoich.get(c, 0) for c in cations)
    total_an = sum(stoich.get(a, 0) for a in anions)

    if total_an == 0:
        return "unknown"

    ratio = total_cat / total_an

    if len(cations) == 1 and total_cat == 2 and total_an == 3:
        return "A2X3"
    if len(cations) == 2 and total_cat == 2 and total_an == 4:
        return "A2B2X4"
    if len(cations) == 1 and total_cat == 1 and total_an == 2:
        return "AX2"
    if len(cations) == 1 and total_cat == 1 and total_an == 1:
        return "AX"
    if len(cations) == 2 and total_cat == 4 and total_an == 7:
        return "A2B4X7"
    if abs(ratio - 1.0) < 0.15:
        return "AX-like"
    if abs(ratio - 0.5) < 0.1:
        return "A2X-like"
    return "custom"
