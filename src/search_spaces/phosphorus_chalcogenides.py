"""Transition-metal phosphorus chalcogenide search space: M-P-X.

Chemical space
--------------
* Transition metals (M): Ti, V, Cr, Mn, Fe, Co, Ni, Cu
* Phosphorus (P): distinct sublattice
* Chalcogen (X): S, Se, Te

Supports multimetal compositions (M, M', M'' ...) and mixed chalcogen anions.
The parser recognises P as a distinct sublattice species, separate from the
chalcogen anions.
"""

from __future__ import annotations

import itertools
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
    generate_charge_balanced_stoichiometries,
    has_mixed_occupancy,
)

# ---------------------------------------------------------------------------
# Element pools
# ---------------------------------------------------------------------------

TRANSITION_METALS: list[str] = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu"]
PHOSPHORUS: list[str] = ["P"]
CHALCOGEN_POOL: list[str] = ["S", "Se", "Te"]
ANION_POOL: list[str] = PHOSPHORUS + CHALCOGEN_POOL

ALL_ALLOWED = set(TRANSITION_METALS + PHOSPHORUS + CHALCOGEN_POOL)

_APPLICATION_DOMAINS = [
    "2D materials",
    "van der Waals materials",
    "spintronics",
    "quantum materials",
    "magnetism",
    "optoelectronics",
    "catalysis",
    "energy conversion",
]

_OXIDATION_STATES = {
    "Ti": [3, 4], "V": [3, 4, 5], "Cr": [3, 6],
    "Mn": [2, 3, 4], "Fe": [2, 3], "Co": [2, 3],
    "Ni": [2, 3], "Cu": [1, 2],
    "P": [3, 4, 5],
    "S": [-2], "Se": [-2], "Te": [-2],
}


class PhosphorusChalcogenideSpace(SearchSpace):
    """Transition-metal phosphorus chalcogenide M-P-X search space."""

    @property
    def name(self) -> str:
        return "phosphorus_chalcogenides"

    @property
    def application_domains(self) -> list[str]:
        return list(_APPLICATION_DOMAINS)

    @property
    def allowed_elements(self) -> set[str]:
        return ALL_ALLOWED

    @property
    def transition_metal_pool(self) -> list[str]:
        return list(TRANSITION_METALS)

    @property
    def chalcogen_pool(self) -> list[str]:
        return list(CHALCOGEN_POOL)

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
                f"Elements outside phosphorus-chalcogenide pool: {sorted(outside)}"
            )

        tm_species = sorted(
            element_set & set(TRANSITION_METALS), key=lambda s: _by_z(s)
        )
        has_p = "P" in element_set
        chalcogen_species = sorted(
            element_set & set(CHALCOGEN_POOL), key=lambda s: _by_z(s)
        )

        if not tm_species:
            raise FormulaError("No transition-metal species found")
        if not has_p:
            raise FormulaError("Phosphorus (P) is required in this search space")
        if not chalcogen_species:
            raise FormulaError("No chalcogen species (S, Se, Te) found")

        stoich = _reduced_amounts(text)
        canonical = canonicalize_search_space_formula(text)

        cation_species = list(tm_species)
        anion_species = ["P"] + chalcogen_species

        warnings: list[str] = []
        n_tm = len(tm_species)
        if n_tm == 1:
            warnings.append("Single transition metal")
        elif n_tm == 2:
            warnings.append("Bimetallic composition")
        elif n_tm >= 3:
            warnings.append(f"Multimetallic ({n_tm} transition metals)")

        if len(chalcogen_species) > 1:
            warnings.append(
                f"Mixed chalcogen: {', '.join(chalcogen_species)}"
            )

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
        cation_fracs = compute_role_fractions(stoich, set(TRANSITION_METALS))
        anion_set = set(ANION_POOL)
        anion_fracs = compute_role_fractions(stoich, anion_set)

        p_stoich = stoich.get("P", 0)
        prototype = _infer_prototype(stoich, tm_species, chalcogen_species)

        return FormulaRecord(
            formula=text,
            canonical_formula=canonical,
            system=self.name,
            elements=elements,
            element_count=len(elements),
            cation_species=cation_species,
            anion_species=anion_species,
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
                "phosphorus_count": p_stoich,
                "n_transition_metals": n_tm,
                "chalcogen_species": chalcogen_species,
                "mixed_chalcogen": len(chalcogen_species) > 1,
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
        min_metals: int = 1,
        max_metals: int = 4,
        max_coeff: int = 8,
        max_total_atoms: int = 20,
        max_results: int = 5000,
    ) -> list[str]:
        """Generate candidate M-P-X formulas.

        Parameters
        ----------
        min_metals, max_metals
            Range for the number of transition-metal species.
        max_coeff
            Maximum stoichiometric coefficient for any element.
        max_total_atoms
            Maximum total atoms in the reduced formula.
        max_results
            Cap on the number of returned formulas.
        """
        results: list[str] = []
        seen: set[str] = set()

        for n_metal in range(min_metals, max_metals + 1):
            for metal_combo in itertools.combinations(TRANSITION_METALS, n_metal):
                metals = list(metal_combo)

                chalcogen_combos = []
                for n_chalc in range(1, len(CHALCOGEN_POOL) + 1):
                    for chalc in itertools.combinations(CHALCOGEN_POOL, n_chalc):
                        chalcogen_combos.append(["P"] + list(chalc))

                for anions in chalcogen_combos:
                    stoichs = generate_charge_balanced_stoichiometries(
                        cations=metals,
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
    tm_species: list[str],
    chalcogen_species: list[str],
) -> str:
    """Assign a prototype label based on stoichiometry pattern."""
    total_m = sum(stoich.get(m, 0) for m in tm_species)
    p_count = stoich.get("P", 0)
    total_chalc = sum(stoich.get(x, 0) for x in chalcogen_species)
    total_anion = p_count + total_chalc

    if total_m == 1 and p_count == 1 and total_chalc == 3:
        return "MPX3"
    if total_m == 2 and p_count == 1 and total_chalc == 4:
        return "M2PX4"
    if total_m == 1 and p_count == 1 and total_chalc == 2:
        return "MPX2"
    if total_m == 3 and p_count == 1 and total_chalc == 4:
        return "M3PX4"
    if total_m == 1 and p_count == 2 and total_chalc == 6:
        return "MP2X6"
    return "custom"
