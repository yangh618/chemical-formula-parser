"""Candidate-composition generator built on SMACT.

Design
------
Candidate compositions are produced by enumerating, for every element
combination of the requested size, all charge-neutral integer ratios allowed
by the chosen SMACT oxidation-state set (``smact.neutral_ratios``), with an
optional Pauling electronegativity test (``smact.screening.pauling_test``).
Results are canonicalised with our formula utilities, so ``AlCrNb2`` and
``Nb2CrAl`` collapse to one row.

This enumerates *chemically plausible* (charge-neutral / SMACT-valid)
candidate compositions. It is **not** a statement about thermodynamic
stability or experimental realizability (scientific safeguards, section 29).
"""

from __future__ import annotations

import itertools
import logging
from typing import Iterable

from ..formula_utils import normalize_formula

logger = logging.getLogger(__name__)

# Attribute name on smact.Element for each named oxidation-state set.
OXIDATION_SET_ATTR = {
    "smact14": "oxidation_states_smact14",
    "icsd16": "oxidation_states_icsd16",
    "icsd24": "oxidation_states_icsd24",
    "pymatgen_sp": "oxidation_states_sp",
    "wiki": "oxidation_states_wiki",
}

# Guard against pathological enumeration blow-ups.
_MAX_OX_PRODUCT = 200_000


def oxidation_states(element_symbol: str, oxidation_set: str = "icsd24") -> list[int]:
    """Return the SMACT oxidation states for *element_symbol* in *oxidation_set*."""
    from smact import Element

    attr = OXIDATION_SET_ATTR.get(str(oxidation_set).lower())
    if attr is None:
        raise ValueError(
            f"Unknown oxidation-state set '{oxidation_set}'. Choose from {sorted(OXIDATION_SET_ATTR)} "
            "or pass a custom value."
        )
    states = getattr(Element(element_symbol), attr)
    return list(states or [])


def _combo_candidates(
    element_symbols: tuple[str, ...],
    max_stoich: int,
    oxidation_set: str,
    use_en: bool,
    include_alloys: bool,
) -> set[str]:
    """Return canonical reduced formulas for one element combination.

    Two regimes are handled:

    * **All-metallic** combos (when ``include_alloys``): SMACT treats any
      all-metal composition as an *alloy* (no integer charge-balancing), so we
      enumerate integer stoichiometry vectors ``(1..max_stoich)^k`` and keep
      their reduced (gcd-collapsed) formulas.
    * **Combos containing a non-metal** (or ``include_alloys=False``): a proper
      charge-balanced enumeration over the SMACT oxidation-state set is used,
      with the optional Pauling electronegativity test.
    """
    from smact import Element, metals, neutral_ratios
    from smact.screening import pauling_test

    result: set[str] = set()
    all_metals = all(sym in metals for sym in element_symbols)

    def _add(vector) -> None:
        parts = []
        for sym, coef in zip(element_symbols, vector, strict=True):
            parts.append(f"{sym}{int(coef) if int(coef) != 1 else ''}")
        try:
            result.add(normalize_formula("".join(parts)))
        except Exception:  # noqa: BLE001 - skip unparseable
            return

    if all_metals and include_alloys:
        for vector in itertools.product(range(1, max_stoich + 1), repeat=len(element_symbols)):
            _add(vector)
        return result

    ox_lists = [oxidation_states(sym, oxidation_set) for sym in element_symbols]
    if any(not states for states in ox_lists):
        logger.debug("No oxidation states in '%s' for %s; skipping combo", oxidation_set, element_symbols)
        return result

    ox_product = 1
    for states in ox_lists:
        ox_product *= len(states)
    if ox_product > _MAX_OX_PRODUCT:
        logger.warning("Oxidation-state product %d too large for %s; skipping.", ox_product, element_symbols)
        return result

    enegs = [Element(sym).pauling_eneg for sym in element_symbols]
    for ox_assignment in itertools.product(*ox_lists):
        ratios = neutral_ratios(list(ox_assignment), threshold=max_stoich)
        for ratio in ratios:
            if use_en and not pauling_test(list(ox_assignment), enegs):
                continue
            _add(ratio)
    return result


def generate_formulas(
    elements: Iterable[str],
    num_elements: Iterable[int] = (2, 3, 4),
    max_stoich: int = 8,
    oxidation_states_set: str = "icsd24",
    electronegativity_filter: bool = False,
    include_alloys: bool = True,
    **_: object,
) -> list[str]:
    """Generate unique canonical reduced formulas over *elements*.

    Parameters
    ----------
    elements : iterable of element symbols used as the combinatorial universe.
    num_elements : sizes of element combinations to consider (e.g. ``(2, 3)``).
    max_stoich : maximum stoichiometric coefficient of any element.
    oxidation_states_set : SMACT oxidation-state set name (default ``icsd24``).
    electronegativity_filter : if True, additionally apply the Pauling test
        (meaningful only for charge-balanced, non-alloy combos).
    include_alloys : if True, treat all-metal combinations as valid alloys and
        enumerate integer stoichiometry vectors for them.

    Returns
    -------
    list[str] : sorted list of canonical reduced formulas (no duplicates).
    """
    elements = sorted({el.strip() for el in elements if el and el.strip()})
    sizes = sorted({int(n) for n in num_elements if int(n) >= 1})
    if len(elements) < 2:
        raise ValueError("Provide at least two elements to generate multi-element formulas.")
    max_stoich = int(max_stoich)
    if max_stoich < 1:
        raise ValueError("max_stoich must be >= 1")

    found: set[str] = set()
    for size in sizes:
        for combo in itertools.combinations(elements, size):
            candidates = _combo_candidates(
                combo, max_stoich, oxidation_states_set, electronegativity_filter, include_alloys)
            found.update(candidates)
    logger.info("Generated %d unique SMACT candidate formulas over %d elements.",
                len(found), len(elements))
    return sorted(found)
