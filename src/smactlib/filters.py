"""SMACT chemical filter predicates (evaluation of a given formula).

Predicates are thin wrappers around :func:`smact.screening.smact_validity`
(and our own oxidation-state availability helper). All return plain booleans
and never raise; invalid formulas simply evaluate to ``False``.
"""

from __future__ import annotations

import logging

from formula_utils import safe_normalize_formula
from .generator import oxidation_states

logger = logging.getLogger(__name__)


def oxidation_states_available(elements, oxidation_set: str = "icsd24") -> bool:
    """True if every element has at least one state in the chosen SMACT set."""
    for sym in elements:
        try:
            if not oxidation_states(sym, oxidation_set):
                return False
        except Exception:  # noqa: BLE001 - unknown element / no data -> not available
            return False
    return True


def _validity(
    formula: str,
    oxidation_set: str,
    use_en: bool,
    include_alloys: bool,
) -> bool:
    """Run ``smact.screening.smact_validity``, never raising."""
    from smact.screening import smact_validity

    canonical = safe_normalize_formula(formula)
    if canonical is None:
        return False
    try:
        return bool(
            smact_validity(
                canonical,
                use_pauling_test=use_en,
                include_alloys=include_alloys,
                oxidation_states_set=oxidation_set,
            )
        )
    except Exception as exc:  # noqa: BLE001 - SMACT can raise on exotic inputs
        logger.debug("smact_validity raised for %s: %s", formula, exc)
        return False


def smact_valid(
    formula: str,
    oxidation_set: str = "icsd24",
    use_en: bool = True,
    include_alloys: bool = True,
) -> bool:
    """Overall SMACT validity: neutral composition under the oxidation-state
    set *and* (if ``use_en``) the Pauling electronegativity criterion, with
    alloys optionally auto-accepted."""
    return _validity(formula, oxidation_set, use_en=use_en, include_alloys=include_alloys)


def charge_neutral(
    formula: str,
    oxidation_set: str = "icsd24",
    include_alloys: bool = True,
) -> bool:
    """Charge-neutrality under the oxidation-state set (Pauling test off)."""
    return _validity(formula, oxidation_set, use_en=False, include_alloys=include_alloys)


def electronegativity_valid(
    formula: str,
    oxidation_set: str = "icsd24",
    include_alloys: bool = True,
) -> bool:
    """Pauling electronegativity condition (interpretable only when the
    composition is also charge neutral)."""
    return _validity(formula, oxidation_set, use_en=True, include_alloys=include_alloys)
