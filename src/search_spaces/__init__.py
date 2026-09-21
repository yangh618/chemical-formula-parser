"""Chemical search spaces for CrystalFNO novel-formula discovery.

Each search space defines a chemically meaningful region of composition space
and provides parsing, validation, and generation capabilities.

Available search spaces:

* :class:`OxychalcogenideSpace` -- multicomponent A-B-O-X oxychalcogenides
* :class:`PhosphorusChalcogenideSpace` -- transition-metal phosphorus chalcogenides
* :class:`DoublePerovskiteSpace` -- double perovskite oxides A2BB'O6

Use :func:`get_space` to look up a space by name, or :func:`all_spaces` to
instantiate every registered space.
"""

from .base import (
    COMMON_OXIDATION_STATES,
    FormulaRecord,
    SearchSpace,
    canonicalize_search_space_formula,
    check_charge_balance,
    compute_fractions,
    expand_mixed_occupancy,
    flag_unusual_oxidation_states,
    generate_charge_balanced_stoichiometries,
    has_mixed_occupancy,
    safe_canonicalize,
)
from .double_perovskites import DoublePerovskiteSpace
from .novelty_filter import (
    NoveltyAssessment,
    assess_novelty,
    assess_novelty_batch,
    filter_by_category,
    filter_novel_only,
    novelty_summary,
)
from .oxychalcogenides import OxychalcogenideSpace
from .phosphorus_chalcogenides import PhosphorusChalcogenideSpace

_REGISTRY: dict[str, type[SearchSpace]] = {
    "oxychalcogenides": OxychalcogenideSpace,
    "phosphorus_chalcogenides": PhosphorusChalcogenideSpace,
    "double_perovskites": DoublePerovskiteSpace,
}


def get_space(name: str) -> SearchSpace:
    """Instantiate a search space by name.

    Raises :class:`KeyError` if the name is not registered.
    """
    cls = _REGISTRY.get(name)
    if cls is None:
        raise KeyError(
            f"Unknown search space {name!r}. "
            f"Available: {sorted(_REGISTRY)}"
        )
    return cls()


def all_spaces() -> list[SearchSpace]:
    """Instantiate all registered search spaces."""
    return [cls() for cls in _REGISTRY.values()]


def register_space(name: str, cls: type[SearchSpace]) -> None:
    """Register a new search space class under *name*."""
    _REGISTRY[name] = cls


def list_spaces() -> list[str]:
    """Return the names of all registered search spaces."""
    return sorted(_REGISTRY)


__all__ = [
    "COMMON_OXIDATION_STATES",
    "DoublePerovskiteSpace",
    "FormulaRecord",
    "NoveltyAssessment",
    "OxychalcogenideSpace",
    "PhosphorusChalcogenideSpace",
    "SearchSpace",
    "all_spaces",
    "assess_novelty",
    "assess_novelty_batch",
    "canonicalize_search_space_formula",
    "check_charge_balance",
    "compute_fractions",
    "expand_mixed_occupancy",
    "filter_by_category",
    "filter_novel_only",
    "flag_unusual_oxidation_states",
    "generate_charge_balanced_stoichiometries",
    "get_space",
    "has_mixed_occupancy",
    "list_spaces",
    "novelty_summary",
    "register_space",
    "safe_canonicalize",
]
