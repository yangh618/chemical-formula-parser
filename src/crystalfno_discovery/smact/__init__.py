"""SMACT-based candidate-composition generation and filtering.

SMACT (>= 3) exposes its rules through :mod:`smact.screening` and the
:mod:`smact` top-level helpers (``neutral_ratios``, ``pauling_test``,
``element_dictionary``). We use those SMACT building blocks here rather than
re-implementing oxidation-state chemistry.

Submodules:

* :mod:`~.generator` - enumerate chemically plausible charge-neutral candidate
  compositions over an element universe;
* :mod:`~.filters`  - evaluate SMACT validity / charge-neutrality /
  electronegativity flags for a given formula.
"""

from __future__ import annotations

from .filters import (
    charge_neutral,
    electronegativity_valid,
    oxidation_states_available,
    smact_valid,
)
from .generator import generate_formulas

__all__ = [
    "generate_formulas",
    "smact_valid",
    "charge_neutral",
    "electronegativity_valid",
    "oxidation_states_available",
]
