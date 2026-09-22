"""Structure-predictor interface + external benchmark/dataset helpers.

This package defines a small abstract interface that downstream structure-
prediction models can implement.  It does not train or ship any model.
"""

from __future__ import annotations

from .dataset import StructurePredictorInterface, load_formulas

__all__ = ["StructurePredictorInterface", "load_formulas"]
