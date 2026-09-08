"""CrystalFNO interface + external benchmark/dataset helpers.

This package only defines the *interface* to a frozen CrystalFNO model. It does
not train CrystalFNO or any other machine-learning model (section 19).
"""

from __future__ import annotations

from .dataset import CrystalFNOInterface, load_formulas

__all__ = ["CrystalFNOInterface", "load_formulas"]
