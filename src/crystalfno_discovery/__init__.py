"""CrystalFNO external validation and novel-composition discovery pipeline.

This package:

* collects crystal structures from external databases (COD, OQMD, JARVIS,
  optionally AFLOW) into a unified representation for the **external
  validation** of a pretrained CrystalFNO model;
* systematically generates chemically plausible candidate formulas with
  SMACT and screens them for **database-unseen** novelty (Task B).

The package does **not** train CrystalFNO or any other machine-learning
model. CrystalFNO is treated as a frozen model connected through the abstract
:class:`~crystalfno_discovery.validation.dataset.CrystalFNOInterface`.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = ["__version__"]
