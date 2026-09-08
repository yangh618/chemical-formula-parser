"""Abstract CrystalFNO interface and formula-loading convenience."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ..formula_utils import FormulaError, normalize_formula


class CrystalFNOInterface(ABC):
    """Abstract interface a pretrained (frozen) CrystalFNO model must satisfy.

    Subclass this with your own model wrapper. This package does **not**
    implement, train or modify the CrystalFNO neural network.
    """

    @abstractmethod
    def predict(self, formula: str) -> Any:
        """Predict a crystal structure (or structural target) for *formula*.

        ``formula`` is a canonical reduced composition string such as
        ``"NbAl2Cr"``. Implementations should return whatever structure-like
        object the model emits (e.g. a pymatgen Structure).
        """
        raise NotImplementedError

    def predict_many(self, formulas) -> list[Any]:
        """Convenience loop over :meth:`predict` (override for batching)."""
        return [self.predict(f) for f in formulas]


def load_formulas(path: str | Path, as_dataframe: bool = False) -> Any:
    """Load candidate formulas from a novelty/benchmark CSV.

    ``path`` points at a CSV with a ``formula`` column (e.g.
    ``data/formulas/novel_formulas.csv``). Returns a pandas DataFrame by
    default so downstream code can select on novelty flags; pass
    ``as_dataframe=False`` to get a plain list of canonical formula strings.

    Rows whose ``formula`` cannot be canonicalised are dropped.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Formula file not found: {path}")
    df = pd.read_csv(path)
    if "formula" not in df.columns:
        raise ValueError(f"'{path}' has no 'formula' column (columns: {list(df.columns)})")

    canonical: list[str] = []
    valid_mask: list[bool] = []
    for raw in df["formula"]:
        try:
            normalized = normalize_formula(str(raw))
        except FormulaError:
            normalized = ""
        canonical.append(normalized)
        valid_mask.append(bool(normalized))

    out = df.copy()
    out["formula"] = canonical
    out = out[valid_mask].reset_index(drop=True)
    if as_dataframe:
        return out
    return out["formula"].tolist()


def to_formula_list(df: Optional[pd.DataFrame]) -> list[str]:
    """Extract the canonical ``formula`` column as a list of strings."""
    if df is None or df.empty:
        return []
    return df["formula"].astype(str).tolist()
