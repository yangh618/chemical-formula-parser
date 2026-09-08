"""Tests for novelty screening and formula loading."""

from __future__ import annotations

import pandas as pd
import pytest

from crystalfno_discovery.novelty import (
    KnownIndex,
    annotate,
    build_known_index,
    canonicalize_formula,
)
from crystalfno_discovery.validation import load_formulas


def _known_df():
    return pd.DataFrame(
        [
            {"database": "COD", "formula": "NbAl2Cr", "anonymous_formula": "ABC2", "elements": "Al Cr Nb"},
            {"database": "OQMD", "formula": "Fe2O3", "anonymous_formula": "A2B3", "elements": "Fe O"},
            {"database": "MP", "formula": "TiO2", "anonymous_formula": "AB2", "elements": "Ti O"},
        ]
    )


def test_known_formula_rejected():
    index = build_known_index(_known_df())
    row = annotate(canonicalize_formula("CrAl2Nb"), index)
    assert row["novel_exact_formula"] is False
    assert row["found_in_cod"] is True
    assert row["found_in_mp"] is False


def test_unknown_formula_retained():
    index = build_known_index(_known_df())
    row = annotate(canonicalize_formula("ZrCoW2"), index)
    assert row["novel_exact_formula"] is True
    assert row["found_in_any"] is False
    assert row["novel_element_combination"] is True
    # ABC2 already appears (NbAl2Cr) -> anonymous not novel despite exact novelty.
    assert row["novel_anonymous_formula"] is False


def test_three_levels_distinct():
    index = build_known_index(_known_df())
    # Level 1 false (present), level 2 false, level 3 false.
    exact = annotate(canonicalize_formula("Al2CrNb"), index)
    assert (exact["novel_exact_formula"], exact["novel_element_combination"], exact["novel_anonymous_formula"]) == (
        False, False, False)
    # Fully unseen composition.
    unseen = annotate(canonicalize_formula("ZrMoW3"), index)
    assert unseen["novel_exact_formula"] is True
    assert unseen["novel_element_combination"] is True
    assert unseen["novel_anonymous_formula"] is True


def test_known_index_mp_membership():
    index = build_known_index(_known_df())
    assert index.formula_in("NbAl2Cr", database="mp") is False  # NbAl2Cr is COD, not MP
    assert index.formula_in("TiO2", database="mp") is True
    assert index.formula_in("Fe2O3") is True


def test_load_formulas(tmp_path):
    path = tmp_path / "n.csv"
    path.write_text("formula,novel_exact_formula\nAl2CrNb,False\nZrCoW2,True\nBogusZZ,False\n")
    formulas = load_formulas(path, as_dataframe=False)
    assert formulas == ["NbAl2Cr", "ZrCoW2"]  # BogusZZ dropped safely
    df = load_formulas(path, as_dataframe=True)
    assert list(df["formula"]) == ["NbAl2Cr", "ZrCoW2"]


def test_empty_known_index_treats_all_as_unseen():
    index = KnownIndex()
    index.finalize()
    row = annotate(canonicalize_formula("Al2CrNb"), index)
    assert row["novel_exact_formula"] is True


@pytest.mark.parametrize("bad", ["", "BogusZZ", "Xx", "NotARealElem2"])
def test_bad_inputs_safe(bad):
    from crystalfno_discovery.novelty import safe_canonicalize_formula

    assert safe_canonicalize_formula(bad) is None


def test_whitespace_formula_accepted():
    from crystalfno_discovery.novelty import safe_canonicalize_formula

    assert safe_canonicalize_formula("Al2CrNb ") is not None
    assert safe_canonicalize_formula("Al") is not None  # single element is a valid formula
