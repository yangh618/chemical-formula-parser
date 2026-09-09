"""Tests for formula parsing / canonicalization / vectors."""

from __future__ import annotations

import numpy as np
import pytest

from formula_utils import (
    FormulaError,
    composition_distance,
    element_set,
    element_set_distance,
    get_anonymous_formula,
    get_elemental_fraction_vector,
    get_elements,
    get_num_elements,
    get_stoichiometric_vector,
    is_valid_formula,
    normalize_formula,
    safe_normalize_formula,
)


def test_equivalent_reduced_formulas():
    forms = ["Al2CrNb", "CrAl2Nb", "NbCrAl2", "Al0.6666667Cr0.3333333Nb0.3333333", "Al0.6667Cr0.3333Nb0.3333"]
    reduced = {normalize_formula(f) for f in forms}
    assert len(reduced) == 1
    assert reduced == {"NbAl2Cr"}


def test_anonymous_same_for_same_stoich():
    assert get_anonymous_formula("Al2CrNb") == "ABC2"
    assert get_anonymous_formula("Ti2VTa") == "ABC2"
    assert get_anonymous_formula("Al2CrNb") == get_anonymous_formula("Ti2VTa")


def test_elements_and_counts():
    els = get_elements("Al2CrNb")
    assert els == ["Al", "Cr", "Nb"]  # sorted by atomic number
    assert get_num_elements("Al2CrNb") == 3
    vec = get_stoichiometric_vector("Al2CrNb")
    assert vec.shape == (118,)
    assert vec.sum() == 4
    assert vec[Element_index("Al")] == 2 and vec[Element_index("Nb")] == 1


def Element_index(sym: str) -> int:
    from pymatgen.core.periodic_table import Element

    return Element(sym).Z - 1


def test_fraction_vector_sums_to_one():
    v = get_elemental_fraction_vector("Fe2O3")
    assert v.shape == (118,)
    assert np.isclose(v.sum(), 1.0)


def test_invalid_formulas_fail_safely():
    for bad in ["", "XxYy", "Al2Zz", "NotAReal", "..", "Fe3+"]:
        assert not is_valid_formula(bad)
        with pytest.raises(FormulaError):
            normalize_formula(bad)
        assert safe_normalize_formula(bad) is None


def test_distances():
    assert composition_distance("Al2CrNb", "CrAl2Nb") == 0.0
    assert composition_distance("Fe2O3", "Fe2O3") == 0.0
    assert composition_distance("Fe2O3", "FeO3") > 0.0
    assert element_set_distance("Fe2O3", "Fe2O3") == 0.0
    assert element_set_distance("Fe2O3", "Al2O3") == pytest.approx(1 - 1 / 3)  # share O only
    assert element_set("Fe2O3") == {"Fe", "O"}
