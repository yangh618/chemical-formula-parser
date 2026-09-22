"""Tests for formula parsing / canonicalization / vectors."""

from __future__ import annotations

import numpy as np
import pytest

from formula_utils import (
    FormulaError,
    composition_distance,
    element_set,
    element_set_distance,
    flatten_formula,
    get_anonymous_formula,
    get_elemental_fraction_vector,
    get_elements,
    get_num_elements,
    get_stoichiometric_vector,
    is_valid_formula,
    normalize_formula,
    safe_flatten_formula,
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


# ---------------------------------------------------------------------------
# flatten_formula: parenthesis-free output
# ---------------------------------------------------------------------------

def test_flatten_removes_parentheses():
    assert flatten_formula("La4(SO)3") == "La4S3O3"
    assert flatten_formula("Ti(PS3)2") == "TiP2S6"
    assert flatten_formula("La(NiO2)3") == "LaNi3O6"
    assert "(" not in flatten_formula("Ag(AuF4)2")
    assert flatten_formula("Ag(AuF4)2") == "AgAu2F8"


def test_flatten_preserves_conventional_element_order():
    """Cations come before anions, matching normalize_formula's ordering."""
    assert flatten_formula("La2O3") == "La2O3"      # not O3La2
    assert flatten_formula("TiO2") == "TiO2"        # not O2Ti
    assert flatten_formula("Ca2La2O5") == "Ca2La2O5"


def test_flatten_is_notation_independent():
    """Any spelling of the same composition flattens identically."""
    for equivalent in ["La4(SO)3", "La4S3O3", "La4O3S3", "O3S3La4", "S3La4O3"]:
        assert flatten_formula(equivalent) == "La4S3O3"
    for equivalent in ["Ti(PS3)2", "TiP2S6", "S6TiP2", "P2TiS6"]:
        assert flatten_formula(equivalent) == "TiP2S6"


def test_flatten_round_trips_to_same_canonical_form():
    """Flattening must not change the underlying composition."""
    for f in ["La4(SO)3", "Ti(PS3)2", "La(NiO2)3", "La6Se(SO3)2", "Ag(AuF4)2", "Fe2O3"]:
        assert normalize_formula(flatten_formula(f)) == normalize_formula(f)


def test_flatten_is_idempotent():
    for f in ["La4(SO)3", "Ti(PS3)2", "La2O3"]:
        once = flatten_formula(f)
        assert flatten_formula(once) == once


def test_flatten_omits_unit_counts():
    assert flatten_formula("FeO") == "FeO"
    assert "1" not in flatten_formula("FeO")
    assert "1" not in flatten_formula("La4(SO)3")
    # pymatgen orders by electronegativity, so S precedes O here
    assert flatten_formula("LaTiO2S") == normalize_formula("LaTiO2S") == "LaTiSO2"


def test_flatten_invalid_input():
    for bad in ["", "XxYy", "Fe3+", "(La2O3", "La2O3)"]:
        with pytest.raises(FormulaError):
            flatten_formula(bad)
        assert safe_flatten_formula(bad) is None
