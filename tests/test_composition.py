"""Tests for the generic composition module."""

from __future__ import annotations

import pytest

from formula_utils import FormulaError
from composition import (
    CompositionSpace,
    FormulaRecord,
    check_charge_balance,
    expand_mixed_occupancy,
    generate_charge_balanced_stoichiometries,
    has_mixed_occupancy,
)
from novelty_filter import (
    NoveltyAssessment,
    assess_novelty,
    filter_novel_only,
    novelty_summary,
)


# ---------------------------------------------------------------------------
# Fixtures: minimal composition-space configs for testing
# ---------------------------------------------------------------------------

_BINARY_OXIDE_CFG = {
    "cations": ["Mg", "Ca", "Ti"],
    "anions": ["O"],
    "oxidation_states": {
        "Mg": [2], "Ca": [2], "Ti": [3, 4],
        "O": [-2],
    },
    "constraints": {"max_stoich": 6, "max_total_atoms": 10},
}

_OXYCHALC_CFG = {
    "cations": ["La", "Ti"],
    "anions": ["O", "S"],
    "oxidation_states": {
        "La": [3], "Ti": [3, 4],
        "O": [-2], "S": [-2],
    },
    "constraints": {
        "required_anions": ["O"],
        "max_stoich": 6,
        "max_total_atoms": 12,
    },
}

_MPX_CFG = {
    "cations": ["Fe", "Ti"],
    "anions": ["P", "S"],
    "oxidation_states": {
        "Fe": [2, 3], "Ti": [3, 4],
        "P": [3, 4, 5],
        "S": [-2],
    },
    "constraints": {"max_stoich": 6, "max_total_atoms": 12},
}


# ---------------------------------------------------------------------------
# CompositionSpace construction
# ---------------------------------------------------------------------------

class TestCompositionSpaceInit:
    def test_basic_properties(self):
        space = CompositionSpace("test", _BINARY_OXIDE_CFG)
        assert space.name == "test"
        assert set(space.cations) == {"Mg", "Ca", "Ti"}
        assert set(space.anions) == {"O"}
        assert space.allowed_elements == {"Mg", "Ca", "Ti", "O"}

    def test_oxidation_states(self):
        space = CompositionSpace("test", _BINARY_OXIDE_CFG)
        assert space.oxidation_states["Ti"] == [3, 4]
        assert space.oxidation_states["O"] == [-2]

    def test_empty_config(self):
        space = CompositionSpace("empty", {"cations": [], "anions": [], "oxidation_states": {}})
        assert space.cations == []
        assert space.anions == []


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

class TestParse:
    def test_valid_binary_oxide(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        rec = space.parse("TiO2")
        assert rec.formula == "TiO2"
        assert rec.canonical_formula == "TiO2"
        assert set(rec.cation_species) == {"Ti"}
        assert rec.anion_species == ["O"]
        assert rec.charge_balanced is True

    def test_valid_multi_cation(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        rec = space.parse("CaTiO3")
        assert set(rec.cation_species) == {"Ca", "Ti"}
        assert rec.charge_balanced is True

    def test_element_outside_pool(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        with pytest.raises(FormulaError, match="outside"):
            space.parse("Fe2O3")

    def test_no_cation(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        with pytest.raises(FormulaError, match="cation"):
            space.parse("O2")

    def test_no_anion(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        with pytest.raises(FormulaError, match="anion"):
            space.parse("TiCa")

    def test_required_anion_enforced(self):
        space = CompositionSpace("oxychalc", _OXYCHALC_CFG)
        # O is required; a formula with only S should fail
        with pytest.raises(FormulaError, match="Required anion"):
            space.parse("La2S3")

    def test_required_anion_present(self):
        space = CompositionSpace("oxychalc", _OXYCHALC_CFG)
        rec = space.parse("LaTiO2S")
        assert rec.charge_balanced is True
        assert set(rec.anion_species) == {"O", "S"}

    def test_empty_formula(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        with pytest.raises(FormulaError):
            space.parse("")

    def test_is_valid_true(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        assert space.is_valid("TiO2") is True

    def test_is_valid_false(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        assert space.is_valid("Fe2O3") is False
        assert space.is_valid("") is False


# ---------------------------------------------------------------------------
# Charge balance
# ---------------------------------------------------------------------------

class TestChargeBalance:
    def test_balanced(self):
        ox = {"Ti": [3, 4], "O": [-2]}
        result = check_charge_balance({"Ti": 1, "O": 2}, ox)
        assert len(result) > 0
        assert result[0]["Ti"] == 4
        assert result[0]["O"] == -2

    def test_unbalanced(self):
        ox = {"Mg": [2], "O": [-2]}
        result = check_charge_balance({"Mg": 1, "O": 3}, ox)
        assert len(result) == 0

    def test_multiple_solutions(self):
        ox = {"Fe": [2, 3], "O": [-2]}
        result = check_charge_balance({"Fe": 2, "O": 3}, ox)
        assert len(result) >= 1


# ---------------------------------------------------------------------------
# Stoichiometry generation
# ---------------------------------------------------------------------------

class TestGeneration:
    def test_generate_binary(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        candidates = space.generate(max_coeff=4, max_total_atoms=8, max_results=20)
        assert len(candidates) > 0
        assert all(space.is_valid(f) for f in candidates)

    def test_generate_with_required_anion(self):
        space = CompositionSpace("oxychalc", _OXYCHALC_CFG)
        candidates = space.generate(max_coeff=4, max_total_atoms=10, max_results=20)
        assert len(candidates) > 0
        for f in candidates:
            from formula_utils import get_elements
            assert "O" in get_elements(f)

    def test_generate_respects_max_results(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        candidates = space.generate(max_results=5)
        assert len(candidates) <= 5

    def test_generate_charge_balanced_stoichiometries(self):
        ox = {"Mg": [2], "O": [-2]}
        stoichs = generate_charge_balanced_stoichiometries(
            cations=["Mg"], anions=["O"], max_coeff=4, max_total_atoms=6,
            oxidation_states=ox,
        )
        assert len(stoichs) > 0
        # MgO should be among them
        assert any(s == {"Mg": 1, "O": 1} for s in stoichs)


# ---------------------------------------------------------------------------
# Formula helpers
# ---------------------------------------------------------------------------

class TestFormulaHelpers:
    def test_has_mixed_occupancy(self):
        assert has_mixed_occupancy("(La,Sr)2FeMoO6") is True
        assert has_mixed_occupancy("La2FeMoO6") is False

    def test_expand_mixed_occupancy(self):
        result = expand_mixed_occupancy("(La,Sr)TiO3")
        assert len(result) == 2
        assert "LaTiO3" in result
        assert "SrTiO3" in result

    def test_expand_no_parens(self):
        result = expand_mixed_occupancy("LaTiO3")
        assert result == ["LaTiO3"]


# ---------------------------------------------------------------------------
# FormulaRecord
# ---------------------------------------------------------------------------

class TestFormulaRecord:
    def test_to_dict(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        rec = space.parse("TiO2")
        d = rec.to_dict()
        assert d["formula"] == "TiO2"
        assert "canonical_formula" in d
        assert "charge_balanced" in d

    def test_fields(self):
        space = CompositionSpace("oxide", _BINARY_OXIDE_CFG)
        rec = space.parse("CaTiO3")
        assert rec.space_name == "oxide"
        assert rec.element_count == 3
        assert len(rec.elements) == 3


# ---------------------------------------------------------------------------
# Novelty filter (basic smoke tests)
# ---------------------------------------------------------------------------

class TestNoveltyFilter:
    def test_novelty_assessment_dataclass(self):
        a = NoveltyAssessment(
            formula="TiO2", canonical_formula="O2Ti",
            known_formula=False, known_element_set=True,
            known_stoichiometry=True, novelty_category="known_stoichiometry_type",
            databases_found_in=[],
        )
        d = a.to_dict()
        assert d["formula"] == "TiO2"
        assert d["novelty_category"] == "known_stoichiometry_type"

    def test_filter_novel_only(self):
        a1 = NoveltyAssessment("A", "A", True, True, True, "known_formula", ["cod"])
        a2 = NoveltyAssessment("B", "B", False, False, False, "fully_novel", [])
        result = filter_novel_only([a1, a2])
        assert len(result) == 1
        assert result[0].formula == "B"

    def test_novelty_summary(self):
        a1 = NoveltyAssessment("A", "A", True, True, True, "known_formula", ["cod"])
        a2 = NoveltyAssessment("B", "B", False, False, False, "fully_novel", [])
        a3 = NoveltyAssessment("C", "C", False, False, False, "fully_novel", [])
        summary = novelty_summary([a1, a2, a3])
        assert summary["known_formula"] == 1
        assert summary["fully_novel"] == 2
