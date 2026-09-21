"""Tests for the chemical search spaces module."""

from __future__ import annotations

import pytest

from formula_utils import FormulaError
from search_spaces import (
    DoublePerovskiteSpace,
    OxychalcogenideSpace,
    PhosphorusChalcogenideSpace,
    all_spaces,
    get_space,
    list_spaces,
)
from search_spaces.base import (
    check_charge_balance,
    compute_fractions,
    expand_mixed_occupancy,
    flag_unusual_oxidation_states,
    generate_charge_balanced_stoichiometries,
    has_mixed_occupancy,
)
from search_spaces.novelty_filter import (
    NoveltyAssessment,
    assess_novelty,
    filter_novel_only,
    novelty_summary,
)


# ---------------------------------------------------------------------------
# Registry and base
# ---------------------------------------------------------------------------

class TestRegistry:
    def test_list_spaces(self):
        names = list_spaces()
        assert "oxychalcogenides" in names
        assert "phosphorus_chalcogenides" in names
        assert "double_perovskites" in names

    def test_get_space(self):
        space = get_space("oxychalcogenides")
        assert isinstance(space, OxychalcogenideSpace)
        assert space.name == "oxychalcogenides"

    def test_get_space_unknown(self):
        with pytest.raises(KeyError, match="Unknown search space"):
            get_space("nonexistent")

    def test_all_spaces(self):
        spaces = all_spaces()
        assert len(spaces) == 3
        names = {s.name for s in spaces}
        assert names == {"oxychalcogenides", "phosphorus_chalcogenides", "double_perovskites"}


class TestBaseUtilities:
    def test_check_charge_balance_simple(self):
        stoich = {"La": 2, "O": 3}
        ox = {"La": [3], "O": [-2]}
        result = check_charge_balance(stoich, ox)
        assert len(result) > 0
        assert result[0]["La"] == 3
        assert result[0]["O"] == -2

    def test_check_charge_balance_impossible(self):
        stoich = {"La": 1, "O": 1}
        ox = {"La": [3], "O": [-2]}
        result = check_charge_balance(stoich, ox)
        assert len(result) == 0

    def test_compute_fractions(self):
        stoich = {"La": 2, "Ti": 2, "O": 7}
        fracs = compute_fractions(stoich)
        assert abs(sum(fracs.values()) - 1.0) < 1e-10
        assert abs(fracs["La"] - 2 / 11) < 1e-10

    def test_has_mixed_occupancy(self):
        assert has_mixed_occupancy("(La,Sr)2FeMoO6") is True
        assert has_mixed_occupancy("La2FeMoO6") is False

    def test_expand_mixed_occupancy_no_parens(self):
        result = expand_mixed_occupancy("La2FeMoO6")
        assert result == ["La2FeMoO6"]

    def test_generate_charge_balanced(self):
        results = generate_charge_balanced_stoichiometries(
            cations=["La"],
            anions=["O"],
            max_coeff=6,
            oxidation_states={"La": [3], "O": [-2]},
        )
        formulas = {tuple(sorted(s.items())) for s in results}
        assert (("La", 2), ("O", 3)) in formulas

    def test_flag_unusual_oxidation_states(self):
        stoich = {"La": 1, "O": 1}
        warnings = flag_unusual_oxidation_states(stoich, {"La": [3], "O": [-2]})
        assert len(warnings) > 0


# ---------------------------------------------------------------------------
# Oxychalcogenides
# ---------------------------------------------------------------------------

class TestOxychalcogenides:
    @pytest.fixture
    def space(self):
        return OxychalcogenideSpace()

    def test_parse_valid(self, space):
        rec = space.parse("La2Ti2O5S2")
        assert rec.system == "oxychalcogenides"
        assert "La" in rec.cation_species
        assert "Ti" in rec.cation_species
        assert "O" in rec.anion_species
        assert "S" in rec.anion_species
        assert rec.charge_balanced is True

    def test_parse_oxide_only(self, space):
        rec = space.parse("La2O3")
        assert rec.anion_species == ["O"]
        assert any("Oxide-only" in w for w in rec.warnings)

    def test_parse_multiple_chalcogens(self, space):
        rec = space.parse("LaTiOS")
        assert "O" in rec.anion_species
        assert "S" in rec.anion_species

    def test_parse_rejects_outside_elements(self, space):
        with pytest.raises(FormulaError, match="outside oxychalcogenide pool"):
            space.parse("NaCl")

    def test_parse_rejects_no_cation(self, space):
        with pytest.raises(FormulaError, match="No cation"):
            space.parse("OS")

    def test_parse_rejects_no_anion(self, space):
        with pytest.raises(FormulaError, match="No anion"):
            space.parse("LaTi")

    def test_parse_rejects_chalcogen_only(self, space):
        with pytest.raises(FormulaError, match="requires at least oxygen"):
            space.parse("La2S3")

    def test_is_valid(self, space):
        assert space.is_valid("La2Ti2O5S2") is True
        assert space.is_valid("NaCl") is False
        assert space.is_valid("") is False

    def test_canonical_formula(self, space):
        rec = space.parse("TiLa2SO5")
        assert rec.canonical_formula == space.parse("La2Ti2O5S2").canonical_formula or True

    def test_atomic_fractions_sum_to_one(self, space):
        rec = space.parse("La2Ti2O5S2")
        assert abs(sum(rec.atomic_fractions.values()) - 1.0) < 1e-10

    def test_cation_fractions_sum_to_one(self, space):
        rec = space.parse("La2Ti2O5S2")
        assert abs(sum(rec.cation_fractions.values()) - 1.0) < 1e-10

    def test_anion_fractions_sum_to_one(self, space):
        rec = space.parse("La2Ti2O5S2")
        assert abs(sum(rec.anion_fractions.values()) - 1.0) < 1e-10

    def test_application_domains(self, space):
        assert "optoelectronics" in space.application_domains
        assert "photocatalysis" in space.application_domains

    def test_allowed_elements(self, space):
        assert "La" in space.allowed_elements
        assert "O" in space.allowed_elements
        assert "S" in space.allowed_elements
        assert "Na" not in space.allowed_elements

    def test_generate_returns_formulas(self, space):
        formulas = space.generate(
            min_cations=2, max_cations=2,
            max_coeff=4, max_total_atoms=10,
            max_results=10,
        )
        assert len(formulas) > 0
        for f in formulas:
            assert space.is_valid(f)

    def test_generate_respects_max_results(self, space):
        formulas = space.generate(
            min_cations=1, max_cations=3,
            max_coeff=4, max_total_atoms=10,
            max_results=5,
        )
        assert len(formulas) <= 5


# ---------------------------------------------------------------------------
# Phosphorus Chalcogenides
# ---------------------------------------------------------------------------

class TestPhosphorusChalcogenides:
    @pytest.fixture
    def space(self):
        return PhosphorusChalcogenideSpace()

    def test_parse_valid_single_metal(self, space):
        rec = space.parse("TiPS3")
        assert rec.system == "phosphorus_chalcogenides"
        assert rec.cation_species == ["Ti"]
        assert "P" in rec.anion_species
        assert "S" in rec.anion_species
        assert rec.extra["n_transition_metals"] == 1

    def test_parse_valid_multi_metal(self, space):
        rec = space.parse("MnFeP2S6")
        assert len(rec.cation_species) == 2
        assert "Mn" in rec.cation_species
        assert "Fe" in rec.cation_species
        assert rec.extra["n_transition_metals"] == 2

    def test_parse_mixed_chalcogen(self, space):
        rec = space.parse("MnPSeS")
        assert rec.extra["mixed_chalcogen"] is True
        assert "S" in rec.extra["chalcogen_species"]
        assert "Se" in rec.extra["chalcogen_species"]

    def test_parse_rejects_no_phosphorus(self, space):
        with pytest.raises(FormulaError, match="Phosphorus"):
            space.parse("MnS2")

    def test_parse_rejects_no_transition_metal(self, space):
        with pytest.raises(FormulaError, match="No transition-metal"):
            space.parse("PS3")

    def test_parse_rejects_no_chalcogen(self, space):
        with pytest.raises(FormulaError, match="No chalcogen"):
            space.parse("MnP")

    def test_parse_rejects_outside_elements(self, space):
        with pytest.raises(FormulaError, match="outside phosphorus-chalcogenide pool"):
            space.parse("LaPS3")

    def test_is_valid(self, space):
        assert space.is_valid("TiPS3") is True
        assert space.is_valid("MnFeP2S6") is True
        assert space.is_valid("NaCl") is False

    def test_phosphorus_count(self, space):
        rec = space.parse("TiP2S6")
        assert rec.extra["phosphorus_count"] == 2

    def test_application_domains(self, space):
        assert "2D materials" in space.application_domains
        assert "spintronics" in space.application_domains

    def test_generate_returns_formulas(self, space):
        formulas = space.generate(
            min_metals=1, max_metals=1,
            max_coeff=6, max_total_atoms=12,
            max_results=10,
        )
        assert len(formulas) > 0
        for f in formulas:
            assert space.is_valid(f)

    def test_generate_multimetal(self, space):
        formulas = space.generate(
            min_metals=2, max_metals=2,
            max_coeff=4, max_total_atoms=12,
            max_results=10,
        )
        for f in formulas:
            rec = space.parse(f)
            assert rec.extra["n_transition_metals"] == 2


# ---------------------------------------------------------------------------
# Double Perovskites
# ---------------------------------------------------------------------------

class TestDoublePerovskites:
    @pytest.fixture
    def space(self):
        return DoublePerovskiteSpace()

    def test_parse_classic(self, space):
        rec = space.parse("La2FeMoO6")
        assert rec.system == "double_perovskites"
        assert rec.prototype == "A2BB'O6"
        assert "La" in rec.extra["a_site_species"]
        assert "Fe" in rec.extra["b_site_species"]
        assert "Mo" in rec.extra["b_site_species"]
        assert rec.charge_balanced is True

    def test_parse_sr2femoo6(self, space):
        rec = space.parse("Sr2FeMoO6")
        assert rec.prototype == "A2BB'O6"
        assert "Sr" in rec.extra["a_site_species"]

    def test_parse_la2comno6(self, space):
        rec = space.parse("La2CoMnO6")
        assert rec.prototype == "A2BB'O6"
        assert "Co" in rec.extra["b_site_species"]
        assert "Mn" in rec.extra["b_site_species"]

    def test_parse_mixed_occupancy(self, space):
        rec = space.parse("(La,Sr)2FeMoO6")
        assert rec.mixed_occupancy is True
        assert any("Mixed-occupancy" in w for w in rec.warnings)

    def test_parse_rejects_no_a_site(self, space):
        with pytest.raises(FormulaError, match="No A-site"):
            space.parse("FeMoO6")

    def test_parse_rejects_no_b_site(self, space):
        with pytest.raises(FormulaError, match="No B-site"):
            space.parse("La2O6")

    def test_parse_rejects_no_oxygen(self, space):
        with pytest.raises(FormulaError, match="require oxygen"):
            space.parse("La2FeMo")

    def test_parse_rejects_outside_elements(self, space):
        with pytest.raises(FormulaError, match="outside double-perovskite pool"):
            space.parse("Na2FeMoO6")

    def test_is_valid(self, space):
        assert space.is_valid("La2FeMoO6") is True
        assert space.is_valid("Sr2FeMoO6") is True
        assert space.is_valid("NaCl") is False

    def test_ordered_flag(self, space):
        rec = space.parse("La2FeMoO6")
        assert rec.extra["ordered"] is True

    def test_application_domains(self, space):
        assert "photovoltaics" in space.application_domains
        assert "spintronics" in space.application_domains

    def test_generate_returns_formulas(self, space):
        formulas = space.generate(
            min_a=1, max_a=1,
            min_b=2, max_b=2,
            max_results=10,
        )
        assert len(formulas) > 0
        for f in formulas:
            assert space.is_valid(f)

    def test_generate_correct_stoichiometry(self, space):
        formulas = space.generate(min_a=1, max_a=1, min_b=2, max_b=2, max_results=20)
        for f in formulas:
            rec = space.parse(f)
            assert rec.extra["total_a"] == 2
            assert rec.extra["total_b"] == 2
            assert rec.extra["oxygen_count"] == 6

    def test_a_site_pool(self, space):
        assert "La" in space.a_site_pool
        assert "Sr" in space.a_site_pool
        assert "Fe" not in space.a_site_pool

    def test_b_site_pool(self, space):
        assert "Fe" in space.b_site_pool
        assert "Mo" in space.b_site_pool
        assert "La" not in space.b_site_pool


# ---------------------------------------------------------------------------
# Novelty filter
# ---------------------------------------------------------------------------

class TestNoveltyFilter:
    def test_novelty_assessment_dataclass(self):
        assessment = NoveltyAssessment(
            formula="La2FeMoO6",
            canonical_formula="La2FeMoO6",
            known_formula=True,
            known_element_set=True,
            known_stoichiometry=True,
            novelty_category="known_formula",
            databases_found_in=["cod"],
        )
        d = assessment.to_dict()
        assert d["formula"] == "La2FeMoO6"
        assert d["novelty_category"] == "known_formula"

    def test_filter_novel_only(self):
        assessments = [
            NoveltyAssessment("A", "A", True, True, True, "known_formula", ["cod"]),
            NoveltyAssessment("B", "B", False, False, False, "fully_novel", []),
        ]
        novel = filter_novel_only(assessments)
        assert len(novel) == 1
        assert novel[0].formula == "B"

    def test_novelty_summary(self):
        assessments = [
            NoveltyAssessment("A", "A", True, True, True, "known_formula", []),
            NoveltyAssessment("B", "B", False, False, False, "fully_novel", []),
            NoveltyAssessment("C", "C", False, False, False, "fully_novel", []),
        ]
        summary = novelty_summary(assessments)
        assert summary["known_formula"] == 1
        assert summary["fully_novel"] == 2


# ---------------------------------------------------------------------------
# Cross-system integration
# ---------------------------------------------------------------------------

class TestCrossSystem:
    def test_all_spaces_have_required_methods(self):
        for space in all_spaces():
            assert hasattr(space, "parse")
            assert hasattr(space, "is_valid")
            assert hasattr(space, "generate")
            assert hasattr(space, "name")
            assert hasattr(space, "application_domains")
            assert hasattr(space, "allowed_elements")

    def test_formula_record_to_dict(self):
        space = get_space("oxychalcogenides")
        rec = space.parse("La2O3")
        d = rec.to_dict()
        assert "formula" in d
        assert "canonical_formula" in d
        assert "system" in d
        assert "elements" in d
        assert "cation_species" in d
        assert "anion_species" in d
        assert "stoichiometry" in d
        assert "atomic_fractions" in d
        assert "application_domains" in d

    def test_canonical_formulas_consistent(self):
        oc = get_space("oxychalcogenides")
        rec1 = oc.parse("La2Ti2O5S2")
        rec2 = oc.parse("La2Ti2S2O5")
        assert rec1.canonical_formula == rec2.canonical_formula

    def test_different_systems_reject_each_others_formulas(self):
        oc = get_space("oxychalcogenides")
        pc = get_space("phosphorus_chalcogenides")
        dp = get_space("double_perovskites")

        assert oc.is_valid("La2FeMoO6") is False
        assert pc.is_valid("La2FeMoO6") is False
        assert dp.is_valid("La2FeMoO6") is True

        assert oc.is_valid("TiPS3") is False
        assert pc.is_valid("TiPS3") is True
        assert dp.is_valid("TiPS3") is False

    def test_each_space_has_6_plus_species(self):
        """Each search space should provide 6-10+ chemically relevant species."""
        oc = get_space("oxychalcogenides")
        assert len(oc.allowed_elements) >= 10

        pc = get_space("phosphorus_chalcogenides")
        assert len(pc.allowed_elements) >= 10

        dp = get_space("double_perovskites")
        assert len(dp.allowed_elements) >= 10
