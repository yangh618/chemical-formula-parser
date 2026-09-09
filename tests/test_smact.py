"""Tests for SMACT generation and filtering."""

from __future__ import annotations

from smactlib import filters as f
from smactlib.generator import generate_formulas, oxidation_states


def test_oxide_generation():
    out = generate_formulas(["Al", "O"], num_elements=(2,), max_stoich=6, oxidation_states_set="icsd24")
    assert "Al2O3" in out
    assert all("Al" in x and "O" in x for x in out)


def test_alloy_generation_with_include_alloys():
    out = generate_formulas(["Al", "Nb"], num_elements=(2,), max_stoich=4, oxidation_states_set="icsd24")
    assert out  # non-empty
    # Charge-neutral metal-only pair is impossible, so outputs rely on alloys.
    assert "NbAl" in out


def test_no_alloys_excludes_metal_pair():
    out = generate_formulas(["Al", "Nb"], num_elements=(2,), max_stoich=4,
                            oxidation_states_set="icsd24", include_alloys=False)
    assert out == []


def test_generator_dedupes_reduced_forms():
    # Al2Nb and NbAl2 (same) should appear once.
    out = generate_formulas(["Al", "Nb"], num_elements=(2,), max_stoich=6)
    assert len(out) == len(set(out))
    assert "NbAl" in out


def test_oxidation_states_lookup():
    assert oxidation_states("Al", "icsd24") == [3]
    assert -2 in oxidation_states("O", "icsd24")


def test_filters():
    assert f.smact_valid("NaCl", oxidation_set="icsd24") is True
    assert f.charge_neutral("Al2O3", oxidation_set="icsd24") is True
    assert f.electronegativity_valid("Al2O3", oxidation_set="icsd24") is True
    assert f.oxidation_states_available(["Al", "O"]) is True
    # Invalid element -> filter fails safely.
    assert f.smact_valid("FeX2", oxidation_set="icsd24") is False
    assert f.oxidation_states_available(["NotARealElem"]) is False
