"""Tests for diversity selection and stratified generation."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pytest

from composition import CompositionSpace
from diversity import build_features, diversity_report, select_diverse
from formula_utils import get_elements


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_OXIDE_CFG = {
    "cations": ["Mg", "Ca", "Ti", "Fe"],
    "anions": ["O"],
    "oxidation_states": {
        "Mg": [2], "Ca": [2], "Ti": [3, 4], "Fe": [2, 3], "O": [-2],
    },
    "constraints": {"max_stoich": 4, "max_total_atoms": 8, "max_results": 500},
}

# Wider pool (~1000 formulas, ~380 element sets) for diversity assertions that
# need more candidates than the small fixture can supply.
_OXYSULFIDE_CFG = {
    "cations": ["Mg", "Ca", "Ti", "Fe", "Mn", "Co", "Ni", "Zn"],
    "anions": ["O", "S"],
    "oxidation_states": {
        "Mg": [2], "Ca": [2], "Ti": [3, 4], "Fe": [2, 3],
        "Mn": [2, 3, 4], "Co": [2, 3], "Ni": [2, 3], "Zn": [2],
        "O": [-2], "S": [-2],
    },
    "constraints": {"max_stoich": 4, "max_total_atoms": 8, "max_results": 2000},
}


def _space(cfg: dict | None = None) -> CompositionSpace:
    return CompositionSpace("test_oxides", cfg or _OXIDE_CFG)


def _big_pool(per_combo_limit: int = 4, max_results: int = 2000) -> list[str]:
    """A pool large enough for meaningful diversity assertions."""
    return _space(_OXYSULFIDE_CFG).generate(
        max_results=max_results, per_combo_limit=per_combo_limit
    )


def _element_set(formula: str) -> frozenset[str]:
    return frozenset(get_elements(formula))


# ---------------------------------------------------------------------------
# Stratified generation (per_combo_limit)
# ---------------------------------------------------------------------------

def test_per_combo_limit_caps_each_element_combination():
    """With a cap of 1, no element combination may contribute twice."""
    candidates = _space().generate(max_results=500, per_combo_limit=1)
    counts = Counter(_element_set(f) for f in candidates)
    assert counts, "expected some candidates"
    assert max(counts.values()) == 1


def test_per_combo_limit_spreads_coverage():
    """Stratification yields more distinct element sets for the same count."""
    space = _space(_OXYSULFIDE_CFG)
    n = 40
    naive = space.generate(max_results=n)
    stratified = space.generate(max_results=n, per_combo_limit=2)

    naive_sets = len({_element_set(f) for f in naive})
    strat_sets = len({_element_set(f) for f in stratified})
    assert strat_sets > naive_sets


def test_per_combo_limit_none_means_unlimited():
    """Omitting the cap reproduces the plain enumeration order."""
    space = _space()
    with_none = space.generate(max_results=100, per_combo_limit=None)
    without = space.generate(max_results=100)
    assert with_none == without


def test_max_cation_species_constraint_respected():
    cfg = dict(_OXIDE_CFG)
    cfg["constraints"] = {**_OXIDE_CFG["constraints"], "max_cation_species": 2}
    space = _space(cfg)

    assert space.constraints["max_cation_species"] == 2
    for f in space.generate(max_results=200):
        cations = _element_set(f) - {"O"}
        assert len(cations) <= 2


def test_max_cation_species_defaults_to_full_pool():
    space = _space()
    assert space.constraints["max_cation_species"] == 0
    # The full four-cation pool is reachable when uncapped.
    sizes = {len(_element_set(f) - {"O"}) for f in space.generate(max_results=500)}
    assert max(sizes) > 2


def test_generate_kwargs_override_config():
    cfg = dict(_OXIDE_CFG)
    cfg["constraints"] = {**_OXIDE_CFG["constraints"], "max_cation_species": 1}
    space = _space(cfg)

    capped = space.generate(max_results=200)
    assert all(len(_element_set(f) - {"O"}) == 1 for f in capped)

    overridden = space.generate(max_results=200, max_cations=3)
    assert any(len(_element_set(f) - {"O"}) > 1 for f in overridden)


# ---------------------------------------------------------------------------
# Feature building
# ---------------------------------------------------------------------------

def test_build_features_shapes_and_valid_filtering():
    formulas = ["Fe2O3", "TiO2", "Xx2Yy3", "MgO"]
    valid, frac, pres = build_features(formulas)

    assert valid == [0, 1, 3]           # the invalid entry is dropped
    assert frac.shape == (3, 118)
    assert pres.shape == (3, 118)
    assert np.allclose(frac.sum(axis=1), 1.0)
    assert pres[0].sum() == 2           # Fe2O3 spans Fe and O
    assert pres[2].sum() == 2           # MgO spans Mg and O


def test_build_features_presence_is_binary():
    _, _, pres = build_features(["Fe2O3"])
    assert set(pres.ravel().tolist()) <= {0.0, 1.0}


# ---------------------------------------------------------------------------
# Diversity selection
# ---------------------------------------------------------------------------

def test_select_diverse_returns_requested_count():
    pool = _big_pool()
    selected = select_diverse(pool, 25)
    assert len(selected) == 25


def test_select_diverse_is_subset_without_duplicates():
    pool = _big_pool()
    selected = select_diverse(pool, 40)
    assert len(selected) == len(set(selected))
    assert set(selected) <= set(pool)


def test_select_diverse_clamps_to_pool_size():
    pool = ["Fe2O3", "TiO2", "MgO"]
    assert len(select_diverse(pool, 100)) == 3


def test_select_diverse_deduplicates_input():
    pool = ["Fe2O3", "Fe2O3", "TiO2", "MgO", "MgO"]
    selected = select_diverse(pool, 10)
    assert len(selected) == 3
    assert len(set(selected)) == 3


def test_select_diverse_handles_empty_and_invalid():
    assert select_diverse([], 10) == []
    assert select_diverse(["Xx2Yy3", "Zz4Qq"], 10) == []


def test_select_diverse_is_deterministic():
    pool = _big_pool()
    assert select_diverse(pool, 30) == select_diverse(pool, 30)


def test_select_diverse_beats_enumeration_head_on_min_distance():
    """Max-min selection must raise the closest-pair distance."""
    pool = _big_pool()
    n = 30
    naive = pool[:n]
    diverse = select_diverse(pool, n)

    assert diversity_report(diverse).min_pairwise_distance > \
           diversity_report(naive).min_pairwise_distance


def test_select_diverse_covers_more_element_sets_than_head():
    pool = _big_pool()
    n = 30
    naive_sets = len({_element_set(f) for f in pool[:n]})
    diverse_sets = len({_element_set(f) for f in select_diverse(pool, n)})
    assert diverse_sets > naive_sets


def test_weights_change_the_selection():
    """Emphasising element identity gives a different subset than composition."""
    pool = _big_pool()
    n = 30
    by_comp = select_diverse(pool, n, w_composition=1.0, w_element_set=0.0)
    by_species = select_diverse(pool, n, w_composition=0.0, w_element_set=1.0)

    assert len(by_species) == n
    # With no composition term, picks are spread over distinct element sets
    # until the pool's element sets run out.
    available = len({_element_set(f) for f in pool})
    assert len({_element_set(f) for f in by_species}) == min(n, available)
    assert by_comp != by_species


def test_zero_weights_do_not_crash():
    pool = _big_pool(max_results=50)
    selected = select_diverse(pool, 5, w_composition=0.0, w_element_set=0.0)
    assert len(selected) == 5


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def test_diversity_report_basic_fields():
    report = diversity_report(["Fe2O3", "TiO2", "MgO", "CaTiO3"]).to_dict()
    assert report["n_formulas"] == 4
    assert report["distinct_element_sets"] == 4
    assert report["elements_covered"] == 5          # Fe, Ti, O, Mg, Ca
    assert report["min_pairwise_distance"] > 0.0
    assert report["mean_pairwise_distance"] >= report["min_pairwise_distance"]


def test_diversity_report_identical_set_has_zero_distance():
    report = diversity_report(["Fe2O3", "Fe2O3", "Fe2O3"])
    assert report.n_formulas == 1
    assert report.min_pairwise_distance == 0.0
    assert report.mean_pairwise_distance == 0.0


def test_diversity_report_empty():
    report = diversity_report([])
    assert report.n_formulas == 0
    assert report.distinct_element_sets == 0


def test_diversity_report_samples_large_sets():
    pool = _big_pool()
    assert len(pool) > 50
    report = diversity_report(pool, max_pairs=50)
    assert report.sampled is True
    assert report.n_formulas == len(set(pool))


def test_diversity_report_most_common_element_set():
    report = diversity_report(["Fe2O3", "FeO", "Fe3O4", "TiO2"])
    assert set(report.most_common_element_set) == {"Fe", "O"}
    assert report.most_common_element_set_count == 3
