"""Tests for novelty screening and formula loading."""

from __future__ import annotations

import pandas as pd
import pytest

from novelty import (
    KnownIndex,
    annotate,
    build_known_index,
    canonicalize_formula,
)
from validation import load_formulas


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
    from novelty import safe_canonicalize_formula

    assert safe_canonicalize_formula(bad) is None


def test_whitespace_formula_accepted():
    from novelty import safe_canonicalize_formula

    assert safe_canonicalize_formula("Al2CrNb ") is not None
    assert safe_canonicalize_formula("Al") is not None  # single element is a valid formula


# ---------------------------------------------------------------------------
# generate_novel_candidates: adaptive limit + optional validator hook
# ---------------------------------------------------------------------------

class _FakeSpace:
    """Minimal stand-in for CompositionSpace: returns a fixed candidate list."""

    def __init__(self, candidates):
        self._candidates = list(candidates)
        self.limits: list[int] = []

    def generate(self, max_results=None, **_):
        self.limits.append(max_results)
        return self._candidates[: max_results or len(self._candidates)]


def test_generate_novel_candidates_scales_limit_until_target_met():
    from novelty_filter import generate_novel_candidates

    candidates = [f"C{i}" for i in range(50)]
    known = {f"C{i}" for i in range(40)}  # only the last 10 are novel
    space = _FakeSpace(candidates)

    result = generate_novel_candidates(space, 10, known, base_limit=20, growth=2.0, verbose=False)

    assert result.target_met is True
    assert result.novel_found == 10
    assert result.known_filtered == 40
    assert space.limits == [20, 40, 80]  # grew until the target was met


def test_generate_novel_candidates_validator_drops_rejected():
    from novelty_filter import generate_novel_candidates

    candidates = ["NaCl", "Al2O3", "FeX2", "TiO2"]
    space = _FakeSpace(candidates)
    validator = lambda f: f != "FeX2"  # noqa: E731

    result = generate_novel_candidates(space, 3, set(), validator=validator, verbose=False)

    assert result.rejected_by_validator == 1
    assert result.novel == ["NaCl", "Al2O3", "TiO2"]
    assert result.total_generated == 4


def test_generate_novel_candidates_without_validator_rejects_nothing():
    from novelty_filter import generate_novel_candidates

    space = _FakeSpace(["NaCl", "Al2O3"])
    result = generate_novel_candidates(space, 2, set(), verbose=False)

    assert result.rejected_by_validator == 0
    assert result.novel == ["NaCl", "Al2O3"]


def test_generate_novel_candidates_reports_exhausted_space():
    from novelty_filter import generate_novel_candidates

    space = _FakeSpace(["NaCl", "Al2O3"])
    result = generate_novel_candidates(space, 100, set(), base_limit=50, verbose=False)

    assert result.target_met is False
    assert result.space_exhausted is True
    assert result.novel_found == 2


def test_generate_novel_candidates_passes_through_generate_kwargs():
    from novelty_filter import generate_novel_candidates

    class _Recording(_FakeSpace):
        def generate(self, max_results=None, **kwargs):
            self.kwargs = kwargs
            return super().generate(max_results=max_results, **kwargs)

    space = _Recording(["NaCl"])
    generate_novel_candidates(space, 1, set(), per_combo_limit=4, max_cations=3, verbose=False)

    assert space.kwargs == {"per_combo_limit": 4, "max_cations": 3}


# ---------------------------------------------------------------------------
# Merged unique-formula CSV (datasets/known_formulas.csv) support
# ---------------------------------------------------------------------------

def _merged_fixture(tmp_path):
    """Three tiny per-database CSVs with deliberate cross-database overlap."""
    (tmp_path / "cod_formulas.csv").write_text(
        "formula\nAc(AgAs)2\nNbAl2Cr\nTiO2\n", encoding="utf-8")
    (tmp_path / "mp_formulas.csv").write_text(
        "formula\nNbAl2Cr\nFe2O3\n", encoding="utf-8")
    (tmp_path / "oqmd_formulas.csv").write_text(
        "formula\nNbAl2Cr\nZrO2\n", encoding="utf-8")


def test_merge_database_csv_files_writes_each_unique_formula_once(tmp_path):
    import csv

    from novelty_filter import merge_database_csv_files

    _merged_fixture(tmp_path)
    result = merge_database_csv_files(datasets_dir=tmp_path)

    assert result.output_path == tmp_path / "known_formulas.csv"
    assert result.per_database == {"cod": 3, "mp": 2, "oqmd": 2}
    assert result.overlap_removed == 2          # NbAl2Cr counted once, not thrice
    assert result.total_unique == 5

    with open(result.output_path, newline="", encoding="utf-8") as f:
        rows = {r["formula"]: r["databases"] for r in csv.DictReader(f)}

    assert sorted(rows) == ["AcAg2As2", "Fe2O3", "NbAl2Cr", "TiO2", "ZrO2"]
    assert rows["AcAg2As2"] == "cod"            # polyatomic group flattened
    assert rows["NbAl2Cr"] == "cod;mp;oqmd"     # provenance preserved


def test_merge_database_csv_files_subset_and_custom_output(tmp_path):
    from novelty_filter import merge_database_csv_files

    _merged_fixture(tmp_path)
    out = tmp_path / "sub" / "cod_mp.csv"
    result = merge_database_csv_files(
        datasets_dir=tmp_path, databases=["cod", "mp"], output_path=out)

    assert set(result.per_database) == {"cod", "mp"}
    assert result.total_unique == 4
    assert out.exists()


def test_merge_database_csv_files_adds_derived_metadata(tmp_path):
    import csv

    from novelty_filter import load_known_index_from_csvs, merge_database_csv_files

    _merged_fixture(tmp_path)
    result = merge_database_csv_files(datasets_dir=tmp_path)
    with open(result.output_path, newline="", encoding="utf-8") as f:
        rows = {r["formula"]: r for r in csv.DictReader(f)}

    assert rows["NbAl2Cr"]["anonymous_formula"] == "ABC2"
    assert rows["NbAl2Cr"]["elements"] == "Al Cr Nb"
    assert rows["NbAl2Cr"]["num_elements"] == "3"

    # the derived columns feed the element/anonymous novelty levels
    index = load_known_index_from_csvs(csv_path=result.output_path)
    assert frozenset({"Al", "Cr", "Nb"}) in index.element_combos
    assert "ABC2" in index.anonymous


def test_merge_database_csv_files_no_metadata_skips_derivation(tmp_path):
    import csv

    from novelty_filter import merge_database_csv_files

    _merged_fixture(tmp_path)
    result = merge_database_csv_files(datasets_dir=tmp_path, with_metadata=False)
    with open(result.output_path, newline="", encoding="utf-8") as f:
        rows = {r["formula"]: r for r in csv.DictReader(f)}

    assert rows["NbAl2Cr"]["anonymous_formula"] == ""
    assert rows["NbAl2Cr"]["elements"] == ""
    assert rows["NbAl2Cr"]["databases"] == "cod;mp;oqmd"


def test_comparison_key_flattens_polyatomic_groups():
    from novelty_filter import comparison_key

    assert comparison_key("Ac(AgAs)2") == "AcAg2As2"
    assert comparison_key(" NbAl2Cr ") == "NbAl2Cr"   # already canonical
    assert comparison_key("") is None
    assert comparison_key("   ") is None


def test_load_known_formulas_set_prefers_merged_csv(tmp_path):
    from novelty_filter import load_known_formulas_set, merge_database_csv_files

    _merged_fixture(tmp_path)
    merge_database_csv_files(datasets_dir=tmp_path)
    with open(tmp_path / "cod_formulas.csv", "a", encoding="utf-8") as f:
        f.write("OnlyInFile\n")

    formulas = load_known_formulas_set(datasets_dir=tmp_path)
    assert {"AcAg2As2", "Fe2O3", "NbAl2Cr", "TiO2", "ZrO2"} <= formulas
    assert "OnlyInFile" not in formulas            # merged CSV is read

    per_db = load_known_formulas_set(datasets_dir=tmp_path, databases=["cod"])
    assert "OnlyInFile" in per_db                  # --db uses the per-database file


def test_load_known_formulas_set_reads_explicit_csv(tmp_path):
    from novelty_filter import load_known_formulas_set

    path = tmp_path / "known.csv"
    path.write_text("formula,databases\nNbAl2Cr,cod;mp\nTiO2,mp\n", encoding="utf-8")

    assert load_known_formulas_set(csv_path=path) == {"NbAl2Cr", "TiO2"}


def test_load_known_formulas_set_missing_explicit_csv_raises(tmp_path):
    from novelty_filter import load_known_formulas_set

    with pytest.raises(FileNotFoundError):
        load_known_formulas_set(csv_path=tmp_path / "does_not_exist.csv")


def test_load_known_index_from_merged_csv_keeps_provenance(tmp_path):
    from novelty_filter import load_known_index_from_csvs

    path = tmp_path / "known.csv"
    path.write_text(
        "formula,databases\nNbAl2Cr,cod;mp\nAc(AgAs)2,cod\n", encoding="utf-8")
    index = load_known_index_from_csvs(csv_path=path)

    assert index.formula_in("NbAl2Cr") is True
    assert index.formula_in("NbAl2Cr", database="cod") is True
    assert index.formula_in("NbAl2Cr", database="oqmd") is False
    assert index.formula_in("AcAg2As2", database="cod") is True   # flattened key
    assert set(index.databases_containing("NbAl2Cr")) == {"cod", "mp"}


def test_load_known_index_from_csvs_legacy_schema(tmp_path):
    from novelty_filter import load_known_index_from_csvs

    path = tmp_path / "legacy.csv"
    path.write_text(
        "formula,database,anonymous_formula,elements\n"
        "NbAl2Cr,COD,ABC2,Al Cr Nb\n"
        "TiO2,MP,AB2,O Ti\n",
        encoding="utf-8",
    )
    index = load_known_index_from_csvs(csv_path=path)

    assert index.formula_in("TiO2", database="mp") is True
    assert index.formula_in("TiO2", database="cod") is False
    assert "ABC2" in index.anonymous
    assert frozenset({"Al", "Cr", "Nb"}) in index.element_combos


def test_candidate_cli_filter_known_uses_merged_csv(tmp_path):
    from argparse import Namespace

    from candidate_cli import _load_known

    path = tmp_path / "known.csv"
    path.write_text("formula,databases\nNbAl2Cr,cod;mp\nTiO2,mp\n", encoding="utf-8")

    known = _load_known(Namespace(known_csv=str(path), db=None))
    assert known == {"NbAl2Cr", "TiO2"}
