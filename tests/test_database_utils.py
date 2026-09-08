"""Offline tests for database helpers (no network required)."""

from __future__ import annotations

import pandas as pd
import pytest

from crystalfno_discovery import database_utils as du
from crystalfno_discovery.databases import cod, record_to_metadata
from crystalfno_discovery.novelty import build_known_index


def test_structure_from_optimade():
    attrs = {
        "lattice_vectors": [[2.8, 0, 0], [0, 2.8, 0], [0, 0, 2.8]],
        "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.4, 1.4, 1.4]],
        "species_at_sites": ["Li", "Cl"],
        "elements": ["Li", "Cl"],
    }
    struct, meta, status = record_to_metadata(attrs, "OQMD", "oq-1")
    assert struct is not None
    assert meta["reduced_formula"] == "LiCl"
    assert status == du.STATUS_VALID


def test_structure_from_optimade_missing_coords():
    struct, meta, status = record_to_metadata({"elements": ["Na"]}, "COD", "c1")
    assert struct is None
    assert status == du.STATUS_MISSING_COORDINATES


def test_cod_metadata_from_attributes():
    attrs = {
        "chemical_formula_reduced": "Al2O3",
        "lattice_vectors": [[4.76, 0, 0], [0, 4.76, 0], [0, 0, 12.99]],
        "space_group_symbol_hermann_mauguin": "R-3c",
    }
    meta = cod.metadata_from_attributes(attrs, "9000000")
    assert meta["database"] == "COD"
    assert meta["reduced_formula"] == "Al2O3"
    assert meta["lattice_a"] > 0
    assert meta["space_group"] == "R-3c"


def test_cod_cif_url_template():
    url = cod.cif_url(9000001)
    assert url.endswith("9000001.cif") and "crystallography.net/cod/" in url


def test_build_known_index_from_known_df():
    frame = pd.DataFrame(
        [
            {"database": "COD", "formula": "Al2O3", "anonymous_formula": "A2B3", "elements": "Al O"},
            {"database": "OQMD", "formula": "Al2O3", "anonymous_formula": "A2B3", "elements": "Al O"},
            {"database": "MP", "formula": "Fe2O3", "anonymous_formula": "A2B3", "elements": "Fe O"},
        ]
    )
    index = build_known_index(frame)
    assert index.formula_in("Al2O3", database="COD") is True
    assert index.formula_in("Al2O3", database="OQMD") is True
    assert index.formula_in("Al2O3", database="MP") is False
    assert index.formula_in("Fe2O3") is True
