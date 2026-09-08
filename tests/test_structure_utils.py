"""Tests for structure utilities (metadata, validation, comparison)."""

from __future__ import annotations

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from crystalfno_discovery.structure_utils import (
    compare_structures,
    lattice_errors,
    radial_distribution,
    rdf_similarity,
    save_cif,
    site_issues,
    structure_is_ordered,
    structure_to_metadata,
    structures_match,
    validate_structure,
)


def _nacl():
    na = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
    cl = [[0.5, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]
    return Structure(Lattice.from_parameters(5.64, 5.64, 5.64, 90, 90, 90), ["Na"] * 4 + ["Cl"] * 4, na + cl)


def test_metadata_columns():
    st = _nacl()
    meta = structure_to_metadata(st, database="COD", database_id="1")
    assert meta["reduced_formula"] == "NaCl"
    assert meta["space_group"] == "Fm-3m"
    assert meta["space_group_number"] == 225
    assert meta["nsites"] == 8
    assert meta["nelements"] == 2
    assert meta["volume"] == pytest.approx(5.64 ** 3)
    assert structure_is_ordered(st)


def test_validation_valid():
    st = _nacl()
    out = validate_structure(st, expected_formula="NaCl")
    assert out["valid"] is True
    assert out["reason"] is None
    assert out["min_distance"] > 2.0


def test_validation_composition_mismatch():
    st = _nacl()
    out = validate_structure(st, expected_formula="Fe2O3")
    assert out["valid"] is False
    assert "composition mismatch" in out["reason"]


def test_validation_overlapping_atoms():
    lat = Lattice.from_parameters(4, 4, 4, 90, 90, 90)
    st = Structure(lat, ["Fe", "Fe"], [[0, 0, 0], [0.001, 0, 0]])
    out = validate_structure(st)
    assert out["valid"] is False
    assert out["overlapping_atoms"] is True


def test_validation_invalid_lattice():
    # Zero volume lattice is rejected.
    bad = Lattice([[1, 0, 0], [0, 0, 0], [0, 0, 1]])
    st = Structure(bad, ["Fe"], [[0, 0, 0]])
    out = validate_structure(st)
    assert out["valid"] is False


def test_disorder_and_partial_occupancy():
    lat = Lattice.from_parameters(3, 3, 3, 90, 90, 90)
    dis = Structure(lat, [{"Fe": 0.5, "Co": 0.5}], [[0, 0, 0]])
    partial = Structure(lat, [{"Fe": 0.75}], [[0, 0, 0]])
    assert site_issues(dis) == (True, True)
    assert site_issues(partial) == (True, False)
    assert structure_is_ordered(dis) is False


def test_cif_roundtrip_and_parse(tmp_path):
    st = _nacl()
    path = tmp_path / "nacl.cif"
    save_cif(st, str(path))
    from crystalfno_discovery.database_utils import parse_structure_file

    _s, meta, status, _r = parse_structure_file(path, database="COD", database_id="t")
    assert status == "valid"
    assert meta["reduced_formula"] == "NaCl"


def test_rdf_and_comparison():
    st = _nacl()
    r, g = radial_distribution(st, r_max=5.0, n_bins=30)
    assert r.shape == g.shape == (30,)
    assert np.isfinite(g).all()
    assert g.max() > 1.0
    assert rdf_similarity(st, st) == pytest.approx(1.0)

    scaled = st.copy()
    scaled.scale_lattice(st.volume * 1.05)
    rep = compare_structures(st, scaled)
    assert rep["composition_match"] is True
    assert rep["lattice_errors"]["volume_error"] == pytest.approx(0.05, abs=0.01)
    assert "rdf_similarity" in rep


def test_structures_match_duplicate():
    st1 = _nacl()
    st2 = _nacl()
    st2.translate_sites(0, [0.0, 0.0, 0.0], frac_coords=True)
    assert structures_match(st1, st2) is True
    # An unrelated structure should not match.
    other = Structure(Lattice.from_parameters(4, 4, 4, 90, 90, 90), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    assert structures_match(st1, other) is False or lattice_errors(st1, other)["volume_error"] > 0.1
