"""Structure inspection, validation, metadata and comparison utilities.

Depends on :mod:`pymatgen.core`. Two pymatgen-2026 notes that matter here:

* space-group symmetry is obtained with ``pymatgen.symmetry.analyzer
  .SpacegroupAnalyzer`` (``SymmetryFinder`` was removed);
* the old ``pymatgen.analysis.radial_distribution_function`` module is gone,
  so a small self-contained periodic RDF is implemented here.

All structure inputs are :class:`pymatgen.core.Structure` objects.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Composition, Structure
from pymatgen.core.periodic_table import Element

from .formula_utils import FormulaError, normalize_formula

# Tolerance (Angstrom) below which two atoms are considered overlapping.
_OVERLAP_TOL = 0.5
# Smallest physically meaningful interatomic distance (Angstrom).
_DEFAULT_MIN_DIST = 0.5

# Metadata columns written for every structure (section 3).
METADATA_COLUMNS = [
    "database",
    "database_id",
    "formula",
    "reduced_formula",
    "anonymous_formula",
    "elements",
    "nelements",
    "nsites",
    "lattice_a",
    "lattice_b",
    "lattice_c",
    "alpha",
    "beta",
    "gamma",
    "volume",
    "space_group",
    "space_group_number",
    "structure_file",
]


def _composition_species_map(site: Any) -> dict[str, float]:
    """Return ``{symbol: fractional_occupancy}`` for one pymatgen site."""
    species = site.species
    if isinstance(species, Element) or not hasattr(species, "items"):
        symbol = getattr(getattr(species, "element", species), "symbol", str(species))
        return {symbol: 1.0}
    out: dict[str, float] = {}
    for sp, amt in species.items():
        symbol = getattr(getattr(sp, "element", sp), "symbol", str(sp))
        out[symbol] = out.get(symbol, 0.0) + float(amt)
    return out


def site_issues(structure: Structure) -> tuple[bool, bool]:
    """Return ``(has_partial_occupancy, has_disorder)`` for *structure*.

    A site is considered *disordered* when more than one chemical element
    shares it and *partially occupied* when any species has occupancy < 1.
    """
    has_partial = False
    has_disorder = False
    for site in structure:
        mapping = _composition_species_map(site)
        if len(mapping) > 1:
            has_disorder = True
        if any(abs(occ - 1.0) > 1e-6 for occ in mapping.values()):
            has_partial = True
    return has_partial, has_disorder


def structure_is_ordered(structure: Structure) -> bool:
    """Return ``True`` if every site is fully occupied by a single element."""
    partial, disordered = site_issues(structure)
    return not partial and not disordered


def get_space_group(structure: Structure) -> tuple[str, int]:
    """Return ``(symbol, international_number)`` for *structure*."""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    finder = SpacegroupAnalyzer(structure, symprec=1e-3)
    try:
        symbol = finder.get_space_group_symbol()
        number = finder.get_space_group_number()
    except Exception as exc:  # noqa: BLE001 - symmetry may be undetermined
        warnings.warn(f"Space-group detection failed: {exc}")
        symbol, number = "P1", 1
    return symbol, int(number)


def structure_to_metadata(
    structure: Structure,
    database: str = "",
    database_id: str = "",
    structure_file: str = "",
    **extra: Any,
) -> dict[str, Any]:
    """Convert a :class:`Structure` to the unified metadata row (section 3).

    ``extra`` key/value pairs are preserved so database-specific information
    is not discarded.
    """
    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    reduced = structure.composition.reduced_formula
    sym, num = get_space_group(structure)
    elements = sorted(
        {getattr(getattr(el, "element", el), "symbol", str(el)) for el in structure.composition.elements},
        key=lambda s: Element(s).Z,
    )
    row: dict[str, Any] = {
        "database": database,
        "database_id": database_id,
        "formula": structure.composition.formula,
        "reduced_formula": reduced,
        "anonymous_formula": structure.composition.anonymized_formula,
        "elements": " ".join(elements),
        "nelements": len(elements),
        "nsites": structure.num_sites,
        "lattice_a": float(a),
        "lattice_b": float(b),
        "lattice_c": float(c),
        "alpha": float(alpha),
        "beta": float(beta),
        "gamma": float(gamma),
        "volume": float(structure.volume),
        "space_group": sym,
        "space_group_number": num,
        "structure_file": structure_file,
    }
    row.update(extra)
    return row


def save_cif(structure: Structure, path: str) -> str:
    """Write *structure* as a CIF to *path* (creating parent dirs)."""
    from pathlib import Path

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    structure.to(filename=str(path))
    return path


def _validate_lattice(structure: Structure) -> tuple[bool, str | None]:
    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    lengths = np.array([a, b, c], dtype=float)
    if not np.all(np.isfinite(lengths)) or np.any(lengths <= 0.0):
        return False, "Non-positive or non-finite lattice lengths"
    if not np.all(np.isfinite([alpha, beta, gamma])):
        return False, "Non-finite lattice angles"
    if not np.all((np.array([alpha, beta, gamma]) > 0.0) & (np.array([alpha, beta, gamma]) < 180.0)):
        return False, "Lattice angles outside (0, 180) degrees"
    if not structure.volume > 0.0:
        return False, "Non-positive cell volume"
    return True, None


def validate_structure(
    structure: Structure,
    expected_formula: str | None = None,
    min_distance: float = _DEFAULT_MIN_DIST,
) -> dict[str, Any]:
    """Run structural sanity checks (section 20).

    Checks: valid lattice, positive volume, finite in-cell fractional
    coordinates, composition consistency, atomic-species/formula match,
    a minimum interatomic distance, atom overlap and periodic consistency.

    Returns::

        {"valid": bool, "reason": str|None, "min_distance": float,
         "volume": float, "overlapping_atoms": bool, "partial_occupancy": bool,
         "disordered": bool}
    """
    result: dict[str, Any] = {
        "valid": False,
        "reason": None,
        "min_distance": float("inf"),
        "volume": float(structure.volume),
        "nsites": structure.num_sites,
        "overlapping_atoms": False,
        "partial_occupancy": False,
        "disordered": False,
    }

    ok, why = _validate_lattice(structure)
    if not ok:
        result["reason"] = f"lattice: {why}"
        return result
    if structure.num_sites < 1:
        result["reason"] = "no atoms"
        return result

    frac = structure.frac_coords
    if not np.all(np.isfinite(frac)):
        result["reason"] = "non-finite fractional coordinates"
        return result
    wrapped = frac % 1.0
    if not np.allclose(frac, wrapped, atol=1e-6):
        structure = Structure(structure.lattice, structure.species, wrapped)

    if expected_formula is not None:
        try:
            expected = normalize_formula(expected_formula)
        except FormulaError:
            expected = None
        if expected is not None:
            actual = structure.composition.reduced_formula
            if actual != expected:
                result["reason"] = f"composition mismatch: {actual} != {expected}"
                return result

    try:
        matrix = structure.distance_matrix
    except Exception as exc:  # noqa: BLE001
        result["reason"] = f"could not compute distances: {exc}"
        return result
    if matrix.size == 0:
        result["reason"] = "empty distance matrix"
        return result

    n = structure.num_sites
    off_diag = matrix[~np.eye(n, dtype=bool)]
    min_dist = float(off_diag.min()) if off_diag.size else float("inf")
    result["min_distance"] = min_dist
    result["overlapping_atoms"] = bool(min_dist < _OVERLAP_TOL)
    if min_dist < min_distance:
        result["reason"] = f"minimum interatomic distance {min_dist:.3f} < {min_distance}"
        return result

    result["partial_occupancy"], result["disordered"] = site_issues(structure)
    result["valid"] = True
    return result


def radial_distribution(structure: Structure, r_max: float = 6.0, n_bins: int = 60) -> tuple[np.ndarray, np.ndarray]:
    """Periodic radial distribution function, returning ``(r, g(r))``.

    A simple baseline used for structure comparison (the pymatgen RDF module
    was removed in recent versions). Uses a translated-supercell neighbour
    list and normalises by the mean number density of the cell.
    """
    from itertools import product

    matrix = np.asarray(structure.lattice.matrix, dtype=float)  # rows = lattice vectors
    volume = float(structure.volume)
    if not np.all(np.isfinite(matrix)) or volume <= 0.0:
        raise ValueError("Cannot compute RDF for an invalid lattice")

    n = structure.num_sites
    if n < 2:
        raise ValueError("RDF needs at least two atoms")

    r_min = 1e-3
    dr = (r_max - r_min) / n_bins
    radii = np.linspace(r_min + dr / 2, r_max - dr / 2, n_bins)
    reps = [int(np.ceil(r_max / max(np.linalg.norm(matrix[d]), 1e-12))) + 1 for d in range(3)]

    base = np.asarray(structure.cart_coords, dtype=float)
    translations = list(product(*[range(-reps[d], reps[d] + 1) for d in range(3)]))
    shifts = np.array(translations, dtype=float)
    shift_cart = shifts @ matrix  # (ntrans, 3)

    ntrans = len(translations)
    super_coords = (base[None, :, :] + shift_cart[:, None, :]).reshape(-1, 3)  # (ntrans*n, 3)

    counts = np.zeros(n_bins, dtype=float)
    zero_t_indices = [ti for ti, t in enumerate(translations) if all(x == 0 for x in t)]

    for i in range(n):
        dists = np.linalg.norm(super_coords - base[i], axis=1)
        drop_idx = [ti * n + i for ti in zero_t_indices]
        dists = np.delete(dists, drop_idx)
        dists = dists[(dists >= r_min) & (dists < r_max)]
        if dists.size:
            hist, _ = np.histogram(dists, bins=n_bins, range=(r_min, r_max))
            counts += hist

    avg = counts / n
    rho = n / volume
    shell_vol = 4.0 * np.pi / 3.0 * ((radii + dr / 2) ** 3 - (radii - dr / 2) ** 3)
    shell_vol = np.maximum(shell_vol, 1e-30)
    g = avg / (rho * shell_vol)
    return radii, g


def rdf_similarity(structure_a: Structure, structure_b: Structure, r_max: float = 6.0) -> float:
    """Pearson correlation of two structures' RDFs (higher = more alike)."""
    _, ga = radial_distribution(structure_a, r_max=r_max)
    _, gb = radial_distribution(structure_b, r_max=r_max)
    if ga.std() == 0.0 or gb.std() == 0.0:
        return float("nan")
    return float(np.corrcoef(ga, gb)[0, 1])


def structures_match(structure_a: Structure, structure_b: Structure) -> bool:
    """StructureMatcher duplicate test (composition + lattice + atomic order)."""
    matcher = StructureMatcher(primitive_cell=False, attempt_supercell=True)
    try:
        return bool(matcher.fit(structure_a, structure_b))
    except Exception:  # noqa: BLE001
        return False


def lattice_errors(structure_ref: Structure, structure_pred: Structure) -> dict[str, float]:
    """Relative (absolute) lattice-parameter, volume and density errors."""
    rp = structure_ref.lattice.parameters
    pp = structure_pred.lattice.parameters
    out = {
        "a_error": abs(pp[0] - rp[0]) / rp[0],
        "b_error": abs(pp[1] - rp[1]) / rp[1],
        "c_error": abs(pp[2] - rp[2]) / rp[2],
        "volume_error": abs(structure_pred.volume - structure_ref.volume) / structure_ref.volume,
    }
    if structure_ref.density > 0:
        out["density_error"] = abs(structure_pred.density - structure_ref.density) / structure_ref.density
    else:
        out["density_error"] = float("nan")
    return out


def compare_structures(structure_ref: Structure, structure_pred: Structure) -> dict[str, Any]:
    """Full structural comparison report (section 21)."""
    ref_meta = structure_to_metadata(structure_ref)
    pred_meta = structure_to_metadata(structure_pred)
    return {
        "composition_match": ref_meta["reduced_formula"] == pred_meta["reduced_formula"],
        "space_group_match": ref_meta["space_group"] == pred_meta["space_group"],
        "ref_space_group": ref_meta["space_group"],
        "pred_space_group": pred_meta["space_group"],
        "lattice_errors": lattice_errors(structure_ref, structure_pred),
        "rdf_similarity": rdf_similarity(structure_ref, structure_pred),
        "structure_matcher_match": structures_match(structure_ref, structure_pred),
    }

