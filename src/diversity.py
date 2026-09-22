"""Diversity-aware selection of composition candidates.

Naive enumeration of a composition space returns candidates in a deterministic
but *clustered* order: every stoichiometry of the first element combination is
emitted before the next combination is reached. The head of such a list is
dominated by a handful of species — hundreds of ``La-S-O`` variants before any
``Ce``-containing formula appears.

This module selects a subset that is spread out over composition space using
greedy **max-min diversity** (farthest-point sampling): repeatedly take the
candidate that is farthest from everything already chosen.

The distance combines two complementary terms:

``composition``
    Euclidean distance between 118-D elemental *fraction* vectors. Separates
    different stoichiometries of the same elements.
``element set``
    ``1 - Jaccard`` similarity of the element sets. Separates different
    chemical species, which is exactly the axis that enumeration order
    collapses.

Both terms are simple, deterministic geometric baselines — no learned model.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np
from pymatgen.core.periodic_table import Element

from formula_utils import (
    FormulaError,
    get_elemental_fraction_vector,
    get_elements,
)

logger = logging.getLogger(__name__)

_VEC_LEN = 118


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

@dataclass
class DiversityReport:
    """Summary statistics describing how spread out a set of formulas is."""

    n_formulas: int
    mean_pairwise_distance: float
    min_pairwise_distance: float
    distinct_element_sets: int
    elements_covered: int
    most_common_element_set: tuple[str, ...]
    most_common_element_set_count: int
    sampled: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_formulas": self.n_formulas,
            "mean_pairwise_distance": round(self.mean_pairwise_distance, 4),
            "min_pairwise_distance": round(self.min_pairwise_distance, 4),
            "distinct_element_sets": self.distinct_element_sets,
            "elements_covered": self.elements_covered,
            "most_common_element_set": "".join(self.most_common_element_set),
            "most_common_element_set_count": self.most_common_element_set_count,
            "sampled": self.sampled,
        }


# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------

def build_features(
    formulas: Sequence[str],
) -> tuple[list[int], np.ndarray, np.ndarray]:
    """Build vectorised features for *formulas*.

    Returns ``(valid_indices, fraction_matrix, presence_matrix)`` where
    ``valid_indices[i]`` is the position of row ``i`` in the input list.
    Unparseable formulas are dropped rather than raising.
    """
    frac = np.zeros((len(formulas), _VEC_LEN), dtype=float)
    pres = np.zeros((len(formulas), _VEC_LEN), dtype=float)
    valid: list[int] = []

    row = 0
    for i, formula in enumerate(formulas):
        try:
            frac[row] = get_elemental_fraction_vector(formula)
            elements = get_elements(formula)
        except (FormulaError, Exception):  # noqa: BLE001
            continue
        for el in elements:
            pres[row, Element(el).Z - 1] = 1.0
        valid.append(i)
        row += 1

    return valid, frac[:row], pres[:row]


def _distance_to(
    frac: np.ndarray,
    pres: np.ndarray,
    set_sizes: np.ndarray,
    j: int,
    w_composition: float,
    w_element_set: float,
) -> np.ndarray:
    """Distance from row *j* to every row (vectorised)."""
    d_comp = np.linalg.norm(frac - frac[j], axis=1)

    intersection = pres @ pres[j]
    union = set_sizes + set_sizes[j] - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard = np.where(union > 0, intersection / np.where(union > 0, union, 1.0), 0.0)
    d_el = 1.0 - jaccard

    return w_composition * d_comp + w_element_set * d_el


def _pairwise_matrix(
    frac: np.ndarray,
    pres: np.ndarray,
    w_composition: float,
    w_element_set: float,
) -> np.ndarray:
    """Full pairwise distance matrix (O(N^2) memory — use on modest N)."""
    sq = np.einsum("ij,ij->i", frac, frac)
    d_comp2 = np.maximum(sq[:, None] + sq[None, :] - 2.0 * (frac @ frac.T), 0.0)
    d_comp = np.sqrt(d_comp2)

    intersection = pres @ pres.T
    sizes = pres.sum(axis=1)
    union = sizes[:, None] + sizes[None, :] - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard = np.where(union > 0, intersection / np.where(union > 0, union, 1.0), 0.0)
    d_el = 1.0 - jaccard

    return w_composition * d_comp + w_element_set * d_el


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def select_diverse(
    formulas: Iterable[str],
    n_select: int,
    w_composition: float = 1.0,
    w_element_set: float = 1.0,
) -> list[str]:
    """Select up to *n_select* formulas that are maximally spread out.

    Greedy farthest-point sampling: the first pick is the candidate farthest
    from the pool centroid, and each subsequent pick is the candidate whose
    distance to its *nearest* already-selected neighbour is largest. This is a
    deterministic 1/2-approximation to the max-min dispersion problem.

    Parameters
    ----------
    formulas
        Candidate pool. Duplicates are collapsed, order is otherwise irrelevant.
    n_select
        How many formulas to return. Clamped to the pool size.
    w_composition, w_element_set
        Weights on the stoichiometry and species-identity distance terms.

    Returns
    -------
    list[str]
        Selected formulas, in the order they were chosen.
    """
    unique: list[str] = list(dict.fromkeys(formulas))
    if not unique:
        return []

    valid_idx, frac, pres = build_features(unique)
    n_pool = frac.shape[0]
    if n_pool == 0:
        return []

    k = min(int(n_select), n_pool)
    if k <= 0:
        return []
    if k == n_pool:
        return [unique[i] for i in valid_idx]

    set_sizes = pres.sum(axis=1)

    # Deterministic seed: the point farthest from the pool centroid.
    centroid = frac.mean(axis=0)
    seed = int(np.argmax(np.linalg.norm(frac - centroid, axis=1)))

    selected: list[int] = [seed]
    min_dist = _distance_to(frac, pres, set_sizes, seed, w_composition, w_element_set)
    min_dist[seed] = -np.inf

    for _ in range(k - 1):
        nxt = int(np.argmax(min_dist))
        selected.append(nxt)
        d = _distance_to(frac, pres, set_sizes, nxt, w_composition, w_element_set)
        min_dist = np.minimum(min_dist, d)
        min_dist[selected] = -np.inf

    return [unique[valid_idx[i]] for i in selected]


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def diversity_report(
    formulas: Sequence[str],
    w_composition: float = 1.0,
    w_element_set: float = 1.0,
    max_pairs: int = 3000,
) -> DiversityReport:
    """Quantify how spread out *formulas* are.

    Pairwise statistics are computed on a deterministic subsample when the set
    is larger than *max_pairs*, to keep memory bounded.
    """
    unique: list[str] = list(dict.fromkeys(formulas))
    if not unique:
        return DiversityReport(
            n_formulas=0,
            mean_pairwise_distance=0.0,
            min_pairwise_distance=0.0,
            distinct_element_sets=0,
            elements_covered=0,
            most_common_element_set=(),
            most_common_element_set_count=0,
        )

    valid_idx, frac, pres = build_features(unique)
    n = frac.shape[0]

    sampled = False
    if n > max_pairs:
        rng = np.random.default_rng(0)
        take = np.sort(rng.choice(n, size=max_pairs, replace=False))
        frac_s, pres_s = frac[take], pres[take]
        sampled = True
    else:
        frac_s, pres_s = frac, pres

    if frac_s.shape[0] > 1:
        dist = _pairwise_matrix(frac_s, pres_s, w_composition, w_element_set)
        iu = np.triu_indices(dist.shape[0], k=1)
        pairs = dist[iu]
        mean_d = float(pairs.mean())
        min_d = float(pairs.min())
    else:
        mean_d = min_d = 0.0

    element_sets = [frozenset(get_elements(unique[i])) for i in valid_idx]
    counter = Counter(element_sets)
    top_set, top_count = counter.most_common(1)[0]
    all_elements: set[str] = set()
    for es in element_sets:
        all_elements |= es

    return DiversityReport(
        n_formulas=len(unique),
        mean_pairwise_distance=mean_d,
        min_pairwise_distance=min_d,
        distinct_element_sets=len(counter),
        elements_covered=len(all_elements),
        most_common_element_set=tuple(sorted(top_set, key=lambda s: Element(s).Z)),
        most_common_element_set_count=top_count,
        sampled=sampled,
    )
