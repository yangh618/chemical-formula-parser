"""Robust formula parsing, canonicalization and composition-vector utilities.

All heavy lifting is delegated to :class:`pymatgen.core.Composition`. In recent
pymatgen versions (>= 2026) ``Composition.reduced_formula`` only reduces
*integer* amounts; fractional inputs such as ``Al0.6667Cr0.3333Nb0.3333`` are
not collapsed. We therefore canonicalize through
``Composition.get_integer_formula_and_factor``, which finds the smallest
integer ratio consistent (within tolerance) with the fractional input.
"""

from __future__ import annotations

import re

import numpy as np
from pymatgen.core import Composition
from pymatgen.core.periodic_table import Element


class FormulaError(ValueError):
    """Raised for malformed/unrepresentable formulas.

    The package never silently coerces an invalid formula into a valid one:
    callers handle :class:`FormulaError` or use the ``safe_*`` helpers that
    return ``None``.
    """


# Length of the elemental fraction vector (118 elements).
_VEC_LEN = 118


def _parse(formula: str | Composition) -> Composition:
    """Parse *formula* into a validated :class:`Composition`.

    Raises :class:`FormulaError` for empty, unparseable, charged, dummy or
    non-positive compositions.
    """
    if isinstance(formula, Composition):
        comp = formula
    else:
        if formula is None or not str(formula).strip():
            raise FormulaError("Empty formula")
        text = str(formula).strip()
        try:
            comp = Composition(text)
        except Exception as exc:  # noqa: BLE001 - pymatgen raises several types
            raise FormulaError(f"Could not parse formula: {text!r}") from exc

    if len(comp) == 0:
        raise FormulaError("Formula contains no elements")

    for el, amount in comp.items():
        charge = getattr(el, "oxi_state", 0.0) or 0.0
        if abs(float(charge)) > 1e-9:
            raise FormulaError(f"Formula contains charged species: {el!s}")
        symbol = getattr(getattr(el, "element", el), "symbol", None)
        if symbol is None:
            raise FormulaError(f"Formula contains a non-element species: {el!s}")
        try:
            Element(symbol)
        except Exception as exc:  # noqa: BLE001
            raise FormulaError(f"Formula contains an unknown element: {symbol!r}") from exc
        if float(amount) <= 0.0:
            raise FormulaError(f"Non-positive amount for element {symbol}: {amount}")
    return comp


def is_valid_formula(formula: str) -> bool:
    """Return ``True`` if *formula* is a neutral, real-element composition."""
    try:
        _parse(formula)
    except FormulaError:
        return False
    return True


def normalize_formula(formula: str, max_denominator: int = 10000) -> str:
    """Return a canonical reduced (integer-ratio) formula string.

    ``Al2CrNb``, ``CrAl2Nb`` and ``Al0.6667Cr0.3333Nb0.3333`` all normalise to
    ``NbAl2Cr``. Raises :class:`FormulaError` for unparseable input.
    """
    comp = _parse(formula)
    reduced, _factor = comp.get_integer_formula_and_factor(max_denominator=max_denominator)
    return reduced


def safe_normalize_formula(formula: str) -> str | None:
    """``normalize_formula`` that returns ``None`` instead of raising."""
    try:
        return normalize_formula(formula)
    except FormulaError:
        return None


def flatten_formula(formula: str, max_denominator: int = 10000) -> str:
    """Return the reduced formula as a plain, parenthesis-free string.

    :func:`normalize_formula` preserves polyatomic groupings, so it renders
    ``La4S3O3`` as ``La4(SO)3``. This expands such groupings into a flat
    element-by-count string, omitting counts of 1::

        flatten_formula("La4(SO)3")  -> "La4S3O3"
        flatten_formula("Ti(PS3)2")  -> "TiP2S6"
        flatten_formula("La2O3")     -> "La2O3"

    Element order follows pymatgen's conventional chemical ordering (the order
    used by :func:`normalize_formula`) rather than atomic number, so cations
    are written before anions.

    The result is canonical and notation-independent: any two formulas with the
    same reduced composition flatten to the same string, which makes it a safe
    comparison key as well as a display form.

    Raises :class:`FormulaError` for unparseable input.
    """
    canonical = normalize_formula(formula, max_denominator=max_denominator)
    amounts = _expand_groups(canonical)
    if not amounts:
        raise FormulaError(f"Formula contains no elements: {formula!r}")

    parts = []
    for symbol, amount in amounts.items():
        count = int(round(float(amount)))
        if count <= 0:
            raise FormulaError(f"Non-positive amount for element {symbol} in {formula!r}")
        parts.append(symbol if count == 1 else f"{symbol}{count}")
    return "".join(parts)


def safe_flatten_formula(formula: str) -> str | None:
    """``flatten_formula`` that returns ``None`` instead of raising."""
    try:
        return flatten_formula(formula)
    except FormulaError:
        return None


_ELEMENT_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def _expand_groups(text: str) -> dict[str, float]:
    """Expand parenthesised groups in *text*, preserving first-appearance order.

    Returns an ordered ``{element: amount}`` mapping. Handles nesting, e.g.
    ``Ag(AuF4)2`` -> ``{'Ag': 1, 'Au': 2, 'F': 8}``.
    """
    result: dict[str, float] = {}
    i = 0
    n = len(text)

    while i < n:
        if text[i] == "(":
            depth = 0
            j = i
            while j < n:
                if text[j] == "(":
                    depth += 1
                elif text[j] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            if depth != 0:
                raise FormulaError(f"Unbalanced parentheses in formula: {text!r}")

            inner = text[i + 1 : j]
            k = j + 1
            while k < n and (text[k].isdigit() or text[k] == "."):
                k += 1
            mult_str = text[j + 1 : k]
            multiplier = float(mult_str) if mult_str else 1.0

            for symbol, amount in _expand_groups(inner).items():
                result[symbol] = result.get(symbol, 0.0) + amount * multiplier
            i = k
        else:
            m = _ELEMENT_TOKEN.match(text, i)
            if not m:
                raise FormulaError(f"Could not parse formula fragment: {text[i:]!r}")
            symbol = m.group(1)
            count_str = m.group(2)
            count = float(count_str) if count_str else 1.0
            result[symbol] = result.get(symbol, 0.0) + count
            i = m.end()

    return result


def get_anonymous_formula(formula: str) -> str:
    """Anonymous stoichiometry (``A``, ``B``, ``C``, ...) of *formula*.

    Depends only on the multiset of coefficients, so e.g. both ``Al2CrNb`` and
    ``Ti2VTa`` map to ``ABC2``. Raises :class:`FormulaError` if invalid.
    """
    reduced = normalize_formula(formula)
    return Composition(reduced).anonymized_formula


def get_elements(formula: str) -> list[str]:
    """Element symbols of *formula*, sorted by atomic number."""
    comp = _parse(formula)
    symbols = {getattr(getattr(el, "element", el), "symbol", str(el)) for el in comp.elements}
    return sorted(symbols, key=lambda sym: Element(sym).Z)


def get_num_elements(formula: str) -> int:
    """Number of distinct chemical elements in *formula*."""
    return len(get_elements(formula))


def _reduced_amounts(formula: str, max_denominator: int = 10000) -> dict[str, int]:
    """Return a symbol -> integer-count map for the canonical reduced formula."""
    reduced = normalize_formula(formula, max_denominator=max_denominator)
    comp = Composition(reduced)
    return {
        getattr(getattr(el, "element", el), "symbol", str(el)): int(round(float(amount)))
        for el, amount in comp.items()
    }


def get_stoichiometric_vector(formula: str, max_denominator: int = 10000) -> np.ndarray:
    """Length-118 integer count vector of the reduced composition.

    Index ``i`` holds the (reduced) atom count of element with Z = ``i + 1``.
    """
    vec = np.zeros(_VEC_LEN, dtype=int)
    for symbol, count in _reduced_amounts(formula, max_denominator=max_denominator).items():
        vec[Element(symbol).Z - 1] = count
    return vec


def get_elemental_fraction_vector(formula: str) -> np.ndarray:
    """Length-118 elemental *fraction* vector (sums to 1)."""
    comp = _parse(formula)
    amounts = {
        getattr(getattr(el, "element", el), "symbol", str(el)): float(amount)
        for el, amount in comp.items()
    }
    total = sum(amounts.values())
    if total <= 0.0:
        raise FormulaError("Formula has non-positive total atom count")
    vec = np.zeros(_VEC_LEN, dtype=float)
    for symbol, amount in amounts.items():
        vec[Element(symbol).Z - 1] = amount / total
    return vec


def element_set(formula: str) -> frozenset[str]:
    """Set of element symbols present in *formula*."""
    return frozenset(get_elements(formula))


# ---------------------------------------------------------------------------
# Composition-space distances (section 18).
# ---------------------------------------------------------------------------

def composition_distance(formula_a: str, formula_b: str) -> float:
    """Euclidean distance between 118-D elemental fraction vectors.

    Deliberately a *simple baseline* metric, not a learned descriptor.
    """
    va = get_elemental_fraction_vector(formula_a)
    vb = get_elemental_fraction_vector(formula_b)
    return float(np.linalg.norm(va - vb))


def element_set_distance(formula_a: str, formula_b: str) -> float:
    """``1 - Jaccard`` index of the element sets (0 identical, 1 disjoint)."""
    sa = element_set(formula_a)
    sb = element_set(formula_b)
    if not sa and not sb:
        return 0.0
    return 1.0 - len(sa & sb) / len(sa | sb)

