"""Rules of an ``energy_decomposition`` payload that JSON Schema cannot express (spec § 5.16).

Shared by the validator (:mod:`qvf.validate`) and the writer
(:meth:`qvf.QvfWriter.add_energy_decomposition`), so that a payload the writer accepts is one
the validator accepts. Standard library only.

License: Apache-2.0.
"""

from __future__ import annotations

import math
from typing import Any, Optional

#: Hartree per unit of every energy unit a decomposition may use (§ 4.3,
#: CODATA 2018 as in § 3.6), legacy spellings of hartree and ev included.
ENERGY_UNITS: dict[str, float] = {
    "hartree": 1.0,
    "Eh": 1.0,
    "Ha": 1.0,
    "Hartree": 1.0,
    "ev": 1.0 / 27.211386245988,
    "eV": 1.0 / 27.211386245988,
    "kcal_per_mol": 1.0 / 627.5094740631,
    "kj_per_mol": 1.0 / 2625.4996394799,
}

#: Default sum-rule tolerance, hartree (§ 5.16).
DEFAULT_SUM_TOLERANCE_EH = 1e-6

SECTION = "spec § 5.16"


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _is_num(value: Any) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value))


def energy_eh(energy: Any) -> Optional[float]:
    """An ``{value, units}`` energy in hartree, or None when it is not one
    (a wrong unit is reported separately)."""
    if not isinstance(energy, dict) or not _is_num(energy.get("value")):
        return None
    factor = ENERGY_UNITS.get(energy.get("units"))
    return None if factor is None else energy["value"] * factor


def _energies(payload: dict) -> list[tuple[str, Any]]:
    """Every ``{value, units}`` energy of a payload with its location."""
    found: list[tuple[str, Any]] = []
    for k, frag in enumerate(_list(payload.get("fragments"))):
        if isinstance(frag, dict) and "energy" in frag:
            found.append((f"fragments[{k}].energy", frag["energy"]))
    energies = payload.get("energies")
    if isinstance(energies, dict) and "system" in energies:
        found.append(("energies.system", energies["system"]))
    for k, term in enumerate(_list(payload.get("terms"))):
        if isinstance(term, dict) and "energy" in term:
            found.append((f"terms[{k}].energy", term["energy"]))
    for k, rule in enumerate(_list(payload.get("sum_rules"))):
        if isinstance(rule, dict) and "tolerance" in rule:
            found.append((f"sum_rules[{k}].tolerance", rule["tolerance"]))
    pairs = payload.get("orbital_pairs")
    if isinstance(pairs, dict):
        if "total" in pairs:
            found.append(("orbital_pairs.total", pairs["total"]))
        for k, pair in enumerate(_list(pairs.get("pairs"))):
            if isinstance(pair, dict) and "energy" in pair:
                found.append((f"orbital_pairs.pairs[{k}].energy", pair["energy"]))
    return found


def _list(value: Any) -> list:
    return value if isinstance(value, list) else []


def pair_keys(payload: Any) -> Optional[set[tuple[int, str]]]:
    """The ``(index, spin)`` of every orbital pair, or None when the payload
    has no pair table (for resolving ``volume.difference`` references)."""
    if not isinstance(payload, dict) or not isinstance(payload.get("orbital_pairs"), dict):
        return None
    keys = set()
    for pair in _list(payload["orbital_pairs"].get("pairs")):
        if isinstance(pair, dict) and _is_int(pair.get("index")):
            keys.add((pair["index"], pair.get("spin", "restricted")))
    return keys


def payload_errors(payload: Any, n_atoms: Optional[int] = None) -> list[tuple[str, str]]:
    """``(code, message)`` for every rule of § 5.16 an ``energy_decomposition`` payload breaks.

    ``n_atoms`` is the atom count of the section's structure when it
    resolves; without it the atom indices are not range-checked. Codes:
    ``E_MEMBER_MALFORMED`` (fragments, terms, indices, pairs, per-k arrays),
    ``E_UNIT_CONTRADICTED`` (an energy not in an energy unit) and
    ``E_SUM_RULE`` (a declared sum rule that does not hold). A payload that
    is not an object, or whose arrays have the wrong type, is left to the
    schema check.
    """
    if not isinstance(payload, dict):
        return []
    errors: list[tuple[str, str]] = []

    def bad(msg: str) -> None:
        errors.append(("E_MEMBER_MALFORMED", f"{msg} ({SECTION})"))

    # --- fragments ---------------------------------------------------------
    fragments = _list(payload.get("fragments"))
    n_frag = len(fragments)
    atoms_of: list[Optional[set[int]]] = []
    for frag in fragments:
        atoms = frag.get("atoms") if isinstance(frag, dict) else None
        atoms_of.append({a for a in atoms if _is_int(a)} if isinstance(atoms, list) else None)
    top = [k for k, f in enumerate(fragments) if isinstance(f, dict) and "parent" not in f]
    if isinstance(payload.get("fragments"), list) and len(top) < 2:
        bad(f"fragments has {len(top)} interacting fragment(s) (no parent); "
            f"a decomposition needs at least two")
    owner: dict[int, int] = {}
    for k in top:
        for a in sorted(atoms_of[k] or ()):
            if a in owner:
                bad(f"fragments[{k}] and fragments[{owner[a]}] both contain atom {a}; "
                    f"interacting fragments are disjoint")
            else:
                owner[a] = k
    for k, frag in enumerate(fragments):
        if not isinstance(frag, dict):
            continue
        if n_atoms is not None:
            for a in sorted(atoms_of[k] or ()):
                if not 0 <= a < n_atoms:
                    bad(f"fragments[{k}].atoms names atom {a}, but the structure has "
                        f"{n_atoms} atoms")
        if "parent" in frag:
            p = frag["parent"]
            if not _is_int(p) or not 0 <= p < n_frag or p == k:
                bad(f"fragments[{k}].parent {p!r} is not the index of another fragment")
            elif p not in top:
                bad(f"fragments[{k}].parent {p} is itself a sub-fragment; a sub-fragment "
                    f"subdivides an interacting fragment")
            elif atoms_of[k] is not None and atoms_of[p] is not None \
                    and not atoms_of[k] <= atoms_of[p]:
                extra = sorted(atoms_of[k] - atoms_of[p])
                bad(f"fragments[{k}] has atoms {extra} that its parent fragments[{p}] "
                    f"does not contain")

    # --- terms -------------------------------------------------------------
    terms = _list(payload.get("terms"))
    n_terms = len(terms)
    keys: dict[tuple, int] = {}
    interactions = []
    for k, term in enumerate(terms):
        if not isinstance(term, dict):
            continue
        frs = term.get("fragments")
        if isinstance(frs, list):
            for f in frs:
                if _is_int(f) and not 0 <= f < n_frag:
                    bad(f"terms[{k}].fragments names fragment {f}, but there are "
                        f"{n_frag} fragments")
            if term.get("directed") and len(frs) != 2:
                bad(f"terms[{k}] is directed but names {len(frs)} fragments; a directed "
                    f"term names [donor, acceptor]")
        if "parent" in term:
            p = term["parent"]
            if not _is_int(p) or not 0 <= p < n_terms or p == k:
                bad(f"terms[{k}].parent {p!r} is not the index of another term")
        scope = None
        if isinstance(frs, list):
            scope = tuple(frs) if term.get("directed") else tuple(sorted(
                f for f in frs if _is_int(f)))
        key = (term.get("name"), scope, term.get("spin"), term.get("parent"))
        if key in keys:
            bad(f"terms[{k}] repeats terms[{keys[key]}] (same name, fragments, spin "
                f"and parent)")
        else:
            keys[key] = k
        if (term.get("name") == "interaction" and "fragments" not in term
                and "spin" not in term and "parent" not in term):
            interactions.append(k)
    # cycles in the parent chain
    for k, term in enumerate(terms):
        seen = {k}
        node = term.get("parent") if isinstance(term, dict) else None
        while _is_int(node) and 0 <= node < n_terms:
            if node in seen:
                bad(f"terms[{k}] is its own ancestor through parent")
                break
            seen.add(node)
            parent_term = terms[node]
            node = parent_term.get("parent") if isinstance(parent_term, dict) else None
    if isinstance(payload.get("terms"), list) and len(interactions) != 1:
        bad(f"terms has {len(interactions)} interaction terms (name 'interaction' "
            f"without fragments, spin or parent); exactly one is required")

    # --- units -------------------------------------------------------------
    for where, energy in _energies(payload):
        if isinstance(energy, dict) and isinstance(energy.get("units"), str) \
                and energy["units"] not in ENERGY_UNITS:
            errors.append(("E_UNIT_CONTRADICTED",
                           f"{where} is in {energy['units']!r}, which is not an energy "
                           f"unit (hartree, ev, kcal_per_mol, kj_per_mol) ({SECTION})"))

    # --- sum rules ---------------------------------------------------------
    for k, rule in enumerate(_list(payload.get("sum_rules"))):
        if not isinstance(rule, dict):
            continue
        total, parts = rule.get("total"), rule.get("parts")
        idx = ([total] if _is_int(total) else []) + \
              [p for p in _list(parts) if _is_int(p)]
        out_of_range = [i for i in idx if not 0 <= i < n_terms]
        if out_of_range:
            bad(f"sum_rules[{k}] names term(s) {out_of_range}, but there are "
                f"{n_terms} terms")
            continue
        if not _is_int(total) or not isinstance(parts, list):
            continue
        if total in parts:
            bad(f"sum_rules[{k}] lists its total terms[{total}] among its parts")
            continue
        values = [energy_eh(terms[i].get("energy")) if isinstance(terms[i], dict) else None
                  for i in [total] + parts]
        tol = DEFAULT_SUM_TOLERANCE_EH
        if "tolerance" in rule:
            tol = energy_eh(rule["tolerance"])
        if tol is None or any(v is None for v in values):
            continue  # an energy that is not one is reported above or by the schema
        residual = values[0] - sum(values[1:])
        if abs(residual) > tol:
            errors.append(("E_SUM_RULE",
                           f"sum_rules[{k}]: terms[{total}] is {values[0]:.10g} hartree "
                           f"but its parts {parts} add up to {sum(values[1:]):.10g} "
                           f"hartree (residual {residual:.3g}, tolerance {tol:.3g}) "
                           f"({SECTION})"))

    # --- orbital pairs -----------------------------------------------------
    pairs = payload.get("orbital_pairs")
    n_k = None
    periodic = payload.get("periodic")
    if isinstance(periodic, dict) and isinstance(periodic.get("kpoints"), dict):
        kps = periodic["kpoints"].get("kpoints")
        if isinstance(kps, list):
            n_k = len(kps)
            weights = periodic["kpoints"].get("weights")
            if isinstance(weights, list) and len(weights) != n_k:
                bad(f"periodic.kpoints has {n_k} k-points but {len(weights)} weights")
    if isinstance(pairs, dict):
        if pairs.get("type") == "nocv" and "eigenvalue_convention" not in pairs:
            bad("orbital_pairs of type 'nocv' must state eigenvalue_convention")
        seen_pairs: dict[tuple, int] = {}
        for k, pair in enumerate(_list(pairs.get("pairs"))):
            if not isinstance(pair, dict):
                continue
            key = (pair.get("index"), pair.get("spin", "restricted"))
            if key in seen_pairs:
                bad(f"orbital_pairs.pairs[{k}] repeats index {key[0]} of spin {key[1]!r} "
                    f"(pairs[{seen_pairs[key]}])")
            else:
                seen_pairs[key] = k
            for field in ("donor", "acceptor"):
                f = pair.get(field)
                if _is_int(f) and not 0 <= f < n_frag:
                    bad(f"orbital_pairs.pairs[{k}].{field} names fragment {f}, but there "
                        f"are {n_frag} fragments")
            for field in ("kpoint_eigenvalues", "kpoint_energies"):
                if field in pair:
                    if n_k is None:
                        bad(f"orbital_pairs.pairs[{k}].{field} needs periodic.kpoints")
                    elif isinstance(pair[field], list) and len(pair[field]) != n_k:
                        bad(f"orbital_pairs.pairs[{k}].{field} has {len(pair[field])} "
                            f"entries for {n_k} k-points")
    return errors
