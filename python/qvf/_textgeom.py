"""Helpers shared by the text geometry converters of _espresso and _aims."""
from __future__ import annotations

import re

from ._legacy_writer import _ELEMENT_SYMBOLS
from .convert import _atom, _lattice_inverse, _number


def element_of_label(label, fmt):
    """The element of a species label such as Si, Fe1, Fe_up or O-d: a one
    or two letter element symbol (any case) followed by nothing or by a
    character that is not a letter."""
    for size in (2, 1):
        head, rest = label[:size], label[size:]
        symbol = head[:1].upper() + head[1:].lower()
        if len(head) == size and symbol in _ELEMENT_SYMBOLS[1:] \
                and (not rest or not rest[0].isalpha()):
            return symbol
    raise ValueError(f"{fmt} species label {label!r} does not start with an element symbol")


def cartesian(fractional, lattice):
    return [sum(fractional[i] * lattice[i][j] for i in range(3)) for j in range(3)]


def write_structure(writer, atoms, lattice, label, frozen):
    """One structure section: periodic on all three axes when ``lattice``
    is given, else a molecule."""
    section = {"id": "structure", "kind": "structure"}
    if label:
        section["label"] = label
    writer.add_section(section)
    payload = {"atoms": [_atom(s, [p + 0.0 for p in x]) for s, x in atoms],
               "pbc": [lattice is not None] * 3}
    if lattice is not None:
        payload["lattice_vectors"] = lattice
    writer.add_json("structure", "structure", payload)
    if frozen:
        writer.set_constraints(frozen_atoms=sorted(frozen))


def volume(lattice):
    a, b, c = lattice
    return abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
               + a[2] * (b[0] * c[1] - b[1] * c[0]))


def export_geometry(reader, selected, fmt, *, periodic_only):
    """(atoms [(symbol, position)], lattice rows or None, frozen indices,
    constraints not written, label) of a structure for a text export.

    ``periodic_only`` refuses a molecule; a structure periodic on only some
    axes is always refused."""
    geometry = reader.json(selected.id, "structure")
    pbc = geometry.get("pbc", [False] * 3)
    lattice = geometry.get("lattice_vectors")
    if all(pbc) and lattice is not None:
        rows = [[_number(v) for v in row] for row in lattice]
        _lattice_inverse(rows)
    elif not any(pbc) and not periodic_only:
        rows = None
    else:
        raise ValueError(f"{fmt} export requires a structure periodic on all three axes"
                         + ("" if periodic_only else " or a molecule"))
    atoms = []
    for atom in geometry["atoms"]:
        if atom["symbol"] not in _ELEMENT_SYMBOLS[1:]:
            raise ValueError(f"{fmt} export requires element symbols, not {atom['symbol']!r}")
        atoms.append((atom["symbol"], [_number(v) for v in atom["position"]]))
    if not atoms:
        raise ValueError(f"{fmt} export requires at least one atom")
    first = next(s.id for s in reader.sections if s.kind == "structure")
    constraints = reader.raw_manifest.get("constraints") or {}
    frozen = set(constraints.get("frozen_atoms", [])) if selected.id == first else set()
    if any(not 0 <= i < len(atoms) for i in frozen):
        raise ValueError("constraints.frozen_atoms names an atom outside the structure")
    omitted = []
    if constraints.get("frozen_atoms") and not frozen:
        omitted.append("frozen_atoms (they index the archive's first structure)")
    if constraints.get("distance_constraints"):
        omitted.append("distance_constraints")
    raw = next(s for s in reader.raw_manifest["sections"] if s["id"] == selected.id)
    label = " ".join(str(raw.get("label") or selected.id).split())
    label = re.sub(r"[^ -~]", "?", label)
    return atoms, rows, frozen, omitted, label


def fixed(value):
    """A fixed-point number that Fortran list-directed and QE card readers
    accept, with 15 decimals."""
    return f"{value + 0.0:.15f}"
