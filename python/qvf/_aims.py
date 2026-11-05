"""FHI-aims geometry.in interchange; see docs/converters.md.

Import reads ``lattice_vector``, ``atom`` and ``atom_frac`` lines (angstrom)
and relaxation constraints of one geometry.in; per-atom moments, charges and
velocities are counted, not stored. Export writes ``lattice_vector`` and
Cartesian ``atom`` lines in archive order.
"""
from __future__ import annotations

import re

from . import __version__
from ._textgeom import cartesian, element_of_label, export_geometry, fixed, volume, \
    write_structure
from .convert import _export_report, _new_output, _number, _select, _text, _validated
from .reader import Reader
from .writer import Writer

_REAL = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?")
# Per-atom keywords that follow an atom line and are reported, not stored.
_PER_ATOM = {"initial_moment", "initial_charge", "velocity", "constraint_region",
             "magnetic_response", "magnetic_moment", "nuclear_spin", "isotope",
             "RT_TDDFT_initial_velocity"}
# Global keywords that are reported, not stored.
_GLOBAL = {"hessian_block", "hessian_block_lv", "hessian_block_lv_atom", "hessian_file",
           "trust_radius", "symmetry_n_params", "symmetry_params", "symmetry_lv",
           "symmetry_frac", "homogeneous_field", "multipole", "set_vacuum_level",
           "esp_constraint", "verbatim_writeout", "calculate_friction"}


def _vector(tokens, label):
    if len(tokens) < 3 or not all(_REAL.fullmatch(t) for t in tokens[:3]):
        raise ValueError(f"geometry.in {label} requires three numbers")
    return [_number(t) for t in tokens[:3]]


def _parse(text):
    if any(ord(c) < 32 and c not in "\t\r\n" for c in text):
        raise ValueError("geometry.in requires text without control characters")
    lattice, sites, omitted, comments = [], [], {}, 0
    for line in text.splitlines():
        content = line.split("#", 1)[0].strip()
        comments += "#" in line
        if not content:
            continue
        tokens = content.split()
        key = tokens[0]
        if key in ("atom", "atom_frac"):
            if len(tokens) != 5:
                raise ValueError(f"geometry.in {key} requires three coordinates and a species")
            sites.append({"fractional": key == "atom_frac", "x": _vector(tokens[1:4], key),
                          "label": tokens[4],
                          "symbol": element_of_label(tokens[4], "geometry.in"),
                          "fixed": set()})
        elif key == "lattice_vector":
            if len(tokens) != 4:
                raise ValueError("geometry.in lattice_vector requires three numbers")
            lattice.append(_vector(tokens[1:], key))
        elif key == "constrain_relaxation":
            if not sites:
                raise ValueError("geometry.in constrain_relaxation must follow an atom")
            if len(tokens) != 2 or tokens[1].lower() not in (".true.", ".false.", "x", "y", "z"):
                raise ValueError("geometry.in constrain_relaxation requires .true., .false., "
                                 "x, y or z")
            flag = tokens[1].lower()
            if flag == ".true.":
                sites[-1]["fixed"] |= {"x", "y", "z"}
            elif flag != ".false.":
                sites[-1]["fixed"].add(flag)
        elif key in _PER_ATOM:
            if not sites:
                raise ValueError(f"geometry.in {key} must follow an atom")
            omitted[key] = omitted.get(key, 0) + 1
        elif key in _GLOBAL:
            omitted[key] = omitted.get(key, 0) + 1
        elif key == "empty":
            if len(tokens) != 5:
                raise ValueError("geometry.in empty requires three coordinates and a species")
            omitted["empty"] = omitted.get("empty", 0) + 1
        else:
            raise ValueError(f"geometry.in keyword {key!r} is not supported")
    if len(lattice) not in (0, 3):
        raise ValueError("geometry.in requires no lattice_vector (a molecule) or three")
    if lattice and volume(lattice) < 1e-12:
        raise ValueError("geometry.in lattice vectors are linearly dependent")
    if not sites:
        raise ValueError("geometry.in requires at least one atom")
    if any(s["fractional"] for s in sites) and not lattice:
        raise ValueError("geometry.in atom_frac requires three lattice_vector lines")
    atoms = [(s["symbol"], cartesian(s["x"], lattice) if s["fractional"] else s["x"])
             for s in sites]
    frozen = [k for k, s in enumerate(sites) if len(s["fixed"]) == 3]
    partial = sum(0 < len(s["fixed"]) < 3 for s in sites)
    labels = {s["label"]: s["symbol"] for s in sites}
    return {"atoms": atoms, "lattice": lattice or None, "frozen": frozen, "partial": partial,
            "omitted": omitted, "labels": labels, "comments": comments,
            "fractional": sum(s["fractional"] for s in sites)}


def import_aims(source, destination):
    head = _parse(_text(source))
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "FHI-aims geometry import") \
                as writer:
            write_structure(writer, head["atoms"], head["lattice"], None, head["frozen"])
        findings = _validated(temporary)
    return {"format": "aims", "output": str(destination), "atoms": len(head["atoms"]),
            "periodic": head["lattice"] is not None,
            "cell_volume_angstrom3": volume(head["lattice"]) if head["lattice"] else None,
            "fractional_atoms": head["fractional"], "species_labels": head["labels"],
            "frozen_atoms": len(head["frozen"]),
            "partially_fixed_atoms_omitted": head["partial"],
            "keywords_omitted": dict(sorted(head["omitted"].items())),
            "comment_lines": head["comments"], "findings": findings,
            "scope": "Lattice, species elements, positions and atoms fixed on all three axes "
                     "(as constraints.frozen_atoms) only; initial moments and charges, "
                     "velocities, empty sites, symmetry and Hessian blocks are not stored. Keep "
                     "the original geometry.in."}


def export_aims(source, destination, *, section=None):
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, {"structure"})
        atoms, lattice, frozen, omitted, label = export_geometry(
            reader, selected, "geometry.in", periodic_only=False)
        lines = [f"# {label}"]
        lines += ["lattice_vector " + " ".join(f"{fixed(v):>22s}" for v in row)
                  for row in lattice or []]
        for k, (symbol, x) in enumerate(atoms):
            lines.append("atom " + " ".join(f"{fixed(v):>22s}" for v in x) + f" {symbol}")
            if k in frozen:
                lines.append("    constrain_relaxation .true.")
        with _new_output(destination) as temporary:
            temporary.write_text("\n".join(lines) + "\n", encoding="ascii")
        report = _export_report(reader, selected, {selected.id}, destination, "aims")
        report.update(atoms=len(atoms), periodic=lattice is not None, frozen_atoms=len(frozen),
                      constraints_not_written=omitted,
                      generated="The section label as a comment, lattice_vector rows and "
                                "Cartesian atom lines in angstrom with 15 decimals in archive "
                                "order, species named by element, constrain_relaxation .true. "
                                "for frozen atoms.")
        return report
