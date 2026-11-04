"""Explicit, standard-library-only XYZ/extXYZ, PDB, CIF, mmCIF, POSCAR, CHGCAR, Quantum ESPRESSO,
FHI-aims, XSF, cube and Molden conversions.

Conversions create a new file and return a report. They do not execute a job
or infer a grid's physical quantity, periodicity, charge or multiplicity.
See docs/converters.md for the supported subsets and loss reporting.
"""
from __future__ import annotations

import math
from pathlib import Path
import re
import shlex
import struct

from . import __version__
from ._legacy_writer import BOHR_PER_ANGSTROM, _ELEMENT_SYMBOLS
from .units import HARTREE_IN_EV
from .operations import _new_output
from .reader import Reader
from .writer import Writer
from .validate import run

__all__ = ["import_xyz", "export_xyz", "import_cube", "export_cube", "import_xsf", "export_xsf",
           "import_pdb", "export_pdb", "import_molden", "export_molden", "import_cif",
           "export_cif", "import_poscar", "export_poscar", "import_chgcar", "export_chgcar",
           "import_mmcif", "export_mmcif", "import_espresso_in", "export_espresso_in",
           "import_espresso_xml", "import_espresso_cube", "import_aims", "export_aims"]
_MAX_TEXT = 64 * 1024 * 1024
_MAX_VALUES = 5_000_000
_SCALAR_KINDS = {"volume." + name for name in (
    "density", "orbital", "spin", "elf", "difference", "potential", "rdg",
    "generic", "laplacian", "transition_density")}


def import_molden(source, destination, *, spin, dialect="standard"):
    """Import molecular real s to g orbitals with an explicit spin and dialect.

    ``spin`` is restricted, restricted_open or unrestricted; it cannot be
    inferred from an Alpha-only Molden file. ``dialect`` names the producer
    convention (standard, renormalize, orca, turbomole, cfour, psi4-1.3 or
    xtb); it is never guessed. See docs/converters.md.
    """
    from ._molden import import_molden as convert
    return convert(source, destination, spin=spin, dialect=dialect)


def export_molden(source, destination, *, section=None):
    """Export one molecular real s/p/d wavefunction; report omitted metadata."""
    from ._molden import export_molden as convert
    return convert(source, destination, section=section)


def import_cif(source, destination, *, block=None):
    """Import one core-CIF data block as the periodic structure of its unit cell.

    The block is the only one with atom sites, or the one named by
    ``block``. The listed symmetry operations (else those of the Hall
    symbol; none are needed for P 1) are applied to the fractional sites,
    copies of a site within 0.1 A are merged, and a structure.symmetry
    section keeps the space group. Elements come only from
    ``_atom_site_type_symbol``. See docs/converters.md.
    """
    from ._cif import import_cif as convert
    return convert(source, destination, block=block)


def export_cif(source, destination, *, section=None):
    """Export one structure periodic on all three axes as a P 1 CIF."""
    from ._cif import export_cif as convert
    return convert(source, destination, section=section)


def import_poscar(source, destination, *, species=None):
    """Import a VASP POSCAR/CONTCAR as one periodic structure.

    Elements come from the VASP 5 species line, or from ``species`` for a
    VASP 4 file. Atoms fixed on all three axes by Selective dynamics become
    the root ``constraints.frozen_atoms``; velocities are reported, not
    stored. See docs/converters.md.
    """
    from ._poscar import import_poscar as convert
    return convert(source, destination, species=species)


def export_poscar(source, destination, *, section=None):
    """Export one structure periodic on all three axes as a VASP 5 POSCAR."""
    from ._poscar import export_poscar as convert
    return convert(source, destination, section=section)


def import_mmcif(source, destination, *, block=None):
    """Import one mmCIF (PDBx) data block's atom sites as structures.

    Atoms keep the PDB atom-site identity from the author columns; models,
    the cell, helices, sheet strands and struct_conn bonds follow the PDB
    import. See docs/converters.md.
    """
    from ._mmcif import import_mmcif as convert
    return convert(source, destination, block=block)


def export_mmcif(source, destination, *, section=None):
    """Export one structure as an mmCIF _atom_site loop (with _cell when periodic)."""
    from ._mmcif import export_mmcif as convert
    return convert(source, destination, section=section)


def import_chgcar(source, destination, *, species=None):
    """Import a VASP CHGCAR, CHG, AECCAR or PARCHG charge density.

    The POSCAR head becomes the periodic structure (as in import_poscar);
    the first grid, rho * V_cell, becomes ``volume.density`` in e / bohr^3 on
    a periodic grid without the endpoint, and a collinear or noncollinear
    magnetization becomes ``volume.spin`` or ``volume.vector``. PAW
    augmentation occupancies are skipped. See docs/converters.md.
    """
    from ._chgcar import import_chgcar as convert
    return convert(source, destination, species=species)


def export_chgcar(source, destination, *, section=None):
    """Export one periodic ``volume.density`` (and a magnetization on its grid) as CHGCAR."""
    from ._chgcar import export_chgcar as convert
    return convert(source, destination, section=section)


def import_espresso_in(source, destination):
    """Import the structure of a Quantum ESPRESSO pw.x input.

    Any ibrav of INPUT_PW with celldm or A, B, C, cosAB, cosAC, cosBC, or
    ibrav = 0 with CELL_PARAMETERS in alat, bohr or angstrom; ATOMIC_POSITIONS
    in alat, bohr, angstrom or crystal. Elements come from the species labels;
    atoms with all three if_pos flags 0 become ``constraints.frozen_atoms``.
    See docs/converters.md.
    """
    from ._espresso import import_espresso_in as convert
    return convert(source, destination)


def export_espresso_in(source, destination, *, section=None):
    """Export one structure periodic on all three axes as a pw.x input with ibrav = 0."""
    from ._espresso import export_espresso_in as convert
    return convert(source, destination, section=section)


def import_espresso_xml(source, destination):
    """Import the final structure, total energy, forces and stress of a pw.x XML data file.

    Lengths are bohr and energies hartree in the file; the stress sign is
    flipped to the QVF convention. Results of an SCF that did not converge
    are omitted. See docs/converters.md.
    """
    from ._espresso import import_espresso_xml as convert
    return convert(source, destination)


def import_espresso_cube(source, destination, *, kind="volume.density"):
    """Import a pp.x cube file as a periodic grid without the endpoint plane."""
    from ._espresso import import_espresso_cube as convert
    return convert(source, destination, kind=kind)


def import_aims(source, destination):
    """Import an FHI-aims geometry.in as one periodic structure or a molecule.

    ``lattice_vector``, ``atom`` and ``atom_frac`` are read in angstrom; atoms
    with ``constrain_relaxation .true.`` (or x, y and z) become
    ``constraints.frozen_atoms``. See docs/converters.md.
    """
    from ._aims import import_aims as convert
    return convert(source, destination)


def export_aims(source, destination, *, section=None):
    """Export one periodic structure or molecule as an FHI-aims geometry.in."""
    from ._aims import export_aims as convert
    return convert(source, destination, section=section)


def _text(path):
    with open(path, "rb") as stream:
        raw = stream.read(_MAX_TEXT + 1)
    if len(raw) > _MAX_TEXT:
        raise ValueError("converter text limit exceeded (64 MiB)")
    return raw.decode("utf-8")


def _number(value):
    result = float(str(value).replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError("conversion requires finite numeric values")
    return result


def _symbol(value):
    if value not in _ELEMENT_SYMBOLS[1:]:
        raise ValueError(f"unsupported element symbol {value!r}")
    return value


def _atom(symbol, position):
    return {"symbol": _symbol(symbol), "atomic_number": _ELEMENT_SYMBOLS.index(symbol),
            "position": [_number(v) for v in position]}


def _binary(writer, sid, role, values, shape):
    with writer.begin_binary(sid, role, dtype="float64", shape=shape) as sink:
        chunk = []
        for value in values:
            chunk.append(_number(value))
            if len(chunk) == 4096:
                sink.write(struct.pack("<4096d", *chunk)); chunk.clear()
        if chunk:
            sink.write(struct.pack("<" + str(len(chunk)) + "d", *chunk))


def _check_cell(pbc, lattice):
    if any(pbc) and lattice is None:
        raise ValueError("periodic geometry needs a cell")
    if lattice is not None:
        active = [row for row, flag in zip(lattice, pbc) if flag]
        if any(not any(v != 0 for v in row) for row in active):
            raise ValueError("periodic cell vectors must be nonzero")
        if len(active) >= 2:
            a, b = active[:2]
            cross = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
            if not any(v != 0 for v in cross) or (len(active) == 3 and sum(x*y for x, y in zip(cross, active[2])) == 0):
                raise ValueError("periodic cell vectors must be independent")


def _structure(writer, atoms, pbc, lattice):
    _check_cell(pbc, lattice)
    writer.add_section({"id": "structure", "kind": "structure"})
    writer.add_json("structure", "structure", {
        "atoms": atoms, "pbc": pbc, "lattice_vectors": lattice})


def _validated(path):
    report = run(path)
    if not report["ok"]:
        errors = [f["code"] for f in report["findings"] if f["severity"] == "error"]
        raise ValueError("converted archive failed validation: " + ", ".join(errors))
    return report["findings"]


_EV_PER_ANGSTROM = 1 / (HARTREE_IN_EV * BOHR_PER_ANGSTROM)       # to hartree/bohr
_EV_PER_ANGSTROM3 = 1 / (HARTREE_IN_EV * BOHR_PER_ANGSTROM ** 3)  # to hartree/bohr^3


def _xyz_properties(value):
    """[(name, type, count, first column)] and the column count of an extXYZ
    Properties value such as species:S:1:pos:R:3:forces:R:3."""
    parts = value.split(":")
    if len(parts) % 3:
        raise ValueError("extXYZ Properties needs name:type:count triplets")
    props, column = [], 0
    for name, kind, count in zip(parts[::3], parts[1::3], parts[2::3]):
        if kind not in ("S", "R", "I", "L") or not re.fullmatch(r"[1-9][0-9]*", count):
            raise ValueError(f"extXYZ property {name!r} needs a type S, R, I or L and a count")
        props.append((name, kind, int(count), column))
        column += int(count)
    names = [p[0] for p in props]
    if len(set(names)) != len(names):
        raise ValueError("extXYZ Properties repeats a name")
    shapes = {p[0]: p[1:3] for p in props}
    for name, shape in (("species", ("S", 1)), ("pos", ("R", 3))):
        if shapes.get(name) != shape:
            raise ValueError("extXYZ requires species:S:1 and pos:R:3")
    for name, shape in (("forces", ("R", 3)), ("Z", ("I", 1))):
        if name in shapes and shapes[name] != shape:
            raise ValueError(f"extXYZ {name} must be {shape[0]}:{shape[1]}")
    return props, column


def _xyz_frames(text):
    lines = text.splitlines(); index = 0; frames = []; count = 0; comments = 0
    while index < len(lines):
        if not lines[index].strip() and not any(line.strip() for line in lines[index:]):
            break
        n = int(lines[index]); index += 1
        count += n * 3
        if n < 1 or count > _MAX_VALUES or index + n >= len(lines):
            raise ValueError("invalid, truncated or oversized XYZ frame")
        comment = lines[index]; index += 1
        lattice = None; pbc = [False] * 3
        props, columns = [("species", "S", 1, 0), ("pos", "R", 3, 1)], 4
        energy = stress = None; omitted_keys = []
        if re.search(r'\b(?:Properties|Lattice|pbc|energy|stress)\s*=', comment):
            attrs = {}
            # key=value pairs, with quoted space-separated vectors.
            try:
                tokens = shlex.split(re.sub(r'\s*=\s*', '=', comment), comments=False)
            except ValueError:
                raise ValueError("extXYZ header has an unclosed quotation") from None
            for token in tokens:
                key, _, value = token.partition('=')
                if key in attrs:
                    raise ValueError(f"repeated extXYZ field {key!r}")
                attrs[key] = value
            omitted_keys = sorted(k for k in attrs
                                  if k not in {"Properties", "Lattice", "pbc", "energy", "stress"})
            if "Properties" in attrs:
                props, columns = _xyz_properties(attrs["Properties"])
            if "Lattice" in attrs:
                values = [_number(x) for x in attrs["Lattice"].split()]
                if len(values) != 9:
                    raise ValueError("extXYZ Lattice must contain nine row-vector components")
                lattice = [values[i:i + 3] for i in (0, 3, 6)]
                pbc = [True] * 3
            if "pbc" in attrs:
                flags = attrs["pbc"].lower().split()
                if len(flags) != 3 or any(x not in {"t", "f", "true", "false"} for x in flags):
                    raise ValueError("extXYZ pbc must contain three booleans")
                pbc = [x in {"t", "true"} for x in flags]
            if any(pbc) and lattice is None:
                raise ValueError("periodic extXYZ needs Lattice")
            if "energy" in attrs:
                energy = _number(attrs["energy"])
            if "stress" in attrs:
                values = [_number(x) for x in attrs["stress"].split()]
                if len(values) != 9:
                    raise ValueError("extXYZ stress must contain nine components")
                stress = [values[i:i + 3] for i in (0, 3, 6)]
        elif comment:
            comments += 1
        column = {p[0]: p for p in props}
        atoms, forces = [], [] if "forces" in column else None
        for line in lines[index:index + n]:
            fields = line.split()
            if len(fields) != columns:
                raise ValueError("extXYZ atom rows must have one field per Properties column")
            species, pos = column["species"][3], column["pos"][3]
            atom = _atom(fields[species], fields[pos:pos + 3])
            if "Z" in column and (not re.fullmatch(r"[0-9]+", fields[column["Z"][3]])
                                  or int(fields[column["Z"][3]]) != atom["atomic_number"]):
                raise ValueError("extXYZ Z disagrees with the species")
            if forces is not None:
                start = column["forces"][3]
                forces.append([_number(v) for v in fields[start:start + 3]])
            atoms.append(atom)
        index += n
        frames.append({"atoms": atoms, "pbc": pbc, "lattice_vectors": lattice,
                       "energy": energy, "forces": forces, "stress": stress,
                       "omitted_keys": omitted_keys,
                       "omitted_properties": [p[0] for p in props
                                              if p[0] not in ("species", "pos", "forces", "Z")]})
    if not frames:
        raise ValueError("XYZ has no frames")
    first = frames[0]
    for frame in frames[1:]:
        if [a["symbol"] for a in frame["atoms"]] != [a["symbol"] for a in first["atoms"]]:
            raise ValueError("trajectory frames must preserve atom count and symbol order")
        if frame["pbc"] != first["pbc"] or (frame["lattice_vectors"] is None) != (first["lattice_vectors"] is None):
            raise ValueError("trajectory frames must preserve periodicity and cell availability")
    return frames, comments


def import_xyz(source, destination):
    """Convert XYZ/extXYZ to a structure, plus a trajectory for several frames.

    extXYZ ``energy`` (eV), ``forces`` (eV/A) and ``stress`` (eV/A^3) become
    the energy, the gradient and the stress in atomic units; other header
    keys and per-atom columns are reported as omitted.
    """
    frames, comments = _xyz_frames(_text(source)); first = frames[0]
    energies = all(f["energy"] is not None for f in frames)
    forces = all(f["forces"] is not None for f in frames)
    single = len(frames) == 1
    stress = single and first["stress"] is not None and any(first["pbc"])
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "XYZ geometry import") as writer:
            _structure(writer, first["atoms"], first["pbc"], first["lattice_vectors"])
            if single and energies:
                writer.set_energies(first["energy"] / HARTREE_IN_EV)
            if single and forces:
                writer.add_section({"id": "forces", "kind": "forces", "structure_ref": "structure"})
                _binary(writer, "forces", "gradient",
                        (-v * _EV_PER_ANGSTROM for row in first["forces"] for v in row),
                        [len(first["atoms"]), 3])
                if stress:
                    _binary(writer, "forces", "stress",
                            (v * _EV_PER_ANGSTROM3 for row in first["stress"] for v in row), [3, 3])
            if len(frames) > 1:
                writer.add_section({"id": "trajectory", "kind": "trajectory", "structure_ref": "structure"})
                metadata = {"symbols": [a["symbol"] for a in first["atoms"]], "pbc": first["pbc"]}
                if energies:
                    metadata["energies"] = [f["energy"] / HARTREE_IN_EV for f in frames]
                writer.add_json("trajectory", "metadata", metadata)
                if forces:
                    _binary(writer, "trajectory", "gradient",
                            (-v * _EV_PER_ANGSTROM for f in frames for row in f["forces"]
                             for v in row), [len(frames), len(first["atoms"]), 3])
                _binary(writer, "trajectory", "coords", (v for f in frames for a in f["atoms"] for v in a["position"]),
                        [len(frames), len(first["atoms"]), 3])
                if first["lattice_vectors"] is not None:
                    _binary(writer, "trajectory", "lattice", (v for f in frames for row in f["lattice_vectors"] for v in row),
                            [len(frames), 3, 3])
        findings = _validated(temporary)
    with_energy = sum(f["energy"] is not None for f in frames)
    with_forces = sum(f["forces"] is not None for f in frames)
    with_stress = sum(f["stress"] is not None for f in frames)
    return {"format": "xyz", "frames": len(frames), "findings": findings,
            "omitted_comments": comments, "output": str(destination),
            "energies": len(frames) if energies else 0, "forces": len(frames) if forces else 0,
            "stress": int(stress),
            "omitted": {"energies": with_energy if not energies else 0,
                        "forces": with_forces if not forces else 0,
                        "stress": with_stress - int(stress),
                        "header_keys": sorted({k for f in frames for k in f["omitted_keys"]}),
                        "properties": sorted({p for f in frames
                                              for p in f["omitted_properties"]})}}


def _assigned_forces(reader, structure_id):
    """The one forces section assigned to ``structure_id`` (section-level or
    meta structure_ref, else the archive's only structure), else None."""
    raw = {x["id"]: x for x in reader.raw_manifest["sections"]}
    structures = [x.id for x in reader.sections if x.kind == "structure"]
    found = []
    for x in reader.sections:
        if x.kind == "forces":
            meta = reader.json(x.id, "meta") if x.member("meta") else {}
            ref = raw[x.id].get("structure_ref") or meta.get("structure_ref")
            if ref == structure_id or (ref is None and len(structures) == 1):
                found.append(x)
    return found[0] if len(found) == 1 else None


def _select(reader, section, kinds):
    if section is not None:
        result = reader.section(section)
        if result.kind not in kinds:
            raise ValueError("selected section has an unsupported conversion kind")
        return result
    choices = [s for s in reader.sections if s.kind in kinds]
    if len(choices) != 1:
        raise ValueError("select a section explicitly when zero or multiple candidates exist")
    return choices[0]


def _export_report(reader, selected, used, destination, fmt):
    return {"format": fmt, "section": selected.id, "output": str(destination),
            "omitted_sections": [s.id for s in reader.sections if s.id not in used],
            "scope": "Selected geometry/grid only; QVF metadata, bonds, identity, jobs and other results are not preserved. Keep the original archive."}


_PDB_METADATA = set("""HEADER OBSLTE TITLE SPLIT CAVEAT COMPND SOURCE KEYWDS
    EXPDTA MDLTYP AUTHOR REVDAT SPRSDE JRNL REMARK DBREF DBREF1 DBREF2 SEQADV
    SEQRES MODRES HET HETNAM HETSYN FORMUL CISPEP SITE
    ORIGX1 ORIGX2 ORIGX3
    ANISOU TER MASTER""".split())


def _pdb_integer(field, label):
    if not re.fullmatch(r"[0-9]+", field.strip()) or int(field) < 1:
        raise ValueError(f"PDB {label} requires a positive decimal integer")
    return int(field)


def _pdb_decimal(field, label, low, high):
    """An optional fixed-column decimal such as occupancy (blank: absent)."""
    if not field.strip():
        return None
    if not re.fullmatch(r"[+-]?[0-9]+\.[0-9]+", field.strip()):
        raise ValueError(f"PDB {label} requires a fixed-point decimal")
    value = float(field)
    if not low <= value <= high:
        raise ValueError(f"PDB {label} must lie in [{low}, {high}]")
    return value


def _pdb_atom(line, record):
    """One ATOM/HETATM record as a structure atom with its atom-site identity
    (spec § 5.1). Blank columns are omitted, never written as placeholders."""
    fields = [line[i:i + 8].strip() for i in (30, 38, 46)]
    if any(not re.fullmatch(r"[+-]?[0-9]+\.[0-9]{3}", v) for v in fields):
        raise ValueError("PDB coordinates require fixed-width Real(8.3) fields")
    # Atom names such as CA are not element declarations. Never guess.
    result = _atom(line[76:78].strip().capitalize(), fields)
    seq = line[22:26].strip()
    if seq and not re.fullmatch(r"-?[0-9]+", seq):
        raise ValueError("PDB residue sequence numbers must be decimal (no hybrid-36)")
    identity = {"atom_name": line[12:16].strip(), "residue_name": line[17:20].strip(),
                "residue_seq": int(seq) if seq else None, "chain_id": line[21].strip(),
                "insertion_code": line[26].strip(), "alt_loc": line[16].strip(),
                "occupancy": _pdb_decimal(line[54:60], "occupancy", 0, 1),
                "b_factor": _pdb_decimal(line[60:66], "temperature factor", -1e6, 1e6)}
    result.update({k: v for k, v in identity.items() if v not in (None, "")})
    result["hetero"] = record == "HETATM"
    return result


def _pdb_seq(field, label):
    if not re.fullmatch(r"-?[0-9]+", field.strip()):
        raise ValueError(f"PDB {label} requires a decimal residue number (no hybrid-36)")
    return int(field)


def _pdb_range(line, record):
    """One HELIX or SHEET record as a ``secondary_structure`` range (spec
    § 5.1, Secondary-structure ranges); blank optional columns are absent."""
    line = line.ljust(80)
    if record == "HELIX":
        start, end = (line[19], line[21:25], line[25]), (line[31], line[33:37], line[37])
        result = {"type": "helix"}
        if line[38:40].strip():
            helix_class = _pdb_integer(line[38:40], "helix class")
            if helix_class > 10:
                raise ValueError("PDB helix class must lie in 1-10")
            result["helix_class"] = helix_class
    else:
        start, end = (line[21], line[22:26], line[26]), (line[32], line[33:37], line[37])
        result = {"type": "sheet", "strand": _pdb_integer(line[7:10], "SHEET strand")}
        if line[11:14].strip():
            result["sheet_id"] = line[11:14].strip()
    if start[0] != end[0]:
        raise ValueError(f"PDB {record} ranges must start and end on one chain")
    result.update(chain=start[0].strip(), start_seq=_pdb_seq(start[1], record),
                  end_seq=_pdb_seq(end[1], record))
    for key, code in (("start_insertion_code", start[2]), ("end_insertion_code", end[2])):
        if code.strip():
            result[key] = code.strip()
    order = ("type", "chain", "start_seq", "start_insertion_code", "end_seq",
             "end_insertion_code", "helix_class", "sheet_id", "strand")
    return {k: result[k] for k in order if k in result}


def _pdb_unresolved(ranges, atoms):
    """How many ranges have an end that is not a residue of ``atoms``."""
    present = {(a.get("chain_id", ""), a["residue_seq"], a.get("insertion_code"))
               for a in atoms if "residue_seq" in a}
    return sum((r["chain"], r["start_seq"], r.get("start_insertion_code")) not in present
               or (r["chain"], r["end_seq"], r.get("end_insertion_code")) not in present
               for r in ranges)


def _pdb_residues(atoms):
    """(chains, residues) in file order, keyed by (chain, seq, insertion
    code); None when an atom lacks a residue name or number."""
    if any("residue_name" not in a or "residue_seq" not in a for a in atoms):
        return None
    chains, residues = [], {}
    for i, a in enumerate(atoms):
        chain = a.get("chain_id", "")
        if chain and chain not in chains:
            chains.append(chain)
        key = (chain, a["residue_seq"], a.get("insertion_code"))
        if key not in residues:
            residues[key] = {"name": a["residue_name"], "seq": a["residue_seq"],
                             "chain": chain, "atom_indices": []}
            if key[2]:
                residues[key]["insertion_code"] = key[2]
        residues[key]["atom_indices"].append(i)
    return chains, list(residues.values())


def _pdb_conect(line):
    """(serial, bonded serials, obsolete hydrogen-bond/salt-bridge fields?)
    of one CONECT record."""
    line = line.ljust(80)
    serial = _pdb_integer(line[6:11], "CONECT serial")
    partners = [_pdb_integer(line[k:k + 5], "CONECT bonded serial")
                for k in (11, 16, 21, 26) if line[k:k + 5].strip()]
    if not partners:
        raise ValueError("PDB CONECT names no bonded atom")
    if serial in partners:
        raise ValueError("PDB CONECT bonds an atom to itself")
    return serial, partners, bool(line[31:61].strip())


def _pdb_bonds(conects, serials, label):
    """Distinct pairs of 0-based atom indices, sorted, and the number of
    repeated listings (a pair listed again from the same atom). Bond order
    is never inferred from repetition."""
    index = {s: i for i, s in enumerate(serials)}
    pairs, listed, repeated = set(), set(), 0
    for serial, partners, _ in conects:
        for partner in partners:
            for s in (serial, partner):
                if s not in index:
                    raise ValueError(f"PDB CONECT names atom serial {s}, which is not an "
                                     f"atom of {label}")
            if (serial, partner) in listed:
                repeated += 1
            listed.add((serial, partner))
            pairs.add(tuple(sorted((index[serial], index[partner]))))
    return sorted(pairs), repeated


def _pdb_link(line, record):
    """(record, two atom ends, symmetry operators) of one SSBOND or LINK
    record. An end is (chain, seq, insertion code, residue name, atom name,
    alternate location); SSBOND names the SG atoms of two cysteines. A blank
    operator means 1555, the identity."""
    line = line.ljust(80)
    if record == "SSBOND":
        ends = [(line[15], line[17:21], line[21], line[11:14], "SG", " "),
                (line[29], line[31:35], line[35], line[25:28], "SG", " ")]
    else:
        ends = [(line[21], line[22:26], line[26], line[17:20], line[12:16], line[16]),
                (line[51], line[52:56], line[56], line[47:50], line[42:46], line[46])]
    resolved = []
    for chain, seq, code, residue, atom, alt in ends:
        if not atom.strip() or not residue.strip():
            raise ValueError(f"PDB {record} requires residue and atom names")
        resolved.append((chain.strip(), _pdb_seq(seq, record), code.strip() or None,
                         residue.strip(), atom.strip(), alt.strip() or None))
    operators = tuple(line[k:k + 6].strip() or "1555" for k in (59, 66))
    if not all(re.fullmatch(r"[0-9]{4,6}", op) for op in operators):
        raise ValueError(f"PDB {record} symmetry operators must be codes such as 1555")
    return record, resolved, operators


def _pdb_link_bonds(links, atoms, label):
    """{record: distinct pairs} of the identity-operator SSBOND and LINK
    records, resolved by atom-site identity (chain, residue number, insertion
    code, residue and atom name). A blank alternate location names every
    alternate of the atom, a given one that alternate or an atom without
    alternates; alternates pair only with the same alternate or an atom
    without one. A missing, ambiguous or self-bonded atom refuses."""
    index = {}
    for i, a in enumerate(atoms):
        index.setdefault((a.get("chain_id", ""), a.get("residue_seq"), a.get("insertion_code"),
                          a.get("residue_name"), a.get("atom_name")), []).append(i)
    found = {"SSBOND": set(), "LINK": set()}
    for record, ends, operators in links:
        if operators != ("1555", "1555"):
            continue
        sides = []
        for chain, seq, code, residue, name, alt in ends:
            where = f"{name} of {residue} {chain or '-'}{seq}{code or ''}"
            candidates = index.get((chain, seq, code, residue, name), [])
            alternates = [atoms[i].get("alt_loc") for i in candidates]
            if len(set(alternates)) != len(alternates):
                raise ValueError(f"PDB {record} atom {where} is ambiguous in {label}")
            if alt is not None:
                candidates = ([i for i in candidates if atoms[i].get("alt_loc") == alt]
                              or [i for i in candidates if atoms[i].get("alt_loc") is None])
            if not candidates:
                raise ValueError(f"PDB {record} names atom {where}"
                                 f"{f' altloc {alt}' if alt else ''}, which is not an atom "
                                 f"of {label}")
            sides.append(candidates)
        pairs = {tuple(sorted((i, j))) for i in sides[0] for j in sides[1]
                 if None in (atoms[i].get("alt_loc"), atoms[j].get("alt_loc"))
                 or atoms[i].get("alt_loc") == atoms[j].get("alt_loc")}
        if not pairs:
            raise ValueError(f"PDB {record} joins atoms of no common alternate location")
        if any(i == j for i, j in pairs):
            raise ValueError(f"PDB {record} bonds an atom to itself")
        found[record] |= pairs
    return found


def _pdb_model_bonds(conects, links, model, label):
    """(sorted pairs, repeated CONECT listings, per-record counts) of one
    model: CONECT pairs merged with resolved SSBOND and LINK pairs."""
    pairs, repeated = _pdb_bonds(conects, model[2], label)
    linked = _pdb_link_bonds(links, model[1], label)
    listed = set(pairs)
    counts = {r.lower(): {"bonds": len(p), "also_conect": len(p & listed)}
              for r, p in linked.items()}
    return sorted(listed.union(*linked.values())), repeated, counts


def _pdb_cell_vectors(a, b, c, alpha, beta, gamma):
    """Lattice rows in the PDB standard orientation: a along x, b in the
    xy-plane, c completing a right-handed cell (wwPDB CRYST1/SCALE)."""
    ca, cb, cg = (math.cos(math.radians(v)) for v in (alpha, beta, gamma))
    sg = math.sin(math.radians(gamma))
    cy = (ca - cb * cg) / sg
    cz2 = 1 - cb * cb - cy * cy
    if sg <= 0 or cz2 <= 0:
        raise ValueError("PDB CRYST1 angles do not describe a cell")
    return [[a, 0.0, 0.0], [b * cg, b * sg, 0.0], [c * cb, c * cy, c * math.sqrt(cz2)]]


def _lattice_inverse(rows):
    """The inverse of a lattice with cell vectors as rows, so fractional
    coordinates are the row vector x times it."""
    det = (rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
           - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
           + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0]))
    if abs(det) < 1e-9:
        raise ValueError("the lattice is singular")
    return [[(rows[(j + 1) % 3][(i + 1) % 3] * rows[(j + 2) % 3][(i + 2) % 3]
              - rows[(j + 1) % 3][(i + 2) % 3] * rows[(j + 2) % 3][(i + 1) % 3]) / det
             for j in range(3)] for i in range(3)]


def _pdb_scale(lattice):
    """The SCALE matrix (fractional = SCALE x) of an upper-triangular-form
    lattice whose rows are the cell vectors."""
    (ax, _, _), (bx, by, _), (cx, cy, cz) = lattice
    return [[1 / ax, -bx / (ax * by), (bx * cy - by * cx) / (ax * by * cz)],
            [0.0, 1 / by, -cy / (by * cz)],
            [0.0, 0.0, 1 / cz]]


def _pdb_cryst1(line):
    """(cell lengths and angles, space group, Z) of a CRYST1 record."""
    line = line.ljust(80)
    values = [_pdb_decimal(line[i:j], "CRYST1 cell parameter", 1e-6, 1e6)
              for i, j in ((6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54))]
    if any(v is None for v in values) or not all(0 < v < 180 for v in values[3:]):
        raise ValueError("PDB CRYST1 requires three positive lengths and angles in (0, 180)")
    group = " ".join(line[55:66].split())
    if not group:
        raise ValueError("PDB CRYST1 requires a space group")
    z = _pdb_integer(line[66:70], "CRYST1 Z") if line[66:70].strip() else None
    return values, group, z


def _pdb_cell(cryst1, scales):
    """(lattice or None, space group, Z, placeholder) from CRYST1 and SCALEn.
    SCALE must be the standard orthogonalization of the cell with zero
    translation; the wwPDB placeholder 1 A, 90 degrees, P 1 means no cell."""
    if cryst1 is None:
        if scales:
            raise ValueError("PDB SCALE records require a CRYST1 cell")
        return None, None, None, False
    values, group, z = cryst1
    lattice = _pdb_cell_vectors(*values)
    if scales:
        if sorted(scales) != [1, 2, 3]:
            raise ValueError("PDB SCALE requires exactly SCALE1, SCALE2 and SCALE3")
        expected = _pdb_scale(lattice)
        for n, (row, shift) in scales.items():
            if abs(shift) > 1e-5 or any(abs(v - e) > 1e-4 * max(abs(v), abs(e)) + 2e-6
                                        for v, e in zip(row, expected[n - 1])):
                raise ValueError("PDB SCALE does not match the standard orthogonalization "
                                 "of the CRYST1 cell; the coordinates use another frame")
    placeholder = values == [1.0, 1.0, 1.0, 90.0, 90.0, 90.0] and group == "P 1" \
        and z in (None, 1)
    return (None if placeholder else lattice), group, z, placeholder


def _pdb_scale_record(line):
    line = line.ljust(80)
    n = _pdb_integer(line[5], "SCALE row")
    row = [_pdb_decimal(line[i:i + 10], "SCALE element", -1e6, 1e6) for i in (10, 20, 30)]
    shift = _pdb_decimal(line[45:55], "SCALE translation", -1e6, 1e6)
    if any(v is None for v in row) or shift is None:
        raise ValueError("PDB SCALE requires three matrix elements and a translation")
    return n, (row, shift)


def _pdb_operator(rows, label):
    """The 3x4 operator [R | t] of rows 1-3; R must be a proper rotation."""
    if sorted(rows) != [1, 2, 3]:
        raise ValueError(f"PDB {label} requires rows 1, 2 and 3, each once")
    matrix = [rows[n] for n in (1, 2, 3)]
    r = [row[:3] for row in matrix]
    det = (r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
           - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
           + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]))
    if det < 0 or any(abs(sum(r[i][k] * r[j][k] for k in range(3)) - (i == j)) > 1e-4
                      for i in range(3) for j in range(3)):
        raise ValueError(f"PDB {label} is not a proper rotation")
    return matrix


def _pdb_mtrix(lines):
    """(operators, given, not given, [(serial, operator)] not given) of
    MTRIX1-3 records. iGiven 1 means the copy's coordinates are in the file;
    blank means they are not."""
    rows = {}
    for line in lines:
        line = line.ljust(80)
        n, serial = int(line[5]), _pdb_integer(line[7:10], "MTRIX serial")
        values = [_pdb_decimal(line[i:i + 10], "MTRIX element", -1e6, 1e6)
                  for i in (10, 20, 30)] + [_pdb_decimal(line[45:55], "MTRIX translation",
                                                         -1e6, 1e6)]
        if None in values or line[59] not in " 1":
            raise ValueError("PDB MTRIX requires three elements, a translation and iGiven "
                             "1 or blank")
        entry = rows.setdefault(serial, ({}, set()))
        if n in entry[0]:
            raise ValueError(f"PDB MTRIX{n} {serial} is repeated")
        entry[0][n] = values
        entry[1].add(line[59] == "1")
    for serial, (matrix, given) in rows.items():
        _pdb_operator(matrix, f"MTRIX {serial}")
        if len(given) != 1:
            raise ValueError(f"PDB MTRIX {serial} rows disagree on iGiven")
    given = sum(next(iter(g)) for _, g in rows.values())
    missing = [(serial, _pdb_operator(matrix, f"MTRIX {serial}"))
               for serial, (matrix, g) in sorted(rows.items()) if not next(iter(g))]
    return len(rows), given, len(rows) - given, missing


def _pdb_assemblies(lines):
    """[(biomolecule, [(chains, {serial: operator})])] of REMARK 350: one
    block per APPLY THE FOLLOWING TO CHAINS, with its BIOMT operators in file
    order. Other REMARK 350 lines (such as the determination) are skipped."""
    assemblies = []
    for line in lines:
        text = line[10:].strip()
        blocks = assemblies[-1][1] if assemblies else None
        if text.startswith("BIOMOLECULE:"):
            ident = text.split(":", 1)[1].strip()
            if not re.fullmatch(r"[0-9]+", ident) or any(a[0] == ident for a in assemblies):
                raise ValueError("PDB REMARK 350 biomolecules need distinct decimal numbers")
            assemblies.append((ident, []))
        elif text.startswith(("APPLY THE FOLLOWING TO CHAINS:", "AND CHAINS:")):
            if blocks is None:
                raise ValueError("PDB REMARK 350 lists chains before BIOMOLECULE")
            if text.startswith("APPLY"):
                blocks.append(([], {}))
            elif not blocks or blocks[-1][1]:
                raise ValueError("PDB REMARK 350 AND CHAINS must continue a chain list")
            blocks[-1][0].extend(c.strip() for c in text.split(":", 1)[1].split(",")
                                 if c.strip())
        elif text.startswith("BIOMT"):
            tokens = text.split()
            if len(tokens) != 6 or tokens[0] not in ("BIOMT1", "BIOMT2", "BIOMT3"):
                raise ValueError("PDB REMARK 350 BIOMT requires a serial, three elements "
                                 "and a translation")
            if not blocks:
                raise ValueError("PDB REMARK 350 BIOMT must follow a chain list")
            serial = _pdb_integer(tokens[1], "BIOMT serial")
            rows = blocks[-1][1].setdefault(serial, {})
            if int(tokens[0][-1]) in rows:
                raise ValueError(f"PDB REMARK 350 {tokens[0]} {serial} is repeated")
            rows[int(tokens[0][-1])] = [_pdb_decimal(t, "BIOMT element", -1e6, 1e6)
                                        for t in tokens[2:]]
    for ident, blocks in assemblies:
        if not blocks or any(not chains or not ops for chains, ops in blocks):
            raise ValueError(f"PDB REMARK 350 biomolecule {ident} needs chains and BIOMT "
                             "operators")
        for k, (chains, ops) in enumerate(blocks):
            if len(set(chains)) != len(chains):
                raise ValueError(f"PDB REMARK 350 biomolecule {ident} repeats a chain")
            blocks[k] = (chains, {s: _pdb_operator(r, f"BIOMT {s} of biomolecule {ident}")
                                  for s, r in ops.items()})
    return assemblies


def _pdb_ncs(atoms, bonds, ranges, operators):
    """(atoms, bonds, ranges, copies) of the asymmetric unit completed by the
    MTRIX operators whose copies are not given: the deposited atoms, then
    each operator (x' = R x + t, identity skipped) applied to all of them,
    with chains named <chain>-<serial>."""
    out_atoms, out_bonds, out_ranges = list(atoms), list(bonds), list(ranges)
    copies = []
    for serial, op in operators:
        if all(abs(op[i][k] - (i == k)) < 1e-6 for i in range(3) for k in range(4)):
            continue
        offset = len(out_atoms)
        for a in atoms:
            x = a["position"]
            atom = dict(a, position=[round(sum(op[r][k] * x[k] for k in range(3)) + op[r][3], 6)
                                     + 0.0 for r in range(3)])
            atom["chain_id"] = f"{a.get('chain_id', '')}-{serial}"
            out_atoms.append(atom)
        out_bonds += [(i + offset, j + offset) for i, j in bonds]
        out_ranges += [{**r, "chain": f"{r['chain']}-{serial}"} for r in ranges]
        copies.append({"operator": serial, "matrix": op})
    return out_atoms, out_bonds, out_ranges, copies


def _pdb_assembly(blocks, atoms, bonds, ranges, label):
    """(atoms, bonds, ranges, copies) of one biological assembly: each BIOMT
    operator applied to the atoms of its chains. An identity copy keeps the
    chain id; another copy is named <chain>-<operator serial>. Bonds and
    secondary-structure ranges are copied within each copy."""
    present = {a.get("chain_id", "") for a in atoms}
    out_atoms, out_bonds, out_ranges, names, copies = [], [], [], set(), []
    for chains, operators in blocks:
        missing = [c for c in chains if c not in present]
        if missing:
            raise ValueError(f"PDB REMARK 350 names chain {missing[0]}, which has no atoms "
                             f"in {label}")
        chosen = [i for i, a in enumerate(atoms) if a.get("chain_id", "") in chains]
        for serial, op in operators.items():
            identity = all(abs(op[i][k] - (i == k)) < 1e-6 for i in range(3) for k in range(4))
            rename = {c: c if identity else f"{c}-{serial}" for c in chains}
            if names & set(rename.values()):
                raise ValueError(f"PDB REMARK 350 copies a chain twice under BIOMT {serial}")
            names |= set(rename.values())
            index = {i: len(out_atoms) + k for k, i in enumerate(chosen)}
            for i in chosen:
                atom = dict(atoms[i]); x = atom["position"]
                atom["position"] = [round(sum(op[r][k] * x[k] for k in range(3)) + op[r][3], 6)
                                    + 0.0 for r in range(3)]
                atom["chain_id"] = rename[atom["chain_id"]]
                out_atoms.append(atom)
            out_bonds += [(index[i], index[j]) for i, j in bonds if i in index and j in index]
            out_ranges += [{**r, "chain": rename[r["chain"]]} for r in ranges
                           if r["chain"] in rename]
            copies.append({"operator": serial, "chains": [rename[c] for c in chains],
                           "matrix": op})
    return out_atoms, sorted(out_bonds), out_ranges, copies


def _pdb_unit_cell(atoms, bonds, ranges, lattice, group):
    """(atoms, bonds, ranges, summary) of the unit cell: every operator of
    ``group`` applied to ``atoms``, each copy shifted by whole lattice vectors
    so its fractional centre lies in [0, 1). The identity copy keeps its
    chains; copy n is named <chain>-<n>. A copy of an atom within 0.1 A of an
    earlier copy of the same atom, modulo the lattice (an atom on a symmetry
    element), is merged into it; a bond to it is kept only when the merged
    atoms coincide in place."""
    from ._spacegroups import PDB_HALL, operators
    if group not in PDB_HALL:
        raise ValueError(f"PDB symmetry expansion does not support space group {group!r}; "
                         "it covers the chiral space groups under their wwPDB names")
    ops = [(r, [float(v) for v in t]) for r, t in operators(PDB_HALL[group])]
    scale = _pdb_scale(lattice)

    def cartesian(f):
        return [sum(f[i] * lattice[i][k] for i in range(3)) for k in range(3)]

    def apart(f, g):
        d = [a - b for a, b in zip(f, g)]
        return math.dist(cartesian([v - round(v) for v in d]), [0, 0, 0])

    fractional = [[sum(scale[i][k] * a["position"][k] for k in range(3)) for i in range(3)]
                  for a in atoms]
    images = [[[sum(r[i][k] * f[k] for k in range(3)) + t[i] for i in range(3)]
               for f in fractional] for r, t in ops]
    shifts = [[-math.floor(sum(f[i] for f in copy) / len(copy)) for i in range(3)]
              for copy in images]
    special = {n for n, f in enumerate(fractional)
               if any(apart(copy[n], f) < 0.1 for copy in images[1:])}
    out_atoms, out_bonds, out_ranges, placed, merged, dropped = [], set(), [], {}, 0, 0
    for k, copy in enumerate(images):
        rename = (lambda c: c) if k == 0 else (lambda c: f"{c}-{k + 1}")
        for n, f in enumerate(copy):
            earlier = next((j for j in range(k) if n in special
                            and apart(images[j][n], f) < 0.1), None)
            position = [round(v, 6) + 0.0 for v in
                        cartesian([a + b for a, b in zip(f, shifts[k])])]
            if earlier is not None:
                merged += 1
                kept = placed[(earlier, n)]
                placed[(k, n)] = kept if math.dist(out_atoms[kept]["position"],
                                                   position) < 0.1 else None
                continue
            atom = dict(atoms[n], position=position)
            if "chain_id" in atom or k:
                atom["chain_id"] = rename(atom.get("chain_id", ""))
            placed[(k, n)] = len(out_atoms)
            out_atoms.append(atom)
        for i, j in bonds:
            a, b = placed[(k, i)], placed[(k, j)]
            if a is None or b is None or a == b:
                dropped += 1
            else:
                out_bonds.add(tuple(sorted((a, b))))
        out_ranges += [{**r, "chain": rename(r["chain"])} for r in ranges]
    summary = {"space_group": group, "hall_symbol": PDB_HALL[group], "operators": len(ops),
               "atoms": len(out_atoms), "special_position_atoms_merged": merged,
               "bonds_dropped": dropped}
    return out_atoms, sorted(out_bonds), out_ranges, summary


def _pdb_structure(writer, sid, atoms, label=None, ranges=(), bonds=(), cell=None):
    section = {"id": sid, "kind": "structure"}
    if label:
        section["label"] = label
    grouped = _pdb_residues(atoms)
    if grouped:
        if grouped[0]:
            section["chains"] = grouped[0]
        section["residues"] = grouped[1]
    if ranges:
        section["secondary_structure"] = list(ranges)
    writer.add_section(section)
    writer.add_json(sid, "structure", {"atoms": atoms, "pbc": [cell is not None] * 3,
                                       "lattice_vectors": cell})
    if bonds:
        writer.add_json(sid, "bonds", {"pairs": [{"i": i, "j": j} for i, j in bonds]})
    return grouped is not None


def import_pdb(source, destination, *, expand_symmetry=False, expand_ncs=False):
    """Import PDB atom records with their atom-site identity.

    Every ATOM/HETATM record becomes an atom in source order, alternate
    locations included, carrying its name, residue, chain, insertion code,
    alternate location, occupancy, temperature factor and record type; the
    section lists chains and residues, HELIX/SHEET records become its
    ``secondary_structure`` ranges, and CONECT records its ``bonds`` member
    (distinct pairs, no inferred order) together with the SSBOND and LINK
    records resolved by atom-site identity; a bond to a symmetry mate (an
    operator other than 1555) is skipped and counted. Nothing is selected by occupancy or
    alternate location. A multi-model ensemble whose models share their
    ordered atoms and metadata becomes model 1 plus a trajectory of every
    model; otherwise each model is its own structure section, with CONECT
    serials resolved within each model. A CRYST1 cell (not the wwPDB
    placeholder) makes the structure periodic in the PDB standard
    orientation, checked against any SCALE records, and a structure.symmetry
    section records its space group; symmetry mates are not generated.
    Each REMARK 350 biomolecule becomes a nonperiodic structure section
    ``assembly_<n>`` built from the first model (identity copies keep their
    chain, other copies are named ``<chain>-<BIOMT serial>``), unless MTRIX
    records say NCS copies are missing from the file. With ``expand_ncs``,
    those copies are generated into an ``asymmetric_unit`` section, from
    which assemblies and the unit cell are then built. With
    ``expand_symmetry``, a periodic ``unit_cell`` section (with a P 1
    structure.symmetry section) also holds every space-group copy of the
    first model, for the chiral space groups under their wwPDB names.
    """
    models = []; serials = set(); omitted = {}; current = None; ranges = []; conects = []
    links = []; mtrix = []; remark350 = []
    cryst1 = None; scales = {}
    in_model = False; ended = False; charged = 0; declared = None
    text = _text(source)
    if any((ord(c) < 32 and c not in '\r\n') or ord(c) > 126 for c in text):
        raise ValueError("PDB requires ASCII records without control characters")
    for line in text.splitlines():
        if not line.strip():
            continue
        if any(ord(c) < 32 or ord(c) > 126 for c in line) or line[80:].strip():
            raise ValueError("PDB requires ASCII fixed-column records of at most 80 columns")
        if ended:
            raise ValueError("PDB has records after END")
        record = line[:6].strip()
        if line[:6].ljust(6) != record.ljust(6):
            raise ValueError("PDB record names must be left-aligned in columns 1-6")
        if record == "MODEL":
            if in_model:
                raise ValueError("PDB MODEL is missing ENDMDL")
            if current is not None and current[0] is None:
                raise ValueError("PDB coordinate records must all lie inside MODEL/ENDMDL")
            if line[6:10].strip() or line[14:].strip():
                raise ValueError("PDB MODEL serial must fit columns 11-14")
            serial = _pdb_integer(line[10:14], "model serial")
            if any(m[0] == serial for m in models):
                raise ValueError("PDB MODEL serials must be unique")
            current = (serial, [], []); models.append(current); serials = set()
            in_model = True
        elif record == "ENDMDL":
            if not in_model or not current[1]:
                raise ValueError("PDB ENDMDL requires a nonempty MODEL")
            in_model = False
        elif record == "NUMMDL":
            if line[6:10].strip() or line[14:].strip():
                raise ValueError("PDB NUMMDL count must fit columns 11-14")
            if declared is not None:
                raise ValueError("PDB allows one NUMMDL record")
            declared = _pdb_integer(line[10:14], "NUMMDL")
        elif record in {"ATOM", "HETATM"}:
            if current is not None and current[0] is not None and not in_model:
                raise ValueError("PDB coordinate record outside MODEL/ENDMDL")
            if current is None:
                current = (None, [], []); models.append(current)
            if len(line) < 78:
                raise ValueError("PDB atom requires coordinates and explicit element columns 77-78")
            if line[11] != ' ' or line[20] != ' ' or line[27:30].strip():
                raise ValueError("PDB atom fields overflow fixed-column separators")
            serial = _pdb_integer(line[6:11], "atom serial")
            if serial in serials:
                raise ValueError("PDB atom serials must be unique within a model")
            serials.add(serial)
            if len(serials) > 99999:
                raise ValueError("PDB model exceeds 99999 atoms")
            current[1].append(_pdb_atom(line, record))
            current[2].append(serial)
            charged += bool(line[78:80].strip())
        elif record == "CONECT":
            conects.append(_pdb_conect(line))
        elif record in {"MTRIX1", "MTRIX2", "MTRIX3"}:
            mtrix.append(line)
        elif record == "REMARK" and line[6:10].strip() == "350":
            remark350.append(line)
        elif record in {"SSBOND", "LINK"}:
            links.append(_pdb_link(line, record))
        elif record == "CRYST1":
            if cryst1 is not None:
                raise ValueError("PDB allows one CRYST1 record")
            cryst1 = _pdb_cryst1(line)
        elif record in {"SCALE1", "SCALE2", "SCALE3"}:
            n, row = _pdb_scale_record(line)
            if record[-1] != str(n) or n in scales:
                raise ValueError("PDB SCALE records must be SCALE1-3, each once")
            scales[n] = row
        elif record in {"HELIX", "SHEET"}:
            ranges.append(_pdb_range(line, record))
        elif record == "END":
            ended = True
        elif record in _PDB_METADATA:
            omitted[record] = omitted.get(record, 0) + 1
        else:
            raise ValueError(f"unsupported PDB record {record!r}")
    if in_model:
        raise ValueError("PDB MODEL is missing ENDMDL")
    if not models:
        raise ValueError("PDB has no atoms")
    if declared is not None and declared != len(models):
        raise ValueError("PDB NUMMDL disagrees with the number of models")

    def identity(atoms):
        return [{k: v for k, v in a.items() if k != "position"} for a in atoms]

    def where(m):
        return "the file" if m[0] is None else f"MODEL {m[0]}"

    cell, group, z, placeholder = _pdb_cell(cryst1, scales)
    ncs = _pdb_mtrix(mtrix)
    if expand_symmetry and cell is None:
        raise ValueError("PDB symmetry expansion requires a CRYST1 cell (not the placeholder)")
    if expand_symmetry and ncs[2] and not expand_ncs:
        raise ValueError("PDB symmetry expansion requires every MTRIX copy to be given "
                         "(or expand_ncs to generate them)")
    completed = None
    expansion = None
    biomolecules = _pdb_assemblies(remark350)
    built = []
    symmetry = []

    def add_symmetry(writer, sid):
        if cell is not None:
            symmetry.append(f"symmetry_{sid}" if sid != "structure" else "symmetry")
            writer.add_section({"id": symmetry[-1], "kind": "structure.symmetry",
                                "structure_ref": sid})
            writer.add_json(symmetry[-1], "data", {"space_group_symbol": group})

    first = models[0][1]
    shared = len(models) > 1 and all(identity(m[1]) == identity(first) for m in models[1:])
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "PDB structure import") as writer:
            if len(models) == 1 or shared:
                sections = ["structure"]
                bonds, repeated, linked = _pdb_model_bonds(conects, links, models[0],
                                                           where(models[0]))
                residues = _pdb_structure(writer, "structure", first, ranges=ranges,
                                          bonds=bonds, cell=cell)
                add_symmetry(writer, "structure")
                if shared:
                    sections.append("trajectory")
                    writer.add_section({"id": "trajectory", "kind": "trajectory",
                                        "structure_ref": "structure"})
                    writer.add_json("trajectory", "metadata", {
                        "symbols": [a["symbol"] for a in first], "pbc": [cell is not None] * 3,
                        "model_serials": [m[0] for m in models]})
                    _binary(writer, "trajectory", "coords",
                            (v for m in models for a in m[1] for v in a["position"]),
                            [len(models), len(first), 3])
                    if cell is not None:
                        _binary(writer, "trajectory", "lattice",
                                (v for _ in models for row in cell for v in row),
                                [len(models), 3, 3])
            else:
                sections = [f"model_{m[0]}" for m in models]
                resolved = [_pdb_model_bonds(conects, links, m, where(m)) for m in models]
                bonds, repeated, linked = resolved[0]
                residues = all([_pdb_structure(writer, sid, m[1], f"MODEL {m[0]}", ranges, b,
                                               cell)
                                for sid, m, (b, _, _) in zip(sections, models, resolved)])
                for sid in sections:
                    add_symmetry(writer, sid)
            base, source_section = (first, bonds, ranges), sections[0]
            if expand_ncs and ncs[2]:
                n_atoms, n_bonds, n_ranges, n_copies = _pdb_ncs(first, bonds, ranges, ncs[3])
                _pdb_structure(writer, "asymmetric_unit", n_atoms,
                               "Asymmetric unit (MTRIX copies generated)", n_ranges, n_bonds,
                               cell)
                add_symmetry(writer, "asymmetric_unit")
                sections.append("asymmetric_unit")
                base, source_section = (n_atoms, n_bonds, n_ranges), "asymmetric_unit"
                completed = {"section": "asymmetric_unit", "atoms": len(n_atoms),
                             "copies": n_copies}
            if not ncs[2] or completed:
                for ident, blocks in biomolecules:
                    sid = f"assembly_{ident}"
                    a_atoms, a_bonds, a_ranges, copies = _pdb_assembly(
                        blocks, *base, where(models[0]))
                    _pdb_structure(writer, sid, a_atoms, f"Biological assembly {ident}",
                                   a_ranges, a_bonds)
                    built.append({"section": sid, "biomolecule": int(ident),
                                  "source_section": source_section, "atoms": len(a_atoms),
                                  "copies": copies})
                sections += [b["section"] for b in built]
            if expand_symmetry:
                u_atoms, u_bonds, u_ranges, expansion = _pdb_unit_cell(*base, cell, group)
                _pdb_structure(writer, "unit_cell", u_atoms, f"Unit cell ({group} expanded)",
                               u_ranges, u_bonds, cell)
                writer.add_section({"id": "symmetry_unit_cell", "kind": "structure.symmetry",
                                    "structure_ref": "unit_cell"})
                writer.add_json("symmetry_unit_cell", "data", {"space_group_symbol": "P 1"})
                expansion.update(section="unit_cell", source_section=source_section)
                sections += ["unit_cell"]
        findings = _validated(temporary)
    atoms = [a for m in models for a in m[1]]
    return {"format": "pdb", "output": str(destination), "atoms": len(first),
            "model": models[0][0] if len(models) == 1 else None,
            "models": [m[0] for m in models] if len(models) > 1 else [],
            "ensemble": "trajectory" if shared else "structures" if len(models) > 1 else None,
            "sections": sections, "findings": findings,
            "alternate_location_atoms": sum("alt_loc" in a for a in atoms),
            "residues_listed": residues,
            "bonds": len(bonds), "conect_repeated_listings": repeated,
            "conect_hydrogen_bond_fields_omitted": sum(c[2] for c in conects),
            "ssbond": {"records": sum(r == "SSBOND" for r, _, _ in links), **linked["ssbond"]},
            "link": {"records": sum(r == "LINK" for r, _, _ in links), **linked["link"]},
            "symmetry_mate_bonds_skipped": {
                "records": sum(ops != ("1555", "1555") for _, _, ops in links),
                "operators": sorted({" ".join(ops) for _, _, ops in links
                                     if ops != ("1555", "1555")})},
            "cell": None if cell is None else {
                "space_group": group, "z": z, "symmetry_sections": symmetry,
                "expanded": expansion is not None,
                "note": "Coordinates are as deposited (usually the asymmetric unit); "
                        + ("the unit_cell section holds every space-group copy."
                           if expansion else "symmetry mates are not generated.")},
            "placeholder_cell": placeholder,
            "symmetry_expansion": expansion,
            "assemblies": built,
            "assemblies_not_built": None if not biomolecules or not ncs[2] or completed else
                "MTRIX copies are not given, so the deposited structure is incomplete",
            "mtrix": {"operators": ncs[0], "given": ncs[1], "not_given": ncs[2]},
            "ncs_expansion": completed,
            "secondary_structure_ranges": len(ranges),
            "secondary_structure_unresolved": _pdb_unresolved(ranges, first),
            "omitted_records": dict(sorted(omitted.items())),
            "atom_charge_fields_omitted": charged,
            "omitted": ["atom serials, segment identifiers and atom charges",
                        "helix identifiers, comments and lengths; strand sense and registration",
                        "SSBOND/LINK serials and distances, and which record listed a bond"],
            "scope": "Elements, coordinates, atom-site identity, HELIX/SHEET ranges, CONECT, SSBOND and LINK bonds and the CRYST1 cell only; REMARK 350 assemblies are expanded as separate structures; no bond order, bond to a symmetry mate, crystal symmetry expansion (unless requested) or NCS copy is reconstructed. No charge, multiplicity or job is inferred, and no alternate location or model is selected. Keep the original PDB."}


_PDB_IDENTITY = ("atom_name", "residue_name", "residue_seq")


def _pdb_fixed(value, width, label):
    text = f"{_number(value):{width}.2f}"
    if len(text) != width:
        raise ValueError(f"PDB {label} overflows its fixed-width field")
    return text


def _pdb_atom_name(name, symbol):
    # The standard alignment: a name shorter than four characters starts in
    # column 14 when its element has one letter.
    return name.ljust(4) if len(name) == 4 or len(symbol) == 2 else (" " + name).ljust(4)


def _pdb_identity_line(serial, atom, symbol, fields):
    """One ATOM/HETATM record from an atom's atom-site identity."""
    name, residue = str(atom["atom_name"]), str(atom["residue_name"])
    chain = str(atom.get("chain_id", "")) or " "
    code = str(atom.get("insertion_code", "")) or " "
    alt = str(atom.get("alt_loc", "")) or " "
    seq = atom["residue_seq"]
    if not (0 < len(name) <= 4 and len(residue) <= 3 and len(chain) == 1
            and len(code) == 1 and len(alt) == 1):
        raise ValueError("PDB atom-site identity overflows its fixed columns")
    if not isinstance(seq, int) or isinstance(seq, bool) or not -999 <= seq <= 9999:
        raise ValueError("PDB residue sequence numbers must be integers in -999..9999")
    name = _pdb_atom_name(name, symbol)
    occupancy = ("      " if atom.get("occupancy") is None
                 else _pdb_fixed(atom["occupancy"], 6, "occupancy"))
    b_factor = ("      " if atom.get("b_factor") is None
                else _pdb_fixed(atom["b_factor"], 6, "temperature factor"))
    record = "HETATM" if atom.get("hetero", False) else "ATOM  "
    return (f"{record}{serial:5d} {name}{alt}{residue:>3s} {chain}{seq:4d}{code}   "
            + ''.join(fields) + occupancy + b_factor + f"          {symbol:>2s}  ")


def _pdb_range_records(ranges, atoms):
    """HELIX then SHEET records for ``ranges``; returns (lines, coil count).
    Residue names come from the atoms; an end the structure does not contain
    gets a blank name. A sheet range without ``sheet_id`` is its own sheet."""
    names = {(a.get("chain_id", ""), a["residue_seq"], a.get("insertion_code")):
             str(a["residue_name"]) for a in atoms}
    helices, sheets, coils = [], {}, 0
    for i, r in enumerate(ranges):
        if r["type"] == "coil":
            coils += 1
        elif r["type"] == "helix":
            helices.append(r)
        else:
            sheets.setdefault(r.get("sheet_id", f"S{i + 1}"), []).append(r)

    def end(r, which):
        seq, code = r[f"{which}_seq"], str(r.get(f"{which}_insertion_code", ""))
        chain = str(r["chain"])
        if not (isinstance(seq, int) and -999 <= seq <= 9999 and len(chain) <= 1
                and len(code) <= 1):
            raise ValueError("PDB secondary-structure range overflows its fixed columns")
        name = names.get((chain, seq, code or None), "")
        return f"{name:>3s} {chain or ' '}", f"{seq:4d}{code or ' '}"

    lines = []
    for serial, r in enumerate(helices, 1):
        if serial > 999:
            raise ValueError("PDB export supports at most 999 helices")
        (res1, seq1), (res2, seq2) = end(r, "start"), end(r, "end")
        cls = f"{r['helix_class']:2d}" if "helix_class" in r else "  "
        lines.append(f"HELIX  {serial:3d} {serial:>3d} {res1} {seq1} {res2} {seq2}{cls}")
    for sheet, strands in sheets.items():
        if len(sheet) > 3 or len(strands) > 99:
            raise ValueError("PDB sheet identifiers have at most 3 characters and 99 strands")
        for k, r in enumerate(strands, 1):
            strand = r.get("strand", k)
            if not 1 <= strand <= 999:
                raise ValueError("PDB strand numbers must lie in 1-999")
            (res1, seq1), (res2, seq2) = end(r, "start"), end(r, "end")
            sense = " 0" if strand == 1 else "  "
            lines.append(f"SHEET  {strand:3d} {sheet:>3s}{len(strands):2d} {res1}{seq1} "
                         f"{res2}{seq2}{sense}")
    return lines, coils


def _pdb_assigned_bonds(reader, selected, n_atoms):
    """(pairs, section ids used, orders not written) of the bonds spec § 5.1
    assigns to ``selected``: its own ``bonds`` member, and every ``bonds``
    section whose ``structure_ref`` names it (or that has none, when it is
    the archive's first structure)."""
    sources = []
    if "bonds" in (selected.members or {}):
        sources.append((selected.id, selected.id))
    first = next(s.id for s in reader.sections if s.kind == "structure")
    raw = {s["id"]: s for s in reader.raw_manifest["sections"]}
    for s in reader.sections:
        if s.kind == "bonds" and raw[s.id].get("structure_ref", first) == selected.id:
            sources.append((s.id, s.id))
    pairs, orders = set(), 0
    for sid, _ in sources:
        for pair in reader.json(sid, "bonds")["pairs"]:
            i, j = pair["i"], pair["j"]
            if i == j or not (0 <= i < n_atoms and 0 <= j < n_atoms):
                raise ValueError("PDB CONECT export requires bonds between two atoms of "
                                 "the structure")
            pairs.add(tuple(sorted((i, j))))
            orders += pair.get("order", 1) != 1
    return sorted(pairs), {sid for sid, _ in sources}, orders


def _pdb_cell_records(reader, selected, lattice):
    """CRYST1 and SCALE1-3 for a lattice in the PDB standard orientation;
    returns (lines, space group, largest cell rounding error in angstrom)."""
    rows = [[_number(v) for v in row] for row in lattice]
    if any(abs(rows[i][k]) > 1e-4 for i, k in ((0, 1), (0, 2), (1, 2))) or \
            not all(rows[i][i] > 0 for i in range(3)):
        raise ValueError("PDB export requires the lattice in PDB orientation (a along x, "
                         "b in the xy-plane, right-handed); coordinates are not rotated")
    lengths = [math.sqrt(sum(v * v for v in row)) for row in rows]

    def angle(u, v):
        dot = sum(x * y for x, y in zip(rows[u], rows[v])) / (lengths[u] * lengths[v])
        return math.degrees(math.acos(max(-1.0, min(1.0, dot))))

    values = lengths + [angle(1, 2), angle(0, 2), angle(0, 1)]
    raw = {s["id"]: s for s in reader.raw_manifest["sections"]}
    first = next(s.id for s in reader.sections if s.kind == "structure")
    group = "P 1"
    for s in reader.sections:
        if s.kind == "structure.symmetry" and raw[s.id].get("structure_ref", first) == selected.id:
            data = reader.json(s.id, "data")
            group = data.get("space_group_symbol") or data.get("international") or group
    fields = [f"{v:9.3f}" for v in values[:3]] + [f"{v:7.2f}" for v in values[3:]]
    if any(len(f) != w for f, w in zip(fields, [9] * 3 + [7] * 3)) or len(group) > 11:
        raise ValueError("PDB CRYST1 cell or space group overflows its fixed columns")
    z = "   1" if group == "P 1" else "    "
    lines = ["CRYST1" + "".join(fields) + f" {group:<11s}{z}"]
    rounded = _pdb_cell_vectors(*(float(f) for f in fields))
    for n, row in enumerate(_pdb_scale(rounded), 1):
        lines.append(f"SCALE{n}    " + "".join(f"{v:10.6f}" for v in row) + f"     {0:10.5f}")
    error = max(abs(a - b) for r, q in zip(rows, rounded) for a, b in zip(r, q))
    return lines, group, error


def _pdb_link_records(pairs, atoms):
    """SSBOND then LINK records for the bonds between different residues:
    SSBOND for SG-SG between two cysteines (one record per residue pair),
    LINK for any other such bond with a hetero atom. Operators are 1555 and
    the length is the bond's distance. Returns (lines, SSBOND, LINK count)."""
    def residue(a):
        return a.get("chain_id", ""), a["residue_seq"], a.get("insertion_code")

    def place(fields):
        line = [" "] * 80
        for start, width, text in fields:
            if len(text) != width:
                raise ValueError("PDB SSBOND/LINK field overflows its fixed columns")
            line[start:start + width] = text
        return "".join(line).rstrip()

    def end(a):
        chain, seq, code = residue(a)
        if not -999 <= seq <= 9999:
            raise ValueError("PDB residue sequence numbers must be integers in -999..9999")
        return str(chain) or " ", f"{seq:4d}", str(code or " "), f"{str(a['residue_name']):>3s}"

    ssbonds, links = {}, []
    for i, j in pairs:
        a, b = atoms[i], atoms[j]
        if residue(a) == residue(b):
            continue
        length = f"{math.dist(*([_number(v) for v in x['position']] for x in (a, b))):5.2f}"
        if all(x["residue_name"] == "CYS" and x["atom_name"] == "SG" for x in (a, b)):
            ssbonds.setdefault(tuple(sorted((residue(a), residue(b)), key=str)), (a, b, length))
        elif a.get("hetero", False) or b.get("hetero", False):
            links.append((a, b, length))
    lines = []
    for serial, (a, b, length) in enumerate(ssbonds.values(), 1):
        if serial > 999:
            raise ValueError("PDB export supports at most 999 SSBOND records")
        (c1, s1, i1, _), (c2, s2, i2, _) = end(a), end(b)
        lines.append(place([(0, 6, "SSBOND"), (7, 3, f"{serial:3d}"), (11, 3, "CYS"), (15, 1, c1),
                            (17, 4, s1), (21, 1, i1), (25, 3, "CYS"), (29, 1, c2), (31, 4, s2),
                            (35, 1, i2), (59, 6, "  1555"), (66, 6, "  1555"),
                            (73, 5, length)]))
    for a, b, length in links:
        fields = [(0, 4, "LINK"), (59, 6, "  1555"), (66, 6, "  1555"), (73, 5, length)]
        for x, (name, alt, res, chain, seq, code) in ((a, (12, 16, 17, 21, 22, 26)),
                                                      (b, (42, 46, 47, 51, 52, 56))):
            c, s, i, r = end(x)
            atom_name = str(x["atom_name"])
            if not 0 < len(atom_name) <= 4:
                raise ValueError("PDB atom-site identity overflows its fixed columns")
            fields += [(name, 4, _pdb_atom_name(atom_name, _symbol(x["symbol"]).upper())),
                       (alt, 1, str(x.get("alt_loc", "")) or " "), (res, 3, r), (chain, 1, c),
                       (seq, 4, s), (code, 1, i)]
        lines.append(place(fields))
    return lines, len(ssbonds), len(links)


def _pdb_conect_records(pairs):
    """CONECT records listing every bond from both atoms, four partners per
    record, serials being 1-based atom positions."""
    partners: dict[int, list[int]] = {}
    for i, j in pairs:
        partners.setdefault(i, []).append(j)
        partners.setdefault(j, []).append(i)
    lines = []
    for atom in sorted(partners):
        bonded = sorted(partners[atom])
        for k in range(0, len(bonded), 4):
            lines.append(f"CONECT{atom + 1:5d}" + ''.join(f"{b + 1:5d}" for b in bonded[k:k + 4]))
    return lines


def export_pdb(source, destination, *, section=None):
    """Export one nonperiodic structure as PDB atom records.

    When every atom carries ``atom_name``, ``residue_name`` and
    ``residue_seq``, the records keep the atom-site identity (record type,
    names, residue, chain, insertion code, alternate location, occupancy and
    temperature factor) and the section's helix and sheet ranges become
    HELIX/SHEET records. Otherwise each atom gets its own synthetic UNL
    residue, a coordinate interchange subset with no biomolecular meaning.
    Either way the structure's bonds become CONECT records, without bond
    orders; with identity, a bond between two cysteine SG atoms also becomes
    an SSBOND record and any other bond between residues that involves a
    hetero atom a LINK record. A structure periodic on all three axes, with
    its lattice in the PDB standard orientation, gets CRYST1 and SCALE
    records (space group from an assigned structure.symmetry section, else
    P 1). Neither output is a deposition-ready PDB file.
    """
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, {"structure", "trajectory"})
        if selected.kind != "structure":
            raise ValueError("PDB export requires one structure, not a trajectory")
        geometry = reader.json(selected.id, "structure")
        pbc = geometry.get("pbc", [False] * 3)
        periodic = all(pbc) and geometry.get("lattice_vectors") is not None
        if not periodic and (geometry.get("lattice_vectors") is not None or any(pbc)):
            raise ValueError("PDB export needs periodicity on all three axes or no cell; "
                             "use extXYZ for partial periodicity")
        atoms = geometry["atoms"]
        identity = all(all(atom.get(k) is not None for k in _PDB_IDENTITY) for atom in atoms)
        limit = 99999 if identity else 9999
        if not 1 <= len(atoms) <= limit:
            raise ValueError("PDB export requires 1-99999 atoms with atom-site identity, or "
                             "1-9999 atoms without it (one synthetic residue per atom)")
        coils = ssbonds = link_count = 0
        cell_lines, group, cell_error = (_pdb_cell_records(reader, selected,
                                                           geometry["lattice_vectors"])
                                         if periodic else ([], None, None))
        bonds, bond_sections, orders = _pdb_assigned_bonds(reader, selected, len(atoms))
        connectivity = "CONECT, NO BOND ORDERS" if bonds else "NO CONNECTIVITY"
        if identity:
            cell_text = "CELL, NO SYMMETRY MATES" if periodic else "NO CELL"
            lines = [f"REMARK   1 QVF STRUCTURE EXPORT; {connectivity}; {cell_text}"]
            raw = next(s for s in reader.raw_manifest["sections"] if s["id"] == selected.id)
            records, coils = _pdb_range_records(raw.get("secondary_structure") or [], atoms)
            linked, ssbonds, link_count = _pdb_link_records(bonds, atoms)
            lines += records + linked + cell_lines
        else:
            lines = [f"REMARK   1 QVF GEOMETRY ONLY; SYNTHETIC UNL RESIDUES, {connectivity}",
                     "REMARK   1 OCCUPANCY 1.00 AND TEMPERATURE FACTOR 0.00 ARE PLACEHOLDERS"]
            lines += cell_lines
        max_error = 0.0
        for i, atom in enumerate(atoms, 1):
            symbol = _symbol(atom["symbol"]).upper()
            position = [_number(v) for v in atom["position"]]
            fields = [f"{v:8.3f}" for v in position]
            if len(position) != 3 or any(len(field) != 8 for field in fields):
                raise ValueError("PDB coordinates overflow fixed-width Real(8.3) fields")
            max_error = max(max_error, *(abs(v - float(field)) for v, field in zip(position, fields)))
            if identity:
                lines.append(_pdb_identity_line(i, atom, symbol, fields))
            else:
                name = (" " + symbol).ljust(4) if len(symbol) == 1 else symbol.ljust(4)
                lines.append(f"HETATM{i:5d} {name} UNL  {i:4d}    " + ''.join(fields)
                             + f"  1.00  0.00          {symbol:>2s}  ")
        lines += _pdb_conect_records(bonds)
        lines.append("END")
        if any(len(line) > 80 for line in lines):
            raise ValueError("PDB export produced a record wider than 80 columns")
        with _new_output(destination) as temporary:
            temporary.write_text('\n'.join(line.ljust(80) for line in lines) + '\n', encoding="ascii")
        report = _export_report(reader, selected, {selected.id} | bond_sections,
                                destination, "pdb")
        report.update(atoms=len(atoms), coordinate_precision_angstrom=0.001,
                      max_rounding_error_angstrom=max_error,
                      atom_site_identity=identity, coil_ranges_omitted=coils,
                      bonds=len(bonds), bond_orders_not_written=orders,
                      ssbond_records=ssbonds, link_records=link_count,
                      space_group=group, max_cell_rounding_error_angstrom=cell_error)
        if identity:
            report["generated"] = ("Sequential serials; ATOM unless the atom is marked hetero. "
                                   "Occupancy and temperature factor are blank where the archive "
                                   "records none; atom charges are not written.")
        else:
            report["generated"] = ("Synthetic serials, element-based atom names, one UNL residue per atom, occupancy 1.00 and temperature factor 0.00; no chemical or experimental meaning.")
        return report

def export_xyz(source, destination, *, section=None, extended=True):
    """Export one structure or trajectory; extXYZ preserves cell and pbc.

    extXYZ also writes the energy (eV), forces (eV/A) and stress (eV/A^3)
    of a structure (root energies when it is the archive's first structure,
    and its forces section) or of a trajectory (metadata energies and the
    gradient member).
    """
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, {"structure", "trajectory"})
        energies = gradients = stresses = None
        used = {selected.id}
        if selected.kind == "structure":
            frames = [reader.json(selected.id, "structure")]
            if extended:
                from .units import convert as convert_units
                structures = [x.id for x in reader.sections if x.kind == "structure"]
                total = (reader.raw_manifest.get("energies") or {}).get("total")
                if total is not None and structures[0] == selected.id:
                    energies = [convert_units(_number(total["value"]), total["units"], "ev")]
                forces = _assigned_forces(reader, selected.id)
                if forces is not None:
                    used.add(forces.id)
                    gradients = [reader.tensor(forces.id, "gradient").tolist()]
                    if forces.member("stress"):
                        stresses = [reader.tensor(forces.id, "stress").tolist()]
        else:
            metadata = reader.json(selected.id, "metadata")
            symbols = metadata.get("symbols")
            if symbols is None:
                numbers = metadata.get("atomic_numbers")
                if numbers is None:
                    atoms = metadata.get("atoms", [])
                    symbols = [a.get("symbol") for a in atoms]
                else:
                    symbols = [_ELEMENT_SYMBOLS[z] for z in numbers]
            coords = reader.tensor(selected.id, "coords").tolist()
            lattice = reader.tensor(selected.id, "lattice") if selected.member("lattice") else None
            cells = lattice.tolist() if lattice else None
            if extended and metadata.get("energies") is not None:
                energies = [_number(e) * HARTREE_IN_EV for e in metadata["energies"]]
            if extended and selected.member("gradient"):
                gradients = reader.tensor(selected.id, "gradient").tolist()
            frames = []
            for i, positions in enumerate(coords):
                if len(symbols) != len(positions):
                    raise ValueError("trajectory atom identities are incomplete")
                frames.append({"atoms": [_atom(s, p) for s, p in zip(symbols, positions)],
                               "pbc": metadata.get("pbc", [False] * 3),
                               "lattice_vectors": cells[i] if lattice and len(lattice.shape) == 3 else cells})
        lines = []
        for k, frame in enumerate(frames):
            cell = frame.get("lattice_vectors"); pbc = frame.get("pbc", [False] * 3)
            _check_cell(pbc, cell)
            if not extended and (cell is not None or any(pbc)):
                raise ValueError("plain XYZ cannot retain a cell or periodicity; use extXYZ")
            header = "QVF geometry export"
            if extended:
                properties = "species:S:1:pos:R:3" + (":forces:R:3" if gradients else "")
                header = f'Properties={properties} pbc="' + ' '.join('T' if x else 'F' for x in pbc) + '"'
                if cell is not None:
                    header += ' Lattice="' + ' '.join(format(_number(v), '.17g') for row in cell for v in row) + '"'
                if energies:
                    header += f" energy={energies[k]:.17g}"
                if stresses:
                    header += ' stress="' + ' '.join(format(_number(v) / _EV_PER_ANGSTROM3, '.17g')
                                                     for row in stresses[k] for v in row) + '"'
            lines.extend([str(len(frame["atoms"])), header])
            for i, atom in enumerate(frame["atoms"]):
                row = _symbol(atom["symbol"]) + ' ' + ' '.join(format(_number(v), '.17g') for v in atom["position"])
                if gradients:
                    row += ' ' + ' '.join(format(-_number(v) / _EV_PER_ANGSTROM, '.17g')
                                          for v in gradients[k][i])
                lines.append(row)
        with _new_output(destination) as temporary:
            temporary.write_text('\n'.join(lines) + '\n', encoding="utf-8")
        report = _export_report(reader, selected, used, destination, "extxyz" if extended else "xyz")
        report.update(frames=len(frames), energy=bool(energies), forces=bool(gradients),
                      stress=bool(stresses))
        return report


def import_cube(source, destination, *, kind="volume.generic", periodic=(False, False, False), includes_endpoint=True):
    """Import a single scalar cube in bohr, with explicit quantity/periodicity.

    A periodic cell is derived from voxel vectors times n-1 (endpoint) or n.
    Scalar values are unchanged; choosing a specific kind asserts its units.
    """
    if kind not in _SCALAR_KINDS:
        raise ValueError("cube conversion requires a scalar volume kind")
    if len(periodic) != 3 or not all(type(x) is bool for x in periodic) or type(includes_endpoint) is not bool:
        raise ValueError("periodic requires three booleans and includes_endpoint a boolean")
    lines = _text(source).splitlines()
    if len(lines) < 6:
        raise ValueError("truncated cube header")
    axes = re.findall(r'(?:OUTER|MIDDLE|INNER) LOOP:\s*([XYZ])', lines[1].upper())
    if axes and axes != ['X', 'Y', 'Z']:
        raise ValueError("cube axis order must be X/Y/Z")
    head = lines[2].split()
    if len(head) not in (4, 5):
        raise ValueError("invalid cube atom/origin line")
    natoms = int(head[0]); origin = [_number(v) for v in head[1:4]]
    if natoms < 1 or (len(head) == 5 and int(head[4]) != 1):
        raise ValueError("only positive-atom, single-scalar cubes are supported")
    shape = []; vectors = []
    for line in lines[3:6]:
        fields = line.split()
        if len(fields) != 4 or int(fields[0]) < 1:
            raise ValueError("cube axes require positive counts and bohr vectors")
        shape.append(int(fields[0])); vectors.append([_number(v) for v in fields[1:]])
    count = math.prod(shape)
    if count > _MAX_VALUES or natoms * 3 > _MAX_VALUES or len(lines) < 6 + natoms:
        raise ValueError("truncated or oversized cube")
    atoms = []; charges = []
    for line in lines[6:6 + natoms]:
        fields = line.split()
        if len(fields) != 5 or not 1 <= int(fields[0]) < len(_ELEMENT_SYMBOLS):
            raise ValueError("unsupported cube atom record")
        charges.append(_number(fields[1]))
        atoms.append(_atom(_ELEMENT_SYMBOLS[int(fields[0])], [_number(v) / BOHR_PER_ANGSTROM for v in fields[2:]]))
    values = [_number(v) for line in lines[6 + natoms:] for v in line.split()]
    if len(values) != count:
        raise ValueError("cube data count differs from grid shape")
    cell = None
    if any(periodic):
        cell = [[v * (shape[i] - int(includes_endpoint)) / BOHR_PER_ANGSTROM for v in row]
                for i, row in enumerate(vectors)]
    grid = {"shape": shape, "origin": origin, "voxel_vectors": vectors,
            "periodic": list(periodic), "includes_endpoint": includes_endpoint, "structure_ref": "structure"}
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "Scalar cube import") as writer:
            _structure(writer, atoms, list(periodic), cell)
            writer.add_section({"id": "volume", "kind": kind})
            writer.add_json("volume", "grid", grid)
            _binary(writer, "volume", "data", values, shape)
        findings = _validated(temporary)
    return {"format": "cube", "output": str(destination), "kind": kind, "findings": findings,
            "omitted": ["two comment lines", "cube atom charge column"],
            "nonzero_atom_charges_omitted": sum(v != 0 for v in charges)}


def export_cube(source, destination, *, section=None):
    """Export a selected scalar grid and its unambiguous structure in bohr."""
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, _SCALAR_KINDS)
        grid = reader.json(selected.id, "grid"); tensor = reader.tensor(selected.id, "data")
        if len(tensor.shape) != 3:
            raise ValueError("cube needs a scalar three-dimensional grid")
        sid = selected.to_json().get("structure_ref") or grid.get("structure_ref")
        structure = _select(reader, sid, {"structure"})
        atoms = reader.json(structure.id, "structure")["atoms"]
        if not atoms:
            raise ValueError("cube requires at least one atom")
        rows = ["QVF scalar grid export (bohr)", "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"]
        rows.append(str(len(atoms)) + ' ' + ' '.join(format(_number(v), '.17g') for v in grid["origin"]))
        for n, vector in zip(tensor.shape, grid["voxel_vectors"]):
            rows.append(str(n) + ' ' + ' '.join(format(_number(v), '.17g') for v in vector))
        for atom in atoms:
            z = _ELEMENT_SYMBOLS.index(_symbol(atom["symbol"]))
            rows.append(f"{z} 0 " + ' '.join(format(_number(v) * BOHR_PER_ANGSTROM, '.17g') for v in atom["position"]))
        with _new_output(destination) as temporary, temporary.open('w', encoding='utf-8') as stream:
            stream.write('\n'.join(rows) + '\n')
            for i, (value,) in enumerate(struct.iter_unpack('<f' if tensor.dtype == 'float32' else '<d', tensor.data)):
                stream.write(format(_number(value), '.17g') + ('\n' if i % 6 == 5 else ' '))
            stream.write('\n')
        report = _export_report(reader, selected, {selected.id, structure.id}, destination, "cube")
        report.update(kind=selected.kind, periodic=grid.get("periodic", [False] * 3),
                      includes_endpoint=grid.get("includes_endpoint", True),
                      note="Cube does not encode QVF kind or periodicity; retain these options when reimporting. Atom charge column is zero, not a population analysis.")
        return report


class _XsfLines:
    """Keep comments at section boundaries; never discard unknown records."""
    def __init__(self, text):
        self.lines = [line.strip() for line in text.splitlines() if line.strip()]
        self.index = 0
        self.comments = 0

    def boundary(self):
        while self.index < len(self.lines) and self.lines[self.index].startswith('#'):
            self.index += 1
            self.comments += 1
        return self.lines[self.index] if self.index < len(self.lines) else None

    def take(self):
        if self.index == len(self.lines):
            raise ValueError("truncated XSF section")
        line = self.lines[self.index]
        self.index += 1
        return line

    def expect(self, expected):
        if self.take() != expected:
            raise ValueError(f"XSF subset requires {expected}")

    def vector(self):
        fields = self.take().split()
        if len(fields) != 3:
            raise ValueError("XSF vector requires three components")
        return [_number(v) for v in fields]


def _xsf_atom(line):
    """(atom, forces in hartree/angstrom or None) of one XSF atom row."""
    fields = line.split()
    if len(fields) not in (4, 7):
        raise ValueError("XSF atom rows require identity, three coordinates and optionally "
                         "three force components")
    symbol = fields[0]
    if re.fullmatch(r'[0-9]+', symbol):
        z = int(symbol)
        if not 1 <= z < len(_ELEMENT_SYMBOLS):
            raise ValueError("unsupported XSF atomic number")
        symbol = _ELEMENT_SYMBOLS[z]
    return _atom(symbol, fields[1:4]), ([_number(v) for v in fields[4:]] if len(fields) == 7
                                        else None)


def import_xsf(source, destination, *, kind=None, periodic=None):
    """Import one XSF geometry and optional scalar datagrid, without resampling.

    Geometry periodicity follows CRYSTAL/SLAB/POLYMER. Grid periodicity is
    independent and must be asserted explicitly; its default is a finite box.
    Grid coordinates convert angstrom to bohr; choosing kind asserts the
    scalar values already have that canonical kind's units.
    """
    if kind is not None and kind not in _SCALAR_KINDS:
        raise ValueError("XSF conversion requires a scalar volume kind")
    if periodic is not None and (len(periodic) != 3 or not all(type(v) is bool for v in periodic)):
        raise ValueError("periodic requires three booleans")
    lines = _XsfLines(_text(source))
    header = lines.boundary()
    dimensions = {"CRYSTAL": 3, "SLAB": 2, "POLYMER": 1}
    lattice = None; conventional = False
    if header in dimensions:
        lines.take()
        pbc = [i < dimensions[header] for i in range(3)]
        lines.boundary(); lines.expect("PRIMVEC")
        lattice = [lines.vector() for _ in range(3)]
        if lines.boundary() == "CONVVEC":
            lines.take(); conventional = True
            _check_cell([True] * 3, [lines.vector() for _ in range(3)])
        lines.boundary()
        if lines.take() not in ("PRIMCOORD", "PRIMCOORD 1"):
            raise ValueError("XSF subset requires one primitive PRIMCOORD section")
        count = lines.take().split()
        if len(count) != 2 or count[1] != '1' or not 0 < int(count[0]) * 3 <= _MAX_VALUES:
            raise ValueError("invalid or oversized XSF primitive atom count")
        rows = [_xsf_atom(lines.take()) for _ in range(int(count[0]))]
        if lines.boundary() == "CONVCOORD":
            lines.take(); conventional = True
            count = lines.take().split()
            if len(count) != 2 or not re.fullmatch(r"[0-9]+", count[0]):
                raise ValueError("XSF CONVCOORD requires an atom count and 1")
            for _ in range(int(count[0])):
                _xsf_atom(lines.take())
    else:
        if header == "MOLECULE":
            lines.take(); lines.boundary()
        lines.expect("ATOMS")
        pbc = [False] * 3
        rows = []
        while lines.index < len(lines.lines):
            line = lines.lines[lines.index]
            if line.startswith('#') or line.startswith('BEGIN_'):
                break
            rows.append(_xsf_atom(lines.take()))
            if len(rows) * 3 > _MAX_VALUES:
                raise ValueError("oversized XSF geometry")
        if not rows:
            raise ValueError("XSF requires at least one atom")
    atoms = [atom for atom, _ in rows]
    forces = [f for _, f in rows]
    if any(f is None for f in forces):
        if any(f is not None for f in forces):
            raise ValueError("XSF forces must be given for every atom or none")
        forces = None
    _check_cell(pbc, lattice)
    flags = list(periodic) if periodic is not None else [False] * 3
    if any(flag and not axis for flag, axis in zip(flags, pbc)):
        raise ValueError("periodic XSF grid axes require matching periodic structure axes")
    grids = []
    while lines.boundary() == "BEGIN_BLOCK_DATAGRID_3D":
        lines.take()
        block = lines.take()
        if len(block.split()) != 1:
            raise ValueError("XSF datagrid block identifier must be one word")
        count = len(grids)
        while lines.boundary() != "END_BLOCK_DATAGRID_3D":
            # Some established writers omit the underscore before the identifier.
            match = re.fullmatch(r'BEGIN_DATAGRID_3D_?(\S+)', lines.take())
            if not match:
                raise ValueError("XSF requires named scalar 3D datagrids")
            fields = lines.take().split()
            if len(fields) != 3:
                raise ValueError("XSF grid requires three counts")
            shape = [int(v) for v in fields]
            total = sum(math.prod(g[1]["shape"]) for g in grids) + math.prod(shape)
            if min(shape) < 2 or total > _MAX_VALUES:
                raise ValueError("XSF grid requires 2 or more points per axis within the value limit")
            origin = lines.vector()
            spans = [lines.vector() for _ in range(3)]
            _check_cell([True] * 3, spans)
            values = []
            while len(values) < math.prod(shape):
                values.extend(_number(v) for v in lines.take().split())
            if len(values) != math.prod(shape):
                raise ValueError("XSF data count differs from grid shape")
            lines.expect("END_DATAGRID_3D")
            grids.append((f"{block}/{match.group(1)}", {
                "shape": shape, "origin": [v * BOHR_PER_ANGSTROM for v in origin],
                "voxel_vectors": [[v * BOHR_PER_ANGSTROM / (shape[i] - 1) for v in row]
                                  for i, row in enumerate(spans)],
                "includes_endpoint": True, "periodic": flags, "structure_ref": "structure"},
                values))
        lines.take()
        if len(grids) == count:
            raise ValueError("XSF datagrid block has no grid")
    if lines.boundary() is not None:
        raise ValueError("unsupported XSF record or additional geometry")
    if not grids and (kind is not None or periodic is not None):
        raise ValueError("XSF kind/periodic options require a scalar datagrid")
    sections = []
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "XSF geometry and scalar import") as writer:
            _structure(writer, atoms, pbc, lattice)
            if forces is not None:
                writer.add_section({"id": "forces", "kind": "forces", "structure_ref": "structure"})
                _binary(writer, "forces", "gradient",
                        (-v / BOHR_PER_ANGSTROM for row in forces for v in row), [len(atoms), 3])
            for n, (name, grid, values) in enumerate(grids, 1):
                sid = "volume" if n == 1 else f"volume_{n}"
                sections.append(sid)
                writer.add_section({"id": sid, "kind": kind or "volume.generic", "label": name})
                writer.add_json(sid, "grid", grid)
                nx, ny, nz = grid["shape"]
                # XSF index x + nx*(y + ny*z); QVF index z + nz*(y + ny*x).
                _binary(writer, sid, "data",
                        (values[x + nx*(y + ny*z)] for x in range(nx) for y in range(ny) for z in range(nz)),
                        grid["shape"])
        findings = _validated(temporary)
    return {"format": "xsf", "output": str(destination), "findings": findings,
            "kind": (kind or "volume.generic") if grids else None,
            "periodic": flags if grids else None,
            "grids": [{"section": sid, "label": g[0]} for sid, g in zip(sections, grids)],
            "forces": forces is not None, "conventional_cell_omitted": conventional,
            "omitted_comments": lines.comments, "omitted": []}


def export_xsf(source, destination, *, section=None):
    """Export one structure or scalar grid with its structure to static XSF.

    Periodic grids without duplicate endpoints are refused. No points are
    synthesized. Scalar kind and field periodicity are returned in the report
    so callers can assert them again when importing the text representation.
    """
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, _SCALAR_KINDS | {"structure"})
        grid = None; tensor = None
        if selected.kind != "structure":
            grid = reader.json(selected.id, "grid")
            tensor = reader.tensor(selected.id, "data")
            if len(tensor.shape) != 3 or min(tensor.shape) < 2 or math.prod(tensor.shape) > _MAX_VALUES:
                raise ValueError("XSF requires a bounded scalar 3D grid with at least two points per axis")
            if any(grid.get("periodic", [])) and not grid.get("includes_endpoint", True):
                raise ValueError("XSF periodic grids require included endpoints; no planes are synthesized")
            sid = selected.to_json().get("structure_ref") or grid.get("structure_ref")
            structure = _select(reader, sid, {"structure"})
        else:
            structure = selected
        payload = reader.json(structure.id, "structure")
        pbc = payload.get("pbc", [False] * 3); cell = payload.get("lattice_vectors")
        _check_cell(pbc, cell)
        headers = {(True, True, True): "CRYSTAL", (True, True, False): "SLAB", (True, False, False): "POLYMER"}
        fmt = lambda row: ' '.join(format(_number(v), '.17g') for v in row)
        atoms = payload["atoms"]
        if not atoms or len(atoms) * 3 > _MAX_VALUES:
            raise ValueError("XSF requires a bounded, nonempty geometry")
        if any(pbc):
            if tuple(pbc) not in headers:
                raise ValueError("XSF periodic axes must be the leading cell axes; axes are never permuted")
            rows = [headers[tuple(pbc)], "PRIMVEC", *(fmt(row) for row in cell), "PRIMCOORD", f"{len(atoms)} 1"]
        else:
            if cell is not None:
                raise ValueError("XSF molecular subset cannot retain a nonperiodic cell")
            rows = ["ATOMS"]
        forces = _assigned_forces(reader, structure.id) if grid is None else None
        gradient = reader.tensor(forces.id, "gradient").tolist() if forces is not None else None
        for k, a in enumerate(atoms):
            row = str(_ELEMENT_SYMBOLS.index(_symbol(a["symbol"]))) + ' ' + fmt(a["position"])
            if gradient is not None:
                # XSF forces are hartree/angstrom.
                row += ' ' + fmt(-_number(v) * BOHR_PER_ANGSTROM for v in gradient[k])
            rows.append(row)
        if grid is not None:
            spans = [[v * (tensor.shape[i] - 1) / BOHR_PER_ANGSTROM for v in row]
                     for i, row in enumerate(grid["voxel_vectors"])]
            _check_cell([True] * 3, spans)
            rows.extend(["BEGIN_BLOCK_DATAGRID_3D", "qvf_grid", "BEGIN_DATAGRID_3D_field",
                         ' '.join(str(n) for n in tensor.shape),
                         fmt(v / BOHR_PER_ANGSTROM for v in grid["origin"]), *(fmt(row) for row in spans)])
        with _new_output(destination) as temporary, temporary.open('w', encoding='utf-8') as stream:
            stream.write('\n'.join(rows) + '\n')
            if tensor is not None:
                nx, ny, nz = tensor.shape
                code = '<f' if tensor.dtype == 'float32' else '<d'
                size = struct.calcsize(code)
                for z in range(nz):
                    for y in range(ny):
                        stream.write(fmt(struct.unpack_from(code, tensor.data, size*(z + nz*(y + ny*x)))[0]
                                         for x in range(nx)) + '\n')
                stream.write("END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D\n")
        report = _export_report(reader, selected, {selected.id, structure.id}
                                | ({forces.id} if forces is not None else set()), destination, "xsf")
        report["forces"] = forces is not None
        if grid is not None:
            report.update(kind=selected.kind, periodic=grid.get("periodic", [False] * 3), includes_endpoint=True,
                          note="XSF retains structure periodicity but not QVF field kind or field periodicity; supply kind and periodic again on import.")
        return report
