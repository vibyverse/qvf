"""Export molecular GFN-xTB results without importing or running xTB.

The input contract is xTB's JSON writer and explicitly selected companion
files, documented at https://xtb-docs.readthedocs.io/en/latest/properties.html.
Field units and the synthetic Gaussian mode layout follow xTB 6.7.1
``src/main/json.F90`` and ``src/freq/io.f90``. No console-output scraping,
directory discovery, convergence inference or normal-mode calculation occurs.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import re
import tempfile

from .. import __version__
from ..convert import (_atom, _binary, _number, _structure, _text, _validated,
                       _xyz_frames, import_cube)
from .._molden import import_molden
from ..operations import _new_output
from ..reader import Reader
from ..units import BOHR_PER_ANGSTROM, HARTREE_IN_EV
from ..writer import Writer

__all__ = ["write_xtb"]
_ENERGY = re.compile(r"^\s*energy:\s*(\S+)\s+gnorm:\s*(\S+)\s+xtb:\s*(\S+)")
_METHODS = {"GFN0-xTB", "GFN1-xTB", "GFN2-xTB"}


def _integer(value, name, minimum=None):
    if type(value) is not int or (minimum is not None and value < minimum):
        raise ValueError(f"{name} must be an integer" + (f" >= {minimum}" if minimum is not None else ""))
    return value


def _vector(value, name, size=None):
    if not isinstance(value, list) or (size is not None and len(value) != size):
        raise ValueError(f"{name} has an invalid length")
    if any(type(x) not in (int, float) for x in value):
        raise ValueError(f"{name} requires numeric values")
    return [_number(x) for x in value]


def _json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate xTB JSON key: {key}")
            result[key] = value
        return result
    def constant(value):
        raise ValueError(f"non-finite xTB JSON constant: {value}")
    data = json.loads(_text(path), object_pairs_hook=pairs, parse_constant=constant)
    if not isinstance(data, dict) or data.get("method") not in _METHODS:
        raise ValueError("xTB JSON must identify GFN0-xTB, GFN1-xTB or GFN2-xTB")
    if not isinstance(data.get("xtb version"), str) or not data["xtb version"].strip():
        raise ValueError("xTB JSON must contain an xtb version")
    if "total energy" not in data or type(data["total energy"]) not in (int, float):
        raise ValueError("xTB JSON must contain a numeric total energy")
    return data


def _geometry(path):
    text = _text(path)
    if not text.lstrip().startswith("$"):
        frames = _xtb_xyz(text)
        if len(frames) != 1 or any(frames[0]["pbc"]) or frames[0]["lattice_vectors"] is not None:
            raise ValueError("geometry must be one molecular XYZ frame")
        return frames[0]["atoms"], {}
    atoms = []; identity = {}; group = None; factor = 1.0; ended = False
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ended:
            raise ValueError("unexpected data after coord $end")
        if line.startswith("$"):
            fields = line.lower().split()
            group = fields[0]
            if group == "$coord" and not atoms and fields in (["$coord"], ["$coord", "angs"]):
                factor = 1 / BOHR_PER_ANGSTROM if len(fields) == 1 else 1.0
            elif group in ("$chrg", "$spin") and len(fields) == 2:
                key = "charge" if group == "$chrg" else "unpaired"
                if key in identity:
                    raise ValueError("repeated coord identity group")
                identity[key] = int(fields[1])
            elif fields == ["$end"]:
                ended = True
            else:
                raise ValueError("coord supports only molecular $coord [angs], $chrg and $spin; use XYZ for other inputs")
        elif group == "$coord":
            fields = line.split()
            if len(fields) != 4:
                raise ValueError("coord atom rows require three coordinates and an element")
            atoms.append(_atom(fields[3].capitalize(), [_number(x) * factor for x in fields[:3]]))
        else:
            raise ValueError("unexpected coord data")
    if not atoms or not ended:
        raise ValueError("coord needs atoms and $end")
    return atoms, identity


def _xtb_xyz(text):
    """xTB preserves lowercase coord symbols in its otherwise standard XYZ."""
    lines = text.splitlines(); index = 0
    while index < len(lines) and lines[index].strip():
        count = int(lines[index])
        if count < 1 or index + count + 2 > len(lines):
            raise ValueError("invalid or truncated xTB XYZ frame")
        for i in range(index + 2, index + count + 2):
            fields = lines[i].split()
            if fields:
                fields[0] = fields[0].capitalize()
                lines[i] = " ".join(fields)
        index += count + 2
    return _xyz_frames("\n".join(lines))[0]


def _same_atoms(atoms, other, tolerance=2e-5):
    if len(atoms) != len(other) or any(
        a["atomic_number"] != b["atomic_number"] or
        any(abs(x - y) > tolerance for x, y in zip(a["position"], b["position"]))
        for a, b in zip(atoms, other)
    ):
        raise ValueError("companion geometry does not match the selected result geometry and atom order")


def _trajectory(path, atoms, energy):
    text = _text(path)
    frames = _xtb_xyz(text)
    energies = []; lines = text.splitlines(); index = 0
    for frame in frames:
        match = _ENERGY.match(lines[index + 1])
        if match is None:
            raise ValueError("trajectory requires xTB energy/gnorm/xtb frame comments")
        energies.append(_number(match[1])); _number(match[2])
        index += len(frame["atoms"]) + 2
        if any(frame["pbc"]) or frame["lattice_vectors"] is not None:
            raise ValueError("only molecular xTB trajectories are supported")
    _same_atoms(atoms, frames[-1]["atoms"])
    if not math.isclose(energies[-1], energy, abs_tol=2e-7, rel_tol=0):
        raise ValueError("last trajectory energy does not match xTB JSON")
    return frames, energies


def _bonds(path, count):
    pairs = []; seen = set()
    for line in _text(path).splitlines():
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError("wbo rows require two 1-based indices and a bond order")
        i, j = int(fields[0]) - 1, int(fields[1]) - 1
        key = tuple(sorted((i, j)))
        if not (0 <= i < count and 0 <= j < count) or i == j or key in seen:
            raise ValueError("wbo contains an invalid or repeated atom pair")
        seen.add(key)
        pairs.append({"i": i, "j": j, "order": _number(fields[2])})
    return pairs


def _modes(path, atoms, frequencies):
    """Read only xTB's g98fake2 layout, never a general Gaussian output.

    xTB excludes abs(frequency) <= 0.1 cm^-1 and rounds coordinates to
    6 decimals and mode vectors to 2. Its dummy force constants and Raman
    activities are deliberately ignored, not exported as calculated zeros.
    """
    lines = _text(path).splitlines()
    if not any("frequency output generated by the xtb code" in line for line in lines[:10]):
        raise ValueError("normal modes require xTB-generated g98.out")
    start = next((i for i, line in enumerate(lines) if "Standard orientation:" in line), None)
    if start is None:
        raise ValueError("g98.out is missing its reference geometry")
    reference = []
    for line in lines[start + 5:start + 5 + len(atoms)]:
        fields = line.split()
        if len(fields) != 6 or int(fields[0]) != len(reference) + 1:
            raise ValueError("invalid g98.out reference geometry")
        reference.append({"atomic_number": int(fields[1]), "position": [_number(x) for x in fields[3:]]})
    _same_atoms(atoms, reference)
    indices = [i for i, value in enumerate(frequencies) if abs(value) > 0.1]
    observed = []; displacements = []
    for i, line in enumerate(lines):
        if not line.strip().startswith("Frequencies --"):
            continue
        block = [_number(x) for x in line.split("--", 1)[1].split()]
        if not 1 <= len(block) <= 3 or i + 7 + len(atoms) > len(lines):
            raise ValueError("truncated g98.out mode block")
        if not lines[i + 6].split()[:2] == ["Atom", "AN"]:
            raise ValueError("invalid g98.out mode header")
        rows = [[] for _ in block]
        for j, row in enumerate(lines[i + 7:i + 7 + len(atoms)]):
            fields = row.split()
            if (len(fields) != 2 + 3 * len(block) or int(fields[0]) != j + 1
                    or int(fields[1]) != atoms[j]["atomic_number"]):
                raise ValueError("g98.out mode atom order or dimensions are invalid")
            values = [_number(x) for x in fields[2:]]
            for k in range(len(block)):
                rows[k].append(values[3*k:3*k + 3])
        observed.extend(block); displacements.extend(rows)
    if len(observed) != len(indices) or any(abs(a - frequencies[b]) > 0.00011 for a, b in zip(observed, indices)):
        raise ValueError("g98.out modes do not match xTB JSON frequencies")
    if not observed:
        raise ValueError("g98.out contains no nonzero modes")
    return indices, displacements


def _copy_companion(writer, archive, atoms, prefix, canonical_orbitals=False):
    with Reader(archive) as reader:
        _same_atoms(atoms, reader.json("structure", "structure")["atoms"])
        for section in reader.sections:
            if section.kind == "structure":
                continue
            spec = section.to_json(); spec.pop("members")
            spec["id"] = prefix; spec["structure_ref"] = "structure"
            writer.add_section(spec)
            for role, member in section.members.items():
                if member.format == "json":
                    payload = reader.json(section.id, role)
                    if role == "mo_metadata" and canonical_orbitals:
                        payload["orbital_kind"] = "canonical"
                    writer.add_json(prefix, role, payload)
                else:
                    with writer.begin_binary(prefix, role, member.dtype, member.shape) as sink:
                        reader.tensor_into(section.id, role, sink)


def write_xtb(results, geometry, target, *, trajectory=None, wbo=None,
              molden=None, spin=None, modes=None, cubes=(), charge=None,
              multiplicity=None, run_status=None, canonical_orbitals=False):
    """Write one new validated QVF from molecular GFN0/1/2-xTB outputs.

    ``results`` is xtbout.json, ``geometry`` is the corresponding final XYZ
    or molecular coord, and ``target`` must not exist. Companion paths are
    opt-in: an xTB XYZ trajectory, wbo, molden.input (explicit ``spin``),
    g98.out modes, and ``cubes`` as (path, QVF volume kind) pairs. Files must
    belong to the same calculation. Geometry is checked where available;
    JSON and wbo do not contain coordinates and cannot establish that link.

    Charge is supplied explicitly or read from coord; JSON's electron count
    is a valence count, not a nuclear-charge sum. Multiplicity comes from
    JSON's unpaired electrons or explicit input; contradictions are errors.
    Status is omitted unless the caller supplies
    ``converged`` or ``failed`` based on independent execution evidence.
    Missing optional data are reported, never manufactured. The report
    includes imported sections, omitted JSON fields and validation findings.
    No xTB, NumPy or viewer dependency is required.

    ``canonical_orbitals=True`` asserts that Molden contains the untransformed
    canonical xTB orbitals. The managed runner establishes this; generic
    file export leaves the orbital kind ``custom`` unless explicitly asserted.
    """
    data = _json(results); atoms, identity = _geometry(geometry)
    n = len(atoms); energy = _number(data["total energy"])
    if charge is not None:
        _integer(charge, "charge")
    if multiplicity is not None:
        _integer(multiplicity, "multiplicity", 1)
    if (molden is None) != (spin is None):
        raise ValueError("molden and its explicit spin interpretation must be supplied together")
    if type(canonical_orbitals) is not bool or (canonical_orbitals and molden is None):
        raise ValueError("canonical_orbitals requires a Molden file and a boolean assertion")
    if run_status not in (None, "converged", "failed"):
        raise ValueError("run_status must be converged, failed or omitted")
    used = {"method", "xtb version", "total energy"}
    electrons = data.get("number of electrons")
    if electrons is not None:
        electrons = _integer(electrons, "number of electrons", 0)
        used.add("number of electrons")
    charge = identity.get("charge") if charge is None else charge
    unpaired = data.get("number of unpaired electrons", identity.get("unpaired"))
    if unpaired is not None:
        unpaired = _integer(unpaired, "number of unpaired electrons", 0)
        if multiplicity is not None and multiplicity != unpaired + 1:
            raise ValueError("multiplicity disagrees with xTB unpaired electron count")
        multiplicity = unpaired + 1; used.add("number of unpaired electrons")
    if multiplicity is not None:
        _integer(multiplicity, "multiplicity", 1)
    if ("charge" in identity and charge != identity["charge"]) or (
        "unpaired" in identity and multiplicity != identity["unpaired"] + 1
    ):
        raise ValueError("coord charge/spin disagrees with xTB results")
    if electrons is not None and multiplicity is not None:
        if multiplicity - 1 > electrons or (electrons - multiplicity + 1) % 2:
            raise ValueError("inconsistent valence electron count and multiplicity")
    if charge is not None:
        total_electrons = sum(a["atomic_number"] for a in atoms) - charge
        if total_electrons < 0 or (multiplicity is not None and
                            (multiplicity - 1 > total_electrons or (total_electrons - multiplicity + 1) % 2)):
            raise ValueError("inconsistent electron count and multiplicity")
    if spin == "restricted" and multiplicity not in (None, 1):
        raise ValueError("restricted Molden spin conflicts with multiplicity")
    frames = _trajectory(trajectory, atoms, energy) if trajectory is not None else None
    pairs = _bonds(wbo, n) if wbo is not None else None
    charges = _vector(data["partial charges"], "partial charges", n) if "partial charges" in data else None
    if charges is not None:
        used.add("partial charges")
        if charge is not None and abs(sum(charges) - charge) > max(1e-5, n * 1e-8):
            raise ValueError("partial charges do not sum to the molecular charge")
    freqkey = "vibrational frequencies / rcm"
    frequencies = _vector(data[freqkey], freqkey) if freqkey in data else None
    vibration = {}
    if frequencies is not None:
        used.add(freqkey); vibration["frequencies"] = frequencies
        for key, name in (("IR intensities / km/mol", "ir_intensities"),
                          ("Raman activities / A^4/amu", "raman_activities"),
                          ("reduced masses", "reduced_masses")):
            if key in data:
                vibration[name] = _vector(data[key], key, len(frequencies)); used.add(key)
    if modes is not None and frequencies is None:
        raise ValueError("normal modes require JSON vibrational frequencies")
    normal_modes = _modes(modes, atoms, frequencies) if modes is not None else None
    extra = {}
    for key in ("orbital energies / eV", "fractional occupation"):
        if key in data:
            extra[key] = _vector(data[key], key); used.add(key)
    if extra:
        if set(extra) != {"orbital energies / eV", "fractional occupation"} or len(extra["orbital energies / eV"]) != len(extra["fractional occupation"]):
            raise ValueError("orbital energies and fractional occupation must have matching lengths")
        count = _integer(data.get("number of molecular orbitals"), "number of molecular orbitals", 1)
        if count != len(extra["fractional occupation"]):
            raise ValueError("orbital arrays disagree with number of molecular orbitals")
        if any(x < 0 or x > 2 for x in extra["fractional occupation"]):
            raise ValueError("fractional occupations must lie between zero and two")
        if electrons is not None and abs(sum(extra["fractional occupation"]) - electrons) > max(1e-5, count * 1e-8):
            raise ValueError("fractional occupations disagree with the valence electron count")
        used.add("number of molecular orbitals")
        extra = {"orbital_energies_hartree": [x / HARTREE_IN_EV for x in extra["orbital energies / eV"]],
                 "orbital_occupations": extra["fractional occupation"]}
    if vibration:
        extra["vibrational_properties"] = vibration
    if electrons is not None:
        extra["valence_electrons"] = electrons
    imported = ["structure"]; companion_reports = []; missing = []
    with _new_output(target) as temporary, tempfile.TemporaryDirectory(prefix="qvf-xtb-") as work:
        with Writer(temporary, "qvf.adapters.xtb", __version__, "xTB result export") as writer:
            _structure(writer, atoms, [False] * 3, None)
            provenance = {"method": data["method"], "software": [{"name": "xtb", "version": data["xtb version"]}],
                          "scf_energy": {"value": energy, "units": "hartree"}, "x_xtb": extra}
            for key, value in (("charge", charge), ("multiplicity", multiplicity),
                               ("run_status", run_status)):
                if value is not None:
                    provenance[key] = value
            writer.set_provenance(**provenance)
            writer.set_energies(energy, method=data["method"])
            if "dipole / a.u." in data:
                dipole = _vector(data["dipole / a.u."], "dipole / a.u.", 3)
                # NIST: https://physics.nist.gov/cgi-bin/cuu/Value?auedm
                vector = [x * 2.541746473 for x in dipole]
                writer.set_dipole_moment(vector_debye=vector, total_debye=math.sqrt(sum(x*x for x in vector)))
                used.add("dipole / a.u.")
            def section(sid, kind):
                writer.add_section({"id": sid, "kind": kind, "structure_ref": "structure"})
                imported.append(sid)
            if charges is not None:
                section("charges", "atom_properties")
                if data["method"] == "GFN0-xTB":
                    # GFN0 uses an electronegativity-equilibration charge model.
                    role = "custom.xtb_charge"
                    writer.add_json("charges", "meta", {role: {
                        "units": "e", "kind": "scalar",
                        "description": "GFN0-xTB partial charge; positive means electron depletion."}})
                else:
                    role = "mulliken_charge"
                _binary(writer, "charges", role, charges, [n])
            if pairs is not None:
                section("bond_orders", "bond_orders")
                writer.add_json("bond_orders", "bond_orders", {"method": "wiberg", "pairs": pairs})
            if frames is not None:
                section("trajectory", "trajectory")
                writer.add_json("trajectory", "metadata", {"symbols": [a["symbol"] for a in atoms],
                                "coordinate_unit": "angstrom", "energies": frames[1]})
                _binary(writer, "trajectory", "coords", (x for f in frames[0] for a in f["atoms"] for x in a["position"]), [len(frames[0]), n, 3])
            for name, kind in (("ir_intensities", "spectra.ir"), ("raman_activities", "spectra.raman")):
                if name in vibration:
                    section(name, kind)
                    writer.add_json(name, "spectrum", {"frequencies": frequencies, "intensities": vibration[name]})
            if normal_modes is not None:
                indices, vectors = normal_modes
                section("vibrations", "vibrations")
                writer.add_json("vibrations", "metadata", {**{k: [v[i] for i in indices] for k, v in vibration.items()},
                                "symbols": [a["symbol"] for a in atoms], "x_displacement_decimals": 2})
                _binary(writer, "vibrations", "displacements", (x for mode in vectors for row in mode for x in row), [len(vectors), n, 3])
            if molden is not None:
                archive = Path(work) / "orbitals.qvf"
                report = import_molden(molden, archive, spin=spin, dialect="xtb")
                _copy_companion(writer, archive, atoms, "wavefunction", canonical_orbitals)
                companion_reports.append(report); imported.append("wavefunction")
            for i, (path, kind) in enumerate(cubes):
                archive = Path(work) / f"cube-{i}.qvf"
                report = import_cube(path, archive, kind=kind, periodic=[False] * 3)
                _copy_companion(writer, archive, atoms, f"volume_{i}")
                companion_reports.append(report); imported.append(f"volume_{i}")
        findings = _validated(temporary)
    if molden is None:
        missing.append("Orbital coefficients/basis not supplied; JSON orbital levels are provenance only.")
    if frequencies is not None and modes is None:
        missing.append("Normal-mode displacements not supplied; frequencies/spectra do not imply animated modes.")
    if run_status is None:
        missing.append("Convergence status not supplied; output-file presence is not success evidence.")
    return {"output": str(target), "sections": imported, "findings": findings,
            "omitted_json_fields": sorted(set(data) - used), "limitations": missing,
            "companions": [{k: v for k, v in report.items() if k != "output"} for report in companion_reports]}
