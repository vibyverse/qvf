"""The unit vocabulary (spec § 4.3, "Unit names" and "Declared units").

The schema's ``$defs/UnitName`` is the single source of truth. These tests
hold the spec table, ``qvf_units``, the C++ ``qvf/units.hpp`` table and the
JavaScript ``units.ts`` table to it, check the conversions, and pin how the
writers and the validator treat declared units.
"""

import hashlib
import json
import os
import re
import sys
import zipfile

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))

import qvf_reader  # noqa: E402
import qvf_units  # noqa: E402
import qvf_writer  # noqa: E402
from qvf_reader import validate_qvf  # noqa: E402
from qvf_units import (LEGACY_UNITS, UNITS, canonical, convert,  # noqa: E402
                       fixed_units, is_unit)
from qvf_writer import QvfWriter  # noqa: E402

_ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

with open(qvf_reader._SCHEMA_PATH, encoding="utf-8") as _fh:
    DEFS = json.load(_fh)["$defs"]

WATER = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]


def _repo_text(*parts):
    path = os.path.join(_ROOT, *parts)
    if not os.path.exists(path):
        pytest.skip(f"{'/'.join(parts)} not present (standalone checkout)")
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _table_after(text, heading):
    """First-column backticked names of the first table after ``heading``."""
    at = text.index(heading)
    names, in_table = [], False
    for line in text[at:].splitlines()[1:]:
        if line.startswith("|"):
            in_table = True
            m = re.match(r"^\| `([^`]+)` \|", line)
            if m:
                names.append(m.group(1))
        elif in_table:
            break
    return names


# --- one vocabulary everywhere ------------------------------------------------

def test_module_reads_the_schema_vocabulary():
    assert UNITS == tuple(DEFS["UnitName"]["enum"])
    assert LEGACY_UNITS == tuple(DEFS["LegacyUnitName"]["enum"])
    assert len(UNITS) == len(set(UNITS)) == 45
    assert not set(UNITS) & set(LEGACY_UNITS)


# The 25 names of the first vocabulary (2026-10-08). Later names are appended,
# so every earlier name keeps its place and its meaning.
_FIRST_VOCABULARY = (
    "angstrom", "bohr", "hartree", "ev", "kcal_per_mol", "kj_per_mol", "cm-1",
    "nm", "debye", "mhz", "hz", "ppm", "km_per_mol", "a3", "bohr3", "gpa", "fs",
    "ps", "kelvin", "atm", "e", "amu", "radian", "degree", "dimensionless")


def test_vocabulary_only_grows_at_the_end():
    assert UNITS[:len(_FIRST_VOCABULARY)] == _FIRST_VOCABULARY


def test_compound_names_follow_the_spelling_rule():
    """Spec § 4.3, "Compound names": unit names joined by `_`, one `_per_`
    before the denominator, a power as a trailing digit, a lone reciprocal as
    `-1`. The rule describes the closed list. It is not a parser."""
    base = {"angstrom", "bohr", "hartree", "ev", "debye", "gpa", "fs", "e",
            "c", "m", "states"}
    factor = r"(?:%s)[2-9]?" % "|".join(sorted(base, key=len, reverse=True))
    product = rf"{factor}(?:_{factor})*"
    rule = re.compile(rf"^(?:{product}(?:_per_{product})?|{factor}-1)$")
    for name in UNITS[len(_FIRST_VOCABULARY):]:
        assert rule.match(name), name
    # Rule-shaped strings outside the list are not unit names.
    for name in ("e_per_angstrom3", "ev_per_angstrom2", "bohr-1", "hartree/bohr"):
        assert not is_unit(name), name


def test_spec_table_is_the_schema_enum():
    """Every unit name of the spec § 4.3 table is in $defs/UnitName and vice
    versa, in the same order; the same for the legacy spellings."""
    spec = _repo_text("spec", "qvf-format-spec.md")
    assert _table_after(spec, "#### Unit names") == list(UNITS)
    assert _table_after(spec, "**Legacy spellings.**") == list(LEGACY_UNITS)


def test_cpp_table_is_the_schema_enum():
    # Generated from the schema by scripts/gen_kinds.py (test_generated_tables
    # checks it is current); qvf/units.hpp re-exports it.
    header = _repo_text("cpp", "include", "qvf", "kinds_generated.hpp")
    body = header[header.index("kUnitNames"):]
    body = body[body.index("{"):body.index("};")]
    assert re.findall(r'"([^"]+)"', body) == list(UNITS)
    legacy = header[header.index("kLegacyUnitSpellings"):]
    legacy = legacy[:legacy.index("};")]
    assert re.findall(r'\{"([^"]+)", Unit::', legacy) == list(LEGACY_UNITS)
    enum = header[header.index("enum class Unit {"):]
    enum = enum[:enum.index("};")]
    assert re.findall(r"//\s*(\S+)\s*$", enum, re.M) == list(UNITS)


def test_js_table_is_the_schema_enum():
    source = _repo_text("js", "src", "units.ts")
    body = source[source.index("UNIT_NAMES = ["):]
    body = body[:body.index("]")]
    assert re.findall(r'"([^"]+)"', body) == list(UNITS)
    legacy = source[source.index("LEGACY_UNIT_NAMES"):]
    legacy = legacy[legacy.index("{"):legacy.index("}")]
    assert re.findall(r"(\w+): \"", legacy) == list(LEGACY_UNITS)


def test_value_with_units_accepts_names_and_legacy_spellings():
    assert DEFS["ValueWithUnits"]["properties"]["units"]["oneOf"] == [
        {"$ref": "#/$defs/UnitName"}, {"$ref": "#/$defs/LegacyUnitName"}]
    assert DEFS["LegacyUnitName"].get("deprecated") is True


def _units_specs(node, path=""):
    """Every `units` property definition of the payload definitions."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "properties" and isinstance(value, dict) and "units" in value:
                yield path, value["units"]
            yield from _units_specs(value, f"{path}/{key}")
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _units_specs(value, f"{path}/{i}")


def test_every_payload_units_object_is_typed():
    seen = set()
    for path, spec in _units_specs(DEFS):
        name = path.split("/")[1]
        if name in ("ValueWithUnits",):
            continue
        if name == "SpectrumPayload":          # broadening.units, a string
            assert spec["$ref"] == "#/$defs/UnitName", path
            continue
        if spec.get("$ref") == "#/$defs/UnitName":   # a single declared unit
            continue
        if "$ref" in spec:                             # a named units object
            spec = DEFS[spec["$ref"].rsplit("/", 1)[1]]
        seen.add(name)
        assert spec.get("type") == "object", path
        if name != "ReactionPathMetadataPayload":   # typed through allOf
            # Either an open map of unit names, or a closed object whose
            # every value is a fixed vocabulary name (forces / hessian meta).
            assert spec["additionalProperties"] in (
                {"$ref": "#/$defs/UnitName"}, False), path
            if spec["additionalProperties"] is False:
                assert spec.get("properties"), path
                assert all("const" in prop for prop in spec["properties"].values()), path
        for key, prop in spec.get("properties", {}).items():
            assert prop.get("const", "angstrom") in UNITS, (path, key)
    assert seen >= {"VolumeGridPayload", "TrajectoryMetadataPayload",
                    "ReactionPathMetadataPayload", "ScanSurfaceMetadataPayload",
                    "KpathPayload", "DosMetaPayload", "FermiMeshPayload",
                    "QpathPayload", "PhononMetaPayload", "QtaimPayload"}
    spectrum = DEFS["SpectrumPayload"]["properties"]
    assert spectrum["x_units"]["$ref"] == spectrum["y_units"]["$ref"] == \
        "#/$defs/UnitName"


def test_fixed_units_follow_the_schema():
    assert fixed_units("VolumeGridPayload") == {"origin": "bohr",
                                                "voxel_vectors": "bohr"}
    # The two trajectory-family payloads share TrajectoryMetadataCommon and
    # differ in the convention of their `lattice` member (spec § 4.3).
    common = {"energies": "hartree", "coords": "angstrom", "times": "fs",
              "timestep_fs": "fs", "velocities": "angstrom_per_fs",
              "gradient": "hartree_per_bohr"}
    assert fixed_units("TrajectoryMetadataCommon") == common
    assert fixed_units("TrajectoryMetadataPayload") == {
        **common, "lattice": "angstrom", "cell_volumes": "a3",
        "temperatures": "kelvin"}
    assert fixed_units("ReactionPathMetadataPayload") == {**common,
                                                          "lattice": "bohr"}
    assert fixed_units("KpathPayload") == {"fermi": "ev", "fermi_energy_ev": "ev",
                                           "eigenvalues": "ev"}
    assert fixed_units("QtaimPayload") == {}
    with pytest.raises(KeyError):
        fixed_units("NoSuchPayload")


# --- names and conversions -----------------------------------------------------

def test_is_unit_and_canonical():
    assert is_unit("ev") and is_unit("cm-1") and is_unit("dimensionless")
    assert not is_unit("eV") and not is_unit("Ev") and not is_unit(None)
    assert canonical("Eh") == canonical("Hartree") == canonical("hartree") == "hartree"
    assert canonical("eV") == "ev"
    for bad in ("rad", "Angstrom", "", 3):
        with pytest.raises(ValueError):
            canonical(bad)


def test_bohr_per_angstrom_is_the_writer_constant():
    assert qvf_units.BOHR_PER_ANGSTROM is qvf_writer.BOHR_PER_ANGSTROM
    assert convert(1.0, "angstrom", "bohr") == qvf_writer.BOHR_PER_ANGSTROM
    assert convert(convert(2.5, "angstrom", "bohr"), "bohr", "angstrom") == \
        pytest.approx(2.5, rel=1e-15)
    # 1 bohr = 0.529177210903 Å (CODATA 2018), to the ten decimals kept.
    assert convert(1.0, "bohr", "angstrom") == pytest.approx(0.529177210903,
                                                             rel=1e-10)


@pytest.mark.parametrize("src,dst,value,want,rel", [
    ("hartree", "ev", 1.0, 27.211386245988, 1e-15),
    ("hartree", "kj_per_mol", 1.0, 2625.4996394799, 1e-15),
    ("hartree", "kcal_per_mol", 1.0, 627.5094740631, 1e-12),
    ("kcal_per_mol", "kj_per_mol", 1.0, 4.184, 1e-14),
    ("hartree", "cm-1", 1.0, 219474.6313632, 1e-15),
    ("ev", "cm-1", 1.0, 8065.543937349, 1e-11),
    ("ev", "kj_per_mol", 1.0, 96.48533212331, 1e-11),
    ("ev", "nm", 1.0, 1239.841984332, 1e-15),
    ("nm", "ev", 300.0, 4.132806614440, 1e-11),
    ("nm", "cm-1", 500.0, 20000.0, 1e-11),
    ("cm-1", "nm", 20000.0, 500.0, 1e-11),
    ("Eh", "eV", 2.0, 54.422772491976, 1e-15),
])
def test_energy_conversions(src, dst, value, want, rel):
    assert convert(value, src, dst) == pytest.approx(want, rel=rel)


def test_conversions_round_trip_and_take_arrays():
    energies = ("hartree", "ev", "kcal_per_mol", "kj_per_mol", "cm-1", "nm")
    for a in energies:
        for b in energies:
            assert convert(convert(0.37, a, b), b, a) == pytest.approx(0.37,
                                                                       rel=1e-13)
    arr = np.array([1.0, 2.0])
    assert np.allclose(convert(arr, "hartree", "ev"), arr * 27.211386245988)
    assert np.allclose(convert(np.array([200.0, 400.0]), "nm", "ev"),
                       [6.19920992166, 3.09960496083])
    assert convert(5.0, "kelvin", "kelvin") == 5.0


@pytest.mark.parametrize("src,dst", [("angstrom", "ev"), ("kelvin", "atm"),
                                     ("bohr", "nm"), ("debye", "e"),
                                     ("ev", "furlong")])
def test_unsupported_conversions_raise(src, dst):
    with pytest.raises(ValueError):
        convert(1.0, src, dst)


# --- writers declare, never convert ------------------------------------------

GRID = ([-1.0, -1.0, -1.0], [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]])


def _doc(path, kind, role):
    with qvf_reader.QvfReader(path) as r:
        sec = next(s for s in r.sections if s["kind"] == kind)
        return r.read_json(sec["members"][role])


def test_provenance_conveniences_write_unit_names():
    w = QvfWriter("t", "1")
    w.set_provenance(scf_energy_eh=-74.9, fermi_energy_ev=-5.1)
    prov = w.build_manifest()["provenance"]
    assert prov["scf_energy"] == {"value": -74.9, "units": "hartree"}
    assert prov["fermi_energy"] == {"value": -5.1, "units": "ev"}


def test_writer_declares_units(tmp_path):
    w = QvfWriter("t", "1")
    w.add_structure(WATER)
    rho = np.ones((2, 2, 2), dtype=np.float32)
    w.add_density(rho, *GRID, units={"origin": "bohr", "voxel_vectors": "bohr"},
                  section_id="rho")
    w.add_trajectory({"symbols": ["O", "H", "H"]}, np.zeros((2, 3, 3)),
                     units={"energies": "hartree", "coords": "angstrom"})
    w.add_scan_surface({"axis_a_label": "r", "axis_b_label": "a"},
                       [1.0, 2.0], [0.1, 0.2], np.zeros((2, 2)),
                       units={"axis_a": "bohr", "axis_b": "degree"})
    w.add_spectrum("spectra.uvvis", frequencies=[200.0], intensities=[0.1],
                   x_units="nm", y_units="dimensionless")
    path = str(tmp_path / "u.qvf")
    w.write(path)
    report = validate_qvf(path)
    assert report["ok"] and not report["warnings"], report["findings"]
    grid = _doc(path, "volume.density", "grid")
    assert grid["units"] == {"origin": "bohr", "voxel_vectors": "bohr"}
    assert grid["origin"] == [-1.0, -1.0, -1.0]          # nothing converted
    assert _doc(path, "trajectory", "metadata")["units"]["coords"] == "angstrom"
    assert _doc(path, "scan.surface", "metadata")["units"]["axis_b"] == "degree"
    assert _doc(path, "spectra.uvvis", "spectrum")["frequencies"] == [200.0]


@pytest.mark.parametrize("call,match", [
    (lambda w: w.add_density(np.ones((2, 2, 2)), *GRID,
                             units={"origin": "Bohr"}), "not a QVF unit name"),
    (lambda w: w.add_density(np.ones((2, 2, 2)), *GRID,
                             units={"origin": "angstrom"}), "must be 'bohr'"),
    (lambda w: w.add_trajectory({"symbols": ["O"], "units": {"energies": "Eh"}},
                                np.zeros((1, 1, 3))), "not a QVF unit name"),
    (lambda w: w.add_trajectory({"symbols": ["O"]}, np.zeros((1, 1, 3)),
                                units={"energies": "ev"}), "must be 'hartree'"),
    (lambda w: w.add_trajectory({"symbols": ["O"], "units": {}},
                                np.zeros((1, 1, 3)), units={}), "both"),
    (lambda w: w.add_reaction_path({"symbols": ["O"],
                                    "units": {"lattice": "angstrom"}},
                                   np.zeros((1, 1, 3))), "must be 'bohr'"),
    (lambda w: w.add_scan_surface({}, [1.0], [1.0], np.zeros((1, 1)),
                                  units={"axis_b": "rad"}), "not a QVF unit name"),
    (lambda w: w.add_spectrum("spectra.ir", frequencies=[1.0], intensities=[1.0],
                              x_units="1/cm"), "not a QVF unit name"),
    (lambda w: w.add_spectrum("spectra.ir", payload={
        "frequencies": [1.0], "intensities": [1.0],
        "broadening": {"type": "gaussian", "units": "wavenumbers"}}),
     "not a QVF unit name"),
    (lambda w: w.add_topology_qtaim({"points": []},
                                    units={"position": "au"}), "not a QVF unit name"),
])
def test_writer_refuses_bad_declarations(call, match):
    w = QvfWriter("t", "1")
    with pytest.raises(ValueError, match=match):
        call(w)
    assert not w.build_manifest()["sections"]


# --- the validator: W_UNIT_UNKNOWN --------------------------------------------

def _codes(report):
    return sorted(f["code"] for f in report["findings"])


def _member(path, data):
    return {"path": path, "format": "json",
            "sha256": hashlib.sha256(data).hexdigest()}


def _archive(tmp_path, provenance=None, kind=None, role=None, payload=None):
    """A water structure plus, optionally, one JSON-member section, written
    byte by byte so that the writer's own checks do not interfere."""
    struct = json.dumps({"atoms": WATER, "pbc": [False] * 3,
                         "lattice_vectors": None}).encode()
    files = {"s/structure.json": struct}
    sections = [{"id": "structure", "kind": "structure",
                 "members": {"structure": _member("s/structure.json", struct)}}]
    if kind is not None:
        data = json.dumps(payload).encode()
        files["p/payload.json"] = data
        sections.append({"id": "p", "kind": kind,
                         "members": {role: _member("p/payload.json", data)}})
    manifest = {"qvf_version": 1,
                "source": {"program": "t", "version": "1", "calculation": ""},
                "sections": sections}
    if provenance is not None:
        manifest["provenance"] = provenance
    path = str(tmp_path / "raw.qvf")
    with zipfile.ZipFile(path, "w") as zf:
        for name, data in files.items():
            zf.writestr(name, data)
        zf.writestr("manifest.json", json.dumps(manifest))
    return path


@pytest.mark.parametrize("units,warn", [
    ("hartree", False), ("Eh", False), ("eV", False), ("Hartree", False),
    ("kJ/mol", True), ("EH", True), (12, True)])
def test_provenance_units(tmp_path, units, warn):
    path = _archive(tmp_path, provenance={
        "scf_energy": {"value": -1.0, "units": units},
        "total_energy": {"value": -1.0, "units": units}})
    report = validate_qvf(path)
    unit = [f for f in report["findings"] if f["code"] == "W_UNIT_UNKNOWN"]
    assert len(unit) == (2 if warn else 0), report["findings"]
    if warn:
        assert {f["where"] for f in unit} == {"/provenance/scf_energy/units",
                                              "/provenance/total_energy/units"}


def test_payload_units_with_and_without_the_schema_engine(tmp_path, monkeypatch):
    grid = {"origin": [0, 0, 0], "voxel_vectors": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "shape": [1, 1, 1], "units": {"data": "e_per_bohr3"}}
    path = _archive(tmp_path, kind="x_test.grid", role="grid", payload=grid)
    assert "W_UNIT_UNKNOWN" not in _codes(validate_qvf(path))   # vendor: free
    spectrum = {"frequencies": [1.0], "intensities": [1.0], "x_units": "nm",
                "y_units": "arb", "broadening": {"type": "none", "units": "1/cm"}}
    path = _archive(tmp_path, kind="spectra.uvvis", role="spectrum",
                    payload=spectrum)
    report = validate_qvf(path)
    unit = [f for f in report["findings"] if f["code"] == "W_UNIT_UNKNOWN"]
    assert len(unit) == 2
    assert all(f["where"] == "/sections/1/members/spectrum" for f in unit)
    pytest.importorskip("jsonschema")
    assert "E_PAYLOAD_SCHEMA" in _codes(report) and not report["ok"]
    monkeypatch.setattr(qvf_reader, "_load_schema_validator", lambda: None)
    monkeypatch.setattr(qvf_reader, "_load_payload_validator", lambda _d: None)
    report = validate_qvf(path)
    assert report["ok"]
    assert _codes(report) == ["W_SCHEMA_SKIPPED", "W_UNIT_UNKNOWN",
                              "W_UNIT_UNKNOWN"]


def test_fixed_unit_contradiction_is_a_schema_error(tmp_path):
    pytest.importorskip("jsonschema")
    meta = {"symbols": ["O"], "units": {"coords": "bohr"}}
    path = _archive(tmp_path, kind="trajectory", role="metadata", payload=meta)
    codes = _codes(validate_qvf(path))
    assert "E_PAYLOAD_SCHEMA" in codes and "W_UNIT_UNKNOWN" not in codes
