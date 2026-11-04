"""The property kinds ``properties.electric``, ``properties.elastic`` and
``properties.dielectric`` (spec § 5.13; P14, P17).

Covers the writers (round trips and refusals), the schema definitions, and the
validator rules that JSON Schema cannot express: symmetric elastic matrices,
``static`` = ``electronic`` + ``ionic``, one Born tensor per atom of the
section's structure, and the agreement of an electric dipole with the root
``dipole_moment`` (``E_ROOT_INCONSISTENT``).
"""

import hashlib
import json
import os
import sys
import zipfile

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))

import qvf_reader  # noqa: E402
from qvf_reader import QvfReader, validate_qvf  # noqa: E402
from qvf_writer import QvfWriter  # noqa: E402

_SCHEMA_CODES = {"E_SCHEMA", "E_PAYLOAD_SCHEMA", "W_SCHEMA_SKIPPED"}

WATER = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]
CELL = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]


def _identity(x):
    return [[x if i == j else 0.0 for j in range(3)] for i in range(3)]


def _cubic(c11=36.4, c12=8.8, c44=8.0):
    m = np.zeros((6, 6))
    for i in range(3):
        m[i, i], m[3 + i, 3 + i] = c11, c44
        m[i, :3] += c12
        m[i, i] = c11
    return m


def _codes(report):
    return {f["code"] for f in report["findings"]} - _SCHEMA_CODES


def _where(report, code):
    return sorted(f["where"] for f in report["findings"] if f["code"] == code)


# --- raw archives, so the validator sees what the writer would refuse ----------

def _json_bytes(obj):
    return json.dumps(obj).encode("utf-8")


class _Raw:
    """A manifest plus member files, written as a ZIP with the manifest last."""

    def __init__(self):
        self.manifest = {"qvf_version": 1,
                         "source": {"program": "t", "version": "1",
                                    "calculation": "properties"},
                         "sections": []}
        self.files = {}

    def structure(self, sid="structure", atoms=WATER):
        data = _json_bytes({"atoms": atoms, "pbc": [False] * 3,
                            "lattice_vectors": None})
        path = f"{sid}/structure.json"
        self.files[path] = data
        self.manifest["sections"].append(
            {"id": sid, "kind": "structure",
             "members": {"structure": {
                 "path": path, "format": "json",
                 "sha256": hashlib.sha256(data).hexdigest()}}})
        return self

    def properties(self, kind, payload, sid=None, **fields):
        sid = sid or kind.split(".")[1]
        data = _json_bytes(payload)
        path = f"{sid}/properties.json"
        self.files[path] = data
        self.manifest["sections"].append(
            {"id": sid, "kind": kind, **fields,
             "members": {"properties": {
                 "path": path, "format": "json",
                 "sha256": hashlib.sha256(data).hexdigest()}}})
        return self

    def validate(self, tmp_path, name="a.qvf"):
        path = str(tmp_path / name)
        with zipfile.ZipFile(path, "w") as zf:
            for p, data in self.files.items():
                zf.writestr(p, data)
            zf.writestr(zipfile.ZipInfo("manifest.json"), _json_bytes(self.manifest),
                        compress_type=zipfile.ZIP_STORED)
        return validate_qvf(path)


@pytest.fixture(params=["with-schema", "without-schema"])
def engine(request, monkeypatch):
    """The semantic rules hold with and without a JSON Schema engine."""
    if request.param == "without-schema":
        monkeypatch.setattr(qvf_reader, "_load_schema_validator", lambda: None)
        monkeypatch.setattr(qvf_reader, "_load_payload_validator", lambda name: None)
    else:
        pytest.importorskip("jsonschema")
    return request.param


# --- writer -------------------------------------------------------------------

def _writer(periodic=False):
    w = QvfWriter(program="t", version="1", calculation="properties")
    if periodic:
        w.add_structure(WATER[:2], pbc=[True] * 3, lattice_vectors=CELL)
    else:
        w.add_structure(WATER)
    return w


ELECTRIC = {
    "origin": "center_of_mass",
    "dipole": {"vector": [0.0, 0.0, -0.75], "total": 0.75, "units": "debye"},
    "quadrupole": {"matrix": np.diag([1.0, -0.4, -0.6]), "traceless": True,
                   "convention": "buckingham"},
    "octupole": {"tensor": np.arange(27.0).reshape(3, 3, 3)},
    "polarizability": {"matrix": np.diag([8.0, 9.0, 10.0]), "isotropic": 9.0},
    "hyperpolarizability": {"tensor": np.zeros(27), "frequencies_ev": [0.0, 0.0]},
    "units": {"dipole": "debye"},
}


def test_writer_round_trip(tmp_path):
    w = _writer(periodic=True)
    w.set_dipole_moment(total_debye=0.75, vector_debye=[0.0, 0.0, -0.75])
    w.add_electric_properties(ELECTRIC)
    stiffness = _cubic()
    w.add_elastic_properties(stiffness, compliance=np.linalg.inv(stiffness),
                             moduli={"bulk_modulus_gpa": 18.0,
                                     "poisson_ratio": 0.3},
                             average="hill", units={"stiffness": "gpa"})
    w.add_dielectric_properties(
        np.diag([2.0, 2.0, 2.0]), ionic=np.diag([1.0, 1.0, 1.0]),
        static=np.diag([3.0, 3.0, 3.0]), born_charges=np.zeros((2, 3, 3)),
        piezoelectric=np.zeros((3, 6)), units={"born_charges": "e"})
    path = str(tmp_path / "p.qvf")
    w.write(path)
    report = validate_qvf(path)
    assert report["ok"] and not _codes(report), report["findings"]
    assert [s["kind"] for s in w.build_manifest()["sections"]] == [
        "structure", "properties.electric", "properties.elastic",
        "properties.dielectric"]
    with QvfReader(path) as r:
        sections = {s["id"]: s for s in r.sections}
        electric = r.read_json(sections["electric"]["members"]["properties"])
        elastic = r.read_json(sections["elastic"]["members"]["properties"])
        dielectric = r.read_json(sections["dielectric"]["members"]["properties"])
    assert electric["octupole"]["tensor"][:4] == [0.0, 1.0, 2.0, 3.0]
    assert len(electric["octupole"]["tensor"]) == 27
    assert electric["quadrupole"]["matrix"][1][1] == -0.4
    assert elastic["stiffness"][0][1] == 8.8
    assert elastic["average"] == "hill" and elastic["poisson_ratio"] == 0.3
    assert "shear_modulus_gpa" not in elastic
    assert len(dielectric["born_charges"]) == 2
    assert dielectric["units"] == {"born_charges": "e"}


def test_writer_accepts_flat_or_nested_tensors_and_nested_lists(tmp_path):
    w = _writer()
    w.add_electric_properties(
        {"octupole": {"tensor": [float(k) for k in range(27)]},
         "polarizability": {"matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}},
        section_id="e1")
    w.add_electric_properties(
        {"octupole": {"tensor": np.arange(27.0).reshape(3, 3, 3)}}, section_id="e2")
    manifest = w.build_manifest()
    path = str(tmp_path / "t.qvf")
    w.write(path)
    with QvfReader(path) as r:
        a, b = (r.read_json(s["members"]["properties"])
                for s in r.sections if s["kind"] == "properties.electric")
    assert a["octupole"] == b["octupole"]
    assert manifest["sections"][1]["members"]["properties"]["path"] == \
        "e1/properties.json"


def test_writer_omits_optional_fields():
    w = _writer()
    w.add_elastic_properties(_cubic())
    w.add_dielectric_properties(_identity(2.0))
    w.add_electric_properties({"dipole": {"vector": [0, 0, 1]}})
    assert [s["id"] for s in w.build_manifest()["sections"]] == [
        "structure", "elastic", "dielectric", "electric"]
    for payload in (w._files["elastic/properties.json"],
                    w._files["dielectric/properties.json"],
                    w._files["electric/properties.json"]):
        assert not {"units", "x_note"} & set(json.loads(payload))


def test_writer_passes_vendor_keys_through():
    w = _writer()
    w.add_electric_properties({
        "dipole": {"vector": [0, 0, 1], "x_method": "ccsd"},
        "x_note": {"anything": [1, 2]}})
    payload = json.loads(w._files["electric/properties.json"])
    assert payload["dipole"]["x_method"] == "ccsd"
    assert payload["x_note"] == {"anything": [1, 2]}


def _bad(**overrides):
    return {**{"dipole": {"vector": [0, 0, 1]}}, **overrides}


@pytest.mark.parametrize("call, match", [
    (lambda w: w.add_electric_properties({}), "at least one of"),
    (lambda w: w.add_electric_properties({"origin": "com"}), "at least one of"),
    (lambda w: w.add_electric_properties(_bad(dipol=1)), "unknown keys"),
    (lambda w: w.add_electric_properties(_bad(origin=[0, 0])), "3 numbers"),
    (lambda w: w.add_electric_properties(_bad(origin="")), "must not be empty"),
    (lambda w: w.add_electric_properties({"dipole": {"total": 1.0}}),
     "requires a vector"),
    (lambda w: w.add_electric_properties({"dipole": {"vector": [0, 0, 1],
                                                     "units": "e_bohr"}}),
     "units must be 'debye'"),
    (lambda w: w.add_electric_properties({"dipole": {"vector": [0, 0, 1],
                                                     "total": -1.0}}), ">= 0"),
    (lambda w: w.add_electric_properties({"dipole": {"vector": [0, 0, 1],
                                                     "magnitude": 1}}),
     "unknown keys"),
    (lambda w: w.add_electric_properties({"dipole": {"vector": [0, 0, float("nan")]}}),
     "finite"),
    (lambda w: w.add_electric_properties({"quadrupole": {"matrix": np.eye(3),
                                                         "traceless": False}}),
     "requires convention"),
    (lambda w: w.add_electric_properties({"quadrupole": {
        "matrix": np.eye(3), "traceless": False, "convention": "buckingham"}}),
     "buckingham convention is traceless"),
    (lambda w: w.add_electric_properties({"quadrupole": {
        "matrix": np.eye(3), "traceless": False, "convention": "spherical"}}),
     "one of buckingham, cartesian"),
    (lambda w: w.add_electric_properties({"quadrupole": {
        "matrix": [[0, 1, 0], [2, 0, 0], [0, 0, 0]], "traceless": True,
        "convention": "cartesian"}}), "must be symmetric"),
    (lambda w: w.add_electric_properties({"quadrupole": {
        "matrix": np.eye(2), "traceless": False, "convention": "cartesian"}}),
     "3 x 3"),
    (lambda w: w.add_electric_properties({"quadrupole": {
        "matrix": np.eye(3), "traceless": "no", "convention": "cartesian"}}),
     "boolean"),
    (lambda w: w.add_electric_properties({"octupole": {"tensor": np.zeros(26)}}),
     "27 numbers"),
    (lambda w: w.add_electric_properties({"octupole": {"tensor": np.zeros(27),
                                                       "convention": "x"}}),
     "convention"),
    (lambda w: w.add_electric_properties({"polarizability": {
        "matrix": np.eye(3), "frequency_ev": -1.0}}), ">= 0"),
    (lambda w: w.add_electric_properties({"polarizability": {
        "matrix": np.eye(3), "isotropic": "1"}}), "must be a number"),
    (lambda w: w.add_electric_properties({"hyperpolarizability": {
        "tensor": np.zeros(27), "frequencies_ev": [1.0]}}), "2 numbers"),
    (lambda w: w.add_electric_properties(_bad(units={"dipole": "e"})),
     "must be 'debye'"),
    (lambda w: w.add_electric_properties(_bad(units={"dipole": "hartree"})),
     "must be 'debye'"),
    (lambda w: w.add_electric_properties(_bad(units={"dipole": "D"})),
     "not a QVF unit name"),
    (lambda w: w.add_electric_properties("dipole"), "must be a mapping"),
    (lambda w: w.add_elastic_properties(np.zeros((5, 6))), "6 x 6"),
    (lambda w: w.add_elastic_properties(
        [[0, 1, 0, 0, 0, 0], [2, 0, 0, 0, 0, 0]] + [[0] * 6] * 4),
     "must be symmetric"),
    (lambda w: w.add_elastic_properties(_cubic(), compliance=np.ones((6, 5))),
     "6 x 6"),
    (lambda w: w.add_elastic_properties(_cubic(), moduli={"bulk": 1.0}),
     "unknown keys"),
    (lambda w: w.add_elastic_properties(_cubic(), moduli=[1.0]), "mapping"),
    (lambda w: w.add_elastic_properties(_cubic(), average="mean"),
     "one of voigt, reuss, hill"),
    (lambda w: w.add_elastic_properties(_cubic(), units={"stiffness": "ev"}),
     "must be 'gpa'"),
    (lambda w: w.add_dielectric_properties(np.ones((3, 2))), "3 x 3"),
    (lambda w: w.add_dielectric_properties(
        _identity(2.0), ionic=_identity(1.0), static=_identity(3.5)),
     "static must equal their sum"),
    (lambda w: w.add_dielectric_properties(_identity(2.0), born_charges=np.zeros((3, 3))),
     r"\[n_atoms, 3, 3\]"),
    (lambda w: w.add_dielectric_properties(_identity(2.0),
                                           born_charges=np.zeros((0, 3, 3))),
     r"\[n_atoms, 3, 3\]"),
    (lambda w: w.add_dielectric_properties(_identity(2.0), piezoelectric=np.zeros((3, 3))),
     "3 x 6"),
    (lambda w: w.add_dielectric_properties(_identity(2.0), units={"born_charges": "hartree"}),
     "must be 'e'"),
])
def test_writer_refuses_before_writing(call, match):
    w = _writer()
    before = (len(w._sections), dict(w._files))
    with pytest.raises(ValueError, match=match):
        call(w)
    assert (len(w._sections), w._files) == before
    assert len(w.build_manifest()["sections"]) == 1


def test_writer_checks_the_section_id():
    w = _writer()
    w.add_elastic_properties(_cubic(), section_id="c")
    with pytest.raises(ValueError, match="duplicate section id"):
        w.add_elastic_properties(_cubic(), section_id="c")
    with pytest.raises(ValueError, match="must match"):
        w.add_elastic_properties(_cubic(), section_id="bad id")


def test_writer_accepts_a_symmetric_matrix_within_tolerance():
    stiff = _cubic()
    stiff[0, 1] += 1e-9          # relative 3e-11 of the largest entry
    _writer().add_elastic_properties(stiff)


# --- schema -------------------------------------------------------------------

def _payload_errors(defname, payload):
    pytest.importorskip("jsonschema")
    return list(qvf_reader._load_payload_validator(defname).iter_errors(payload))


def _schema_errors(manifest):
    jsonschema = pytest.importorskip("jsonschema")
    with open(qvf_reader._SCHEMA_PATH, encoding="utf-8") as fh:
        validator = jsonschema.Draft202012Validator(json.load(fh))
    return list(validator.iter_errors(manifest))


M3 = _identity(1.0)
T27 = [0.0] * 27
C6 = _cubic().tolist()


@pytest.mark.parametrize("payload, ok", [
    ({"dipole": {"vector": [0, 0, 1]}}, True),
    ({"quadrupole": {"matrix": M3, "traceless": False,
                     "convention": "cartesian"}}, True),
    ({"octupole": {"tensor": T27}}, True),
    ({"polarizability": {"matrix": M3, "frequency_ev": 0}}, True),
    ({"hyperpolarizability": {"tensor": T27}}, True),
    ({"dipole": {"vector": [0, 0, 1]}, "origin": [0.0, 0.0, 1.5],
      "x_note": "n", "units": {"dipole": "debye", "origin": "bohr"}}, True),
    ({"dipole": {"vector": [0, 0, 1], "x_basis": "cc-pvdz"}}, True),
    ({}, False),                                           # nothing to carry
    ({"origin": "origin"}, False),
    ({"dipole": {"total": 1.0}}, False),                   # no vector
    ({"dipole": {"vector": [0, 1]}}, False),
    ({"dipole": {"vector": [0, 0, 1], "total": -1}}, False),
    ({"dipole": {"vector": [0, 0, 1], "units": "e"}}, False),
    ({"dipole": {"vector": [0, 0, 1], "extra": 1}}, False),
    ({"dipole": {"vector": [0, 0, 1]}, "quadropole": {}}, False),
    ({"dipole": {"vector": [0, 0, 1]}, "units": {"dipole": "e"}}, False),
    ({"dipole": {"vector": [0, 0, 1]}, "units": {"quadrupole": "debye"}}, False),
    ({"dipole": {"vector": [0, 0, 1]}, "origin": ""}, False),
    ({"dipole": {"vector": [0, 0, 1]}, "origin": [0, 0]}, False),
    ({"quadrupole": {"matrix": M3, "traceless": True}}, False),
    ({"quadrupole": {"matrix": M3, "convention": "cartesian"}}, False),
    ({"quadrupole": {"matrix": M3, "traceless": True,
                     "convention": "traceless"}}, False),
    ({"quadrupole": {"matrix": M3, "traceless": False,
                     "convention": "buckingham"}}, False),
    ({"quadrupole": {"matrix": M3[:2], "traceless": True,
                     "convention": "cartesian"}}, False),
    ({"octupole": {"tensor": T27[:26]}}, False),
    ({"octupole": {"tensor": T27, "convention": "spherical"}}, False),
    ({"polarizability": {"matrix": M3, "frequency_ev": -0.1}}, False),
    ({"polarizability": {"frequency_ev": 0.0}}, False),
    ({"hyperpolarizability": {"tensor": T27, "frequencies_ev": [1.0]}}, False),
    ({"hyperpolarizability": {"tensor": T27, "frequencies_ev": [1, 2, 3]}}, False),
])
def test_electric_payload_definition(payload, ok):
    assert (not _payload_errors("ElectricPropertiesPayload", payload)) is ok


@pytest.mark.parametrize("payload, ok", [
    ({"stiffness": C6}, True),
    ({"stiffness": C6, "compliance": C6, "bulk_modulus_gpa": 1, "shear_modulus_gpa": 1,
      "youngs_modulus_gpa": 1, "poisson_ratio": 0.3, "average": "reuss",
      "units": {"stiffness": "gpa", "poisson_ratio": "dimensionless"},
      "x_code": 1}, True),
    ({}, False),
    ({"compliance": C6}, False),
    ({"stiffness": C6[:5]}, False),
    ({"stiffness": [row[:5] for row in C6]}, False),
    ({"stiffness": C6, "average": "mean"}, False),
    ({"stiffness": C6, "bulk": 1}, False),
    ({"stiffness": C6, "bulk_modulus_gpa": "1"}, False),
    ({"stiffness": C6, "units": {"stiffness": "ev"}}, False),
    ({"stiffness": C6, "units": {"stiffness": "GPa"}}, False),
    ({"stiffness": C6, "units": {"poisson_ratio": "gpa"}}, False),
    ({"stiffness": C6, "units": {"compliance": "gpa"}}, False),
])
def test_elastic_payload_definition(payload, ok):
    assert (not _payload_errors("ElasticPropertiesPayload", payload)) is ok


@pytest.mark.parametrize("payload, ok", [
    ({"electronic": M3}, True),
    ({"electronic": M3, "ionic": M3, "static": M3, "born_charges": [M3, M3],
      "piezoelectric": [[0.0] * 6] * 3,
      "units": {"electronic": "dimensionless", "static": "dimensionless",
                "born_charges": "e"}, "x_k": 0}, True),
    ({}, False),
    ({"ionic": M3}, False),
    ({"electronic": M3[:2]}, False),
    ({"electronic": M3, "born_charges": []}, False),
    ({"electronic": M3, "born_charges": [M3[:2]]}, False),
    ({"electronic": M3, "piezoelectric": M3}, False),
    ({"electronic": M3, "piezoelectric": [[0.0] * 6] * 2}, False),
    ({"electronic": M3, "dielectric": M3}, False),
    ({"electronic": M3, "units": {"born_charges": "hartree"}}, False),
    ({"electronic": M3, "units": {"electronic": "e"}}, False),
    ({"electronic": M3, "units": {"piezoelectric": "dimensionless"}}, False),
])
def test_dielectric_payload_definition(payload, ok):
    assert (not _payload_errors("DielectricPropertiesPayload", payload)) is ok


def test_section_definitions_in_the_manifest():
    w = _writer()
    w.add_electric_properties({"dipole": {"vector": [0, 0, 1]}})
    manifest = json.loads(json.dumps(w.build_manifest()))
    assert not _schema_errors(manifest)
    sec = manifest["sections"][1]
    for mutate in (lambda s: s["members"].update(extra=s["members"]["properties"]),
                   lambda s: s["members"].pop("properties"),
                   lambda s: s["members"]["properties"].update(format="binary")):
        broken = json.loads(json.dumps(manifest))
        mutate(broken["sections"][1])
        assert _schema_errors(broken)
    assert sec["kind"] == "properties.electric"


def test_payload_is_checked_by_the_validator_with_a_schema_engine(tmp_path):
    pytest.importorskip("jsonschema")
    raw = _Raw().structure().properties(
        "properties.elastic", {"stiffness": C6, "average": "mean"})
    report = raw.validate(tmp_path)
    assert "E_PAYLOAD_SCHEMA" in {f["code"] for f in report["findings"]}
    assert not report["ok"]


# --- validator: elastic symmetry ----------------------------------------------

def test_valid_property_archive_is_clean(tmp_path, engine):
    raw = (_Raw().structure()
           .properties("properties.elastic", {"stiffness": C6, "compliance": C6})
           .properties("properties.dielectric", {
               "electronic": M3, "ionic": M3, "static": _identity(2.0),
               "born_charges": [M3] * 3})
           .properties("properties.electric", {"dipole": {"vector": [0, 0, 1]}}))
    report = raw.validate(tmp_path)
    assert report["ok"] and not _codes(report), report["findings"]


@pytest.mark.parametrize("name", ["stiffness", "compliance"])
def test_elastic_matrices_must_be_symmetric(tmp_path, engine, name):
    bad = _cubic().tolist()
    bad[2][4] = 3.0
    payload = {"stiffness": C6, name: bad}
    report = _Raw().structure().properties("properties.elastic",
                                           payload).validate(tmp_path)
    assert _codes(report) == {"E_MEMBER_MALFORMED"}
    assert _where(report, "E_MEMBER_MALFORMED") == ["/sections/1/members/properties"]
    assert name in report["errors"][0]


def test_elastic_symmetry_tolerance_is_relative_to_the_largest_entry(tmp_path, engine):
    near = _cubic().tolist()
    near[0][1] += 3.0e-5                       # 8.3e-7 of 36.4: within tolerance
    assert not _codes(_Raw().structure().properties(
        "properties.elastic", {"stiffness": near}).validate(tmp_path, "near.qvf"))
    far = _cubic().tolist()
    far[0][1] += 5.0e-5                        # 1.4e-6 of 36.4: outside
    assert _codes(_Raw().structure().properties(
        "properties.elastic", {"stiffness": far}).validate(tmp_path, "far.qvf")
    ) == {"E_MEMBER_MALFORMED"}
    zero = [[0.0] * 6 for _ in range(6)]       # an all-zero matrix is symmetric
    assert not _codes(_Raw().structure().properties(
        "properties.elastic", {"stiffness": zero}).validate(tmp_path, "zero.qvf"))


@pytest.mark.parametrize("stiffness", [
    [[1.0] * 6] * 5, [[1.0] * 5] * 6, "stiff", None, [[1.0, "x"] + [0.0] * 4] * 6])
def test_malformed_elastic_matrices_are_left_to_the_schema(tmp_path, stiffness):
    report = _Raw().structure().properties(
        "properties.elastic", {"stiffness": stiffness}).validate(tmp_path)
    assert not _codes(report)


# --- validator: dielectric ----------------------------------------------------

def test_static_must_equal_electronic_plus_ionic(tmp_path, engine):
    payload = {"electronic": _identity(2.62), "ionic": _identity(4.58),
               "static": _identity(7.5)}
    report = _Raw().structure().properties("properties.dielectric",
                                           payload).validate(tmp_path)
    assert _codes(report) == {"E_MEMBER_MALFORMED"}
    assert _where(report, "E_MEMBER_MALFORMED") == ["/sections/1/members/properties"]


@pytest.mark.parametrize("delta, ok", [(5e-7, True), (2e-6, False)])
def test_static_tolerance_is_one_millionth(tmp_path, delta, ok):
    static = _identity(7.2)
    static[1][2] = delta
    payload = {"electronic": _identity(2.62), "ionic": _identity(4.58),
               "static": static}
    report = _Raw().structure().properties("properties.dielectric",
                                           payload).validate(tmp_path)
    assert (not _codes(report)) is ok


@pytest.mark.parametrize("present", [("electronic", "static"), ("electronic", "ionic"),
                                     ("electronic",)])
def test_static_is_compared_only_when_all_three_are_present(tmp_path, present):
    full = {"electronic": _identity(1.0), "ionic": _identity(1.0),
            "static": _identity(9.0)}
    payload = {k: full[k] for k in present}
    report = _Raw().structure().properties("properties.dielectric",
                                           payload).validate(tmp_path)
    assert not _codes(report)


def test_born_charges_must_number_the_atoms(tmp_path, engine):
    def report(n, name):
        return (_Raw().structure().properties(
            "properties.dielectric",
            {"electronic": M3, "born_charges": [M3] * n}).validate(tmp_path, name))
    assert not _codes(report(3, "three.qvf"))
    bad = report(2, "two.qvf")
    assert _codes(bad) == {"E_MEMBER_MALFORMED"}
    assert _where(bad, "E_MEMBER_MALFORMED") == ["/sections/1/members/properties"]
    assert "2 tensors" in bad["errors"][0] and "3 atoms" in bad["errors"][0]


def test_born_charges_use_the_section_structure_ref(tmp_path):
    def raw(n, **fields):
        return (_Raw().structure().structure("dimer", WATER[:2]).properties(
            "properties.dielectric", {"electronic": M3, "born_charges": [M3] * n},
            **fields))
    assert not _codes(raw(2, structure_ref="dimer").validate(tmp_path, "a.qvf"))
    assert _codes(raw(3, structure_ref="dimer").validate(tmp_path, "b.qvf")) == {
        "E_MEMBER_MALFORMED"}
    # Without a reference, two structures leave the count unchecked.
    assert not _codes(raw(7).validate(tmp_path, "c.qvf"))
    # A reference to something that is not a structure is a reference error.
    report = raw(3, structure_ref="nowhere").validate(tmp_path, "d.qvf")
    assert _codes(report) == {"E_REF_UNRESOLVED"}


def test_born_charges_before_the_structure_are_still_checked(tmp_path):
    raw = _Raw()
    raw.properties("properties.dielectric",
                   {"electronic": M3, "born_charges": [M3] * 5})
    raw.structure()
    report = raw.validate(tmp_path)
    assert _codes(report) == {"E_MEMBER_MALFORMED"}
    assert _where(report, "E_MEMBER_MALFORMED") == ["/sections/0/members/properties"]


def test_born_charges_without_any_structure_are_not_counted(tmp_path):
    report = _Raw().properties("properties.dielectric",
                               {"electronic": M3, "born_charges": [M3] * 5}
                               ).validate(tmp_path)
    assert not _codes(report)


# --- validator: the electric dipole against the root block ---------------------

def _electric_report(tmp_path, dipole, root, name="e.qvf"):
    raw = _Raw().structure().properties("properties.electric", {"dipole": dipole})
    if root is not None:
        raw.manifest["dipole_moment"] = root
    return raw.validate(tmp_path, name)


ROOT = {"total_debye": 1.85, "vector_debye": [0.0, 0.0, -1.85],
        "origin": "center_of_mass"}


def test_matching_dipole_is_consistent(tmp_path, engine):
    report = _electric_report(
        tmp_path, {"vector": [0.0, 0.0, -1.85 + 5e-7], "total": 1.85 - 5e-7}, ROOT)
    assert report["ok"] and not _codes(report), report["findings"]


@pytest.mark.parametrize("dipole", [
    {"vector": [0.0, 0.0, -1.5]},
    {"vector": [0.0, 0.1, -1.85]},
    {"vector": [0.0, 0.0, -1.85 + 2e-6]},
    {"vector": [0.0, 0.0, -1.85], "total": 1.5},
])
def test_disagreeing_dipole_is_root_inconsistent(tmp_path, engine, dipole):
    report = _electric_report(tmp_path, dipole, ROOT)
    assert _codes(report) == {"E_ROOT_INCONSISTENT"}, report["findings"]
    assert _where(report, "E_ROOT_INCONSISTENT") == ["/sections/1/members/properties"]
    assert not report["ok"]


def test_both_vector_and_total_can_disagree(tmp_path):
    report = _electric_report(tmp_path, {"vector": [0, 0, -1.5], "total": 1.5}, ROOT)
    assert len(_where(report, "E_ROOT_INCONSISTENT")) == 2


@pytest.mark.parametrize("root", [
    None,
    {},
    {"origin": "origin"},
    {"total_debye": 1.85},                        # only the total, not compared
    {"vector_debye": [0.0, 0.0]},                 # malformed: left to the schema
    {"vector_debye": [0.0, 0.0, "x"]},
    "dipole",
])
def test_nothing_to_compare_means_no_finding(tmp_path, root):
    report = _electric_report(tmp_path, {"vector": [0.0, 0.0, -1.5]}, root)
    assert "E_ROOT_INCONSISTENT" not in {f["code"] for f in report["findings"]}


def test_a_dipole_without_a_total_compares_only_the_vector(tmp_path):
    report = _electric_report(tmp_path, {"vector": [0.0, 0.0, -1.85]},
                              {"total_debye": 99.0})
    assert not _codes(report)


def test_every_electric_section_is_compared(tmp_path):
    raw = (_Raw().structure()
           .properties("properties.electric", {"dipole": {"vector": [0, 0, -1.85]}},
                       sid="static")
           .properties("properties.electric", {"dipole": {"vector": [0, 0, 1.0]}},
                       sid="dynamic"))
    raw.manifest["dipole_moment"] = ROOT
    report = raw.validate(tmp_path)
    assert _where(report, "E_ROOT_INCONSISTENT") == ["/sections/2/members/properties"]


# --- units --------------------------------------------------------------------

def test_unknown_declared_unit_is_a_warning_everywhere(tmp_path, monkeypatch):
    monkeypatch.setattr(qvf_reader, "_load_schema_validator", lambda: None)
    monkeypatch.setattr(qvf_reader, "_load_payload_validator", lambda name: None)
    report = _Raw().structure().properties(
        "properties.elastic",
        {"stiffness": C6, "units": {"stiffness": "GPa"}}).validate(tmp_path)
    assert report["ok"]
    assert [f["code"] for f in report["findings"]
            if f["code"] == "W_UNIT_UNKNOWN"] == ["W_UNIT_UNKNOWN"]


def test_schema_fixes_the_declarable_units():
    fixed = {
        "ElectricPropertiesPayload": {"dipole": "debye", "origin": "bohr",
                                      "quadrupole": "debye_angstrom",
                                      "octupole": "debye_angstrom2",
                                      "polarizability": "e2_bohr2_per_hartree",
                                      "hyperpolarizability": "e3_bohr3_per_hartree2"},
        "ElasticPropertiesPayload": {"stiffness": "gpa", "compliance": "gpa-1",
                                     "bulk_modulus_gpa": "gpa",
                                     "shear_modulus_gpa": "gpa",
                                     "youngs_modulus_gpa": "gpa",
                                     "poisson_ratio": "dimensionless"},
        "DielectricPropertiesPayload": {"electronic": "dimensionless",
                                        "ionic": "dimensionless",
                                        "static": "dimensionless",
                                        "born_charges": "e",
                                        "piezoelectric": "c_per_m2"},
    }
    from qvf_units import fixed_units
    for defname, expected in fixed.items():
        assert fixed_units(defname) == expected


# --- registry, kinds and corpus -----------------------------------------------

def test_kinds_are_canonical_and_registered():
    root = os.path.join(HERE, "..", "..")
    with open(os.path.join(root, "registry.json"), encoding="utf-8") as fh:
        registry = json.load(fh)
    kinds = {"properties.electric", "properties.elastic", "properties.dielectric"}
    assert kinds <= set(registry["canonical_kinds"])
    assert kinds <= qvf_reader._CANONICAL_KINDS
    promoted = {p["kind"]: p for p in registry["promotions"]}
    for kind in kinds:
        assert promoted[kind]["from"] is None
    assert {k: v for k, v in qvf_reader._PAYLOAD_DEFS.items()
            if k[0] in kinds} == {
        ("properties.electric", "properties"): "ElectricPropertiesPayload",
        ("properties.elastic", "properties"): "ElasticPropertiesPayload",
        ("properties.dielectric", "properties"): "DielectricPropertiesPayload"}
