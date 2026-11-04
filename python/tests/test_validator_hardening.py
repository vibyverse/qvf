"""Validator hardening (spec § 2.1, § 6) and stable finding codes (Appendix B).

Covers what the negative conformance corpus cannot: the report shape and its
backward-compatible string lists, the linear JSON depth scanner, limits as
parameters, caps checked *before* decompression, the reader's refusals, the
CLI, and drift guards that keep Appendix B, the Python code table, and the C++
validator in lock-step.
"""

import json
import os
import re
import subprocess
import sys
import zipfile

import pytest

HERE = os.path.dirname(__file__)
PYDIR = os.path.normpath(os.path.join(HERE, ".."))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PYDIR)

import qvf_reader  # noqa: E402
from qvf_reader import (  # noqa: E402
    FINDING_CODES, Limits, QvfReader, _json_depth, _path_problem, validate_qvf)
from qvf_writer import QvfWriter  # noqa: E402

WATER = {
    "atoms": [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8}],
    "pbc": [False, False, False],
    "lattice_vectors": None,
}


def _member(path, data, fmt="json"):
    import hashlib
    return {"path": path, "format": fmt, "sha256": hashlib.sha256(data).hexdigest()}


def _write_zip(path, entries):
    """entries: list of (name, bytes[, compress_type]); manifest STOREd."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")      # duplicate names are the point
        with zipfile.ZipFile(path, "w") as zf:
            for entry in entries:
                name, data = entry[0], entry[1]
                method = entry[2] if len(entry) > 2 else zipfile.ZIP_DEFLATED
                zf.writestr(zipfile.ZipInfo(name), data, compress_type=method)
    return str(path)


def _archive(tmp_path, *, sections=None, files=None, extra_entries=(),
             manifest_extra=None, name="a.qvf"):
    """A minimal valid archive (one water structure) with optional changes."""
    sdata = json.dumps(WATER).encode()
    if sections is None:
        sections = [{"id": "structure", "kind": "structure",
                     "members": {"structure": _member("structure/structure.json",
                                                      sdata)}}]
        files = {"structure/structure.json": sdata, **(files or {})}
    manifest = {"qvf_version": 1,
                "source": {"program": "t", "version": "1", "calculation": "c"},
                "sections": sections, **(manifest_extra or {})}
    entries = [(n, d) for n, d in (files or {}).items()] + list(extra_entries)
    entries.append(("manifest.json", json.dumps(manifest).encode(),
                    zipfile.ZIP_STORED))
    return _write_zip(tmp_path / name, entries)


def _codes(report):
    return {f["code"] for f in report["findings"]}


# --- report shape --------------------------------------------------------------

def test_report_keeps_legacy_keys_and_adds_findings(tmp_path):
    report = validate_qvf(_archive(tmp_path))
    assert set(report) == {"ok", "errors", "warnings", "n_sections", "findings", "validation_context"}
    assert report["ok"] and report["errors"] == [] and report["n_sections"] == 1
    for f in report["findings"]:
        assert set(f) == {"severity", "code", "where", "message"}


def test_string_lists_are_derived_from_findings(tmp_path):
    sdata = json.dumps(WATER).encode()
    bad = _member("structure/structure.json", b"other bytes")
    path = _archive(tmp_path, sections=[{"id": "structure", "kind": "structure",
                                         "members": {"structure": bad}}],
                    files={"structure/structure.json": sdata},
                    extra_entries=[("stale.json", b"{}")])
    report = validate_qvf(path)
    assert not report["ok"]
    errors = [f for f in report["findings"] if f["severity"] == "error"]
    warnings = [f for f in report["findings"] if f["severity"] == "warning"]
    assert report["errors"] == [f"{f['code']}: {f['message']}" for f in errors]
    assert report["warnings"] == [f"{f['code']}: {f['message']}" for f in warnings]
    assert "E_SHA256_MISMATCH" in _codes(report)
    assert any(w.startswith("W_ORPHAN_MEMBER: ") for w in report["warnings"])
    sha = next(f for f in errors if f["code"] == "E_SHA256_MISMATCH")
    assert sha["where"] == "/sections/0/members/structure"


def test_every_emitted_code_is_registered(tmp_path):
    neg = os.path.join(REPO, "conformance", "negative")
    if not os.path.isdir(neg):
        pytest.skip("negative corpus not present")
    for name in sorted(os.listdir(neg)):
        if name.endswith(".qvf"):
            for f in validate_qvf(os.path.join(neg, name))["findings"]:
                assert FINDING_CODES.get(f["code"]) == f["severity"], (name, f)


def test_schema_skipped_is_a_warning(tmp_path, monkeypatch):
    monkeypatch.setattr(qvf_reader, "_load_schema_validator", lambda: None)
    monkeypatch.setattr(qvf_reader, "_load_payload_validator", lambda name: None)
    report = validate_qvf(_archive(tmp_path))
    assert report["ok"]
    assert _codes(report) == {"W_SCHEMA_SKIPPED"}


def test_schema_errors_carry_a_pointer(tmp_path):
    pytest.importorskip("jsonschema")
    report = validate_qvf(_archive(tmp_path, manifest_extra={"unknown_root": 1}))
    schema = [f for f in report["findings"] if f["code"] == "E_SCHEMA"]
    assert schema and all(f["where"] == "manifest.json" for f in schema)


# --- JSON depth scanner ------------------------------------------------------------

@pytest.mark.parametrize("raw, depth", [
    (b"1", 0),
    (b"[]", 1),
    (b'{"a": [1, {"b": []}]}', 4),
    (b'{"k": "[[[[{{{{"}', 1),               # brackets inside a string
    (b'["\\"[[["]', 1),                      # escaped quote keeps the string open
    (b'["\\\\", [[]]]', 3),                  # escaped backslash closes it
    (b'"\\u005b\\u005b"', 0),
])
def test_json_depth_scanner(raw, depth):
    assert _json_depth(raw, 64) == depth


def test_json_depth_scanner_stops_past_the_limit():
    assert _json_depth(b"[" * 100000, 64) == 65


def test_json_depth_boundary_is_inclusive(tmp_path):
    def with_blob(depth):
        blob = b"[" * depth + b"]" * depth
        sdata = json.dumps(WATER).encode()
        return _archive(tmp_path, name=f"d{depth}.qvf", sections=[
            {"id": "structure", "kind": "structure",
             "members": {"structure": _member("structure/structure.json", sdata)}},
            {"id": "blob", "kind": "x_test.blob",
             "members": {"data": _member("x/blob.json", blob)}}],
            files={"structure/structure.json": sdata, "x/blob.json": blob})
    assert validate_qvf(with_blob(64))["ok"]
    assert "E_JSON_TOO_DEEP" in _codes(validate_qvf(with_blob(65)))
    assert "E_JSON_TOO_DEEP" in _codes(
        validate_qvf(with_blob(10), limits=Limits(max_json_depth=9)))


# --- path rules and limits -----------------------------------------------------------

@pytest.mark.parametrize("path, bad", [
    ("structure/structure.json", False),
    ("a..b/c.json", False),
    ("dir/.../c.json", False),
    ("dir\\c.json", True),
    ("C:/c.json", True),
    ("z:c.json", True),
    ("/c.json", True),
    ("../c.json", True),
    ("a/../c.json", True),
    ("a/..", True),
    ("a\0b", True),
])
def test_path_rules(path, bad):
    assert (_path_problem(path) is not None) == bad


def test_limits_from_dict():
    assert Limits.from_dict(None) == Limits()
    assert Limits.from_dict({"max_ratio": 5}).max_ratio == 5
    with pytest.raises(TypeError):
        Limits.from_dict({"max_bogus": 1})


def test_caps_are_checked_before_decompression(tmp_path, monkeypatch):
    zeros = bytes(1 << 16)
    sdata = json.dumps(WATER).encode()
    path = _archive(tmp_path, sections=[
        {"id": "structure", "kind": "structure",
         "members": {"structure": _member("structure/structure.json", sdata)}},
        {"id": "zeros", "kind": "x_test.zeros",
         "members": {"data": dict(_member("x/zeros.bin", zeros, "binary"),
                                  dtype="uint8", shape=[len(zeros)])}}],
        files={"structure/structure.json": sdata, "x/zeros.bin": zeros})
    read = []
    real_open = zipfile.ZipFile.open

    def spy(self, name, mode="r", pwd=None, **kwargs):
        read.append(getattr(name, "filename", name))
        return real_open(self, name, mode, pwd, **kwargs)

    # Every entry is decompressed through ZipFile.open (bounded reads).
    monkeypatch.setattr(zipfile.ZipFile, "open", spy)
    for lim, code in ((Limits(max_ratio=10, ratio_exempt_bytes=1024),
                       "E_COMPRESSION_RATIO"),
                      (Limits(max_member_bytes=1024), "E_MEMBER_TOO_LARGE"),
                      (Limits(max_total_bytes=1024), "E_TOTAL_TOO_LARGE")):
        read.clear()
        report = validate_qvf(path, limits=lim)
        assert code in _codes(report), report["findings"]
        assert "x/zeros.bin" not in read, (code, read)
    read.clear()
    assert validate_qvf(path)["ok"]          # within the default caps
    assert "x/zeros.bin" in read


def test_json_size_cap(tmp_path):
    path = _archive(tmp_path)
    report = validate_qvf(path, limits=Limits(max_json_bytes=16))
    assert [f["code"] for f in report["findings"]] == ["E_JSON_TOO_LARGE"]
    assert report["findings"][0]["where"] == "manifest.json"


def test_reserved_ids(tmp_path):
    sdata = json.dumps(WATER).encode()
    for sid in ("auto_open", "bookmarks"):
        path = _archive(tmp_path, name=f"{sid}.qvf", sections=[
            {"id": sid, "kind": "structure",
             "members": {"structure": _member("structure/structure.json", sdata)}}],
            files={"structure/structure.json": sdata})
        assert "E_RESERVED_ID" in _codes(validate_qvf(path))


def test_orphans_ignore_directory_entries(tmp_path):
    path = _archive(tmp_path, extra_entries=[("stale/", b""), ("stale/x.json", b"{}")])
    report = validate_qvf(path)
    assert report["ok"]
    orphans = [f["where"] for f in report["findings"] if f["code"] == "W_ORPHAN_MEMBER"]
    assert orphans == ["entry:stale/x.json"]


def test_malformed_manifest_shapes_do_not_crash(tmp_path):
    cases = {
        "list.qvf": b"[1, 2]",
        "sections.qvf": b'{"qvf_version": 1, "source": [], "sections": {}}',
        "section.qvf": b'{"qvf_version": 1, "source": {}, "sections": [1, {"id": [1]}]}',
        "refs.qvf": (b'{"qvf_version": true, "source": {"program": "p", '
                     b'"version": "v", "calculation": ""}, "sections": [{"id": "w", '
                     b'"kind": "reaction.waypoints", "trajectory_ref": [1], '
                     b'"members": {"m": {"path": 5}}}]}'),
    }
    for name, manifest in cases.items():
        path = _write_zip(tmp_path / name, [("manifest.json", manifest,
                                             zipfile.ZIP_STORED)])
        report = validate_qvf(path)
        assert not report["ok"], name
        assert _codes(report) <= set(FINDING_CODES), name


# --- reader refusals -------------------------------------------------------------

def test_reader_refuses_duplicate_entries(tmp_path):
    sdata = json.dumps(WATER).encode()
    path = _archive(tmp_path, extra_entries=[("structure/structure.json", sdata)])
    assert "E_ZIP_DUPLICATE_ENTRY" in _codes(validate_qvf(path))
    with pytest.raises(ValueError, match="E_ZIP_DUPLICATE_ENTRY"):
        QvfReader(path)


def test_reader_applies_manifest_caps(tmp_path):
    path = _archive(tmp_path)
    with pytest.raises(ValueError, match="E_JSON_TOO_LARGE"):
        QvfReader(path, limits=Limits(max_json_bytes=16))
    with pytest.raises(ValueError, match="E_JSON_TOO_DEEP"):
        QvfReader(path, limits=Limits(max_json_depth=2))
    with QvfReader(path) as r:
        assert r.sections[0]["id"] == "structure"


# --- writers refuse what validators reject ------------------------------------------

def test_writer_refuses_reserved_ids_and_bad_paths():
    w = QvfWriter(program="t", version="1")
    for sid in ("auto_open", "bookmarks"):
        with pytest.raises(ValueError, match="reserved"):
            w.add_citations("@misc{x}", section_id=sid)
    with pytest.raises(ValueError, match="illegal member path"):
        w.add_vendor_section("x_test.thing", json_members={"a\\b": {}})


# --- CLI ------------------------------------------------------------------------

def _cli(*args):
    return subprocess.run([sys.executable, os.path.join(PYDIR, "qvf_reader.py"),
                           *args], capture_output=True, text=True)


def test_cli_prints_codes_and_keeps_exit_codes(tmp_path):
    good = _archive(tmp_path, name="good.qvf")
    bad = _archive(tmp_path, name="bad.qvf", sections=[
        {"id": "bookmarks", "kind": "structure",
         "members": {"structure": _member("structure/structure.json",
                                          json.dumps(WATER).encode())}}],
        files={"structure/structure.json": json.dumps(WATER).encode()})
    assert _cli(good).returncode == 0
    r = _cli(bad)
    assert r.returncode == 1
    assert re.search(r"^\s+E_RESERVED_ID: ", r.stdout, re.M), r.stdout
    assert _cli().returncode == 2
    j = _cli("--json", bad)
    assert j.returncode == 1
    doc = json.loads(j.stdout)
    assert doc["ok"] is False
    assert any(f["code"] == "E_RESERVED_ID" for f in doc["findings"])


# --- drift guards: Appendix B ≡ Python ≡ C++ ---------------------------------------

def _repo_file(*parts):
    path = os.path.join(REPO, *parts)
    if not os.path.exists(path):
        pytest.skip(f"{os.path.join(*parts)} not present (standalone checkout)")
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def test_appendix_b_matches_finding_codes():
    spec = _repo_file("spec", "qvf-format-spec.md")
    appendix = spec.split("## Appendix B.", 1)[1]
    rows = dict(re.findall(r"^\| `([EW]_[A-Z0-9_]+)` \| (error|warning) \|",
                           appendix, re.M))
    assert rows == FINDING_CODES
    for code, severity in rows.items():
        assert code.startswith("E_" if severity == "error" else "W_"), code


def test_cpp_validator_emits_the_same_codes():
    src = _repo_file("cpp", "src", "qvf_reader.cpp")
    cpp_codes = set(re.findall(r'"([EW]_[A-Z0-9_]+)', src))
    # The C++ validator runs the JSON Schema through the checker generated
    # from it (E_SCHEMA, E_PAYLOAD_SCHEMA), or reports W_SCHEMA_SKIPPED when
    # the check is disabled.
    assert cpp_codes == set(FINDING_CODES)
    generated = _repo_file("cpp", "include", "qvf", "kinds_generated.hpp")
    block = generated.split("kFindingCodes[] = {", 1)[1].split("};", 1)[0]
    rows = re.findall(r'\{"([EW]_[A-Z0-9_]+)", (true|false)\}', block)
    assert {c: "error" if e == "true" else "warning" for c, e in rows} == FINDING_CODES


def test_cpp_canonical_kinds_match_python_and_schema():
    src = _repo_file("cpp", "include", "qvf", "kinds_generated.hpp")
    block = src.split("kCanonicalKinds[] = {", 1)[1].split("};", 1)[0]
    cpp_kinds = set(re.findall(r'"([^"]+)"', block))
    schema = json.loads(_repo_file("spec", "qvf_manifest.schema.json"))
    defs = schema["$defs"]
    schema_kinds = set()
    for branch in defs["Section"]["oneOf"]:
        kind = defs[branch["$ref"].split("/")[-1]].get("properties", {}).get("kind", {})
        if "const" in kind:
            schema_kinds.add(kind["const"])
    assert cpp_kinds == set(qvf_reader._CANONICAL_KINDS) == schema_kinds
