#!/usr/bin/env python3
"""Generate the QVF negative (must-fail) conformance corpus.

Writes one small ``.qvf`` per validation rule next to this file, each paired
with a ``.expected.json`` naming the finding code(s) a conforming validator must
report (spec § 6, Appendix B). Every archive is a minimal *valid* archive with
exactly one defect, so it isolates one check.

The archives are built byte by byte with :mod:`zipfile` — deliberately **not**
with the reference writer, which refuses to produce them. Like the positive
corpus they are frozen test vectors: commit them, and regenerate only to extend
the set (zlib versions may change the compressed bytes, never the verdict).

Archives that exercise a resource cap stay tiny by naming a deliberately
lowered cap in their ``expected.json`` (``"limits": {...}``); the runner passes
those limits to the validator. With the default limits such an archive is
valid, which is the point: it proves the cap, not the archive, is what fails.

Run:  python generate_negative.py
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import os
import struct
import warnings
import zipfile
from typing import Any, Callable, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE

_PROGRAM = "qvf-conformance"
_VERSION = "0.1.0"
_STORED = zipfile.ZIP_STORED
_DEFLATED = zipfile.ZIP_DEFLATED

WATER = {
    "atoms": [
        {"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
        {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
        {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1},
    ],
    "pbc": [False, False, False],
    "lattice_vectors": None,
}


# --- building blocks ---------------------------------------------------------

def _json(obj: Any) -> bytes:
    return json.dumps(obj, indent=1).encode("utf-8")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_member(path: str, data: bytes) -> dict:
    return {"path": path, "format": "json", "sha256": _sha(data)}


def binary_member(path: str, data: bytes, dtype: str, shape: list) -> dict:
    return {"path": path, "format": "binary", "dtype": dtype, "shape": shape,
            "sha256": _sha(data)}


def opaque_member(path: str, data: bytes) -> dict:
    return {"path": path, "format": "binary", "sha256": _sha(data)}


class Archive:
    """A manifest plus ZIP entries, written in order (manifest STOREd last)."""

    def __init__(self, calculation: str) -> None:
        self.manifest: dict = {
            "qvf_version": 1,
            "source": {"program": _PROGRAM, "version": _VERSION,
                       "calculation": calculation},
            "sections": [],
        }
        self.entries: list[tuple[str, bytes, int]] = []
        self.manifest_bytes: Optional[bytes] = None   # override (broken JSON)
        self.include_manifest = True

    def add_entry(self, name: str, data: bytes, method: int = _DEFLATED) -> None:
        self.entries.append((name, data, method))

    def add_section(self, section: dict, files: dict[str, bytes]) -> dict:
        self.manifest["sections"].append(section)
        for name, data in files.items():
            self.add_entry(name, data)
        return section

    def add_structure(self, sid: str = "structure",
                      path: str = "structure/structure.json",
                      payload: Any = None) -> dict:
        data = _json(WATER if payload is None else payload)
        return self.add_section(
            {"id": sid, "kind": "structure",
             "members": {"structure": json_member(path, data)}},
            {path: data})

    def to_bytes(self) -> bytes:
        entries = list(self.entries)
        if self.include_manifest:
            mf = (self.manifest_bytes if self.manifest_bytes is not None
                  else _json(self.manifest))
            entries.append(("manifest.json", mf, _STORED))
        buf = io.BytesIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")      # zipfile warns on duplicates
            with zipfile.ZipFile(buf, "w") as zf:
                for name, data, method in entries:
                    info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                    info.compress_type = method
                    info.create_system = 3
                    info.external_attr = 0o644 << 16
                    zf.writestr(info, data, compress_type=method, compresslevel=9)
        return buf.getvalue()


def base_archive() -> Archive:
    """The valid archive every negative case starts from (one water structure)."""
    a = Archive("negative-base")
    a.add_structure()
    return a


def _vendor(a: Archive, sid: str, role: str, path: str, data: bytes,
            fmt: str = "json", **section_extra: Any) -> dict:
    if fmt == "json":
        member = json_member(path, data)
    else:
        member = binary_member(path, data, "uint8", [len(data)])
    return a.add_section({"id": sid, "kind": f"x_test.{sid}",
                          "members": {role: member}, **section_extra},
                         {path: data})


def _nested(depth: int) -> bytes:
    return b"[" * depth + b"]" * depth


def _patch_declared_size(zip_bytes: bytes, name: str, new_size: int) -> bytes:
    """Rewrite one entry's uncompressed size in its local header and central
    directory record (the CRC is left alone)."""
    out = bytearray(zip_bytes)
    raw_name = name.encode("utf-8")
    patched = 0
    for sig, size_off, name_len_off, name_off in (
            (b"PK\x03\x04", 22, 26, 30), (b"PK\x01\x02", 24, 28, 46)):
        pos = out.find(sig)
        while pos != -1:
            n = struct.unpack_from("<H", out, pos + name_len_off)[0]
            if bytes(out[pos + name_off:pos + name_off + n]) == raw_name:
                struct.pack_into("<I", out, pos + size_off, new_size)
                patched += 1
            pos = out.find(sig, pos + 4)
    assert patched == 2, patched
    return bytes(out)


# --- cases -------------------------------------------------------------------
# Each case returns (archive bytes, description, expected findings, limits).

Case = Callable[[], tuple]
CASES: dict[str, Case] = {}


def case(name: str) -> Callable[[Case], Case]:
    def register(fn: Case) -> Case:
        CASES[name] = fn
        return fn
    return register


def _err(*codes: str) -> list[dict]:
    return [{"code": c, "severity": "error"} for c in codes]


# Container rules (new in this revision of spec § 6).

@case("zip_duplicate_entry")
def _():
    a = base_archive()
    a.entries.insert(0, ("structure/structure.json", _json({"atoms": []}),
                         _DEFLATED))
    return (a.to_bytes(),
            "structure/structure.json occurs twice in the ZIP with different "
            "bytes; readers disagree on which entry wins (spec § 2.1).",
            _err("E_ZIP_DUPLICATE_ENTRY"), None)


def _bad_path_case(path: str, what: str):
    a = Archive("bad-path")
    a.add_structure(path=path)
    return (a.to_bytes(), f"Member path {path!r} {what} (spec § 2.1).",
            _err("E_MEMBER_PATH"), None)


@case("path_backslash")
def _():
    return _bad_path_case("structure\\structure.json", "uses a backslash")


@case("path_drive_letter")
def _():
    return _bad_path_case("C:/structure/structure.json",
                          "begins with a drive letter")


@case("path_absolute")
def _():
    return _bad_path_case("/structure/structure.json", "begins with '/'")


@case("path_dotdot")
def _():
    return _bad_path_case("structure/../../structure.json",
                          "contains a '..' segment")


@case("json_too_large")
def _():
    a = base_archive()
    _vendor(a, "blob", "data", "x_test/blob.json",
            _json({"values": list(range(600))}))
    return (a.to_bytes(),
            "A JSON member larger than the JSON size cap (lowered to 2 KiB "
            "here; the reference value is 256 MiB).",
            _err("E_JSON_TOO_LARGE"), {"max_json_bytes": 2048})


@case("json_too_deep")
def _():
    a = base_archive()
    _vendor(a, "deep", "data", "x_test/deep.json", _nested(100))
    return (a.to_bytes(),
            "A JSON member nesting 100 arrays, over the depth cap of 64.",
            _err("E_JSON_TOO_DEEP"), None)


@case("manifest_too_deep")
def _():
    a = base_archive()
    a.manifest["provenance"] = {"x_nested": json.loads(_nested(70))}
    return (a.to_bytes(),
            "manifest.json nests 72 levels deep (70 arrays inside provenance), "
            "over the depth cap of 64.",
            _err("E_JSON_TOO_DEEP"), None)


@case("member_too_large")
def _():
    a = base_archive()
    _vendor(a, "big", "data", "x_test/big.bin",
            bytes(range(256)) * 16, fmt="binary")
    return (a.to_bytes(),
            "A 4 KiB member over the per-member size cap (lowered to 2 KiB "
            "here; the reference value is 8 GiB).",
            _err("E_MEMBER_TOO_LARGE"), {"max_member_bytes": 2048})


@case("total_too_large")
def _():
    a = base_archive()
    for k in range(3):
        _vendor(a, f"part{k}", "data", f"x_test/part{k}.bin",
                bytes(range(256)) * 8, fmt="binary")
    return (a.to_bytes(),
            "Members that are each under the member cap but together expand "
            "past the total cap (lowered to 4 KiB here; the reference value "
            "is 64 GiB).",
            _err("E_TOTAL_TOO_LARGE"), {"max_total_bytes": 4096})


@case("compression_ratio")
def _():
    a = base_archive()
    _vendor(a, "zeros", "data", "x_test/zeros.bin", bytes(65536), fmt="binary")
    return (a.to_bytes(),
            "A 64 KiB run of zeros deflated to well under 1 KiB: a ratio over "
            "the explicit application cap (100 with a 4 KiB exemption floor here). "
            "The default ratio policy is disabled; this is valid QVF.",
            _err("E_COMPRESSION_RATIO"),
            {"max_ratio": 100, "ratio_exempt_bytes": 4096})


def _crc_mismatch(method):
    a = Archive("crc-mismatch")
    data = bytes(4096)
    _vendor(a, "data", "data", "x_test/data.bin", data, fmt="binary")
    del a.manifest["sections"][0]["members"]["data"]["sha256"]
    a.entries = [(name, body, method) for name, body, _ in a.entries]
    raw = bytearray(a.to_bytes())
    for signature, offset in ((b"PK\x03\x04", 14), (b"PK\x01\x02", 16)):
        at = raw.index(signature)
        crc = struct.unpack_from("<I", raw, at + offset)[0]
        struct.pack_into("<I", raw, at + offset, crc ^ 1)
    return (bytes(raw),
            "A member without SHA-256 whose valid payload has a wrong ZIP CRC-32.",
            _err("E_ZIP_UNREADABLE"), None)


@case("crc_mismatch_stored")
def _():
    return _crc_mismatch(_STORED)


@case("crc_mismatch_deflated")
def _():
    return _crc_mismatch(_DEFLATED)


@case("entry_size_understated")
def _():
    a = base_archive()
    data = bytes((i * 7) % 251 for i in range(8192))
    _vendor(a, "liar", "data", "x_test/liar.bin", data, fmt="binary")
    return (_patch_declared_size(a.to_bytes(), "x_test/liar.bin", 1024),
            "A member whose ZIP headers declare 1 KiB uncompressed while its "
            "data inflates to 8 KiB: a reader must not inflate past the "
            "declared size, so the entry is unreadable.",
            _err("E_ZIP_UNREADABLE"), None)


@case("orphan_member")
def _():
    a = base_archive()
    a.add_entry("stale/old_structure.json", _json(WATER))
    return (a.to_bytes(),
            "A ZIP entry no manifest member references (the by-product of an "
            "in-place update): a warning only, the archive stays valid.",
            [{"code": "W_ORPHAN_MEMBER", "severity": "warning"}], None)


# Member size and byte order (spec § 2.2, § 4.2). A duplicated entry name
# other than manifest.json is zip_duplicate_entry above.

@case("member_size_mismatch")
def _():
    a = base_archive()
    member = a.manifest["sections"][0]["members"]["structure"]
    data = _json(WATER)
    member["size"] = len(data) + 8
    return (a.to_bytes(),
            "The structure member declares a size 8 bytes larger than the "
            "uncompressed size its ZIP entry records; sha256 and payload are "
            "otherwise correct. A validator compares size with the central "
            "directory before reading the member (spec § 4.2).",
            _err("E_BINARY_SIZE"), None)


@case("byte_order_big")
def _():
    a = base_archive()
    a.manifest["byte_order"] = "big"
    return (a.to_bytes(),
            "A root byte_order of 'big': 'little' is the only byte order of "
            "QVF v1 (spec § 2.2), and the JSON Schema rejects any other value. "
            "Only a validator with a JSON Schema engine detects this; one "
            "without (it reports W_SCHEMA_SKIPPED) accepts the archive, while "
            "the reference readers refuse to open it.",
            _err("E_SCHEMA"), None)


@case("reserved_id_auto_open")
def _():
    a = Archive("reserved-id")
    a.add_structure(sid="auto_open")
    return (a.to_bytes(),
            "A section id 'auto_open', a reserved viewer_defaults key "
            "(spec § 3.1).",
            _err("E_RESERVED_ID"), None)


@case("reserved_id_bookmarks")
def _():
    a = Archive("reserved-id")
    a.add_structure(sid="bookmarks")
    return (a.to_bytes(),
            "A section id 'bookmarks', a reserved viewer_defaults key "
            "(spec § 3.1).",
            _err("E_RESERVED_ID"), None)


# Pre-existing checks, each with a stable code now.

@case("not_a_zip")
def _():
    return (b"This file is not a ZIP archive.\n",
            "Not a ZIP archive at all.", _err("E_ZIP_UNREADABLE"), None)


@case("manifest_missing")
def _():
    a = base_archive()
    a.include_manifest = False
    return (a.to_bytes(), "A ZIP archive without manifest.json.",
            _err("E_MANIFEST_MISSING"), None)


@case("manifest_json_invalid")
def _():
    a = base_archive()
    a.manifest_bytes = b'{"qvf_version": 1, "source": '
    return (a.to_bytes(), "manifest.json is truncated, invalid JSON.",
            _err("E_MANIFEST_JSON"), None)


@case("version_unsupported")
def _():
    a = base_archive()
    a.manifest["qvf_version"] = 3
    return (a.to_bytes(), "qvf_version 3, an unknown major version.",
            _err("E_VERSION"), None)


@case("source_empty")
def _():
    a = base_archive()
    a.manifest["source"]["program"] = ""
    return (a.to_bytes(), "source.program is empty.", _err("E_SOURCE"), None)


@case("section_malformed")
def _():
    a = base_archive()
    a.manifest["sections"].append(42)
    return (a.to_bytes(), "A sections entry that is not an object.",
            _err("E_SECTION_MALFORMED"), None)


@case("section_id_duplicate")
def _():
    a = base_archive()
    a.add_structure(path="structure/second.json")
    return (a.to_bytes(), "Two sections share the id 'structure'.",
            _err("E_SECTION_ID_DUPLICATE"), None)


@case("kind_unknown")
def _():
    a = base_archive()
    data = _json({"x": 1})
    a.add_section({"id": "bogus", "kind": "volume.bogus",
                   "members": {"data": json_member("bogus/data.json", data)}},
                  {"bogus/data.json": data})
    return (a.to_bytes(), "A non-vendor kind that is not canonical.",
            _err("E_KIND_UNKNOWN"), None)


@case("critical_undeclared")
def _():
    a = base_archive()
    _vendor(a, "ecp", "data", "x_test/ecp.json", _json({"core": 2}),
            critical=True)
    return (a.to_bytes(),
            "A critical vendor section whose namespace x_test is not declared "
            "in root extensions (spec § 7.4).",
            _err("E_CRITICAL_UNDECLARED"), None)


@case("member_malformed")
def _():
    a = base_archive()
    a.manifest["sections"][0]["members"]["bonds"] = {"format": "json"}
    return (a.to_bytes(), "A member spec without a path.",
            _err("E_MEMBER_MALFORMED"), None)


@case("member_missing")
def _():
    a = base_archive()
    a.manifest["sections"][0]["members"]["bonds"] = json_member(
        "structure/bonds.json", _json({"pairs": []}))
    return (a.to_bytes(), "A declared member path absent from the ZIP.",
            _err("E_MEMBER_MISSING"), None)


@case("sha256_mismatch")
def _():
    a = base_archive()
    a.manifest["sections"][0]["members"]["structure"]["sha256"] = _sha(b"other")
    return (a.to_bytes(), "A member whose bytes do not match its sha256.",
            _err("E_SHA256_MISMATCH"), None)


@case("json_invalid")
def _():
    a = Archive("json-invalid")
    data = b'{"atoms": ['
    a.add_section({"id": "structure", "kind": "structure",
                   "members": {"structure": json_member(
                       "structure/structure.json", data)}},
                  {"structure/structure.json": data})
    return (a.to_bytes(), "A JSON member that is not valid JSON.",
            _err("E_JSON_INVALID"), None)


@case("dtype_unknown")
def _():
    a = base_archive()
    data = bytes(6)
    a.add_section({"id": "atom_properties", "kind": "atom_properties",
                   "members": {"mulliken_charge": binary_member(
                       "atom_properties/q.bin", data, "float16", [3])}},
                  {"atom_properties/q.bin": data})
    return (a.to_bytes(), "A binary member with an unknown dtype.",
            _err("E_DTYPE"), None)


@case("binary_size")
def _():
    a = base_archive()
    data = struct.pack("<2d", -0.68, 0.34)
    a.add_section({"id": "atom_properties", "kind": "atom_properties",
                   "members": {"mulliken_charge": binary_member(
                       "atom_properties/q.bin", data, "float64", [3])}},
                  {"atom_properties/q.bin": data})
    return (a.to_bytes(),
            "A float64 [3] member that stores 16 bytes instead of 24.",
            _err("E_BINARY_SIZE"), None)


@case("run_record_utf8")
def _():
    a = base_archive()
    log = b"\xff\xfe not UTF-8\n"
    a.add_section({"id": "run_record", "kind": "run.record", "program": "t",
                   "members": {"log": opaque_member("run_record/log.txt", log)}},
                  {"run_record/log.txt": log})
    return (a.to_bytes(), "A run.record log that does not decode as UTF-8.",
            _err("E_UTF8"), None)


@case("files_index")
def _():
    a = base_archive()
    log = b"done\n"
    files = _json({"attachment.restart": {"filename": "restart.chk"}})
    a.add_section({"id": "run_record", "kind": "run.record", "program": "t",
                   "members": {
                       "log": opaque_member("run_record/log.txt", log),
                       "files": json_member("run_record/files.json", files)}},
                  {"run_record/log.txt": log, "run_record/files.json": files})
    return (a.to_bytes(),
            "A run.record files index naming a role with no member.",
            _err("E_FILES_INDEX"), None)


@case("ref_unresolved")
def _():
    a = base_archive()
    wp = _json({"waypoints": [{"frame": 0, "label": "reactant"}]})
    a.add_section({"id": "waypoints", "kind": "reaction.waypoints",
                   "trajectory_ref": "trajectory",
                   "members": {"waypoints": json_member("reaction/wp.json", wp)}},
                  {"reaction/wp.json": wp})
    return (a.to_bytes(),
            "reaction.waypoints.trajectory_ref names a section that does not "
            "exist.", _err("E_REF_UNRESOLVED"), None)


@case("structure_periodicity")
def _():
    payload = dict(WATER)
    del payload["pbc"]
    payload["lattice_vectors"] = [[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]]
    a = Archive("periodicity")
    a.add_structure(payload=payload)
    return (a.to_bytes(),
            "A structure payload with lattice_vectors but no pbc (spec § 5.1).",
            _err("E_STRUCTURE_PERIODICITY"), None)


@case("structure_occupancy_range")
def _():
    payload = json.loads(json.dumps(WATER))
    payload["atoms"][0].update({"alt_loc": "A", "occupancy": 1.5})
    a = Archive("occupancy")
    a.add_structure(payload=payload)
    return (a.to_bytes(),
            "A structure atom with alt_loc A and occupancy 1.5; occupancy is a "
            "fraction from 0 to 1 (spec § 5.1, Biomolecular atom identity).",
            _err("E_PAYLOAD_SCHEMA"), None)



@case("structure_helix_class_on_sheet")
def _():
    a = Archive("secondary structure")
    section = a.add_structure()
    section["secondary_structure"] = [
        {"type": "sheet", "chain": "A", "start_seq": 1, "end_seq": 3, "helix_class": 1}]
    return (a.to_bytes(),
            "A sheet secondary-structure range with a helix_class, which applies "
            "only to helix ranges (spec § 5.1, Secondary-structure ranges). Only a "
            "validator with a JSON Schema engine detects this; one without (it "
            "reports W_SCHEMA_SKIPPED) accepts the archive.",
            _err("E_SCHEMA"), None)

# Member payloads that break their $defs definition (spec § 5.4, § 5.7). The
# expected E_PAYLOAD_SCHEMA is a schema-engine code (Appendix B): a validator
# without a JSON Schema engine reports W_SCHEMA_SKIPPED instead and need not
# reject these archives, so the runners require the code only from a validator
# that ran the schema check.

def _payload_case(kind: str, sid: str, role: str, payload: Any) -> Archive:
    a = base_archive()
    path = f"{sid}/{role}.json"
    data = _json(payload)
    a.add_section({"id": sid, "kind": kind,
                   "members": {role: json_member(path, data)}}, {path: data})
    return a


@case("spectrum_length_mismatch")
def _():
    a = _payload_case("spectra.ir", "spectra_ir", "spectrum",
                      {"frequencies": [1595.0, 3657.0, 3756.0],
                       "intensities": [67.0, 4.0]})
    return (a.to_bytes(),
            "A spectra.ir payload with three frequencies but two intensities; "
            "the arrays must have one length (spec § 5.4).",
            _err("E_PAYLOAD_SCHEMA"), None)


@case("qtaim_point_type")
def _():
    a = _payload_case("topology.qtaim", "qtaim", "critical_points",
                      {"points": [{"type": "bond",
                                   "position": [0.0, 0.5788, -0.3149],
                                   "rho": 0.3653, "atom_pair": [0, 1]}]})
    return (a.to_bytes(),
            "A topology.qtaim critical point of type 'bond'; the type is one "
            "of bcp, rcp, ccp, ncp (spec § 5.7).",
            _err("E_PAYLOAD_SCHEMA"), None)


@case("scf_history_unknown_key")
def _():
    a = _payload_case("scf_history", "scf_history", "iterations",
                      {"iterations": [{"iter": 1, "energy": -75.9}]})
    return (a.to_bytes(),
            "An scf_history iteration that names its energy 'energy' instead "
            "of 'energy_eh': the strict payload rejects the unknown key and "
            "misses the required one (spec § 5.7).",
            _err("E_PAYLOAD_SCHEMA"), None)


# Payload definitions (spec § 5, § 6). Only a validator with a JSON Schema
# engine reports E_PAYLOAD_SCHEMA; one without (it reports W_SCHEMA_SKIPPED)
# must find nothing else wrong with these archives.

def _volume_density(a: Archive, grid: Any) -> dict:
    g = _json(grid)
    data = struct.pack("<64f", *([0.125] * 64))
    return a.add_section(
        {"id": "rho", "kind": "volume.density",
         "members": {"grid": json_member("volumes/rho_grid.json", g),
                     "data": binary_member("volumes/rho.dat", data, "float32",
                                           [4, 4, 4])}},
        {"volumes/rho_grid.json": g, "volumes/rho.dat": data})


def _wavefunction(a: Archive, basis: Any, mo_metadata: Any,
                  unrestricted: bool = False) -> dict:
    b, m = _json(basis), _json(mo_metadata)
    c = struct.pack("<4d", 1.0, 0.0, 0.0, 1.0)
    members = {"basis": json_member("wf/basis.json", b),
               "mo_metadata": json_member("wf/mo_metadata.json", m)}
    files = {"wf/basis.json": b, "wf/mo_metadata.json": m}
    roles = (("mo_coefficients_alpha", "mo_coefficients_beta") if unrestricted
             else ("mo_coefficients",))
    for role in roles:
        members[role] = binary_member(f"wf/{role}.dat", c, "float64", [2, 2])
        files[f"wf/{role}.dat"] = c
    return a.add_section({"id": "wf", "kind": "wavefunction.gto",
                          "members": members}, files)


_H2_BASIS = {"structure_ref": "structure", "pure": True, "n_ao": 2,
             "shells": [{"center": 1, "l": 0, "exponents": [1.24],
                         "coefficients": [1.0]},
                        {"center": 2, "l": 0, "exponents": [1.24],
                         "coefficients": [1.0]}]}


@case("payload_grid_shape")
def _():
    a = base_archive()
    _volume_density(a, {"origin": [-1.0, -1.0, -1.0],
                        "voxel_vectors": [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0],
                                          [0.0, 0.0, 0.5]],
                        "shape": [4, 4]})
    return (a.to_bytes(),
            "A volume.density grid whose shape has two entries instead of "
            "three ($defs/VolumeGridPayload, spec § 5.2).",
            _err("E_PAYLOAD_SCHEMA"), None)


@case("payload_basis_exponent")
def _():
    a = base_archive()
    basis = json.loads(json.dumps(_H2_BASIS))
    basis["shells"][1]["exponents"] = [-1.24]
    _wavefunction(a, basis, {"n_mo": 2, "n_ao": 2, "spin": "restricted",
                             "orbital_kind": "canonical",
                             "energies": [-0.58, 0.67],
                             "occupations": [2.0, 0.0]})
    return (a.to_bytes(),
            "A wavefunction.gto basis shell with a negative primitive exponent "
            "($defs/BasisPayload, spec § 5.3).",
            _err("E_PAYLOAD_SCHEMA"), None)


@case("payload_mo_metadata_shape")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, {"n_ao": 2, "spin": "unrestricted",
                                 "orbital_kind": "canonical",
                                 "energies": [-0.58, 0.67],
                                 "occupations": [1.0, 0.0]},
                  unrestricted=True)
    return (a.to_bytes(),
            "An unrestricted mo_metadata that carries the restricted shape "
            "(top-level energies and occupations) instead of alpha and beta "
            "objects ($defs/MoMetadataPayload, spec § 5.3).",
            _err("E_PAYLOAD_SCHEMA"), None)


# --- driver ------------------------------------------------------------------

@case("bands_spin_mismatch")
def _():
    a = base_archive()
    kpath = _json({"n_spin": 1, "n_kpoints": 2, "n_bands": 2})
    eig = struct.pack("<4d", -1.0, 1.0, -0.5, 1.5)          # [1, 2, 2]
    proj = struct.pack("<8d", *([0.5] * 8))                 # [2, 2, 2, 1]
    a.add_section({"id": "bands", "kind": "bands", "members": {
        "kpath": json_member("bands/kpath.json", kpath),
        "eigenvalues": binary_member("bands/eigenvalues.bin", eig, "float64",
                                     [1, 2, 2]),
        "projections": binary_member("bands/projections.bin", proj, "float64",
                                     [2, 2, 2, 1])}},
        {"bands/kpath.json": kpath, "bands/eigenvalues.bin": eig,
         "bands/projections.bin": proj})
    return (a.to_bytes(),
            "A bands section whose rank-4 projections carry n_spin 2 while its "
            "eigenvalues carry n_spin 1 (spec § 5.6).",
            _err("E_MEMBER_MALFORMED"), None)


@case("eos_model_unknown")
def _():
    a = base_archive()
    vols = struct.pack("<3d", 90.0, 100.0, 110.0)
    energies = struct.pack("<3d", -5.1, -5.2, -5.15)
    fit = _json({"model": "spline", "V0": 100.0, "E0": -5.2, "B0": 74.0,
                 "B0_prime": 4.0})
    a.add_section({"id": "eos", "kind": "equation_of_state", "members": {
        "volumes": binary_member("eos/volumes.bin", vols, "float64", [3]),
        "energies": binary_member("eos/energies.bin", energies, "float64", [3]),
        "fit": json_member("eos/fit.json", fit)}},
        {"eos/volumes.bin": vols, "eos/energies.bin": energies,
         "eos/fit.json": fit})
    return (a.to_bytes(),
            "An equation_of_state fit payload whose model is not one of the "
            "models of spec § 5.6. Only a validator with a JSON Schema engine "
            "detects this; one without (it reports W_SCHEMA_SKIPPED) accepts "
            "the archive.",
            _err("E_PAYLOAD_SCHEMA"), None)


# Declared units (spec § 4.3). An unknown unit name is the warning
# W_UNIT_UNKNOWN in every validator; the schema also types it, so a validator
# with a JSON Schema engine adds E_PAYLOAD_SCHEMA and rejects the archive.

_GRID = {"origin": [-1.0, -1.0, -1.0],
         "voxel_vectors": [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]],
         "shape": [4, 4, 4]}


@case("unit_unknown")
def _():
    a = base_archive()
    _volume_density(a, {**_GRID, "units": {"origin": "bohr",
                                           "data": "electrons_per_bohr3"}})
    return (a.to_bytes(),
            "A volume.density grid whose units object declares "
            "'electrons_per_bohr3', which is not a unit name of spec § 4.3. "
            "Every validator warns (W_UNIT_UNKNOWN); one with a JSON Schema "
            "engine also rejects the payload (E_PAYLOAD_SCHEMA), while one "
            "without accepts the archive with the warning.",
            [{"code": "W_UNIT_UNKNOWN", "severity": "warning"}]
            + _err("E_PAYLOAD_SCHEMA"), None)


@case("unit_contradicts_fixed")
def _():
    a = base_archive()
    _volume_density(a, {**_GRID, "units": {"origin": "angstrom"}})
    return (a.to_bytes(),
            "A volume.density grid that declares units.origin 'angstrom', a "
            "unit name that contradicts the bohr spec § 4.3 fixes for the grid "
            "origin. Only a validator with a JSON Schema engine detects this; "
            "one without accepts the archive.",
            _err("E_PAYLOAD_SCHEMA"), None)

# Energy derivatives and the root energies block (spec § 3.6, § 5.10, § 5.11).

@case("forces_atom_count")
def _():
    a = base_archive()
    grad = struct.pack("<6d", 0.0, 0.0, 0.0025, 0.0, 0.0016, -0.0013)  # [2, 3]
    meta = _json({"structure_ref": "structure", "convention": "gradient"})
    a.add_section({"id": "forces", "kind": "forces", "members": {
        "gradient": binary_member("forces/gradient.bin", grad, "float64", [2, 3]),
        "meta": json_member("forces/meta.json", meta)}},
        {"forces/gradient.bin": grad, "forces/meta.json": meta})
    return (a.to_bytes(),
            "A forces gradient with 2 rows for a 3-atom structure named by "
            "meta.structure_ref; gradient must be [n_atoms, 3] (spec § 5.10).",
            _err("E_MEMBER_MALFORMED"), None)


@case("hessian_not_square")
def _():
    a = base_archive()
    matrix = struct.pack("<54d", *([0.0] * 54))                 # [9, 6]
    a.add_section({"id": "hessian", "kind": "hessian", "members": {
        "matrix": binary_member("hessian/matrix.bin", matrix, "float64", [9, 6])}},
        {"hessian/matrix.bin": matrix})
    return (a.to_bytes(),
            "A hessian matrix of shape [9, 6]; it must be the square [3N, 3N] "
            "(spec § 5.11).",
            _err("E_MEMBER_MALFORMED"), None)


@case("energies_inconsistent")
def _():
    a = base_archive()
    a.manifest["provenance"] = {
        "method": "rhf", "scf_energy": {"value": -74.9659011, "units": "hartree"}}
    a.manifest["energies"] = {
        "total": {"value": -74.9600000, "units": "hartree"}, "method": "rhf"}
    return (a.to_bytes(),
            "energies.total, stated for the method provenance.method names, "
            "differs from provenance.scf_energy by 5.9e-3 hartree; the two "
            "must agree within 1e-8 hartree (spec § 3.6).",
            _err("E_ROOT_INCONSISTENT"), None)

# Volumes (spec § 5.2): periodic grids against the cell of their structure, the
# data shape of the 4-D kinds, and a series axis against its slabs.

_CUBIC_CELL = {"atoms": [{"symbol": "Na", "position": [0.0, 0.0, 0.0],
                          "atomic_number": 11}],
               "pbc": [True, True, True],
               "lattice_vectors": [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0],
                                   [0.0, 0.0, 4.0]]}
# One lattice step of the 4 Å cell over four points, in bohr (4 x 1.8897261246 / 4).
_CHGCAR_STEP = 1.8897261246


def _periodic_grid(voxel_vectors: list, includes_endpoint: bool) -> Archive:
    a = Archive("negative-periodic-grid")
    a.add_structure(payload=_CUBIC_CELL)
    _volume_density(a, {"origin": [0.0, 0.0, 0.0], "voxel_vectors": voxel_vectors,
                        "shape": [4, 4, 4], "periodic": [True, True, True],
                        "includes_endpoint": includes_endpoint})
    return a


def _volume_4d(a: Archive, kind: str, shape: list,
               axis: Optional[dict] = None) -> dict:
    n = 1
    for d in shape:
        n *= d
    data = struct.pack(f"<{n}f", *([0.25] * n))
    spatial = shape[1:] if kind == "volume.series" else shape[:3]
    g = _json({"origin": [-1.0, -1.0, -1.0],
               "voxel_vectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                 [0.0, 0.0, 1.0]],
               "shape": spatial})
    members = {"grid": json_member("volumes/v_grid.json", g),
               "data": binary_member("volumes/v.dat", data, "float32", shape)}
    files = {"volumes/v_grid.json": g, "volumes/v.dat": data}
    if axis is not None:
        ax = _json(axis)
        members["axis"] = json_member("volumes/v_axis.json", ax)
        files["volumes/v_axis.json"] = ax
    return a.add_section({"id": "v", "kind": kind, "members": members}, files)


@case("volume_periodic_span")
def _():
    a = _periodic_grid([[_CHGCAR_STEP, 0.0, 0.0], [0.0, _CHGCAR_STEP, 0.0],
                        [0.0, 0.0, 1.9]], False)
    return (a.to_bytes(),
            "A volume.density grid periodic on all three axes of a 4 Å cubic "
            "cell, CHGCAR style (includes_endpoint false), whose third voxel "
            "vector times shape[2] is 7.6 bohr instead of the 7.5589 bohr of "
            "lattice row 2: a periodic grid must span its cell within 1e-6 "
            "relative (spec § 5.2).",
            _err("E_MEMBER_MALFORMED"), None)


@case("volume_endpoint_count")
def _():
    step = [[_CHGCAR_STEP, 0.0, 0.0], [0.0, _CHGCAR_STEP, 0.0],
            [0.0, 0.0, _CHGCAR_STEP]]
    a = _periodic_grid(step, True)
    return (a.to_bytes(),
            "A periodic volume.density grid with the CHGCAR-style point count "
            "(voxel_vectors[a] x shape[a] is the 4 Å cell) that declares "
            "includes_endpoint true, under which voxel_vectors[a] x "
            "(shape[a] - 1) must span the cell: the silent cube/CHGCAR mix-up "
            "the endpoint rule exists to catch (spec § 5.2).",
            _err("E_MEMBER_MALFORMED"), None)


@case("volume_series_axis_length")
def _():
    a = base_archive()
    _volume_4d(a, "volume.series", [3, 2, 2, 2],
               {"kind": "time", "values": [0.0, 10.0], "units": "fs"})
    return (a.to_bytes(),
            "A volume.series whose data holds 3 slabs but whose axis has 2 "
            "values; the axis needs one value per slab (spec § 5.2).",
            _err("E_MEMBER_MALFORMED"), None)


@case("volume_vector_trailing_axis")
def _():
    a = base_archive()
    _volume_4d(a, "volume.vector", [2, 2, 2, 2])
    return (a.to_bytes(),
            "A volume.vector whose data is [2, 2, 2, 2]; a vector field is "
            "[nx, ny, nz, 3] (spec § 5.2). A validator with a JSON Schema "
            "engine also reports E_SCHEMA.",
            _err("E_MEMBER_MALFORMED", "E_SCHEMA"), None)

# atom_properties roles (spec § 5.7). The role list, the vector shape and the
# custom-role rule are checked by every validator except the vector shape and
# the closed role list, which only a JSON Schema engine can judge.

def _atom_properties(a: Archive, roles: dict, meta: Any = None) -> dict:
    """An atom_properties section for the two atoms of a structure; `roles`
    maps a role to (values, shape), `meta` (when not None) is its meta member."""
    members: dict = {}
    files: dict[str, bytes] = {}
    for role, (values, shape) in roles.items():
        data = struct.pack(f"<{len(values)}d", *values)
        members[role] = binary_member(f"atom_properties/{role}.bin", data,
                                      "float64", shape)
        files[f"atom_properties/{role}.bin"] = data
    if meta is not None:
        data = _json(meta)
        members["meta"] = json_member("atom_properties/meta.json", data)
        files["atom_properties/meta.json"] = data
    return a.add_section({"id": "atom_properties", "kind": "atom_properties",
                          "members": members}, files)


_CUSTOM_ENTRY = {"units": "e", "kind": "scalar",
                 "description": "Illustrative net flux."}


def _two_atoms() -> Archive:
    a = Archive("negative-atom-properties")
    a.add_structure(payload={"atoms": WATER["atoms"][:2], "pbc": [False] * 3,
                             "lattice_vectors": None})
    return a


@case("atom_properties_custom_no_meta")
def _():
    a = _two_atoms()
    _atom_properties(a, {"custom.flux": ([0.5, -0.5], [2])})
    return (a.to_bytes(),
            "An atom_properties section with a custom.flux role and no meta "
            "member: a custom role is allowed only with a meta member that "
            "describes it (spec § 5.7). A validator with a JSON Schema engine "
            "also reports E_SCHEMA.",
            _err("E_MEMBER_MALFORMED", "E_SCHEMA"), None)


@case("atom_properties_custom_meta_entry_missing")
def _():
    a = _two_atoms()
    _atom_properties(a, {"custom.flux": ([0.5, -0.5], [2]),
                         "custom.other": ([1.0, 2.0], [2])},
                     {"custom.flux": _CUSTOM_ENTRY})
    return (a.to_bytes(),
            "An atom_properties section whose meta member describes "
            "custom.flux but not custom.other: every custom role needs its own "
            "entry (spec § 5.7).",
            _err("E_MEMBER_MALFORMED"), None)


@case("atom_properties_vector_shape")
def _():
    a = _two_atoms()
    _atom_properties(a, {"magnetic_moment": ([0.9, 0.1], [2])})
    return (a.to_bytes(),
            "An atom_properties magnetic_moment of shape [n_atoms] instead of "
            "[n_atoms, 3]; a vector role is float64 [n_atoms, 3] "
            "($defs/AtomVectorBinary, spec § 5.7). Only a validator with a "
            "JSON Schema engine detects this; one without (it reports "
            "W_SCHEMA_SKIPPED) accepts the archive.",
            _err("E_SCHEMA"), None)


@case("atom_properties_role_unknown")
def _():
    a = _two_atoms()
    _atom_properties(a, {"bogus_charge": ([-0.3, 0.3], [2])})
    return (a.to_bytes(),
            "An atom_properties role that is neither in the role list nor a "
            "custom.<name> role (spec § 5.7). Only a validator with a JSON "
            "Schema engine detects this; one without (it reports "
            "W_SCHEMA_SKIPPED) accepts the archive.",
            _err("E_SCHEMA"), None)


# Section references (spec § 4.1, § 5.5) and the root metadata block (§ 3.5).

def _bonds(a: Archive, **section_extra: Any) -> dict:
    data = _json({"pairs": [{"i": 0, "j": 1, "order": 1.0}]})
    return a.add_section({"id": "bonds", "kind": "bonds",
                          "members": {"bonds": json_member("bonds/bonds.json",
                                                           data)},
                          **section_extra},
                         {"bonds/bonds.json": data})


@case("structure_ref_wrong_kind")
def _():
    a = base_archive()
    log = b"SCF converged\n"
    a.add_section({"id": "run", "kind": "run.record", "program": "t",
                   "members": {"log": opaque_member("run/log.txt", log)}},
                  {"run/log.txt": log})
    _bonds(a, structure_ref="run")
    return (a.to_bytes(),
            "A section whose structure_ref names an existing section of kind "
            "run.record instead of a structure section (spec § 4.1).",
            _err("E_REF_UNRESOLVED"), None)


@case("run_ref_unresolved")
def _():
    a = base_archive()
    _bonds(a, structure_ref="structure", run_ref="run")
    return (a.to_bytes(),
            "A section whose run_ref names a section that does not exist "
            "(spec § 4.1).",
            _err("E_REF_UNRESOLVED"), None)


@case("trajectory_ref_wrong_kind")
def _():
    a = base_archive()
    wp = _json({"waypoints": [{"frame_index": 0, "label": "R",
                               "kind": "reactant"}]})
    a.add_section({"id": "waypoints", "kind": "reaction.waypoints",
                   "trajectory_ref": "structure",
                   "members": {"waypoints": json_member("reaction/wp.json", wp)}},
                  {"reaction/wp.json": wp})
    return (a.to_bytes(),
            "reaction.waypoints.trajectory_ref names an existing section of "
            "kind structure instead of a trajectory section (spec § 5.5).",
            _err("E_REF_UNRESOLVED"), None)


@case("metadata_uuid_invalid")
def _():
    a = base_archive()
    a.manifest["metadata"] = {"title": "uuid without hyphens",
                              "uuid": "3f2504e04f8941d39a0c0305e82c3301"}
    return (a.to_bytes(),
            "A root metadata uuid of 32 hex digits without the hyphens of the "
            "RFC 4122 8-4-4-4-12 form (spec § 3.5). Only a validator with a "
            "JSON Schema engine detects this; one without (it reports "
            "W_SCHEMA_SKIPPED) accepts the archive.",
            _err("E_SCHEMA"), None)


# Wavefunction extensions (spec § 5.3, Appendix A.1).

def _primitive_norm(alpha: float, l: int) -> float:
    """N_i of spec Appendix A.1: the norm of the axial component (l, 0, 0)."""
    df = 1.0
    for k in range(1, 2 * l, 2):
        df *= k
    return ((2.0 * alpha / math.pi) ** 0.75 * (4.0 * alpha) ** (l / 2.0)
            / math.sqrt(df))


def _one_centre_overlap(aos: list) -> list:
    """Overlap matrix of single N_i-normalized Cartesian primitives on one
    centre; ``aos`` is a list of ``(alpha, l, (l_x, l_y, l_z))``. Per axis,
    int x^n exp(-p x^2) dx = (n - 1)!! / (2p)^(n/2) sqrt(pi / p) for even n."""
    rows = []
    for a, la, pa in aos:
        row = []
        for b, lb, pb in aos:
            p = a + b
            value = _primitive_norm(a, la) * _primitive_norm(b, lb)
            for k in range(3):
                n = pa[k] + pb[k]
                if n % 2:
                    value = 0.0
                    break
                df = 1.0
                for j in range(1, n, 2):
                    df *= j
                value *= df / (2.0 * p) ** (n // 2) * math.sqrt(math.pi / p)
            row.append(value)
        rows.append(row)
    return rows


@case("wavefunction_normalization")
def _():
    a = base_archive()
    # An s and a Cartesian d shell on the oxygen, one primitive each with
    # coefficient 1, so every axial self-overlap is 1 and d_xy's is 1/3.
    basis = _json({"structure_ref": "structure", "pure": False, "n_ao": 7,
                   "shells": [{"center": 0, "l": 0, "exponents": [1.2],
                               "coefficients": [1.0]},
                              {"center": 0, "l": 2, "exponents": [0.8],
                               "coefficients": [1.0]}]})
    meta = _json({"n_mo": 7, "n_ao": 7, "spin": "restricted",
                  "orbital_kind": "canonical", "energies": [0.0] * 7,
                  "occupations": [2.0] + [0.0] * 6})
    powers = [(0, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0),
              (0, 1, 1), (0, 0, 2)]
    overlap = _one_centre_overlap(
        [(1.2, 0, powers[0])] + [(0.8, 2, p) for p in powers[1:]])
    overlap[2][2] *= 3.0          # d_xy: 1/3 becomes 1, the per-component bug
    s = struct.pack("<49d", *[x for row in overlap for x in row])
    c = struct.pack("<49d", *[1.0 if i == j else 0.0
                              for i in range(7) for j in range(7)])
    a.add_section({"id": "wf", "kind": "wavefunction.gto", "members": {
        "basis": json_member("wf/basis.json", basis),
        "mo_metadata": json_member("wf/mo_metadata.json", meta),
        "mo_coefficients": binary_member("wf/mo_coefficients.dat", c,
                                         "float64", [7, 7]),
        "overlap": binary_member("wf/overlap.dat", s, "float64", [7, 7])}},
        {"wf/basis.json": basis, "wf/mo_metadata.json": meta,
         "wf/mo_coefficients.dat": c, "wf/overlap.dat": s})
    return (a.to_bytes(),
            "A wavefunction.gto overlap member whose d_xy diagonal is 1, three "
            "times the 1/3 that Appendix A.1 gives for the stored Cartesian d "
            "shell: the overlap of a producer that normalized each Cartesian "
            "component on its own (spec § 5.3, Appendix A.1).",
            _err("E_NORMALIZATION"), None)


@case("wavefunction_n_electrons")
def _():
    a = base_archive()
    basis = json.loads(json.dumps(_H2_BASIS))
    basis["n_electrons"] = 3
    _wavefunction(a, basis, {"n_mo": 2, "n_ao": 2, "spin": "restricted",
                             "orbital_kind": "canonical",
                             "energies": [-0.58, 0.67],
                             "occupations": [2.0, 0.0]})
    return (a.to_bytes(),
            "A wavefunction.gto basis whose n_electrons is 3 while the "
            "mo_metadata occupations sum to 2; the two must agree (spec § 5.3).",
            _err("E_MEMBER_MALFORMED"), None)


@case("wavefunction_complex_rank2")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, {"n_mo": 2, "n_ao": 2, "spin": "restricted",
                                 "orbital_kind": "canonical", "complex": True,
                                 "energies": [-0.58, 0.67],
                                 "occupations": [2.0, 0.0]})
    return (a.to_bytes(),
            "A wavefunction.gto whose mo_metadata says complex: true while its "
            "mo_coefficients member is the real rank-2 [2, 2], without the "
            "trailing [re, im] axis (spec § 5.3).",
            _err("E_MEMBER_MALFORMED"), None)


# Excited states (spec § 5.12). The state table's index order, its references
# and its orbital rows are checked by every validator, with or without a JSON
# Schema engine.

def _excited_states(a: Archive, states: list, **section_extra: Any) -> dict:
    data = _json({"states": states, "meta": {"method": "tddft",
                                              "n_states": len(states)}})
    return a.add_section({"id": "excited_states", "kind": "excited_states",
                          **section_extra,
                          "members": {"states": json_member(
                              "excited_states/states.json", data)}},
                         {"excited_states/states.json": data})


def _state(index: int, energy_ev: float, **fields: Any) -> dict:
    return {"index": index,
            "excitation_energy": {"value": energy_ev, "units": "ev"}, **fields}


_H2_MO = {"n_mo": 2, "n_ao": 2, "spin": "restricted", "orbital_kind": "canonical",
          "energies": [-0.58, 0.67], "occupations": [2.0, 0.0]}


@case("excited_states_index_duplicate")
def _():
    a = base_archive()
    _excited_states(a, [_state(1, 7.4), _state(2, 9.2), _state(2, 9.85)])
    return (a.to_bytes(),
            "An excited_states payload whose third state repeats index 2; "
            "state indices are unique and ascending (spec § 5.12).",
            _err("E_MEMBER_MALFORMED"), None)


@case("excited_states_nto_ref_wrong_kind")
def _():
    a = base_archive()
    _excited_states(a, [_state(1, 7.4, nto_ref="structure")])
    return (a.to_bytes(),
            "An excited state whose nto_ref names an existing section of kind "
            "structure instead of a wavefunction.gto section holding its "
            "natural transition orbitals (spec § 5.12).",
            _err("E_REF_UNRESOLVED"), None)


@case("excited_states_to_mo_out_of_range")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, _H2_MO)
    _excited_states(a, [_state(1, 7.4, contributions=[
        {"from_mo": 0, "to_mo": 2, "weight": 0.98}])], wavefunction_ref="wf")
    return (a.to_bytes(),
            "An excited state whose contribution has to_mo 2, but the "
            "wavefunction.gto that wavefunction_ref names has 2 orbital rows "
            "(0 and 1) in mo_coefficients (spec § 5.12).",
            _err("E_MEMBER_MALFORMED"), None)



# Volumes that name their wavefunction and orbital (spec § 5.2, Volumes of a
# wavefunction).

_H2_MO_UNRESTRICTED = {"n_ao": 2, "spin": "unrestricted", "orbital_kind": "canonical",
                       "alpha": {"energies": [-0.61, 0.31], "occupations": [1.0, 0.0]},
                       "beta": {"energies": [-0.42, 0.55], "occupations": [1.0, 0.0]}}


def _linked_volume(a: Archive, kind: str, sid: str, **fields: Any) -> dict:
    grid = _json({"origin": [-1.0, -1.0, -1.0],
                  "voxel_vectors": [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0],
                                    [0.0, 0.0, 0.5]],
                  "shape": [4, 4, 4]})
    data = struct.pack("<64f", *([0.125] * 64))
    return a.add_section(
        {"id": sid, "kind": kind,
         "members": {"grid": json_member(f"volumes/{sid}_grid.json", grid),
                     "data": binary_member(f"volumes/{sid}.dat", data, "float32",
                                           [4, 4, 4])},
         **fields},
        {f"volumes/{sid}_grid.json": grid, f"volumes/{sid}.dat": data})


@case("volume_wavefunction_ref_wrong_kind")
def _():
    a = base_archive()
    _linked_volume(a, "volume.density", "rho", wavefunction_ref="structure")
    return (a.to_bytes(),
            "A volume.density whose wavefunction_ref names an existing section of "
            "kind structure instead of a wavefunction.gto section (spec § 5.2).",
            _err("E_REF_UNRESOLVED"), None)


@case("volume_orbital_mo_index_out_of_range")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, _H2_MO)
    _linked_volume(a, "volume.orbital", "mo_2", wavefunction_ref="wf", mo_index=2)
    return (a.to_bytes(),
            "A volume.orbital with mo_index 2, but the wavefunction.gto that "
            "wavefunction_ref names has 2 rows (0 and 1) in mo_coefficients "
            "(spec § 5.2).",
            _err("E_REF_UNRESOLVED"), None)


@case("volume_orbital_spin_missing")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, _H2_MO_UNRESTRICTED, unrestricted=True)
    _linked_volume(a, "volume.orbital", "mo_0", wavefunction_ref="wf", mo_index=0)
    return (a.to_bytes(),
            "A volume.orbital with mo_index 0 and no spin, but the wavefunction.gto "
            "that wavefunction_ref names is unrestricted: row 0 of "
            "mo_coefficients_alpha or of mo_coefficients_beta (spec § 5.2).",
            _err("E_REF_UNRESOLVED"), None)


@case("volume_orbital_spin_restricted")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, _H2_MO)
    _linked_volume(a, "volume.orbital", "mo_0", wavefunction_ref="wf", mo_index=0,
                   spin="alpha")
    return (a.to_bytes(),
            "A volume.orbital with spin 'alpha', but the wavefunction.gto that "
            "wavefunction_ref names carries one set of orbitals in "
            "mo_coefficients (spec § 5.2).",
            _err("E_REF_UNRESOLVED"), None)


@case("volume_orbital_mo_index_without_ref")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, _H2_MO)
    _linked_volume(a, "volume.orbital", "mo_0", mo_index=0)
    return (a.to_bytes(),
            "A volume.orbital with mo_index but no wavefunction_ref, which "
            "mo_index requires (spec § 5.2). Only a validator with a JSON Schema "
            "engine detects this; one without (it reports W_SCHEMA_SKIPPED) "
            "accepts the archive.",
            _err("E_SCHEMA"), None)

# Typed provenance fields (spec § 3.2) and the grown job.spec payload (§ 5.9).

@case("job_spec_structure_ref_wrong_kind")
def _():
    a = base_archive()
    _bonds(a)
    spec = _json({"job_type": "molecular", "method": "rhf",
                  "structure_ref": "bonds"})
    a.add_section({"id": "job_spec", "kind": "job.spec",
                   "members": {"spec": json_member("job_spec/spec.json", spec)}},
                  {"job_spec/spec.json": spec})
    return (a.to_bytes(),
            "A job.spec payload whose structure_ref names an existing section "
            "of kind bonds instead of a structure section (spec § 5.9).",
            _err("E_REF_UNRESOLVED"), None)


@case("provenance_git_commit_not_hex")
def _():
    a = base_archive()
    a.manifest["provenance"] = {"method": "rhf", "git_commit": "9f2b6c1d4e7g"}
    return (a.to_bytes(),
            "A provenance git_commit of 12 characters, one of which ('g') is "
            "not a hexadecimal digit (spec § 3.2). Only a validator with a "
            "JSON Schema engine detects this; one without (it reports "
            "W_SCHEMA_SKIPPED) accepts the archive.",
            _err("E_SCHEMA"), None)


@case("provenance_started_utc_not_utc")
def _():
    a = base_archive()
    a.manifest["provenance"] = {"method": "rhf",
                                "started_utc": "2026-10-07T14:00:00+02:00"}
    return (a.to_bytes(),
            "A provenance started_utc with a +02:00 offset instead of the "
            "trailing Z of the UTC form YYYY-MM-DDThh:mm:ssZ (spec § 3.2). "
            "Only a validator with a JSON Schema engine detects this; one "
            "without (it reports W_SCHEMA_SKIPPED) accepts the archive.",
            _err("E_SCHEMA"), None)

# Trajectory members (spec § 5.5): rules every validator checks, schema or not.

_TRAJ_COORDS = [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
                0.0, 0.0, 0.1180, 0.0, 0.7580, -0.4700, 0.0, -0.7580, -0.4700]


def _trajectory(a: Archive, meta: dict, coords_dtype: str = "float64",
                extra: Optional[dict] = None) -> dict:
    """A 2-frame water trajectory with metadata `meta`, coords of
    `coords_dtype` and float64 `extra` members (role -> (values, shape))."""
    m = _json(meta)
    fmt = "<18f" if coords_dtype == "float32" else "<18d"
    coords = struct.pack(fmt, *_TRAJ_COORDS)
    members = {"metadata": json_member("trajectory/metadata.json", m),
               "coords": binary_member("trajectory/coords.bin", coords,
                                       coords_dtype, [2, 3, 3])}
    files = {"trajectory/metadata.json": m, "trajectory/coords.bin": coords}
    for role, (values, shape) in (extra or {}).items():
        data = struct.pack(f"<{len(values)}d", *values)
        members[role] = binary_member(f"trajectory/{role}.bin", data, "float64",
                                      shape)
        files[f"trajectory/{role}.bin"] = data
    return a.add_section({"id": "md", "kind": "trajectory", "members": members},
                         files)


_TRAJ_META = {"symbols": ["O", "H", "H"], "kind": "md", "timestep_fs": 0.5}


@case("trajectory_lattice_without_pbc")
def _():
    a = base_archive()
    cell = [9.0, 0.0, 0.0, 0.3, 9.0, 0.0, 0.2, 0.4, 9.0]
    _trajectory(a, _TRAJ_META, extra={"lattice": (cell, [3, 3])})
    return (a.to_bytes(),
            "A trajectory with a lattice member (rows, Angstrom) whose metadata "
            "has no pbc; pbc is REQUIRED when lattice is present (spec § 5.5).",
            _err("E_MEMBER_MALFORMED"), None)


@case("trajectory_float32_without_precision")
def _():
    a = base_archive()
    _trajectory(a, _TRAJ_META, coords_dtype="float32")
    return (a.to_bytes(),
            "Trajectory coords stored as float32 while the metadata sets no "
            "precision; float32 coords require metadata.precision 'float32' "
            "(spec § 5.5).",
            _err("E_MEMBER_MALFORMED"), None)


@case("trajectory_velocities_frame_count")
def _():
    a = base_archive()
    _trajectory(a, _TRAJ_META,
                extra={"velocities": ([0.001] * 27, [3, 3, 3])})
    return (a.to_bytes(),
            "Trajectory velocities of shape [3, 3, 3] for coords of 2 frames; "
            "every per-frame member has the leading extent of coords "
            "(spec § 5.5).",
            _err("E_MEMBER_MALFORMED"), None)


@case("trajectory_energies_disagree")
def _():
    a = base_archive()
    _trajectory(a, {**_TRAJ_META, "energies": [-76.4120, -76.4118]},
                extra={"energies": ([-76.4120, -76.4128], [2])})
    return (a.to_bytes(),
            "A trajectory whose binary energies member differs from "
            "metadata.energies by 1e-3 hartree in frame 1; the two must agree "
            "within 1e-10 hartree (spec § 5.5).",
            _err("E_MEMBER_MALFORMED"), None)


# ==== WP-29: property kinds (begin) ===========================================
# The rules of spec § 5.13 that JSON Schema cannot express: symmetric elastic
# matrices, static = electronic + ionic, one Born tensor per atom, and a dipole
# that agrees with the root dipole_moment. Every validator checks them, so each
# archive is invalid to the Python, C++ and JavaScript validators alike.

def _properties(a: Archive, kind: str, sid: str, payload: Any,
                **section_extra: Any) -> dict:
    """A property section: one JSON member `properties` (spec § 5.13)."""
    data = _json(payload)
    path = f"{sid}/properties.json"
    return a.add_section({"id": sid, "kind": kind, **section_extra,
                          "members": {"properties": json_member(path, data)}},
                         {path: data})


def _cubic_stiffness(c11: float = 36.4, c12: float = 8.8,
                     c44: float = 8.0) -> list:
    m = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        m[i][i], m[3 + i][3 + i] = c11, c44
        for j in range(3):
            if i != j:
                m[i][j] = c12
    return m


def _scaled_identity(x: float) -> list:
    return [[x if i == j else 0.0 for j in range(3)] for i in range(3)]


@case("elastic_not_symmetric")
def _():
    a = base_archive()
    stiffness = _cubic_stiffness()
    stiffness[1][0] = 9.9                      # entry (0, 1) is 8.8
    _properties(a, "properties.elastic", "elastic", {"stiffness": stiffness})
    return (a.to_bytes(),
            "A properties.elastic stiffness whose entry (0, 1) is 8.8 GPa and "
            "entry (1, 0) 9.9 GPa; the matrix must be symmetric to within 1e-6 "
            "of its largest entry (spec § 5.13).",
            _err("E_MEMBER_MALFORMED"), None)


@case("dielectric_static_mismatch")
def _():
    a = base_archive()
    _properties(a, "properties.dielectric", "dielectric", {
        "electronic": _scaled_identity(2.62), "ionic": _scaled_identity(4.58),
        "static": _scaled_identity(7.5)})     # electronic + ionic is 7.2
    return (a.to_bytes(),
            "A properties.dielectric whose static tensor has 7.5 on the "
            "diagonal while electronic + ionic is 7.2; static must equal "
            "their sum within 1e-6 (spec § 5.13).",
            _err("E_MEMBER_MALFORMED"), None)


@case("dielectric_born_atom_count")
def _():
    a = base_archive()
    _properties(a, "properties.dielectric", "dielectric", {
        "electronic": _scaled_identity(1.8),
        "born_charges": [_scaled_identity(-0.8), _scaled_identity(0.4)]})
    return (a.to_bytes(),
            "A properties.dielectric with 2 Born tensors for the 3-atom "
            "structure of the archive; there must be one per atom "
            "(spec § 5.13).",
            _err("E_MEMBER_MALFORMED"), None)


@case("electric_dipole_inconsistent")
def _():
    a = base_archive()
    a.manifest["dipole_moment"] = {"total_debye": 1.85,
                                   "vector_debye": [0.0, 0.0, -1.85],
                                   "origin": "center_of_mass"}
    _properties(a, "properties.electric", "electric", {
        "origin": "center_of_mass",
        "dipole": {"vector": [0.0, 0.0, -1.5], "total": 1.5, "units": "debye"}})
    return (a.to_bytes(),
            "A properties.electric dipole of 1.5 debye next to a root "
            "dipole_moment of 1.85 debye; the two state one dipole and must "
            "agree within 1e-6 debye (spec § 5.13).",
            _err("E_ROOT_INCONSISTENT"), None)
# ==== WP-29: property kinds (end) =============================================


# k-resolved Bloch orbitals (spec § 5.3.1). A minimal valid wavefunction.bloch
# section on a periodic hydrogen chain, with one defect per archive.

_CHAIN = {
    "atoms": [{"symbol": "H", "position": [0.0, 0.0, 0.0], "atomic_number": 1},
              {"symbol": "H", "position": [0.8, 0.0, 0.0], "atomic_number": 1}],
    "pbc": [True, False, False],
    "lattice_vectors": [[2.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
}


def _f64(shape: tuple, value: float = 0.5) -> bytes:
    n = 1
    for d in shape:
        n *= d
    return struct.pack(f"<{n}d", *([value] * n))


def _bloch_archive(*, structure: Any = None, weights: tuple = (0.5, 0.5),
                   coeff_shape: tuple = (2, 2, 2, 2, 2),
                   occ_shape: tuple = (2, 2, 2)) -> Archive:
    a = Archive("bloch")
    a.add_structure(payload=_CHAIN if structure is None else structure)
    docs = {
        "wf_bloch/basis.json": _json({
            "structure_ref": "structure", "pure": True, "n_ao": 2,
            "shells": [{"center": c, "l": 0, "exponents": [0.4166],
                        "coefficients": [1.0]} for c in (0, 1)]}),
        "wf_bloch/kpoints.json": _json({
            "kpoints": [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
            "weights": list(weights), "labels": ["G", "X"]}),
        "wf_bloch/mo_metadata.json": _json({
            "spin": "unrestricted", "n_spin": 2, "n_mo": 2, "n_ao": 2,
            "complex": True, "k_point_convention": "fractional"}),
    }
    arrays = {"mo_coefficients": coeff_shape, "energies": (2, 2, 2),
              "occupations": occ_shape}
    members = {role: json_member(path, docs[path]) for role, path in (
        ("basis", "wf_bloch/basis.json"), ("kpoints", "wf_bloch/kpoints.json"),
        ("mo_metadata", "wf_bloch/mo_metadata.json"))}
    files = dict(docs)
    for role, shape in arrays.items():
        data = _f64(shape)
        path = f"wf_bloch/{role}.bin"
        members[role] = binary_member(path, data, "float64", list(shape))
        files[path] = data
    a.add_section({"id": "wf_bloch", "kind": "wavefunction.bloch",
                   "members": members}, files)
    return a


@case("bloch_weights_sum")
def _():
    a = _bloch_archive(weights=(0.5, 0.25))
    return (a.to_bytes(),
            "A wavefunction.bloch kpoints payload whose two weights sum to "
            "0.75; the weights must sum to 1 within 1e-6 (spec § 5.3.1).",
            _err("E_MEMBER_MALFORMED"), None)


@case("bloch_coefficients_trailing_axis")
def _():
    a = _bloch_archive(coeff_shape=(2, 2, 2, 2, 1))
    return (a.to_bytes(),
            "wavefunction.bloch mo_coefficients of shape [2, 2, 2, 2, 1]: the "
            "last axis must be the (real, imaginary) pair, of length 2 (spec "
            "§ 5.3.1). Every validator reports the shape; one with a JSON "
            "Schema engine also reports the manifest schema violation.",
            _err("E_MEMBER_MALFORMED", "E_SCHEMA"), None)


@case("bloch_occupations_shape")
def _():
    a = _bloch_archive(occ_shape=(2, 2, 3))
    return (a.to_bytes(),
            "wavefunction.bloch occupations of shape [2, 2, 3] beside "
            "mo_coefficients [2, 2, 2, 2, 2]; occupations must be [n_k, "
            "n_spin, n_mo] = [2, 2, 2] (spec § 5.3.1).",
            _err("E_MEMBER_MALFORMED"), None)


@case("bloch_structure_molecular")
def _():
    a = _bloch_archive(structure=WATER)
    return (a.to_bytes(),
            "A wavefunction.bloch section whose basis.structure_ref names a "
            "molecular structure (no periodic axis in pbc); Bloch orbitals "
            "need a periodic structure (spec § 5.3.1).",
            _err("E_MEMBER_MALFORMED"), None)


# Sketch documents (spec § 5.1 chemical identity, § 5.8 job_ref, § 5.14
# reaction.scheme, § 5.15 preview).

def _png(width: int = 1, height: int = 1) -> bytes:
    """A valid grey PNG image."""
    import zlib

    def chunk(tag: bytes, data: bytes) -> bytes:
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))
    raw = b"".join(b"\x00" + b"\x80" * width for _ in range(height))
    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


_PNG = _png()


def _scheme(a: Archive, steps: list) -> dict:
    data = _json({"steps": steps})
    return a.add_section(
        {"id": "scheme", "kind": "reaction.scheme",
         "members": {"scheme": json_member("scheme/scheme.json", data)}},
        {"scheme/scheme.json": data})


def _preview(a: Archive, thumbnail: bytes, thumbnail_type: Optional[str],
             image: Optional[bytes] = None,
             image_type: Optional[str] = None) -> dict:
    members = {"thumbnail": opaque_member("preview/thumbnail.png", thumbnail)}
    files = {"preview/thumbnail.png": thumbnail}
    if thumbnail_type is not None:
        members["thumbnail"]["media_type"] = thumbnail_type
    if image is not None:
        members["image"] = opaque_member("preview/image.svg", image)
        members["image"]["media_type"] = image_type
        files["preview/image.svg"] = image
    return a.add_section({"id": "preview", "kind": "preview",
                          "members": members}, files)


def _two_waters() -> Archive:
    a = Archive("sketch")
    a.add_structure("water_a", "water_a/structure.json")
    a.add_structure("water_b", "water_b/structure.json")
    return a


@case("structure_charge_impossible")
def _():
    a = base_archive()
    a.manifest["sections"][0].update({"charge": 0, "multiplicity": 2})
    return (a.to_bytes(),
            "A water structure section with charge 0 and multiplicity 2: ten "
            "electrons cannot leave one unpaired, because the unpaired "
            "electrons must have the parity of the electron count (spec § 5.1, "
            "Chemical identity).",
            _err("E_CHARGE_MULTIPLICITY"), None)


@case("structure_charge_provenance")
def _():
    a = base_archive()
    a.manifest["sections"][0].update({"charge": 0, "multiplicity": 1})
    a.manifest["provenance"] = {"charge": 1, "multiplicity": 2}
    return (a.to_bytes(),
            "The archive's only structure section states charge 0 and "
            "multiplicity 1 while provenance states charge 1 and multiplicity "
            "2; with one structure the two describe the same species and must "
            "agree (spec § 5.1).",
            _err("E_ROOT_INCONSISTENT"), None)


@case("reaction_scheme_ref_wrong_kind")
def _():
    a = base_archive()
    _bonds(a)
    _scheme(a, [{"id": "r1", "reactants": [{"structure_ref": "structure"}],
                 "products": [{"structure_ref": "bonds"}]}])
    return (a.to_bytes(),
            "A reaction.scheme step whose product names an existing section of "
            "kind bonds; every participant must name a structure section "
            "(spec § 5.14).",
            _err("E_REF_UNRESOLVED"), None)


@case("reaction_scheme_step_id_duplicate")
def _():
    a = _two_waters()
    step = {"id": "r1", "reactants": [{"structure_ref": "water_a"}],
            "products": [{"structure_ref": "water_b"}]}
    _scheme(a, [step, dict(step)])
    return (a.to_bytes(),
            "Two steps of a reaction.scheme share the id r1; step ids are "
            "unique within the payload (spec § 5.14).",
            _err("E_MEMBER_MALFORMED"), None)


@case("reaction_scheme_participant_repeated")
def _():
    a = _two_waters()
    _scheme(a, [{"id": "r1",
                 "reactants": [{"structure_ref": "water_a"},
                               {"structure_ref": "water_a"}],
                 "products": [{"structure_ref": "water_b"}]}])
    return (a.to_bytes(),
            "A reaction.scheme step that lists water_a twice among its "
            "reactants; a structure occurs once per role of a step, and a "
            "coefficient says how many (spec § 5.14).",
            _err("E_MEMBER_MALFORMED"), None)


@case("reaction_scheme_energies_reserved")
def _():
    a = _two_waters()
    _scheme(a, [{"id": "r1", "reactants": [{"structure_ref": "water_a"}],
                 "products": [{"structure_ref": "water_b"}],
                 "energies": {"delta_e": -0.01}}])
    return (a.to_bytes(),
            "A reaction.scheme step whose energies object has an ad-hoc "
            "shape ({delta_e}) instead of $defs/ReactionStepEnergies: no "
            "version, deltas or terms, and an unknown key (spec § 5.14, "
            "Reaction energies). Only a validator that runs the JSON Schema "
            "reports it.",
            _err("E_PAYLOAD_SCHEMA"), None)


def _energies_step(**changes: Any) -> dict:
    """water_a -> water_b with a valid energies object; ``changes`` edits it.

    ``changes`` keys are ``arrow`` (on the step), ``delta`` (energy value_eh),
    ``coefficient_a``, ``method_b``, ``run_ref_b`` and ``thermal`` (adds a
    Gibbs delta; ``True`` with, ``"bare"`` without temperature/pressure).
    """
    method = {"program": "qvf-conformance", "program_version": "0.1.0",
              "specification": {"task": "energy", "method": "HF", "basis": "STO-3G"}}

    def term(ref, coefficient, energy, gibbs):
        t = {"structure_ref": ref, "coefficient": coefficient, "energy_eh": energy,
             "method": dict(method),
             "source": {"name": ref + ".out", "sha256": _sha(ref.encode()),
                        "parser": "conformance/1",
                        "locations": {"energy_eh": "line 10"}}}
        if changes.get("thermal"):
            t["gibbs_free_energy_eh"] = gibbs
        return t
    a = term("water_a", changes.get("coefficient_a", -1), -74.96, -74.95)
    b = term("water_b", 1, -74.95, -74.945)
    if "method_b" in changes:
        b["method"] = changes["method_b"]
    if "run_ref_b" in changes:
        b["run_ref"] = changes["run_ref_b"]
    energies: dict = {"version": 1,
                      "deltas": {"energy": {"value_eh": changes.get("delta", 0.01)}},
                      "terms": [a, b]}
    if changes.get("thermal"):
        energies["deltas"]["gibbs_free_energy"] = {"value_eh": 0.005}
        if changes["thermal"] != "bare":
            energies.update(temperature_k=298.15, pressure_atm=1.0)
    step = {"id": "r1", "reactants": [{"structure_ref": "water_a"}],
            "products": [{"structure_ref": "water_b"}], "energies": energies}
    if "arrow" in changes:
        step["arrow"] = changes["arrow"]
    return step


def _energies_case(description: str, codes: list, **changes: Any):
    a = _two_waters()
    _scheme(a, [_energies_step(**changes)])
    return a.to_bytes(), description + " (spec § 5.14, Reaction energies).", codes, None


@case("reaction_scheme_energies_sum")
def _():
    return _energies_case(
        "Step r1 water_a -> water_b states an energy delta of 0.02 hartree, "
        "but its terms -1 x -74.96 + 1 x -74.95 sum to 0.01",
        _err("E_MEMBER_MALFORMED"), delta=0.02)


@case("reaction_scheme_energies_sign")
def _():
    return _energies_case(
        "The reactant water_a has the term coefficient +1 instead of -1, so "
        "the delta, though equal to the sum of its terms, adds a reactant",
        _err("E_MEMBER_MALFORMED"), coefficient_a=1, delta=-149.91)


@case("reaction_scheme_energies_method")
def _():
    return _energies_case(
        "The two terms of step r1 come from different methods (HF/STO-3G and "
        "HF/6-31G); a reaction energy needs one method",
        _err("E_MEMBER_MALFORMED"),
        method_b={"program": "qvf-conformance", "program_version": "0.1.0",
                  "specification": {"task": "energy", "method": "HF", "basis": "6-31G"}})


@case("reaction_scheme_energies_resonance")
def _():
    return _energies_case(
        "A resonance step with energies: two resonance structures of one "
        "species have no reaction energy",
        _err("E_MEMBER_MALFORMED"), arrow="resonance")


@case("reaction_scheme_energies_run_ref")
def _():
    return _energies_case(
        "A term's run_ref names the structure section water_b instead of a "
        "run.record section",
        _err("E_REF_UNRESOLVED"), run_ref_b="water_b")


@case("reaction_scheme_energies_thermal_conditions")
def _():
    return _energies_case(
        "A Gibbs free energy delta without temperature_k and pressure_atm. "
        "Only a validator that runs the JSON Schema reports it",
        _err("E_PAYLOAD_SCHEMA"), thermal="bare")


@case("preview_signature_mismatch")
def _():
    a = base_archive()
    _preview(a, b"\xff\xd8\xff\xe0 not a png", "image/png")
    return (a.to_bytes(),
            "A preview thumbnail whose media_type is image/png but whose bytes "
            "begin with the JPEG marker FF D8 FF instead of the PNG signature "
            "(spec § 5.15).",
            _err("E_MEMBER_MALFORMED"), None)


@case("preview_svg_not_utf8")
def _():
    a = base_archive()
    _preview(a, _PNG, "image/png", b"<svg>\xff\xfe</svg>", "image/svg+xml")
    return (a.to_bytes(),
            "A preview image of media_type image/svg+xml whose bytes do not "
            "decode as UTF-8 (spec § 5.15).",
            _err("E_UTF8"), None)


@case("preview_media_type_missing")
def _():
    a = base_archive()
    _preview(a, _PNG, None)
    return (a.to_bytes(),
            "A preview thumbnail without media_type; preview members must say "
            "what their bytes are (spec § 5.15). Only a validator that runs "
            "the JSON Schema reports it.",
            _err("E_SCHEMA"), None)


@case("run_record_job_ref_wrong_kind")
def _():
    a = base_archive()
    log = b"done\n"
    a.add_section({"id": "run", "kind": "run.record", "program": "test-code",
                   "job_ref": "structure",
                   "members": {"log": opaque_member("run/log.txt", log)}},
                  {"run/log.txt": log})
    return (a.to_bytes(),
            "A run.record whose job_ref names the structure section instead of "
            "a job.spec section (spec § 5.8).",
            _err("E_REF_UNRESOLVED"), None)


# Toolkit 0.3.0: structure_ref on a standalone bonds section (spec § 5.1), the
# IAO reference basis of localized orbitals (§ 5.3) and the reserved
# viewer_defaults members camera and crossfade (§ 3.4).

_H2 = {
    "atoms": [
        {"symbol": "H", "position": [0.0, 0.0, 0.0], "atomic_number": 1},
        {"symbol": "H", "position": [0.0, 0.0, 0.74], "atomic_number": 1},
    ],
    "pbc": [False, False, False],
    "lattice_vectors": None,
}


def _water_and_h2() -> Archive:
    a = Archive("bonds")
    a.add_structure("structure_water", "structure_water/structure.json")
    a.add_structure("structure_h2", "structure_h2/structure.json", payload=_H2)
    return a


def _standalone_bonds(a: Archive, pairs: list, structure_ref: str) -> dict:
    data = _json({"pairs": pairs})
    return a.add_section({"id": "bonds", "kind": "bonds",
                          "structure_ref": structure_ref,
                          "members": {"bonds": json_member(
                              "bonds/connectivity.json", data)}},
                         {"bonds/connectivity.json": data})


@case("bonds_structure_ref_unresolved")
def _():
    a = _water_and_h2()
    _standalone_bonds(a, [{"i": 0, "j": 1}], "structure_ethanol")
    return (a.to_bytes(),
            "A standalone bonds section whose structure_ref names no section "
            "of the archive (spec § 4.1, § 5.1).",
            _err("E_REF_UNRESOLVED"), None)


@case("bonds_structure_ref_wrong_kind")
def _():
    a = _water_and_h2()
    spec = _json({"job_type": "molecular", "structure_ref": "structure_h2"})
    a.add_section({"id": "job", "kind": "job.spec",
                   "structure_ref": "structure_h2",
                   "members": {"spec": json_member("job/spec.json", spec)}},
                  {"job/spec.json": spec})
    _standalone_bonds(a, [{"i": 0, "j": 1}], "job")
    return (a.to_bytes(),
            "A standalone bonds section whose structure_ref names a job.spec "
            "section instead of a structure section (spec § 4.1, § 5.1).",
            _err("E_REF_UNRESOLVED"), None)


@case("bonds_atom_out_of_range")
def _():
    a = _water_and_h2()
    _standalone_bonds(a, [{"i": 0, "j": 1}, {"i": 1, "j": 2, "order": 1.0}],
                      "structure_h2")
    return (a.to_bytes(),
            "A standalone bonds section whose structure_ref names the second "
            "structure, H2 with atoms 0 and 1, and whose second pair joins "
            "atoms 1 and 2. Atom 2 exists in the first structure, water, so "
            "a validator that counted the atoms of the first structure "
            "instead of the referenced one would accept it (spec § 5.1).",
            _err("E_MEMBER_MALFORMED"), None)


@case("payload_localization_reference_basis_empty")
def _():
    a = base_archive()
    _wavefunction(a, _H2_BASIS, {"n_mo": 2, "n_ao": 2, "spin": "restricted",
                                 "orbital_kind": "ibo",
                                 "energies": [-0.58, 0.67],
                                 "occupations": [2.0, 0.0],
                                 "localization_method": "ibo",
                                 "localization_reference_basis": ""})
    return (a.to_bytes(),
            "Localized orbitals whose mo_metadata names the IAO reference "
            "basis with an empty string; localization_reference_basis is a "
            "non-empty name ($defs/MoMetadataPayload, spec § 5.3).",
            _err("E_PAYLOAD_SCHEMA"), None)


@case("reserved_id_camera")
def _():
    a = Archive("reserved-id")
    a.add_structure(sid="camera")
    return (a.to_bytes(),
            "A section id 'camera', a reserved viewer_defaults member "
            "(spec § 3.1, § 3.4).",
            _err("E_RESERVED_ID"), None)


@case("reserved_id_crossfade")
def _():
    a = Archive("reserved-id")
    a.add_structure(sid="crossfade")
    return (a.to_bytes(),
            "A section id 'crossfade', a reserved viewer_defaults member "
            "(spec § 3.1, § 3.4).",
            _err("E_RESERVED_ID"), None)


def _unrestricted_nto_archive(defect):
    """Two raw NTO sections; isolate one issue-12 contract defect."""
    import copy
    a = base_archive()
    basis = copy.deepcopy(_H2_BASIS)
    for role in ("hole", "particle"):
        if defect == "partner_missing" and role == "particle":
            continue
        sid = "nto_" + role
        weights = {"alpha": [0.55, 0.15], "beta": [0.3]}
        rows = {"alpha": 2, "beta": 1}
        if role == "particle":
            if defect == "rows": rows["alpha"] = 1; weights["alpha"] = [0.7]
            if defect == "weights": weights["beta"] = [0.4]
        if defect == "order": weights["alpha"] = [0.15, 0.55]
        if defect == "negative": weights["beta"] = [-0.3]
        meta = {"spin": "unrestricted", "orbital_kind": "nto",
                "occupation_semantics": "transition_weight", "nto_role": role,
                "state_ref": "es", "state_index": 2 if defect == "state" else 1,
                "n_ao": 2}
        for ch in ("alpha", "beta"):
            meta[ch] = {"n_mo": rows[ch], "energies": [0.0] * rows[ch], "occupations": weights[ch]}
        if defect == "metadata_rows" and role == "hole": meta["alpha"]["n_mo"] = 3
        if defect == "energies_rows" and role == "hole": meta["beta"]["energies"] = []
        if defect == "missing_role" and role == "hole": del meta["nto_role"]
        if defect == "missing_state" and role == "hole": del meta["state_index"]
        members, files = {}, {}
        for key, value in (("basis", basis), ("mo_metadata", meta)):
            path = sid + "/" + key + ".json"; data = _json(value)
            members[key] = json_member(path, data); files[path] = data
        for ch in ("alpha", "beta"):
            key = "mo_coefficients_" + ch; path = sid + "/" + key + ".dat"
            data = struct.pack("<" + "d" * (rows[ch] * 2), *([1.0, 0.0] * rows[ch]))
            members[key] = binary_member(path, data, "float64", [rows[ch], 2]); files[path] = data
        a.add_section({"id": sid, "kind": "wavefunction.gto", "members": members}, files)
    path = "es/states.json"
    data = _json({"states": [{"index": 1, "excitation_energy": {"value": 4.0, "units": "ev"}, "nto_ref": "nto_hole"}]})
    a.add_section({"id": "es", "kind": "excited_states", "members": {"states": json_member(path, data)}}, {path: data})
    return a


for _defect in ("rows", "weights", "order", "negative", "metadata_rows", "energies_rows", "partner_missing", "state", "missing_role", "missing_state"):
    def _make_nto_case(defect=_defect):
        code = "E_REF_UNRESOLVED" if defect in ("partner_missing", "state", "missing_role") else "E_PAYLOAD_SCHEMA" if defect == "missing_state" else "E_MEMBER_MALFORMED"
        expected = _err(code)
        if defect == "missing_state":
            expected += _err("E_REF_UNRESOLVED")  # The particle cannot resolve its pair.
        return (_unrestricted_nto_archive(defect).to_bytes(),
                "Unrestricted NTO defect: " + defect + "; spec § 5.3.", expected, None)
    case("nto_unrestricted_" + _defect)(_make_nto_case)


def build(out_dir: str = OUT) -> list[str]:
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for name, fn in CASES.items():
        data, description, expect, limits = fn()
        assert len(data) < 100 * 1024, (name, len(data))
        with open(os.path.join(out_dir, f"{name}.qvf"), "wb") as fh:
            fh.write(data)
        expected: dict[str, Any] = {"archive": f"{name}.qvf",
                                    "description": description}
        if limits:
            expected["limits"] = limits
        expected["expect"] = expect
        with open(os.path.join(out_dir, f"{name}.expected.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(expected, fh, indent=2)
            fh.write("\n")
        written.append(name)
        print(f"  {name}.qvf  ({len(data)} bytes) -> "
              + ", ".join(e["code"] for e in expect))
    return written


if __name__ == "__main__":
    build()
