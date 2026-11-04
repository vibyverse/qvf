"""Regression tests for defects found by fuzzing the Python readers and validator.

Each archive is built in memory. The fuzz harnesses themselves are opt-in
(python/tests/fuzz, QVF_FUZZ=1); these tests pin what they found.
"""
from __future__ import annotations

import hashlib
import io
import json
import struct
import sys
import tracemalloc
import zipfile
import zlib
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from qvf import validate as V  # noqa: E402
from qvf.reader import Reader  # noqa: E402

WATER = {"atoms": [{"symbol": "O", "position": [0, 0, 0], "atomic_number": 8}],
         "pbc": [False, False, False], "lattice_vectors": None}


def _manifest(members):
    return {"qvf_version": 1,
            "source": {"program": "t", "version": "1", "calculation": "t"},
            "sections": [{"id": "x", "kind": "x_test.blob", "members": members}]}


def _archive(payload: bytes, *, method=zipfile.ZIP_DEFLATED, manifest=None,
             member=None) -> bytes:
    buf = io.BytesIO()
    member = member or {"path": "x/blob.json", "format": "json"}
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(zipfile.ZipInfo("x/blob.json"), payload, compress_type=method)
        zf.writestr("manifest.json", manifest if manifest is not None
                    else json.dumps(_manifest({"blob": member})))
    return buf.getvalue()


def _patch_entry(raw: bytes, name: bytes, *, usize=None, crc=None) -> bytes:
    """Rewrite one entry's declared size and CRC-32 in both headers."""
    out = bytearray(raw)
    for sig, crc_off, size_off, nlen_off, name_off in (
            (b"PK\x03\x04", 14, 22, 26, 30), (b"PK\x01\x02", 16, 24, 28, 46)):
        pos = out.find(sig)
        while pos != -1:
            n = struct.unpack_from("<H", out, pos + nlen_off)[0]
            if bytes(out[pos + name_off:pos + name_off + n]) == name:
                if usize is not None:
                    struct.pack_into("<I", out, pos + size_off, usize)
                if crc is not None:
                    struct.pack_into("<I", out, pos + crc_off, crc)
            pos = out.find(sig, pos + 4)
    return bytes(out)


def _codes(report):
    return {f["code"] for f in report["findings"]}


def _lying_bomb() -> bytes:
    """Declares 1 KiB with the CRC-32 of that 1 KiB; the stream holds 64 MiB."""
    payload = b"[" + b" " * ((64 << 20) - 2) + b"]"
    raw = _archive(payload)
    return _patch_entry(raw, b"x/blob.json", usize=1024, crc=zlib.crc32(payload[:1024]))


def test_lying_stream_is_not_inflated_past_its_declared_size():
    raw = _lying_bomb()
    assert len(raw) < 200_000
    tracemalloc.start()
    try:
        report = V.validate_qvf(io.BytesIO(raw))
        with Reader.from_bytes(raw) as reader:
            with pytest.raises(ValueError, match=r"^E_ZIP_UNREADABLE"):
                reader.json("x", "blob")
        with V.QvfReader(io.BytesIO(raw)) as legacy:
            with pytest.raises(ValueError, match=r"^E_ZIP_UNREADABLE"):
                legacy.read_member(legacy.sections[0]["members"]["blob"])
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert "E_ZIP_UNREADABLE" in _codes(report)
    # zipfile.read() used to inflate the whole 64 MiB stream (twice over).
    assert peak < 24 << 20, peak


def test_stream_shorter_than_declared_is_unreadable():
    payload = json.dumps({"a": "x" * 1000}).encode()
    raw = _patch_entry(_archive(payload), b"x/blob.json", usize=2 * len(payload))
    assert "E_ZIP_UNREADABLE" in _codes(V.validate_qvf(io.BytesIO(raw)))
    with Reader.from_bytes(raw) as reader, \
            pytest.raises(ValueError, match=r"^E_ZIP_UNREADABLE"):
        reader.json("x", "blob")


@pytest.mark.parametrize("method", [zipfile.ZIP_BZIP2, zipfile.ZIP_LZMA])
def test_only_store_and_deflate_are_read(method):
    raw = _archive(json.dumps(WATER).encode(), method=method)
    report = V.validate_qvf(io.BytesIO(raw))
    assert "E_ZIP_UNREADABLE" in _codes(report), report["findings"]
    with Reader.from_bytes(raw) as reader, \
            pytest.raises(ValueError, match=r"^E_ZIP_UNREADABLE.*compression method"):
        reader.json("x", "blob")


@pytest.mark.parametrize("damage", ["deflate", "method", "encrypted", "patched"])
def test_damaged_entries_raise_coded_value_errors(damage):
    payload = json.dumps(WATER).encode()
    raw = bytearray(_archive(payload))
    local = raw.find(b"PK\x03\x04")
    central = raw.find(b"PK\x01\x02")
    if damage == "deflate":
        start = local + 30 + len(b"x/blob.json")
        raw[start] ^= 0xFF            # invalid block type: a zlib.error before
    elif damage == "method":
        struct.pack_into("<H", raw, central + 10, 99)   # NotImplementedError
    elif damage == "encrypted":
        for at in (local + 6, central + 8):
            raw[at] |= 0x1                                # RuntimeError
    else:
        for at in (local + 6, central + 8):
            raw[at] |= 0x20                               # NotImplementedError
    raw = bytes(raw)
    assert "E_ZIP_UNREADABLE" in _codes(V.validate_qvf(io.BytesIO(raw)))
    with Reader.from_bytes(raw) as reader, \
            pytest.raises(ValueError, match=r"^E_ZIP_UNREADABLE"):
        reader.json("x", "blob")


def test_not_a_zip_is_a_coded_value_error_and_still_a_bad_zip_file():
    for raw in (b"not a zip", b"PK\x05\x06" + b"\xff" * 18):
        with pytest.raises(ValueError, match=r"^E_ZIP_UNREADABLE") as info:
            Reader.from_bytes(raw)
        assert isinstance(info.value, zipfile.BadZipFile)
        with pytest.raises(zipfile.BadZipFile):
            V.QvfReader(io.BytesIO(raw))


@pytest.mark.parametrize("text", [
    b'{"qvf_version": NaN}',
    b'{"qvf_version": 1, "x": -Infinity}',
    json.dumps(_manifest({})).encode("utf-16"),
    json.dumps(_manifest({})).encode("utf-32-le"),
])
def test_manifest_must_be_rfc8259_utf8_json(text):
    raw = _archive(b"{}", manifest=text)
    assert "E_MANIFEST_JSON" in _codes(V.validate_qvf(io.BytesIO(raw)))
    with pytest.raises(ValueError, match=r"^E_MANIFEST_JSON"):
        Reader.from_bytes(raw)
    with pytest.raises(ValueError, match=r"^E_MANIFEST_JSON"):
        V.QvfReader(io.BytesIO(raw))


def test_member_json_must_be_rfc8259_utf8_json():
    raw = _archive(b'{"a": NaN}')
    assert "E_JSON_INVALID" in _codes(V.validate_qvf(io.BytesIO(raw)))
    with Reader.from_bytes(raw) as reader, \
            pytest.raises(ValueError, match=r"^E_JSON_INVALID"):
        reader.json("x", "blob")


def test_utf8_bom_is_still_accepted():
    raw = _archive(b"{}", manifest=b"\xef\xbb\xbf" + json.dumps(_manifest({})).encode())
    assert "E_MANIFEST_JSON" not in _codes(V.validate_qvf(io.BytesIO(raw)))
    Reader.from_bytes(raw).close()


def test_legacy_reader_errors_carry_codes_and_keep_their_types():
    raw = _archive(b'{"a": ', member={"path": "x/blob.json", "format": "json",
                                       "sha256": "0" * 64})
    with V.QvfReader(io.BytesIO(raw)) as legacy:
        member = legacy.sections[0]["members"]["blob"]
        with pytest.raises(ValueError, match=r"^E_SHA256_MISMATCH: sha256 mismatch"):
            legacy.read_member(member)
        with pytest.raises(json.JSONDecodeError, match=r"^E_JSON_INVALID"):
            legacy.read_json({"path": "x/blob.json", "format": "json"})
    huge = b'{"qvf_version": 1, "x": 1' + b"0" * 5000 + b"}"
    with pytest.raises(ValueError, match=r"^E_MANIFEST_JSON"):
        V.QvfReader(io.BytesIO(_archive(b"{}", manifest=huge)))


def test_tolist_of_a_huge_zero_extent_shape_is_refused():
    member = {"path": "x/blob.json", "format": "binary", "dtype": "float64",
              "shape": [2 ** 40, 0], "sha256": hashlib.sha256(b"").hexdigest()}
    raw = _archive(b"", member=member)
    assert V.validate_qvf(io.BytesIO(raw))["ok"]
    with Reader.from_bytes(raw) as reader:
        tensor = reader.tensor("x", "blob")
        assert len(tensor.data) == 0
        with pytest.raises(ValueError, match=r"tolist\(\).*nested empty lists"):
            tensor.tolist()
    small = Reader.from_bytes(_archive(b"", member=dict(member, shape=[3, 0, 2])))
    with small:
        assert small.tensor("x", "blob").tolist() == [[], [], []]
