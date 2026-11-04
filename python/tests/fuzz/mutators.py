"""Structure-aware mutations of QVF archives, at the ZIP and the JSON level.

Every operator takes a `random.Random` and an `Archive` model and edits the
model in place. `mutate(raw, ops)` applies a list of `(name, seed)` pairs to a
seed archive and returns the new bytes, so Hypothesis can draw and shrink the
operator list and atheris can drive the same operators from its seed.
"""
from __future__ import annotations

import copy
import functools
import hashlib
import json
import random
import struct
import zlib
from pathlib import Path

from zipmodel import Archive, Entry, parse, store, SAT32

ROOT = Path(__file__).resolve().parents[3]
CORPUS_DIRS = [ROOT / "conformance" / d for d in ("corpus", "negative", "legacy")]

INTERESTING_INTS = [0, 1, -1, 2, 3, 7, 255, 256, 0x7FFF, 0xFFFF, 0x10000, 2 ** 31 - 1,
                    2 ** 31, 2 ** 32 - 1, 2 ** 32, 2 ** 53, 2 ** 53 + 1, 2 ** 63 - 1,
                    2 ** 63, 2 ** 64 - 1, 2 ** 64, 10 ** 30, -(2 ** 63)]
INTERESTING_SIZES = [0, 1, 2, 30, 0xFFFF, 0x10000, 2 ** 24, 2 ** 30, 2 ** 31 - 1, 2 ** 31,
                     0xFFFFFFFE, SAT32, 2 ** 33, 8 * 1024 ** 3, 8 * 1024 ** 3 + 1,
                     2 ** 53 - 1, 2 ** 53, 2 ** 62, 2 ** 64 - 1]
EVIL_NAMES = ["../escape", "a/../../b", "/abs/path", "C:/drive", "c:relative", "a\\b",
              "back\\..\\slash", "nul\x00byte", "dir/", "", ".", "..", "./manifest.json",
              "manifest.json", "MANIFEST.JSON", "manifest.json/", "a//b", "\u00e9t\u00e9",
              "\ud800".encode("utf-8", "surrogatepass").decode("latin-1"),
              "x" * 300, "x" * 0xFFFF]
EVIL_STRINGS = ["", " ", "x", "\u0000", "\ud83d\ude00", "float128", "complex64", "<f8",
                "bool", "float64 ", "Float64", "json", "binary", "opaque", "../x", "/x",
                "0" * 64, "A" * 64, "g" * 64, "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b"
                "934ca495991b7852b855", "x" * 100000]


def seeds() -> list[tuple[str, bytes]]:
    out = []
    for d in CORPUS_DIRS:
        for p in sorted(d.glob("*.qvf")):
            out.append((f"{d.name}/{p.name}", p.read_bytes()))
    return out


def _entry(rng: random.Random, a: Archive) -> Entry | None:
    return rng.choice(a.entries) if a.entries else None


def _copy(rng) -> str:
    return rng.choice(["local", "central", "both"])


def _set(e: Entry, where: str, key: str, value) -> None:
    if where in ("local", "both"):
        e.local[key] = value
    if where in ("central", "both"):
        e.central[key] = value


# --- ZIP-level operators -----------------------------------------------------

def op_method(rng, a):
    e = _entry(rng, a)
    if e:
        _set(e, _copy(rng), "method", rng.choice([0, 8, 1, 6, 9, 12, 14, 93, 95, 98, 99,
                                                   rng.randrange(0x10000)]))


def op_flags(rng, a):
    e = _entry(rng, a)
    if e:
        bit = rng.choice([0x1, 0x8, 0x40, 0x800, 0x2000, 0x20, 1 << rng.randrange(16)])
        _set(e, _copy(rng), "flags", e.flags ^ bit)


def op_crc(rng, a):
    e = _entry(rng, a)
    if e:
        _set(e, _copy(rng), "crc", rng.choice([0, SAT32, e.crc ^ (1 << rng.randrange(32)),
                                               rng.randrange(2 ** 32)]))


def _size_value(rng, true: int) -> int:
    r = rng.random()
    if r < 0.4:
        return max(0, true + rng.choice([-1, 1, -2, 2, -1024, 1024, rng.randint(-true, 4096)]))
    if r < 0.8:
        return rng.choice(INTERESTING_SIZES)
    return rng.randrange(2 ** 32)


def op_csize(rng, a):
    e = _entry(rng, a)
    if e:
        _set(e, _copy(rng), "csize", _size_value(rng, len(e.data)) & SAT32)


def op_usize(rng, a):
    e = _entry(rng, a)
    if e:
        _set(e, _copy(rng), "usize", _size_value(rng, e.usize) & SAT32)


def op_zip64(rng, a):
    """Saturate fields and add a ZIP64 extra field, with true or lying values."""
    e = _entry(rng, a)
    if not e:
        return
    keys = [k for k in ("usize", "csize", "offset") if rng.random() < 0.6] or ["usize"]
    e.zip64 = tuple(keys)
    if rng.random() < 0.5:
        k = rng.choice(keys)
        e.zip64_values[k] = rng.choice(INTERESTING_SIZES + [2 ** 63, 2 ** 64 - 1])
    if rng.random() < 0.15:
        e.zip64_truncate = rng.randint(1, 8 * len(keys))
    if rng.random() < 0.1:
        e.zip64_values["length"] = rng.choice([0, 4, 7, 9, 0xFFFF])
    if rng.random() < 0.3:
        a.zip64_end = True


def op_zip64_end(rng, a):
    a.zip64_end = True
    r = rng.random()
    if r < 0.2:
        a.end64["count"] = rng.choice(INTERESTING_SIZES)
    elif r < 0.4:
        a.end64["cd_offset"] = rng.choice(INTERESTING_SIZES)
    elif r < 0.5:
        a.end64["cd_size"] = rng.choice(INTERESTING_SIZES)
    elif r < 0.6:
        a.end64["locator_offset"] = rng.choice(INTERESTING_SIZES)
    elif r < 0.7:
        a.end64["disks"] = rng.choice([0, 2, SAT32])
    elif r < 0.8:
        a.end64["disk"] = rng.choice([1, SAT32])
    elif r < 0.9:
        a.end64["record_size"] = rng.choice([0, 12, 2 ** 63])


def op_descriptor(rng, a):
    e = _entry(rng, a)
    if not e:
        return
    e.flags |= 0x8
    e.descriptor = rng.choice(["sig32", "nosig32", "sig64"])
    e.zero_local_sizes = rng.random() < 0.7
    if rng.random() < 0.4:
        key = rng.choice(["crc", "csize", "usize"])
        e.descriptor_values[key] = rng.choice(INTERESTING_SIZES)
    if rng.random() < 0.1:
        e.flags &= ~0x8   # a descriptor the flags do not announce


def op_duplicate(rng, a):
    e = _entry(rng, a)
    if not e:
        return
    d = copy.deepcopy(e)
    if rng.random() < 0.3:
        d.data = d.data[: len(d.data) // 2]
    a.entries.insert(rng.randrange(len(a.entries) + 1), d)


def op_duplicate_manifest(rng, a):
    m = a.find(b"manifest.json")
    if m:
        d = copy.deepcopy(m)
        a.entries.insert(rng.randrange(len(a.entries) + 1), d)


def op_rename(rng, a):
    e = _entry(rng, a)
    if not e:
        return
    r = rng.random()
    if r < 0.6:
        name = rng.choice(EVIL_NAMES).encode("utf-8", "surrogatepass")
    elif r < 0.8 and len(a.entries) > 1:
        name = rng.choice(a.entries).name
    else:
        name = bytes(rng.randrange(256) for _ in range(rng.randint(0, 12)))
    where = rng.choice(["both", "both", "local", "central"])
    if where == "both":
        e.name = name
    else:
        _set(e, where, "name", name)
    if rng.random() < 0.3:
        e.flags ^= 0x800


def op_offset(rng, a):
    e = _entry(rng, a)
    if e:
        e.central["offset"] = rng.choice(INTERESTING_SIZES + [rng.randrange(4096)])


def op_corrupt_data(rng, a):
    e = _entry(rng, a)
    if not e or not e.data:
        return
    d = bytearray(e.data)
    r = rng.random()
    if r < 0.5:
        for _ in range(rng.randint(1, 4)):
            d[rng.randrange(len(d))] ^= 1 << rng.randrange(8)
    elif r < 0.75:
        del d[rng.randrange(len(d)):]
    else:
        d += bytes(rng.randrange(256) for _ in range(rng.randint(1, 64)))
    e.data = bytes(d)


def _payload(e: Entry) -> bytes | None:
    try:
        if e.method == 0:
            return e.data
        if e.method == 8:
            return zlib.decompress(e.data, -15)
    except zlib.error:
        return None
    return None


def op_recompress(rng, a):
    """Switch between STORE and DEFLATE, keeping the entry valid."""
    e = _entry(rng, a)
    if not e:
        return
    payload = _payload(e)
    if payload is None:
        return
    fresh = store("x", payload, deflate=e.method == 0 or rng.random() < 0.3)
    e.data, e.method = fresh.data, fresh.method


def op_drop(rng, a):
    if a.entries:
        a.entries.pop(rng.randrange(len(a.entries)))


def op_reorder(rng, a):
    rng.shuffle(a.entries)


def op_extra(rng, a):
    e = _entry(rng, a)
    if not e:
        return
    blob = bytearray()
    for _ in range(rng.randint(1, 3)):
        tag = rng.choice([0x0001, 0x5455, 0x7875, 0x000A, 0x9901, rng.randrange(0x10000)])
        body = bytes(rng.randrange(256) for _ in range(rng.randint(0, 24)))
        length = len(body) if rng.random() < 0.7 else rng.randrange(0x10000)
        blob += struct.pack("<HH", tag, length) + body
    if rng.random() < 0.5:
        e.central_extra = bytes(blob)
    else:
        e.local_extra = bytes(blob)
    if rng.random() < 0.1:
        _set(e, _copy(rng), "extra_length", rng.choice([0, 1, 3, 0xFFFF]))


def op_end(rng, a):
    r = rng.random()
    n = len(a.entries)
    if r < 0.2:
        a.end["count"] = rng.choice([0, max(0, n - 1), n + 1, 0xFFFF])
    elif r < 0.35:
        a.end["cd_size"] = rng.choice([0, 1, 45, 46, SAT32, rng.randrange(1 << 20)])
    elif r < 0.5:
        a.end["cd_offset"] = rng.choice([0, 1, SAT32, rng.randrange(1 << 20)])
    elif r < 0.6:
        a.end["comment_length"] = rng.choice([1, 100, 0xFFFF])
    elif r < 0.7:
        a.end["disk"] = rng.choice([1, 0xFFFF])
    elif r < 0.8:
        a.comment = b"PK\x05\x06" + bytes(rng.randrange(256) for _ in range(18))
    elif r < 0.9:
        a.prefix = bytes(rng.randrange(256) for _ in range(rng.randint(1, 64)))
    else:
        a.suffix = bytes(rng.randrange(256) for _ in range(rng.randint(1, 64)))


@functools.lru_cache(maxsize=None)
def _zeros(n: int) -> tuple[bytes, int]:
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    zeros = bytes(n)
    return c.compress(zeros) + c.flush(), zlib.crc32(zeros) & SAT32


def op_bomb(rng, a):
    """Replace an entry with a highly compressible stream, honest or lying."""
    e = _entry(rng, a)
    if not e:
        return
    n = rng.choice([1 << 20, 1 << 24, 1 << 26])
    e.data, e.crc = _zeros(n)
    e.method = 8
    r = rng.random()
    if r < 0.4:
        e.usize = n
    elif r < 0.7:
        e.usize = rng.choice([n // 2, 1024, 0])      # stream longer than declared
    else:
        e.usize = rng.choice([n * 2, 2 ** 31, 0xFFFFFFFE])   # stream shorter


def op_declare_huge(rng, a):
    """Declare a huge uncompressed size for a tiny entry, ZIP64 or not."""
    e = _entry(rng, a)
    if not e:
        return
    size = rng.choice([2 ** 30, 2 ** 31, 0xFFFFFFFE, 2 ** 33, 7 * 1024 ** 3, 8 * 1024 ** 3])
    if size > 0xFFFFFFFE or rng.random() < 0.3:
        e.zip64 = ("usize",)
        e.zip64_values["usize"] = size
        e.usize = size if size <= SAT32 else e.usize
    else:
        e.usize = size
    e.local.pop("usize", None)


GROUP_ZIP = [op_method, op_flags, op_crc, op_csize, op_usize, op_zip64, op_zip64_end,
             op_descriptor, op_duplicate, op_duplicate_manifest, op_rename, op_offset,
             op_corrupt_data, op_recompress, op_drop, op_reorder, op_extra, op_end,
             op_bomb, op_declare_huge]


# --- manifest / JSON-level operators ------------------------------------------

def _deep(depth: int, leaf=0):
    v = leaf
    for i in range(depth):
        v = [v] if i % 2 else {"k": v}
    return v


def _evil_value(rng):
    r = rng.random()
    if r < 0.25:
        return rng.choice(INTERESTING_INTS)
    if r < 0.4:
        return rng.choice([None, True, False, 0.5, -0.0, 1e308, -1e308, 5e-324, 1.5,
                           float(2 ** 53), 1e300 * 10])
    if r < 0.55:
        return rng.choice(EVIL_STRINGS)
    if r < 0.7:
        return rng.choice([[], {}, [None], [[]], {"": None}, [0] * 1000])
    if r < 0.8:
        return _deep(rng.choice([10, 62, 63, 64, 65, 200]))
    if r < 0.9:
        return [rng.choice(INTERESTING_INTS) for _ in range(rng.randint(0, 6))]
    return {"path": rng.choice(EVIL_STRINGS), "format": rng.choice(["json", "binary"])}


def _nodes(doc, path=()):
    yield path, doc
    if isinstance(doc, dict):
        for k, v in doc.items():
            yield from _nodes(v, path + (k,))
    elif isinstance(doc, list):
        for i, v in enumerate(doc):
            yield from _nodes(v, path + (i,))


def _replace(doc, path, value):
    if not path:
        return value
    parent = doc
    for p in path[:-1]:
        parent = parent[p]
    parent[path[-1]] = value
    return doc


def _manifest(a: Archive):
    """(entry, parsed manifest) of the last manifest.json, or (None, None)."""
    e = None
    for x in a.entries:
        if x.name == b"manifest.json":
            e = x
    if e is None:
        return None, None
    payload = _payload(e)
    if payload is None:
        return e, None
    try:
        return e, json.loads(payload)
    except ValueError:
        return e, None


def _write_json(e: Entry, doc, rng) -> None:
    text = json.dumps(doc, allow_nan=True, ensure_ascii=rng.random() < 0.5)
    fresh = store("x", text.encode("utf-8", "surrogatepass"), deflate=e.method == 8)
    e.data, e.method, e.crc, e.usize = fresh.data, fresh.method, fresh.crc, fresh.usize
    e.local.clear()
    e.central.clear()


def _members(doc):
    if not isinstance(doc, dict) or not isinstance(doc.get("sections"), list):
        return []
    out = []
    for s in doc["sections"]:
        if isinstance(s, dict) and isinstance(s.get("members"), dict):
            for role, m in s["members"].items():
                if isinstance(m, dict):
                    out.append((s, role, m))
    return out


def op_json_node(rng, a):
    """Replace, delete or add a random node of manifest.json."""
    e, doc = _manifest(a)
    if doc is None:
        return
    nodes = list(_nodes(doc))
    path, node = rng.choice(nodes)
    r = rng.random()
    if r < 0.6 or not path:
        doc = _replace(doc, path, _evil_value(rng))
    elif r < 0.8:
        parent = doc
        for p in path[:-1]:
            parent = parent[p]
        del parent[path[-1]]
    else:
        if isinstance(node, dict):
            node[rng.choice(["x", "critical", "sha256", "size", "dtype", "shape", "kind",
                             "id", "structure_ref", "members", "qvf_version"])] = _evil_value(rng)
        elif isinstance(node, list):
            node.insert(rng.randrange(len(node) + 1), _evil_value(rng))
    _write_json(e, doc, rng)


def op_member_field(rng, a):
    """Edit dtype, shape, size, sha256, path or format of one member spec."""
    e, doc = _manifest(a)
    ms = _members(doc)
    if not ms:
        return
    s, role, m = rng.choice(ms)
    field = rng.choice(["dtype", "shape", "size", "sha256", "path", "format"])
    if field == "dtype":
        m["dtype"] = rng.choice(["int8", "uint8", "int16", "uint16", "int32", "uint32",
                                 "int64", "uint64", "float32", "float64"] + EVIL_STRINGS[:10]
                                + [None, 8, ["float64"]])
    elif field == "shape":
        m["shape"] = rng.choice([[], [0], [-1], [1.5], ["3"], [True], 3, "3", None,
                                 [2 ** 31, 2 ** 31], [2 ** 40, 2 ** 40, 2 ** 40],
                                 [2 ** 64], [10 ** 30], [2 ** 40, 0], [2 ** 62, 2 ** 62, 0], [1] * 100, [0, 2 ** 62],
                                 [rng.randint(0, 100) for _ in range(rng.randint(1, 5))]])
    elif field == "size":
        m["size"] = rng.choice(INTERESTING_INTS + [None, True, 1.0, "12", -0.0, 1e20])
    elif field == "sha256":
        h = m.get("sha256") if isinstance(m.get("sha256"), str) else "0" * 64
        m["sha256"] = rng.choice([h.upper(), h[:-1], h + "0", "", None, 0, "0" * 64,
                                  h[:-1] + ("0" if h[-1:] != "0" else "1")])
    elif field == "path":
        others = [x["path"] for _, _, x in ms if isinstance(x.get("path"), str)]
        m["path"] = rng.choice(EVIL_NAMES[:-1] + others + [None, 7, ["a"]])
    else:
        m["format"] = rng.choice(["json", "binary", "opaque", None, 1, ""])
    _write_json(e, doc, rng)


def op_payload(rng, a):
    """Mutate a JSON member's payload, keeping or refreshing its sha256."""
    e, doc = _manifest(a)
    ms = [(s, r, m) for s, r, m in _members(doc) if m.get("format") == "json"
          and isinstance(m.get("path"), str)]
    if not ms:
        return
    s, role, m = rng.choice(ms)
    target = a.find(m["path"].encode("utf-8", "surrogatepass"))
    if target is None:
        return
    raw = _payload(target)
    try:
        payload = json.loads(raw)
    except (TypeError, ValueError):
        return
    nodes = list(_nodes(payload))
    path, _ = rng.choice(nodes)
    payload = _replace(payload, path, _evil_value(rng))
    text = json.dumps(payload, allow_nan=True).encode("utf-8", "surrogatepass")
    fresh = store("x", text, deflate=target.method == 8)
    target.data, target.method, target.crc, target.usize = (
        fresh.data, fresh.method, fresh.crc, fresh.usize)
    if rng.random() < 0.8:
        m["sha256"] = hashlib.sha256(text).hexdigest()
    if "size" in m and rng.random() < 0.8:
        m["size"] = len(text)
    _write_json(e, doc, rng)


def op_binary_payload(rng, a):
    """Resize or rewrite a binary member, refreshing sha256 and size."""
    e, doc = _manifest(a)
    ms = [(s, r, m) for s, r, m in _members(doc) if m.get("format") == "binary"
          and isinstance(m.get("path"), str)]
    if not ms:
        return
    s, role, m = rng.choice(ms)
    target = a.find(m["path"].encode("utf-8", "surrogatepass"))
    if target is None:
        return
    raw = _payload(target) or b""
    r = rng.random()
    if r < 0.3:
        raw = raw[: rng.randrange(len(raw) + 1)]
    elif r < 0.5:
        raw = raw + bytes(rng.randrange(256) for _ in range(rng.randint(1, 16)))
    elif r < 0.75:
        # NaN, infinities and denormals in a float64 member.
        vals = [float("nan"), float("inf"), -float("inf"), 5e-324, 1e308, -0.0]
        raw = b"".join(struct.pack("<d", rng.choice(vals)) for _ in range(len(raw) // 8))
    else:
        raw = bytes(rng.randrange(256) for _ in range(len(raw)))
    fresh = store("x", raw, deflate=target.method == 8)
    target.data, target.method, target.crc, target.usize = (
        fresh.data, fresh.method, fresh.crc, fresh.usize)
    m["sha256"] = hashlib.sha256(raw).hexdigest()
    if "size" in m:
        m["size"] = len(raw)
    _write_json(e, doc, rng)


def op_json_text(rng, a):
    """Rewrite manifest.json at the text level: what json.dumps cannot produce."""
    e = a.find(b"manifest.json")
    if e is None:
        return
    raw = _payload(e)
    if raw is None:
        return
    r = rng.random()
    if r < 0.1:
        raw = raw.replace(b"1", b"NaN", 1)
    elif r < 0.2:
        raw = raw.replace(b": 1", b": Infinity", 1) if b": 1" in raw else b"-Infinity"
    elif r < 0.3:
        raw = raw.replace(b"1", b"1" + b"0" * rng.choice([400, 5000]), 1)
    elif r < 0.4:
        raw = raw.decode("utf-8", "replace").encode(rng.choice(["utf-16", "utf-16-le",
                                                                "utf-32"]))
    elif r < 0.45:
        raw = b"\xef\xbb\xbf" + raw
    elif r < 0.55:
        raw = raw + rng.choice([b" x", b"}", b"\x00", b"\n\n"])
    elif r < 0.65:
        i = rng.randrange(len(raw) + 1)
        raw = raw[:i] + rng.choice([b"\xff", b"\xc0\xaf", b"\xed\xa0\x80", b"\x80"]) + raw[i:]
    elif r < 0.75:
        d = rng.choice([63, 64, 65, 1000, 100000])
        raw = b"[" * d + b"]" * d
    elif r < 0.85:
        d = rng.choice([63, 64, 65, 3000])
        raw = raw[:-1] + b', "deep": ' + b"[" * d + b"]" * d + b"}" if raw.endswith(b"}") \
            else raw
    elif r < 0.9:
        raw = raw.replace(b'"', b'"\\ud800', 1)
    else:
        i = rng.randrange(len(raw) + 1)
        raw = raw[:i] + raw[i + rng.randint(1, 8):]
    fresh = store("x", raw, deflate=e.method == 8)
    e.data, e.method, e.crc, e.usize = fresh.data, fresh.method, fresh.crc, fresh.usize
    e.local.clear()
    e.central.clear()


GROUP_JSON = [op_json_node, op_member_field, op_payload, op_binary_payload, op_json_text]
OPS = {f.__name__[3:]: f for f in GROUP_ZIP + GROUP_JSON}
ZIP_OPS = [f.__name__[3:] for f in GROUP_ZIP]
JSON_OPS = [f.__name__[3:] for f in GROUP_JSON]


# --- byte level ----------------------------------------------------------------

def byte_mutate(rng: random.Random, raw: bytes) -> bytes:
    d = bytearray(raw)
    r = rng.random()
    if not d:
        return bytes(rng.randrange(256) for _ in range(rng.randint(1, 64)))
    if r < 0.4:
        for _ in range(rng.randint(1, 8)):
            d[rng.randrange(len(d))] = rng.randrange(256)
    elif r < 0.6:
        del d[rng.randrange(len(d)):]
    elif r < 0.8:
        # Overwrite a little-endian integer at a header-looking position.
        sig = rng.choice([b"PK\x03\x04", b"PK\x01\x02", b"PK\x05\x06", b"PK\x06\x06"])
        pos = [i for i in range(len(d) - 3) if d[i:i + 4] == sig]
        if pos:
            p = rng.choice(pos) + rng.randrange(4, 46)
            v = rng.choice(INTERESTING_SIZES) & SAT32
            if p + 4 <= len(d):
                d[p:p + 4] = struct.pack("<I", v)
    else:
        i = rng.randrange(len(d))
        j = min(len(d), i + rng.randint(1, 256))
        d[i:i] = d[i:j]
    return bytes(d)


def mutate(raw: bytes, ops: list[tuple[str, int]], byte_seed: int | None = None) -> bytes:
    """Apply structure-aware operators to a seed archive, then optional byte noise."""
    try:
        a = parse(raw)
    except Exception:  # noqa: BLE001 - a non-ZIP seed only gets byte noise
        a = None
    if a is not None:
        for name, seed in ops:
            OPS[name](random.Random(seed), a)
        raw = a.build()
    if byte_seed is not None:
        raw = byte_mutate(random.Random(byte_seed), raw)
    return raw
