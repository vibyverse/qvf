"""An editable model of a ZIP archive for structure-aware fuzzing.

`parse` turns archive bytes into `Archive` (one `Entry` per central-directory
record, with the compressed bytes sliced from the file), and `Archive.build`
writes it back. Every header field can be set separately for the local and
the central copy, so a mutation can make the two disagree, saturate a field
and add a ZIP64 extra field, or add a data descriptor. Standard library only.
"""
from __future__ import annotations

import dataclasses
import struct
import zipfile
import zlib
from typing import Optional

SIG_LOCAL = 0x04034B50
SIG_CENTRAL = 0x02014B50
SIG_EOCD = 0x06054B50
SIG_EOCD64 = 0x06064B50
SIG_LOCATOR64 = 0x07064B50
SIG_DESCRIPTOR = 0x08074B50
SAT16 = 0xFFFF
SAT32 = 0xFFFFFFFF


@dataclasses.dataclass
class Entry:
    name: bytes
    data: bytes            # the stored (compressed) bytes
    method: int = 0
    flags: int = 0
    crc: int = 0
    usize: int = 0
    # Field overrides for one copy only: keys of the header fields
    # (method, flags, crc, csize, usize, name, offset, version).
    local: dict = dataclasses.field(default_factory=dict)
    central: dict = dataclasses.field(default_factory=dict)
    local_extra: bytes = b""
    central_extra: bytes = b""
    comment: bytes = b""
    # ZIP64: saturate the 32-bit fields named here and add a 0x0001 extra
    # field carrying `zip64_values` (default: the true values).
    zip64: tuple = ()
    zip64_values: dict = dataclasses.field(default_factory=dict)
    zip64_truncate: int = 0     # drop this many bytes from the ZIP64 field
    descriptor: Optional[str] = None   # None, "sig32", "nosig32", "sig64"
    descriptor_values: dict = dataclasses.field(default_factory=dict)
    zero_local_sizes: bool = False


@dataclasses.dataclass
class Archive:
    entries: list
    comment: bytes = b""
    prefix: bytes = b""          # bytes before the first local header
    suffix: bytes = b""          # bytes after the end record
    zip64_end: bool = False      # write ZIP64 end records
    end: dict = dataclasses.field(default_factory=dict)      # EOCD overrides
    end64: dict = dataclasses.field(default_factory=dict)    # ZIP64 EOCD overrides

    def build(self) -> bytes:
        out = bytearray(self.prefix)
        central = bytearray()
        for e in self.entries:
            offset = len(out)
            out += _local(e)
            out += e.data
            out += _descriptor(e)
            central += _central(e, offset)
        cd_offset = len(out)
        out += central
        n = len(self.entries)
        count = self.end.get("count", n)
        cd_size = self.end.get("cd_size", len(central))
        cd_off = self.end.get("cd_offset", cd_offset)
        if self.zip64_end:
            record = len(out)
            out += struct.pack(
                "<IQHHIIQQQQ", SIG_EOCD64, self.end64.get("record_size", 44),
                45, 45, self.end64.get("disk", 0), self.end64.get("cd_disk", 0),
                self.end64.get("count", n), self.end64.get("count", n),
                self.end64.get("cd_size", len(central)),
                self.end64.get("cd_offset", cd_offset))
            out += struct.pack("<IIQI", SIG_LOCATOR64, self.end64.get("locator_disk", 0),
                               self.end64.get("locator_offset", record),
                               self.end64.get("disks", 1))
            count, cd_size, cd_off = SAT16, SAT32, SAT32
        out += struct.pack(
            "<IHHHHIIH", SIG_EOCD, self.end.get("disk", 0) & SAT16,
            self.end.get("cd_disk", 0) & SAT16, count & SAT16, count & SAT16,
            cd_size & SAT32, cd_off & SAT32,
            self.end.get("comment_length", len(self.comment)) & SAT16)
        out += self.comment
        out += self.suffix
        return bytes(out)

    def find(self, name: bytes) -> Optional[Entry]:
        for e in self.entries:
            if e.name == name:
                return e
        return None


def _zip64_extra(e: Entry, offset: int) -> bytes:
    if not e.zip64:
        return b""
    true = {"usize": e.usize, "csize": len(e.data), "offset": offset}
    body = b""
    for key in ("usize", "csize", "offset"):
        if key in e.zip64:
            body += struct.pack("<Q", e.zip64_values.get(key, true[key]) & (2 ** 64 - 1))
    if e.zip64_truncate:
        body = body[: max(0, len(body) - e.zip64_truncate)]
    return struct.pack("<HH", 1, e.zip64_values.get("length", len(body)) & SAT16) + body


def _local(e: Entry) -> bytes:
    c = e.local
    usize = c.get("usize", e.usize)
    csize = c.get("csize", len(e.data))
    crc = c.get("crc", e.crc)
    if e.zero_local_sizes:
        usize = csize = crc = 0
    extra = e.local_extra
    if e.zip64:
        # The local ZIP64 field always carries both sizes (APPNOTE 4.5.3).
        v = e.zip64_values
        body = struct.pack("<QQ", v.get("usize", e.usize) & (2 ** 64 - 1),
                           v.get("csize", len(e.data)) & (2 ** 64 - 1))
        if e.zip64_truncate:
            body = body[: max(0, len(body) - e.zip64_truncate)]
        extra = struct.pack("<HH", 1, len(body)) + body + extra
        usize = csize = SAT32
    name = c.get("name", e.name)
    return struct.pack(
        "<IHHHHHIIIHH", SIG_LOCAL, c.get("version", 20) & SAT16,
        c.get("flags", e.flags) & SAT16, c.get("method", e.method) & SAT16, 0, 0x21,
        crc & SAT32, csize & SAT32, usize & SAT32, len(name) & SAT16,
        c.get("extra_length", len(extra)) & SAT16) + name + extra


def _central(e: Entry, offset: int) -> bytes:
    c = e.central
    offset = c.get("offset", offset)
    usize = c.get("usize", e.usize)
    csize = c.get("csize", len(e.data))
    extra = _zip64_extra(e, offset) + e.central_extra
    if "usize" in e.zip64:
        usize = SAT32
    if "csize" in e.zip64:
        csize = SAT32
    if "offset" in e.zip64:
        offset = SAT32
    name = c.get("name", e.name)
    return struct.pack(
        "<IHHHHHHIIIHHHHHII", SIG_CENTRAL, 45 if e.zip64 else 20,
        c.get("version", 45 if e.zip64 else 20) & SAT16, c.get("flags", e.flags) & SAT16,
        c.get("method", e.method) & SAT16, 0, 0x21, c.get("crc", e.crc) & SAT32,
        csize & SAT32, usize & SAT32, len(name) & SAT16,
        c.get("extra_length", len(extra)) & SAT16, len(e.comment) & SAT16, 0, 0, 0,
        offset & SAT32) + name + extra + e.comment


def _descriptor(e: Entry) -> bytes:
    if e.descriptor is None:
        return b""
    v = e.descriptor_values
    crc = v.get("crc", e.crc) & SAT32
    csize = v.get("csize", len(e.data))
    usize = v.get("usize", e.usize)
    if e.descriptor == "sig64":
        return struct.pack("<IIQQ", SIG_DESCRIPTOR, crc, csize & (2 ** 64 - 1),
                           usize & (2 ** 64 - 1))
    body = struct.pack("<III", crc, csize & SAT32, usize & SAT32)
    return body if e.descriptor == "nosig32" else struct.pack("<I", SIG_DESCRIPTOR) + body


def parse(raw: bytes) -> Archive:
    """The model of a well-formed archive (raises for anything zipfile rejects)."""
    import io
    entries = []
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        for info in zf.infolist():
            off = info.header_offset
            nlen, xlen = struct.unpack_from("<HH", raw, off + 26)
            start = off + 30 + nlen + xlen
            entries.append(Entry(
                name=info.orig_filename.encode("utf-8") if info.flag_bits & 0x800
                else info.orig_filename.encode("cp437", "replace"),
                data=raw[start:start + info.compress_size], method=info.compress_type,
                flags=info.flag_bits & ~0x8, crc=info.CRC, usize=info.file_size))
    return Archive(entries=entries)


def store(name: str, payload: bytes, *, deflate: bool = True) -> Entry:
    """A fresh, valid entry for `payload`."""
    if deflate:
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        data = c.compress(payload) + c.flush()
        method = 8
    else:
        data, method = payload, 0
    return Entry(name=name.encode("utf-8"), data=data, method=method,
                 flags=0x800 if not name.isascii() else 0,
                 crc=zlib.crc32(payload) & SAT32, usize=len(payload))
