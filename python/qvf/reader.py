"""Lazy, bounded readers with typed manifests and optional NumPy arrays."""
from __future__ import annotations

import dataclasses
import io
import math
import struct
from typing import Any

from . import model, validate as _validation
from ._container import DTYPES, Limits, check_archive, checked_member, stream_member


@dataclasses.dataclass(frozen=True)
class Tensor:
    """Owning little-endian tensor with a read-only byte ``memoryview``.

    ``data`` is always bytes, including on big-endian hosts. ``numpy()`` is an
    explicit optional dependency; ``tolist()`` decodes with the standard
    library and preserves scalar and zero-extent shapes.
    """
    data: memoryview
    dtype: str
    shape: tuple[int, ...]

    def numpy(self):
        import numpy as np
        return np.frombuffer(self.data, dtype=np.dtype(self.dtype).newbyteorder("<")).reshape(self.shape)

    #: Most empty lists ``tolist()`` builds for a zero-extent shape, whose
    #: nested lists cost memory that no stored byte pays for.
    MAX_EMPTY_LISTS = 1 << 20

    def tolist(self):
        if 0 in self.shape:
            lists = math.prod(self.shape[:self.shape.index(0)])
            if lists > self.MAX_EMPTY_LISTS:
                raise ValueError(f"tolist(): shape {self.shape} holds no elements but "
                                 f"would need {lists} nested empty lists; use data")
        fmt = "<" + DTYPES[self.dtype][1]
        values = iter(x[0] for x in struct.iter_unpack(fmt, self.data))
        def unpack(shape):
            if not shape:
                return next(values)
            return [unpack(shape[1:]) for _ in range(shape[0])]
        return unpack(self.shape)


class Reader:
    """Open a path or seekable binary stream, reading only the manifest.

    JSON and tensor payloads are loaded on demand. All entry paths and central
    directory resource caps are checked at open. Keep the context open while
    using lazy kind views returned by :meth:`get`.
    """
    def __init__(self, source, *, limits=None):
        self.limits = limits or Limits()
        self._zip = _validation._open_zip(source)
        try:
            check_archive(self._zip, self.limits)
            try:
                info = self._zip.getinfo("manifest.json")
            except KeyError:
                raise ValueError("E_MANIFEST_MISSING: manifest.json") from None
            raw = self._json_bytes({"path": "manifest.json"}, verify=False)
            try:
                self.raw_manifest = _validation._parse_json(raw)
            except (ValueError, UnicodeError) as exc:
                raise ValueError(f"E_MANIFEST_JSON: {exc}") from exc
            if not isinstance(self.raw_manifest, dict):
                raise ValueError("E_MANIFEST_JSON: manifest must be an object")
            version = self.raw_manifest.get("qvf_version")
            if isinstance(version, bool) or version not in (1, 2):
                raise ValueError(f"E_VERSION: unsupported qvf_version {version!r}")
            if self.raw_manifest.get("byte_order", "little") != "little":
                raise ValueError("E_SCHEMA: byte_order must be little")
            self.manifest = model.Manifest.from_json(self.raw_manifest)
        except BaseException:
            self._zip.close()
            raise

    @classmethod
    def open(cls, path, *, limits=None):
        return cls(path, limits=limits)

    @classmethod
    def from_bytes(cls, data, *, limits=None):
        return cls(io.BytesIO(data), limits=limits)

    @property
    def typed_manifest(self):
        return self.manifest

    @property
    def sections(self):
        return self.manifest.sections or []

    def sections_of(self, kind):
        return [section for section in self.sections if section.kind == kind]

    def section(self, section_id):
        return self.manifest.section(section_id)

    def _member(self, section, role):
        if isinstance(section, str):
            section = self.section(section)
        if not isinstance(section, model.Section):
            raise KeyError("unknown section")
        member = section.member(role)
        if member is None:
            raise KeyError(f"section {section.id!r} has no member {role!r}")
        return member.to_json()

    def read_member(self, member, *, verify=True):
        if isinstance(member, model.MemberSpec):
            member = member.to_json()
        chunks = []
        stream_member(self._zip, member, self.limits, chunks.append, verify=verify)
        return b"".join(chunks)

    def _json_bytes(self, member, *, verify=True):
        info = checked_member(self._zip, member, self.limits)
        if info.file_size > self.limits.max_json_bytes:
            raise ValueError(f"E_JSON_TOO_LARGE: {info.filename!r}")
        raw = self.read_member(member, verify=verify)
        if _validation._json_depth(raw, self.limits.max_json_depth) > self.limits.max_json_depth:
            raise ValueError(f"E_JSON_TOO_DEEP: {info.filename!r}")
        return raw

    def json(self, section, role, *, verify=True):
        member = self._member(section, role)
        if member.get("format") != "json":
            raise TypeError(f"member {role!r} is not JSON")
        try:
            return _validation._parse_json(self._json_bytes(member, verify=verify))
        except (ValueError, UnicodeError) as exc:
            if str(exc).startswith("E_"):
                raise
            raise ValueError(f"E_JSON_INVALID: {member['path']!r}: {exc}") from exc

    def tensor(self, section, role, *, verify=True):
        member = self._member(section, role)
        if member.get("format") != "binary" or "dtype" not in member or "shape" not in member:
            raise TypeError(f"member {role!r} is not a tensor")
        raw = self.read_member(member, verify=verify)
        return Tensor(memoryview(raw), member["dtype"], tuple(member["shape"]))

    def tensor_into(self, section, role, sink, *, verify=True, chunk_size=1 << 20):
        """Stream little-endian byte chunks to a callback or binary file.

        Shape/size is checked first; checksum/CRC verification finishes at EOF.
        A consumer must discard the result on an exception.
        """
        member = self._member(section, role)
        if member.get("format") != "binary" or "dtype" not in member or "shape" not in member:
            raise TypeError(f"member {role!r} is not a tensor")
        callback = sink if callable(sink) else sink.write
        return stream_member(self._zip, member, self.limits, callback,
                             verify=verify, chunk_size=chunk_size)

    def get(self, section_id, kind_type=None):
        """Return a lazy kind view (or a requested class with ``from_reader``)."""
        from .kinds import from_reader
        section = self.section(section_id)
        if section is None:
            raise KeyError(section_id)
        return (kind_type.from_reader(self, section) if kind_type is not None
                else from_reader(self, section))

    def validate(self, *, level=_validation.Level.All):
        findings = _validation._Findings()
        count = _validation._validate_open_zip(self._zip, self.limits, findings, level)
        return findings.report(count)

    def close(self):
        self._zip.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
