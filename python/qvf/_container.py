"""Bounded ZIP and tensor primitives used by the public facades (stdlib only)."""
from __future__ import annotations

import hashlib
import math
import zipfile
from . import validate as _validation

Limits = _validation.Limits
DTYPES = {"int8": (1, "b"), "uint8": (1, "B"), "int16": (2, "h"),
          "uint16": (2, "H"), "int32": (4, "i"), "uint32": (4, "I"),
          "int64": (8, "q"), "uint64": (8, "Q"), "float32": (4, "f"),
          "float64": (8, "d")}


def tensor_size(dtype, shape):
    if not isinstance(dtype, str) or dtype not in DTYPES:
        raise ValueError(f"E_DTYPE: unsupported dtype {dtype!r}")
    try:
        shape = tuple(shape)
    except TypeError as exc:
        raise ValueError("E_BINARY_SIZE: shape must be an iterable of dimensions") from exc
    if any(not isinstance(n, int) or isinstance(n, bool) or n < 0 for n in shape):
        raise ValueError("E_BINARY_SIZE: shape must contain non-negative integers")
    return DTYPES[dtype][0] * math.prod(shape)


def check_archive(zf, limits):
    """Check central-directory caps before decompressing even the manifest."""
    duplicates = _validation._duplicate_entries(zf)
    if duplicates:
        raise ValueError(f"E_ZIP_DUPLICATE_ENTRY: {duplicates!r}")
    total = 0
    for info in zf.infolist():
        why = _validation._path_problem(info.orig_filename)
        if why:
            raise ValueError(f"E_MEMBER_PATH: {info.orig_filename!r}: {why}")
        if info.file_size > limits.max_member_bytes:
            raise ValueError(f"E_MEMBER_TOO_LARGE: {info.filename!r}")
        total += info.file_size
        if (limits.max_ratio > 0 and info.compress_type != zipfile.ZIP_STORED
                and info.file_size >= limits.ratio_exempt_bytes
                and (info.compress_size == 0 or
                     info.file_size > limits.max_ratio * info.compress_size)):
            raise ValueError(f"E_COMPRESSION_RATIO: {info.filename!r}")
    if total > limits.max_total_bytes:
        raise ValueError("E_TOTAL_TOO_LARGE: archive exceeds total expansion cap")


def checked_member(zf, member, limits):
    path = member.get("path")
    if not isinstance(path, str) or not path or _validation._path_problem(path):
        raise ValueError("E_MEMBER_PATH: illegal or absent member path")
    try:
        info = zf.getinfo(path)
    except KeyError:
        raise ValueError(f"E_MEMBER_MISSING: {path!r}") from None
    if info.file_size > limits.max_member_bytes:
        raise ValueError(f"E_MEMBER_TOO_LARGE: {path!r}")
    problem = _validation._size_problem(member, info)
    if problem:
        raise ValueError(f"{problem[0]}: {problem[1]}")
    if "dtype" in member:
        expected = tensor_size(member["dtype"], member.get("shape", ()))
        if expected != info.file_size:
            raise ValueError(f"E_BINARY_SIZE: {path!r}: shape expects {expected} bytes, "
                             f"entry holds {info.file_size}")
    return info


def iter_member(zf, member, limits, *, verify=True, chunk_size=1 << 20):
    """Deliver bounded byte chunks; verification completes after the last chunk.

    Exhaust the iterator to complete verification, and discard consumed data
    if it raises. A partially consumed stream is not authenticated.
    """
    if not isinstance(chunk_size, int) or chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    info = checked_member(zf, member, limits)
    digest = hashlib.sha256()
    # Bounded by the declared size; a short, long or damaged stream and a
    # CRC-32 mismatch raise E_ZIP_UNREADABLE (spec § 6.1).
    for chunk in _validation._iter_entry(zf, info, chunk_size):
        digest.update(chunk)
        yield memoryview(chunk)
    if verify and "sha256" in member and digest.hexdigest() != member["sha256"]:
        raise ValueError(f"E_SHA256_MISMATCH: {info.filename!r}")


def stream_member(zf, member, limits, sink, *, verify=True, chunk_size=1 << 20):
    """Stream bounded chunks to a callback; discard output on verification error."""
    size = 0
    for chunk in iter_member(zf, member, limits, verify=verify, chunk_size=chunk_size):
        sink(chunk)
        size += len(chunk)
    return size
