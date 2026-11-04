#!/usr/bin/env python3
# License: Apache-2.0.
"""Check archives written by the C++ library on another host (CONTRIBUTING.md).

    python scripts/check_endian_roundtrip.py DIR [--reference REFDIR]

Every ``*.qvf`` in DIR must pass the Python validator. Archives written by
``qvf_endian_probe write`` (cpp/tests/endian_probe.cpp) are also compared
bit for bit with values this script computes itself with NumPy. With
``--reference``, every typed member must equal the same member of the
archive of the same name in REFDIR (written on a host of the other byte
order), as raw bytes and as decoded NumPy arrays. Exit status 1 on any
mismatch. Needs ``pip install -e "python/[numpy,validate]"``.
"""
import argparse
import math
import pathlib
import sys

import numpy as np

import qvf

N = 24
SHAPE = (2, 3, 4)
INTS = ["int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64"]


def u64(i):
    return (0x0102030405060708 * (i + 1) + 0x8000000000000081 * (i & 1)) % (1 << 64)


def f64(i):
    if i == 0:
        return -0.0
    if i == 1:
        return math.ldexp(1.0, -1074)
    x = math.ldexp((i + 1) * 0.1, i - 12)
    return -x if i & 1 else x


def expected(dtype):
    """The probe values of one dtype, as a little-endian NumPy array."""
    if dtype in INTS:
        width = np.dtype(dtype).itemsize
        raw = b"".join((u64(i) % (1 << (8 * width))).to_bytes(width, "little") for i in range(N))
        return np.frombuffer(raw, dtype=np.dtype(dtype).newbyteorder("<")).reshape(SHAPE)
    values = np.array([f64(i) for i in range(N)], dtype="<f8")
    if dtype == "float32":
        values = values.astype("<f4")
        values[1] = np.float32(math.ldexp(1.0, -149))
    return values.reshape(SHAPE)


def probe_expectation(section, role):
    if role == "convert_int16_to_float64":
        return expected("int16").astype("<f8")
    if role == "convert_float64_to_float32":
        return expected("float64").astype("<f4")
    name = role if section in ("stream", "legacy") else role.split("_", 1)[1]
    if section == "legacy" and name == "float32":
        return expected("float64").astype("<f4")
    return expected(name)


def same_bits(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def tensors(path):
    with qvf.Reader.open(path) as r:
        out = {}
        for s in r.manifest.sections or []:
            for role, spec in (s.members or {}).items():
                if spec.dtype is not None:
                    t = r.tensor(s, role)
                    out[(s.id, role)] = (bytes(t.data), t.numpy())
        return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dir", type=pathlib.Path)
    ap.add_argument("--reference", type=pathlib.Path)
    args = ap.parse_args()
    paths = sorted(args.dir.glob("*.qvf"))
    if not paths:
        sys.exit(f"no .qvf archives in {args.dir}")
    bad = 0
    for path in paths:
        report = qvf.validate_qvf(str(path))
        problems = list(report["errors"])
        got = tensors(path)
        if path.stem.startswith("endian_probe"):
            for (sid, role), (_, arr) in got.items():
                if not same_bits(arr, probe_expectation(sid, role)):
                    problems.append(f"{sid}/{role}: values differ from the expected bits")
        if args.reference:
            ref = tensors(args.reference / path.name)
            if ref.keys() != got.keys():
                problems.append("typed members differ from the reference archive")
            for key in ref.keys() & got.keys():
                if ref[key][0] != got[key][0] or not same_bits(ref[key][1], got[key][1]):
                    problems.append(f"{key[0]}/{key[1]}: differs from the reference archive")
        status = "ok" if not problems else "FAIL"
        print(f"{status:4} {path.name}: {len(got)} typed members, "
              f"{len(report['warnings'])} warnings")
        for p in problems:
            print(f"     {p}")
        bad += bool(problems)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
