#!/usr/bin/env python3
"""Generate the QVF conformance golden corpus.

Writes a set of fixed reference ``.qvf`` archives under ``corpus/``, each paired
with a ``.expected.json`` describing the values a conforming consumer must decode
from it. The archives are **frozen test vectors**: they define the format, not
the writer, so they are committed and only regenerated when the corpus is
intentionally extended (a producer proves conformance by *decoding to the same
expected values*, not by reproducing bytes).

Run:  python generate_corpus.py

This uses the reference Python writer as the source of truth. The runner
(``run_conformance.py``) is what every consumer implements to self-certify.
"""

from __future__ import annotations

import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = os.path.join(HERE, "corpus")
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "python")))
from qvf_writer import (  # noqa: E402
    BOHR_PER_ANGSTROM, QvfWriter as _QvfWriter, primitive_norm)




def QvfWriter(*args, **kwargs):  # noqa: N802 - stands in for the class
    """The reference writer without automatic metadata (spec § 3.5): no time
    stamp and no random uuid, so a regenerated archive has the same manifest.
    An archive that wants metadata sets it explicitly with set_metadata."""
    kwargs.setdefault("auto_metadata", False)
    return _QvfWriter(*args, **kwargs)


_PROGRAM = "qvf-conformance"
_VERSION = "0.1.0"


def _write(name, writer, checks, description):
    """Write one archive + its expected.json."""
    path = os.path.join(CORPUS, f"{name}.qvf")
    writer.write(path)
    manifest = writer.build_manifest()
    expected = {
        "archive": f"{name}.qvf",
        "description": description,
        "source_program": _PROGRAM,
        "kinds": sorted({s["kind"] for s in manifest["sections"]}),
        "checks": checks,
    }
    with open(os.path.join(CORPUS, f"{name}.expected.json"), "w",
              encoding="utf-8") as fh:
        json.dump(expected, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"  {name}.qvf  ({len(manifest['sections'])} sections, "
          f"{len(checks)} checks)")


# =============================================================================
# WP-14: conformance corpus v2. Archives for the canonical kinds the corpus did
# not exercise, and the edge cases (unrestricted and Cartesian wavefunctions, a
# 1-D wire, a variable-cell reaction path, ZIP64 and data-descriptor containers,
# a legacy qvf_version 2 archive). The code lives here, apart from build(), so
# that everything above this banner stays exactly what produced the frozen
# vectors; build() calls the two functions at its end.
# =============================================================================
import io  # noqa: E402
import math  # noqa: E402
import struct  # noqa: E402
import zipfile  # noqa: E402
import zlib  # noqa: E402

from qvf_writer import BOHR_PER_ANGSTROM, Lattice  # noqa: E402

_EPOCH = (1980, 1, 1, 0, 0, 0)       # fixed ZIP timestamp: no wall-clock bytes


def _legacy_dir() -> str:
    """``legacy/``, beside the corpus directory in use."""
    return os.path.join(os.path.dirname(os.path.abspath(CORPUS)), "legacy")


# --- container helpers: re-pack an archive the reference writer produced -----

def _entries_of(zip_bytes):
    """``(name, data, compress_type)`` of every entry, in archive order."""
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        return [(i.filename, zf.read(i), i.compress_type)
                for i in zf.infolist()]


def _zip_info(name, method):
    info = zipfile.ZipInfo(name, date_time=_EPOCH)
    info.compress_type = method
    info.create_system = 3
    info.external_attr = 0o644 << 16
    return info


def _plain_zip(entries):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, data, method in entries:
            zf.writestr(_zip_info(name, method), data)
    return buf.getvalue()


class _StreamSink(io.RawIOBase):
    """A write-only stream with no ``tell`` or ``seek``. ``zipfile`` cannot go
    back to patch a local header on such a stream, so it sets the
    data-descriptor flag (general-purpose bit 3) and writes the CRC and the
    sizes after each member's data."""

    def __init__(self):
        super().__init__()
        self._parts = []

    def writable(self):
        return True

    def write(self, data):
        self._parts.append(bytes(data))
        return len(data)

    def getvalue(self):
        return b"".join(self._parts)


def _data_descriptor_zip(entries):
    """Every member streamed: local headers carry bit 3 with zero CRC and zero
    sizes, a data descriptor follows the data, and only the central directory
    knows the sizes. All members are DEFLATEd (a streaming reader cannot find
    the end of a STOREd member that has a descriptor)."""
    sink = _StreamSink()
    with zipfile.ZipFile(sink, "w") as zf:
        for name, data, _method in entries:
            zf.writestr(_zip_info(name, zipfile.ZIP_DEFLATED), data)
    return sink.getvalue()


def _zip64_zip(entries):
    """An archive with every ZIP64 structure forced although every member is
    tiny: a ZIP64 extra field in each local header, a ZIP64 extra field in each
    central-directory record with all three sizes and the offset saturated to
    0xFFFFFFFF, a ZIP64 end-of-central-directory record and locator, and an
    end record whose count, size and offset are saturated too.

    ``zipfile`` cannot produce this for small members: ``force_zip64`` writes
    the extra field into the local header only, and the central directory and
    end record stay in their 32-bit forms until a size or offset really exceeds
    the 32-bit limit. So the archive is serialized here."""
    out = bytearray()
    central = bytearray()
    dos_date = (1 << 5) | 1                 # 1980-01-01, 00:00:00
    for name, data, method in entries:
        raw_name = name.encode("utf-8")
        crc = zlib.crc32(data) & 0xFFFFFFFF
        if method == zipfile.ZIP_STORED:
            body = data
        else:
            deflater = zlib.compressobj(9, zlib.DEFLATED, -15)
            body = deflater.compress(data) + deflater.flush()
        offset = len(out)
        # Local header: both sizes live in the extra field (APPNOTE 4.5.3).
        extra = struct.pack("<HHQQ", 0x0001, 16, len(data), len(body))
        out += struct.pack("<IHHHHHIIIHH", 0x04034B50, 45, 0, method, 0,
                           dos_date, crc, 0xFFFFFFFF, 0xFFFFFFFF,
                           len(raw_name), len(extra))
        out += raw_name + extra + body
        # Central record: uncompressed size, compressed size, local offset.
        c_extra = struct.pack("<HHQQQ", 0x0001, 24, len(data), len(body),
                              offset)
        central += struct.pack("<IHHHHHHIIIHHHHHII", 0x02014B50,
                               (3 << 8) | 45, 45, 0, method, 0, dos_date, crc,
                               0xFFFFFFFF, 0xFFFFFFFF, len(raw_name),
                               len(c_extra), 0, 0, 0, 0o100644 << 16,
                               0xFFFFFFFF)
        central += raw_name + c_extra
    cd_offset = len(out)
    out += central
    end64_offset = len(out)
    n = len(entries)
    out += struct.pack("<IQHHIIQQQQ", 0x06064B50, 44, (3 << 8) | 45, 45, 0, 0,
                       n, n, len(central), cd_offset)
    out += struct.pack("<IIQI", 0x07064B50, 0, end64_offset, 1)
    out += struct.pack("<IHHHHIIH", 0x06054B50, 0, 0, 0xFFFF, 0xFFFF,
                       0xFFFFFFFF, 0xFFFFFFFF, 0)
    return bytes(out)


class _Rebuilt:
    """Adapts a finished ``QvfWriter`` to ``_write``, re-packing the ZIP the
    writer produced (``transform`` maps its entries to archive bytes). The
    member bytes, and so every sha256 in the manifest, are unchanged."""

    def __init__(self, writer, transform):
        self._writer = writer
        self._transform = transform

    def write(self, path):
        data = self._transform(_entries_of(self._writer.to_bytes()))
        with open(path, "wb") as fh:
            fh.write(data)
        return path

    def build_manifest(self):
        return self._writer.build_manifest()


def _write_legacy(name, data, description, expect):
    """Write one archive of ``legacy/`` and its negative-format expected.json."""
    directory = _legacy_dir()
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, f"{name}.qvf"), "wb") as fh:
        fh.write(data)
    expected = {"archive": f"{name}.qvf", "description": description,
                "expect": expect}
    with open(os.path.join(directory, f"{name}.expected.json"), "w",
              encoding="utf-8") as fh:
        json.dump(expected, fh, indent=2)
        fh.write("\n")
    print(f"  legacy/{name}.qvf  ({len(data)} bytes) -> "
          + ", ".join(e["code"] for e in expect))


# --- Appendix A.1 reference: self-overlap of a contracted Cartesian function --

def _odd_double_factorial(n):
    """(2n - 1)!!, with (-1)!! = 1."""
    out = 1
    for k in range(1, 2 * n, 2):
        out *= k
    return out


def _cartesian_self_overlap(exponents, coefficients, powers):
    """<chi|chi> for chi = sum_i c_i N_i x^lx y^ly z^lz exp(-a_i r^2) with the
    format's N_i (spec Appendix A.1): one factor per primitive from the total
    l = lx + ly + lz, not a per-component unit norm.

    Uses  int x^(2n) exp(-p x^2) dx = (2n-1)!! / (2p)^n * sqrt(pi / p)  per
    axis, so a primitive pair contributes
    (2lx-1)!! (2ly-1)!! (2lz-1)!! / (2p)^l * (pi / p)^(3/2)."""
    lx, ly, lz = powers
    l = lx + ly + lz
    angular = (_odd_double_factorial(lx) * _odd_double_factorial(ly)
               * _odd_double_factorial(lz))
    total = 0.0
    for a_i, c_i in zip(exponents, coefficients):
        for a_j, c_j in zip(exponents, coefficients):
            p = a_i + a_j
            total += (c_i * primitive_norm(a_i, l) * c_j * primitive_norm(a_j, l)
                      * angular / (2.0 * p) ** l * (math.pi / p) ** 1.5)
    return total


def _axial_normalized(exponents, coefficients, l):
    """``coefficients`` scaled so that the axial component ``(l, 0, 0)`` of the
    contracted shell has self-overlap 1 (the usual contraction normalization)."""
    scale = 1.0 / math.sqrt(_cartesian_self_overlap(exponents, coefficients,
                                                    (l, 0, 0)))
    return [c * scale for c in coefficients]


# --- volume helpers ----------------------------------------------------------

def _volume_grid(shape, step):
    """Origin, voxel vectors and the three coordinate arrays (bohr) of a
    point-centred grid of ``shape`` points centred on the origin."""
    origin = [-(n - 1) * step / 2.0 for n in shape]
    axes = [origin[a] + step * np.arange(shape[a]) for a in range(3)]
    coords = np.meshgrid(*axes, indexing="ij")
    voxels = [[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, step]]
    return origin, voxels, coords


def _r2(coords, center):
    return sum((c - x0) ** 2 for c, x0 in zip(coords, center))


def _volume_checks(sid, data, shape, origin, step, indices, tol=1e-6):
    """Grid checks (shape on every axis, origin, a voxel step) plus the stored
    value at each of ``indices``."""
    checks = [{"section": sid, "member": "grid",
               "json_pointer": ["shape", a], "equals": shape[a]}
              for a in range(3)]
    checks.append({"section": sid, "member": "grid",
                   "json_pointer": ["origin", 0], "equals": origin[0],
                   "tol": 1e-12})
    checks.append({"section": sid, "member": "grid",
                   "json_pointer": ["voxel_vectors", 1, 1], "equals": step,
                   "tol": 1e-12})
    for idx in indices:
        checks.append({"section": sid, "member": "data",
                       "binary_index": list(idx), "equals": float(data[idx]),
                       "tol": tol})
    return checks


def _wp14_coverage_archives():
    """One archive for each canonical kind the corpus did not yet exercise:
    ``bonds``, ``volume.orbital`` / ``spin`` / ``elf`` / ``potential`` / ``rdg``
    / ``generic`` / ``difference``, ``basis.ao``, ``dos.coop`` and
    ``citations``."""
    ang = BOHR_PER_ANGSTROM

    # bonds: the standalone connectivity section (the water archive carries its
    # bonds as a member of the structure section instead).
    w = QvfWriter(_PROGRAM, _VERSION, "formaldehyde-bonds")
    w.add_structure(
        [{"symbol": "C", "position": [0.0, 0.0, 0.0], "atomic_number": 6},
         {"symbol": "O", "position": [0.0, 0.0, 1.205], "atomic_number": 8},
         {"symbol": "H", "position": [0.943, 0.0, -0.587], "atomic_number": 1},
         {"symbol": "H", "position": [-0.943, 0.0, -0.587],
          "atomic_number": 1}])
    w.add_bonds([{"i": 0, "j": 1, "order": 2.0}, {"i": 0, "j": 2, "order": 1.0},
                 {"i": 0, "j": 3}])
    _write("bonds", w, [
        {"section": "bonds", "member": "bonds",
         "json_pointer": ["pairs", 0, "order"], "equals": 2.0, "tol": 1e-12},
        {"section": "bonds", "member": "bonds",
         "json_pointer": ["pairs", 1, "j"], "equals": 2},
        # A pair without `order` means a single bond (spec § 5.7).
        {"section": "bonds", "member": "bonds",
         "json_pointer": ["pairs", 2], "equals": {"i": 0, "j": 3}},
        {"section": "bonds", "member": "bonds",
         "json_pointer": ["pairs"], "equals": [
             {"i": 0, "j": 1, "order": 2.0}, {"i": 0, "j": 2, "order": 1.0},
             {"i": 0, "j": 3}]},
    ], "Standalone bonds section ($defs/BondsPayload) for formaldehyde: a "
       "double C=O bond, a single C-H bond, and a pair that omits `order`, "
       "which means a single bond.")

    # volume.orbital: signed amplitude (an antibonding sigma_u orbital of H2).
    shape, step = (5, 4, 6), 0.8
    origin, voxels, xyz = _volume_grid(shape, step)
    h2_a, h2_b = (0.0, 0.0, -0.7), (0.0, 0.0, 0.7)
    orb = (np.exp(-0.6 * _r2(xyz, h2_a))
           - np.exp(-0.6 * _r2(xyz, h2_b))).astype(np.float32)
    assert orb.min() < 0.0 < orb.max()
    w = QvfWriter(_PROGRAM, _VERSION, "h2-orbital")
    w.add_structure(
        [{"symbol": "H", "position": [0.0, 0.0, z / ang], "atomic_number": 1}
         for z in (-0.7, 0.7)])
    w.add_orbital(orb, origin, voxels, label="LUMO (sigma_u)",
                  section_id="lumo")
    _write("volume_orbital", w, [
        {"manifest_pointer": ["sections", 1, "label"],
         "equals": "LUMO (sigma_u)"},
        {"manifest_pointer": ["sections", 1, "members", "data", "dtype"],
         "equals": "float32"},
    ] + _volume_checks("lumo", orb, shape, origin, step,
                       [(2, 1, 0), (2, 1, 5), (4, 3, 2)]),
        "volume.orbital: a signed orbital amplitude (the data change sign "
        "across the molecule) on an anisotropic 5x4x6 float32 grid in bohr; "
        "origin and voxel_vectors in bohr.")

    # volume.spin: signed spin density, float64 (the opt-in dtype).
    shape, step = (4, 5, 5), 0.9
    origin, voxels, xyz = _volume_grid(shape, step)
    o_pos, h_pos = (0.0, 0.0, 0.0), (0.0, 0.0, 1.84)
    spin = (np.exp(-0.8 * _r2(xyz, o_pos))
            - 0.12 * np.exp(-0.9 * _r2(xyz, h_pos)))
    assert spin.min() < 0.0 < spin.max()
    w = QvfWriter(_PROGRAM, _VERSION, "oh-spin")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.0], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.0, 1.84 / ang],
          "atomic_number": 1}])
    w.add_volume("volume.spin", spin, origin, voxels, dtype="float64",
                 label="Spin density", section_id="spin")
    _write("volume_spin", w, [
        {"manifest_pointer": ["sections", 1, "members", "data", "dtype"],
         "equals": "float64"},
    ] + _volume_checks("spin", spin, shape, origin, step,
                       [(1, 2, 2), (1, 2, 4), (3, 4, 0)], tol=1e-12),
        "volume.spin: a signed spin density of the OH radical (positive on O, "
        "negative on H) stored as float64 [nx, ny, nz] on a 4x5x5 grid in "
        "bohr.")

    # volume.elf: dimensionless, within [0, 1].
    shape, step = (5, 5, 5), 0.7
    origin, voxels, xyz = _volume_grid(shape, step)
    elf = (0.05 + 0.9 * np.exp(-0.15 * _r2(xyz, (0.0, 0.0, 0.0)))
           ).astype(np.float32)
    assert 0.0 <= elf.min() and elf.max() <= 1.0
    w = QvfWriter(_PROGRAM, _VERSION, "h2o-elf")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692],
          "atomic_number": 1}])
    w.add_volume("volume.elf", elf, origin, voxels,
                 label="Electron localization function", section_id="elf")
    _write("volume_elf", w, _volume_checks(
        "elf", elf, shape, origin, step, [(2, 2, 2), (0, 1, 4), (4, 3, 1)]),
        "volume.elf: a dimensionless localization function with values in "
        "[0, 1] (float32, 5x5x5 grid in bohr).")

    # volume.potential: hartree / e, positive and negative.
    shape, step = (6, 5, 4), 0.8
    origin, voxels, xyz = _volume_grid(shape, step)
    hf_h, hf_f = (0.0, 0.0, 0.9), (0.0, 0.0, -0.9)
    esp = (0.42 / np.sqrt(_r2(xyz, hf_h) + 0.4)
           - 0.42 / np.sqrt(_r2(xyz, hf_f) + 0.4)).astype(np.float32)
    assert esp.min() < 0.0 < esp.max()
    w = QvfWriter(_PROGRAM, _VERSION, "hf-esp")
    w.add_structure(
        [{"symbol": "H", "position": [0.0, 0.0, 0.9 / ang], "atomic_number": 1},
         {"symbol": "F", "position": [0.0, 0.0, -0.9 / ang],
          "atomic_number": 9}])
    w.add_volume("volume.potential", esp, origin, voxels,
                 label="Electrostatic potential", section_id="esp")
    _write("volume_potential", w, _volume_checks(
        "esp", esp, shape, origin, step, [(1, 2, 0), (3, 2, 3), (5, 4, 3)]),
        "volume.potential: an electrostatic potential in hartree / e, "
        "positive near the hydrogen and negative near the fluorine "
        "(float32, 6x5x4 grid in bohr).")

    # volume.rdg: dimensionless, >= 0 (a van der Waals pair).
    shape, step = (5, 5, 6), 0.9
    origin, voxels, xyz = _volume_grid(shape, step)
    rdg = (0.05 + 0.9 * (1.0 - np.exp(-0.25 * _r2(xyz, (0.0, 0.0, 0.0))))
           ).astype(np.float32)
    assert rdg.min() > 0.0
    w = QvfWriter(_PROGRAM, _VERSION, "ne2-rdg")
    w.add_structure(
        [{"symbol": "Ne", "position": [0.0, 0.0, z], "atomic_number": 10}
         for z in (-1.545, 1.545)])
    w.add_volume("volume.rdg", rdg, origin, voxels,
                 label="Reduced density gradient", section_id="rdg")
    _write("volume_rdg", w, _volume_checks(
        "rdg", rdg, shape, origin, step, [(2, 2, 2), (0, 4, 5), (4, 1, 3)]),
        "volume.rdg: the dimensionless reduced density gradient of a neon "
        "pair (float32, 5x5x6 grid in bohr); small between the atoms.")

    # volume.generic: the escape hatch, with a label that says what it is.
    shape, step = (4, 6, 5), 0.8
    origin, voxels, xyz = _volume_grid(shape, step)
    lap = (4.0 * np.exp(-1.1 * _r2(xyz, (0.0, 0.0, -1.1)))
           - 3.0 * np.exp(-0.9 * _r2(xyz, (0.0, 0.0, 1.1)))).astype(np.float32)
    w = QvfWriter(_PROGRAM, _VERSION, "co-laplacian")
    w.add_structure(
        [{"symbol": "C", "position": [0.0, 0.0, -1.1 / ang],
          "atomic_number": 6},
         {"symbol": "O", "position": [0.0, 0.0, 1.1 / ang],
          "atomic_number": 8}])
    w.add_volume("volume.generic", lap, origin, voxels,
                 label="Laplacian of the electron density",
                 section_id="laplacian")
    _write("volume_generic", w, [
        {"manifest_pointer": ["sections", 1, "label"],
         "equals": "Laplacian of the electron density"},
    ] + _volume_checks("laplacian", lap, shape, origin, step,
                       [(1, 2, 1), (2, 3, 4), (3, 5, 0)]),
        "volume.generic: a scalar field that no specific kind describes, "
        "identified by its section label (float32, 4x6x5 grid in bohr).")

    # volume.difference: operand_a - operand_b, both operands in the archive.
    shape, step = (4, 4, 5), 0.8
    origin, voxels, xyz = _volume_grid(shape, step)
    rho_a = (0.30 * np.exp(-0.7 * _r2(xyz, (0.0, 0.0, -0.6)))
             + 0.20 * np.exp(-0.7 * _r2(xyz, (0.0, 0.0, 0.6)))
             ).astype(np.float32)
    rho_b = (0.25 * np.exp(-0.7 * _r2(xyz, (0.0, 0.0, -0.6)))
             + 0.25 * np.exp(-0.7 * _r2(xyz, (0.0, 0.0, 0.6)))
             ).astype(np.float32)
    diff = (rho_a - rho_b).astype(np.float32)
    assert diff.min() < 0.0 < diff.max()
    w = QvfWriter(_PROGRAM, _VERSION, "h2-difference")
    w.add_structure(
        [{"symbol": "H", "position": [0.0, 0.0, z / ang], "atomic_number": 1}
         for z in (-0.6, 0.6)])
    w.add_density(rho_a, origin, voxels, label="rho(A)", section_id="rho_a")
    w.add_density(rho_b, origin, voxels, label="rho(B)", section_id="rho_b")
    w.add_volume("volume.difference", diff, origin, voxels,
                 label="rho(A) - rho(B)", section_id="rho_diff",
                 operand_a="rho_a", operand_b="rho_b",
                 description="rho(A) - rho(B)")
    _write("volume_difference", w, [
        {"manifest_pointer": ["sections", 3, "operand_a"], "equals": "rho_a"},
        {"manifest_pointer": ["sections", 3, "operand_b"], "equals": "rho_b"},
        {"manifest_pointer": ["sections", 3, "description"],
         "equals": "rho(A) - rho(B)"},
    ] + _volume_checks("rho_diff", diff, shape, origin, step,
                       [(1, 1, 1), (2, 1, 3), (0, 3, 4)])
      + [{"section": "rho_a", "member": "data", "binary_index": [1, 1, 1],
          "equals": float(rho_a[1, 1, 1]), "tol": 1e-6},
         {"section": "rho_b", "member": "data", "binary_index": [1, 1, 1],
          "equals": float(rho_b[1, 1, 1]), "tol": 1e-6}],
        "volume.difference: data = operand_a - operand_b, with both operands "
        "(two volume.density sections on the same grid) in the archive and "
        "named by operand_a / operand_b section ids.")

    # basis.ao: Cartesian components, named by cartesian_powers (Appendix A.3).
    alpha = 1.2
    norm = primitive_norm(alpha, 2)
    shape, step = (7, 7, 7), 0.6
    origin, voxels, xyz = _volume_grid(shape, step)
    x, y, z = xyz
    r2 = _r2(xyz, (0.0, 0.0, 0.0))
    ao_xx = (norm * x * x * np.exp(-alpha * r2)).astype(np.float32)
    ao_xy = (norm * x * y * np.exp(-alpha * r2)).astype(np.float32)
    w = QvfWriter(_PROGRAM, _VERSION, "ne-cartesian-ao")
    w.add_structure([{"symbol": "Ne", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 10}])
    common = {"atom_index": 0, "atom_symbol": "Ne", "shell_index": 0,
              "primitive_index": 0, "angular_momentum": [2, 0],
              "shell_type": "d", "exponent": alpha, "coefficient": 1.0,
              "is_primitive": True, "is_contracted": False,
              "basis_label": "model-cartesian"}
    w.add_basis_ao(dict(common, ao_index=0, cartesian_powers=[2, 0, 0]),
                   ao_xx, origin, voxels, label="d_xx", section_id="ao_xx")
    w.add_basis_ao(dict(common, ao_index=1, cartesian_powers=[1, 1, 0]),
                   ao_xy, origin, voxels, label="d_xy", section_id="ao_xy")
    _write("basis_ao", w, [
        {"manifest_pointer": ["sections", 1, "ao_metadata",
                              "cartesian_powers"], "equals": [2, 0, 0]},
        {"manifest_pointer": ["sections", 2, "ao_metadata",
                              "cartesian_powers"], "equals": [1, 1, 0]},
        # The m slot of a Cartesian component is 0 by convention and means
        # nothing (Appendix A.3).
        {"manifest_pointer": ["sections", 2, "ao_metadata",
                              "angular_momentum"], "equals": [2, 0]},
        {"manifest_pointer": ["sections", 2, "ao_metadata", "ao_index"],
         "equals": 1},
        {"manifest_pointer": ["sections", 2, "ao_metadata", "exponent"],
         "equals": alpha, "tol": 1e-12},
        {"section": "ao_xy", "member": "grid",
         "json_pointer": ["shape", 2], "equals": 7},
        # N_i of the axial component, applied to d_xy as well (Appendix A.1):
        # the stored values are N * x * y * exp(-alpha r^2), with no extra
        # per-component factor.
        {"section": "ao_xy", "member": "data", "binary_index": [4, 4, 3],
         "equals": float(ao_xy[4, 4, 3]), "tol": 1e-6},
        {"section": "ao_xy", "member": "data", "binary_index": [2, 4, 3],
         "equals": -float(ao_xy[4, 4, 3]), "tol": 1e-6},
        {"section": "ao_xx", "member": "data", "binary_index": [5, 3, 3],
         "equals": float(ao_xx[5, 3, 3]), "tol": 1e-6},
        {"section": "ao_xx", "member": "data", "binary_index": [3, 5, 3],
         "equals": float(ao_xx[3, 5, 3]), "tol": 1e-6},
    ], "basis.ao for two components of a Cartesian (pure: false) d primitive, "
       "alpha = 1.2 bohr^-2: d_xx and d_xy, each named by cartesian_powers "
       "with angular_momentum [2, 0] (Appendix A.3). Data are "
       "N * x^lx y^ly z^lz exp(-alpha r^2) with the axial N of Appendix A.1.")

    # dos.coop: bonding analysis without a spin axis (the lower-rank shapes).
    e_coop = np.linspace(-4.0, 4.0, 9)
    coop = np.stack([0.30 * np.exp(-(e_coop + 1.0) ** 2)
                     - 0.15 * np.exp(-(e_coop - 1.5) ** 2),
                     0.05 * np.exp(-(e_coop + 0.5) ** 2)
                     - 0.10 * np.exp(-(e_coop - 2.0) ** 2)])     # [2, 9]
    icoop = np.array([0.18, -0.02])                              # [2]
    a_fe = 2.87
    w = QvfWriter(_PROGRAM, _VERSION, "fe-coop")
    w.add_structure(
        [{"symbol": "Fe", "position": [0.0, 0.0, 0.0], "atomic_number": 26},
         {"symbol": "Fe", "position": [a_fe / 2] * 3, "atomic_number": 26}],
        pbc=[True, True, True],
        lattice_vectors=[[a_fe, 0, 0], [0, a_fe, 0], [0, 0, a_fe]])
    w.add_dos_coop_cohp("dos.coop", e_coop, coop, icoop, {
        "method": "coop", "n_spin": 1, "fermi_energy_ev": 5.2,
        "smearing": 0.1, "smearing_type": "gaussian",
        "pairs": [{"i": 0, "j": 1, "label": "Fe1-Fe2", "distance_ang": 2.4855},
                  {"i": 0, "j": 0, "label": "Fe1-Fe1'", "distance_ang": a_fe}],
        "pair_labels": ["Fe1-Fe2", "Fe1-Fe1'"]})
    _write("dos_coop", w, [
        {"manifest_pointer": ["sections", 1, "members", "projections",
                              "shape"], "equals": [2, 9]},
        {"manifest_pointer": ["sections", 1, "members", "integrated",
                              "shape"], "equals": [2]},
        {"section": "dos_coop", "member": "meta",
         "json_pointer": ["method"], "equals": "coop"},
        {"section": "dos_coop", "member": "meta",
         "json_pointer": ["n_spin"], "equals": 1},
        {"section": "dos_coop", "member": "meta",
         "json_pointer": ["pairs", 1, "distance_ang"], "equals": a_fe,
         "tol": 1e-12},
        {"section": "dos_coop", "member": "meta",
         "json_pointer": ["pair_labels", 1], "equals": "Fe1-Fe1'"},
        {"section": "dos_coop", "member": "energies", "binary_index": [8],
         "equals": 4.0, "tol": 1e-12},
        {"section": "dos_coop", "member": "projections",
         "binary_index": [1, 3], "equals": float(coop[1, 3]), "tol": 1e-12},
        {"section": "dos_coop", "member": "projections",
         "binary_index": [0, 7], "equals": float(coop[0, 7]), "tol": 1e-12},
        {"section": "dos_coop", "member": "integrated", "binary_index": [1],
         "equals": -0.02, "tol": 1e-12},
    ], "dos.coop without a spin axis: projections [n_pairs, n_points] and "
       "integrated [n_pairs] (n_spin = 1), energies in eV referenced to E_F, "
       "pairs named by the meta payload.")

    # citations: UTF-8 BibTeX bytes in an opaque binary member.
    bib = ("@article{Hohenberg1964,\n"
           "  author  = {Hohenberg, P. and Kohn, W.},\n"
           "  title   = {Inhomogeneous Electron Gas},\n"
           "  journal = {Phys. Rev.},\n"
           "  volume  = {136},\n"
           "  pages   = {B864--B871},\n"
           "  year    = {1964}\n"
           "}\n\n"
           "@article{Schroedinger1926,\n"
           "  author  = {Schrödinger, E.},\n"
           "  title   = {Quantisierung als Eigenwertproblem},\n"
           "  journal = {Annalen der Physik},\n"
           "  volume  = {79},\n"
           "  pages   = {361--376},\n"
           "  year    = {1926}\n"
           "}\n")
    w = QvfWriter(_PROGRAM, _VERSION, "citations")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692],
          "atomic_number": 1}])
    w.add_citations(bib)
    _write("citations", w, [
        {"section": "citations", "member": "references", "decode": "utf8",
         "contains": "@article{Hohenberg1964,"},
        # A non-ASCII character: the member is UTF-8, not Latin-1.
        {"section": "citations", "member": "references", "decode": "utf8",
         "contains": "Schrödinger, E."},
        {"section": "citations", "member": "references", "decode": "utf8",
         "equals": bib},
    ], "citations: BibTeX as an opaque binary member (no dtype or shape) "
       "that decodes as UTF-8, including a non-ASCII author name.")


def _wp14_edge_case_archives():
    """Edge cases: an unrestricted and a Cartesian wavefunction, a 1-D wire, a
    variable-cell reaction path, the same logical archive in a ZIP64 and a
    data-descriptor container, and a legacy ``qvf_version: 2`` archive."""
    ang = BOHR_PER_ANGSTROM

    # Unrestricted wavefunction: a lithium atom, alpha 2 / beta 1. The alpha
    # and the beta members are separate arrays that differ.
    li_s1_e = [16.12, 2.36, 0.51]
    li_s2_e = [0.63, 0.075]
    li_p_e = [0.55]
    shells = [
        {"center": 0, "l": 0, "exponents": li_s1_e,
         "coefficients": _axial_normalized(li_s1_e, [0.15, 0.53, 0.52], 0)},
        {"center": 0, "l": 0, "exponents": li_s2_e,
         "coefficients": _axial_normalized(li_s2_e, [0.60, 0.50], 0)},
        {"center": 0, "l": 1, "exponents": li_p_e,
         "coefficients": _axial_normalized(li_p_e, [1.0], 1)}]

    def rotation(theta):
        m = np.eye(5)
        m[0, 0], m[0, 1] = math.cos(theta), math.sin(theta)
        m[1, 0], m[1, 1] = -math.sin(theta), math.cos(theta)
        return m

    theta_a, theta_b = 0.10, 0.22
    w = QvfWriter(_PROGRAM, _VERSION, "li-uhf")
    w.add_structure([{"symbol": "Li", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 3}])
    w.add_wavefunction_gto(
        shells, mo_coefficients_alpha=rotation(theta_a),
        mo_coefficients_beta=rotation(theta_b),
        energies_alpha=[-2.4777, -0.1963, 0.0640, 0.0640, 0.0640],
        occupations_alpha=[1.0, 1.0, 0.0, 0.0, 0.0],
        energies_beta=[-2.4325, -0.0970, 0.0720, 0.0720, 0.0720],
        occupations_beta=[1.0, 0.0, 0.0, 0.0, 0.0])
    _write("wavefunction_uhf", w, [
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["spin"], "equals": "unrestricted"},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["n_ao"], "equals": 5},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["alpha", "n_mo"], "equals": 5},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["alpha", "occupations", 1], "equals": 1.0},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["beta", "occupations", 1], "equals": 0.0},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["alpha", "energies", 0], "equals": -2.4777,
         "tol": 1e-12},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["beta", "energies", 0], "equals": -2.4325,
         "tol": 1e-12},
        {"manifest_pointer": ["sections", 1, "members",
                              "mo_coefficients_alpha", "shape"],
         "equals": [5, 5]},
        {"manifest_pointer": ["sections", 1, "members",
                              "mo_coefficients_beta", "shape"],
         "equals": [5, 5]},
        {"section": "wf", "member": "mo_coefficients_alpha",
         "binary_index": [0, 1], "equals": math.sin(theta_a), "tol": 1e-12},
        {"section": "wf", "member": "mo_coefficients_alpha",
         "binary_index": [1, 0], "equals": -math.sin(theta_a), "tol": 1e-12},
        {"section": "wf", "member": "mo_coefficients_beta",
         "binary_index": [0, 1], "equals": math.sin(theta_b), "tol": 1e-12},
        {"section": "wf", "member": "mo_coefficients_beta",
         "binary_index": [0, 0], "equals": math.cos(theta_b), "tol": 1e-12},
    ], "Unrestricted wavefunction.gto (spin unrestricted): a lithium atom with "
       "two alpha and one beta electron. mo_metadata carries alpha and beta "
       "objects (no top-level energies or occupations) and the coefficients "
       "are two separate members, mo_coefficients_alpha and "
       "mo_coefficients_beta. Illustrative numbers, not a converged UHF.")

    # Cartesian d shell: the Appendix A.1 normalization contract, checkable.
    s_exp, d_exp = [2.2, 0.55], [2.4, 0.7]
    s_c = _axial_normalized(s_exp, [0.7, 0.4], 0)
    d_c = _axial_normalized(d_exp, [0.55, 0.60], 2)
    powers = [(0, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0),
              (0, 1, 1), (0, 0, 2)]
    ao_labels = ["s", "d_xx", "d_xy", "d_xz", "d_yy", "d_yz", "d_zz"]
    overlap = [_cartesian_self_overlap(s_exp, s_c, powers[0])] + [
        _cartesian_self_overlap(d_exp, d_c, p) for p in powers[1:]]
    w = QvfWriter(_PROGRAM, _VERSION, "ne-cartesian")
    w.add_structure([{"symbol": "Ne", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 10}])
    w.add_wavefunction_gto(
        [{"center": 0, "l": 0, "exponents": s_exp, "coefficients": s_c},
         {"center": 0, "l": 2, "exponents": d_exp, "coefficients": d_c}],
        mo_coefficients=np.eye(7),
        energies=[-1.0, -0.4, -0.2, 0.1, 0.3, 0.5, 0.8],
        occupations=[2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], pure=False)
    w.set_extensions({"x_conformance": {"version": "1.0", "critical": False}})
    w.add_vendor_section(
        "x_conformance.overlap_diagonal",
        json_members={"payload": {
            "convention": "Appendix A.1: one N_i per primitive from the total "
                          "l; a mixed Cartesian component is not individually "
                          "unit-normalized",
            "basis_section": "wf",
            "ao_labels": ao_labels,
            "cartesian_powers": [list(p) for p in powers],
            "self_overlap": overlap}},
        section_id="overlap_check")
    _write("wavefunction_cartesian", w, [
        {"section": "wf", "member": "basis", "json_pointer": ["pure"],
         "equals": False},
        {"section": "wf", "member": "basis", "json_pointer": ["n_ao"],
         "equals": 7},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 1, "l"], "equals": 2},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 1, "pure"], "equals": False},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 1, "exponents", 0], "equals": d_exp[0],
         "tol": 1e-12},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 1, "exponents", 1], "equals": d_exp[1],
         "tol": 1e-12},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 1, "coefficients", 0], "equals": d_c[0],
         "tol": 1e-12},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 1, "coefficients", 1], "equals": d_c[1],
         "tol": 1e-12},
        {"section": "wf", "member": "mo_coefficients", "binary_index": [2, 2],
         "equals": 1.0, "tol": 1e-12},
        # Appendix A.1: the self-overlap of the axial components is 1 and of a
        # mixed component 1/3 (for d, <xy|xy> = (1 * 1) / 3!!). A consumer
        # that evaluates the stored basis and gets 1 for d_xy has applied a
        # per-component correction the format forbids.
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["self_overlap", 0], "equals": 1.0, "tol": 1e-12},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["self_overlap", 1], "equals": 1.0, "tol": 1e-12},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["self_overlap", 2], "equals": 1.0 / 3.0,
         "tol": 1e-12},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["self_overlap", 3], "equals": 1.0 / 3.0,
         "tol": 1e-12},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["self_overlap", 4], "equals": 1.0, "tol": 1e-12},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["self_overlap", 5], "equals": 1.0 / 3.0,
         "tol": 1e-12},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["self_overlap", 6], "equals": 1.0, "tol": 1e-12},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["ao_labels", 2], "equals": "d_xy"},
        {"section": "overlap_check", "member": "payload",
         "json_pointer": ["cartesian_powers", 2], "equals": [1, 1, 0]},
    ], "wavefunction.gto with a Cartesian (pure: false) d shell: the "
       "Appendix A.1 normalization contract in numbers. The vendor section "
       "x_conformance.overlap_diagonal stores the self-overlap of each AO, "
       "computed from the stored exponents and coefficients with the format's "
       "N_i: 1 for the axial components and 1/3 for d_xy, d_xz and d_yz. A "
       "consumer recomputes them from basis.json to verify its normalization.")

    # 1-D wire: the periodic axis is y, not the leading one.
    w = QvfWriter(_PROGRAM, _VERSION, "carbon-wire")
    w.add_structure(
        [{"symbol": "C", "position": [0.0, 0.0, 0.0], "atomic_number": 6},
         {"symbol": "C", "position": [0.0, 1.28, 0.0], "atomic_number": 6}],
        pbc=[False, True, False],
        lattice_vectors=[[14.0, 0.0, 0.0], [0.0, 2.46, 0.0],
                         [0.0, 0.0, 17.5]],
        dimensionality=1)
    _write("structure_wire_1d", w, [
        {"section": "structure", "member": "structure",
         "json_pointer": ["pbc"], "equals": [False, True, False]},
        {"section": "structure", "member": "structure",
         "json_pointer": ["dimensionality"], "equals": 1},
        {"section": "structure", "member": "structure",
         "json_pointer": ["lattice_vectors", 1, 1], "equals": 2.46,
         "tol": 1e-12},
        # Rows 0 and 2 belong to aperiodic axes: bookkeeping, not cell edges.
        {"section": "structure", "member": "structure",
         "json_pointer": ["lattice_vectors", 0, 0], "equals": 14.0,
         "tol": 1e-12},
        {"section": "structure", "member": "structure",
         "json_pointer": ["lattice_vectors", 2, 2], "equals": 17.5,
         "tol": 1e-12},
    ], "1-D periodic wire (a carbon chain): pbc [false, true, false], so the "
       "periodic axis is y and is not the leading one; dimensionality 1. "
       "Rows 0 and 2 of lattice_vectors are synthesized normals of different "
       "lengths (14.0 and 17.5 Angstrom), not cell edges.")

    # Variable-cell reaction path: lattice [n_frames, 3, 3], columns, bohr.
    n_frames = 4
    rows = [[[4.20 - 0.03 * f, 0.0, 0.0],
             [0.40 + 0.10 * f, 4.10 - 0.02 * f, 0.0],
             [0.15, 0.25 + 0.05 * f, 4.30 + 0.04 * f]]
            for f in range(n_frames)]
    coords = np.array([[[0.0, 0.0, 0.0],
                        [2.10 - 0.03 * f, 2.05 - 0.02 * f, 2.15 - 0.03 * f]]
                       for f in range(n_frames)])
    w = QvfWriter(_PROGRAM, _VERSION, "mgo-variable-cell")
    w.add_structure(
        [{"symbol": "Mg", "position": list(coords[0, 0]), "atomic_number": 12},
         {"symbol": "O", "position": list(coords[0, 1]), "atomic_number": 8}],
        pbc=[True, True, True],
        lattice_vectors=Lattice.from_rows_angstrom(rows[0]))
    w.add_reaction_path(
        {"symbols": ["Mg", "O"], "atomic_numbers": [12, 8],
         "energies": [-274.50, -274.53, -274.55, -274.56],
         "kind": "optimization", "pbc": [True, True, True], "dim": 3,
         "n_frames": n_frames},
        coords, lattice=[Lattice.from_rows_angstrom(r) for r in rows],
        section_id="path")

    def lat(f, i, j):
        # L[f][i][j] = f_bohr * R[f][j][i]: column j of L is row j of R.
        return ang * rows[f][j][i]
    _write("reaction_path_variable_cell", w, [
        {"manifest_pointer": ["sections", 1, "members", "lattice", "shape"],
         "equals": [4, 3, 3]},
        {"manifest_pointer": ["sections", 1, "members", "lattice", "dtype"],
         "equals": "float64"},
        {"section": "path", "member": "metadata",
         "json_pointer": ["dim"], "equals": 3},
        {"section": "path", "member": "lattice", "binary_index": [0, 0, 0],
         "equals": lat(0, 0, 0), "tol": 1e-12},
        # b_x of the last frame: row 0 of column 1, so R[1][0], not R[0][1].
        {"section": "path", "member": "lattice", "binary_index": [3, 0, 1],
         "equals": lat(3, 0, 1), "tol": 1e-12},
        {"section": "path", "member": "lattice", "binary_index": [3, 1, 1],
         "equals": lat(3, 1, 1), "tol": 1e-12},
        {"section": "path", "member": "lattice", "binary_index": [2, 1, 2],
         "equals": lat(2, 1, 2), "tol": 1e-12},
        {"section": "path", "member": "lattice", "binary_index": [3, 2, 2],
         "equals": lat(3, 2, 2), "tol": 1e-12},
        {"section": "path", "member": "lattice", "binary_index": [1, 2, 1],
         "equals": 0.0, "tol": 1e-12},
        {"section": "path", "member": "coords", "binary_index": [3, 1, 2],
         "equals": float(coords[3, 1, 2]), "tol": 1e-12},
    ], "Periodic reaction.path with a variable cell: lattice float64 "
       "[n_frames, 3, 3], columns a, b, c in bohr, a different non-diagonal "
       "cell in every frame. The cells are not symmetric matrices, so a "
       "consumer that reads the columns as rows gets different numbers.")

    # The same logical archive in three containers. The archive is small so
    # that only the container is under test.
    def container_base():
        w = QvfWriter(_PROGRAM, _VERSION, "container-water")
        w.add_structure(
            [{"symbol": "O", "position": [0.0, 0.0, 0.1173],
              "atomic_number": 8},
             {"symbol": "H", "position": [0.0, 0.7572, -0.4692],
              "atomic_number": 1},
             {"symbol": "H", "position": [0.0, -0.7572, -0.4692],
              "atomic_number": 1}],
            bonds=[{"i": 0, "j": 1, "order": 1.0},
                   {"i": 0, "j": 2, "order": 1.0}])
        w.add_atom_properties(mulliken_charge=[-0.70, 0.35, 0.35])
        w.add_spectrum("spectra.ir", frequencies=[1595.0, 3657.0, 3756.0],
                       intensities=[67.0, 3.0, 36.0])
        return w

    container_checks = [
        {"manifest_pointer": ["qvf_version"], "equals": 1},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 1, "position", 1], "equals": 0.7572,
         "tol": 1e-12},
        {"section": "structure", "member": "bonds",
         "json_pointer": ["pairs", 1, "j"], "equals": 2},
        {"section": "atom_properties", "member": "mulliken_charge",
         "binary_index": [0], "equals": -0.70, "tol": 1e-12},
        {"section": "atom_properties", "member": "mulliken_charge",
         "binary_index": [2], "equals": 0.35, "tol": 1e-12},
        {"section": "spectra_ir", "member": "spectrum",
         "json_pointer": ["frequencies", 2], "equals": 3756.0, "tol": 1e-9},
    ]

    zip64_bytes = _zip64_zip(_entries_of(container_base().to_bytes()))
    assert zip64_bytes.count(b"PK\x06\x06") == 1       # ZIP64 end record
    assert zip64_bytes.count(b"PK\x06\x07") == 1       # ZIP64 locator
    _write("zip64", _Rebuilt(container_base(), _zip64_zip), container_checks,
           "ZIP64 forced on small members: a ZIP64 extra field in every local "
           "header and every central-directory record (sizes and local-header "
           "offset saturated), a ZIP64 end-of-central-directory record and "
           "locator, and a saturated end record. The member bytes are those "
           "of an ordinary archive of the same calculation (water: structure "
           "with bonds, atom_properties, IR spectrum).")

    descriptor_bytes = _data_descriptor_zip(
        _entries_of(container_base().to_bytes()))
    with zipfile.ZipFile(io.BytesIO(descriptor_bytes)) as zf:
        for info in zf.infolist():
            flags = struct.unpack_from("<H", descriptor_bytes,
                                       info.header_offset + 6)[0]
            assert flags & 0x08, info.filename          # data descriptor
    _write("data_descriptor", _Rebuilt(container_base(), _data_descriptor_zip),
           container_checks,
           "Members written in streaming mode: every local header carries the "
           "data-descriptor flag (bit 3) with zero CRC and zero sizes, a data "
           "descriptor follows each member's data, and only the central "
           "directory knows the sizes, so a reader must use the central "
           "directory. The same calculation as zip64.qvf; all members are "
           "DEFLATEd, manifest.json included.")

    # Legacy: qvf_version 2, as vibe-qc v0.10 - v0.15 stamped a periodic
    # reaction.path. Producers must not emit it; validators read it as
    # version 1 and report the warning W_VERSION_WITHDRAWN (spec § 3.1,
    # GOVERNANCE.md 2026-07-10 and 2026-10-14). The archive is otherwise a
    # valid v1 archive.
    slab = Lattice.from_rows_angstrom([[2.46, 0.0, 0.0], [1.23, 2.1304, 0.0],
                                       [0.0, 0.0, 15.875]])
    w = QvfWriter(_PROGRAM, _VERSION, "graphene-path-v2")
    w.add_structure(
        [{"symbol": "C", "position": [0.0, 0.0, 0.0], "atomic_number": 6},
         {"symbol": "C", "position": [1.23, 0.7101, 0.0], "atomic_number": 6}],
        pbc=[True, True, False], lattice_vectors=slab, dimensionality=2)
    w.add_reaction_path(
        {"symbols": ["C", "C"], "atomic_numbers": [6, 6],
         "energies": [-76.30, -76.31, -76.30], "kind": "scan",
         "pbc": [True, True, False], "dim": 2, "n_frames": 3},
        np.array([[[0.0, 0.0, 0.0], [1.23, 0.7101, dz]]
                  for dz in (-0.05, 0.0, 0.05)]),
        lattice=slab, section_id="path")

    def restamp(entries):
        out = []
        for name, data, method in entries:
            if name == "manifest.json":
                manifest = json.loads(data)
                manifest["qvf_version"] = 2
                manifest["schema_uri"] = (
                    "https://vibe-qc.org/spec/qvf/2/manifest.schema.json")
                data = json.dumps(manifest, ensure_ascii=False,
                                  indent=2).encode("utf-8")
            out.append((name, data, method))
        return _plain_zip(out)

    _write_legacy(
        "reaction_path_v2", restamp(_entries_of(w.to_bytes())),
        "A periodic reaction.path (fixed [3, 3] lattice member) stamped "
        "qvf_version 2, as vibe-qc v0.10 - v0.15 wrote it before the bump was "
        "withdrawn (GOVERNANCE.md, 2026-07-10). Apart from the version and "
        "the v2 schema_uri it is a valid v1 archive. Consumers SHOULD read it "
        "as version 1, and validators check it as version 1 and report the "
        "warning W_VERSION_WITHDRAWN, so the archive is valid: the warning "
        "records that producers must not write version 2 (spec § 3.1).",
        [{"code": "W_VERSION_WITHDRAWN", "severity": "warning"}])


def build() -> None:
    os.makedirs(CORPUS, exist_ok=True)

    # 1) minimal molecular structure + explicit bonds
    w = QvfWriter(_PROGRAM, _VERSION, "water")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}],
        bonds=[{"i": 0, "j": 1, "order": 1.0}, {"i": 0, "j": 2, "order": 1.0}])
    _write("structure_water", w, [
        {"manifest_pointer": ["qvf_version"], "equals": 1},
        {"manifest_pointer": ["source", "program"], "equals": _PROGRAM},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 0, "symbol"], "equals": "O"},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 1, "position", 1], "equals": 0.7572,
         "tol": 1e-9},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 0, "atomic_number"], "equals": 8},
        {"section": "structure", "member": "bonds",
         "json_pointer": ["pairs", 0, "order"], "equals": 1.0, "tol": 1e-9},
    ], "Minimal molecular structure (H2O) with explicit bonds; Angstrom.")

    # 2) periodic structure with a cubic lattice
    w = QvfWriter(_PROGRAM, _VERSION, "mgo")
    w.add_structure(
        [{"symbol": "Mg", "position": [0.0, 0.0, 0.0], "atomic_number": 12},
         {"symbol": "O", "position": [2.1, 2.1, 2.1], "atomic_number": 8}],
        pbc=[True, True, True],
        lattice_vectors=[[4.2, 0, 0], [0, 4.2, 0], [0, 0, 4.2]])
    _write("structure_periodic", w, [
        {"section": "structure", "member": "structure",
         "json_pointer": ["pbc", 0], "equals": True},
        {"section": "structure", "member": "structure",
         "json_pointer": ["lattice_vectors", 0, 0], "equals": 4.2, "tol": 1e-9},
    ], "Periodic structure (MgO) with a cubic lattice; row vectors, Angstrom.")

    # 2b) 2-D slab. pbc marks the third axis aperiodic, so lattice row 2 is a
    # synthesized normal that a consumer must not draw, tile, or measure. Its
    # length is deliberately vacuum-gap-sized (15.875 A ~ 30 bohr): a reader
    # keying off lattice geometry instead of off pbc gets this archive wrong.
    w = QvfWriter(_PROGRAM, _VERSION, "graphene")
    w.add_structure(
        [{"symbol": "C", "position": [0.0, 0.0, 0.0], "atomic_number": 6},
         {"symbol": "C", "position": [1.23, 0.7101, 0.0], "atomic_number": 6}],
        pbc=[True, True, False],
        lattice_vectors=[[2.46, 0.0, 0.0], [1.23, 2.1304, 0.0],
                         [0.0, 0.0, 15.875]],
        dimensionality=2)
    _write("structure_slab_2d", w, [
        {"section": "structure", "member": "structure",
         "json_pointer": ["pbc", 0], "equals": True},
        {"section": "structure", "member": "structure",
         "json_pointer": ["pbc", 1], "equals": True},
        {"section": "structure", "member": "structure",
         "json_pointer": ["pbc", 2], "equals": False},
        {"section": "structure", "member": "structure",
         "json_pointer": ["dimensionality"], "equals": 2},
        {"section": "structure", "member": "structure",
         "json_pointer": ["lattice_vectors", 2, 2], "equals": 15.875,
         "tol": 1e-9},
    ], "2-D slab (graphene): pbc [true, true, false], dimensionality 2. Lattice "
       "row 2 is a non-physical synthesized normal, present only so the matrix "
       "stays full-rank 3x3; it is not a cell edge.")

    # 3) restricted wavefunction.gto with the primitive-norm divide
    shells = [{"center": 0, "l": 0, "exponents": [1.0, 0.5],
               "coefficients": [0.6, 0.4]},
              {"center": 0, "l": 1, "exponents": [1.0], "coefficients": [1.0]}]
    C = np.eye(4, dtype=np.float64)
    w = QvfWriter(_PROGRAM, _VERSION, "wf_rhf")
    w.add_structure([{"symbol": "Ne", "position": [0, 0, 0], "atomic_number": 10}])
    w.add_wavefunction_gto(shells, mo_coefficients=C, energies=[-1, -0.5, 0.5, 1],
                           occupations=[2, 2, 0, 0],
                           coeffs_are_libint_normalized=True)
    # The stored coefficient must be the input divided by the primitive norm.
    expect_c0 = 0.6 / primitive_norm(1.0, 0)
    _write("wavefunction_rhf", w, [
        {"section": "wf", "member": "basis", "json_pointer": ["n_ao"], "equals": 4},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 0, "coefficients", 0], "equals": expect_c0,
         "tol": 1e-9},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 0, "pure"], "equals": True},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["spin"], "equals": "restricted"},
        {"section": "wf", "member": "mo_coefficients", "binary_index": [0, 0],
         "equals": 1.0, "tol": 1e-12},
    ], "Restricted wavefunction.gto; coefficients on N_i-normalized primitives "
       "(the N_i divide is the key contract; see spec Appendix A.1).")

    # 4) IR spectrum (frequency/intensity) + provenance + dipole (root metadata)
    w = QvfWriter(_PROGRAM, _VERSION, "ir")
    w.add_structure([{"symbol": "C", "position": [0, 0, 0], "atomic_number": 6},
                     {"symbol": "O", "position": [0, 0, 1.13], "atomic_number": 8}])
    w.add_spectrum("spectra.ir", frequencies=[2143.0], intensities=[500.0])
    # The frozen vector keeps the legacy spelling "Eh" that it was written
    # with (accepted but deprecated, spec § 4.3); the writer's scf_energy_eh
    # convenience now writes the unit name "hartree".
    w.set_provenance(scf_energy={"value": -113.3, "units": "Eh"}, method="RKS",
                     functional="PBE", basis="def2-SVP", scf_converged=True)
    w.set_dipole_moment(total_debye=0.11, vector_debye=[0.0, 0.0, 0.11],
                        origin="center_of_mass")
    _write("spectra_ir", w, [
        {"section": "spectra_ir", "member": "spectrum",
         "json_pointer": ["frequencies", 0], "equals": 2143.0, "tol": 1e-9},
        {"section": "spectra_ir", "member": "spectrum",
         "json_pointer": ["intensities", 0], "equals": 500.0, "tol": 1e-9},
        {"manifest_pointer": ["provenance", "functional"], "equals": "PBE"},
        {"manifest_pointer": ["provenance", "scf_energy", "value"],
         "equals": -113.3, "tol": 1e-9},
        {"manifest_pointer": ["provenance", "scf_energy", "units"], "equals": "Eh"},
        {"manifest_pointer": ["dipole_moment", "total_debye"], "equals": 0.11,
         "tol": 1e-9},
    ], "IR spectrum (cm^-1 / km*mol^-1) with root provenance + dipole metadata.")

    # 5) canonical EPR (object-shaped spectrum)
    w = QvfWriter(_PROGRAM, _VERSION, "epr")
    w.add_structure([{"symbol": "N", "position": [0, 0, 0], "atomic_number": 7}])
    w.add_spectrum("spectra.epr", payload={
        "g_tensor": {"principal": [2.0023, 2.0021, 2.0089], "isotropic": 2.0044},
        "hyperfine": [{"atom_index": 0, "symbol": "N", "a_iso_mhz": 45.2}],
        "zero_field_splitting": {"d_mhz": 1200.0, "e_mhz": 30.0}})
    _write("spectra_epr", w, [
        {"section": "spectra_epr", "member": "spectrum",
         "json_pointer": ["g_tensor", "principal", 2], "equals": 2.0089,
         "tol": 1e-9},
        {"section": "spectra_epr", "member": "spectrum",
         "json_pointer": ["zero_field_splitting", "d_mhz"], "equals": 1200.0,
         "tol": 1e-9},
    ], "Canonical spectra.epr: g-tensor / hyperfine / ZFS; g dimensionless, "
       "A and D in MHz.")

    # 6) volume.density grid (bohr) — binary + grid JSON contract
    rng = np.random.default_rng(0)
    dens = rng.random((6, 6, 6)).astype(np.float32)
    w = QvfWriter(_PROGRAM, _VERSION, "density")
    w.add_structure([{"symbol": "H", "position": [0, 0, 0], "atomic_number": 1}])
    w.add_density(dens, origin=[-1.5, -1.5, -1.5],
                  voxel_vectors=[[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
                  label="rho", section_id="rho")
    _write("volume_density", w, [
        {"section": "rho", "member": "grid",
         "json_pointer": ["shape", 0], "equals": 6},
        {"section": "rho", "member": "grid",
         "json_pointer": ["voxel_vectors", 0, 0], "equals": 0.5, "tol": 1e-9},
        {"section": "rho", "member": "data", "binary_index": [1, 2, 3],
         "equals": float(dens[1, 2, 3]), "tol": 1e-6},
    ], "volume.density grid; origin + voxel_vectors in bohr, float32 [nx,ny,nz].")

    # 7) band structure — rank-3 eigenvalue tensor (eV)
    eig = np.linspace(-5.0, 5.0, 1 * 8 * 3).reshape(1, 8, 3).astype(np.float64)
    w = QvfWriter(_PROGRAM, _VERSION, "bands")
    w.add_structure([{"symbol": "Si", "position": [0, 0, 0], "atomic_number": 14}],
                    pbc=[True, True, True],
                    lattice_vectors=[[3.8, 0, 0], [0, 3.8, 0], [0, 0, 3.8]])
    w.add_bands({"n_spin": 1, "n_kpoints": 8, "n_bands": 3, "fermi": 0.0,
                 "segments": [{"label_start": "G", "label_end": "X",
                               "k_start": [0, 0, 0], "k_end": [0.5, 0, 0],
                               "n_points": 8}]}, eig, section_id="bands0")
    _write("bands", w, [
        {"section": "bands0", "member": "kpath",
         "json_pointer": ["n_bands"], "equals": 3},
        {"section": "bands0", "member": "eigenvalues", "binary_index": [0, 0, 0],
         "equals": float(eig[0, 0, 0]), "tol": 1e-9},
        {"section": "bands0", "member": "eigenvalues", "binary_index": [0, 7, 2],
         "equals": float(eig[0, 7, 2]), "tol": 1e-9},
    ], "Band structure; eigenvalues [n_spin, n_kpoints, n_bands] in eV.")

    # 8) analysis: atom_properties + bond_orders + scf_history
    w = QvfWriter(_PROGRAM, _VERSION, "analysis")
    w.add_structure([{"symbol": "O", "position": [0, 0, 0], "atomic_number": 8},
                     {"symbol": "H", "position": [0, 0, 0.96], "atomic_number": 1}])
    w.add_atom_properties(mulliken_charge=[-0.34, 0.34],
                          loewdin_charge=[-0.20, 0.20])
    w.add_bond_orders("mayer", [{"i": 0, "j": 1, "order": 0.98,
                                 "distance_ang": 0.96}])
    w.add_scf_history([{"iter": 1, "energy_eh": -75.0, "delta_e": 1.0,
                        "diis_error": 0.1},
                       {"iter": 2, "energy_eh": -75.9, "delta_e": -0.9,
                        "diis_error": 1e-4}])
    _write("analysis", w, [
        {"section": "atom_properties", "member": "mulliken_charge",
         "binary_index": [0], "equals": -0.34, "tol": 1e-9},
        {"section": "bond_orders", "member": "bond_orders",
         "json_pointer": ["pairs", 0, "order"], "equals": 0.98, "tol": 1e-9},
        {"section": "scf_history", "member": "iterations",
         "json_pointer": ["iterations", 1, "energy_eh"], "equals": -75.9,
         "tol": 1e-9},
    ], "Analysis sections: atom_properties (float64 [n_atoms]), bond_orders, "
       "scf_history.")

    # 9) vendor extension + extensions governance
    w = QvfWriter(_PROGRAM, _VERSION, "vendor")
    w.add_structure([{"symbol": "H", "position": [0, 0, 0], "atomic_number": 1}])
    w.set_extensions({"x_demo": {"version": "1.0", "critical": False}})
    w.add_vendor_section("x_demo.custom",
                         json_members={"payload": {"answer": 42}},
                         critical=False)
    _write("vendor", w, [
        {"manifest_pointer": ["extensions", "x_demo", "version"], "equals": "1.0"},
        {"section": "x_demo_custom", "member": "payload",
         "json_pointer": ["answer"], "equals": 42},
    ], "Vendor namespace (x_demo.*) with an extensions governance declaration.")

    # 10) run.record — self-contained input + log of one program invocation
    input_text = ("! RHF STO-3G TightSCF\n"
                  "* xyz 0 1\n"
                  "O  0.0000  0.0000  0.1173\n"
                  "H  0.0000  0.7572 -0.4692\n"
                  "H  0.0000 -0.7572 -0.4692\n"
                  "*\n")
    log_text = ("democode 1.0.0 — water.inp\n"
                "SCF iteration   1  E = -74.880\n"
                "SCF iteration   8  E = -74.965901\n"
                "SCF CONVERGED\n"
                "FINAL SINGLE POINT ENERGY  -74.965901\n")
    w = QvfWriter(_PROGRAM, _VERSION, "run-record")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692],
          "atomic_number": 1}])
    w.add_run_record(
        program="democode", program_version="1.0.0",
        command="democode water.inp", exit_status=0,
        started_utc="2026-07-24T12:00:00Z", finished_utc="2026-07-24T12:00:41Z",
        input_text=input_text, log_text=log_text,
        files={"input": {"filename": "water.inp"},
               "log": {"filename": "water.out"},
               "attachment.geom": {"filename": "water.xyz",
                                   "media_type": "chemical/x-xyz"}},
        attachments={"geom": b"3\nwater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\n"
                             b"H 0 -0.7572 -0.4692\n"})
    _write("run_record", w, [
        {"manifest_pointer": ["sections", 1, "program"], "equals": "democode"},
        {"manifest_pointer": ["sections", 1, "program_version"],
         "equals": "1.0.0"},
        {"manifest_pointer": ["sections", 1, "exit_status"], "equals": 0},
        {"section": "run_record", "member": "input", "decode": "utf8",
         "equals": input_text},
        {"section": "run_record", "member": "log", "decode": "utf8",
         "contains": "FINAL SINGLE POINT ENERGY  -74.965901"},
        {"section": "run_record", "member": "files",
         "json_pointer": ["input", "filename"], "equals": "water.inp"},
        {"section": "run_record", "member": "attachment.geom",
         "decode": "utf8", "contains": "water"},
    ], "run.record: verbatim program input + full log + attachment, making "
       "the archive a self-contained record of the calculation.")

    # 11) job.spec — a pending archive: structure + declarative job request,
    # provenance.run_status = "pending" (nothing has run yet).
    w = QvfWriter(_PROGRAM, _VERSION, "pending-job")
    w.add_structure(
        [{"symbol": "Mg", "position": [0.0, 0.0, 0.0], "atomic_number": 12},
         {"symbol": "O", "position": [2.1, 2.1, 2.1], "atomic_number": 8}],
        pbc=[True, True, True],
        lattice_vectors=[[4.2, 0, 0], [0, 4.2, 0], [0, 0, 4.2]])
    w.add_job_spec(job_type="periodic", method="rks", basis="pob-tzvp-rev2",
                   functional="pbe", charge=0, multiplicity=1,
                   kpoints=[4, 4, 4], tasks=["single_point"],
                   options={"cutoff_ha": 300.0})
    w.set_provenance(run_status="pending")
    _write("job_spec", w, [
        {"manifest_pointer": ["provenance", "run_status"], "equals": "pending"},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["job_type"], "equals": "periodic"},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["method"], "equals": "rks"},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["kpoints", 1], "equals": 4},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["tasks", 0], "equals": "single_point"},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["options", "cutoff_ha"], "equals": 300.0,
         "tol": 1e-9},
    ], "job.spec: a pending container — structure + declarative job request "
       "(provenance.run_status = 'pending'), the executable-input half of "
       "run.record.")

    # 12) every 1-D spectrum kind (spec § 5.4, $defs/SpectrumPayload): the x
    # axis stays in `frequencies` in the kind's § 4.3 unit unless `x_units`
    # says otherwise; peaks, broadening and the kind blocks exercise the
    # optional fields. H2O2 (C2, chiral) so ECD and VCD are not zero by
    # symmetry; all values are synthetic.
    w = QvfWriter(_PROGRAM, _VERSION, "h2o2-spectra")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.725, 0.0], "atomic_number": 8},
         {"symbol": "O", "position": [0.0, -0.725, 0.0], "atomic_number": 8},
         {"symbol": "H", "position": [0.5342, 0.8934, 0.7919],
          "atomic_number": 1},
         {"symbol": "H", "position": [0.5342, -0.8934, -0.7919],
          "atomic_number": 1}])
    vib = [371.0, 877.0, 1266.0, 1402.0, 3599.0, 3608.0]
    w.add_spectrum(
        "spectra.ir", frequencies=vib,
        intensities=[160.0, 0.1, 90.0, 1.0, 20.0, 45.0],
        x_units="cm-1", y_units="km_per_mol",
        peaks=[{"x": 371.0, "y": 160.0, "label": "OO torsion"},
               {"x": 1266.0, "y": 90.0, "width": 12.0,
                "label": "OOH antisymmetric bend"}],
        broadening={"type": "lorentzian", "fwhm": 10.0, "units": "cm-1"})
    w.add_spectrum(
        "spectra.raman", payload={
            "frequencies": vib,
            "intensities": [2.1, 18.5, 1.3, 6.2, 95.0, 30.0],
            "temperature_k": 298.15,
            "raman": {"laser_nm": 532.0,
                      "depolarization_ratios": [0.75, 0.12, 0.75, 0.31,
                                                0.08, 0.75]}})
    w.add_spectrum(
        "spectra.uvvis", payload={
            "frequencies": [5.92, 7.21, 8.05],
            "intensities": [0.0012, 0.021, 0.054],
            "x": [5.92, 7.21, 8.05],
            "x_units": "ev", "y_units": "dimensionless",
            "peaks": [{"x": 8.05, "y": 0.054,
                       "assignment": {"state_index": 3}}],
            "broadening": {"type": "gaussian", "fwhm": 0.3, "units": "ev"},
            "uvvis": {"x_axis": "wavelength_nm"}})
    w.add_spectrum("spectra.ecd", frequencies=[5.92, 7.21, 8.05],
                   intensities=[0.8, -2.4, 1.1])
    w.add_spectrum("spectra.vcd", frequencies=vib,
                   intensities=[0.02, -0.001, 0.015, -0.004, 0.003, -0.006])
    w.add_spectrum("spectra.generic", frequencies=[11.7, 12.9, 15.5],
                   intensities=[1.0, 0.82, 0.35], x_units="ev",
                   y_units="dimensionless", label="Photoelectron spectrum")
    _write("spectra_all", w, [
        {"section": "spectra_ir", "member": "spectrum",
         "json_pointer": ["frequencies", 5], "equals": 3608.0, "tol": 1e-9},
        {"section": "spectra_ir", "member": "spectrum",
         "json_pointer": ["x_units"], "equals": "cm-1"},
        {"section": "spectra_ir", "member": "spectrum",
         "json_pointer": ["peaks", 1, "width"], "equals": 12.0, "tol": 1e-9},
        {"section": "spectra_ir", "member": "spectrum",
         "json_pointer": ["broadening", "type"], "equals": "lorentzian"},
        {"section": "spectra_raman", "member": "spectrum",
         "json_pointer": ["raman", "laser_nm"], "equals": 532.0, "tol": 1e-9},
        {"section": "spectra_raman", "member": "spectrum",
         "json_pointer": ["raman", "depolarization_ratios", 1],
         "equals": 0.12, "tol": 1e-9},
        {"section": "spectra_raman", "member": "spectrum",
         "json_pointer": ["temperature_k"], "equals": 298.15, "tol": 1e-9},
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["x", 2], "equals": 8.05, "tol": 1e-9},
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["peaks", 0, "assignment", "state_index"],
         "equals": 3},
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["broadening", "fwhm"], "equals": 0.3, "tol": 1e-9},
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["uvvis", "x_axis"], "equals": "wavelength_nm"},
        {"section": "spectra_ecd", "member": "spectrum",
         "json_pointer": ["intensities", 1], "equals": -2.4, "tol": 1e-9},
        {"section": "spectra_vcd", "member": "spectrum",
         "json_pointer": ["frequencies", 0], "equals": 371.0, "tol": 1e-9},
        {"section": "spectra_generic", "member": "spectrum",
         "json_pointer": ["y_units"], "equals": "dimensionless"},
        {"manifest_pointer": ["sections", 6, "label"],
         "equals": "Photoelectron spectrum"},
    ], "All six 1-D spectrum kinds ($defs/SpectrumPayload) of a synthetic "
       "H2O2: frequencies carry the x axis in the kind's spec section 4.3 unit (cm^-1 "
       "vibrational, eV electronic) unless x_units declares one; IR with peaks "
       "and Lorentzian broadening, Raman and UV-Vis kind blocks, the x alias, "
       "and a generic spectrum with declared units.")

    # 13) typed spectra.nmr (spec § 5.4, $defs/NmrSpectrumPayload): shifts in
    # ppm relative to `reference`, shielding tensors in ppm, J in Hz.
    w = QvfWriter(_PROGRAM, _VERSION, "methanol-nmr")
    w.add_structure(
        [{"symbol": "C", "position": [-0.047, 0.665, 0.0], "atomic_number": 6},
         {"symbol": "O", "position": [-0.047, -0.757, 0.0], "atomic_number": 8},
         {"symbol": "H", "position": [0.866, -1.068, 0.0], "atomic_number": 1},
         {"symbol": "H", "position": [-1.086, 0.977, 0.0], "atomic_number": 1},
         {"symbol": "H", "position": [0.437, 1.078, 0.89], "atomic_number": 1},
         {"symbol": "H", "position": [0.437, 1.078, -0.89],
          "atomic_number": 1}])
    w.add_spectrum("spectra.nmr", payload={
        "isotope": "1H", "reference": "TMS", "solvent": "chloroform",
        "chemical_shifts": [
            {"atom_index": 2, "symbol": "H", "isotropic_shift_ppm": 1.52,
             "anisotropy_ppm": 17.9, "asymmetry": 0.21},
            {"atom_index": 3, "symbol": "H", "isotropic_shift_ppm": 3.41,
             "anisotropy_ppm": 8.6, "asymmetry": 0.35},
            {"atom_index": 4, "symbol": "H", "isotropic_shift_ppm": 3.39},
            {"atom_index": 5, "symbol": "H", "isotropic_shift_ppm": 3.39}],
        "shielding_tensors": [
            {"atom_index": 2, "isotropic_ppm": 30.28,
             "tensor": [[24.31, 1.02, 0.0], [1.02, 42.21, 0.0],
                        [0.0, 0.0, 24.32]]}],
        "j_couplings": [
            {"i": 2, "j": 3, "value_hz": 5.1, "isotope_i": "1H",
             "isotope_j": "1H"},
            {"i": 4, "j": 5, "value_hz": -10.8, "isotope_i": "1H",
             "isotope_j": "1H"}]})
    _write("nmr", w, [
        {"section": "spectra_nmr", "member": "spectrum",
         "json_pointer": ["isotope"], "equals": "1H"},
        {"section": "spectra_nmr", "member": "spectrum",
         "json_pointer": ["solvent"], "equals": "chloroform"},
        {"section": "spectra_nmr", "member": "spectrum",
         "json_pointer": ["chemical_shifts", 1, "isotropic_shift_ppm"],
         "equals": 3.41, "tol": 1e-9},
        {"section": "spectra_nmr", "member": "spectrum",
         "json_pointer": ["chemical_shifts", 0, "asymmetry"], "equals": 0.21,
         "tol": 1e-9},
        {"section": "spectra_nmr", "member": "spectrum",
         "json_pointer": ["shielding_tensors", 0, "tensor", 1, 1],
         "equals": 42.21, "tol": 1e-9},
        {"section": "spectra_nmr", "member": "spectrum",
         "json_pointer": ["j_couplings", 1, "value_hz"], "equals": -10.8,
         "tol": 1e-9},
    ], "Typed spectra.nmr ($defs/NmrSpectrumPayload) of methanol: 1H "
       "chemical shifts (ppm vs TMS), one shielding tensor (ppm, Cartesian "
       "rows) and J couplings (Hz); synthetic values.")

    # 14) topology.qtaim with a declared position unit (spec § 5.7,
    # $defs/QtaimPayload): nuclear and bond critical points of water plus the
    # bond paths through the BCPs; positions in Angstrom, the frame of the
    # structure section. rho and laplacian carry no declared unit.
    water = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
             {"symbol": "H", "position": [0.0, 0.7572, -0.4692],
              "atomic_number": 1},
             {"symbol": "H", "position": [0.0, -0.7572, -0.4692],
              "atomic_number": 1}]
    w = QvfWriter(_PROGRAM, _VERSION, "water-qtaim")
    w.add_structure(water)
    w.add_topology_qtaim({
        "points": [
            {"type": "ncp", "position": [0.0, 0.0, 0.1173], "atoms": [0]},
            {"type": "ncp", "position": [0.0, 0.7572, -0.4692], "atoms": [1]},
            {"type": "ncp", "position": [0.0, -0.7572, -0.4692], "atoms": [2]},
            {"type": "bcp", "position": [0.0, 0.5788, -0.3149],
             "rho": 0.3653, "laplacian": -2.4521, "ellipticity": 0.0226,
             "atom_pair": [0, 1]},
            {"type": "bcp", "position": [0.0, -0.5788, -0.3149],
             "rho": 0.3653, "laplacian": -2.4521, "ellipticity": 0.0226,
             "atom_pair": [0, 2]}],
        "bond_paths": [
            {"atoms": [0, 1], "path": [[0.0, 0.0, 0.1173],
                                       [0.0, 0.5788, -0.3149],
                                       [0.0, 0.7572, -0.4692]]},
            {"atoms": [0, 2], "path": [[0.0, 0.0, 0.1173],
                                       [0.0, -0.5788, -0.3149],
                                       [0.0, -0.7572, -0.4692]]}]},
        units={"position": "angstrom"})
    _write("qtaim", w, [
        {"section": "qtaim", "member": "critical_points",
         "json_pointer": ["units", "position"], "equals": "angstrom"},
        {"section": "qtaim", "member": "critical_points",
         "json_pointer": ["points", 0, "type"], "equals": "ncp"},
        {"section": "qtaim", "member": "critical_points",
         "json_pointer": ["points", 3, "type"], "equals": "bcp"},
        {"section": "qtaim", "member": "critical_points",
         "json_pointer": ["points", 3, "position", 1], "equals": 0.5788,
         "tol": 1e-9},
        {"section": "qtaim", "member": "critical_points",
         "json_pointer": ["points", 4, "atom_pair", 1], "equals": 2},
        {"section": "qtaim", "member": "critical_points",
         "json_pointer": ["points", 3, "laplacian"], "equals": -2.4521,
         "tol": 1e-9},
        {"section": "qtaim", "member": "critical_points",
         "json_pointer": ["bond_paths", 1, "path", 2, 1], "equals": -0.7572,
         "tol": 1e-9},
    ], "topology.qtaim ($defs/QtaimPayload) of water: three nuclear and two "
       "bond critical points and the bond paths through the BCPs, positions "
       "declared in Angstrom (units.position), overlaying the structure; "
       "synthetic values.")

    # 15) structure.symmetry (spec § 5.7, $defs/SymmetryPayload) of CsCl,
    # Pm-3m: the 48 operations of m-3m act on fractional coordinates of the
    # cubic lattice, so each rotation is a signed permutation matrix and every
    # translation is zero. Identity first.
    ops = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            rot = [[0, 0, 0] for _ in range(3)]
            for row, col in enumerate(perm):
                rot[row][col] = signs[row]
            ops.append({"rotation": rot, "translation": [0.0, 0.0, 0.0]})
    w = QvfWriter(_PROGRAM, _VERSION, "cscl-symmetry")
    w.add_structure(
        [{"symbol": "Cs", "position": [0.0, 0.0, 0.0], "atomic_number": 55},
         {"symbol": "Cl", "position": [2.06, 2.06, 2.06], "atomic_number": 17}],
        pbc=[True, True, True],
        lattice_vectors=[[4.12, 0, 0], [0, 4.12, 0], [0, 0, 4.12]])
    w.add_structure_symmetry({
        "spacegroup_number": 221, "international": "Pm-3m",
        "hall": "-P 4 2 3", "point_group": "Oh", "hermann_mauguin": "m-3m",
        "operations": ops, "wyckoffs": ["a", "b"], "equivalent_atoms": [0, 1],
        "symprec": 1e-05})
    _write("symmetry", w, [
        {"section": "symmetry", "member": "data",
         "json_pointer": ["spacegroup_number"], "equals": 221},
        {"section": "symmetry", "member": "data",
         "json_pointer": ["international"], "equals": "Pm-3m"},
        {"section": "symmetry", "member": "data",
         "json_pointer": ["hall"], "equals": "-P 4 2 3"},
        {"section": "symmetry", "member": "data",
         "json_pointer": ["hermann_mauguin"], "equals": "m-3m"},
        {"section": "symmetry", "member": "data",
         "json_pointer": ["operations", 0, "rotation", 1, 1], "equals": 1},
        {"section": "symmetry", "member": "data",
         "json_pointer": ["operations", 47, "translation", 2], "equals": 0.0,
         "tol": 1e-12},
        {"section": "symmetry", "member": "data",
         "json_pointer": ["wyckoffs", 1], "equals": "b"},
        {"section": "symmetry", "member": "data",
         "json_pointer": ["symprec"], "equals": 1e-05, "tol": 1e-12},
    ], "structure.symmetry ($defs/SymmetryPayload) of CsCl in Pm-3m (No. "
       "221): Hall symbol, point group in Schoenflies and Hermann-Mauguin "
       "notation, all 48 operations on fractional coordinates, Wyckoff "
       "letters and equivalent atoms.")
    # 12) trajectory: a geometry optimization with typed metadata, plus a
    # reaction.waypoints annotation layered over it (spec § 5.5).
    h2o_symbols = ["O", "H", "H"]
    base = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692],
                     [0.0, -0.7572, -0.4692]])
    stretch = np.array([[0.0, 0.0, 0.0], [0.0, 0.01, -0.005],
                        [0.0, -0.01, -0.005]])
    opt_coords = np.stack([base + (3 - k) * stretch for k in range(4)])
    opt_energies = [-74.9580, -74.9631, -74.9652, -74.9659]
    w = QvfWriter(_PROGRAM, _VERSION, "h2o-opt")
    w.add_structure([{"symbol": z, "position": list(map(float, xyz)),
                      "atomic_number": n}
                     for z, xyz, n in zip(h2o_symbols, base, [8, 1, 1])])
    w.add_trajectory({"symbols": h2o_symbols, "atomic_numbers": [8, 1, 1],
                      "n_frames": 4, "energies": opt_energies,
                      "kind": "optimization", "pbc": [False, False, False],
                      "precision": "float64"},
                     opt_coords, section_id="opt")
    w.add_reaction_waypoints("opt", {"waypoints": [
        {"frame_index": 0, "label": "start", "kind": "point",
         "energy_eh": opt_energies[0]},
        {"frame_index": 3, "label": "minimum", "kind": "product",
         "energy_eh": opt_energies[3]}]}, section_id="opt_marks")
    _write("trajectory", w, [
        {"section": "opt", "member": "metadata",
         "json_pointer": ["kind"], "equals": "optimization"},
        {"section": "opt", "member": "metadata",
         "json_pointer": ["symbols", 1], "equals": "H"},
        {"section": "opt", "member": "metadata",
         "json_pointer": ["atomic_numbers", 0], "equals": 8},
        {"section": "opt", "member": "metadata",
         "json_pointer": ["energies", 3], "equals": -74.9659, "tol": 1e-9},
        {"section": "opt", "member": "metadata",
         "json_pointer": ["precision"], "equals": "float64"},
        {"section": "opt", "member": "coords", "binary_index": [0, 1, 1],
         "equals": float(opt_coords[0, 1, 1]), "tol": 1e-12},
        {"manifest_pointer": ["sections", 2, "trajectory_ref"], "equals": "opt"},
        {"section": "opt_marks", "member": "waypoints",
         "json_pointer": ["waypoints", 1, "frame_index"], "equals": 3},
        {"section": "opt_marks", "member": "waypoints",
         "json_pointer": ["waypoints", 1, "kind"], "equals": "product"},
    ], "trajectory: a 4-frame H2O optimization; coords [n_frames, n_atoms, 3] "
       "in Angstrom, metadata energies in Hartree, plus reaction.waypoints "
       "annotations referencing it by trajectory_ref.")

    # 13) vibrations: planar NH3 (the inversion transition state), whose one
    # imaginary mode is written as a negative wavenumber (spec § 5.5).
    nh3 = [{"symbol": "N", "position": [0.0, 0.0, 0.0], "atomic_number": 7},
           {"symbol": "H", "position": [1.0, 0.0, 0.0], "atomic_number": 1},
           {"symbol": "H", "position": [-0.5, 0.866, 0.0], "atomic_number": 1},
           {"symbol": "H", "position": [-0.5, -0.866, 0.0], "atomic_number": 1}]
    freqs = [-893.0, 1630.0, 1630.0, 3524.0, 3720.0, 3720.0]
    disp = np.zeros((6, 4, 3))
    disp[0, :, 2] = [0.27, -0.55, -0.55, -0.55]    # umbrella (out of plane)
    disp[1:] = np.random.default_rng(13).standard_normal((5, 4, 3)) * 0.1
    w = QvfWriter(_PROGRAM, _VERSION, "nh3-ts-freq")
    w.add_structure(nh3)
    w.add_vibrations({"frequencies": freqs,
                      "ir_intensities": [512.0, 18.5, 18.5, 0.0, 6.2, 6.2],
                      "reduced_masses": [1.17, 1.08, 1.08, 1.01, 1.10, 1.10],
                      "symmetry_labels": ["A2''", "E'", "E'", "A1'", "E'", "E'"],
                      "imaginary": [True, False, False, False, False, False],
                      "n_modes": 6, "symbols": ["N", "H", "H", "H"]},
                     disp, section_id="vib")
    _write("vibrations", w, [
        {"section": "vib", "member": "metadata",
         "json_pointer": ["frequencies", 0], "equals": -893.0, "tol": 1e-9},
        {"section": "vib", "member": "metadata",
         "json_pointer": ["imaginary", 0], "equals": True},
        {"section": "vib", "member": "metadata",
         "json_pointer": ["frequencies", 5], "equals": 3720.0, "tol": 1e-9},
        {"section": "vib", "member": "metadata",
         "json_pointer": ["ir_intensities", 0], "equals": 512.0, "tol": 1e-9},
        {"section": "vib", "member": "metadata",
         "json_pointer": ["symmetry_labels", 0], "equals": "A2''"},
        {"section": "vib", "member": "metadata",
         "json_pointer": ["n_modes"], "equals": 6},
        {"section": "vib", "member": "displacements", "binary_index": [0, 1, 2],
         "equals": -0.55, "tol": 1e-12},
    ], "vibrations: planar NH3 with one imaginary mode, written as a negative "
       "wavenumber (cm^-1); displacements [n_modes, n_atoms, 3].")

    # 14) scan.surface: a relaxed 2-D scan of r(O-H1) x angle(H1-O-H2), with
    # structured axis descriptions and the relaxed geometry at each node.
    r_ax = [1.70, 1.80, 1.90]                  # bohr
    a_ax = [1.75, 1.85]                        # rad
    scan_e = np.array([[-74.9601, -74.9605], [-74.9652, -74.9659],
                       [-74.9633, -74.9641]])
    geoms = np.stack([base] * 6)
    w = QvfWriter(_PROGRAM, _VERSION, "h2o-scan")
    w.add_structure([{"symbol": z, "position": list(map(float, xyz)),
                      "atomic_number": n}
                     for z, xyz, n in zip(h2o_symbols, base, [8, 1, 1])])
    w.add_scan_surface(
        {"axis_a": {"label": "r(O-H1)", "unit": "bohr", "atoms": [0, 1]},
         "axis_b": {"label": "a(H1-O-H2)", "unit": "rad", "atoms": [1, 0, 2]},
         "atoms": [{"symbol": z, "atomic_number": n}
                   for z, n in zip(h2o_symbols, [8, 1, 1])]},
        r_ax, a_ax, scan_e, geometries=geoms, section_id="scan")
    _write("scan_surface", w, [
        {"section": "scan", "member": "metadata",
         "json_pointer": ["axis_a", "label"], "equals": "r(O-H1)"},
        {"section": "scan", "member": "metadata",
         "json_pointer": ["axis_a", "unit"], "equals": "bohr"},
        {"section": "scan", "member": "metadata",
         "json_pointer": ["axis_b", "atoms", 1], "equals": 0},
        {"section": "scan", "member": "metadata",
         "json_pointer": ["atoms", 0, "symbol"], "equals": "O"},
        {"section": "scan", "member": "axis_a", "binary_index": [2],
         "equals": 1.90, "tol": 1e-12},
        {"section": "scan", "member": "energies", "binary_index": [1, 1],
         "equals": -74.9659, "tol": 1e-12},
        {"section": "scan", "member": "geometries", "binary_index": [5, 2, 1],
         "equals": -0.7572, "tol": 1e-12},
    ], "scan.surface: a relaxed 3x2 scan; axis values in bohr / rad, energies "
       "[nA, nB] in Hartree, relaxed geometries [nA*nB, n_atoms, 3] in "
       "Angstrom.")

    # 15) reaction.path: the molecular HCN -> HNC isomerization with
    # reactant / transition-state / product waypoints and a reaction
    # coordinate. No `lattice` member, so the path is molecular (spec § 5.5).
    hcn = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 1.156], [0.0, 0.0, -1.066]],
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.160], [0.55, 0.0, -0.85]],
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.170], [1.10, 0.0, 0.40]],
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.172], [0.60, 0.0, 1.85]],
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.172], [0.0, 0.0, 2.166]]])
    path_e = [-93.1900, -93.1500, -93.1150, -93.1520, -93.1660]
    w = QvfWriter(_PROGRAM, _VERSION, "hcn-hnc")
    w.add_structure([{"symbol": "C", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 6},
                     {"symbol": "N", "position": [0.0, 0.0, 1.156],
                      "atomic_number": 7},
                     {"symbol": "H", "position": [0.0, 0.0, -1.066],
                      "atomic_number": 1}])
    w.add_reaction_path(
        {"symbols": ["C", "N", "H"], "atomic_numbers": [6, 7, 1],
         "energies": path_e, "kind": "irc",
         "reaction_coordinate": [-2.0, -1.0, 0.0, 1.0, 2.0],
         "reaction_coordinate_label": "IRC arc length",
         "reaction_coordinate_unit": "bohr",
         "waypoints": [
             {"frame_index": 0, "label": "HCN", "kind": "reactant",
              "energy_eh": path_e[0]},
             {"frame_index": 2, "label": "TS", "kind": "transition_state",
              "energy_eh": path_e[2]},
             {"frame_index": 4, "label": "HNC", "kind": "product",
              "energy_eh": path_e[4]}]},
        hcn, section_id="rxn")
    _write("reaction_path", w, [
        {"section": "rxn", "member": "metadata",
         "json_pointer": ["kind"], "equals": "irc"},
        {"section": "rxn", "member": "metadata",
         "json_pointer": ["waypoints", 1, "kind"], "equals": "transition_state"},
        {"section": "rxn", "member": "metadata",
         "json_pointer": ["waypoints", 1, "frame_index"], "equals": 2},
        {"section": "rxn", "member": "metadata",
         "json_pointer": ["waypoints", 1, "energy_eh"], "equals": -93.1150,
         "tol": 1e-9},
        {"section": "rxn", "member": "metadata",
         "json_pointer": ["reaction_coordinate", 4], "equals": 2.0,
         "tol": 1e-12},
        {"section": "rxn", "member": "metadata",
         "json_pointer": ["energies", 0], "equals": -93.19, "tol": 1e-9},
        {"section": "rxn", "member": "coords", "binary_index": [2, 2, 0],
         "equals": 1.10, "tol": 1e-12},
    ], "reaction.path (molecular, no lattice member): HCN -> HNC with "
       "reactant / transition_state / product waypoints (energy_eh in "
       "Hartree) and a reaction coordinate.")

    # 12) phonons: a band structure with complex eigenvectors, a Γ-only set of
    # real ones, and a phonon DOS with its per-atom projection (spec § 5.6).
    n_q, n_atoms = 4, 2
    n_modes = 3 * n_atoms
    ph_freqs = np.linspace(0.0, 700.0, n_q * n_modes).reshape(n_q, n_modes)
    ev_size = n_q * n_modes * n_atoms * 3
    ev_re = np.linspace(-0.5, 0.5, ev_size).reshape(n_q, n_modes, n_atoms, 3)
    ev_im = np.linspace(0.25, -0.25, ev_size).reshape(n_q, n_modes, n_atoms, 3)
    # Mode i of the Γ set moves atom i // 3 along axis i % 3.
    ev_gamma = np.eye(n_modes).reshape(1, n_modes, n_atoms, 3)
    qpath = {"n_atoms": n_atoms, "n_modes": n_modes, "has_eigenvectors": True,
             "segments": [{"label_start": "G", "label_end": "X",
                           "k_start": [0, 0, 0], "k_end": [0.5, 0, 0.5],
                           "n_points": n_q}]}
    dos_grid = np.linspace(0.0, 700.0, 8)
    dos_mg = np.array([0.0, 0.1, 0.4, 0.9, 0.5, 0.2, 0.1, 0.0])
    dos_o = np.array([0.0, 0.0, 0.1, 0.2, 0.6, 1.1, 0.7, 0.0])
    w = QvfWriter(_PROGRAM, _VERSION, "mgo-phonons")
    w.add_structure(
        [{"symbol": "Mg", "position": [0.0, 0.0, 0.0], "atomic_number": 12},
         {"symbol": "O", "position": [2.1, 2.1, 2.1], "atomic_number": 8}],
        pbc=[True, True, True],
        lattice_vectors=[[0.0, 2.1, 2.1], [2.1, 0.0, 2.1], [2.1, 2.1, 0.0]])
    w.add_phonon_bands(qpath, ph_freqs, eigenvectors=ev_re + 1j * ev_im)
    gamma_qpath = dict(qpath, segments=[{"label_start": "G", "label_end": "G",
                                         "k_start": [0, 0, 0], "k_end": [0, 0, 0],
                                         "n_points": 1}])
    w.add_phonon_bands(gamma_qpath, ph_freqs[:1], eigenvectors=ev_gamma,
                       section_id="phonon_gamma")
    w.add_phonon_dos({"n_atoms": n_atoms, "n_modes": n_modes,
                      "channels": [{"atom": 0, "label": "Mg"},
                                   {"atom": 1, "label": "O"}]},
                     dos_grid, dos_mg + dos_o,
                     projected=np.stack([dos_mg, dos_o]))
    _write("phonons", w, [
        {"manifest_pointer": ["sections", 1, "members", "eigenvectors", "dtype"],
         "equals": "float64"},
        {"manifest_pointer": ["sections", 1, "members", "eigenvectors", "shape"],
         "equals": [n_q, n_modes, n_atoms, 3, 2]},
        {"section": "phonon_bands", "member": "frequencies",
         "binary_index": [3, 5], "equals": float(ph_freqs[3, 5]), "tol": 1e-9},
        {"section": "phonon_bands", "member": "eigenvectors",
         "binary_index": [1, 2, 0, 1, 0], "equals": float(ev_re[1, 2, 0, 1]),
         "tol": 1e-12},
        {"section": "phonon_bands", "member": "eigenvectors",
         "binary_index": [1, 2, 0, 1, 1], "equals": float(ev_im[1, 2, 0, 1]),
         "tol": 1e-12},
        {"manifest_pointer": ["sections", 2, "members", "eigenvectors", "shape"],
         "equals": [1, n_modes, n_atoms, 3]},
        {"section": "phonon_gamma", "member": "eigenvectors",
         "binary_index": [0, 4, 1, 1], "equals": 1.0, "tol": 1e-12},
        {"section": "phonon_gamma", "member": "eigenvectors",
         "binary_index": [0, 4, 0, 1], "equals": 0.0, "tol": 1e-12},
        {"manifest_pointer": ["sections", 3, "members", "projected", "shape"],
         "equals": [2, 8]},
        {"section": "phonon_dos", "member": "projected", "binary_index": [1, 5],
         "equals": 1.1, "tol": 1e-12},
        {"section": "phonon_dos", "member": "dos", "binary_index": [5],
         "equals": 1.3, "tol": 1e-12},
        {"section": "phonon_dos", "member": "meta",
         "json_pointer": ["channels", 1, "label"], "equals": "O"},
        {"section": "phonon_dos", "member": "meta",
         "json_pointer": ["channels", 1, "atom"], "equals": 1},
    ], "Phonons (MgO): phonon_bands with complex eigenvectors "
       "[n_q, n_modes, n_atoms, 3, 2] (last axis [re, im]) and a Γ-only set of "
       "real ones [n_q, n_modes, n_atoms, 3]; phonon_dos with a projected DOS "
       "[n_channels, n_points] labelled by meta.channels. Frequencies in cm^-1.")

    # 13) atom_properties carrying every role of the schema's role list.
    w = QvfWriter(_PROGRAM, _VERSION, "oh-radical")
    w.add_structure([{"symbol": "O", "position": [0, 0, 0], "atomic_number": 8},
                     {"symbol": "H", "position": [0, 0, 0.97], "atomic_number": 1}])
    w.add_atom_properties(mulliken_charge=[-0.29, 0.29],
                          loewdin_charge=[-0.18, 0.18],
                          hirshfeld_charge=[-0.12, 0.12],
                          iao_charge=[-0.35, 0.35],
                          spin_population=[1.05, -0.05])
    _write("atom_properties_all", w, [
        {"manifest_pointer": ["sections", 1, "members", "iao_charge", "dtype"],
         "equals": "float64"},
        {"manifest_pointer": ["sections", 1, "members", "hirshfeld_charge",
                              "shape"], "equals": [2]},
        {"section": "atom_properties", "member": "mulliken_charge",
         "binary_index": [0], "equals": -0.29, "tol": 1e-12},
        {"section": "atom_properties", "member": "loewdin_charge",
         "binary_index": [0], "equals": -0.18, "tol": 1e-12},
        {"section": "atom_properties", "member": "hirshfeld_charge",
         "binary_index": [1], "equals": 0.12, "tol": 1e-12},
        {"section": "atom_properties", "member": "iao_charge",
         "binary_index": [0], "equals": -0.35, "tol": 1e-12},
        {"section": "atom_properties", "member": "spin_population",
         "binary_index": [0], "equals": 1.05, "tol": 1e-12},
    ], "atom_properties with all five roles of the schema's role list "
       "(mulliken_charge, loewdin_charge, hirshfeld_charge, iao_charge, "
       "spin_population), each float64 [n_atoms]; illustrative values for the "
       "OH radical.")

    print(f"Corpus written to {CORPUS}")


    # 12) spin-polarized band structure (spec § 5.6): eigenvalues
    # [n_spin, n_k, n_bands], fat-band projections with the spin axis
    # [n_spin, n_k, n_bands, n_channels], occupations [n_spin, n_k, n_bands],
    # and a kpath carrying k-points, per-point labels, explicit segments and
    # channel objects.
    fe_bcc = dict(
        atoms=[{"symbol": "Fe", "position": [0.0, 0.0, 0.0], "atomic_number": 26},
               {"symbol": "Fe", "position": [1.435, 1.435, 1.435],
                "atomic_number": 26}],
        pbc=[True, True, True],
        lattice_vectors=[[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]])
    n_k, n_b = 5, 3
    eig_up = np.linspace(-6.0, 4.0, n_k * n_b).reshape(n_k, n_b)
    eig_s = np.stack([eig_up, eig_up + 1.5])                 # [2, 5, 3], eV
    fermi = 0.5
    occ_s = (eig_s < fermi).astype(np.float64)               # [2, 5, 3]
    weight = np.linspace(0.1, 0.9, 2 * n_k * n_b).reshape(2, n_k, n_b)
    proj_s = np.stack([weight, 1.0 - weight], axis=-1)       # [2, 5, 3, 2]
    w = QvfWriter(_PROGRAM, _VERSION, "bands-spin")
    w.add_structure(**fe_bcc)
    w.add_bands({"kpoints": [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.5, 0.0, 0.0],
                             [0.5, 0.25, 0.0], [0.5, 0.5, 0.0]],
                 "labels": ["G", "", "X", "", "M"],
                 "segments": [{"start": 0, "end": 2,
                               "label_start": "G", "label_end": "X"},
                              {"start": 2, "end": 4,
                               "label_start": "X", "label_end": "M"}],
                 "fermi": fermi, "spinor": False,
                 "n_spin": 2, "n_kpoints": n_k, "n_bands": n_b,
                 "channels": [{"l": 0, "label": "Fe-s"},
                              {"l": 2, "label": "Fe-d"}]},
                eig_s, projections=proj_s, occupations=occ_s,
                section_id="bands_spin")
    _write("bands_spin", w, [
        {"section": "bands_spin", "member": "kpath",
         "json_pointer": ["n_spin"], "equals": 2},
        {"section": "bands_spin", "member": "kpath",
         "json_pointer": ["kpoints", 2, 0], "equals": 0.5, "tol": 1e-12},
        {"section": "bands_spin", "member": "kpath",
         "json_pointer": ["labels", 4], "equals": "M"},
        {"section": "bands_spin", "member": "kpath",
         "json_pointer": ["segments", 1, "start"], "equals": 2},
        {"section": "bands_spin", "member": "kpath",
         "json_pointer": ["channels", 1, "l"], "equals": 2},
        {"section": "bands_spin", "member": "kpath",
         "json_pointer": ["spinor"], "equals": False},
        {"manifest_pointer": ["sections", 1, "members", "projections", "shape"],
         "equals": [2, n_k, n_b, 2]},
        {"manifest_pointer": ["sections", 1, "members", "occupations", "shape"],
         "equals": [2, n_k, n_b]},
        {"section": "bands_spin", "member": "eigenvalues",
         "binary_index": [1, 4, 2], "equals": float(eig_s[1, 4, 2]),
         "tol": 1e-12},
        {"section": "bands_spin", "member": "projections",
         "binary_index": [1, 4, 2, 1], "equals": float(proj_s[1, 4, 2, 1]),
         "tol": 1e-12},
        {"section": "bands_spin", "member": "projections",
         "binary_index": [0, 0, 0, 0], "equals": float(proj_s[0, 0, 0, 0]),
         "tol": 1e-12},
        {"section": "bands_spin", "member": "occupations",
         "binary_index": [0, 0, 0], "equals": 1.0},
        {"section": "bands_spin", "member": "occupations",
         "binary_index": [1, 4, 2], "equals": 0.0},
    ], "Spin-polarized band structure: eigenvalues [n_spin, n_k, n_bands], "
       "fat-band projections with a spin axis [n_spin, n_k, n_bands, "
       "n_channels], occupations [n_spin, n_k, n_bands]; kpath with k-points, "
       "per-point labels, explicit segments and channel objects.")

    # 13) DOS family, spin-polarized: dos.total [n_spin, n_points],
    # dos.projected [n_spin, n_channels, n_points], and a dos.cohp whose pair
    # axis is described by its meta payload.
    e_dos = np.linspace(-3.0, 3.0, 7)
    up = np.exp(-e_dos ** 2)
    down = 0.8 * np.exp(-(e_dos - 0.5) ** 2)
    dos_s = np.stack([up, down])                              # [2, 7]
    pdos_s = np.stack([np.stack([0.3 * up, 0.7 * up]),
                       np.stack([0.3 * down, 0.7 * down])])   # [2, 2, 7]
    cohp = np.stack([[-0.2 * up], [-0.2 * down]])             # [2, 1, 7]
    icohp = np.array([[-0.85], [-0.61]])                      # [2, 1]
    w = QvfWriter(_PROGRAM, _VERSION, "dos-spin")
    w.add_structure(**fe_bcc)
    w.add_dos_total(e_dos, dos_s, smearing=0.1, smearing_type="gaussian",
                    fermi_energy_ev=5.2, n_electrons=16.0, n_spin=2)
    w.add_dos_projected(e_dos, pdos_s,
                        channels=[{"atom": 0, "l": 0, "label": "Fe1-s"},
                                  {"atom": 0, "l": 2, "label": "Fe1-d"}],
                        n_spin=2, fermi_energy_ev=5.2)
    w.add_dos_coop_cohp("dos.cohp", e_dos, cohp, icohp, {
        "method": "cohp", "n_spin": 2, "fermi_energy_ev": 5.2,
        "smearing": 0.1, "smearing_type": "gaussian",
        "pairs": [{"i": 0, "j": 1, "label": "Fe1-Fe2", "distance_ang": 2.4855}],
        "pair_labels": ["Fe1-Fe2"]})
    _write("dos", w, [
        {"manifest_pointer": ["sections", 1, "members", "dos", "shape"],
         "equals": [2, 7]},
        {"section": "dos_total", "member": "dos", "binary_index": [1, 3],
         "equals": float(dos_s[1, 3]), "tol": 1e-12},
        {"manifest_pointer": ["sections", 2, "members", "projections", "shape"],
         "equals": [2, 2, 7]},
        {"section": "dos_projected", "member": "projections",
         "binary_index": [1, 1, 3], "equals": float(pdos_s[1, 1, 3]),
         "tol": 1e-12},
        {"manifest_pointer": ["sections", 2, "channels", 1, "label"],
         "equals": "Fe1-d"},
        {"section": "dos_cohp", "member": "meta",
         "json_pointer": ["n_spin"], "equals": 2},
        {"section": "dos_cohp", "member": "meta",
         "json_pointer": ["pairs", 0, "j"], "equals": 1},
        {"section": "dos_cohp", "member": "meta",
         "json_pointer": ["pair_labels", 0], "equals": "Fe1-Fe2"},
        {"section": "dos_cohp", "member": "meta",
         "json_pointer": ["smearing_type"], "equals": "gaussian"},
        {"section": "dos_cohp", "member": "projections",
         "binary_index": [1, 0, 3], "equals": float(cohp[1, 0, 3]),
         "tol": 1e-12},
        {"section": "dos_cohp", "member": "integrated", "binary_index": [1, 0],
         "equals": -0.61, "tol": 1e-12},
    ], "Spin-polarized DOS family: dos.total [n_spin, n_points], "
       "dos.projected [n_spin, n_channels, n_points] with channel objects, and "
       "dos.cohp [n_spin, n_pairs, n_points] whose meta payload names the "
       "pairs; energies in eV referenced to E_F.")

    # 14) Fermi surface: a 4x4x4 Gamma-centred mesh of E(k) - E_F (eV) for two
    # bands, with the cell, the reciprocal lattice (1/Angstrom, 2*pi included)
    # and the band indices in the mesh payload.
    a_po = 3.359
    b_po = 2.0 * np.pi / a_po
    frac = (np.arange(4) - 4 * (np.arange(4) >= 2)) / 4.0   # 0, 0.25, -0.5, -0.25
    kx, ky, kz = np.meshgrid(frac * b_po, frac * b_po, frac * b_po,
                             indexing="ij")
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    fs = np.stack([3.80998 * k2 - 3.0, 3.80998 * k2 - 0.5],
                  axis=-1).astype(np.float32)               # [4, 4, 4, 2]
    w = QvfWriter(_PROGRAM, _VERSION, "fermi-surface")
    w.add_structure([{"symbol": "Po", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 84}],
                    pbc=[True, True, True],
                    lattice_vectors=[[a_po, 0, 0], [0, a_po, 0], [0, 0, a_po]])
    w.add_fermi_surface({"shape": [4, 4, 4], "shift": [0.0, 0.0, 0.0],
                         "symmetry_expanded": True, "fermi_energy_ev": 6.25,
                         "n_spin": 1, "band_indices": [2, 3],
                         "lattice_vectors": [[a_po, 0, 0], [0, a_po, 0],
                                             [0, 0, a_po]],
                         "reciprocal_lattice": [[b_po, 0, 0], [0, b_po, 0],
                                                [0, 0, b_po]]}, fs)
    _write("fermi_surface", w, [
        {"section": "fermi_surface", "member": "mesh",
         "json_pointer": ["shape", 0], "equals": 4},
        {"section": "fermi_surface", "member": "mesh",
         "json_pointer": ["shift", 2], "equals": 0.0, "tol": 1e-12},
        {"section": "fermi_surface", "member": "mesh",
         "json_pointer": ["symmetry_expanded"], "equals": True},
        {"section": "fermi_surface", "member": "mesh",
         "json_pointer": ["reciprocal_lattice", 1, 1],
         "equals": b_po, "tol": 1e-12},
        {"section": "fermi_surface", "member": "mesh",
         "json_pointer": ["lattice_vectors", 2, 2], "equals": a_po,
         "tol": 1e-12},
        {"section": "fermi_surface", "member": "mesh",
         "json_pointer": ["band_indices", 1], "equals": 3},
        {"section": "fermi_surface", "member": "energies",
         "binary_index": [0, 0, 0, 0], "equals": -3.0, "tol": 1e-6},
        {"section": "fermi_surface", "member": "energies",
         "binary_index": [1, 2, 3, 1], "equals": float(fs[1, 2, 3, 1]),
         "tol": 1e-6},
    ], "Fermi surface: E(k) - E_F in eV on a 4x4x4 Gamma-centred mesh "
       "(float32 [nk1, nk2, nk3, n_bands]); the mesh payload carries shape, "
       "shift, the real-space cell (rows, Angstrom) and the reciprocal lattice "
       "(rows, 1/Angstrom, 2*pi included).")

    # 15) equation of state: E(V) points generated from a third-order
    # Birch-Murnaghan curve, with the fit parameters in the fit payload.
    v0, e0, b0_gpa, b0p = 74.7, -12.3, 160.0, 4.1
    b0_ev = b0_gpa / 160.21766208                             # eV / Angstrom^3
    vols = np.linspace(66.0, 84.0, 7)
    x = (v0 / vols) ** (2.0 / 3.0) - 1.0
    e_eos = e0 + 9.0 * v0 * b0_ev / 16.0 * (
        x ** 3 * b0p + x ** 2 * (6.0 - 4.0 * (v0 / vols) ** (2.0 / 3.0)))
    w = QvfWriter(_PROGRAM, _VERSION, "eos")
    w.add_structure(
        [{"symbol": "Mg", "position": [0.0, 0.0, 0.0], "atomic_number": 12},
         {"symbol": "O", "position": [2.1, 2.1, 2.1], "atomic_number": 8}],
        pbc=[True, True, True],
        lattice_vectors=[[4.2, 0, 0], [0, 4.2, 0], [0, 0, 4.2]])
    w.add_equation_of_state(vols, e_eos, {
        "model": "birch_murnaghan", "V0": v0, "E0": e0, "B0": b0_gpa,
        "B0_prime": b0p, "residual_rms": 0.0, "n_points": 7})
    _write("eos", w, [
        {"section": "eos", "member": "fit", "json_pointer": ["model"],
         "equals": "birch_murnaghan"},
        {"section": "eos", "member": "fit", "json_pointer": ["V0"],
         "equals": v0, "tol": 1e-12},
        {"section": "eos", "member": "fit", "json_pointer": ["B0"],
         "equals": b0_gpa, "tol": 1e-12},
        {"section": "eos", "member": "fit", "json_pointer": ["B0_prime"],
         "equals": b0p, "tol": 1e-12},
        {"section": "eos", "member": "fit", "json_pointer": ["n_points"],
         "equals": 7},
        {"section": "eos", "member": "volumes", "binary_index": [3],
         "equals": 75.0, "tol": 1e-12},
        {"section": "eos", "member": "energies", "binary_index": [3],
         "equals": float(e_eos[3]), "tol": 1e-12},
        {"section": "eos", "member": "energies", "binary_index": [0],
         "equals": float(e_eos[0]), "tol": 1e-12},
    ], "equation_of_state: volumes (Angstrom^3) and energies (eV) on a "
       "third-order Birch-Murnaghan curve; the fit payload carries model, V0, "
       "E0, B0 (GPa) and B0_prime.")

    # units_declared: declared units (spec § 4.3, "Declared units") on a grid,
    # a trajectory metadata, a UV-Vis spectrum stored in nm instead of the
    # kind's eV, and a QTAIM payload in bohr; provenance in the unit name.
    w = QvfWriter(_PROGRAM, _VERSION, "units-declared")
    w.add_structure([{"symbol": z, "position": list(map(float, xyz)),
                      "atomic_number": n}
                     for z, xyz, n in zip(h2o_symbols, base, [8, 1, 1])])
    w.set_provenance(method="RHF", basis="STO-3G", scf_energy_eh=-74.9659)
    rho = np.full((3, 3, 3), 0.05, dtype=np.float32)
    rho[1, 1, 1] = 0.35
    w.add_density(rho, [-1.0, -1.0, -1.0],
                  [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]],
                  units={"origin": "bohr", "voxel_vectors": "bohr"},
                  section_id="rho")
    w.add_trajectory({"symbols": h2o_symbols, "energies": opt_energies[:2],
                      "kind": "optimization"}, opt_coords[:2],
                     units={"energies": "hartree", "coords": "angstrom"},
                     section_id="opt")
    w.add_spectrum("spectra.uvvis", frequencies=[167.6, 128.3, 112.9],
                   intensities=[0.052, 0.110, 0.301], x_units="nm",
                   y_units="dimensionless", section_id="uvvis_nm")
    w.add_topology_qtaim({"points": [
        {"type": "ncp", "position": [0.0, 0.0, 0.2217], "atoms": [0]},
        {"type": "bcp", "position": [0.0, 1.0911, -0.6762], "rho": 0.3653,
         "atom_pair": [0, 1]}]}, units={"position": "bohr"}, section_id="aim")
    _write("units_declared", w, [
        {"manifest_pointer": ["provenance", "scf_energy", "units"],
         "equals": "hartree"},
        {"section": "rho", "member": "grid",
         "json_pointer": ["units", "origin"], "equals": "bohr"},
        {"section": "rho", "member": "grid",
         "json_pointer": ["units", "voxel_vectors"], "equals": "bohr"},
        {"section": "opt", "member": "metadata",
         "json_pointer": ["units", "energies"], "equals": "hartree"},
        {"section": "opt", "member": "metadata",
         "json_pointer": ["units", "coords"], "equals": "angstrom"},
        {"section": "uvvis_nm", "member": "spectrum",
         "json_pointer": ["x_units"], "equals": "nm"},
        {"section": "uvvis_nm", "member": "spectrum",
         "json_pointer": ["frequencies", 0], "equals": 167.6, "tol": 1e-9},
        {"section": "uvvis_nm", "member": "spectrum",
         "json_pointer": ["y_units"], "equals": "dimensionless"},
        {"section": "aim", "member": "critical_points",
         "json_pointer": ["units", "position"], "equals": "bohr"},
        {"section": "aim", "member": "critical_points",
         "json_pointer": ["points", 1, "position", 1], "equals": 1.0911,
         "tol": 1e-9},
    ], "Declared units (spec § 4.3): a volume grid restating bohr, a "
       "trajectory metadata restating hartree and angstrom, a spectra.uvvis "
       "whose x values are stored in nm (x_units) rather than the kind's eV, "
       "QTAIM positions in bohr, and provenance.scf_energy in the unit name "
       "'hartree'. Consumers read every value in its declared unit.")

    # 16) atom_properties beyond the first five roles (spec § 5.7): charges of
    # three more schemes, valences, the condensed Fukui functions and the dual
    # descriptor (consistent with each other), atomic energy and volume, the
    # isotropic NMR shielding, the two vector roles [n_atoms, 3] and one
    # custom.<name> role described by the `meta` member.
    f_plus = np.array([0.62, 0.38])
    f_minus = np.array([0.81, 0.19])
    w = QvfWriter(_PROGRAM, _VERSION, "oh-radical-extended")
    w.add_structure([{"symbol": "O", "position": [0, 0, 0], "atomic_number": 8},
                     {"symbol": "H", "position": [0, 0, 0.97], "atomic_number": 1}])
    w.add_atom_properties(
        cm5_charge=[-0.31, 0.31],
        chelpg_charge=[-0.38, 0.38],
        resp_charge=[-0.41, 0.41],
        valence=[1.92, 0.98],
        bonded_valence=[0.97, 0.97],
        fukui_plus=f_plus,
        fukui_minus=f_minus,
        fukui_zero=(f_plus + f_minus) / 2.0,
        dual_descriptor=f_plus - f_minus,
        atomic_energy_eh=[-75.18, -0.52],
        atomic_volume_bohr3=[64.3, 21.8],
        nmr_shielding_iso_ppm=[370.1, 29.6],
        magnetic_moment=[[0.0, 0.0, 0.93], [0.0, 0.0, 0.07]],
        force_eh_bohr=[[0.0, 0.0, 0.004], [0.0, 0.0, -0.004]],
        custom={"born_charge_mean": (
            [-0.52, 0.52],
            {"units": "e", "kind": "scalar",
             "description": "Trace of the Born effective-charge tensor divided "
                            "by three; illustrative values."})})
    _write("atom_properties_extended", w, [
        {"manifest_pointer": ["sections", 1, "members", "magnetic_moment",
                              "shape"], "equals": [2, 3]},
        {"manifest_pointer": ["sections", 1, "members", "force_eh_bohr",
                              "dtype"], "equals": "float64"},
        {"manifest_pointer": ["sections", 1, "members", "cm5_charge", "shape"],
         "equals": [2]},
        {"manifest_pointer": ["sections", 1, "members",
                              "custom.born_charge_mean", "shape"],
         "equals": [2]},
        {"section": "atom_properties", "member": "cm5_charge",
         "binary_index": [0], "equals": -0.31, "tol": 1e-12},
        {"section": "atom_properties", "member": "resp_charge",
         "binary_index": [1], "equals": 0.41, "tol": 1e-12},
        {"section": "atom_properties", "member": "valence",
         "binary_index": [0], "equals": 1.92, "tol": 1e-12},
        {"section": "atom_properties", "member": "bonded_valence",
         "binary_index": [1], "equals": 0.97, "tol": 1e-12},
        {"section": "atom_properties", "member": "fukui_plus",
         "binary_index": [0], "equals": 0.62, "tol": 1e-12},
        {"section": "atom_properties", "member": "fukui_zero",
         "binary_index": [0], "equals": 0.715, "tol": 1e-12},
        {"section": "atom_properties", "member": "dual_descriptor",
         "binary_index": [0], "equals": -0.19, "tol": 1e-12},
        {"section": "atom_properties", "member": "dual_descriptor",
         "binary_index": [1], "equals": 0.19, "tol": 1e-12},
        {"section": "atom_properties", "member": "atomic_energy_eh",
         "binary_index": [0], "equals": -75.18, "tol": 1e-12},
        {"section": "atom_properties", "member": "atomic_volume_bohr3",
         "binary_index": [1], "equals": 21.8, "tol": 1e-12},
        {"section": "atom_properties", "member": "nmr_shielding_iso_ppm",
         "binary_index": [0], "equals": 370.1, "tol": 1e-12},
        {"section": "atom_properties", "member": "magnetic_moment",
         "binary_index": [0, 2], "equals": 0.93, "tol": 1e-12},
        {"section": "atom_properties", "member": "magnetic_moment",
         "binary_index": [1, 0], "equals": 0.0, "tol": 1e-12},
        {"section": "atom_properties", "member": "force_eh_bohr",
         "binary_index": [1, 2], "equals": -0.004, "tol": 1e-12},
        {"section": "atom_properties", "member": "custom.born_charge_mean",
         "binary_index": [1], "equals": 0.52, "tol": 1e-12},
        {"section": "atom_properties", "member": "meta",
         "json_pointer": ["custom.born_charge_mean", "units"], "equals": "e"},
        {"section": "atom_properties", "member": "meta",
         "json_pointer": ["custom.born_charge_mean", "kind"],
         "equals": "scalar"},
    ], "atom_properties beyond the first five roles: three more charge "
       "schemes (cm5_charge, chelpg_charge, resp_charge), valence and "
       "bonded_valence, fukui_plus / fukui_minus / fukui_zero with "
       "dual_descriptor = fukui_plus - fukui_minus, atomic_energy_eh, "
       "atomic_volume_bohr3, nmr_shielding_iso_ppm, the vector roles "
       "magnetic_moment (collinear, on z) and force_eh_bohr "
       "(float64 [n_atoms, 3]), and the custom role custom.born_charge_mean "
       "described by the meta member; illustrative values for the OH radical.")

    # 16) two structures in one archive, told apart by section references
    # (spec § 4.1), with an explicit root metadata block (spec § 3.5). The
    # density describes the water molecule and the wavefunction the ammonia
    # molecule, so a consumer that takes the first or the only structure for
    # every section draws one of them on the wrong atoms.
    w = QvfWriter(_PROGRAM, _VERSION, "water-ammonia")
    h2o = w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}],
        section_id="h2o")
    nh3 = w.add_structure(
        [{"symbol": "N", "position": [0.0, 0.0, 0.1162], "atomic_number": 7},
         {"symbol": "H", "position": [0.0, 0.9397, -0.2711], "atomic_number": 1},
         {"symbol": "H", "position": [0.8138, -0.4699, -0.2711],
          "atomic_number": 1},
         {"symbol": "H", "position": [-0.8138, -0.4699, -0.2711],
          "atomic_number": 1}],
        section_id="nh3")
    run = w.add_run_record(
        program=_PROGRAM, program_version=_VERSION,
        input_text="water: rhf/sto-3g single point\n"
                   "ammonia: rhf/sto-3g single point\n",
        log_text="water: SCF converged\nammonia: SCF converged\n",
        section_id="run")
    pts = np.arange(4) - 1.5                                  # bohr
    xx, yy, zz = np.meshgrid(pts, pts, pts, indexing="ij")
    rho = w.add_density(np.exp(-(xx ** 2 + yy ** 2 + zz ** 2)),
                        [-1.5, -1.5, -1.5],
                        [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]],
                        label="water density", section_id="rho_h2o")
    shells = [{"center": c, "l": 0, "exponents": [e], "coefficients": [1.0]}
              for c, e in ((0, 0.8), (1, 0.4), (2, 0.4), (3, 0.4))]
    wf = w.add_wavefunction_gto(shells, mo_coefficients=np.eye(4),
                                energies=[-15.5, -1.1, -0.6, 0.6],
                                occupations=[2, 2, 0, 0], structure_ref=nh3,
                                section_id="wf_nh3")
    w.annotate_section(h2o, group="water")
    w.annotate_section(nh3, group="ammonia")
    w.annotate_section(rho, structure_ref=h2o, run_ref=run, group="water")
    w.annotate_section(wf, structure_ref=nh3, run_ref=run, group="ammonia")
    w.set_metadata(
        title="Water and ammonia: two structures in one archive",
        description="A density of water and a toy wavefunction of ammonia, "
                    "each naming its structure through structure_ref.",
        created_utc="2026-10-07T12:00:00Z",
        uuid="6f1c2a4e-8b3d-4c5e-9f70-1a2b3c4d5e6f",
        authors=[{"name": "QVF conformance suite"}],
        license="Apache-2.0",
        keywords=["conformance", "multi-structure", "section references"],
        tags=["test-vector"])
    _write("multi_structure", w, [
        {"manifest_pointer": ["metadata", "license"], "equals": "Apache-2.0"},
        {"manifest_pointer": ["metadata", "uuid"],
         "equals": "6f1c2a4e-8b3d-4c5e-9f70-1a2b3c4d5e6f"},
        {"manifest_pointer": ["metadata", "created_utc"],
         "equals": "2026-10-07T12:00:00Z"},
        {"manifest_pointer": ["metadata", "title"],
         "equals": "Water and ammonia: two structures in one archive"},
        {"manifest_pointer": ["metadata", "keywords", 1],
         "equals": "multi-structure"},
        {"manifest_pointer": ["metadata", "authors", 0, "name"],
         "equals": "QVF conformance suite"},
        {"manifest_pointer": ["sections", 0, "id"], "equals": "h2o"},
        {"manifest_pointer": ["sections", 1, "id"], "equals": "nh3"},
        {"manifest_pointer": ["sections", 2, "id"], "equals": "run"},
        {"manifest_pointer": ["sections", 3, "id"], "equals": "rho_h2o"},
        {"manifest_pointer": ["sections", 4, "id"], "equals": "wf_nh3"},
        {"manifest_pointer": ["sections", 0, "group"], "equals": "water"},
        {"manifest_pointer": ["sections", 1, "group"], "equals": "ammonia"},
        {"manifest_pointer": ["sections", 3, "structure_ref"], "equals": "h2o"},
        {"manifest_pointer": ["sections", 3, "run_ref"], "equals": "run"},
        {"manifest_pointer": ["sections", 3, "group"], "equals": "water"},
        {"manifest_pointer": ["sections", 4, "structure_ref"], "equals": "nh3"},
        {"manifest_pointer": ["sections", 4, "run_ref"], "equals": "run"},
        {"manifest_pointer": ["sections", 4, "group"], "equals": "ammonia"},
        {"section": "h2o", "member": "structure",
         "json_pointer": ["atoms", 0, "symbol"], "equals": "O"},
        {"section": "nh3", "member": "structure",
         "json_pointer": ["atoms", 0, "symbol"], "equals": "N"},
        {"section": "nh3", "member": "structure",
         "json_pointer": ["atoms", 3, "position", 0], "equals": -0.8138,
         "tol": 1e-9},
        {"section": "wf_nh3", "member": "basis",
         "json_pointer": ["structure_ref"], "equals": "nh3"},
        {"section": "wf_nh3", "member": "basis", "json_pointer": ["n_ao"],
         "equals": 4},
        {"section": "wf_nh3", "member": "mo_coefficients",
         "binary_index": [0, 0], "equals": 1.0, "tol": 1e-12},
        {"section": "rho_h2o", "member": "data", "binary_index": [1, 1, 1],
         "equals": float(np.exp(-0.75)), "tol": 1e-6},
        {"section": "run", "member": "log", "decode": "utf8",
         "contains": "ammonia: SCF converged"},
    ], "Two structures (water, ammonia) in one archive. The volume.density "
       "and the wavefunction.gto each name their structure with a "
       "section-level structure_ref and the run.record that produced them "
       "with run_ref; group labels pair each structure with its data. The "
       "root metadata block carries a title, a fixed created_utc and uuid, "
       "the Apache-2.0 licence and keywords (spec § 3.5, § 4.1).")

    # ==== WP-14: conformance corpus v2 (begin) ================================
    _wp14_coverage_archives()
    _wp14_edge_case_archives()
    # ==== WP-14: conformance corpus v2 (end) ==================================
    _forces_hessian()
    _provenance_job_spec()
    _container_rules()
    _md_periodic()
    # ==== WP-22: excited states (begin) =======================================
    _tddft()
    # ==== WP-22: excited states (end) =========================================
    _wavefunction_extensions()
    _wp25_volume_archives()
    # ==== WP-29: property kinds (begin) =======================================
    _wp29_properties()
    # ==== WP-29: property kinds (end) =========================================
    _bloch()
    _sketch_documents()
    _family_drift_archives()


# CODATA 2018: hartree in J, bohr in m, unified atomic mass unit in kg, the
# speed of light in cm/s, and hartree in cm^-1.
_HARTREE_J = 4.3597447222071e-18
_BOHR_M = 5.29177210903e-11
_AMU_KG = 1.66053906660e-27
_C_CM_S = 2.99792458e10
_HARTREE_CM = 219474.6313632


def _forces_hessian() -> None:
    """16) forces + hessian + root energies (spec § 3.6, § 5.10, § 5.11).

    Water with a synthetic but internally consistent derivative set: the
    Hessian (hartree/bohr^2) of three harmonic springs (two O-H stretches and
    an H-H spring standing in for the bend), the gradient (hartree/bohr) of
    those springs with both O-H bonds stretched 0.004 bohr past their rest
    length, harmonic frequencies from the mass-weighted Hessian, and a root
    energies block whose RHF total restates provenance.scf_energy and whose
    zpve restates thermochemistry.zpve_eh.
    """
    atoms = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
             {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
             {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]
    xyz = np.array([a["position"] for a in atoms]) * BOHR_PER_ANGSTROM
    charges = [8.0, 1.0, 1.0]
    masses = [15.999, 1.008, 1.008]
    springs = [(0, 1, 0.52, 0.004), (0, 2, 0.52, 0.004), (1, 2, 0.06, 0.0)]
    hess = np.zeros((9, 9))
    grad = np.zeros((3, 3))
    for a, b, k, stretch in springs:
        u = xyz[b] - xyz[a]
        u /= np.linalg.norm(u)
        block = k * np.outer(u, u)
        hess[3 * a:3 * a + 3, 3 * a:3 * a + 3] += block
        hess[3 * b:3 * b + 3, 3 * b:3 * b + 3] += block
        hess[3 * a:3 * a + 3, 3 * b:3 * b + 3] -= block
        hess[3 * b:3 * b + 3, 3 * a:3 * a + 3] -= block
        grad[a] -= k * stretch * u          # dE/dR_a = -k (r - r0) u
        grad[b] += k * stretch * u
    # Harmonic frequencies: the three non-zero eigenvalues of the mass-weighted
    # Hessian (the six rigid-body ones vanish for a spring model).
    inv_sqrt_m = np.repeat(1.0 / np.sqrt(masses), 3)
    eigval, eigvec = np.linalg.eigh(hess * np.outer(inv_sqrt_m, inv_sqrt_m))
    to_cm = np.sqrt(_HARTREE_J / (_BOHR_M ** 2 * _AMU_KG)) / (2.0 * np.pi * _C_CM_S)
    freqs = [float(to_cm * np.sqrt(lam)) for lam in eigval[-3:]]
    modes = []
    for q in eigvec[:, -3:].T:
        d = (q * inv_sqrt_m).reshape(3, 3)
        d /= np.linalg.norm(d)
        modes.append(d if d.flat[np.argmax(np.abs(d))] > 0 else -d)
    zpve = 0.5 * sum(freqs) / _HARTREE_CM
    e_nn = sum(charges[i] * charges[j] / float(np.linalg.norm(xyz[i] - xyz[j]))
               for i in range(3) for j in range(i + 1, 3))
    e_scf = -74.9659011

    units = {"energy": "hartree", "length": "bohr"}
    w = QvfWriter(_PROGRAM, _VERSION, "water-derivatives")
    w.add_structure(atoms)
    w.add_forces(grad, meta={"structure_ref": "structure", "method": "rhf",
                             "convention": "gradient", "units": units})
    w.add_hessian(hess, meta={"structure_ref": "structure", "mass_weighted": False,
                              "projected": False, "atom_masses_amu": masses,
                              "units": units})
    w.add_vibrations({"symbols": ["O", "H", "H"], "frequencies": freqs,
                      "n_modes": 3}, np.array(modes), hessian_ref="hessian")
    w.set_provenance(method="rhf", basis="sto-3g", scf_converged=True,
                     scf_energy={"value": e_scf, "units": "hartree"})
    w.set_thermochemistry(zpve_eh=zpve, temperature_k=298.15)
    w.set_energies(e_scf, method="rhf",
                   components={"nuclear_repulsion": e_nn,
                               "electronic": e_scf - e_nn, "zpve": zpve},
                   by_method=[{"method": "rhf", "total": e_scf}])
    _write("forces_hessian", w, [
        {"manifest_pointer": ["sections", 1, "members", "gradient", "shape"],
         "equals": [3, 3]},
        {"section": "forces", "member": "gradient", "binary_index": [1, 1],
         "equals": float(grad[1, 1]), "tol": 1e-12},
        {"section": "forces", "member": "gradient", "binary_index": [0, 2],
         "equals": float(grad[0, 2]), "tol": 1e-12},
        {"section": "forces", "member": "meta", "json_pointer": ["convention"],
         "equals": "gradient"},
        {"section": "forces", "member": "meta", "json_pointer": ["structure_ref"],
         "equals": "structure"},
        {"manifest_pointer": ["sections", 2, "members", "matrix", "shape"],
         "equals": [9, 9]},
        {"section": "hessian", "member": "matrix", "binary_index": [4, 4],
         "equals": float(hess[4, 4]), "tol": 1e-12},
        {"section": "hessian", "member": "matrix", "binary_index": [1, 4],
         "equals": float(hess[1, 4]), "tol": 1e-12},
        {"section": "hessian", "member": "matrix", "binary_index": [4, 1],
         "equals": float(hess[1, 4]), "tol": 1e-12},
        {"section": "hessian", "member": "meta",
         "json_pointer": ["atom_masses_amu", 0], "equals": 15.999, "tol": 1e-12},
        {"section": "hessian", "member": "meta", "json_pointer": ["mass_weighted"],
         "equals": False},
        {"manifest_pointer": ["sections", 3, "hessian_ref"], "equals": "hessian"},
        {"section": "vibrations", "member": "metadata",
         "json_pointer": ["frequencies", 2], "equals": freqs[2], "tol": 1e-6},
        {"manifest_pointer": ["energies", "total", "value"], "equals": e_scf,
         "tol": 1e-12},
        {"manifest_pointer": ["energies", "total", "units"], "equals": "hartree"},
        {"manifest_pointer": ["energies", "method"], "equals": "rhf"},
        {"manifest_pointer": ["energies", "by_method", 0, "total", "value"],
         "equals": e_scf, "tol": 1e-12},
        {"manifest_pointer": ["energies", "components", "nuclear_repulsion",
                              "value"], "equals": e_nn, "tol": 1e-9},
        {"manifest_pointer": ["provenance", "scf_energy", "value"],
         "equals": e_scf, "tol": 1e-12},
        {"manifest_pointer": ["thermochemistry", "zpve_eh"], "equals": zpve,
         "tol": 1e-12},
        {"manifest_pointer": ["energies", "components", "zpve", "value"],
         "equals": zpve, "tol": 1e-12},
    ], "Water with energy derivatives: forces.gradient is the energy gradient "
       "dE/dR (hartree/bohr; forces are its negative) as float64 [n_atoms, 3] "
       "with no stress; hessian.matrix is the [3N, 3N] Cartesian Hessian "
       "(hartree/bohr^2, index 3*atom + xyz, symmetric) with atom masses in "
       "its meta; the vibrations section names it through hessian_ref. The "
       "root energies block restates provenance.scf_energy (RHF) and "
       "thermochemistry.zpve_eh, as spec § 3.6 requires.")


# =============================================================================
# WP-22: excited states. tddft.qvf: water with a ground-state wavefunction, an
# excited_states table, the natural transition orbitals and the transition
# density of its first state, and the UV-Vis stick spectrum the table implies.
# =============================================================================

_HARTREE_EV = 27.211386245988      # CODATA 2018


def _tddft() -> None:
    """17) excited_states + NTOs + volume.transition_density + spectra.uvvis
    (spec § 5.2, § 5.3, § 5.4, § 5.12).

    Water (C2v, molecular plane yz, x out of plane) with an STO-3G basis and
    a model set of seven symmetry-adapted orbitals (rows of an orthonormal
    matrix: O 1s, O 2s, O 2p_y, O 2p_z, O 2p_x, and the in- and out-of-phase
    H 1s pairs), five of them doubly occupied. Three synthetic but internally
    consistent singlet states: HOMO (1b1, MO 4) -> 4a1 (MO 5) is B1 and
    x-polarized; HOMO -> 2b2 (MO 6) is A2 and dipole-forbidden (f = 0, its
    magnetic transition dipole along z); 3a1 -> 4a1 is A1 and z-polarized.
    Each transition dipole has the length that f = (2/3) dE |mu|^2 gives,
    the molecule is achiral, so every rotatory strength is zero. State 1 has
    a hole/particle NTO pair and a transition density on a grid; the UV-Vis
    spectrum has one line per state at its excitation energy with the
    oscillator strength as intensity, and names the table through states_ref.
    """
    atoms = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
             {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
             {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]
    # STO-3G (contraction coefficients on normalized primitives).
    o_s1 = ([130.70932, 23.808861, 6.4436083], [0.15432897, 0.53532814, 0.44463454])
    o_sp = [5.0331513, 1.1695961, 0.380389]
    o_s2 = [-0.09996723, 0.39951283, 0.70011547]
    o_p = [0.15591627, 0.60768372, 0.39195739]
    h_s = ([3.42525091, 0.62391373, 0.1688554], [0.15432897, 0.53532814, 0.44463454])
    shells = [{"center": 0, "l": 0, "exponents": o_s1[0], "coefficients": o_s1[1]},
              {"center": 0, "l": 0, "exponents": o_sp, "coefficients": o_s2},
              {"center": 0, "l": 1, "exponents": o_sp, "coefficients": o_p},
              {"center": 1, "l": 0, "exponents": h_s[0], "coefficients": h_s[1]},
              {"center": 2, "l": 0, "exponents": h_s[0], "coefficients": h_s[1]}]
    # AO order: O 1s, O 2s, O 2p (m = -1, 0, +1: y, z, x), H_a 1s, H_b 1s.
    r = 1.0 / math.sqrt(2.0)
    mos = np.zeros((7, 7))
    for mo, ao in enumerate((0, 1, 2, 3, 4)):
        mos[mo, ao] = 1.0          # 1a1, 2a1, 1b2 (y), 3a1 (z), 1b1 (x)
    mos[5, 5], mos[5, 6] = r, r    # 4a1: H_a + H_b
    mos[6, 5], mos[6, 6] = r, -r   # 2b2: H_a - H_b
    eps = [-19.0302, -1.0159, -0.4999, -0.3243, -0.2531, 0.2251, 0.3046]
    e_ground = -75.3129834

    def dipole_length(f, de_ev):
        return math.sqrt(1.5 * f / (de_ev / _HARTREE_EV))

    de = [7.40, 9.20, 9.85]
    f = [0.0521, 0.0, 0.0846]
    mu1 = dipole_length(f[0], de[0])
    mu3 = dipole_length(f[2], de[2])
    states = [
        {"index": 1, "excitation_energy": {"value": de[0], "units": "ev"},
         "total_energy": {"value": e_ground + de[0] / _HARTREE_EV,
                          "units": "hartree"},
         "multiplicity": 1, "symmetry": "B1", "oscillator_strength": f[0],
         "rotatory_strength_cgs": 0.0,
         "transition_dipole": [mu1, 0.0, 0.0],
         "transition_magnetic_dipole": [0.0, 0.2134, 0.0],
         "contributions": [{"from_mo": 4, "to_mo": 5, "weight": 0.99}],
         "nto_ref": "nto_s1_hole", "density_ref": "tden_s1"},
        {"index": 2, "excitation_energy": {"value": de[1], "units": "ev"},
         "total_energy": {"value": e_ground + de[1] / _HARTREE_EV,
                          "units": "hartree"},
         "multiplicity": 1, "symmetry": "A2", "oscillator_strength": f[1],
         "rotatory_strength_cgs": 0.0,
         "transition_dipole": [0.0, 0.0, 0.0],
         "transition_magnetic_dipole": [0.0, 0.0, 0.3517],
         "contributions": [{"from_mo": 4, "to_mo": 6, "weight": 0.98}]},
        {"index": 3, "excitation_energy": {"value": de[2], "units": "ev"},
         "total_energy": {"value": e_ground + de[2] / _HARTREE_EV,
                          "units": "hartree"},
         "multiplicity": 1, "symmetry": "A1", "oscillator_strength": f[2],
         "rotatory_strength_cgs": 0.0,
         "transition_dipole": [0.0, 0.0, mu3],
         "transition_magnetic_dipole": [0.0, 0.0, 0.0],
         "contributions": [{"from_mo": 3, "to_mo": 5, "weight": 0.90},
                           {"from_mo": 2, "to_mo": 6, "weight": 0.08}]},
    ]

    # Transition density of state 1 on a grid: a model of the 1b1 -> 4a1 pair,
    # x exp(-|r - r_O|^2) (exp(-|r - r_Ha|^2) + exp(-|r - r_Hb|^2)), odd in x
    # (so it integrates to zero and its first moment lies along x, as mu does).
    xyz = np.array([a["position"] for a in atoms]) * BOHR_PER_ANGSTROM
    shape, step = (8, 10, 8), 0.5
    origin = [-1.75, -2.25, -2.25]
    ii, jj, kk = np.meshgrid(*(np.arange(n) for n in shape), indexing="ij")
    pts = np.stack([origin[0] + step * ii, origin[1] + step * jj,
                    origin[2] + step * kk], axis=-1)
    def gauss(center):
        return np.exp(-((pts - center) ** 2).sum(axis=-1))

    tden = (pts[..., 0] * gauss(xyz[0]) * (gauss(xyz[1]) + gauss(xyz[2]))
            ).astype(np.float32)
    voxels = [[step, 0, 0], [0, step, 0], [0, 0, step]]

    w = QvfWriter(_PROGRAM, _VERSION, "water-tddft")
    w.add_structure(atoms)
    w.set_provenance(method="rks", functional="b3lyp", basis="sto-3g",
                     scf_converged=True,
                     scf_energy={"value": e_ground, "units": "hartree"})
    w.add_wavefunction_gto(shells, mo_coefficients=mos, energies=eps,
                           occupations=[2, 2, 2, 2, 2, 0, 0], section_id="wf")
    # NTOs of state 1: two hole/particle pairs with weights 0.99 and 0.01.
    hole = np.array([mos[4], mos[2]])
    particle = np.array([mos[5], mos[6]])
    for role, coeffs in (("hole", hole), ("particle", particle)):
        w.add_wavefunction_gto(shells, mo_coefficients=coeffs,
                               occupations=[0.99, 0.01], orbital_kind="nto",
                               state_ref="excited_states", state_index=1,
                               nto_role=role, section_id=f"nto_s1_{role}")
    w.add_volume("volume.transition_density", tden, origin, voxels,
                 label="S1 transition density", excited_states_ref="excited_states",
                 state_index=1, section_id="tden_s1")
    w.add_excited_states(states, meta={
        "method": "tddft", "reference_state": "S0", "n_states": 3, "tda": False,
        "spin_flip": False,
        "units": {"excitation_energy": "ev", "total_energy": "hartree",
                  "oscillator_strength": "dimensionless"}},
        wavefunction_ref="wf")
    w.annotate_section("excited_states", structure_ref="structure")
    labels = [f"S{s['index']} ({s['symmetry']})" for s in states]
    w.add_spectrum("spectra.uvvis", payload={
        "frequencies": de, "intensities": f, "x_units": "ev",
        "y_units": "dimensionless",
        "peaks": [{"x": de[k], "y": f[k], "label": labels[k],
                   "assignment": {"state_index": states[k]["index"]}}
                  for k in range(3)],
        "broadening": {"type": "gaussian", "fwhm": 0.3, "units": "ev"},
        "uvvis": {"x_axis": "energy_ev"}, "states_ref": "excited_states"},
        label="UV-Vis (TD-DFT sticks)")
    probe = [(3, 5, 4), (5, 2, 3), (0, 7, 2)]
    _write("tddft", w, [
        {"manifest_pointer": ["sections", 5, "kind"], "equals": "excited_states"},
        {"manifest_pointer": ["sections", 5, "wavefunction_ref"], "equals": "wf"},
        {"manifest_pointer": ["sections", 5, "structure_ref"],
         "equals": "structure"},
        {"manifest_pointer": ["sections", 4, "excited_states_ref"],
         "equals": "excited_states"},
        {"manifest_pointer": ["sections", 4, "state_index"], "equals": 1},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["meta", "n_states"], "equals": 3},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["meta", "method"], "equals": "tddft"},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "index"], "equals": 1},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "excitation_energy", "value"],
         "equals": de[0], "tol": 1e-12},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "excitation_energy", "units"],
         "equals": "ev"},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 2, "total_energy", "value"],
         "equals": e_ground + de[2] / _HARTREE_EV, "tol": 1e-10},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "oscillator_strength"], "equals": f[0],
         "tol": 1e-12},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "transition_dipole", 0], "equals": mu1,
         "tol": 1e-12},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 2, "transition_dipole", 2], "equals": mu3,
         "tol": 1e-12},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 1, "oscillator_strength"], "equals": 0.0,
         "tol": 1e-12},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 1, "transition_magnetic_dipole", 2],
         "equals": 0.3517, "tol": 1e-12},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 1, "symmetry"], "equals": "A2"},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "contributions", 0],
         "equals": {"from_mo": 4, "to_mo": 5, "weight": 0.99}},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 2, "contributions", 1, "to_mo"], "equals": 6},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "nto_ref"], "equals": "nto_s1_hole"},
        {"section": "excited_states", "member": "states",
         "json_pointer": ["states", 0, "density_ref"], "equals": "tden_s1"},
        {"section": "wf", "member": "mo_coefficients", "binary_index": [6, 6],
         "equals": -r, "tol": 1e-12},
        {"section": "nto_s1_hole", "member": "mo_metadata",
         "json_pointer": ["orbital_kind"], "equals": "nto"},
        {"section": "nto_s1_hole", "member": "mo_metadata",
         "json_pointer": ["nto_role"], "equals": "hole"},
        {"section": "nto_s1_particle", "member": "mo_metadata",
         "json_pointer": ["nto_role"], "equals": "particle"},
        {"section": "nto_s1_particle", "member": "mo_metadata",
         "json_pointer": ["state_ref"], "equals": "excited_states"},
        {"section": "nto_s1_particle", "member": "mo_metadata",
         "json_pointer": ["state_index"], "equals": 1},
        {"section": "nto_s1_particle", "member": "mo_metadata",
         "json_pointer": ["occupation_semantics"], "equals": "transition_weight"},
        {"section": "nto_s1_hole", "member": "mo_metadata",
         "json_pointer": ["occupations", 0], "equals": 0.99, "tol": 1e-12},
        {"section": "nto_s1_hole", "member": "mo_coefficients",
         "binary_index": [0, 4], "equals": 1.0, "tol": 1e-12},
        {"section": "nto_s1_particle", "member": "mo_coefficients",
         "binary_index": [0, 6], "equals": r, "tol": 1e-12},
        *({"section": "tden_s1", "member": "data", "binary_index": list(p),
           "equals": float(tden[p]), "tol": 1e-6} for p in probe),
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["states_ref"], "equals": "excited_states"},
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["frequencies"], "equals": de},
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["intensities", 2], "equals": f[2], "tol": 1e-12},
        {"section": "spectra_uvvis", "member": "spectrum",
         "json_pointer": ["peaks", 2, "assignment", "state_index"], "equals": 3},
    ], "Water TD-DFT: an excited_states table of three singlet states (B1, A2, "
       "A1) with excitation energies in eV and total energies in hartree, "
       "oscillator strengths whose transition dipoles (atomic units) satisfy "
       "f = (2/3) dE |mu|^2, zero rotatory strengths, magnetic transition "
       "dipoles, and orbital contributions indexing the 7-orbital ground-state "
       "wavefunction that wavefunction_ref names. State 1 names its hole NTO "
       "section (nto_ref; the particle section carries the same state_ref and "
       "state_index with nto_role particle, occupations are the pair weights) "
       "and its volume.transition_density (density_ref; the volume names the "
       "state back with excited_states_ref and state_index). The spectra.uvvis "
       "stick spectrum has one line per state at its excitation energy with "
       "the oscillator strength as intensity, and its peaks name states by "
       "state_index through states_ref (spec § 5.2, § 5.3, § 5.4, § 5.12).")


def _provenance_job_spec() -> None:
    """17) the typed provenance fields and the grown job.spec (WP-27; spec
    § 3.2, § 5.9). ``job_spec.qvf`` stays as it was: these are new archives."""
    # job_spec_full: a pending periodic job that carries every optional
    # field of the job.spec payload. The archive holds two structures, so
    # `structure_ref` is what says which one the job runs on (the slab, the
    # second one). The provenance of a pending archive describes the step that
    # prepared the job: the program that wrote the container, its revision,
    # the resources of that step and when it ran.
    w = QvfWriter(_PROGRAM, _VERSION, "pending-job-full")
    a = 4.212
    w.add_structure(
        [{"symbol": "Mg", "position": [0.0, 0.0, 0.0], "atomic_number": 12},
         {"symbol": "O", "position": [a / 2, a / 2, a / 2], "atomic_number": 8}],
        pbc=[True, True, True],
        lattice_vectors=[[a, 0, 0], [0, a, 0], [0, 0, a]],
        section_id="mgo_reference")
    s = 4.3
    w.add_structure(
        [{"symbol": "Mg", "position": [0.0, 0.0, 0.0], "atomic_number": 12},
         {"symbol": "O", "position": [s / 2, s / 2, 0.0], "atomic_number": 8},
         {"symbol": "O", "position": [0.0, 0.0, s / 2], "atomic_number": 8},
         {"symbol": "Mg", "position": [s / 2, s / 2, s / 2],
          "atomic_number": 12}],
        pbc=[True, True, False],
        lattice_vectors=[[s, 0, 0], [0, s, 0], [0, 0, 20.0]],
        section_id="mgo_slab")
    w.add_job_spec(
        job_type="periodic", method="rks", functional="pbe",
        basis="pob-tzvp-rev2", charge=0, multiplicity=1, kpoints=[4, 4, 1],
        tasks=["optimize"], structure_ref="mgo_slab",
        solvation={"model": "cpcm", "solvent": "water", "epsilon": 78.4},
        dispersion="d3bj", relativistic="x2c",
        scf={"convergence": 1e-8, "max_iterations": 100, "guess": "sad",
             "damping": 0.3},
        geometry={"convergence": 4.5e-4, "max_steps": 50,
                  "coordinates": "redundant"},
        resources={"n_cores": 16, "memory_mb": 32000, "walltime_s": 86400},
        options={"cutoff_ha": 300.0})
    w.set_provenance(
        run_status="pending",
        software=[{"name": "job-builder", "version": "0.3.0",
                   "url": "https://example.org/job-builder"},
                  {"name": "integral-lib", "version": "2.9.0"}],
        git_commit="3e8a1c5f7b2d", n_threads=8, n_mpi_ranks=4,
        started_utc="2026-10-07T09:00:00Z",
        finished_utc="2026-10-07T09:00:02Z")

    def job(*pointer):
        return {"section": "job_spec", "member": "spec",
                "json_pointer": list(pointer)}

    _write("job_spec_full", w, [
        {"manifest_pointer": ["provenance", "run_status"], "equals": "pending"},
        {"manifest_pointer": ["provenance", "software", 0, "name"],
         "equals": "job-builder"},
        {"manifest_pointer": ["provenance", "software", 0, "url"],
         "equals": "https://example.org/job-builder"},
        {"manifest_pointer": ["provenance", "software", 1, "version"],
         "equals": "2.9.0"},
        {"manifest_pointer": ["provenance", "git_commit"],
         "equals": "3e8a1c5f7b2d"},
        {"manifest_pointer": ["provenance", "n_threads"], "equals": 8},
        {"manifest_pointer": ["provenance", "n_mpi_ranks"], "equals": 4},
        {"manifest_pointer": ["provenance", "started_utc"],
         "equals": "2026-10-07T09:00:00Z"},
        {"manifest_pointer": ["provenance", "finished_utc"],
         "equals": "2026-10-07T09:00:02Z"},
        {**job("job_type"), "equals": "periodic"},
        {**job("structure_ref"), "equals": "mgo_slab"},
        {"section": "mgo_slab", "member": "structure",
         "json_pointer": ["lattice_vectors", 0, 0], "equals": s, "tol": 1e-12},
        {"section": "mgo_slab", "member": "structure",
         "json_pointer": ["pbc", 2], "equals": False},
        {"section": "mgo_reference", "member": "structure",
         "json_pointer": ["lattice_vectors", 0, 0], "equals": a, "tol": 1e-12},
        {**job("solvation", "model"), "equals": "cpcm"},
        {**job("solvation", "solvent"), "equals": "water"},
        {**job("solvation", "epsilon"), "equals": 78.4, "tol": 1e-12},
        {**job("dispersion"), "equals": "d3bj"},
        {**job("relativistic"), "equals": "x2c"},
        {**job("scf", "convergence"), "equals": 1e-8, "tol": 1e-20},
        {**job("scf", "max_iterations"), "equals": 100},
        {**job("scf", "guess"), "equals": "sad"},
        {**job("scf", "damping"), "equals": 0.3, "tol": 1e-12},
        {**job("geometry", "convergence"), "equals": 4.5e-4, "tol": 1e-15},
        {**job("geometry", "max_steps"), "equals": 50},
        {**job("geometry", "coordinates"), "equals": "redundant"},
        {**job("resources", "n_cores"), "equals": 16},
        {**job("resources", "memory_mb"), "equals": 32000},
        {**job("resources", "walltime_s"), "equals": 86400.0, "tol": 1e-9},
        {**job("kpoints", 2), "equals": 1},
        {**job("tasks", 0), "equals": "optimize"},
        {**job("options", "cutoff_ha"), "equals": 300.0, "tol": 1e-9},
    ], "job.spec with every optional field of the payload: a periodic job "
       "whose structure_ref names the second of two structures (the slab, "
       "not the bulk reference), with solvation, dispersion, relativistic, "
       "scf, geometry and resources blocks (resources are hints). The "
       "provenance of this pending archive describes the step that prepared "
       "it: software (with a url), git_commit, thread and rank counts, and "
       "the two time stamps (spec § 3.2, § 5.9).")

    # provenance_full: a settled (converged) archive whose provenance uses
    # every typed field, beside the older ones. The job.spec that asked for the
    # run is preserved, and its resource hint differs from what the run used.
    atoms = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
             {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
             {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]
    w = QvfWriter(_PROGRAM, _VERSION, "settled-run")
    w.add_structure(atoms)
    w.add_job_spec(job_type="molecular", method="rks", functional="pbe",
                   basis="def2-svp", tasks=["single_point"],
                   structure_ref="structure", dispersion="d3bj",
                   resources={"n_cores": 64, "walltime_s": 3600.0})
    w.add_run_record(program="my-code", program_version="1.2.3",
                     command="my-code water.inp", exit_status=0,
                     started_utc="2026-10-07T12:00:00Z",
                     finished_utc="2026-10-07T12:07:31.25Z", sequence=0,
                     log_text="my-code 1.2.3\nSCF converged in 9 iterations\n")
    w.set_provenance(
        method="rks", functional="pbe", basis="def2-svp", charge=0,
        multiplicity=1, n_electrons=10, n_scf_iterations=9, scf_converged=True,
        scf_energy={"value": -76.3, "units": "hartree"}, wall_seconds=451.25,
        run_status="converged",
        software=[{"name": "my-code", "version": "1.2.3"},
                  {"name": "integral-lib", "version": "2.9.0",
                   "url": "https://example.org/integral-lib"},
                  {"name": "xc-lib", "version": "6.2.2", "x_build": "release"}],
        git_commit="3e8a1c5f7b2d49608a1b3c5d7e9f0a2b4c6d8e1f",
        n_threads=8, n_mpi_ranks=4, n_gpus=2,
        started_utc="2026-10-07T12:00:00Z",
        finished_utc="2026-10-07T12:07:31.25Z")
    _write("provenance_full", w, [
        {"manifest_pointer": ["provenance", "run_status"], "equals": "converged"},
        {"manifest_pointer": ["provenance", "software", 0, "name"],
         "equals": "my-code"},
        {"manifest_pointer": ["provenance", "software", 1, "url"],
         "equals": "https://example.org/integral-lib"},
        {"manifest_pointer": ["provenance", "software", 2, "x_build"],
         "equals": "release"},
        {"manifest_pointer": ["provenance", "git_commit"],
         "equals": "3e8a1c5f7b2d49608a1b3c5d7e9f0a2b4c6d8e1f"},
        {"manifest_pointer": ["provenance", "n_threads"], "equals": 8},
        {"manifest_pointer": ["provenance", "n_mpi_ranks"], "equals": 4},
        {"manifest_pointer": ["provenance", "n_gpus"], "equals": 2},
        {"manifest_pointer": ["provenance", "started_utc"],
         "equals": "2026-10-07T12:00:00Z"},
        {"manifest_pointer": ["provenance", "finished_utc"],
         "equals": "2026-10-07T12:07:31.25Z"},
        {"manifest_pointer": ["provenance", "wall_seconds"], "equals": 451.25,
         "tol": 1e-12},
        {"manifest_pointer": ["provenance", "n_electrons"], "equals": 10},
        {"manifest_pointer": ["provenance", "scf_energy", "value"],
         "equals": -76.3, "tol": 1e-12},
        {"manifest_pointer": ["sections", 2, "started_utc"],
         "equals": "2026-10-07T12:00:00Z"},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["structure_ref"], "equals": "structure"},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["resources", "n_cores"], "equals": 64},
        {"section": "run_record", "member": "log", "decode": "utf8",
         "contains": "SCF converged in 9 iterations"},
    ], "A settled (converged) archive whose provenance uses every typed "
       "field: software (one item with a url, one with a vendor key), a "
       "40-character git_commit, n_threads, n_mpi_ranks, n_gpus, and the two "
       "UTC time stamps (one with a fractional second) that agree with "
       "wall_seconds and with the run.record. The job.spec that asked for the "
       "run is preserved; its n_cores is a hint and differs from what the "
       "run used (spec § 3.2, § 5.9).")

# --- container rules: member size, byte_order, updating in place -------------
# (spec § 2.2, § 2.3, § 4.2)
import warnings  # noqa: E402

_WATER_ATOMS = [
    {"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
    {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
    {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]


def _container_rules() -> None:
    """17) size_declared and 18) updated_in_place."""
    _size_declared()
    _updated_in_place()


def _size_declared() -> None:
    """17) A `size` on every kind of member spec, and the root `byte_order`."""
    bibtex = ("@article{sizes2026,\n"
              "  author = {Ångström, A. and Müller, B.},\n"
              "  title = {Member sizes},\n"
              "  year = {2026}\n"
              "}\n")
    input_text = "! RHF STO-3G\n* xyz 0 1\nO 0 0 0.1173\nH 0 0.7572 -0.4692\n" \
                 "H 0 -0.7572 -0.4692\n*\n"
    log_text = "SCF CONVERGED\nFINAL SINGLE POINT ENERGY  -74.965901\n"
    geom = b"3\nwater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692\n"
    density = np.arange(24, dtype=np.float32).reshape(2, 3, 4) / 24.0
    counts = np.array([8, 1, 1], dtype=np.int64)

    w = QvfWriter(_PROGRAM, _VERSION, "member-sizes")
    w.add_structure(_WATER_ATOMS, bonds=[{"i": 0, "j": 1, "order": 1.0},
                                        {"i": 0, "j": 2, "order": 1.0}])
    w.add_volume("volume.density", density, [-1.0, -1.0, -1.0],
                 [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]],
                 section_id="rho")
    w.add_citations(bibtex)
    w.add_run_record(program="democode", input_text=input_text,
                     log_text=log_text,
                     files={"input": {"filename": "water.inp"},
                            "log": {"filename": "water.out"},
                            "attachment.geom": {"filename": "water.xyz"}},
                     attachments={"geom": geom}, section_id="run")
    w.add_vendor_section("x_conformance.sizes",
                         binary_members={"counts": (counts, "int64")},
                         section_id="x_sizes")
    manifest = w.build_manifest()
    sizes = {s["id"]: {role: m["size"] for role, m in s["members"].items()}
             for s in manifest["sections"]}
    # Binary and opaque members have producer-independent bytes, so their
    # sizes are checked exactly; a JSON member's size depends on how its
    # producer formats the JSON, so only its presence is.
    assert sizes["rho"]["data"] == density.nbytes == 96
    assert sizes["citations"]["references"] == len(bibtex.encode("utf-8"))
    assert sizes["x_sizes"]["counts"] == 24

    def size_of(i, role):
        return ["sections", i, "members", role, "size"]

    _write("size_declared", w, [
        {"manifest_pointer": ["byte_order"], "equals": "little"},
        {"manifest_pointer": size_of(1, "data"), "equals": 96},
        {"manifest_pointer": size_of(2, "references"),
         "equals": len(bibtex.encode("utf-8"))},
        {"manifest_pointer": size_of(3, "input"),
         "equals": len(input_text.encode("utf-8"))},
        {"manifest_pointer": size_of(3, "log"),
         "equals": len(log_text.encode("utf-8"))},
        {"manifest_pointer": size_of(3, "attachment.geom"), "equals": len(geom)},
        {"manifest_pointer": size_of(4, "counts"), "equals": 24},
        {"section": "rho", "member": "data", "binary_index": [1, 2, 3],
         "equals": 23.0 / 24.0, "tol": 1e-7},
        {"section": "citations", "member": "references", "decode": "utf8",
         "contains": "Ångström"},
        {"section": "run", "member": "files",
         "json_pointer": ["log", "filename"], "equals": "water.out"},
        {"section": "x_sizes", "member": "counts", "binary_index": [0],
         "equals": 8},
    ], "Every member spec carries the optional size, the uncompressed byte "
       "count of its ZIP entry (spec § 4.2): JSON members (structure, bonds, "
       "the volume grid, the run.record files index), binary members (the "
       "volume data, a vendor int64 array), the opaque citations.references "
       "and the run.record input, log and attachment. A validator checks "
       "each size against the ZIP central directory before reading the "
       "member. The root declares byte_order 'little' (spec § 2.2). The "
       "checks compare the sizes of the binary and opaque members, whose "
       "bytes are fixed; a JSON member's size depends on its producer's "
       "formatting.")


def _updated_in_place() -> None:
    """18) An archive updated in place (spec § 2.3): a pending job.spec
    archive to which a runner appended its results and a new manifest.

    The original archive is written first, manifest last. The update then
    appends with ``zipfile`` mode ``"a"`` — which writes the new entries where
    the old central directory began and a new central directory that lists
    every entry — the members the original did not hold, the rewritten
    ``run.record`` ``files`` index under a new name, and a new
    ``manifest.json`` last. The first ``manifest.json`` and the old
    ``files.json`` are superseded: two orphans, each a ``W_ORPHAN_MEMBER``.
    """
    created, modified = "2026-10-08T09:00:00Z", "2026-10-08T09:05:41Z"
    uuid = "0b5e7c1a-2d4f-4e6a-9b8c-7d1e2f3a4b5c"
    input_text = ("! RHF STO-3G\n* xyz 0 1\nO 0 0 0.1173\n"
                  "H 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692\n*\n")
    log_text = ("democode 1.0.0 — water.inp\n"
                "SCF iteration   8  E = -74.965901\n"
                "SCF CONVERGED\n"
                "FINAL SINGLE POINT ENERGY  -74.965901\n")
    charges = [-0.3302, 0.1651, 0.1651]

    def job(run_status, files, **record):
        w = QvfWriter(_PROGRAM, _VERSION, "water-job")
        w.add_structure(_WATER_ATOMS)
        w.add_job_spec(job_type="molecular", method="rhf", basis="sto-3g",
                       charge=0, multiplicity=1, tasks=["single_point"])
        w.add_run_record(program="democode", program_version="1.0.0",
                         input_text=input_text, files=files, sequence=0,
                         section_id="run", **record)
        w.set_provenance(run_status=run_status)
        return w

    # The pending archive: structure, job.spec and the prepared input.
    before = job("pending", {"input": {"filename": "water.inp"}})
    before.set_metadata(created_utc=created, uuid=uuid,
                        title="Water RHF/STO-3G: a job updated in place")
    # The settled archive: the same sections, the log, Mulliken charges.
    after = job("converged", {"input": {"filename": "water.inp"},
                              "log": {"filename": "water.out"}},
                log_text=log_text, exit_status=0,
                started_utc="2026-10-08T09:05:00Z",
                finished_utc="2026-10-08T09:05:41Z")
    after.add_atom_properties(mulliken_charge=charges, section_id="charges")
    after.annotate_section("charges", structure_ref="structure", run_ref="run")
    after.set_provenance(run_status="converged", method="rhf", basis="sto-3g",
                         scf_converged=True,
                         scf_energy={"value": -74.965901, "units": "hartree"})
    after.set_metadata(created_utc=created, modified_utc=modified, uuid=uuid,
                       title="Water RHF/STO-3G: a job updated in place")

    old = _entries_of(before.to_bytes())
    old_data = {name: data for name, data, _ in old}
    original = _plain_zip([e for e in old if e[0] != "manifest.json"]
                          + [e for e in old if e[0] == "manifest.json"])

    # The update: reuse every entry whose bytes are unchanged, write a changed
    # member under a new name (§ 2.3: never reuse an entry name), append the
    # new members, and the new manifest last.
    manifest = after.build_manifest()
    new_data = {name: data for name, data, _ in _entries_of(after.to_bytes())}
    appended = []
    for section in manifest["sections"]:
        for spec in section["members"].values():
            path, data = spec["path"], new_data[spec["path"]]
            if old_data.get(path) == data:
                continue                     # an earlier entry, kept
            if path in old_data:             # changed: a new name
                stem, dot, ext = path.rpartition(".")
                spec["path"] = path = f"{stem}.2{dot}{ext}"
            appended.append((path, data, zipfile.ZIP_DEFLATED))
    assert [p for p, _, _ in appended] == [
        "run/log.txt", "run/files.2.json", "atom_properties/mulliken_charge.bin"]
    manifest_bytes = json.dumps(manifest, ensure_ascii=False,
                                indent=2).encode("utf-8")

    buf = io.BytesIO(original)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")          # "Duplicate name: manifest.json"
        with zipfile.ZipFile(buf, "a") as zf:
            for name, data, method in appended:
                zf.writestr(_zip_info(name, method), data)
            zf.writestr(_zip_info("manifest.json", zipfile.ZIP_STORED),
                        manifest_bytes)
    data = buf.getvalue()
    names = [i.filename for i in zipfile.ZipFile(io.BytesIO(data)).infolist()]
    assert names.count("manifest.json") == 2 and names[-1] == "manifest.json"

    path = os.path.join(CORPUS, "updated_in_place.qvf")
    with open(path, "wb") as fh:
        fh.write(data)
    sections = {s["id"]: i for i, s in enumerate(manifest["sections"])}
    expected = {
        "archive": "updated_in_place.qvf",
        "description": (
            "An archive updated in place (spec § 2.3). A pending job "
            "(structure, job.spec, a run.record holding the prepared input, "
            "provenance.run_status 'pending') was written with its manifest "
            "last; a runner then appended the run's log, a rewritten "
            "run.record files index under a new name (run/files.2.json), "
            "Mulliken charges and a new manifest.json, and wrote a new "
            "central directory listing every entry. The archive therefore "
            "holds two manifest.json entries: the last is the manifest "
            "(run_status 'converged', metadata.modified_utc set, created_utc "
            "and uuid kept), the first is superseded. The superseded "
            "manifest and run/files.json are orphans, so a validator reports "
            "W_ORPHAN_MEMBER twice and no error; E_ZIP_DUPLICATE_ENTRY does "
            "not apply to manifest.json. The new manifest references the "
            "kept entries (structure, job.spec, the input) and the appended "
            "ones. A consumer that reads the first manifest.json sees "
            "run_status 'pending' and fails the checks."),
        "source_program": _PROGRAM,
        "kinds": sorted({s["kind"] for s in manifest["sections"]}),
        "checks": [
            {"manifest_pointer": ["provenance", "run_status"],
             "equals": "converged"},
            {"manifest_pointer": ["metadata", "modified_utc"],
             "equals": modified},
            {"manifest_pointer": ["metadata", "created_utc"], "equals": created},
            {"manifest_pointer": ["metadata", "uuid"], "equals": uuid},
            {"manifest_pointer": ["sections", sections["run"], "sequence"],
             "equals": 0},
            {"manifest_pointer": ["sections", sections["run"], "members",
                                  "files", "path"],
             "equals": "run/files.2.json"},
            {"section": "job_spec", "member": "spec",
             "json_pointer": ["method"], "equals": "rhf"},
            {"section": "run", "member": "input", "decode": "utf8",
             "equals": input_text},
            {"section": "run", "member": "log", "decode": "utf8",
             "contains": "SCF CONVERGED"},
            {"section": "run", "member": "files",
             "json_pointer": ["log", "filename"], "equals": "water.out"},
            {"section": "charges", "member": "mulliken_charge",
             "binary_index": [0], "equals": charges[0], "tol": 1e-12},
            {"section": "structure", "member": "structure",
             "json_pointer": ["atoms", 0, "symbol"], "equals": "O"},
        ],
    }
    with open(os.path.join(CORPUS, "updated_in_place.expected.json"), "w",
              encoding="utf-8") as fh:
        json.dump(expected, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"  updated_in_place.qvf  ({len(manifest['sections'])} sections, "
          f"{len(expected['checks'])} checks)")

def _md_periodic() -> None:
    """17) md_periodic: a toy NPT molecular-dynamics run of one water molecule
    in a breathing, sheared cell (spec § 4.3, § 5.5; P12).

    Five frames 0.5 fs apart. The atoms move at constant velocity (Å/fs), so
    coords[k] = coords[0] + velocities * times[k]; the cell of frame k has rows
    a = (9.00 (1 + 0.002 k), 0, 0), b = (0.30, 9.10 (1 + 0.002 k), 0) and
    c = (0.20, 0.40, 9.20 (1 + 0.002 k)) in Å, so it is lower triangular and
    reading the rows as columns gives different numbers; cell_volumes is the
    determinant a_x b_y c_z. coords are float32 (metadata.precision); every
    other member is float64. The binary energies restate metadata.energies.
    """
    n_frames, dt = 5, 0.5
    base = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692],
                     [0.0, -0.7572, -0.4692]]) + 4.5
    velocities = np.array([[0.0012, -0.0004, 0.0003],
                           [-0.0150, 0.0080, 0.0110],
                           [0.0130, -0.0090, -0.0120]])
    times = [dt * k for k in range(n_frames)]
    coords = np.stack([base + velocities * t for t in times])
    vel = np.stack([velocities] * n_frames)
    cells = [Lattice.from_rows_angstrom([[9.00 * (1 + 0.002 * k), 0.0, 0.0],
                                         [0.30, 9.10 * (1 + 0.002 * k), 0.0],
                                         [0.20, 0.40, 9.20 * (1 + 0.002 * k)]])
             for k in range(n_frames)]
    volumes = [cell.volume_angstrom3() for cell in cells]
    grad = np.stack([np.array([[0.0, 0.0, 0.0100 - 0.001 * k],
                               [0.0, 0.0080 + 0.0005 * k, -0.0050],
                               [0.0, -0.0080 - 0.0005 * k, -0.0050]])
                     for k in range(n_frames)])
    energies = [-76.41200, -76.41185, -76.41210, -76.41195, -76.41220]
    temperatures = [300.0, 304.2, 298.7, 301.5, 299.9]
    f32 = coords.astype(np.float32)

    w = QvfWriter(_PROGRAM, _VERSION, "h2o-npt-md")
    w.add_structure([{"symbol": z, "position": list(map(float, xyz)),
                      "atomic_number": n}
                     for z, xyz, n in zip(["O", "H", "H"], coords[0], [8, 1, 1])],
                    pbc=[True, True, True], lattice_vectors=cells[0])
    w.add_trajectory(
        {"symbols": ["O", "H", "H"], "atomic_numbers": [8, 1, 1],
         "n_frames": n_frames, "kind": "md", "ensemble": "npt",
         "timestep_fs": dt, "pbc": [True, True, True], "energies": energies},
        coords, lattice=cells, velocities=vel, gradient=grad, times=times,
        energies=energies, cell_volumes=volumes, temperatures=temperatures,
        precision="float32",
        units={"coords": "angstrom", "lattice": "angstrom", "times": "fs",
               "timestep_fs": "fs", "energies": "hartree",
               "cell_volumes": "a3", "temperatures": "kelvin"},
        section_id="md")
    _write("md_periodic", w, [
        {"manifest_pointer": ["sections", 1, "members", "coords", "dtype"],
         "equals": "float32"},
        {"section": "md", "member": "metadata", "json_pointer": ["precision"],
         "equals": "float32"},
        {"section": "md", "member": "metadata", "json_pointer": ["kind"],
         "equals": "md"},
        {"section": "md", "member": "metadata", "json_pointer": ["ensemble"],
         "equals": "npt"},
        {"section": "md", "member": "metadata", "json_pointer": ["timestep_fs"],
         "equals": dt, "tol": 1e-12},
        {"section": "md", "member": "metadata", "json_pointer": ["pbc"],
         "equals": [True, True, True]},
        {"section": "md", "member": "coords", "binary_index": [4, 1, 2],
         "equals": float(f32[4, 1, 2]), "tol": 1e-6},
        {"section": "md", "member": "coords", "binary_index": [2, 2, 0],
         "equals": float(f32[2, 2, 0]), "tol": 1e-6},
        {"manifest_pointer": ["sections", 1, "members", "lattice", "shape"],
         "equals": [n_frames, 3, 3]},
        # Rows in Å: row b of frame 3 starts with 0.30; read as columns the
        # same element would be a_y = 0.
        {"section": "md", "member": "lattice", "binary_index": [3, 1, 0],
         "equals": 0.30, "tol": 1e-12},
        {"section": "md", "member": "lattice", "binary_index": [3, 0, 1],
         "equals": 0.0, "tol": 1e-12},
        {"section": "md", "member": "lattice", "binary_index": [4, 2, 2],
         "equals": cells[4].rows_angstrom()[2][2], "tol": 1e-12},
        {"section": "md", "member": "velocities", "binary_index": [2, 1, 0],
         "equals": float(velocities[1, 0]), "tol": 1e-12},
        {"section": "md", "member": "gradient", "binary_index": [3, 0, 2],
         "equals": float(grad[3, 0, 2]), "tol": 1e-12},
        {"section": "md", "member": "times", "binary_index": [3],
         "equals": 1.5, "tol": 1e-12},
        {"section": "md", "member": "energies", "binary_index": [2],
         "equals": energies[2], "tol": 1e-12},
        {"section": "md", "member": "metadata", "json_pointer": ["energies", 2],
         "equals": energies[2], "tol": 1e-12},
        {"section": "md", "member": "cell_volumes", "binary_index": [4],
         "equals": volumes[4], "tol": 1e-9},
        {"section": "md", "member": "temperatures", "binary_index": [1],
         "equals": 304.2, "tol": 1e-12},
    ], "A toy NPT molecular-dynamics run of one water molecule: coords float32 "
       "[n_frames, n_atoms, 3] in Angstrom with metadata.precision 'float32'; "
       "a per-frame lattice [n_frames, 3, 3] in the structure convention "
       "(rows, Angstrom; not the columns in bohr of reaction.path) with pbc; "
       "velocities in Angstrom/fs, the energy gradient in hartree/bohr, times "
       "in fs, binary energies in hartree that restate metadata.energies, "
       "cell volumes in Angstrom^3 and temperatures in kelvin; ensemble npt "
       "and timestep_fs 0.5.")


# WP-23: wavefunction extensions (spec § 5.3, Appendix A.1): an effective core
# potential with the electron count it implies, and AO overlap matrices that
# make the normalization contract checkable.
# =============================================================================

# Real solid harmonics as Cartesian monomials, (weight, (l_x, l_y, l_z)), in the
# Appendix A.2 order m = -l..+l. On N_i-normalized primitives (Appendix A.1)
# each is normalized like the shell's axial component.
_SQRT3 = math.sqrt(3.0)
_SOLID_HARMONICS = {
    0: [[(1.0, (0, 0, 0))]],
    1: [[(1.0, (0, 1, 0))], [(1.0, (0, 0, 1))], [(1.0, (1, 0, 0))]],
    2: [[(_SQRT3, (1, 1, 0))],
        [(_SQRT3, (0, 1, 1))],
        [(1.0, (0, 0, 2)), (-0.5, (2, 0, 0)), (-0.5, (0, 2, 0))],
        [(_SQRT3, (1, 0, 1))],
        [(0.5 * _SQRT3, (2, 0, 0)), (-0.5 * _SQRT3, (0, 2, 0))]],
}


def _ao_terms(l, pure):
    """Each AO of a shell as a list of (weight, powers) monomials, in the
    Appendix A.2 order."""
    if pure:
        return _SOLID_HARMONICS[l]
    return [[(1.0, (lx, ly, l - lx - ly))]
            for lx in range(l, -1, -1) for ly in range(l - lx, -1, -1)]


def _overlap_1d(a, b, xa, xb, alpha, beta):
    """int (x - xa)^a (x - xb)^b exp(-alpha (x - xa)^2 - beta (x - xb)^2) dx,
    by expanding both monomials about the Gaussian product centre."""
    p = alpha + beta
    xp = (alpha * xa + beta * xb) / p
    total = 0.0
    for i in range(a + 1):
        for j in range(b + 1):
            if (i + j) % 2:
                continue
            total += (math.comb(a, i) * math.comb(b, j) * (xp - xa) ** (a - i)
                      * (xp - xb) ** (b - j)
                      * _odd_double_factorial((i + j) // 2)
                      / (2.0 * p) ** ((i + j) // 2) * math.sqrt(math.pi / p))
    return math.exp(-alpha * beta / p * (xa - xb) ** 2) * total


def _ao_overlap(shells, centres_bohr):
    """The AO overlap matrix of basis.json ``shells`` (coefficients on
    N_i-normalized primitives, an explicit ``pure`` on every shell), with the
    atoms at ``centres_bohr``: S[mu, nu] = <chi_mu | chi_nu>."""
    aos = [(centres_bohr[sh["center"]], sh["l"], sh["exponents"],
            sh["coefficients"], terms)
           for sh in shells for terms in _ao_terms(sh["l"], sh["pure"])]
    s = np.zeros((len(aos), len(aos)))
    for mu, (ra, la, ea, ca, ta) in enumerate(aos):
        for nu, (rb, lb, eb, cb, tb) in enumerate(aos):
            for a, c1 in zip(ea, ca):
                for b, c2 in zip(eb, cb):
                    prim = sum(w1 * w2 * math.prod(
                        _overlap_1d(p1[k], p2[k], ra[k], rb[k], a, b)
                        for k in range(3)) for w1, p1 in ta for w2, p2 in tb)
                    s[mu, nu] += (c1 * c2 * primitive_norm(a, la)
                                  * primitive_norm(b, lb) * prim)
    return s


def _lowdin(s):
    """S^(-1/2): its rows are MOs that are orthonormal in the metric S."""
    eigval, eigvec = np.linalg.eigh(s)
    return eigvec @ np.diag(eigval ** -0.5) @ eigvec.T


def _wavefunction_extensions() -> None:
    """17) An effective core potential; 18) an overlap member (spec § 5.3)."""
    ang = BOHR_PER_ANGSTROM

    # 17) NaH with a 10-electron ECP on Na: Z(Na) + Z(H) - 10 = 2 electrons.
    na_s, na_p, h_s = [0.85, 0.11], [0.25], [1.24, 0.27]
    shells = [
        {"center": 0, "l": 0, "exponents": na_s,
         "coefficients": _axial_normalized(na_s, [0.3, 0.8], 0), "pure": True},
        {"center": 0, "l": 1, "exponents": na_p,
         "coefficients": _axial_normalized(na_p, [1.0], 1), "pure": True},
        {"center": 1, "l": 0, "exponents": h_s,
         "coefficients": _axial_normalized(h_s, [0.4, 0.7], 0), "pure": True}]
    r_nah = 1.887
    overlap = _ao_overlap(shells, [(0.0, 0.0, 0.0), (0.0, 0.0, r_nah * ang)])
    C = _lowdin(overlap)
    density = 2.0 * np.outer(C[0], C[0])
    ecp = [{"center": 0, "n_core_electrons": 10, "label": "toy-ECP10",
            "projectors": [
                {"l": 2, "exponents": [1.05, 0.31], "coefficients": [-1.20, -0.25],
                 "powers": [2, 2], "local": True},
                {"l": 0, "exponents": [2.40, 0.74], "coefficients": [11.5, -2.1],
                 "powers": [0, 2]},
                {"l": 1, "exponents": [1.70], "coefficients": [5.6],
                 "powers": [1]}]}]
    w = QvfWriter(_PROGRAM, _VERSION, "nah-ecp")
    w.add_structure([{"symbol": "Na", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 11},
                     {"symbol": "H", "position": [0.0, 0.0, r_nah],
                      "atomic_number": 1}])
    w.add_wavefunction_gto(
        shells, mo_coefficients=C, energies=[-0.21, 0.03, 0.06, 0.03, 0.34],
        occupations=[2.0, 0.0, 0.0, 0.0, 0.0], ecp=ecp, n_electrons=2,
        symmetry_labels=["sigma", "pi", "sigma", "pi", "sigma"],
        overlap=overlap, density_matrix=density)
    w.set_provenance(method="rhf", basis="toy", charge=0, multiplicity=1,
                     n_electrons=2)
    _write("wavefunction_ecp", w, [
        {"section": "wf", "member": "basis", "json_pointer": ["n_electrons"],
         "equals": 2},
        {"manifest_pointer": ["provenance", "n_electrons"], "equals": 2},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "center"], "equals": 0},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "n_core_electrons"], "equals": 10},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "label"], "equals": "toy-ECP10"},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "projectors", 0, "local"], "equals": True},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "projectors", 0, "l"], "equals": 2},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "projectors", 1, "l"], "equals": 0},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "projectors", 1, "powers"], "equals": [0, 2]},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "projectors", 1, "coefficients", 0],
         "equals": 11.5, "tol": 1e-12},
        {"section": "wf", "member": "basis",
         "json_pointer": ["ecp", 0, "projectors", 2, "exponents", 0],
         "equals": 1.70, "tol": 1e-12},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["occupations", 0], "equals": 2.0},
        {"section": "wf", "member": "mo_metadata",
         "json_pointer": ["symmetry_labels", 1], "equals": "pi"},
        {"manifest_pointer": ["sections", 1, "members", "overlap", "shape"],
         "equals": [5, 5]},
        {"section": "wf", "member": "overlap", "binary_index": [0, 0],
         "equals": 1.0, "tol": 1e-12},
        {"section": "wf", "member": "overlap", "binary_index": [0, 4],
         "equals": float(overlap[0, 4]), "tol": 1e-12},
        {"section": "wf", "member": "overlap", "binary_index": [1, 2],
         "equals": 0.0, "tol": 1e-12},
        {"section": "wf", "member": "density_matrix", "binary_index": [0, 0],
         "equals": float(density[0, 0]), "tol": 1e-12},
    ], "wavefunction.gto with an effective core potential: NaH, the ECP on Na "
       "replacing 10 core electrons, so basis.n_electrons = 11 + 1 - 10 = 2, "
       "which the occupations and provenance.n_electrons restate. The ECP has "
       "a local part (local: true, l = 2, the highest l of the centre) and s "
       "and p semi-local channels, each U(r) = sum c r^(n - 2) exp(-a r^2) "
       "with the radial powers n in powers. The section also carries "
       "symmetry_labels, the AO overlap matrix (whose diagonal is 1: every "
       "shell is contraction-normalized) and the density matrix of the doubly "
       "occupied orbital. The MOs are the Lowdin orbitals S^(-1/2); every "
       "number is illustrative, not a converged calculation.")

    # 18) An AO overlap matrix whose diagonal is the Appendix A.1 numbers:
    # spherical s, p, d on C, a spherical s and a Cartesian d on O.
    c_s, c_p, c_d = [3.2, 0.62], [2.9, 0.55], [0.8]
    o_s, o_d = [5.8, 1.1], [2.3, 0.65]
    shells = [
        {"center": 0, "l": 0, "exponents": c_s,
         "coefficients": _axial_normalized(c_s, [0.35, 0.75], 0), "pure": True},
        {"center": 0, "l": 1, "exponents": c_p,
         "coefficients": _axial_normalized(c_p, [0.45, 0.70], 1), "pure": True},
        {"center": 0, "l": 2, "exponents": c_d,
         "coefficients": _axial_normalized(c_d, [1.0], 2), "pure": True},
        # Deliberately not contraction-normalized: its diagonal is S_axial.
        {"center": 1, "l": 0, "exponents": o_s, "coefficients": [0.4, 0.5],
         "pure": True},
        {"center": 1, "l": 2, "exponents": o_d,
         "coefficients": _axial_normalized(o_d, [0.55, 0.62], 2), "pure": False}]
    r_co = 1.128
    overlap = _ao_overlap(shells, [(0.0, 0.0, 0.0), (0.0, 0.0, r_co * ang)])
    C = _lowdin(overlap)
    density = 2.0 * C[:7].T @ C[:7]
    ao_labels = ["C s", "C p(m=-1)", "C p(m=0)", "C p(m=+1)",
                 "C d(m=-2)", "C d(m=-1)", "C d(m=0)", "C d(m=+1)", "C d(m=+2)",
                 "O s", "O d_xx", "O d_xy", "O d_xz", "O d_yy", "O d_yz",
                 "O d_zz"]
    w = QvfWriter(_PROGRAM, _VERSION, "co-overlap")
    w.add_structure([{"symbol": "C", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 6},
                     {"symbol": "O", "position": [0.0, 0.0, r_co],
                      "atomic_number": 8}])
    w.add_wavefunction_gto(
        shells, mo_coefficients=C, energies=[0.0] * 16,
        occupations=[2.0] * 7 + [0.0] * 9, n_electrons=14, overlap=overlap,
        density_matrix=density)
    checks = [
        {"section": "wf", "member": "basis", "json_pointer": ["n_ao"],
         "equals": 16},
        {"section": "wf", "member": "basis", "json_pointer": ["pure"],
         "equals": True},
        {"section": "wf", "member": "basis",
         "json_pointer": ["shells", 4, "pure"], "equals": False},
        {"manifest_pointer": ["sections", 1, "members", "overlap", "shape"],
         "equals": [16, 16]},
        {"manifest_pointer": ["sections", 1, "members", "density_matrix",
                              "shape"], "equals": [16, 16]},
    ]
    # The diagonal, AO by AO: 1 for every spherical AO of a normalized shell
    # and for the axial Cartesian ones, 1/3 for d_xy, d_xz and d_yz, and
    # S_axial of the unnormalized O s shell.
    checks += [{"section": "wf", "member": "overlap", "binary_index": [k, k],
                "equals": float(overlap[k, k]), "tol": 1e-12}
               for k in range(16)]
    checks += [
        {"section": "wf", "member": "overlap", "binary_index": [0, 9],
         "equals": float(overlap[0, 9]), "tol": 1e-12},
        {"section": "wf", "member": "overlap", "binary_index": [9, 0],
         "equals": float(overlap[0, 9]), "tol": 1e-12},
        {"section": "wf", "member": "mo_coefficients", "binary_index": [0, 0],
         "equals": float(C[0, 0]), "tol": 1e-12},
        {"section": "wf", "member": "density_matrix", "binary_index": [9, 9],
         "equals": float(density[9, 9]), "tol": 1e-12},
    ]
    _write("wavefunction_overlap", w, checks,
           "wavefunction.gto with an AO overlap member, so that the Appendix "
           "A.1 normalization contract is checked by every validator: CO with "
           "spherical s, p and d shells on C and a spherical s and a "
           "Cartesian (pure: false) d shell on O. The overlap diagonal is 1 "
           "for every spherical AO and for d_xx, d_yy and d_zz, 1/3 for d_xy, "
           "d_xz and d_yz, and for the O s shell, whose contraction is "
           "deliberately not normalized, its axial self-overlap "
           f"{overlap[9, 9]:.12f}. AO order: " + ", ".join(ao_labels) + ". "
           "The MOs are the Lowdin orbitals S^(-1/2), orthonormal in this "
           "metric, and the density matrix is that of the seven lowest, doubly "
           "occupied (tr(P S) = 14 = basis.n_electrons); illustrative "
           "numbers, not a converged calculation.")


# WP-25: the vector, complex, series and Laplacian volume kinds, and a
# periodic grid in both endpoint conventions (spec § 5.2).
# =============================================================================

def _wp25_volume_archives() -> None:
    """17)–21) ``volume_periodic``, ``volume_vector``, ``volume_complex``,
    ``volume_series`` and ``volume_laplacian``."""
    ang = BOHR_PER_ANGSTROM

    # volume_periodic: one periodic density on a hexagonal (non-orthogonal)
    # cell, sampled CHGCAR style (includes_endpoint false: shape[a] planes per
    # period) and cube style (includes_endpoint true: the same planes plus the
    # image of plane 0). Both grids have the same voxel vectors, cell row / n.
    rows = [[2.46, 0.0, 0.0], [-1.23, round(2.46 * math.sqrt(3) / 2, 6), 0.0],
            [0.0, 0.0, 3.35]]
    cell = np.array(rows) * ang                              # rows, bohr
    second = [round(x, 6) for x in (np.array(rows[0]) + 2 * np.array(rows[1])) / 3]
    atoms = [{"symbol": "C", "position": [0.0, 0.0, 0.0], "atomic_number": 6},
             {"symbol": "C", "position": second, "atomic_number": 6}]
    sites = np.array([a["position"] for a in atoms]) * ang
    n = (6, 6, 8)
    frac = np.stack(np.meshgrid(*[np.arange(k) / k for k in n], indexing="ij"),
                    axis=-1)
    r = frac @ cell
    rho = np.zeros(n)
    for site in sites:
        for image in itertools.product((-1, 0, 1), repeat=3):
            d = r - site - np.array(image, dtype=float) @ cell
            rho += 0.4 * np.exp(-0.6 * np.sum(d * d, axis=-1))
    rho = rho.astype(np.float32)
    cube = np.pad(rho, ((0, 1), (0, 1), (0, 1)), mode="wrap")  # + plane 0 again
    vox = [(cell[a] / n[a]).tolist() for a in range(3)]
    w = QvfWriter(_PROGRAM, _VERSION, "hexagonal-cell-density")
    w.add_structure(atoms, pbc=[True, True, True], lattice_vectors=rows)
    w.add_density(rho, [0.0, 0.0, 0.0], vox, label="CHGCAR-style grid",
                  periodic=[True, True, True], includes_endpoint=False,
                  structure_ref="structure", section_id="rho_chgcar")
    w.add_density(cube, [0.0, 0.0, 0.0], vox, label="cube-style grid",
                  periodic=[True, True, True], includes_endpoint=True,
                  section_id="rho_cube")
    _write("volume_periodic", w, [
        {"section": "rho_chgcar", "member": "grid", "json_pointer": ["periodic"],
         "equals": [True, True, True]},
        {"section": "rho_chgcar", "member": "grid",
         "json_pointer": ["includes_endpoint"], "equals": False},
        {"section": "rho_chgcar", "member": "grid", "json_pointer": ["shape"],
         "equals": list(n)},
        {"section": "rho_chgcar", "member": "grid",
         "json_pointer": ["structure_ref"], "equals": "structure"},
        {"section": "rho_chgcar", "member": "grid",
         "json_pointer": ["voxel_vectors", 1, 0], "equals": vox[1][0],
         "tol": 1e-12},
        {"section": "rho_chgcar", "member": "grid",
         "json_pointer": ["voxel_vectors", 2, 2], "equals": vox[2][2],
         "tol": 1e-12},
        {"section": "rho_cube", "member": "grid",
         "json_pointer": ["includes_endpoint"], "equals": True},
        {"section": "rho_cube", "member": "grid", "json_pointer": ["shape"],
         "equals": [k + 1 for k in n]},
        {"section": "rho_chgcar", "member": "data", "binary_index": [1, 2, 3],
         "equals": float(rho[1, 2, 3]), "tol": 1e-6},
        {"section": "rho_chgcar", "member": "data", "binary_index": [5, 0, 7],
         "equals": float(rho[5, 0, 7]), "tol": 1e-6},
        # The cube grid holds the same points, plus the image of plane 0 as
        # its last plane along every axis.
        {"section": "rho_cube", "member": "data", "binary_index": [1, 2, 3],
         "equals": float(rho[1, 2, 3]), "tol": 1e-6},
        {"section": "rho_cube", "member": "data", "binary_index": [6, 2, 3],
         "equals": float(rho[0, 2, 3]), "tol": 1e-6},
        {"section": "rho_cube", "member": "data", "binary_index": [3, 6, 8],
         "equals": float(rho[3, 0, 0]), "tol": 1e-6},
    ], "Periodic grids (spec § 5.2): a density on a hexagonal cell (pbc all "
       "true, non-orthogonal a and b) sampled CHGCAR style, includes_endpoint "
       "false, shape 6x6x8, so voxel_vectors[a] x shape[a] is lattice row a "
       "in bohr; and cube style, includes_endpoint true, shape 7x7x9 with the "
       "same voxel vectors, so voxel_vectors[a] x (shape[a] - 1) is the row "
       "and the last plane on every axis repeats plane 0. A consumer that "
       "tiles the cube-style grid drops that plane.")

    # volume_vector: the gradient of a model H2 density, [nx, ny, nz, 3].
    shape, step = (4, 5, 6), 0.8
    origin, voxels, xyz = _volume_grid(shape, step)
    grad = np.zeros(shape + (3,))
    for z0 in (-0.7, 0.7):
        g = np.exp(-0.8 * _r2(xyz, (0.0, 0.0, z0)))
        for c, x0 in enumerate((0.0, 0.0, z0)):
            grad[..., c] += -1.6 * (xyz[c] - x0) * g
    grad = grad.astype(np.float32)
    w = QvfWriter(_PROGRAM, _VERSION, "h2-density-gradient")
    w.add_structure(
        [{"symbol": "H", "position": [0.0, 0.0, z / ang], "atomic_number": 1}
         for z in (-0.7, 0.7)])
    w.add_volume("volume.vector", grad, origin, voxels,
                 label="Density gradient", section_id="grad_rho")
    _write("volume_vector", w, [
        {"manifest_pointer": ["sections", 1, "members", "data", "shape"],
         "equals": [4, 5, 6, 3]},
    ] + _volume_checks("grad_rho", grad, shape, origin, step,
                       [(1, 1, 2, 0), (2, 3, 4, 1), (0, 4, 1, 2), (3, 0, 5, 2)]),
        "volume.vector: the gradient of a model H2 electron density, float32 "
        "[nx, ny, nz, 3] on an anisotropic 4x5x6 grid in bohr; the last axis "
        "holds the x, y and z components, and the grid shape is the leading "
        "three extents of data.")

    # volume_complex: the hydrogen 2p orbital with m = +1, (x + i y) e^(-r/2)
    # / (8 sqrt(pi)) (no Condon-Shortley phase), float64 [nx, ny, nz, 2].
    shape, step = (6, 5, 4), 1.5
    origin, voxels, xyz = _volume_grid(shape, step)
    radial = np.exp(-0.5 * np.sqrt(_r2(xyz, (0.0, 0.0, 0.0)))) / (8 * math.sqrt(math.pi))
    psi = np.stack([xyz[0] * radial, xyz[1] * radial], axis=-1)
    w = QvfWriter(_PROGRAM, _VERSION, "h-2p-plus")
    w.add_structure([{"symbol": "H", "position": [0.0, 0.0, 0.0],
                      "atomic_number": 1}])
    w.add_volume("volume.complex", psi[..., 0] + 1j * psi[..., 1], origin,
                 voxels, dtype="float64", label="2p (m = +1)",
                 section_id="psi_2p1")
    _write("volume_complex", w, [
        {"manifest_pointer": ["sections", 1, "members", "data", "shape"],
         "equals": [6, 5, 4, 2]},
        {"manifest_pointer": ["sections", 1, "members", "data", "dtype"],
         "equals": "float64"},
    ] + _volume_checks("psi_2p1", psi, shape, origin, step,
                       [(1, 3, 2, 0), (1, 3, 2, 1), (4, 0, 1, 0), (4, 0, 1, 1)],
                       tol=1e-12),
        "volume.complex: the hydrogen 2p orbital with m = +1, (x + i y) "
        "exp(-r/2) / (8 sqrt(pi)), float64 [nx, ny, nz, 2] with the last axis "
        "[re, im], on a 6x5x4 grid in bohr.")

    # volume_series: the density of an H2 stretch at frames 0, 2 and 4 of a
    # five-frame trajectory, [n, nx, ny, nz] with a frame axis.
    bonds = [0.70, 0.80, 0.90, 1.00, 1.10]                     # Å
    coords = np.array([[[0.0, 0.0, -r / 2], [0.0, 0.0, r / 2]] for r in bonds])
    frames = [0, 2, 4]
    shape, step = (4, 4, 6), 0.9
    origin, voxels, xyz = _volume_grid(shape, step)
    slabs = np.zeros((len(frames),) + shape)
    for s, f in enumerate(frames):
        for z in coords[f, :, 2] * ang:
            slabs[s] += 0.3 * np.exp(-0.9 * _r2(xyz, (0.0, 0.0, z)))
    slabs = slabs.astype(np.float32)
    labels = [f"r = {bonds[f]:.2f} Å" for f in frames]
    w = QvfWriter(_PROGRAM, _VERSION, "h2-stretch-densities")
    w.add_structure([{"symbol": "H", "position": list(p), "atomic_number": 1}
                     for p in coords[0]])
    w.add_trajectory({"symbols": ["H", "H"], "n_frames": len(bonds)}, coords)
    w.add_volume("volume.series", slabs, origin, voxels,
                 label="Density along the stretch", section_id="rho_frames",
                 axis={"kind": "frame", "values": frames, "labels": labels,
                       "trajectory_ref": "trajectory"})
    _write("volume_series", w, [
        {"manifest_pointer": ["sections", 2, "members", "data", "shape"],
         "equals": [3, 4, 4, 6]},
        {"section": "rho_frames", "member": "axis", "json_pointer": ["kind"],
         "equals": "frame"},
        {"section": "rho_frames", "member": "axis", "json_pointer": ["values"],
         "equals": frames},
        {"section": "rho_frames", "member": "axis",
         "json_pointer": ["labels", 2], "equals": labels[2]},
        {"section": "rho_frames", "member": "axis",
         "json_pointer": ["trajectory_ref"], "equals": "trajectory"},
        {"section": "trajectory", "member": "coords", "binary_index": [4, 1, 2],
         "equals": 0.55, "tol": 1e-12},
    ] + _volume_checks("rho_frames", slabs, shape, origin, step,
                       [(0, 1, 2, 2), (1, 1, 2, 2), (2, 1, 2, 2), (2, 3, 0, 5)]),
        "volume.series: the density of a stretching H2 molecule at frames 0, "
        "2 and 4 of a five-frame trajectory, float32 [n, nx, ny, nz]; the "
        "axis member names the frames (kind frame, values 0-based frame "
        "indices, not slab indices), labels each slab, and names the "
        "trajectory through trajectory_ref. The grid shape is the trailing "
        "three extents of data.")

    # volume_laplacian: del^2 rho of a two-centre Gaussian density, e/bohr^5.
    shape, step = (5, 4, 6), 0.7
    origin, voxels, xyz = _volume_grid(shape, step)
    lap = np.zeros(shape)
    for c0, alpha, z0 in ((0.8, 1.2, -1.07), (1.1, 1.5, 1.07)):
        r2 = _r2(xyz, (0.0, 0.0, z0))
        lap += c0 * (4 * alpha ** 2 * r2 - 6 * alpha) * np.exp(-alpha * r2)
    lap = lap.astype(np.float32)
    assert lap.min() < 0.0 < lap.max()
    w = QvfWriter(_PROGRAM, _VERSION, "n2-laplacian")
    w.add_structure(
        [{"symbol": "N", "position": [0.0, 0.0, z / ang], "atomic_number": 7}
         for z in (-1.07, 1.07)])
    w.add_volume("volume.laplacian", lap, origin, voxels,
                 label="Laplacian of the electron density", section_id="lap")
    _write("volume_laplacian", w, _volume_checks(
        "lap", lap, shape, origin, step, [(2, 1, 1), (2, 2, 4), (0, 0, 0)]),
        "volume.laplacian: the Laplacian of a model two-centre density in "
        "e/bohr^5 (float32, 5x4x6 grid in bohr), negative near the nuclei "
        "and positive farther out; the grid declares no units.data, "
        "because e/bohr^5 has no unit name.")


# ==== WP-29: property kinds (begin) ===========================================
# One elementary charge times one ångström, in debye (CODATA 2018:
# e a0 = 2.541746473 D, a0 = 0.529177210903 Å).
_E_ANGSTROM_DEBYE = 4.80320471


def _wp29_properties() -> None:
    """17) properties.electric, properties.elastic, properties.dielectric
    (spec § 5.13).

    ``properties_molecular``: water with a three-point-charge model, so that
    every number is derived from the geometry. The root ``dipole_moment`` and
    the ``dipole`` of ``properties.electric`` state the same dipole, as § 5.13
    requires; the Buckingham quadrupole is taken about the centre of mass.
    ``properties_periodic``: a cubic CsCl-type cell with the cubic stiffness
    matrix (C11, C12, C44), its compliance and Voigt moduli, and the dielectric
    tensors (static = electronic + ionic) with one Born tensor per atom.
    """
    # --- molecular -----------------------------------------------------------
    atoms = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
             {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
             {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]
    xyz = np.array([a["position"] for a in atoms])
    masses = np.array([15.999, 1.008, 1.008])
    q = np.array([-0.8, 0.4, 0.4])                    # point charges, in e
    r = xyz - (masses[:, None] * xyz).sum(axis=0) / masses.sum()   # about the COM
    mu = _E_ANGSTROM_DEBYE * (q[:, None] * r).sum(axis=0)           # debye
    r2 = (r * r).sum(axis=1)
    quad = 0.5 * _E_ANGSTROM_DEBYE * (
        3.0 * np.einsum("i,ia,ib->ab", q, r, r)
        - np.eye(3) * float((q * r2).sum()))                         # debye Å
    alpha = np.diag([8.64, 9.15, 9.45])               # static, atomic units
    total = float(np.linalg.norm(mu))

    w = QvfWriter(_PROGRAM, _VERSION, "water-electric")
    w.add_structure(atoms)
    w.set_dipole_moment(total_debye=total, vector_debye=mu.tolist(),
                        origin="center_of_mass")
    w.add_electric_properties({
        "origin": "center_of_mass",
        "dipole": {"vector": mu, "total": total, "units": "debye"},
        "quadrupole": {"matrix": quad, "traceless": True,
                       "convention": "buckingham"},
        "polarizability": {"matrix": alpha, "frequency_ev": 0.0,
                           "isotropic": float(np.trace(alpha) / 3.0)},
        "units": {"dipole": "debye", "origin": "bohr"}})
    _write("properties_molecular", w, [
        {"manifest_pointer": ["dipole_moment", "total_debye"], "equals": total,
         "tol": 1e-12},
        {"manifest_pointer": ["dipole_moment", "vector_debye", 2],
         "equals": float(mu[2]), "tol": 1e-12},
        {"manifest_pointer": ["sections", 1, "kind"],
         "equals": "properties.electric"},
        {"section": "electric", "member": "properties", "json_pointer": ["origin"],
         "equals": "center_of_mass"},
        {"section": "electric", "member": "properties",
         "json_pointer": ["dipole", "vector", 2], "equals": float(mu[2]),
         "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["dipole", "total"], "equals": total, "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["dipole", "units"], "equals": "debye"},
        {"section": "electric", "member": "properties",
         "json_pointer": ["quadrupole", "convention"], "equals": "buckingham"},
        {"section": "electric", "member": "properties",
         "json_pointer": ["quadrupole", "traceless"], "equals": True},
        {"section": "electric", "member": "properties",
         "json_pointer": ["quadrupole", "matrix", 0, 0],
         "equals": float(quad[0, 0]), "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["quadrupole", "matrix", 2, 2],
         "equals": float(quad[2, 2]), "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["quadrupole", "matrix", 1, 2],
         "equals": float(quad[1, 2]), "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["quadrupole", "matrix", 2, 1],
         "equals": float(quad[1, 2]), "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["polarizability", "matrix", 1, 1], "equals": 9.15,
         "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["polarizability", "frequency_ev"], "equals": 0.0,
         "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["polarizability", "isotropic"],
         "equals": float(np.trace(alpha) / 3.0), "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["units", "dipole"], "equals": "debye"},
    ], "Water with the electric properties of a three-point-charge model "
       "(O -0.8 e, H +0.4 e): the root dipole_moment and the dipole of the "
       "properties.electric section state one dipole (spec § 5.13 requires "
       "them to agree within 1e-6 debye); the Buckingham quadrupole "
       "(debye·Å, traceless, symmetric) is taken about the centre of mass, and "
       "the static polarizability (atomic units) is diagonal.")

    # --- periodic ------------------------------------------------------------
    a = 4.12
    c11, c12, c44 = 36.4, 8.8, 8.0                    # GPa
    stiffness = np.zeros((6, 6))
    for i in range(3):
        stiffness[i, i] = c11
        stiffness[3 + i, 3 + i] = c44
        for j in range(3):
            if i != j:
                stiffness[i, j] = c12
    compliance = np.linalg.inv(stiffness)
    bulk = (c11 + 2.0 * c12) / 3.0                    # Voigt averages, cubic
    shear = (c11 - c12 + 3.0 * c44) / 5.0
    young = 9.0 * bulk * shear / (3.0 * bulk + shear)
    poisson = (3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear))
    electronic = np.eye(3) * 2.62
    ionic = np.eye(3) * 4.58
    static = electronic + ionic
    born = np.stack([np.eye(3) * 1.1, np.eye(3) * -1.1])      # sum rule: 0

    w = QvfWriter(_PROGRAM, _VERSION, "cscl-elastic-dielectric")
    w.add_structure(
        [{"symbol": "Cs", "position": [0.0, 0.0, 0.0], "atomic_number": 55},
         {"symbol": "Cl", "position": [a / 2, a / 2, a / 2], "atomic_number": 17}],
        pbc=[True, True, True],
        lattice_vectors=[[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    w.add_elastic_properties(
        stiffness, compliance=compliance,
        moduli={"bulk_modulus_gpa": bulk, "shear_modulus_gpa": shear,
                "youngs_modulus_gpa": young, "poisson_ratio": poisson},
        average="voigt",
        units={"stiffness": "gpa", "bulk_modulus_gpa": "gpa",
               "poisson_ratio": "dimensionless"})
    w.add_dielectric_properties(
        electronic, ionic=ionic, static=static, born_charges=born,
        piezoelectric=np.zeros((3, 6)),
        units={"electronic": "dimensionless", "born_charges": "e"})
    _write("properties_periodic", w, [
        {"manifest_pointer": ["sections", 1, "kind"],
         "equals": "properties.elastic"},
        {"manifest_pointer": ["sections", 2, "kind"],
         "equals": "properties.dielectric"},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["stiffness", 0, 0], "equals": c11, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["stiffness", 0, 1], "equals": c12, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["stiffness", 1, 0], "equals": c12, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["stiffness", 5, 5], "equals": c44, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["stiffness", 0, 3], "equals": 0.0, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["compliance", 0, 0], "equals": float(compliance[0, 0]),
         "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["bulk_modulus_gpa"], "equals": bulk, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["shear_modulus_gpa"], "equals": shear, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["youngs_modulus_gpa"], "equals": young, "tol": 1e-9},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["poisson_ratio"], "equals": poisson, "tol": 1e-12},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["average"], "equals": "voigt"},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["units", "stiffness"], "equals": "gpa"},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["electronic", 0, 0], "equals": 2.62, "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["ionic", 1, 1], "equals": 4.58, "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["static", 2, 2], "equals": float(static[2, 2]),
         "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["static", 0, 1], "equals": 0.0, "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["born_charges", 0, 0, 0], "equals": 1.1, "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["born_charges", 1, 2, 2], "equals": -1.1, "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["born_charges", 1, 0, 1], "equals": 0.0, "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["piezoelectric", 1, 4], "equals": 0.0, "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["units", "born_charges"], "equals": "e"},
    ], "A cubic CsCl-type cell (a = 4.12 Å, two atoms). properties.elastic "
       "carries the cubic stiffness matrix in Voigt notation (C11 = 36.4, "
       "C12 = 8.8, C44 = 8.0 GPa), its compliance and the Voigt bulk, shear "
       "and Young moduli and Poisson ratio. properties.dielectric carries "
       "the electronic and ionic tensors and the static one, which equals "
       "their sum, two Born tensors (+1.1 e on Cs and -1.1 e on Cl, one per "
       "atom) and a zero piezoelectric tensor (the cell is centrosymmetric).")
# ==== WP-29: property kinds (end) =============================================


def build_zero_spin():
    """A valid high-ratio reference-writer archive (2 MiB uncompressed)."""
    w = QvfWriter(_PROGRAM, _VERSION, "zero-spin")
    w.add_structure([{"symbol": "He", "position": [0, 0, 0], "atomic_number": 2}])
    w.add_volume("volume.spin", np.zeros((64, 64, 64), dtype=np.float64),
                 [0, 0, 0], np.eye(3), dtype="float64", section_id="spin")
    _write("volume_spin_zero", w, [
        {"section": "spin", "member": "grid", "json_pointer": ["shape"],
         "equals": [64, 64, 64]},
        {"section": "spin", "member": "data", "binary_index": [63, 63, 63],
         "equals": 0.0, "tol": 0.0},
    ], "64 cubed float64 zero spin density: valid despite DEFLATE ratio >1000.")


def _bloch_model():
    """The model behind ``bloch.qvf`` (spec § 5.3.1), shared with the tests.

    A dimerised chain of hydrogen atoms along x (a = 2.0 Å, H-H 0.8 Å inside
    the cell), one normalised s Gaussian per atom (α = 0.4166 bohr⁻², the
    STO-1G H 1s exponent) and an extended-Hückel Hamiltonian with
    spin-dependent on-site energies, so the two spin channels differ in
    energies and coefficients. At k (fractional) the overlap and the
    Hamiltonian are the Bloch sums S_μν(k) = Σ_n e^{2πi k n} ⟨χ_μ|χ_ν(r − n a)⟩
    and H_μν(k) likewise, with H_μν(T) = K (h_μ + h_ν)/2 S_μν(T) off the
    on-site diagonal (K = 1.75). The orbitals solve H C = S C ε and are
    normalised as C†S(k)C = 1; each column's largest component is made real
    and positive, and the X-point orbitals are then multiplied by the gauge
    phase e^{iπ/3}, which changes no observable but makes their coefficients
    complex. One electron per spin and cell fills the lower band.
    """
    a_ang, d_ang, alpha, k_wh = 2.0, 0.8, 0.4166, 1.75
    x = np.array([0.0, d_ang]) * BOHR_PER_ANGSTROM
    a = a_ang * BOHR_PER_ANGSTROM
    onsite = [(-0.53, -0.51), (-0.49, -0.47)]      # [spin][atom], hartree
    kpoints = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]
    gauge = [1.0, complex(math.cos(math.pi / 3), math.sin(math.pi / 3))]
    cells = range(-12, 13)

    def overlap(mu, nu, n):
        r = x[nu] + n * a - x[mu]
        return math.exp(-alpha * r * r / 2.0)     # equal-exponent s Gaussians

    coeff = np.zeros((2, 2, 2, 2), dtype=complex)  # [k, spin, mo, ao]
    energies = np.zeros((2, 2, 2))
    overlaps = []
    for ik, kp in enumerate(kpoints):
        phase = {n: complex(math.cos(2 * math.pi * kp[0] * n),
                            math.sin(2 * math.pi * kp[0] * n)) for n in cells}
        s_k = np.array([[sum(phase[n] * overlap(m, v, n) for n in cells)
                         for v in range(2)] for m in range(2)])
        overlaps.append(s_k)
        sval, svec = np.linalg.eigh(s_k)
        x_half = svec @ np.diag(sval ** -0.5) @ svec.conj().T
        for spin in range(2):
            h = onsite[spin]
            h_k = np.array([[sum(phase[n] * (h[m] if (m == v and n == 0) else
                                             k_wh * 0.5 * (h[m] + h[v])
                                             * overlap(m, v, n))
                                 for n in cells) for v in range(2)]
                            for m in range(2)])
            eps, u = np.linalg.eigh(x_half.conj().T @ h_k @ x_half)
            c = x_half @ u                         # columns are orbitals
            for m in range(2):
                col = c[:, m]
                big = col[np.argmax(np.abs(col))]
                c[:, m] = col * (abs(big) / big) * gauge[ik]
            coeff[ik, spin] = c.T
            energies[ik, spin] = eps
    occupations = np.zeros((2, 2, 2))
    occupations[:, :, 0] = 1.0
    return {"a_ang": a_ang, "d_ang": d_ang, "alpha": alpha, "kpoints": kpoints,
            "coeff": coeff, "energies": energies, "occupations": occupations,
            "overlaps": overlaps}


def _bloch() -> None:
    """17) wavefunction.bloch + bands (spec § 5.3.1, § 5.6).

    The dimerised hydrogen chain of ``_bloch_model`` at Γ and X, with k-point
    weights 1/2 and 1/2, two spin channels, real coefficients at Γ and complex
    ones at X, and a bands section whose eigenvalues at the same two k-points
    are the orbital energies in eV, as spec § 5.3.1 recommends.
    """
    m = _bloch_model()
    energies, occ = m["energies"], m["occupations"]
    homo = float(energies[:, :, 0].max())
    lumo = float(energies[:, :, 1].min())
    fermi_ev = 0.5 * (homo + lumo) * _HARTREE_EV
    w = QvfWriter(_PROGRAM, _VERSION, "h2-chain-bloch")
    w.add_structure([{"symbol": "H", "position": [0.0, 0.0, 0.0], "atomic_number": 1},
                     {"symbol": "H", "position": [m["d_ang"], 0.0, 0.0],
                      "atomic_number": 1}],
                    pbc=[True, False, False],
                    lattice_vectors=[[m["a_ang"], 0.0, 0.0], [0.0, 10.0, 0.0],
                                     [0.0, 0.0, 10.0]])
    shells = [{"center": c, "l": 0, "exponents": [m["alpha"]], "coefficients": [1.0]}
              for c in (0, 1)]
    w.add_wavefunction_bloch(shells, m["kpoints"], m["coeff"], energies, occ,
                             weights=[0.5, 0.5], spin="unrestricted",
                             labels=["G", "X"], mesh=[2, 1, 1],
                             time_reversal=False, fermi_energy_ev=fermi_ev)
    eig_ev = np.transpose(energies, (1, 0, 2)) * _HARTREE_EV   # [spin, k, band]
    w.add_bands({"kpoints": m["kpoints"], "labels": ["G", "X"],
                 "segments": [{"label_start": "G", "label_end": "X",
                               "k_start": [0.0, 0.0, 0.0], "k_end": [0.5, 0.0, 0.0],
                               "n_points": 2}],
                 "fermi": fermi_ev, "n_spin": 2, "n_kpoints": 2, "n_bands": 2},
                eig_ev, occupations=np.transpose(occ, (1, 0, 2)))
    c = m["coeff"]
    _write("bloch", w, [
        {"manifest_pointer": ["sections", 1, "kind"], "equals": "wavefunction.bloch"},
        {"manifest_pointer": ["sections", 1, "members", "mo_coefficients", "shape"],
         "equals": [2, 2, 2, 2, 2]},
        {"manifest_pointer": ["sections", 1, "members", "energies", "shape"],
         "equals": [2, 2, 2]},
        {"section": "wf_bloch", "member": "kpoints", "json_pointer": ["kpoints", 1, 0],
         "equals": 0.5, "tol": 1e-12},
        {"section": "wf_bloch", "member": "kpoints", "json_pointer": ["weights", 1],
         "equals": 0.5, "tol": 1e-12},
        {"section": "wf_bloch", "member": "kpoints", "json_pointer": ["labels", 1],
         "equals": "X"},
        {"section": "wf_bloch", "member": "mo_metadata", "json_pointer": ["spin"],
         "equals": "unrestricted"},
        {"section": "wf_bloch", "member": "mo_metadata", "json_pointer": ["n_spin"],
         "equals": 2},
        {"section": "wf_bloch", "member": "mo_metadata", "json_pointer": ["complex"],
         "equals": True},
        {"section": "wf_bloch", "member": "mo_metadata",
         "json_pointer": ["k_point_convention"], "equals": "fractional"},
        {"section": "wf_bloch", "member": "basis", "json_pointer": ["n_ao"],
         "equals": 2},
        {"section": "wf_bloch", "member": "basis",
         "json_pointer": ["shells", 1, "exponents", 0], "equals": 0.4166, "tol": 1e-12},
        {"section": "wf_bloch", "member": "mo_coefficients",
         "binary_index": [0, 1, 0, 1, 0], "equals": float(c[0, 1, 0, 1].real),
         "tol": 1e-12},
        {"section": "wf_bloch", "member": "mo_coefficients",
         "binary_index": [0, 0, 1, 0, 1], "equals": 0.0, "tol": 1e-12},
        {"section": "wf_bloch", "member": "mo_coefficients",
         "binary_index": [1, 0, 0, 0, 0], "equals": float(c[1, 0, 0, 0].real),
         "tol": 1e-12},
        {"section": "wf_bloch", "member": "mo_coefficients",
         "binary_index": [1, 0, 0, 0, 1], "equals": float(c[1, 0, 0, 0].imag),
         "tol": 1e-12},
        {"section": "wf_bloch", "member": "mo_coefficients",
         "binary_index": [1, 1, 1, 1, 1], "equals": float(c[1, 1, 1, 1].imag),
         "tol": 1e-12},
        {"section": "wf_bloch", "member": "energies", "binary_index": [0, 0, 0],
         "equals": float(energies[0, 0, 0]), "tol": 1e-12},
        {"section": "wf_bloch", "member": "energies", "binary_index": [1, 1, 1],
         "equals": float(energies[1, 1, 1]), "tol": 1e-12},
        {"section": "wf_bloch", "member": "occupations", "binary_index": [1, 1, 0],
         "equals": 1.0, "tol": 1e-12},
        {"section": "wf_bloch", "member": "occupations", "binary_index": [1, 1, 1],
         "equals": 0.0, "tol": 1e-12},
        {"section": "bands", "member": "eigenvalues", "binary_index": [1, 1, 0],
         "equals": float(energies[1, 1, 0]) * _HARTREE_EV, "tol": 1e-9},
        {"section": "bands", "member": "eigenvalues", "binary_index": [0, 0, 1],
         "equals": float(energies[0, 0, 1]) * _HARTREE_EV, "tol": 1e-9},
        {"section": "bands", "member": "kpath", "json_pointer": ["labels", 1],
         "equals": "X"},
    ], "A dimerised hydrogen chain (pbc along x only) with wavefunction.bloch "
       "orbitals at Γ and X: weights 1/2 and 1/2, two spin channels with "
       "different on-site energies, coefficients float64 [n_k, n_spin, n_mo, "
       "n_ao, 2] that are real at Γ and carry the gauge phase e^{iπ/3} at X, "
       "energies in hartree and per-k occupations [n_k, n_spin, n_mo]. The "
       "orbitals satisfy C†S(k)C = 1 for the Bloch sum of spec § 5.3.1. The "
       "bands section lists the same two k-points, and its eigenvalues are "
       "the orbital energies in eV (spec § 5.3.1, § 5.6).")



# =============================================================================
# Sketch documents: the chemical identity of structure sections (spec § 5.1),
# reaction.scheme (§ 5.14), preview (§ 5.15), opaque vendor members with a
# media_type (§ 4.2, § 7.3) and job sets (§ 5.8, § 5.9).
# =============================================================================

def _png_image(width, height):
    """A deterministic grey-gradient PNG (8-bit greyscale)."""
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))
    raw = b"".join(b"\x00" + bytes((x * 255) // max(width - 1, 1)
                                   for x in range(width))
                   for _ in range(height))
    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


#: The species of the scheme: id -> (label, name, charge, multiplicity,
#: SMILES, InChI, InChIKey, atoms as (symbol, x, y, z) in angstrom). The
#: identifiers were computed with RDKit from the SMILES.
_SKETCH_SPECIES = {
    "s1": ("1", "methane", 0, 1, "C", "InChI=1S/CH4/h1H4",
           "VNWKTOKETHGBQD-UHFFFAOYSA-N",
           [("C", 0.0, 0.0, 0.0), ("H", 0.6291, 0.6291, 0.6291),
            ("H", -0.6291, -0.6291, 0.6291), ("H", -0.6291, 0.6291, -0.6291),
            ("H", 0.6291, -0.6291, -0.6291)]),
    "s2": ("2", "chlorine", 0, 1, "ClCl", "InChI=1S/Cl2/c1-2",
           "KZBUYRJDOAKODT-UHFFFAOYSA-N",
           [("Cl", 0.0, 0.0, 0.0), ("Cl", 0.0, 0.0, 1.988)]),
    "s3": ("3", "chloromethane", 0, 1, "CCl", "InChI=1S/CH3Cl/c1-2/h1H3",
           "NEHMKBQYUWJMIP-UHFFFAOYSA-N",
           [("C", 0.0, 0.0, 0.0), ("Cl", 0.0, 0.0, 1.781),
            ("H", 1.027, 0.0, -0.363), ("H", -0.5135, 0.8894, -0.363),
            ("H", -0.5135, -0.8894, -0.363)]),
    "s4": ("4", "hydrogen chloride", 0, 1, "Cl", "InChI=1S/ClH/h1H",
           "VEXZGXHMUGYJMC-UHFFFAOYSA-N",
           [("H", 0.0, 0.0, 0.0), ("Cl", 0.0, 0.0, 1.275)]),
    "s5": ("5", "chlorine atom", 0, 2, "[Cl]", "InChI=1S/Cl",
           "ZAMOUSCENKQFHK-UHFFFAOYSA-N", [("Cl", 0.0, 0.0, 0.0)]),
    "s6": ("6", "water", 0, 1, "O", "InChI=1S/H2O/h1H2",
           "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
           [("O", 0.0, 0.0, 0.1173), ("H", 0.0, 0.7572, -0.4692),
            ("H", 0.0, -0.7572, -0.4692)]),
    "s7": ("7", "oxonium", 1, 1, "[OH3+]", "InChI=1S/H2O/h1H2/p+1",
           "XLYOFNOQVPJJNP-UHFFFAOYSA-O",
           [("O", 0.0, 0.0, 0.0), ("H", 0.94, 0.0, -0.29),
            ("H", -0.47, 0.814, -0.29), ("H", -0.47, -0.814, -0.29)]),
    "s8": ("8", "chloride", -1, 1, "[Cl-]", "InChI=1S/ClH/h1H/p-1",
           "VEXZGXHMUGYJMC-UHFFFAOYSA-M", [("Cl", 0.0, 0.0, 0.0)]),
}

_SKETCH_SVG = (
    "<svg xmlns='http://www.w3.org/2000/svg' width='320' height='60' "
    "viewBox='0 0 320 60'><text x='10' y='35' font-family='sans-serif' "
    "font-size='16'>CH\u2084 + Cl\u2082 \u2192 CH\u2083Cl + HCl</text></svg>")


def _add_species(w, sid):
    label, name, charge, mult, smiles, inchi, key, atoms = _SKETCH_SPECIES[sid]
    return w.add_structure(
        [{"symbol": a, "position": [x, y, z]} for a, x, y, z in atoms],
        charge=charge, multiplicity=mult, name=name, smiles=smiles,
        inchi=inchi, inchikey=key, label=label, section_id=f"structure_{sid}")


def _sketch_documents() -> None:
    """reaction_scheme.qvf and job_set.qvf."""
    w = QvfWriter(_PROGRAM, _VERSION, "chlorination scheme")
    png = _png_image(32, 16)
    w.add_preview(png, image=_SKETCH_SVG, label="scheme", section_id="preview")
    for sid in _SKETCH_SPECIES:
        _add_species(w, sid)
    w.add_reaction_scheme([
        {"id": "r1", "label": "a",
         "reactants": ["structure_s1", "structure_s2"],
         "reagents": ["structure_s5"],
         "products": ["structure_s3", "structure_s4"],
         "conditions": ["h\u03bd", "25 \u00b0C"],
         "yield_percent": 61.5, "yield_text": "61.5%"},
        {"id": "r2", "label": "b", "arrow": "equilibrium",
         "reactants": ["structure_s4", "structure_s6"],
         "products": ["structure_s7", "structure_s8"],
         "conditions": ["water"]},
        {"id": "r3", "label": "c", "arrow": "retrosynthetic",
         "products": ["structure_s3"],
         "reactants": ["structure_s1", {"structure_ref": "structure_s2",
                                        "coefficient": 0.5}]},
    ], label="Radical chlorination of methane", section_id="scheme")
    w.set_extensions({"x_conformance": {"version": "1", "critical": False}})
    w.add_vendor_section(
        "x_conformance.drawing",
        json_members={"document": {"format": "x-conformance-drawing/1",
                                   "species": sorted(_SKETCH_SPECIES)}},
        opaque_members={"drawing": (_SKETCH_SVG.encode("utf-8"),
                                    "image/svg+xml")},
        label="the editable drawing", section_id="x_conformance_drawing")
    checks = [
        {"manifest_pointer": ["sections", 0, "kind"], "equals": "preview"},
        {"manifest_pointer": ["sections", 0, "members", "thumbnail", "media_type"],
         "equals": "image/png"},
        {"manifest_pointer": ["sections", 0, "members", "thumbnail", "size"],
         "equals": len(png)},
        {"manifest_pointer": ["sections", 0, "members", "image", "media_type"],
         "equals": "image/svg+xml"},
        {"section": "preview", "member": "image", "decode": "utf8",
         "contains": "<svg"},
        {"manifest_pointer": ["sections", 1, "id"], "equals": "structure_s1"},
        {"manifest_pointer": ["sections", 1, "label"], "equals": "1"},
        {"manifest_pointer": ["sections", 1, "name"], "equals": "methane"},
        {"manifest_pointer": ["sections", 1, "inchikey"],
         "equals": "VNWKTOKETHGBQD-UHFFFAOYSA-N"},
        {"manifest_pointer": ["sections", 5, "multiplicity"], "equals": 2},
        {"manifest_pointer": ["sections", 5, "smiles"], "equals": "[Cl]"},
        {"manifest_pointer": ["sections", 7, "charge"], "equals": 1},
        {"manifest_pointer": ["sections", 7, "inchi"],
         "equals": "InChI=1S/H2O/h1H2/p+1"},
        {"manifest_pointer": ["sections", 8, "charge"], "equals": -1},
        {"section": "structure_s3", "member": "structure",
         "json_pointer": ["atoms", 1, "symbol"], "equals": "Cl"},
        {"section": "structure_s3", "member": "structure",
         "json_pointer": ["atoms", 1, "position", 2], "equals": 1.781, "tol": 1e-12},
        {"manifest_pointer": ["sections", 9, "kind"], "equals": "reaction.scheme"},
        {"manifest_pointer": ["sections", 9, "label"],
         "equals": "Radical chlorination of methane"},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 0, "reactants", 1, "structure_ref"],
         "equals": "structure_s2"},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 0, "reagents", 0, "structure_ref"],
         "equals": "structure_s5"},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 0, "conditions"],
         "equals": ["h\u03bd", "25 \u00b0C"]},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 0, "yield_percent"], "equals": 61.5, "tol": 1e-12},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 1, "arrow"], "equals": "equilibrium"},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 2, "arrow"], "equals": "retrosynthetic"},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 2, "products", 0, "structure_ref"],
         "equals": "structure_s3"},
        {"section": "scheme", "member": "scheme",
         "json_pointer": ["steps", 2, "reactants", 1, "coefficient"],
         "equals": 0.5, "tol": 1e-12},
        {"manifest_pointer": ["extensions", "x_conformance", "version"],
         "equals": "1"},
        {"manifest_pointer": ["sections", 10, "members", "drawing", "media_type"],
         "equals": "image/svg+xml"},
        {"section": "x_conformance_drawing", "member": "drawing",
         "decode": "utf8", "contains": "CH\u2084 + Cl\u2082"},
    ]
    _write("reaction_scheme", w, checks,
           "A sketch saved as QVF. Eight species, each a structure section "
           "with its chemical identity on the section object (label as the "
           "compound number, name, charge, multiplicity, SMILES, InChI, "
           "InChIKey): a radical (the chlorine atom, multiplicity 2) and two "
           "ions (oxonium +1, chloride -1). A reaction.scheme relates them in "
           "three steps: a forward step with a reagent drawn at the arrow, "
           "conditions and a yield; an equilibrium; and a retrosynthetic step "
           "whose target is its product and whose precursors are its "
           "reactants, one with coefficient 0.5. A preview section, first in "
           "the manifest, carries a PNG thumbnail and an SVG image with their "
           "media_type, and a vendor section stores its drawing as an opaque "
           "member with a media_type (spec § 4.2, § 5.1, § 5.14, § 5.15, "
           "§ 7.3).")

    # A job set halfway through: two species, one job each; the first has
    # run, the second is pending (spec § 3.2, § 5.8, § 5.9).
    w = QvfWriter(_PROGRAM, _VERSION, "job set")
    w.set_provenance(run_status="running")
    _add_species(w, "s6")
    w.add_structure([{"symbol": "O", "position": [0.0, 0.0, 0.0]},
                     {"symbol": "H", "position": [0.0, 0.0, 0.971]}],
                    charge=0, multiplicity=2, name="hydroxyl radical",
                    smiles="[OH]", inchi="InChI=1S/HO/h1H",
                    inchikey="TUJKJAMUKRIRHC-UHFFFAOYSA-N", label="9",
                    section_id="structure_oh")
    w.add_job_spec(job_type="molecular", method="rks", functional="pbe0",
                   basis="def2-svp", charge=0, multiplicity=1,
                   tasks=["optimize", "hessian"], structure_ref="structure_s6",
                   label="water: opt + freq", section_id="job_water")
    w.add_job_spec(job_type="molecular", method="uks", functional="pbe0",
                   basis="def2-svp", charge=0, multiplicity=2,
                   tasks=["optimize"], structure_ref="structure_oh",
                   label="hydroxyl: opt", section_id="job_oh")
    w.add_run_record(program="qvf-conformance-runner", program_version="1.0",
                     input_text="opt freq pbe0 def2-svp\n",
                     log_text="water: optimisation converged\n",
                     exit_status=0, sequence=0, job_ref="job_water",
                     run_status="converged", section_id="run_water")
    w.add_atom_properties(mulliken_charge=[-0.66, 0.33, 0.33],
                          label="Mulliken charges", section_id="charges_water")
    w.annotate_section("charges_water", structure_ref="structure_s6",
                       run_ref="run_water", group="water")
    _write("job_set", w, [
        {"manifest_pointer": ["provenance", "run_status"], "equals": "running"},
        {"manifest_pointer": ["sections", 0, "multiplicity"], "equals": 1},
        {"manifest_pointer": ["sections", 1, "multiplicity"], "equals": 2},
        {"manifest_pointer": ["sections", 2, "kind"], "equals": "job.spec"},
        {"manifest_pointer": ["sections", 2, "structure_ref"],
         "equals": "structure_s6"},
        {"manifest_pointer": ["sections", 3, "structure_ref"],
         "equals": "structure_oh"},
        {"section": "job_water", "member": "spec",
         "json_pointer": ["structure_ref"], "equals": "structure_s6"},
        {"section": "job_oh", "member": "spec",
         "json_pointer": ["multiplicity"], "equals": 2},
        {"manifest_pointer": ["sections", 4, "job_ref"], "equals": "job_water"},
        {"manifest_pointer": ["sections", 4, "run_status"], "equals": "converged"},
        {"section": "run_water", "member": "log", "decode": "utf8",
         "contains": "converged"},
        {"manifest_pointer": ["sections", 5, "run_ref"], "equals": "run_water"},
        {"manifest_pointer": ["sections", 5, "structure_ref"],
         "equals": "structure_s6"},
        {"section": "charges_water", "member": "mulliken_charge",
         "binary_index": [0], "equals": -0.66, "tol": 1e-12},
    ], "A job set halfway through (spec § 5.9): two structures with their "
       "identity (water, and the hydroxyl radical with multiplicity 2), one "
       "job.spec for each naming its structure with a section-level "
       "structure_ref (and the same in the payload). The water job has run: "
       "its run.record names it with job_ref and states run_status "
       "converged, and the Mulliken charges it produced carry run_ref and "
       "structure_ref. No run.record names the hydroxyl job, so it is "
       "pending, and provenance.run_status of the set is running.")


# =============================================================================
# Toolkit 0.3.0: what the family's programs already write or read. Localized
# orbitals with their IAO reference basis (spec § 5.3), job.spec functional
# and basis null (§ 5.9), the deprecated spellings producers wrote before the
# canonical carriers (§ 5.3, § 5.9), structure_ref on standalone bonds
# sections (§ 5.1) and the viewer hints camera, crossfade and energy_window
# (§ 3.4).
# =============================================================================
import hashlib  # noqa: E402


class _Patched:
    """Adapts a finished ``QvfWriter`` to ``_write`` with some JSON payloads
    rewritten: fields the reference writer has no parameter for, and the
    deprecated spellings it never writes. ``updates`` maps ``(section id,
    role)`` to a function of the payload; the member's sha256 and size in the
    manifest follow, and the archive is re-packed in the writer's entry order
    with fixed time stamps."""

    def __init__(self, writer, updates):
        self._writer = writer
        self._updates = updates

    def _build(self):
        entries = _entries_of(self._writer.to_bytes())
        data = {name: body for name, body, _method in entries}
        manifest = json.loads(data["manifest.json"])
        for (sid, role), update in self._updates.items():
            section = next(s for s in manifest["sections"] if s["id"] == sid)
            member = section["members"][role]
            raw = (json.dumps(update(json.loads(data[member["path"]])),
                              indent=2, ensure_ascii=False) + "\n").encode("utf-8")
            data[member["path"]] = raw
            member["sha256"] = hashlib.sha256(raw).hexdigest()
            member["size"] = len(raw)
        data["manifest.json"] = (json.dumps(manifest, indent=2,
                                            ensure_ascii=False)
                                 + "\n").encode("utf-8")
        return manifest, _plain_zip([(name, data[name], method)
                                     for name, _body, method in entries])

    def write(self, path):
        with open(path, "wb") as fh:
            fh.write(self._build()[1])
        return path

    def build_manifest(self):
        return self._build()[0]


def _s_shells(n_atoms, exponent=1.0):
    return [{"center": a, "l": 0, "exponents": [exponent], "coefficients": [1.0]}
            for a in range(n_atoms)]


def _pair_coefficients(r_angstrom, exponent=1.0):
    """The coefficients 1 / sqrt(2 (1 + S)) and 1 / sqrt(2 (1 - S)) of the
    normalized bonding and antibonding combinations of two normalized s
    Gaussians of one exponent r_angstrom apart, whose overlap is
    S = exp(-exponent r^2 / 2), r in bohr."""
    r = r_angstrom * BOHR_PER_ANGSTROM
    overlap = math.exp(-0.5 * exponent * r * r)
    return (1.0 / math.sqrt(2.0 * (1.0 + overlap)),
            1.0 / math.sqrt(2.0 * (1.0 - overlap)))


#: Two hydrogen molecules 3 Å apart, each 0.74 Å long (Å).
_H2_PAIR = [("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.74),
            ("H", 3.0, 0.0, 0.0), ("H", 3.0, 0.0, 0.74)]


def _family_drift_archives() -> None:
    """wavefunction_localized, job_set_methods, deprecated_spellings,
    bonds_structure_ref and viewer_hints."""
    atoms = [{"symbol": s, "position": [x, y, z]} for s, x, y, z in _H2_PAIR]
    c, c_anti = _pair_coefficients(0.74)
    z_mid = 0.37 * BOHR_PER_ANGSTROM
    x_far = 3.0 * BOHR_PER_ANGSTROM

    # 1) Intrinsic bond orbitals with the IAO reference basis (spec § 5.3).
    w = QvfWriter(_PROGRAM, _VERSION, "two H2: intrinsic bond orbitals")
    w.add_structure(atoms)
    w.add_wavefunction_gto(
        _s_shells(4), mo_coefficients=[[c, c, 0.0, 0.0], [0.0, 0.0, c, c]],
        energies=[-0.59, -0.59], occupations=[2.0, 2.0], n_electrons=4,
        orbital_kind="ibo", label="IBOs", section_id="wf_ibo")
    w.add_atom_properties(iao_charge=[0.0, 0.0, 0.0, 0.0],
                          label="IAO charges", section_id="iao_charges")
    localized = _Patched(w, {("wf_ibo", "mo_metadata"): lambda m: dict(
        m, localization_method="ibo", localization_reference_basis="minao",
        atom_populations=[[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, 0.5]],
        centroids_bohr=[[0.0, 0.0, z_mid], [x_far, 0.0, z_mid]],
        n_centres=[2, 2])})
    _write("wavefunction_localized", localized, [
        {"section": "wf_ibo", "member": "mo_metadata",
         "json_pointer": ["orbital_kind"], "equals": "ibo"},
        {"section": "wf_ibo", "member": "mo_metadata",
         "json_pointer": ["localization_method"], "equals": "ibo"},
        {"section": "wf_ibo", "member": "mo_metadata",
         "json_pointer": ["localization_reference_basis"], "equals": "minao"},
        {"section": "wf_ibo", "member": "mo_metadata",
         "json_pointer": ["atom_populations", 1, 2], "equals": 0.5, "tol": 1e-12},
        {"section": "wf_ibo", "member": "mo_metadata",
         "json_pointer": ["centroids_bohr", 1, 0], "equals": x_far, "tol": 1e-9},
        {"section": "wf_ibo", "member": "mo_metadata",
         "json_pointer": ["n_centres"], "equals": [2, 2]},
        {"section": "wf_ibo", "member": "mo_coefficients",
         "binary_index": [1, 3], "equals": c, "tol": 1e-12},
        {"section": "wf_ibo", "member": "basis",
         "json_pointer": ["n_electrons"], "equals": 4},
        {"section": "iao_charges", "member": "iao_charge",
         "binary_index": [2], "equals": 0.0, "tol": 1e-12},
    ], "Two hydrogen molecules 3 Å apart with their intrinsic bond orbitals "
       "(orbital_kind ibo), one per molecule: each the normalized bonding "
       "combination of the molecule's two s functions, so atom_populations "
       "puts half of each orbital on each of its atoms, n_centres is 2 and "
       "the centroids are the bond midpoints in bohr. mo_metadata names the "
       "localization_method and the localization_reference_basis, the "
       "minimal basis the IAOs were built from (minao), and the IAO charges "
       "are the canonical iao_charge member of an atom_properties section "
       "(spec § 5.3, § 5.7). The orbital energies are illustrative.")

    # 2) A job set at three levels of theory, functional and basis null where
    # the method takes none (spec § 5.9).
    w = QvfWriter(_PROGRAM, _VERSION, "water: three methods")
    w.set_provenance(run_status="pending")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.1173]},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692]},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692]}],
        charge=0, multiplicity=1, name="water")
    w.add_job_spec(job_type="molecular", method="rhf", basis="def2-svp",
                   structure_ref="structure", label="Hartree-Fock",
                   section_id="job_rhf")
    w.add_job_spec(job_type="molecular", method="gfn2_xtb",
                   structure_ref="structure", label="GFN2-xTB",
                   section_id="job_xtb")
    w.add_job_spec(job_type="molecular", method="rks", functional="pbe0",
                   basis="def2-svp", structure_ref="structure",
                   label="PBE0", section_id="job_pbe0")
    methods = _Patched(w, {
        ("job_rhf", "spec"): lambda p: dict(p, functional=None),
        ("job_xtb", "spec"): lambda p: dict(p, functional=None, basis=None),
    })
    _write("job_set_methods", methods, [
        {"manifest_pointer": ["provenance", "run_status"], "equals": "pending"},
        {"section": "job_rhf", "member": "spec",
         "json_pointer": ["functional"], "equals": None},
        {"section": "job_rhf", "member": "spec",
         "json_pointer": ["basis"], "equals": "def2-svp"},
        {"section": "job_xtb", "member": "spec",
         "json_pointer": ["method"], "equals": "gfn2_xtb"},
        {"section": "job_xtb", "member": "spec",
         "json_pointer": ["functional"], "equals": None},
        {"section": "job_xtb", "member": "spec",
         "json_pointer": ["basis"], "equals": None},
        {"section": "job_pbe0", "member": "spec",
         "json_pointer": ["functional"], "equals": "pbe0"},
        {"manifest_pointer": ["sections", 3, "structure_ref"],
         "equals": "structure"},
    ], "A pending job set on water at three levels of theory (spec § 5.9): "
       "Hartree-Fock with functional null, GFN2-xTB with functional and "
       "basis null, the method taking neither, and PBE0/def2-SVP. null "
       "states that the method takes no such input; an absent key would "
       "leave the choice to the engine.")

    # 3) The deprecated spellings a producer wrote before the canonical
    # carriers: valid, with no finding (spec § 5.3, § 5.9).
    w = QvfWriter(_PROGRAM, _VERSION, "deprecated spellings")
    w.add_structure(atoms)
    w.add_wavefunction_gto(
        _s_shells(4), mo_coefficients=[[c, c, 0.0, 0.0], [0.0, 0.0, c, c]],
        energies=[-0.59, -0.59], occupations=[2.0, 2.0],
        orbital_kind="localized", section_id="wf_localized_ibo")
    for role, rows, weight in (("hole", [[c, c, 0.0, 0.0]], 1.0),
                               ("electron", [[c_anti, -c_anti, 0.0, 0.0]],
                                1.0)):
        w.add_wavefunction_gto(
            _s_shells(4), mo_coefficients=rows, energies=[0.0],
            occupations=[weight], orbital_kind="natural",
            section_id=f"wf_nto_S1_{role}")
    w.add_job_spec(job_type="molecular", method="gfn2_xtb",
                   structure_ref="structure", section_id="job_spec")

    def nto(meta):
        return dict(meta, occupation_semantics="transition_weight",
                    excitation_index=1, excitation_energy_ev=12.85)

    deprecated = _Patched(w, {
        ("wf_localized_ibo", "mo_metadata"): lambda m: dict(
            m, localization_method="ibo", localization_reference_basis="mini",
            atom_populations=[[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, 0.5]],
            centroids_bohr=[[0.0, 0.0, z_mid], [x_far, 0.0, z_mid]],
            n_centres=[2, 2], iao_charges=[0.0, 0.0, 0.0, 0.0]),
        ("wf_nto_S1_hole", "mo_metadata"): nto,
        ("wf_nto_S1_electron", "mo_metadata"): nto,
        ("job_spec", "spec"): lambda p: dict(p, basis="", functional=None),
    })
    _write("deprecated_spellings", deprecated, [
        {"section": "wf_localized_ibo", "member": "mo_metadata",
         "json_pointer": ["iao_charges"], "equals": [0.0, 0.0, 0.0, 0.0]},
        {"section": "wf_localized_ibo", "member": "mo_metadata",
         "json_pointer": ["localization_reference_basis"], "equals": "mini"},
        {"section": "wf_nto_S1_hole", "member": "mo_metadata",
         "json_pointer": ["excitation_index"], "equals": 1},
        {"section": "wf_nto_S1_electron", "member": "mo_metadata",
         "json_pointer": ["excitation_energy_ev"], "equals": 12.85,
         "tol": 1e-12},
        {"section": "wf_nto_S1_electron", "member": "mo_metadata",
         "json_pointer": ["occupation_semantics"],
         "equals": "transition_weight"},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["basis"], "equals": ""},
        {"section": "job_spec", "member": "spec",
         "json_pointer": ["functional"], "equals": None},
    ], "The deprecated spellings a producer wrote before the canonical "
       "carriers existed, which remain valid with no finding: localized "
       "orbitals whose mo_metadata carries iao_charges (canonical: the "
       "iao_charge member of atom_properties) beside "
       "localization_reference_basis; the hole and electron natural "
       "transition orbitals of a state with excitation_index and "
       "excitation_energy_ev (canonical: state_ref with state_index, and the "
       "state's excitation_energy in excited_states); and a semiempirical "
       "job.spec with basis \"\", the deprecated spelling of null "
       "(spec § 5.3, § 5.9).")

    # 4) Two structures and a standalone bonds section for each, the second
    # naming its structure (spec § 5.1).
    w = QvfWriter(_PROGRAM, _VERSION, "water and hydrogen: bonds")
    w.add_structure(
        [{"symbol": "O", "position": [0.0, 0.0, 0.1173]},
         {"symbol": "H", "position": [0.0, 0.7572, -0.4692]},
         {"symbol": "H", "position": [0.0, -0.7572, -0.4692]}],
        name="water", section_id="structure_water")
    w.add_structure([{"symbol": "H", "position": [0.0, 0.0, 0.0]},
                     {"symbol": "H", "position": [0.0, 0.0, 0.74]}],
                    name="hydrogen", section_id="structure_h2")
    w.add_bonds([{"i": 0, "j": 1}, {"i": 0, "j": 2}], label="O-H bonds",
                section_id="bonds_water")
    w.add_bonds([{"i": 0, "j": 1, "order": 1.0}], structure_ref="structure_h2",
                label="H-H bond", section_id="bonds_h2")
    _write("bonds_structure_ref", w, [
        {"manifest_pointer": ["sections", 2, "id"], "equals": "bonds_water"},
        {"manifest_pointer": ["sections", 3, "structure_ref"],
         "equals": "structure_h2"},
        {"section": "bonds_water", "member": "bonds",
         "json_pointer": ["pairs", 1, "j"], "equals": 2},
        {"section": "bonds_h2", "member": "bonds",
         "json_pointer": ["pairs", 0, "j"], "equals": 1},
        {"section": "structure_h2", "member": "structure",
         "json_pointer": ["atoms", 1, "position", 2], "equals": 0.74,
         "tol": 1e-12},
    ], "Two structures, water and then hydrogen, and two standalone bonds "
       "sections. bonds_h2 names its structure with structure_ref, so its "
       "indices address the two hydrogen atoms; bonds_water has none, so its "
       "indices address the archive's first structure, water, whose atom 2 "
       "the second pair names (spec § 4.1, § 5.1).")

    # 5) The viewer hints camera, crossfade and energy_window (spec § 3.4).
    lattice = [[3.8, 0.0, 0.0], [0.0, 3.8, 0.0], [0.0, 0.0, 3.8]]
    fermi = -4.5
    eig = np.array([[[-1839.2, -9.0, -2.5], [-1839.2, -7.5, -1.0],
                     [-1839.2, -6.0, 0.5], [-1839.2, -5.0, 1.5]]],
                   dtype=np.float64)                         # [1, 4, 3], eV
    e_dos = np.array([-1836.0, -1834.7, -10.0, -5.0, 0.0, 5.0, 10.0])
    dos = np.array([0.0, 2.0, 0.1, 0.8, 0.2, 0.6, 0.1])
    shape, step = (4, 4, 4), 0.5
    origin = [-0.75, -0.75, -0.75]
    voxels = [[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, step]]
    xyz = np.stack(np.meshgrid(*[origin[a] + step * np.arange(shape[a])
                                 for a in range(3)], indexing="ij"), axis=-1)
    rho_a = np.exp(-0.8 * np.sum(xyz ** 2, axis=-1)).astype(np.float32)
    rho_b = (0.5 * np.exp(-0.4 * np.sum(xyz ** 2, axis=-1))).astype(np.float32)
    camera = {"position": [0.0, -12.0, 4.0], "focal_point": [1.9, 1.9, 1.9],
              "view_up": [0.0, 0.0, 1.0], "view_angle": 30.0}
    w = QvfWriter(_PROGRAM, _VERSION, "silicon: viewer hints")
    w.add_structure([{"symbol": "Si", "position": [0.0, 0.0, 0.0]}],
                    pbc=[True, True, True], lattice_vectors=lattice)
    w.add_bands({"kpoints": [[0.0, 0.0, 0.0], [1 / 6, 0.0, 0.0],
                             [1 / 3, 0.0, 0.0], [0.5, 0.0, 0.0]],
                 "labels": ["G", "", "", "X"], "fermi": fermi,
                 "n_spin": 1, "n_kpoints": 4}, eig, section_id="bands")
    w.add_dos_total(e_dos, dos, fermi_energy_ev=fermi, section_id="dos")
    w.add_density(rho_a, origin, voxels, label="rho(A)", section_id="rho_a")
    w.add_density(rho_b, origin, voxels, label="rho(B)", section_id="rho_b")
    w.set_viewer_defaults({
        "auto_open": ["structure", "bands"],
        "camera": camera,
        "crossfade": {"a": "rho_a", "b": "rho_b", "blend": 0.25},
        "bands": {"energy_window": [-15.0, 10.0]},
        "dos": {"energy_window": [-15.0, 10.0]},
        "rho_a": {"isovalue": 0.05},
    })
    _write("viewer_hints", w, [
        {"manifest_pointer": ["viewer_defaults", "camera"], "equals": camera},
        {"manifest_pointer": ["viewer_defaults", "crossfade", "a"],
         "equals": "rho_a"},
        {"manifest_pointer": ["viewer_defaults", "crossfade", "blend"],
         "equals": 0.25, "tol": 1e-12},
        {"manifest_pointer": ["viewer_defaults", "bands", "energy_window"],
         "equals": [-15.0, 10.0]},
        {"manifest_pointer": ["viewer_defaults", "dos", "energy_window", 0],
         "equals": -15.0, "tol": 1e-12},
        {"section": "bands", "member": "kpath", "json_pointer": ["fermi"],
         "equals": fermi, "tol": 1e-12},
        {"section": "bands", "member": "eigenvalues", "binary_index": [0, 3, 0],
         "equals": -1839.2, "tol": 1e-9},
        {"section": "dos", "member": "energies", "binary_index": [1],
         "equals": -1834.7, "tol": 1e-9},
        {"section": "rho_b", "member": "data", "binary_index": [1, 2, 3],
         "equals": float(rho_b[1, 2, 3]), "tol": 1e-6},
    ], "An all-electron silicon cell that opens at a saved view (spec § 3.4): "
       "the camera hint (the shape of a bookmark's camera), a crossfade "
       "between two densities starting at blend 0.25 (0 shows a, 1 shows "
       "b), and an energy_window of [-15, 10] eV relative to the Fermi "
       "level on the bands and the DOS, whose Si 1s states 1835 eV below it "
       "would otherwise squeeze the valence region into a few pixels. The "
       "hints are specified in the spec only, so a malformed one never fails "
       "an archive; these are well-formed.")


def build_unrestricted_ntos():
    """Issue 12: two sections, per-spin pairing, conditional normalization."""
    shells = [{"center": 0, "l": 0, "exponents": [1.0], "coefficients": [1.0]},
              {"center": 1, "l": 0, "exponents": [1.0], "coefficients": [1.0]}]
    cases = [("nto_unrestricted", [0.55, 0.15], [0.30]),
             ("nto_unrestricted_empty_beta", [0.8, 0.2], []),
             ("nto_unrestricted_tddft", [0.7, 0.2], [0.3]),
             ("nto_unrestricted_legacy", [0.55, 0.15], [0.3])]
    for name, alpha, beta in cases:
        w = QvfWriter(_PROGRAM, _VERSION, "unrestricted transition orbitals")
        w.add_structure([{"symbol": "H", "atomic_number": 1, "position": [0, 0, 0]},
                         {"symbol": "H", "atomic_number": 1, "position": [0, 0, 0.74]}])
        # Deliberately write the particle before the hole and the state last.
        for role in ("particle", "hole"):
            w.add_wavefunction_gto(shells, mo_coefficients_alpha=np.eye(2),
                mo_coefficients_beta=np.eye(2)[:len(beta)],
                occupations_alpha=alpha, occupations_beta=beta, n_electrons=2,
                orbital_kind="nto", state_ref="es", state_index=1,
                nto_role=role, section_id="nto_" + role)
        w.add_excited_states([{"index": 1, "excitation_energy": 4.0, "nto_ref": "nto_hole"}], section_id="es")
        checks = []
        for role in ("hole", "particle"):
            sid = "nto_" + role
            for pointer, value in [(["spin"], "unrestricted"), (["nto_role"], role),
                                   (["occupation_semantics"], "transition_weight"),
                                   (["alpha", "occupations"], alpha), (["beta", "occupations"], beta),
                                   (["alpha", "n_mo"], 2), (["beta", "n_mo"], len(beta))]:
                checks.append({"section": sid, "member": "mo_metadata", "json_pointer": pointer, "equals": value})
            checks.append({"section": sid, "member": "mo_coefficients_alpha", "binary_index": [0, 0], "equals": 1.0})
        output = w
        if name == "nto_unrestricted_legacy":
            def legacy(meta):
                return {**{k: v for k, v in meta.items() if k not in ("state_ref", "state_index", "nto_role")},
                        "excitation_index": 1, "excitation_energy_ev": 4.0}
            output = _Patched(w, {("nto_" + role, "mo_metadata"): legacy for role in ("hole", "particle")})
            checks = [c for c in checks if c.get("json_pointer") != ["nto_role"]]
            checks.append({"section": "nto_hole", "member": "mo_metadata", "json_pointer": ["excitation_index"], "equals": 1})
        _write(name, output, checks, "Unrestricted NTO hole/particle sections with two alpha rows and "
               + str(len(beta)) + " beta rows. Weights pair per spin; their combined sum is "
               + str(round(sum(alpha) + sum(beta), 4)) + ". Full TDDFT/RPA sums need not be one (spec § 5.3).")


def build_biomolecule_identity():
    """structure_biomolecule: residues 10 and 10A, alternate locations (#18)."""
    def atom(symbol, position, name, residue, seq, **identity):
        return {"symbol": symbol, "position": position, "atom_name": name,
                "residue_name": residue, "residue_seq": seq, "chain_id": "A",
                **identity}
    atoms = [
        atom("C", [10.512, 4.210, 2.033], "CA", "SER", 10, occupancy=1.0, hetero=False),
        atom("O", [11.402, 5.318, 1.801], "OG", "SER", 10, alt_loc="A", occupancy=0.6,
             hetero=False, b_factor=18.5),
        atom("O", [9.730, 5.104, 2.987], "OG", "SER", 10, alt_loc="B", occupancy=0.4,
             hetero=False, b_factor=21.0),
        atom("C", [13.318, 3.904, 3.115], "CA", "GLY", 10, insertion_code="A",
             occupancy=1.0, hetero=False),
        atom("O", [15.020, 7.441, 0.512], "O", "HOH", 101, occupancy=0.5, hetero=True),
    ]
    w = QvfWriter(_PROGRAM, _VERSION, "biomolecular atom-site identity")
    w.add_structure(atoms, biomolecule={
        "chains": ["A"],
        "residues": [
            {"name": "SER", "seq": 10, "chain": "A", "atom_indices": [0, 1, 2]},
            {"name": "GLY", "seq": 10, "chain": "A", "insertion_code": "A",
             "atom_indices": [3]},
            {"name": "HOH", "seq": 101, "chain": "A", "atom_indices": [4]},
        ]}, name="Ser10-Gly10A fragment with water")
    _write("structure_biomolecule", w, [
        {"manifest_pointer": ["sections", 0, "residues", 0, "seq"], "equals": 10},
        {"manifest_pointer": ["sections", 0, "residues", 1, "seq"], "equals": 10},
        {"manifest_pointer": ["sections", 0, "residues", 1, "insertion_code"],
         "equals": "A"},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 3, "insertion_code"], "equals": "A"},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 1, "alt_loc"], "equals": "A"},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 2, "alt_loc"], "equals": "B"},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 1, "occupancy"], "equals": 0.6, "tol": 1e-12},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 2, "occupancy"], "equals": 0.4, "tol": 1e-12},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 2, "hetero"], "equals": False},
        {"section": "structure", "member": "structure",
         "json_pointer": ["atoms", 4, "hetero"], "equals": True},
    ], "Chain A residues SER 10 and GLY 10A, distinct because their insertion "
       "codes differ, and a HETATM water. SER 10 OG has alternate locations A "
       "(occupancy 0.6) and B (0.4) as two atoms; the CA without alt_loc "
       "belongs to both conformers (spec § 5.1, Biomolecular atom identity).")



def build_secondary_structure_ranges():
    """structure_secondary_structure: typed range ends and classes (#28)."""
    def ca(x, seq, **identity):
        return {"symbol": "C", "position": [x, 0.0, 0.0], "atom_name": "CA",
                "residue_name": "ALA", "residue_seq": seq, "chain_id": "A",
                "hetero": False, **identity}
    atoms = [ca(0.0, 10), ca(3.8, 10, insertion_code="A"), ca(7.6, 11),
             ca(11.4, 12), ca(15.2, 13)]
    w = QvfWriter(_PROGRAM, _VERSION, "secondary-structure ranges")
    w.add_structure(atoms, biomolecule={
        "chains": ["A"],
        "residues": [
            {"name": "ALA", "seq": 10, "chain": "A", "atom_indices": [0]},
            {"name": "ALA", "seq": 10, "chain": "A", "insertion_code": "A",
             "atom_indices": [1]},
            {"name": "ALA", "seq": 11, "chain": "A", "atom_indices": [2]},
            {"name": "ALA", "seq": 12, "chain": "A", "atom_indices": [3]},
            {"name": "ALA", "seq": 13, "chain": "A", "atom_indices": [4]},
        ],
        "secondary_structure": [
            {"type": "helix", "chain": "A", "start_seq": 10, "start_insertion_code": "A",
             "end_seq": 11, "helix_class": 5},
            {"type": "sheet", "chain": "A", "start_seq": 12, "end_seq": 14,
             "sheet_id": "S1", "strand": 1},
            {"type": "coil", "chain": "A", "start_seq": 10, "end_seq": 10},
        ]}, name="alanine CA trace with typed ranges")
    _write("structure_secondary_structure", w, [
        {"manifest_pointer": ["sections", 0, "secondary_structure", 0,
                              "start_insertion_code"], "equals": "A"},
        {"manifest_pointer": ["sections", 0, "secondary_structure", 0, "helix_class"],
         "equals": 5},
        {"manifest_pointer": ["sections", 0, "secondary_structure", 1, "sheet_id"],
         "equals": "S1"},
        {"manifest_pointer": ["sections", 0, "secondary_structure", 1, "strand"],
         "equals": 1},
        {"manifest_pointer": ["sections", 0, "secondary_structure", 1, "end_seq"],
         "equals": 14},
    ], "A 3-10 helix (helix_class 5) from residue 10A to 11, distinct from residue "
       "10 by its start insertion code; strand 1 of sheet S1 from 12 to 14, whose "
       "end residue 14 the structure does not contain, which is allowed; and a "
       "coil at residue 10 (spec § 5.1, Secondary-structure ranges).")

def build_reaction_energies():
    """reaction_scheme_energies: 2 H2 + O2 -> 2 H2O with its energies (#16)."""
    import hashlib
    w = QvfWriter(_PROGRAM, _VERSION, "reaction energies")
    species = {
        "h2": [("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.7414)],
        "o2": [("O", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.2075)],
        "water": [("O", 0.0, 0.0, 0.1173), ("H", 0.0, 0.7572, -0.4692),
                  ("H", 0.0, -0.7572, -0.4692)],
    }
    multiplicity = {"h2": 1, "o2": 3, "water": 1}
    for sid, atoms in species.items():
        w.add_structure([{"symbol": a, "position": [x, y, z]} for a, x, y, z in atoms],
                        charge=0, multiplicity=multiplicity[sid], name=sid,
                        section_id=f"structure_{sid}")
    w.add_run_record(program="democode", program_version="1.0.0",
                     command="democode water.inp", exit_status=0,
                     log_text="FINAL ENERGY -76.0107\nFINAL G -76.0007\n",
                     section_id="run_water")
    method = {"program": "democode", "program_version": "1.0.0",
              "specification": {"task": "opt+freq", "method": "HF",
                                "basis": "6-31G*"}}
    absolute = {"h2": (-1.1269, -1.1279), "o2": (-149.6131, -149.6361),
                "water": (-76.0107, -76.0007)}
    coefficients = {"h2": -2, "o2": -1, "water": 2}
    terms = []
    for sid, c in coefficients.items():
        energy, gibbs = absolute[sid]
        term = {"structure_ref": f"structure_{sid}", "coefficient": c,
                "energy_eh": energy, "gibbs_free_energy_eh": gibbs,
                "method": method, "job_id": f"job_{sid}", "run_id": f"run_{sid}"}
        if sid == "water":
            term["run_ref"] = "run_water"
        else:
            term["source"] = {
                "name": f"{sid}.out",
                "sha256": hashlib.sha256(f"{sid} output".encode()).hexdigest(),
                "parser": "democode-text/1",
                "input_sha256": hashlib.sha256(f"{sid} input".encode()).hexdigest(),
                "locations": {"energy_eh": "line 120", "gibbs_free_energy_eh": "line 388",
                              "temperature_k": "line 371", "pressure_atm": "line 372"}}
        terms.append(term)
    deltas = {}
    for key, field in (("energy", "energy_eh"), ("gibbs_free_energy", "gibbs_free_energy_eh")):
        value = math.fsum(t["coefficient"] * t[field] for t in terms)
        deltas[key] = {"value_eh": value, "kcal_mol": value * 627.509474}
    w.add_reaction_scheme([{
        "id": "r1", "label": "combustion",
        "reactants": [{"structure_ref": "structure_h2", "coefficient": 2},
                      "structure_o2"],
        "products": [{"structure_ref": "structure_water", "coefficient": 2}],
        "energies": {"version": 1, "deltas": deltas, "temperature_k": 298.15,
                     "pressure_atm": 1.0, "terms": terms}}],
        label="hydrogen combustion")
    sid = "reaction_scheme"
    _write("reaction_scheme_energies", w, [
        {"section": sid, "member": "scheme",
         "json_pointer": ["steps", 0, "energies", "version"], "equals": 1},
        {"section": sid, "member": "scheme",
         "json_pointer": ["steps", 0, "energies", "deltas", "energy", "value_eh"],
         "equals": deltas["energy"]["value_eh"], "tol": 1e-12},
        {"section": sid, "member": "scheme",
         "json_pointer": ["steps", 0, "energies", "deltas", "gibbs_free_energy",
                          "kcal_mol"],
         "equals": deltas["gibbs_free_energy"]["kcal_mol"], "tol": 1e-9},
        {"section": sid, "member": "scheme",
         "json_pointer": ["steps", 0, "energies", "temperature_k"], "equals": 298.15},
        {"section": sid, "member": "scheme",
         "json_pointer": ["steps", 0, "energies", "terms", 0, "coefficient"],
         "equals": -2},
        {"section": sid, "member": "scheme",
         "json_pointer": ["steps", 0, "energies", "terms", 2, "run_ref"],
         "equals": "run_water"},
        {"section": sid, "member": "scheme",
         "json_pointer": ["steps", 0, "energies", "terms", 1, "source", "parser"],
         "equals": "democode-text/1"},
    ], "2 H2 + O2 -> 2 H2O with energy and Gibbs free energy deltas at 298.15 K "
       "and 1 atm: each delta is the sum of the signed coefficients (-2, -1, "
       "+2) times the absolute values of its terms, all from one method. The "
       "water term names its run.record with run_ref; the others name the "
       "parsed output files with source (spec § 5.14, Reaction energies).")


def build_compound_units():
    """Compound unit names of spec § 4.3, declared where the format fixes the
    unit (they restate it) and where a payload leaves it open (they choose
    it). Each declaration changes no value: the numbers are what the same
    archive without the declarations would hold."""
    water = [{"symbol": "O", "position": [0.0, 0.0, 0.1173], "atomic_number": 8},
             {"symbol": "H", "position": [0.0, 0.7572, -0.4692], "atomic_number": 1},
             {"symbol": "H", "position": [0.0, -0.7572, -0.4692], "atomic_number": 1}]
    origin = [-1.0, -1.0, -1.0]
    step = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    rho = np.full((3, 3, 3), 0.05)
    rho[1, 1, 1] = 0.35
    lap = np.full((3, 3, 3), 0.02)
    lap[1, 1, 1] = -1.2
    pot = np.full((3, 3, 3), 0.1)
    pot[1, 1, 1] = 0.9
    field = np.zeros((3, 3, 3, 3))
    field[2, 1, 1] = [0.04, 0.0, 0.0]
    coords = np.array([[[a["position"][c] for c in range(3)] for a in water]] * 2)
    coords[1, 0, 2] += 0.002
    velocities = np.zeros((2, 3, 3))
    velocities[1, 0, 2] = 0.001                       # Angstrom/fs
    gradient = np.zeros((2, 3, 3))
    gradient[0, 0, 2] = 0.0123                        # hartree/bohr
    gradient[1, 0, 2] = -0.0047
    quad = np.diag([-1.86, 1.96, -0.10])              # debye Angstrom
    alpha = np.diag([8.64, 9.15, 9.45])               # atomic units

    w = QvfWriter(_PROGRAM, _VERSION, "units-compound-molecular")
    w.add_structure(water)
    w.add_density(rho, origin, step, dtype="float64",
                  units={"origin": "bohr", "data": "e_per_bohr3"},
                  section_id="rho")
    w.add_volume("volume.laplacian", lap, origin, step, dtype="float64",
                 units={"data": "e_per_bohr5"}, section_id="lap")
    w.add_volume("volume.potential", pot, origin, step, dtype="float64",
                 units={"data": "hartree_per_e"}, section_id="esp")
    w.add_volume("volume.vector", field, origin, step, dtype="float64",
                 units={"data": "hartree_per_e_bohr"}, section_id="efield")
    w.add_trajectory({"symbols": ["O", "H", "H"], "kind": "md",
                      "timestep_fs": 0.5}, coords,
                     velocities=velocities, gradient=gradient,
                     times=[0.0, 0.5],
                     units={"coords": "angstrom", "times": "fs",
                            "velocities": "angstrom_per_fs",
                            "gradient": "hartree_per_bohr"},
                     section_id="md")
    w.add_topology_qtaim({"points": [
        {"type": "bcp", "position": [0.0, 1.0911, -0.6762], "rho": 0.3653,
         "laplacian": -2.5361, "atom_pair": [0, 1]}]},
        units={"position": "bohr", "rho": "e_per_bohr3",
               "laplacian": "e_per_bohr5"}, section_id="aim")
    w.add_electric_properties({
        "origin": "center_of_mass",
        "quadrupole": {"matrix": quad, "traceless": True,
                       "convention": "buckingham"},
        "polarizability": {"matrix": alpha, "frequency_ev": 0.0},
        "units": {"quadrupole": "debye_angstrom",
                  "polarizability": "e2_bohr2_per_hartree"}})
    w.add_atom_properties(custom={"gradient_ev_angstrom": (
        [[0.0, 0.0, 0.25], [0.0, 0.12, -0.125], [0.0, -0.12, -0.125]],
        {"units": "ev_per_angstrom", "kind": "vector",
         "description": "Energy gradient dE/dR from a producer that works in "
                        "eV and Angstrom, with illustrative values."})})
    _write("units_compound_molecular", w, [
        {"section": "rho", "member": "grid",
         "json_pointer": ["units", "data"], "equals": "e_per_bohr3"},
        {"section": "rho", "member": "data", "binary_index": [1, 1, 1],
         "equals": 0.35, "tol": 1e-12},
        {"section": "lap", "member": "grid",
         "json_pointer": ["units", "data"], "equals": "e_per_bohr5"},
        {"section": "lap", "member": "data", "binary_index": [1, 1, 1],
         "equals": -1.2, "tol": 1e-12},
        {"section": "esp", "member": "grid",
         "json_pointer": ["units", "data"], "equals": "hartree_per_e"},
        {"section": "efield", "member": "grid",
         "json_pointer": ["units", "data"], "equals": "hartree_per_e_bohr"},
        {"section": "efield", "member": "data", "binary_index": [2, 1, 1, 0],
         "equals": 0.04, "tol": 1e-12},
        {"section": "md", "member": "metadata",
         "json_pointer": ["units", "velocities"], "equals": "angstrom_per_fs"},
        {"section": "md", "member": "metadata",
         "json_pointer": ["units", "gradient"], "equals": "hartree_per_bohr"},
        {"section": "md", "member": "gradient", "binary_index": [0, 0, 2],
         "equals": 0.0123, "tol": 1e-12},
        {"section": "md", "member": "velocities", "binary_index": [1, 0, 2],
         "equals": 0.001, "tol": 1e-12},
        {"section": "aim", "member": "critical_points",
         "json_pointer": ["units", "rho"], "equals": "e_per_bohr3"},
        {"section": "aim", "member": "critical_points",
         "json_pointer": ["units", "laplacian"], "equals": "e_per_bohr5"},
        {"section": "aim", "member": "critical_points",
         "json_pointer": ["points", 0, "laplacian"], "equals": -2.5361,
         "tol": 1e-12},
        {"section": "electric", "member": "properties",
         "json_pointer": ["units", "quadrupole"], "equals": "debye_angstrom"},
        {"section": "electric", "member": "properties",
         "json_pointer": ["units", "polarizability"],
         "equals": "e2_bohr2_per_hartree"},
        {"section": "electric", "member": "properties",
         "json_pointer": ["quadrupole", "matrix", 1, 1], "equals": 1.96,
         "tol": 1e-12},
        {"section": "atom_properties", "member": "meta",
         "json_pointer": ["custom.gradient_ev_angstrom", "units"],
         "equals": "ev_per_angstrom"},
        {"section": "atom_properties", "member": "custom.gradient_ev_angstrom",
         "binary_index": [0, 2], "equals": 0.25, "tol": 1e-12},
    ], "Compound unit names (spec § 4.3) on a molecule: a density grid in "
       "e_per_bohr3, a Laplacian grid in e_per_bohr5, a potential in "
       "hartree_per_e and an electric field (volume.vector) in "
       "hartree_per_e_bohr. Trajectory velocities and gradient restate "
       "angstrom_per_fs and hartree_per_bohr, QTAIM rho and laplacian are "
       "declared in e_per_bohr3 and e_per_bohr5, a quadrupole restates "
       "debye_angstrom and a polarizability e2_bohr2_per_hartree, "
       "and a custom atom_properties role in ev_per_angstrom. Every value is "
       "stored in its declared unit and nothing is converted.")

    a = 4.12
    c11, c12, c44 = 36.4, 8.8, 8.0                    # GPa
    stiffness = np.zeros((6, 6))
    for i in range(3):
        stiffness[i, i] = c11
        stiffness[3 + i, 3 + i] = c44
        for j in range(3):
            if i != j:
                stiffness[i, j] = c12
    compliance = np.linalg.inv(stiffness)             # 1/GPa
    piezo = np.zeros((3, 6))
    piezo[0, 4] = piezo[1, 3] = 0.15                  # C/m^2
    w = QvfWriter(_PROGRAM, _VERSION, "units-compound-periodic")
    w.add_structure(
        [{"symbol": "Cs", "position": [0.0, 0.0, 0.0], "atomic_number": 55},
         {"symbol": "Cl", "position": [a / 2, a / 2, a / 2], "atomic_number": 17}],
        pbc=[True, True, True],
        lattice_vectors=[[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    w.add_elastic_properties(stiffness, compliance=compliance,
                             units={"stiffness": "gpa", "compliance": "gpa-1"})
    w.add_dielectric_properties(np.eye(3) * 2.62, piezoelectric=piezo,
                                units={"electronic": "dimensionless",
                                       "piezoelectric": "c_per_m2"})
    _write("units_compound_periodic", w, [
        {"section": "elastic", "member": "properties",
         "json_pointer": ["units", "compliance"], "equals": "gpa-1"},
        {"section": "elastic", "member": "properties",
         "json_pointer": ["compliance", 0, 0], "equals": float(compliance[0, 0]),
         "tol": 1e-12},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["units", "piezoelectric"], "equals": "c_per_m2"},
        {"section": "dielectric", "member": "properties",
         "json_pointer": ["piezoelectric", 0, 4], "equals": 0.15, "tol": 1e-12},
    ], "Compound unit names (spec § 4.3) on a crystal: an elastic compliance "
       "restating gpa-1 and a piezoelectric tensor restating c_per_m2, the "
       "units spec § 5.13 fixes for them.")


def build_xray_mossbauer():
    """spectra_xas_tddft, spectra_xas_dscf, spectra_mossbauer: the X-ray
    absorption and emission and the Mössbauer kinds (spec § 5.4, format
    review P13). The numbers are synthetic but have the magnitudes of real
    calculations: an Fe K pre-edge of tetrahedral [FeCl4]2- from three TDDFT
    core excitations named through states_ref, its Kβ2,5 valence-to-core
    emission, Cl K-edge ΔSCF sticks of two absorbers as cross sections, an Fe
    L2,3 edge in hartree (PySCF style), and two Fe sites of a mixed-valence
    pair with their isomer shifts, quadrupole splittings, EFG tensors and
    hyperfine field."""
    d = 2.29 / math.sqrt(3.0)
    fecl4 = [{"symbol": "Fe", "position": [0.0, 0.0, 0.0], "atomic_number": 26}]
    for sx, sy, sz in ((1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)):
        fecl4.append({"symbol": "Cl", "position": [sx * d, sy * d, sz * d],
                      "atomic_number": 17})

    # 1) Fe K pre-edge from TDDFT (ORCA OrbWin / Q-Chem TRNSS style).
    w = QvfWriter(_PROGRAM, _VERSION, "FeCl4 K-edge TDDFT")
    w.add_structure(fecl4, charge=-2, multiplicity=5)
    energies = [7111.62, 7112.48, 7113.05]
    fosc = [2.1e-5, 6.4e-5, 1.3e-5]
    w.add_excited_states(
        [{"index": i + 1, "excitation_energy": {"value": e, "units": "ev"},
          "multiplicity": 5, "oscillator_strength": f,
          "contributions": [{"from_mo": 0, "to_mo": 69 + i, "spin": "beta",
                             "weight": 0.97}]}
         for i, (e, f) in enumerate(zip(energies, fosc))],
        meta={"method": "TDDFT (TDA), Fe 1s donor window", "n_states": 3,
              "tda": True},
        section_id="core_states")
    w.add_xray_spectrum(
        "spectra.xas", energies, fosc, edge="K", element="Fe",
        atom_indices=[0], core_orbital="1s", intensity="oscillator_strength",
        y_units="dimensionless", energy_shift_ev=151.3,
        peaks=[{"x": e, "y": f, "label": f"1s -> state {i + 1}",
                "assignment": {"state_index": i + 1, "atom_index": 0}}
               for i, (e, f) in enumerate(zip(energies, fosc))],
        broadening={"type": "lorentzian", "fwhm": 1.25, "units": "ev"},
        states_ref="core_states", label="Fe K pre-edge")
    w.add_xray_spectrum(
        "spectra.xes", [7095.4, 7101.8, 7106.9, 7108.3],
        [0.08, 0.21, 0.64, 1.0], edge="K", element="Fe", atom_indices=[0],
        core_orbital="1s", line="Kβ2,5", intensity="normalized",
        broadening={"type": "gaussian", "fwhm": 2.5},
        label="Fe Kβ valence-to-core")
    _write("spectra_xas_tddft", w, [
        {"section": "spectra_xas", "member": "spectrum",
         "json_pointer": ["xray", "edge"], "equals": "K"},
        {"section": "spectra_xas", "member": "spectrum",
         "json_pointer": ["xray", "intensity"], "equals": "oscillator_strength"},
        {"section": "spectra_xas", "member": "spectrum",
         "json_pointer": ["xray", "energy_shift_ev"], "equals": 151.3,
         "tol": 1e-12},
        {"section": "spectra_xas", "member": "spectrum",
         "json_pointer": ["frequencies", 1], "equals": 7112.48, "tol": 1e-9},
        {"section": "spectra_xas", "member": "spectrum",
         "json_pointer": ["states_ref"], "equals": "core_states"},
        {"section": "spectra_xas", "member": "spectrum",
         "json_pointer": ["peaks", 2, "assignment", "state_index"], "equals": 3},
        {"section": "spectra_xes", "member": "spectrum",
         "json_pointer": ["xray", "line"], "equals": "Kβ2,5"},
        {"section": "spectra_xes", "member": "spectrum",
         "json_pointer": ["intensities", 3], "equals": 1.0, "tol": 1e-12},
    ], "Fe K pre-edge of tetrahedral [FeCl4]2- from three TDDFT core "
       "excitations (eV, oscillator strengths, 151.3 eV shift applied), each "
       "peak naming its state of the excited_states section through "
       "states_ref and its absorber through atom_index, and the Fe Kβ2,5 "
       "valence-to-core emission, normalized (spec § 5.4, X-ray absorption "
       "and emission).")

    # 2) Cl K-edge ΔSCF sticks of two absorbers; Fe L2,3 edge in hartree.
    w = QvfWriter(_PROGRAM, _VERSION, "FeCl4 Cl K-edge dSCF")
    w.add_structure(fecl4, charge=-2, multiplicity=5)
    cl = [2822.31, 2822.35, 2824.90]
    w.add_xray_spectrum(
        "spectra.xas", cl, [0.42, 0.40, 0.11], edge="K", element="Cl",
        atom_indices=[1, 2], core_orbital="1s", intensity="cross_section_mb",
        x_units="ev", energy_shift_ev=0.0,
        peaks=[{"x": cl[0], "label": "Cl1 1s -> Fe 3d (ΔSCF)",
                "assignment": {"atom_index": 1}},
               {"x": cl[1], "label": "Cl2 1s -> Fe 3d (ΔSCF)",
                "assignment": {"atom_index": 2}},
               {"x": cl[2], "label": "Cl1 1s -> Cl 4p (ΔSCF)",
                "assignment": {"atom_index": 1}}],
        broadening={"type": "gaussian", "fwhm": 0.6, "units": "ev"},
        label="Cl K-edge, ΔSCF", section_id="xas_cl_k")
    l23 = [26.0261, 26.0365, 26.4746]
    w.add_xray_spectrum(
        "spectra.xas", l23, [3.1, 4.4, 1.7], edge="L2,3", element="Fe",
        atom_indices=[0], core_orbital="2p", intensity="arbitrary",
        x_units="hartree", label="Fe L2,3 edge", section_id="xas_fe_l23")
    _write("spectra_xas_dscf", w, [
        {"section": "xas_cl_k", "member": "spectrum",
         "json_pointer": ["xray", "intensity"], "equals": "cross_section_mb"},
        {"section": "xas_cl_k", "member": "spectrum",
         "json_pointer": ["xray", "atom_indices", 1], "equals": 2},
        {"section": "xas_cl_k", "member": "spectrum",
         "json_pointer": ["peaks", 1, "assignment", "atom_index"], "equals": 2},
        {"section": "xas_cl_k", "member": "spectrum",
         "json_pointer": ["intensities", 0], "equals": 0.42, "tol": 1e-12},
        {"section": "xas_fe_l23", "member": "spectrum",
         "json_pointer": ["xray", "edge"], "equals": "L2,3"},
        {"section": "xas_fe_l23", "member": "spectrum",
         "json_pointer": ["x_units"], "equals": "hartree"},
        {"section": "xas_fe_l23", "member": "spectrum",
         "json_pointer": ["frequencies", 0], "equals": 26.0261, "tol": 1e-12},
    ], "Cl K-edge ΔSCF sticks of two chloride absorbers as cross sections in "
       "megabarn (no y_units), each peak naming its absorber, and an Fe L2,3 "
       "edge whose energies are declared in hartree, as a PySCF producer "
       "writes them without converting (spec § 5.4).")

    # 3) Mössbauer: two Fe sites, calibrated isomer shifts.
    alpha, beta, c = -0.31, 2.82, 11810.0
    def rho(delta):
        return c + (delta - beta) / alpha
    w = QvfWriter(_PROGRAM, _VERSION, "Fe(II)Fe(III) Mossbauer")
    w.add_structure([
        {"symbol": "Fe", "position": [0.0, 0.0, 0.0], "atomic_number": 26},
        {"symbol": "O", "position": [1.75, 0.95, 0.0], "atomic_number": 8},
        {"symbol": "Fe", "position": [3.5, 0.0, 0.0], "atomic_number": 26}])
    w.add_mossbauer([
        {"atom_index": 0, "symbol": "Fe", "label": "Fe(II) high spin",
         "isomer_shift_mm_s": 0.90, "quadrupole_splitting_mm_s": 2.80,
         "asymmetry": 0.25,
         "efg_au": [[-0.6, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.6]],
         "contact_density_au": rho(0.90), "linewidth_mm_s": 0.30,
         "population": 0.6},
        {"atom_index": 2, "symbol": "Fe", "label": "Fe(III) high spin",
         "isomer_shift_mm_s": 0.45, "quadrupole_splitting_mm_s": -0.72,
         "asymmetry": 0.6, "contact_density_au": rho(0.45),
         "hyperfine_field_t": -48.0, "linewidth_mm_s": 0.28,
         "population": 0.4}],
        isotope="57Fe", reference="alpha-Fe, 298 K", temperature_k=4.2,
        calibration={"alpha": alpha, "beta": beta, "c": c,
                     "label": "synthetic linear calibration"})
    sid = "spectra_mossbauer"
    _write("spectra_mossbauer", w, [
        {"section": sid, "member": "spectrum",
         "json_pointer": ["isotope"], "equals": "57Fe"},
        {"section": sid, "member": "spectrum",
         "json_pointer": ["sites", 0, "isomer_shift_mm_s"], "equals": 0.90,
         "tol": 1e-12},
        {"section": sid, "member": "spectrum",
         "json_pointer": ["sites", 1, "quadrupole_splitting_mm_s"],
         "equals": -0.72, "tol": 1e-12},
        {"section": sid, "member": "spectrum",
         "json_pointer": ["sites", 0, "efg_au", 2, 2], "equals": 1.6,
         "tol": 1e-12},
        {"section": sid, "member": "spectrum",
         "json_pointer": ["sites", 1, "hyperfine_field_t"], "equals": -48.0,
         "tol": 1e-12},
        {"section": sid, "member": "spectrum",
         "json_pointer": ["sites", 0, "contact_density_au"],
         "equals": rho(0.90), "tol": 1e-9},
        {"section": sid, "member": "spectrum",
         "json_pointer": ["calibration", "alpha"], "equals": alpha,
         "tol": 1e-12},
    ], "57Fe Mössbauer parameters of a mixed-valence Fe(II)/Fe(III) pair at "
       "4.2 K: isomer shifts against alpha-Fe and quadrupole splittings in "
       "mm/s, asymmetry parameters, an EFG tensor in atomic units, contact "
       "densities that the stated linear calibration maps onto the shifts, a "
       "hyperfine field in tesla, linewidths and populations (spec § 5.4, "
       "spectra.mossbauer).")

if __name__ == "__main__":
    build()
    build_zero_spin()
    build_unrestricted_ntos()
    build_biomolecule_identity()
    build_secondary_structure_ranges()
    build_reaction_energies()
    build_compound_units()
    build_xray_mossbauer()
