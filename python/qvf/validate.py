"""qvf_reader — a standalone reference reader / validator for QVF archives.

Depends only on the Python standard library and the generated ``qvf_kinds``
module shipped beside it (NumPy is optional and only used for convenience array
decoding). Mirrors the semantic checks that vibe-qc's
``validate_qvf`` performs, so an adopter can self-check archives without
installing vibe-qc. If ``jsonschema`` happens to be importable, the bundled
``spec/qvf_manifest.schema.json`` is additionally enforced for full structural
validation — both of ``manifest.json`` and of every JSON member payload that
has a payload definition (``_PAYLOAD_DEFS``: ``$defs/StructurePayload``,
``$defs/VolumeGridPayload``, …). The canonical kinds, the payload definitions,
the unit vocabulary and the finding codes come from ``qvf_kinds``, the tables
``scripts/gen_kinds.py`` generates from the schema.

Every finding carries a stable code from the specification's Appendix B
("Validation finding codes"), e.g. ``E_SHA256_MISMATCH`` or ``W_ORPHAN_MEMBER``.
The C++ reference validator emits the same codes.

CLI::

    python qvf_reader.py file.qvf          # prints a report; exits non-zero on error
    python qvf_reader.py --json file.qvf   # the same report as JSON (findings list)

Library::

    from qvf_reader import validate_qvf, QvfReader, Limits
    report = validate_qvf("file.qvf")
    if not report["ok"]:
        for f in report["findings"]:
            print(f["code"], f["where"], f["message"])

License: Apache-2.0.
"""

from __future__ import annotations

import dataclasses
import functools
import hashlib
from enum import IntEnum
import json
import math
import os
import re
import struct
import sys
import zipfile
from collections import Counter
from typing import Any, Optional

from . import _kinds as qvf_kinds


def _find_schema() -> str:
    """The schema is installed beside this module as package data."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "schema", "qvf_manifest.schema.json")


_SCHEMA_PATH = _find_schema()

# Canonical kinds (schema Section.oneOf branches), generated from the schema
# (qvf_kinds); a vendor kind (x_<...>) is always accepted.
_CANONICAL_KINDS = frozenset(qvf_kinds.CANONICAL_KINDS)

# JSON member payload definitions (spec § 5, § 6): a JSON member whose section
# kind and role appear here must conform to that ``$defs`` entry of the schema.
# Checked when jsonschema is installed. The table is generated
# (qvf_kinds.PAYLOAD_DEFS); scripts/gen_kinds.py is where a new payload
# definition is wired to its (kind, role).
_PAYLOAD_DEFS: dict[tuple[str, str], str] = dict(qvf_kinds.PAYLOAD_DEFS)

_DTYPE_ITEMSIZE = {
    "int8": 1, "int16": 2, "int32": 4, "int64": 8,
    "uint8": 1, "uint16": 2, "uint32": 4, "uint64": 8,
    "float32": 4, "float64": 8,
}

# Keys of the root `viewer_defaults` object that are not section ids (spec
# § 3.1, § 3.4).
_RESERVED_SECTION_IDS = frozenset({"auto_open", "bookmarks", "camera",
                                   "crossfade"})

# Hartree per unit, for the energy units the root consistency check converts
# (spec § 3.6; CODATA 2018). Keys are ASCII lower case; an energy in any other
# unit is not compared.
_HARTREE_PER_UNIT = {
    "hartree": 1.0, "eh": 1.0, "ha": 1.0,
    "ev": 1.0 / 27.211386245988,
    "kcal_per_mol": 1.0 / 627.5094740631,
    "kj_per_mol": 1.0 / 2625.4996394799,
}
#: Largest allowed difference, in hartree, between two root statements of one
#: energy (spec § 3.6).
_ENERGY_TOLERANCE_EH = 1e-8
#: Largest allowed difference between 1 and the sum of the k-point weights of a
#: ``wavefunction.bloch`` section (spec § 5.3.1).
_BLOCH_WEIGHT_TOLERANCE = 1e-6
_ASCII_LOWER = {c: c + 32 for c in range(ord("A"), ord("Z") + 1)}
# Generic references any section may carry, and the kind each must name
# (spec § 4.1).
_SECTION_REFS = (("structure_ref", "structure"), ("run_ref", "run.record"))
# Per-frame members of the trajectory family (spec § 5.5): the [n_frames,
# n_atoms, 3] roles both kinds carry, and the [n_frames] roles of each kind.
_TRAJECTORY_KINDS = ("trajectory", "reaction.path")
_TRAJECTORY_FRAME_VECTORS = ("velocities", "gradient")
_TRAJECTORY_FRAME_SCALARS = {
    "trajectory": ("times", "energies", "cell_volumes", "temperatures"),
    "reaction.path": ("times", "energies"),
}
#: Largest allowed difference, in hartree, between the binary ``energies``
#: member of a trajectory or reaction path and its ``metadata.energies``
#: (spec § 5.5).
_TRAJECTORY_ENERGY_TOLERANCE_EH = 1e-10

#: Bohr per Ångström (CODATA 2018, ten decimals), the factor ``f`` of spec
#: § 4.3, "Lattice conventions". It is ``qvf_writer.BOHR_PER_ANGSTROM``,
#: repeated here because this module needs no NumPy; a test keeps them equal.
_BOHR_PER_ANGSTROM = 1.8897261246
#: Largest relative difference between the span of a periodic grid axis and
#: the lattice vector it samples (spec § 5.2).
_GRID_SPAN_TOLERANCE = 1e-6
# The volume kinds whose `data` shape validators check (spec § 5.2): the rank,
# the extent of the trailing axis where the kind fixes one, and the form.
_VOLUME_DATA_SHAPES: dict[str, tuple[int, Optional[int], str]] = {
    "volume.laplacian": (3, None, "[nx, ny, nz]"),
    "volume.vector": (4, 3, "[nx, ny, nz, 3]"),
    "volume.complex": (4, 2, "[nx, ny, nz, 2]"),
    "volume.series": (4, None, "[n, nx, ny, nz]"),
}
# The unit name spec § 4.3 fixes for the grid `units.data` of a volume kind.
# One grid payload serves every `volume.*` kind, so the schema cannot fix it
# per kind; another unit name there is E_UNIT_CONTRADICTED.
_VOLUME_DATA_UNITS: dict[str, str] = {
    "volume.density": "e_per_bohr3",
    "volume.transition_density": "e_per_bohr3",
    "volume.potential": "hartree_per_e",
    "volume.rdg": "dimensionless",
    "volume.laplacian": "e_per_bohr5",
}
# The kinds a `volume.series` axis `trajectory_ref` may name (spec § 5.2).
_SERIES_REF_KINDS = ("trajectory", "reaction.path")

#: Stable validation finding codes and their severity (spec Appendix B),
#: generated from the appendix (qvf_kinds.FINDING_CODES).
FINDING_CODES: dict[str, str] = dict(qvf_kinds.FINDING_CODES)

# wavefunction.gto (spec § 5.3): the coefficient members, the optional AO-basis
# matrices, and the tolerances of the electron count and of the Appendix A.1
# normalization check.
_WF_COEFFICIENT_ROLES = ("mo_coefficients", "mo_coefficients_alpha",
                         "mo_coefficients_beta")
_WF_AO_MATRIX_ROLES = ("overlap", "density_matrix", "density_matrix_alpha",
                       "density_matrix_beta")
_N_ELECTRONS_TOLERANCE = 1e-6
_NORMALIZATION_RTOL = 1e-6


@dataclasses.dataclass(frozen=True)
class Limits:
    """Resource caps applied to untrusted archives (spec § 6).

    The defaults are the reference values the specification recommends. A
    validator that exceeds a cap reports an error and does not decompress or
    parse the offending member.
    """

    #: Per ZIP entry, uncompressed bytes (spec § 2.1; 1 gigavoxel of float64).
    max_member_bytes: int = 8 * 1024 ** 3
    #: Sum over all ZIP entries, uncompressed bytes.
    max_total_bytes: int = 64 * 1024 ** 3
    #: Optional application ratio policy; 0 disables it (the default).
    max_ratio: float = 0.0
    #: Entries smaller than this (uncompressed) are exempt from the ratio check.
    ratio_exempt_bytes: int = 1024 ** 2
    #: Largest allowed manifest.json or JSON member, uncompressed bytes.
    max_json_bytes: int = 256 * 1024 ** 2
    #: Deepest allowed nesting of JSON arrays/objects.
    max_json_depth: int = 64

    @classmethod
    def from_dict(cls, overrides: Optional[dict[str, Any]] = None) -> "Limits":
        """The default limits with ``overrides`` applied (unknown keys raise)."""
        return dataclasses.replace(cls(), **(overrides or {}))


# Guard against zip-bomb inputs: 8 GiB per member uncompressed (matches vibe-qc).
# Kept for backward compatibility; ``Limits`` is the source of truth.
_MAX_MEMBER_BYTES = Limits().max_member_bytes

_DRIVE_LETTER = re.compile(r"^[A-Za-z]:")
_JSON_SPECIAL = re.compile(rb'[\[\]{}"\\]')


def _path_problem(path: str) -> Optional[str]:
    """Why a member path / ZIP entry name is illegal (spec § 2.1), or None."""
    if "\\" in path:
        return "contains a backslash"
    if "\0" in path:
        return "contains a NUL character"
    if path.startswith("/"):
        return "begins with '/'"
    if _DRIVE_LETTER.match(path):
        return "begins with a drive letter"
    if ".." in path.split("/"):
        return "contains a '..' segment"
    return None


def _json_depth(raw: bytes, limit: int) -> int:
    """Nesting depth of the JSON arrays/objects in ``raw``, capped at limit + 1.

    A linear scan that visits only brackets, quotes and backslashes, skipping
    string contents (with escapes), so it never recurses and stops as soon as
    the limit is exceeded. Malformed JSON is left to the parser to report.
    """
    depth = deepest = 0
    in_string = False
    escaped_at = -1
    for match in _JSON_SPECIAL.finditer(raw):
        pos = match.start()
        c = raw[pos]
        if in_string:
            if pos == escaped_at:
                continue                    # the escaped character itself
            if c == 0x5C:                   # backslash: escapes the next byte
                escaped_at = pos + 1
            elif c == 0x22:                 # closing quote
                in_string = False
            continue
        if c == 0x22:
            in_string = True
        elif c == 0x5B or c == 0x7B:        # [ {
            depth += 1
            if depth > deepest:
                deepest = depth
                if deepest > limit:
                    return deepest
        elif c == 0x5D or c == 0x7D:        # ] }
            depth -= 1
    return deepest


def _pointer(*parts: Any) -> str:
    """A JSON Pointer (RFC 6901) into manifest.json."""
    return "".join(
        "/" + str(p).replace("~", "~0").replace("/", "~1") for p in parts)


#: The one ZIP entry name that may occur more than once: an archive updated in
#: place keeps its superseded manifests until it is repacked (spec § 2.3).
_MANIFEST = "manifest.json"


def _duplicate_entries(zf: zipfile.ZipFile) -> list[str]:
    """Entry names that occur more than once, other than ``manifest.json``."""
    counts = Counter(info.orig_filename for info in zf.infolist())
    return [name for name, n in counts.items() if n > 1 and name != _MANIFEST]


def _size_problem(member: dict[str, Any], info: zipfile.ZipInfo) -> Optional[tuple[str, str]]:
    """``(code, message)`` when a member spec's ``size`` is malformed or
    disagrees with the uncompressed size of its ZIP entry (spec § 4.2)."""
    if "size" not in member:
        return None
    size = member["size"]
    if not isinstance(size, int) or isinstance(size, bool) or size < 0:
        return ("E_MEMBER_MALFORMED",
                f"member {info.orig_filename!r}: size must be a non-negative "
                f"integer, got {size!r}")
    if size != info.file_size:
        return ("E_BINARY_SIZE",
                f"member {info.orig_filename!r}: declared size {size} != "
                f"{info.file_size} bytes recorded in the ZIP central directory")
    return None


class ZipUnreadable(ValueError, zipfile.BadZipFile):
    """An archive or entry that cannot be read (``E_ZIP_UNREADABLE``).

    A ``ValueError`` like every other refusal of the readers, and a
    ``zipfile.BadZipFile`` so that code catching the ``zipfile`` error keeps
    working.
    """


def _end_record_problem(zf: zipfile.ZipFile) -> Optional[str]:
    """Why the end of central directory record disagrees with the archive.

    ``zipfile`` reads every record that fits in the central directory and
    ignores the disk numbers and the entry count of the end record; the C++
    and JavaScript readers refuse a multi-volume archive and an entry count
    that differs from the records. Uses ``zipfile._EndRecData`` and answers
    ``None`` when it is absent.
    """
    end_record = getattr(zipfile, "_EndRecData", None)
    if end_record is None or zf.fp is None:
        return None
    endrec = end_record(zf.fp)
    if not endrec:
        return None
    if endrec[1] != 0 or endrec[2] != 0:
        return "multi-volume ZIP archives are not supported"
    count = len(zf.infolist())
    if endrec[3] != count or endrec[4] != count:
        return (f"the end record counts {endrec[4]} entries, the central "
                f"directory holds {count}")
    return None


def _open_zip(source: Any) -> zipfile.ZipFile:
    """``zipfile.ZipFile(source)``, refusing any damage as ``ZipUnreadable``."""
    try:
        zf = zipfile.ZipFile(source, "r")
    except MemoryError:
        raise
    except Exception as exc:  # noqa: BLE001 - any parse failure: unreadable
        raise ZipUnreadable(f"E_ZIP_UNREADABLE: not a readable zip: {exc}") from exc
    try:
        problem = _end_record_problem(zf)
    except MemoryError:
        zf.close()
        raise
    except Exception as exc:  # noqa: BLE001
        problem = str(exc)
    if problem:
        zf.close()
        raise ZipUnreadable(f"E_ZIP_UNREADABLE: not a readable zip: {problem}")
    return zf


def _inflates_past_end(source: Any, info: zipfile.ZipInfo) -> bool:
    """Whether a DEFLATE entry read up to its declared size has more output.

    ``zipfile`` stops at the declared size without looking further, so this
    finishes the stream on the reader's own decompressor (bounded: one byte of
    output at a time) and also reports a stream that ends before its final
    block. It relies on ``ZipExtFile`` internals and answers ``False`` when
    they are absent.
    """
    if info.compress_type != zipfile.ZIP_DEFLATED:
        return False
    d = getattr(source, "_decompressor", None)
    fileobj = getattr(source, "_fileobj", None)
    left = getattr(source, "_compress_left", None)
    if d is None or fileobj is None or not isinstance(left, int) \
            or not hasattr(d, "unconsumed_tail"):
        return False
    data = d.unconsumed_tail
    while not d.eof:
        if data:
            if d.decompress(data, 1):
                return True
            data = d.unconsumed_tail
            continue
        if left <= 0:
            return True                     # the stream has no final block
        data = fileobj.read(min(left, 1 << 16))
        if not data:
            return True
        left -= len(data)
    source._compress_left = left
    return False


def _iter_entry(zf: zipfile.ZipFile, info: zipfile.ZipInfo, chunk_size: int = 1 << 20):
    """Yield the uncompressed bytes of an entry in bounded chunks.

    Never inflates past the size the central directory records: a stream that
    is shorter, longer or damaged, a CRC-32 mismatch and an entry ``zipfile``
    cannot decode raise ``ZipUnreadable`` (spec § 2.1, § 6.1). The last chunk
    has been verified only when the iterator is exhausted.
    """
    label = f"ZIP entry {info.orig_filename!r}"
    if info.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
        # zipfile would also decode bzip2, LZMA and Zstandard; QVF allows
        # only STORE and DEFLATE (spec § 2.1), like the other readers.
        raise ZipUnreadable(f"E_ZIP_UNREADABLE: {label} uses compression method "
                            f"{info.compress_type}, not STORE (0) or DEFLATE (8)")
    if (info.compress_type == zipfile.ZIP_STORED
            and info.compress_size != info.file_size):
        raise ZipUnreadable(f"E_ZIP_UNREADABLE: {label} is stored with "
                            f"{info.compress_size} bytes but declares {info.file_size}")
    try:
        with zf.open(info) as source:
            left = info.file_size
            while left > 0:
                chunk = source.read(min(chunk_size, left))
                if not chunk:
                    raise ZipUnreadable(
                        f"E_ZIP_UNREADABLE: {label} inflates to "
                        f"{info.file_size - left} bytes, not its declared "
                        f"{info.file_size}")
                left -= len(chunk)
                yield chunk
            # One more read finishes the entry: zipfile checks the CRC-32 there.
            source.read(1)
            if _inflates_past_end(source, info):
                raise ZipUnreadable(
                    f"E_ZIP_UNREADABLE: {label} does not end at its declared "
                    f"{info.file_size} bytes")
    except (ZipUnreadable, MemoryError):
        raise
    except Exception as exc:  # noqa: BLE001 - any decoder failure: unreadable
        raise ZipUnreadable(f"E_ZIP_UNREADABLE: {label} cannot be read: {exc}") from exc


def _reason(exc: BaseException) -> str:
    """An exception's message without the code and label ``_iter_entry`` adds."""
    text = str(exc).removeprefix("E_ZIP_UNREADABLE: ")
    return text.split(" cannot be read: ", 1)[-1] if text.startswith("ZIP entry ") else text


def _read_entry(zf: zipfile.ZipFile, info: zipfile.ZipInfo) -> bytes:
    """The uncompressed bytes of an entry, read through ``_iter_entry``."""
    return b"".join(_iter_entry(zf, info))


def _reject_constant(name: str) -> Any:
    raise ValueError(f"{name} is not a JSON value (RFC 8259)")


def _parse_json(raw: bytes) -> Any:
    """Parse a JSON member or manifest: UTF-8 JSON as RFC 8259 defines it.

    Unlike ``json.loads`` on bytes, this does not detect UTF-16 or UTF-32
    and refuses the ``NaN``, ``Infinity`` and ``-Infinity`` extensions, as
    the C++ and JavaScript validators do. One leading UTF-8 byte order mark
    is skipped, as before. Raises ``ValueError`` (``UnicodeDecodeError`` for
    bytes that are not UTF-8).
    """
    return json.loads(raw.decode("utf-8-sig"), parse_constant=_reject_constant)


def _loads_coded(raw: bytes, code: str) -> Any:
    """``json.loads`` whose refusals carry a finding code.

    A ``json.JSONDecodeError`` keeps its type (its message gains the code);
    any other parser refusal, such as an integer too long to convert, becomes
    a ``ValueError`` with the code.
    """
    try:
        return _parse_json(raw)
    except json.JSONDecodeError as exc:
        raise json.JSONDecodeError(f"{code}: {exc.msg}", exc.doc, exc.pos) from exc
    except (ValueError, UnicodeError, RecursionError) as exc:
        raise ValueError(f"{code}: {exc}") from exc


class QvfReader:
    """Lazy reader for a ``.qvf`` archive. Verifies member checksums on read.

    Refuses (``ValueError``) an archive whose ZIP entry names are not unique
    (several ``manifest.json`` entries are allowed: the last one is the
    manifest, spec § 2.3), or whose ``byte_order`` is not ``"little"`` (spec
    § 2.2), and applies the ``limits`` JSON size and depth caps to
    ``manifest.json`` before parsing it and to JSON members before decoding
    them. A member whose ``size`` disagrees with its ZIP entry is refused when
    read (``E_BINARY_SIZE``).
    """

    def __init__(self, path: str, *, limits: Optional[Limits] = None) -> None:
        self.path = path
        self.limits = limits or Limits()
        self._zip = _open_zip(path)
        try:
            dups = _duplicate_entries(self._zip)
            if dups:
                raise ValueError(
                    f"E_ZIP_DUPLICATE_ENTRY: duplicate ZIP entry name(s) "
                    f"{sorted(dups)!r}; refusing an ambiguous archive")
            # getinfo() resolves a duplicated name to its last central-directory
            # record, which is the current manifest (spec § 2.3).
            raw = self._read_json_bytes(self._zip.getinfo(_MANIFEST), _MANIFEST)
            self.manifest: dict[str, Any] = _loads_coded(raw, "E_MANIFEST_JSON")
            order = (self.manifest.get("byte_order", "little")
                     if isinstance(self.manifest, dict) else "little")
            if order != "little":
                raise ValueError(
                    f"E_SCHEMA: byte_order {order!r} is not 'little', the only "
                    f"byte order of QVF v1 (spec § 2.2); refusing to guess")
        except BaseException:
            self._zip.close()
            raise

    def close(self) -> None:
        self._zip.close()

    def __enter__(self) -> "QvfReader":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    @property
    def sections(self) -> list[dict[str, Any]]:
        return self.manifest.get("sections", [])

    @property
    def typed_manifest(self) -> Any:
        """The manifest as the typed, lossless model of ``qvf_model``
        (``qvf_model.Manifest.from_json(self.manifest)``); TypeError when
        manifest.json is not a JSON object."""
        from . import model as qvf_model
        return qvf_model.Manifest.from_json(self.manifest)

    def _read_json_bytes(self, info: zipfile.ZipInfo, label: str) -> bytes:
        cap = self.limits.max_json_bytes
        if info.file_size > cap:
            raise ValueError(f"E_JSON_TOO_LARGE: {label} is {info.file_size} "
                             f"bytes, over the {cap}-byte JSON cap")
        raw = _read_entry(self._zip, info)
        limit = self.limits.max_json_depth
        if _json_depth(raw, limit) > limit:
            raise ValueError(f"E_JSON_TOO_DEEP: {label} nests deeper than "
                             f"{limit} levels")
        return raw

    def read_member(self, member: dict[str, Any], *, verify: bool = True) -> bytes:
        info = self._zip.getinfo(member["path"])
        if info.file_size > self.limits.max_member_bytes:
            raise ValueError(f"E_MEMBER_TOO_LARGE: {member['path']} is "
                             f"{info.file_size} bytes uncompressed")
        problem = _size_problem(member, info)
        if problem:
            raise ValueError(f"{problem[0]}: {problem[1]}")
        raw = _read_entry(self._zip, info)
        if verify and "sha256" in member:
            got = hashlib.sha256(raw).hexdigest()
            if got != member["sha256"]:
                raise ValueError(f"E_SHA256_MISMATCH: sha256 mismatch for {member['path']}")
        return raw

    def read_json(self, member: dict[str, Any]) -> Any:
        info = self._zip.getinfo(member["path"])
        if info.file_size > self.limits.max_json_bytes:
            raise ValueError(f"E_JSON_TOO_LARGE: {member['path']} is "
                             f"{info.file_size} bytes")
        raw = self.read_member(member)
        limit = self.limits.max_json_depth
        if _json_depth(raw, limit) > limit:
            raise ValueError(f"E_JSON_TOO_DEEP: {member['path']} nests deeper "
                             f"than {limit} levels")
        return _loads_coded(raw, "E_JSON_INVALID")

    def read_array(self, member: dict[str, Any]):
        import numpy as np
        raw = self.read_member(member)
        dt = np.dtype(member["dtype"]).newbyteorder("<")
        return np.frombuffer(raw, dtype=dt).reshape(member["shape"])


def _load_schema_validator():
    """Return a jsonschema validator for the bundled schema, or None."""
    try:
        import jsonschema  # type: ignore
    except Exception:
        return None
    try:
        with open(_SCHEMA_PATH, "r", encoding="utf-8") as fh:
            schema = json.load(fh)
        return jsonschema.Draft202012Validator(schema)
    except Exception:
        return None


def _load_payload_validator(defname: str):
    """Return a jsonschema validator for one ``$defs`` entry, or None.

    Payload definitions describe archive members rather than manifest content, so
    they are unreachable from the schema root. Re-rooting a document at the
    definition — carrying ``$defs`` along so its internal ``$ref``s still
    resolve — is what makes them usable.
    """
    try:
        import jsonschema  # type: ignore
    except Exception:
        return None
    try:
        with open(_SCHEMA_PATH, "r", encoding="utf-8") as fh:
            schema = json.load(fh)
        return jsonschema.Draft202012Validator({
            "$schema": schema["$schema"],
            "$ref": f"#/$defs/{defname}",
            "$defs": schema["$defs"],
        })
    except Exception:
        return None


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _structure_payload_errors(sid: str, payload: Any) -> list[str]:
    """Semantic checks on a ``structure`` payload (QVF spec § 5.1).

    Enforced without ``jsonschema`` so the periodicity contract holds for every
    validator, not only those with the optional dependency installed. ``pbc`` is
    the normative carrier; ``dimensionality`` is derived from it.
    """
    if not isinstance(payload, dict):
        return [f"section {sid!r}: structure payload must be a JSON object"]

    errors: list[str] = []
    pbc = payload.get("pbc")
    lattice = payload.get("lattice_vectors")
    dim = payload.get("dimensionality")

    if dim is not None and not _is_number(dim):
        errors.append(f"section {sid!r}: dimensionality must be an integer, "
                      f"got {dim!r} (spec § 5.1)")
        dim = None

    if pbc is None:
        if lattice is not None:
            errors.append(
                f"section {sid!r}: structure payload has lattice_vectors but no "
                f"'pbc'; pbc is REQUIRED for periodic structures (spec § 5.1)")
        if dim not in (None, 0):
            errors.append(
                f"section {sid!r}: dimensionality {dim!r} without 'pbc'; absent "
                f"pbc means a molecule, so dimensionality must be 0 (spec § 5.1)")
        return errors

    if (not isinstance(pbc, list) or len(pbc) != 3
            or not all(isinstance(x, bool) for x in pbc)):
        return errors + [f"section {sid!r}: pbc must be three booleans, "
                         f"got {pbc!r}"]

    n_periodic = sum(pbc)
    if n_periodic and not isinstance(lattice, list):
        errors.append(
            f"section {sid!r}: pbc marks {n_periodic} periodic axis/axes but "
            f"lattice_vectors is {'null' if lattice is None else 'absent'}")
    if dim is not None and dim != n_periodic:
        errors.append(
            f"section {sid!r}: dimensionality {dim} != sum(pbc) {n_periodic} "
            f"(spec § 5.1 invariant)")
    return errors


def _spectrum_payload_errors(sid: str, payload: Any) -> list[str]:
    """The equal-length rule of a 1-D ``spectrum`` payload (QVF spec § 5.4).

    ``frequencies``, ``intensities`` and, when present, ``x`` must have one
    shared length. JSON Schema cannot say so, so this runs with or without
    ``jsonschema``; types and required keys are left to the schema check.
    """
    if not isinstance(payload, dict):
        return []
    freqs = payload.get("frequencies")
    if not isinstance(freqs, list):
        return []
    errors = []
    for key in ("intensities", "x"):
        other = payload.get(key)
        if isinstance(other, list) and len(other) != len(freqs):
            errors.append(
                f"section {sid!r}: spectrum {key} has {len(other)} entries but "
                f"frequencies has {len(freqs)}; they must have the same "
                f"length (spec § 5.4)")
    return errors


@functools.lru_cache(maxsize=1)
def _unit_vocabulary() -> Optional[tuple[frozenset, frozenset]]:
    """The unit names and the legacy spellings of spec § 4.3
    (``$defs/UnitName``, ``$defs/LegacyUnitName``), from the tables generated
    from the schema (``qvf_kinds``)."""
    return (frozenset(qvf_kinds.UNITS), frozenset(qvf_kinds.LEGACY_UNITS))


def _declared_units(kind: str, payload: Any,
                    role: Optional[str] = None) -> list[tuple[str, Any]]:
    """``(location, value)`` of every declared unit of a member payload (spec
    § 4.3, "Declared units"): each value of its ``units`` object, the
    ``units`` string of a ``volume.series`` ``axis`` (``role``) and, for a
    one-dimensional spectrum, ``x_units``, ``y_units`` and
    ``broadening.units``, and for ``excited_states`` each value of
    ``meta.units``. The location is a path inside the payload."""
    if not isinstance(payload, dict):
        return []
    out: list[tuple[str, Any]] = []
    units = payload.get("units")
    if isinstance(units, dict):
        out.extend((f"units/{key}", value) for key, value in units.items())
    elif ("units" in payload
          and _PAYLOAD_DEFS.get((kind, role)) == "VolumeSeriesAxisPayload"):
        out.append(("units", units))
    if kind == "excited_states":
        meta = payload.get("meta")
        units = meta.get("units") if isinstance(meta, dict) else None
        if isinstance(units, dict):
            out.extend((f"meta/units/{key}", value)
                       for key, value in units.items())
    if _PAYLOAD_DEFS.get((kind, "spectrum")) == "SpectrumPayload":
        out.extend((key, payload[key]) for key in ("x_units", "y_units")
                   if key in payload)
        broadening = payload.get("broadening")
        if isinstance(broadening, dict) and "units" in broadening:
            out.append(("broadening/units", broadening["units"]))
    return out

def _ascii_lower(text: str) -> str:
    """``text`` with A-Z lowered and everything else unchanged (spec § 3.6)."""
    return text.translate(_ASCII_LOWER)


def _energy_eh(energy: Any) -> Optional[float]:
    """A ``{value, units}`` energy in hartree, or None when it is malformed or
    its unit is not one the root consistency check converts (spec § 3.6)."""
    if not isinstance(energy, dict):
        return None
    value, units = energy.get("value"), energy.get("units")
    if not _is_number(value) or not isinstance(units, str):
        return None
    factor = _HARTREE_PER_UNIT.get(_ascii_lower(units))
    return None if factor is None else float(value) * factor


def _root_consistency_errors(manifest: dict) -> list[tuple[str, str]]:
    """Root blocks that state one energy twice must agree (spec § 3.6).

    Compares every ``energies`` total stated for the method that
    ``provenance.method`` names (``total`` when ``energies.method`` is that
    method, and each such ``by_method`` entry) with ``provenance.scf_energy``,
    and ``energies.components.zpve`` with ``thermochemistry.zpve_eh``. Method
    names match case-insensitively in ASCII. Returns ``(where, message)``
    pairs; a value that is malformed or in an unrecognised unit is skipped.
    """
    energies = manifest.get("energies")
    if not isinstance(energies, dict):
        return []
    errors: list[tuple[str, str]] = []

    def compare(where: str, label: str, got: Optional[float], ref: float,
                ref_label: str) -> None:
        if got is not None and abs(got - ref) > _ENERGY_TOLERANCE_EH:
            errors.append((where, f"{label} is {got:.10f} hartree but "
                                  f"{ref_label} is {ref:.10f} hartree; they "
                                  f"state one energy and must agree within "
                                  f"1e-8 hartree (spec § 3.6)"))

    prov = manifest.get("provenance")
    method = prov.get("method") if isinstance(prov, dict) else None
    scf = _energy_eh(prov.get("scf_energy")) if isinstance(prov, dict) else None
    if isinstance(method, str) and scf is not None:
        key = _ascii_lower(method)
        ref_label = f"provenance.scf_energy (method {method!r})"
        emethod = energies.get("method")
        if isinstance(emethod, str) and _ascii_lower(emethod) == key:
            compare("/energies/total", "energies.total",
                    _energy_eh(energies.get("total")), scf, ref_label)
        by_method = energies.get("by_method")
        for k, entry in enumerate(by_method if isinstance(by_method, list) else []):
            if (isinstance(entry, dict) and isinstance(entry.get("method"), str)
                    and _ascii_lower(entry["method"]) == key):
                compare(_pointer("energies", "by_method", k, "total"),
                        f"energies.by_method[{k}].total",
                        _energy_eh(entry.get("total")), scf, ref_label)
    thermo = manifest.get("thermochemistry")
    components = energies.get("components")
    if (isinstance(thermo, dict) and _is_number(thermo.get("zpve_eh"))
            and isinstance(components, dict)):
        compare("/energies/components/zpve", "energies.components.zpve",
                _energy_eh(components.get("zpve")), float(thermo["zpve_eh"]),
                "thermochemistry.zpve_eh")
    return errors


def _shape_of(member: Any) -> Optional[list[int]]:
    """A member spec's ``shape`` if it is a list of non-negative integers."""
    shape = member.get("shape") if isinstance(member, dict) else None
    if not isinstance(shape, list) or not all(
            isinstance(d, int) and not isinstance(d, bool) and d >= 0
            for d in shape):
        return None
    return shape


def _bands_shape_errors(sid: str, members: Any) -> list[tuple[str, str]]:
    """Spin-axis agreement of a ``bands`` section's binary members (spec § 5.6).

    ``eigenvalues`` is ``[n_spin, n_k, n_bands]``. A rank-4 ``projections``
    (``[n_spin, n_k, n_bands, n_channels]``) and ``occupations``
    (``[n_spin, n_k, n_bands]``) must repeat those three extents; the legacy
    rank-3 ``projections`` (``[n_k, n_bands, n_channels]``) must repeat n_k
    and n_bands. Shapes of any other rank are left to the schema. Returns
    ``(role, message)`` pairs, checked from the manifest alone.
    """
    if not isinstance(members, dict):
        return []
    eig = _shape_of(members.get("eigenvalues"))
    if eig is None or len(eig) != 3:
        return []
    errors: list[tuple[str, str]] = []
    proj = _shape_of(members.get("projections"))
    if proj is not None and len(proj) == 4 and proj[:3] != eig:
        errors.append(("projections",
                       f"section {sid!r}: projections shape {proj} has "
                       f"[n_spin, n_k, n_bands] = {proj[:3]}, but eigenvalues "
                       f"is {eig} (spec § 5.6)"))
    elif proj is not None and len(proj) == 3 and proj[:2] != eig[1:]:
        errors.append(("projections",
                       f"section {sid!r}: rank-3 projections shape {proj} has "
                       f"[n_k, n_bands] = {proj[:2]}, but eigenvalues has "
                       f"{eig[1:]} (spec § 5.6)"))
    occ = _shape_of(members.get("occupations"))
    if occ is not None and len(occ) == 3 and occ != eig:
        errors.append(("occupations",
                       f"section {sid!r}: occupations shape {occ} differs "
                       f"from eigenvalues {eig}; both are "
                       f"[n_spin, n_k, n_bands] (spec § 5.6)"))
    return errors


def _float64_values(raw: Optional[bytes], member: Any) -> Optional[list[float]]:
    """The values of a rank-1 float64 binary member whose byte length matches
    its shape, decoded little-endian; None for anything else."""
    shape = _shape_of(member)
    if (raw is None or not isinstance(member, dict)
            or member.get("dtype") != "float64" or shape is None
            or len(shape) != 1 or len(raw) != 8 * shape[0]):
        return None
    return list(struct.unpack(f"<{shape[0]}d", raw))


def _trajectory_errors(sid: str, kind: str, members: Any, metadata: Any,
                       energies_raw: Optional[bytes]) -> list[tuple[str, str]]:
    """Shape and consistency rules of a ``trajectory`` or ``reaction.path``
    section (spec § 5.5).

    ``coords`` is float32 only when ``metadata.precision`` is ``"float32"``
    (and ``precision``, when present, names the dtype). Against a ``coords``
    of shape ``[n_frames, n_atoms, 3]`` (another rank is left to the schema),
    ``velocities`` and ``gradient`` are ``[n_frames, n_atoms, 3]``; ``times`` and ``energies`` (and, for a trajectory,
    ``cell_volumes`` and ``temperatures``) are ``[n_frames]``; a trajectory
    ``lattice`` is ``[3, 3]`` or ``[n_frames, 3, 3]``, and so is the
    ``lattice`` of a reaction path, each read in its own convention (§ 4.3)
    but checked here by shape only. The binary ``energies`` agree with
    ``metadata.energies`` within 1e-10 hartree. A trajectory with a
    ``lattice`` member carries ``metadata.pbc``, and a true ``pbc`` entry
    needs the ``lattice`` member. Checks that need the metadata are skipped
    when it could not be read (its own finding says why). Returns ``(role,
    message)`` pairs, located at the member of that role.
    """
    if not isinstance(members, dict):
        return []
    errors: list[tuple[str, str]] = []
    sec = f"section {sid!r}: "
    meta = metadata if isinstance(metadata, dict) else None
    coords_spec = members.get("coords")
    coords = _shape_of(coords_spec)
    dtype = coords_spec.get("dtype") if isinstance(coords_spec, dict) else None
    if meta is not None and dtype in ("float32", "float64"):
        precision = meta.get("precision", "float64")
        if dtype == "float32" and precision != "float32":
            errors.append(("coords",
                           sec + f"coords are float32 but metadata.precision is "
                                 f"{precision!r}; float32 coords require "
                                 f"precision 'float32' (spec § 5.5)"))
        elif dtype == "float64" and precision == "float32":
            errors.append(("coords",
                           sec + "metadata.precision is 'float32' but coords are "
                                 "float64; precision must name the dtype of "
                                 "coords (spec § 5.5)"))
    # A coords shape of another rank is left to the schema; the extents below
    # are only compared against a [n_frames, n_atoms, 3] coords.
    if coords is not None and len(coords) == 3 and coords[2] == 3:
        n_frames, n_atoms = coords[0], coords[1]
        for role in _TRAJECTORY_FRAME_VECTORS:
            shape = _shape_of(members.get(role))
            if shape is not None and shape != [n_frames, n_atoms, 3]:
                errors.append((role, sec + f"{role} shape {shape} is not "
                                           f"[n_frames, n_atoms, 3] = "
                                           f"{[n_frames, n_atoms, 3]}, the "
                                           f"extents of coords (spec § 5.5)"))
        for role in _TRAJECTORY_FRAME_SCALARS[kind]:
            shape = _shape_of(members.get(role))
            if shape is not None and shape != [n_frames]:
                errors.append((role, sec + f"{role} shape {shape} is not "
                                           f"[n_frames] = {[n_frames]}, the "
                                           f"frames of coords (spec § 5.5)"))
        lattice = _shape_of(members.get("lattice"))
        if lattice is not None and lattice not in ([3, 3], [n_frames, 3, 3]):
            errors.append(("lattice", sec + f"lattice shape {lattice} is neither "
                                            f"[3, 3] nor [n_frames, 3, 3] = "
                                            f"{[n_frames, 3, 3]} (spec § 5.5)"))
    if meta is not None:
        binary = _float64_values(energies_raw, members.get("energies"))
        listed = meta.get("energies")
        if (binary is not None and isinstance(listed, list)
                and all(_is_number(x) for x in listed)):
            if len(listed) != len(binary):
                errors.append(("energies",
                               sec + f"the energies member has {len(binary)} "
                                     f"values but metadata.energies has "
                                     f"{len(listed)}; they must agree "
                                     f"(spec § 5.5)"))
            else:
                bad = [k for k, (a, b) in enumerate(zip(binary, listed))
                       if not abs(a - float(b)) <= _TRAJECTORY_ENERGY_TOLERANCE_EH]
                if bad:
                    k = bad[0]
                    errors.append((
                        "energies",
                        sec + f"the energies member disagrees with "
                              f"metadata.energies at {len(bad)} frame(s), first "
                              f"frame {k}: {binary[k]!r} vs {float(listed[k])!r} "
                              f"hartree; they must agree within 1e-10 hartree "
                              f"(spec § 5.5)"))
    if kind == "trajectory" and meta is not None:
        pbc = meta.get("pbc")
        if "lattice" in members and pbc is None:
            errors.append(("metadata",
                           sec + "the section has a lattice member but its "
                                 "metadata has no pbc; pbc is REQUIRED when "
                                 "lattice is present (spec § 5.5)"))
        elif (isinstance(pbc, list) and any(x is True for x in pbc)
              and "lattice" not in members):
            errors.append(("metadata",
                           sec + f"metadata.pbc {pbc} marks a periodic axis but "
                                 f"the section has no lattice member "
                                 f"(spec § 5.5)"))
    return errors


def _derivative_errors(sid: str, kind: str, section: dict, members: Any,
                       meta: Any, n_atoms_of: dict[str, int],
                       kinds_by_id: dict[str, str],
                       only_structure: Optional[str]
                       ) -> list[tuple[str, str, str]]:
    """Shape and reference rules of a ``forces`` or ``hessian`` section (spec
    § 5.10, § 5.11).

    ``gradient`` is ``[n_atoms, 3]``, ``stress`` ``[3, 3]``, ``matrix``
    square with an extent 3N divisible by 3, and ``meta.atom_masses_amu`` has
    N entries. The section's structure is the section-level ``structure_ref``
    when the section carries one, else ``meta.structure_ref``, else the
    archive's only ``structure`` section; when it resolves to a structure
    whose atom count is known (``n_atoms_of``), n_atoms and N must equal it.
    A ``meta.structure_ref`` must name a ``structure`` section. Returns
    ``(code, role, message)`` triples, located at the member of that role.
    """
    if not isinstance(members, dict):
        return []
    errors: list[tuple[str, str, str]] = []
    sec = f"section {sid!r}: "
    meta = meta if isinstance(meta, dict) else {}
    meta_ref = meta.get("structure_ref")
    if "structure_ref" in meta and (not isinstance(meta_ref, str)
                                    or kinds_by_id.get(meta_ref) != "structure"):
        errors.append(("E_REF_UNRESOLVED", "meta",
                       sec + f"meta.structure_ref {meta_ref!r} does not name a "
                             f"structure section (spec § 5.10)"))
    if "structure_ref" in section:
        ref = section["structure_ref"]
        structure = ref if isinstance(ref, str) else None
    elif isinstance(meta_ref, str):
        structure = meta_ref
    else:
        structure = only_structure
    n_struct = n_atoms_of.get(structure) if structure is not None else None

    if kind == "forces":
        grad = _shape_of(members.get("gradient"))
        if grad is not None and (len(grad) != 2 or grad[1] != 3):
            errors.append(("E_MEMBER_MALFORMED", "gradient",
                           sec + f"gradient shape {grad} is not [n_atoms, 3] "
                                 f"(spec § 5.10)"))
        elif grad is not None and n_struct is not None and grad[0] != n_struct:
            errors.append(("E_MEMBER_MALFORMED", "gradient",
                           sec + f"gradient has {grad[0]} rows but structure "
                                 f"{structure!r} has {n_struct} atoms "
                                 f"(spec § 5.10)"))
        stress = _shape_of(members.get("stress"))
        if stress is not None and stress != [3, 3]:
            errors.append(("E_MEMBER_MALFORMED", "stress",
                           sec + f"stress shape {stress} is not [3, 3] "
                                 f"(spec § 5.10)"))
        return errors

    mat = _shape_of(members.get("matrix"))
    if mat is None:
        return errors
    if len(mat) != 2 or mat[0] != mat[1]:
        errors.append(("E_MEMBER_MALFORMED", "matrix",
                       sec + f"matrix shape {mat} is not square (spec § 5.11)"))
        return errors
    if mat[0] % 3:
        errors.append(("E_MEMBER_MALFORMED", "matrix",
                       sec + f"matrix extent {mat[0]} is not divisible by 3, so "
                             f"the matrix is not [3N, 3N] (spec § 5.11)"))
        return errors
    n = mat[0] // 3
    if n_struct is not None and n != n_struct:
        errors.append(("E_MEMBER_MALFORMED", "matrix",
                       sec + f"matrix {mat} is for {n} atoms but structure "
                             f"{structure!r} has {n_struct} atoms (spec § 5.11)"))
    masses = meta.get("atom_masses_amu")
    if isinstance(masses, list) and len(masses) != n:
        errors.append(("E_MEMBER_MALFORMED", "meta",
                       sec + f"meta.atom_masses_amu has {len(masses)} entries "
                             f"but the matrix is for {n} atoms (spec § 5.11)"))
    return errors

def _bonds_errors(sid: str, section: dict, payload: Any,
                  kinds_by_id: dict[str, str],
                  n_atoms_of: dict[str, int]) -> list[str]:
    """Atom indices of a standalone ``bonds`` section against the structure
    its ``structure_ref`` names (spec § 5.1).

    Only a section whose ``structure_ref`` names a ``structure`` section with
    a known atom count (``n_atoms_of``) is checked; a reference that does not
    resolve is the generic reference check's finding, and a section without
    one is not checked. Every integer ``i`` and ``j`` of the payload's
    ``pairs`` must be an atom of that structure. Returns messages, located at
    the ``bonds`` member.
    """
    ref = section.get("structure_ref")
    if not isinstance(ref, str) or kinds_by_id.get(ref) != "structure":
        return []
    n = n_atoms_of.get(ref)
    pairs = payload.get("pairs") if isinstance(payload, dict) else None
    if n is None or not isinstance(pairs, list):
        return []
    errors: list[str] = []
    for k, pair in enumerate(pairs):
        if not isinstance(pair, dict):
            continue
        for key in ("i", "j"):
            atom = pair.get(key)
            if _is_int(atom) and not 0 <= atom < n:
                errors.append(f"section {sid!r}: pairs[{k}].{key} {atom} is "
                              f"not one of the {n} atoms of structure "
                              f"{ref!r}, which structure_ref names "
                              f"(spec § 5.1)")
    return errors


# --- excited states, NTOs and transition densities (spec § 5.12) -------------

#: Kinds a state's ``density_ref`` may name (spec § 5.12).
_STATE_DENSITY_KINDS = ("volume.transition_density", "volume.difference")


def _is_int(x: Any) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


def _mo_rows(members: Any) -> dict[str, int]:
    """Row counts of a ``wavefunction.gto`` section's coefficient members,
    read from their manifest shapes: key ``""`` for ``mo_coefficients``,
    ``"alpha"`` and ``"beta"`` for the unrestricted pair. A member without a
    usable shape is left out."""
    rows: dict[str, int] = {}
    if not isinstance(members, dict):
        return rows
    for role, key in (("mo_coefficients", ""), ("mo_coefficients_alpha", "alpha"),
                      ("mo_coefficients_beta", "beta")):
        shape = _shape_of(members.get(role))
        if shape:
            rows[key] = shape[0]
    return rows


def _contribution_rows(rows: dict[str, int], spin: Any) -> Optional[int]:
    """The number of orbital rows a contribution with ``spin`` may index:
    those of ``mo_coefficients``, else of the channel ``spin`` selects, else
    (no ``spin``) of both channels. None when unknown."""
    if "" in rows:
        return rows[""]
    if spin in ("alpha", "beta"):
        return rows.get(spin)
    known = [rows[c] for c in ("alpha", "beta") if c in rows]
    return min(known) if known else None


def _state_indices(payload: Any) -> Optional[set[int]]:
    """The ``index`` values of an ``excited_states`` payload's states, or None
    when the payload has no ``states`` array."""
    if not isinstance(payload, dict) or not isinstance(payload.get("states"), list):
        return None
    return {s["index"] for s in payload["states"]
            if isinstance(s, dict) and _is_int(s.get("index"))}


def _excited_states_errors(sid: str, section: dict, payload: Any,
                           kinds_by_id: dict[str, str],
                           mo_rows_of: dict[str, dict[str, int]]
                           ) -> list[tuple[str, str, str]]:
    """Rules of an ``excited_states`` section that JSON Schema cannot express
    (spec § 5.12).

    The ``index`` values of ``states`` are unique and ascending, and
    ``meta.n_states`` equals the number of states. The section's
    ``wavefunction_ref`` and every state's ``nto_ref`` name
    ``wavefunction.gto`` sections, and every ``density_ref`` a
    ``volume.transition_density`` or ``volume.difference`` section. When
    ``wavefunction_ref`` names a ``wavefunction.gto`` section, every
    contribution's ``from_mo`` and ``to_mo`` is a row of its coefficient
    member (``mo_rows_of``, from the manifest shapes). Returns ``(code,
    field, message)`` triples; ``field`` is the section field a finding is
    about, or ``""`` for the ``states`` member.
    """
    errors: list[tuple[str, str, str]] = []
    sec = f"section {sid!r}: "
    wf_ref: Optional[str] = None
    if "wavefunction_ref" in section:
        ref = section["wavefunction_ref"]
        if not isinstance(ref, str) or kinds_by_id.get(ref) != "wavefunction.gto":
            errors.append(("E_REF_UNRESOLVED", "wavefunction_ref",
                           sec + f"wavefunction_ref {ref!r} does not name a "
                                 f"wavefunction.gto section (spec § 5.12)"))
        else:
            wf_ref = ref
    states = payload.get("states") if isinstance(payload, dict) else None
    if not isinstance(states, list):
        return errors
    rows = mo_rows_of.get(wf_ref, {}) if wf_ref is not None else {}
    seen: set[int] = set()
    prev: Optional[int] = None
    for k, state in enumerate(states):
        if not isinstance(state, dict):
            continue
        idx = state.get("index")
        if _is_int(idx):
            if idx in seen:
                errors.append(("E_MEMBER_MALFORMED", "",
                               sec + f"states[{k}].index {idx} repeats an "
                                     f"earlier index; state indices are unique "
                                     f"(spec § 5.12)"))
            elif prev is not None and idx < prev:
                errors.append(("E_MEMBER_MALFORMED", "",
                               sec + f"states[{k}].index {idx} follows index "
                                     f"{prev}; states are listed in ascending "
                                     f"order of index (spec § 5.12)"))
            seen.add(idx)
            prev = idx
        for field, targets, what in (
                ("nto_ref", ("wavefunction.gto",), "a wavefunction.gto section"),
                ("density_ref", _STATE_DENSITY_KINDS,
                 "a volume.transition_density or volume.difference section")):
            if field in state:
                ref = state[field]
                if not isinstance(ref, str) or kinds_by_id.get(ref) not in targets:
                    errors.append(("E_REF_UNRESOLVED", "",
                                   sec + f"states[{k}].{field} {ref!r} does not "
                                         f"name {what} (spec § 5.12)"))
        contributions = state.get("contributions")
        if wf_ref is None or not isinstance(contributions, list):
            continue
        for c, contrib in enumerate(contributions):
            if not isinstance(contrib, dict):
                continue
            n = _contribution_rows(rows, contrib.get("spin"))
            if n is None:
                continue
            for key in ("from_mo", "to_mo"):
                mo = contrib.get(key)
                if _is_int(mo) and not 0 <= mo < n:
                    errors.append(("E_MEMBER_MALFORMED", "",
                                   sec + f"states[{k}].contributions[{c}].{key} "
                                         f"{mo} is not one of the {n} orbital "
                                         f"rows of wavefunction {wf_ref!r} "
                                         f"(spec § 5.12)"))
    meta = payload.get("meta")
    if (isinstance(meta, dict) and _is_int(meta.get("n_states"))
            and meta["n_states"] != len(states)):
        errors.append(("E_MEMBER_MALFORMED", "",
                       sec + f"meta.n_states is {meta['n_states']} but states "
                             f"has {len(states)} entries (spec § 5.12)"))
    return errors


def _spectrum_states_errors(sid: str, payload: Any, kinds_by_id: dict[str, str],
                            indices_of: dict[str, Optional[set[int]]]
                            ) -> list[str]:
    """A 1-D spectrum's ``states_ref`` names an ``excited_states`` section,
    and every ``peaks[].assignment.state_index`` is the ``index`` of one of
    its states (spec § 5.4). Returns messages."""
    if not isinstance(payload, dict) or "states_ref" not in payload:
        return []
    ref = payload["states_ref"]
    sec = f"section {sid!r}: "
    if not isinstance(ref, str) or kinds_by_id.get(ref) != "excited_states":
        return [sec + f"spectrum states_ref {ref!r} does not name an "
                      f"excited_states section (spec § 5.4)"]
    indices = indices_of.get(ref)
    peaks = payload.get("peaks")
    if indices is None or not isinstance(peaks, list):
        return []
    errors = []
    for k, peak in enumerate(peaks):
        assignment = peak.get("assignment") if isinstance(peak, dict) else None
        state = (assignment.get("state_index")
                 if isinstance(assignment, dict) else None)
        if _is_int(state) and state not in indices:
            errors.append(sec + f"peaks[{k}].assignment.state_index {state} "
                                f"names no state of excited_states section "
                                f"{ref!r} (spec § 5.4)")
    return errors


def _transition_density_errors(sid: str, section: dict,
                               kinds_by_id: dict[str, str],
                               indices_of: dict[str, Optional[set[int]]]
                               ) -> list[tuple[str, str]]:
    """A ``volume.transition_density`` section's ``excited_states_ref`` names
    an ``excited_states`` section, and its ``state_index`` the ``index`` of
    one of its states (spec § 5.2). Returns ``(field, message)`` pairs."""
    if "excited_states_ref" not in section:
        return []
    ref = section["excited_states_ref"]
    sec = f"section {sid!r}: "
    if not isinstance(ref, str) or kinds_by_id.get(ref) != "excited_states":
        return [("excited_states_ref",
                 sec + f"excited_states_ref {ref!r} does not name an "
                       f"excited_states section (spec § 5.2)")]
    state = section.get("state_index")
    indices = indices_of.get(ref)
    if _is_int(state) and indices is not None and state not in indices:
        return [("state_index",
                 sec + f"state_index {state} names no state of excited_states "
                       f"section {ref!r} (spec § 5.2)")]
    return []


#: Largest allowed disagreement between two statements of one property
#: (spec § 5.13): the dipole of ``properties.electric`` against the root
#: ``dipole_moment`` (debye), and ``static`` against ``electronic + ionic``.
_PROPERTY_TOLERANCE = 1e-6
_PROPERTY_KINDS = ("properties.electric", "properties.elastic",
                   "properties.dielectric")


def _matrix_of(value: Any, rows: int, cols: int) -> Optional[list[list[float]]]:
    """``value`` as a ``rows`` x ``cols`` matrix of finite numbers, or None
    when it is anything else (the payload schema reports the malformation)."""
    if not isinstance(value, list) or len(value) != rows:
        return None
    out = []
    for row in value:
        if (not isinstance(row, list) or len(row) != cols
                or not all(_is_number(x) and math.isfinite(x) for x in row)):
            return None
        out.append([float(x) for x in row])
    return out


def _asymmetry(matrix: list[list[float]]) -> Optional[tuple[int, int]]:
    """The first ``(i, j)``, i < j, whose two entries differ by more than
    1e-6 of the largest absolute entry of the matrix; None when symmetric."""
    scale = max(abs(x) for row in matrix for x in row)
    n = len(matrix)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(matrix[i][j] - matrix[j][i]) > _PROPERTY_TOLERANCE * scale:
                return i, j
    return None


def _property_errors(sid: str, kind: str, section: dict, payload: Any,
                     manifest: dict, n_atoms_of: dict[str, int],
                     only_structure: Optional[str]) -> list[tuple[str, str]]:
    """Semantic rules of the property kinds that JSON Schema cannot express
    (spec § 5.13). Returns ``(code, message)`` pairs, located at the section's
    ``properties`` member.

    ``properties.elastic``: ``stiffness`` and ``compliance`` are symmetric to
    within 1e-6 of their largest entry (``E_MEMBER_MALFORMED``).
    ``properties.dielectric``: ``static`` equals ``electronic + ionic`` within
    1e-6 when all three are present, and ``born_charges`` has one tensor per
    atom of the section's structure (the section-level ``structure_ref``, else
    the archive's only ``structure``) when that resolves
    (``E_MEMBER_MALFORMED``). ``properties.electric``: a ``dipole`` agrees
    with the root ``dipole_moment`` within 1e-6 debye
    (``E_ROOT_INCONSISTENT``). A malformed field is left to the schema check.
    """
    if not isinstance(payload, dict):
        return []
    sec = f"section {sid!r}: "
    errors: list[tuple[str, str]] = []
    if kind == "properties.elastic":
        for name in ("stiffness", "compliance"):
            m = _matrix_of(payload.get(name), 6, 6)
            bad = _asymmetry(m) if m is not None else None
            if bad is not None:
                i, j = bad
                errors.append((
                    "E_MEMBER_MALFORMED",
                    sec + f"{name} is not symmetric: entry ({i}, {j}) is "
                          f"{m[i][j]!r} but entry ({j}, {i}) is {m[j][i]!r}; "
                          f"they must agree within 1e-6 of the largest entry "
                          f"(spec § 5.13)"))
    elif kind == "properties.dielectric":
        tensors = {name: _matrix_of(payload.get(name), 3, 3)
                   for name in ("electronic", "ionic", "static")}
        if all(m is not None for m in tensors.values()):
            e, i_, s = tensors["electronic"], tensors["ionic"], tensors["static"]
            off = next(((a, b) for a in range(3) for b in range(3)
                        if abs(s[a][b] - (e[a][b] + i_[a][b])) > _PROPERTY_TOLERANCE),
                       None)
            if off is not None:
                a, b = off
                errors.append((
                    "E_MEMBER_MALFORMED",
                    sec + f"static[{a}][{b}] is {s[a][b]!r} but electronic + "
                          f"ionic is {e[a][b] + i_[a][b]!r}; static must equal "
                          f"their sum within 1e-6 (spec § 5.13)"))
        born = payload.get("born_charges")
        if isinstance(born, list):
            ref = section.get("structure_ref", only_structure)
            n_struct = n_atoms_of.get(ref) if isinstance(ref, str) else None
            if n_struct is not None and len(born) != n_struct:
                errors.append((
                    "E_MEMBER_MALFORMED",
                    sec + f"born_charges has {len(born)} tensors but structure "
                          f"{ref!r} has {n_struct} atoms; there must be one per "
                          f"atom (spec § 5.13)"))
    else:   # properties.electric
        dipole = payload.get("dipole")
        root = manifest.get("dipole_moment")
        if isinstance(dipole, dict) and isinstance(root, dict):
            vec, root_vec = dipole.get("vector"), root.get("vector_debye")
            if (isinstance(vec, list) and isinstance(root_vec, list)
                    and len(vec) == len(root_vec) == 3
                    and all(_is_number(x) and math.isfinite(x)
                            for x in vec + root_vec)
                    and any(abs(a - b) > _PROPERTY_TOLERANCE
                            for a, b in zip(vec, root_vec))):
                errors.append((
                    "E_ROOT_INCONSISTENT",
                    sec + f"dipole.vector {vec!r} differs from the root "
                          f"dipole_moment.vector_debye {root_vec!r}; the two "
                          f"state one dipole and must agree within 1e-6 debye "
                          f"(spec § 5.13)"))
            total, root_total = dipole.get("total"), root.get("total_debye")
            if (_is_number(total) and _is_number(root_total)
                    and math.isfinite(total) and math.isfinite(root_total)
                    and abs(total - root_total) > _PROPERTY_TOLERANCE):
                errors.append((
                    "E_ROOT_INCONSISTENT",
                    sec + f"dipole.total {total!r} differs from the root "
                          f"dipole_moment.total_debye {root_total!r}; the two "
                          f"state one dipole and must agree within 1e-6 debye "
                          f"(spec § 5.13)"))
    return errors



def _is_integral(x: Any) -> bool:
    """A JSON number with an integral value (``2`` or ``2.0``), not a bool."""
    return _is_number(x) and float(x).is_integer()


def _bloch_errors(sid: str, section: dict, members: Any, docs: dict[str, Any],
                  kinds_by_id: dict[str, str], periodic_of: dict[str, bool]
                  ) -> list[tuple[str, str, str]]:
    """Shape and structure rules of a ``wavefunction.bloch`` section that JSON
    Schema cannot express (spec § 5.3.1).

    ``mo_coefficients`` is ``[n_k, n_spin, n_mo, n_ao, 2]`` with n_spin 1 or 2;
    ``energies`` and ``occupations`` repeat its ``[n_k, n_spin, n_mo]`` (each
    other's, when the coefficients are malformed); the ``kpoints`` payload
    lists n_k k-points and, when present, n_k ``labels`` and n_k ``weights``
    that sum to 1 within 1e-6; ``mo_metadata`` states the same ``n_spin`` (1
    for ``restricted``, 2 for ``unrestricted``), ``n_mo`` and ``n_ao``, and
    ``n_mo_per_k`` has n_k entries of at most n_mo; ``basis.n_ao`` equals n_ao.
    The section's structure, named by the section-level ``structure_ref``
    (whose resolution the generic rule of spec § 4.1 checks) or else by
    ``basis.structure_ref``, must be a ``structure`` section with at least one
    periodic axis. Returns ``(code, location, message)`` triples; the location
    is relative to the section, ``members/<role>`` or ``structure_ref``.
    """
    if not isinstance(members, dict):
        return []
    errors: list[tuple[str, str, str]] = []

    def err(code: str, where: str, msg: str) -> None:
        errors.append((code, where, f"section {sid!r}: {msg} (spec § 5.3.1)"))

    coeff = _shape_of(members.get("mo_coefficients"))
    if coeff is not None and (len(coeff) != 5 or coeff[4] != 2
                              or coeff[1] not in (1, 2)):
        err("E_MEMBER_MALFORMED", "members/mo_coefficients",
            f"mo_coefficients shape {coeff} is not [n_k, n_spin, n_mo, n_ao, 2] "
            f"with n_spin 1 or 2")
        coeff = None
    ref = coeff[:3] if coeff is not None else None   # [n_k, n_spin, n_mo]
    ref_role = "mo_coefficients"
    for role in ("energies", "occupations"):
        shape = _shape_of(members.get(role))
        if shape is None:
            continue
        if len(shape) != 3:
            err("E_MEMBER_MALFORMED", f"members/{role}",
                f"{role} shape {shape} is not [n_k, n_spin, n_mo]")
        elif ref is None:
            ref, ref_role = shape, role
        elif shape != ref:
            err("E_MEMBER_MALFORMED", f"members/{role}",
                f"{role} shape {shape} differs from the [n_k, n_spin, n_mo] = "
                f"{ref} of {ref_role}")

    kp = docs.get("kpoints")
    if isinstance(kp, dict):
        points = kp.get("kpoints")
        n_k = len(points) if isinstance(points, list) else None
        if n_k is not None and ref is not None and n_k != ref[0]:
            err("E_MEMBER_MALFORMED", "members/kpoints",
                f"kpoints lists {n_k} k-points but axis 0 of {ref_role} is "
                f"{ref[0]}")
        for key in ("weights", "labels"):
            arr = kp.get(key)
            if isinstance(arr, list) and n_k is not None and len(arr) != n_k:
                err("E_MEMBER_MALFORMED", "members/kpoints",
                    f"kpoints.{key} has {len(arr)} entries for {n_k} k-points")
        weights = kp.get("weights")
        if isinstance(weights, list) and all(_is_number(w) for w in weights):
            total = math.fsum(weights)
            if abs(total - 1.0) > _BLOCH_WEIGHT_TOLERANCE:
                err("E_MEMBER_MALFORMED", "members/kpoints",
                    f"kpoints.weights sum to {total:.9g}, not to 1 within 1e-6")

    meta = docs.get("mo_metadata")
    if isinstance(meta, dict):
        spin, n_spin = meta.get("spin"), meta.get("n_spin")
        want = {"restricted": 1, "unrestricted": 2}.get(spin) \
            if isinstance(spin, str) else None
        if want is not None and _is_integral(n_spin) and n_spin != want:
            err("E_MEMBER_MALFORMED", "members/mo_metadata",
                f"mo_metadata.spin {spin!r} has n_spin {want}, not {n_spin}")
        for key, axis, dims, name in (("n_spin", 1, ref, ref_role),
                                      ("n_mo", 2, ref, ref_role),
                                      ("n_ao", 3, coeff, "mo_coefficients")):
            value = meta.get(key)
            if _is_integral(value) and dims is not None and value != dims[axis]:
                err("E_MEMBER_MALFORMED", "members/mo_metadata",
                    f"mo_metadata.{key} is {value} but axis {axis} of {name} "
                    f"is {dims[axis]}")
        per_k = meta.get("n_mo_per_k")
        if isinstance(per_k, list):
            if ref is not None and len(per_k) != ref[0]:
                err("E_MEMBER_MALFORMED", "members/mo_metadata",
                    f"mo_metadata.n_mo_per_k has {len(per_k)} entries but "
                    f"axis 0 of {ref_role} is {ref[0]}")
            n_mo = ref[2] if ref is not None else meta.get("n_mo")
            if _is_integral(n_mo) and any(_is_number(v) and v > n_mo
                                          for v in per_k):
                err("E_MEMBER_MALFORMED", "members/mo_metadata",
                    f"mo_metadata.n_mo_per_k has an entry above n_mo = {n_mo}")

    basis = docs.get("basis")
    if isinstance(basis, dict) and coeff is not None:
        n_ao = basis.get("n_ao")
        if _is_integral(n_ao) and n_ao != coeff[3]:
            err("E_MEMBER_MALFORMED", "members/basis",
                f"basis.n_ao is {n_ao} but axis 3 of mo_coefficients is "
                f"{coeff[3]}")

    if "structure_ref" in section:
        structure, where = section["structure_ref"], "structure_ref"
        if not isinstance(structure, str) or kinds_by_id.get(structure) != "structure":
            return errors       # reported by the generic rule (spec § 4.1)
    elif isinstance(basis, dict):
        structure, where = basis.get("structure_ref"), "members/basis"
        if not isinstance(structure, str) or kinds_by_id.get(structure) != "structure":
            err("E_REF_UNRESOLVED", where,
                f"basis.structure_ref {structure!r} does not name a structure "
                f"section")
            return errors
    else:
        return errors           # no basis payload: its own finding says why
    if periodic_of.get(structure) is False:
        err("E_MEMBER_MALFORMED", where,
            f"structure {structure!r} has no periodic axis; Bloch orbitals need "
            f"a periodic structure")
    return errors

def _atom_properties_meta_errors(sid: str, members: Any,
                                 docs: dict[str, Any]) -> list[tuple[str, str]]:
    """Custom roles of an ``atom_properties`` section need a ``meta`` entry.

    Every ``custom.<name>`` member must have an entry named ``custom.<name>``,
    a JSON object, in the section's ``meta`` member (spec § 5.7). A section
    without a ``meta`` member has none, so each custom role is reported. A
    ``meta`` member that is declared but could not be read (its own finding
    says why) is not judged. JSON Schema cannot compare two members, so this
    runs with or without ``jsonschema``. Returns ``(role, message)`` pairs.
    """
    if not isinstance(members, dict):
        return []
    custom = [r for r in members if isinstance(r, str) and r.startswith("custom.")]
    if not custom:
        return []
    sec = f"section {sid!r}: custom role "
    if "meta" not in members:
        return [(r, f"{sec}{r!r} has no meta member to describe it (spec § 5.7)")
                for r in custom]
    if "meta" not in docs:
        return []
    meta = docs["meta"]
    if not isinstance(meta, dict):
        return [(r, f"{sec}{r!r}: the meta member is not a JSON object "
                    f"(spec § 5.7)") for r in custom]
    return [(r, f"{sec}{r!r} has no entry in the meta member (spec § 5.7)")
            for r in custom if not isinstance(meta.get(r), dict)]


# --- wavefunction.gto (spec § 5.3, Appendix A.1) ---------------------------------

def _odd_double_factorial(n: int) -> int:
    """(2n - 1)!!, with (-1)!! = 1."""
    out = 1
    for k in range(1, 2 * n, 2):
        out *= k
    return out


def _basis_shells(basis: Any) -> Optional[list[tuple]]:
    """``(l, pure, exponents, coefficients)`` of every shell of a ``basis``
    payload, or None when the AO layout cannot be read (no shell list, or a
    shell whose ``l`` or ``pure`` is malformed). ``exponents`` and
    ``coefficients`` are None for a shell whose primitive arrays are unusable
    (not numbers, an exponent not above 0, or lengths that differ); the
    schema check reports those."""
    if not isinstance(basis, dict) or not isinstance(basis.get("shells"), list):
        return None
    default_pure = basis.get("pure", True)
    if not isinstance(default_pure, bool):
        return None
    out = []
    for sh in basis["shells"]:
        if not isinstance(sh, dict):
            return None
        l, pure = sh.get("l"), sh.get("pure", default_pure)
        if (not isinstance(l, int) or isinstance(l, bool) or not 0 <= l <= 7
                or not isinstance(pure, bool)):
            return None
        exps, coefs = sh.get("exponents"), sh.get("coefficients")
        usable = (isinstance(exps, list) and isinstance(coefs, list)
                  and len(exps) == len(coefs) and len(exps) > 0
                  and all(_is_number(a) and a > 0 for a in exps)
                  and all(_is_number(c) for c in coefs))
        out.append((l, pure, [float(a) for a in exps] if usable else None,
                    [float(c) for c in coefs] if usable else None))
    return out


def _shell_size(l: int, pure: bool) -> int:
    return 2 * l + 1 if pure else (l + 1) * (l + 2) // 2


def _basis_n_ao(basis: Any) -> Optional[int]:
    """n_ao of a ``basis`` payload: counted from its shells, or its ``n_ao``
    when the shells cannot be read (spec § 5.3); None when neither works."""
    shells = _basis_shells(basis)
    if shells is not None:
        return sum(_shell_size(l, pure) for l, pure, _e, _c in shells)
    n_ao = basis.get("n_ao") if isinstance(basis, dict) else None
    if isinstance(n_ao, int) and not isinstance(n_ao, bool) and n_ao >= 0:
        return n_ao
    return None


def _axial_self_overlap(l: int, exponents: list[float],
                        coefficients: list[float]) -> float:
    """Self-overlap of the axial component ``(l, 0, 0)`` of a contracted shell
    on the ``N_i``-normalized primitives of Appendix A.1:
    ``sum_ij c_i c_j (2 sqrt(a_i a_j) / (a_i + a_j))^(l + 3/2)``."""
    total = 0.0
    for a, ca in zip(exponents, coefficients):
        for b, cb in zip(exponents, coefficients):
            total += ca * cb * (2.0 * math.sqrt(a * b) / (a + b)) ** (l + 1.5)
    return total


def _shell_components(l: int, pure: bool) -> list[tuple[str, float]]:
    """``(name, factor)`` of each component of a shell, in Appendix A.2 order:
    the diagonal overlap of a component is ``factor * S_axial`` (Appendix
    A.1, "Verifying with overlap")."""
    if pure:
        return [(f"m = {m:+d}" if m else "m = 0", 1.0) for m in range(-l, l + 1)]
    out = []
    for lx in range(l, -1, -1):
        for ly in range(l - lx, -1, -1):
            lz = l - lx - ly
            factor = (_odd_double_factorial(lx) * _odd_double_factorial(ly)
                      * _odd_double_factorial(lz)) / _odd_double_factorial(l)
            name = "x" * lx + "y" * ly + "z" * lz
            out.append((f"Cartesian {name or 's'} [{lx}, {ly}, {lz}]", factor))
    return out


def _normalization_errors(sid: str, basis: Any, overlap: Any,
                          raw: Optional[bytes]) -> list[str]:
    """The Appendix A.1 check of a ``wavefunction.gto`` ``overlap`` member.

    The self-overlap of every AO follows from ``basis.json``: ``S_axial`` of
    its shell for a spherical component, ``S_axial`` times
    ``(2l_x-1)!!(2l_y-1)!!(2l_z-1)!! / (2l-1)!!`` for a Cartesian one. Each must
    agree with the stored diagonal at that AO's index to a relative tolerance
    of 1e-6. Runs only on a well-formed float64 ``[n_ao, n_ao]`` overlap of a
    basis whose AO layout can be read; a shell with unusable primitive arrays
    is skipped. Returns at most one message, naming the first disagreement.
    """
    shells = _basis_shells(basis)
    shape = _shape_of(overlap)
    if shells is None or shape is None or raw is None:
        return []
    n_ao = sum(_shell_size(l, pure) for l, pure, _e, _c in shells)
    if (shape != [n_ao, n_ao] or overlap.get("dtype") != "float64"
            or len(raw) != 8 * n_ao * n_ao):
        return []
    bad: list[tuple] = []
    ao = 0
    for k, (l, pure, exps, coefs) in enumerate(shells):
        components = _shell_components(l, pure)
        if exps is not None:
            s_axial = _axial_self_overlap(l, exps, coefs)
            for j, (name, factor) in enumerate(components):
                i = ao + j
                want = s_axial * factor
                got = struct.unpack_from("<d", raw, 8 * (i * n_ao + i))[0]
                if not abs(got - want) <= _NORMALIZATION_RTOL * abs(want):
                    bad.append((i, k, l, name, got, want))
        ao += len(components)
    if not bad:
        return []
    i, k, l, name, got, want = bad[0]
    return [f"section {sid!r}: {len(bad)} of {n_ao} overlap diagonal entries "
            f"disagree with the self-overlaps Appendix A.1 gives for the basis; "
            f"first: AO {i} (shell {k}, l = {l}, {name}) is {got:.10g} but the "
            f"basis gives {want:.10g}. The basis.json coefficients are not on "
            f"N_i-normalized primitives, or the AO order disagrees with "
            f"Appendix A.2 (spec § 5.3)"]


def _coefficients_complex(meta: Any) -> Optional[bool]:
    """Whether ``mo_metadata`` declares complex coefficients (spec § 5.3):
    ``complex``, else ``coefficient_encoding``, else real. None when the
    declaration is malformed (the schema check reports it)."""
    if not isinstance(meta, dict):
        return None
    if "complex" in meta:
        value = meta["complex"]
        return value if isinstance(value, bool) else None
    if "coefficient_encoding" in meta:
        encoding = meta["coefficient_encoding"]
        if encoding == "complex_split_last_axis":
            return True
        return False if encoding == "real" else None
    return False


def _unrestricted_nto_errors(meta: Any, members: Any) -> list[str]:
    """Local unrestricted NTO rules; malformed schema inputs are skipped."""
    if (not isinstance(meta, dict) or "nto_role" not in meta
            or meta.get("orbital_kind") != "nto" or meta.get("spin") != "unrestricted"):
        return []
    errors = []
    rows = _mo_rows(members)
    for ch in ("alpha", "beta"):
        block = meta.get(ch)
        if not isinstance(block, dict):
            continue
        n = rows.get(ch)
        if n is None:
            errors.append(f"unrestricted NTOs need mo_coefficients_{ch} (spec § 5.3)")
        else:
            if any(isinstance(block.get(key), list) and len(block[key]) != n
                   for key in ("energies", "occupations")) or (
                       _is_integral(block.get("n_mo")) and block["n_mo"] != n):
                errors.append(f"NTO {ch} metadata must have one entry per coefficient row (spec § 5.3)")
        occ = block.get("occupations")
        if isinstance(occ, list) and all(_is_number(x) for x in occ):
            if any(x < 0 for x in occ) or any(a < b for a, b in zip(occ, occ[1:])):
                errors.append(f"NTO {ch} weights must be nonnegative and non-increasing (spec § 5.3)")
    return errors


def _nto_pair_errors(meta: Any, members: Any, ntos: list[tuple],
                     kinds: dict, indices: dict) -> list[tuple[str, str]]:
    """Forward references and equal per-spin rows/weights, scoped to unrestricted NTOs."""
    if (not isinstance(meta, dict) or "nto_role" not in meta
            or meta.get("orbital_kind") != "nto" or meta.get("spin") != "unrestricted"):
        return []
    ref, index, role = meta.get("state_ref"), meta.get("state_index"), meta.get("nto_role")
    if not isinstance(ref, str) or not _is_integral(index) or role not in ("hole", "particle"):
        return []  # payload schema reports these
    errors = []
    if kinds.get(ref) != "excited_states" or (indices.get(ref) is not None and index not in indices[ref]):
        errors.append(("E_REF_UNRESOLVED", "NTO state_ref/state_index must name an excited state (spec § 5.3)"))
    opposite = "particle" if role == "hole" else "hole"
    partners = [(m, ms) for _, _, m, ms in ntos if isinstance(m, dict)
                and m.get("orbital_kind") == "nto" and m.get("state_ref") == ref
                and m.get("state_index") == index and m.get("nto_role") == opposite]
    if len(partners) != 1:
        errors.append(("E_REF_UNRESOLVED", "NTOs must have one opposite-role section for the same state (spec § 5.3)"))
        return errors
    other, other_members = partners[0]
    if other.get("spin") != "unrestricted":
        errors.append(("E_MEMBER_MALFORMED", "NTO partners must have the same spin layout (spec § 5.3)"))
        return errors
    rows, other_rows = _mo_rows(members), _mo_rows(other_members)
    for ch in ("alpha", "beta"):
        a, b = meta.get(ch), other.get(ch)
        if ch in rows and ch in other_rows and rows[ch] != other_rows[ch]:
            errors.append(("E_MEMBER_MALFORMED", f"NTO partners must have the same {ch} row count (spec § 5.3)"))
        if isinstance(a, dict) and isinstance(b, dict):
            av, bv = a.get("occupations"), b.get("occupations")
            if (isinstance(av, list) and isinstance(bv, list)
                    and all(_is_number(x) for x in av + bv) and av != bv):
                errors.append(("E_MEMBER_MALFORMED", f"NTO partners must have identical {ch} weights (spec § 5.3)"))
    return errors


def _wavefunction_errors(sid: str, members: Any, docs: dict[str, Any],
                         provenance: Any) -> list[tuple[str, str, str]]:
    """Cross-member rules of a ``wavefunction.gto`` section (spec § 5.3).

    The rank of every coefficient member against ``complex``; the length of
    every ``symmetry_labels`` against its ``energies``; the ``[n_ao, n_ao]``
    shape of ``overlap`` and of the density members, and the pairing of the
    density members with the coefficient members; ``basis.n_electrons``
    against the occupations (``E_MEMBER_MALFORMED``) and against
    ``provenance.n_electrons`` (``E_ROOT_INCONSISTENT``). A rule whose inputs
    are missing or malformed is skipped; the schema check reports those.
    Returns ``(code, role, message)`` triples.
    """
    if not isinstance(members, dict):
        return []
    errors: list[tuple[str, str, str]] = []
    sec = f"section {sid!r}: "
    basis, meta = docs.get("basis"), docs.get("mo_metadata")

    cplx = _coefficients_complex(meta)
    if cplx is not None:
        for role in _WF_COEFFICIENT_ROLES:
            shape = _shape_of(members.get(role)) if role in members else None
            if shape is None:
                continue
            if cplx and (len(shape) != 3 or shape[2] != 2):
                errors.append(("E_MEMBER_MALFORMED", role,
                               sec + f"{role} shape {shape} has no trailing "
                                     f"[re, im] axis of 2, but mo_metadata "
                                     f"declares complex coefficients "
                                     f"(spec § 5.3)"))
            elif not cplx and len(shape) != 2:
                errors.append(("E_MEMBER_MALFORMED", role,
                               sec + f"{role} shape {shape} is not "
                                     f"[n_mo, n_ao], and mo_metadata does not "
                                     f"declare complex coefficients "
                                     f"(spec § 5.3)"))

    if isinstance(meta, dict):
        sets = [("", meta)] + [(f"{ch}.", meta[ch]) for ch in ("alpha", "beta")
                               if isinstance(meta.get(ch), dict)]
        for prefix, block in sets:
            labels, energies = block.get("symmetry_labels"), block.get("energies")
            if (isinstance(labels, list) and isinstance(energies, list)
                    and len(labels) != len(energies)):
                errors.append(("E_MEMBER_MALFORMED", "mo_metadata",
                               sec + f"mo_metadata {prefix}symmetry_labels has "
                                     f"{len(labels)} entries but "
                                     f"{prefix}energies has {len(energies)}; "
                                     f"there is one label per MO (spec § 5.3)"))

    n_ao = _basis_n_ao(basis)
    for role in _WF_AO_MATRIX_ROLES:
        shape = _shape_of(members.get(role)) if role in members else None
        if shape is not None and n_ao is not None and shape != [n_ao, n_ao]:
            errors.append(("E_MEMBER_MALFORMED", role,
                           sec + f"{role} shape {shape} is not [n_ao, n_ao] = "
                                 f"[{n_ao}, {n_ao}] for the basis "
                                 f"(spec § 5.3)"))
    unrestricted = ("mo_coefficients_alpha" in members
                    and "mo_coefficients_beta" in members)
    if "density_matrix" in members and "mo_coefficients" not in members:
        errors.append(("E_MEMBER_MALFORMED", "density_matrix",
                       sec + "density_matrix accompanies mo_coefficients only; "
                             "an unrestricted wavefunction carries "
                             "density_matrix_alpha and density_matrix_beta "
                             "(spec § 5.3)"))
    for role, other in (("density_matrix_alpha", "density_matrix_beta"),
                        ("density_matrix_beta", "density_matrix_alpha")):
        if role not in members:
            continue
        if not unrestricted:
            errors.append(("E_MEMBER_MALFORMED", role,
                           sec + f"{role} accompanies mo_coefficients_alpha "
                                 f"and mo_coefficients_beta only; a single-set "
                                 f"wavefunction carries density_matrix "
                                 f"(spec § 5.3)"))
        elif other not in members:
            errors.append(("E_MEMBER_MALFORMED", role,
                           sec + f"{role} without {other}; the two come "
                                 f"together (spec § 5.3)"))

    n_el = basis.get("n_electrons") if isinstance(basis, dict) else None
    if not _is_number(n_el):
        return errors
    if (isinstance(meta, dict)
            and meta.get("occupation_semantics") != "transition_weight"):
        lists = [meta["occupations"]] if "occupations" in meta else []
        lists += [meta[ch]["occupations"] for ch in ("alpha", "beta")
                  if isinstance(meta.get(ch), dict) and "occupations" in meta[ch]]
        if lists and all(isinstance(occ, list) and all(map(_is_number, occ))
                         for occ in lists):
            total = math.fsum(x for occ in lists for x in occ)
            if not abs(total - n_el) <= _N_ELECTRONS_TOLERANCE:
                errors.append(("E_MEMBER_MALFORMED", "basis",
                               sec + f"basis.n_electrons is {n_el} but the "
                                     f"mo_metadata occupations sum to "
                                     f"{total:.10g}; they must agree within "
                                     f"1e-6 (spec § 5.3)"))
    prov_n = (provenance.get("n_electrons") if isinstance(provenance, dict)
              else None)
    if _is_number(prov_n) and prov_n != n_el:
        errors.append(("E_ROOT_INCONSISTENT", "basis",
                       sec + f"basis.n_electrons is {n_el} but "
                             f"provenance.n_electrons is {prov_n}; they state "
                             f"one electron count and must be equal "
                             f"(spec § 3.2, § 5.3)"))
    return errors


def _volume_errors(sid: str, kind: str, members: Any, docs: dict[str, Any],
                   kinds_by_id: dict[str, str]) -> list[tuple[str, str, str]]:
    """Data shape and series axis of a volume section (spec § 5.2).

    The ``data`` of ``volume.laplacian`` is ``[nx, ny, nz]``, of
    ``volume.vector`` ``[nx, ny, nz, 3]``, of ``volume.complex``
    ``[nx, ny, nz, 2]`` and of ``volume.series`` ``[n, nx, ny, nz]``. A
    series ``axis`` has one ``values`` entry, and one ``labels`` entry when it
    has labels, per slab, and its ``trajectory_ref`` names a ``trajectory`` or
    ``reaction.path`` section. JSON Schema cannot compare two members, so this
    runs with or without ``jsonschema``. Returns ``(code, role, message)``.
    """
    rule = _VOLUME_DATA_SHAPES.get(kind)
    if rule is None or not isinstance(members, dict):
        return []
    rank, trailing, form = rule
    sec = f"section {sid!r}: "
    errors: list[tuple[str, str, str]] = []
    shape = _shape_of(members.get("data"))
    if shape is not None and (len(shape) != rank or (
            trailing is not None and shape[-1] != trailing)):
        errors.append(("E_MEMBER_MALFORMED", "data",
                       sec + f"{kind} data shape {shape} is not {form} "
                             f"(spec § 5.2)"))
        shape = None
    axis = docs.get("axis")
    if kind != "volume.series" or not isinstance(axis, dict):
        return errors
    for key in ("values", "labels"):
        seq = axis.get(key)
        if shape is not None and isinstance(seq, list) and len(seq) != shape[0]:
            errors.append(("E_MEMBER_MALFORMED", "axis",
                           sec + f"axis.{key} has {len(seq)} entries but data "
                                 f"{shape} has {shape[0]} slabs; it needs one "
                                 f"per slab (spec § 5.2)"))
    if "trajectory_ref" in axis:
        ref = axis["trajectory_ref"]
        if not isinstance(ref, str) or kinds_by_id.get(ref) not in _SERIES_REF_KINDS:
            errors.append(("E_REF_UNRESOLVED", "axis",
                           sec + f"axis.trajectory_ref {ref!r} does not name a "
                                 f"trajectory or reaction.path section "
                                 f"(spec § 5.2)"))
    return errors


def _grid_data_unit_error(sid: str, kind: Any, grid: Any) -> Optional[str]:
    """A grid ``units.data`` that contradicts the unit spec § 4.3 fixes for
    the volume kind, or None. Only a unit name contradicts: a string outside
    the vocabulary is ``W_UNIT_UNKNOWN`` alone."""
    fixed = _VOLUME_DATA_UNITS.get(kind) if isinstance(kind, str) else None
    units = grid.get("units") if isinstance(grid, dict) else None
    if fixed is None or not isinstance(units, dict):
        return None
    unit = units.get("data")
    if not isinstance(unit, str) or unit == fixed or unit not in qvf_kinds.UNITS:
        return None
    return (f"section {sid!r}: grid units.data {unit!r} contradicts the unit "
            f"{fixed!r} that spec § 4.3 fixes for {kind} data")


def _is_bool3(value: Any) -> bool:
    return (isinstance(value, list) and len(value) == 3
            and all(isinstance(x, bool) for x in value))


def _vec3(value: Any) -> Optional[list[float]]:
    """``value`` as three floats when it is a list of three numbers."""
    if (isinstance(value, list) and len(value) == 3
            and all(_is_number(x) for x in value)):
        return [float(x) for x in value]
    return None


def _show_vec(v: list[float]) -> str:
    return "[" + ", ".join(f"{x:.6g}" for x in v) + "]"


def _grid_errors(sid: str, section: dict, grid: Any,
                 kinds_by_id: dict[str, str],
                 cells_of: dict[str, tuple[Any, Any]],
                 only_structure: Optional[str]) -> list[tuple[str, str]]:
    """Reference and periodic-span rules of a grid payload (spec § 5.2).

    A grid ``structure_ref`` names a ``structure`` section. The grid's
    structure is the section-level ``structure_ref`` when the section carries
    one, else the grid's ``structure_ref``, else the archive's only
    ``structure`` section. When it resolves to a structure whose payload was
    read (``cells_of``: id -> ``(pbc, lattice_vectors)``), every axis a with
    ``periodic[a]`` true must have ``pbc[a]`` true, and the span
    ``voxel_vectors[a] × (shape[a] − 1)`` (``includes_endpoint`` true, the
    default) or ``voxel_vectors[a] × shape[a]`` (false) must equal lattice
    row a converted to bohr within 1e-6 relative. Malformed fields are left
    to the schema. Returns ``(code, message)`` pairs.
    """
    if not isinstance(grid, dict):
        return []
    sec = f"section {sid!r}: "
    errors: list[tuple[str, str]] = []
    ref = grid.get("structure_ref")
    if "structure_ref" in grid and (not isinstance(ref, str)
                                    or kinds_by_id.get(ref) != "structure"):
        errors.append(("E_REF_UNRESOLVED",
                       sec + f"grid structure_ref {ref!r} does not name a "
                             f"structure section (spec § 5.2)"))
    periodic = grid.get("periodic")
    if not _is_bool3(periodic) or not any(periodic):
        return errors
    if "structure_ref" in section:
        structure = (section["structure_ref"]
                     if isinstance(section["structure_ref"], str) else None)
    elif isinstance(ref, str):
        structure = ref
    else:
        structure = only_structure
    if structure is None or structure not in cells_of:
        return errors
    pbc, lattice = cells_of[structure]
    if pbc is None:
        pbc = [False, False, False]          # absent pbc: a molecule (§ 5.1)
    endpoint = grid.get("includes_endpoint", True)
    shape = grid.get("shape")
    vox = grid.get("voxel_vectors")
    if (not _is_bool3(pbc) or not isinstance(endpoint, bool)
            or not isinstance(shape, list) or len(shape) != 3
            or not all(isinstance(n, int) and not isinstance(n, bool)
                       for n in shape)
            or not isinstance(vox, list) or len(vox) != 3):
        return errors
    for a in range(3):
        if not periodic[a]:
            continue
        if not pbc[a]:
            errors.append(("E_MEMBER_MALFORMED",
                           sec + f"grid axis {a} is periodic but structure "
                                 f"{structure!r} has pbc[{a}] false; a grid "
                                 f"axis can be periodic only along a periodic "
                                 f"axis of its structure (spec § 5.2)"))
            continue
        step = _vec3(vox[a])
        row = (_vec3(lattice[a]) if isinstance(lattice, list)
               and len(lattice) == 3 else None)
        if step is None or row is None:
            continue
        n = shape[a] - (1 if endpoint else 0)
        span = [x * n for x in step]
        target = [x * _BOHR_PER_ANGSTROM for x in row]
        diff = math.sqrt(sum((s - t) ** 2 for s, t in zip(span, target)))
        if diff > _GRID_SPAN_TOLERANCE * math.sqrt(sum(t * t for t in target)):
            errors.append((
                "E_MEMBER_MALFORMED",
                sec + f"periodic grid axis {a} spans {_show_vec(span)} bohr "
                      f"(voxel_vectors[{a}] × {n}, includes_endpoint "
                      f"{str(endpoint).lower()}) but lattice vector {a} of "
                      f"structure {structure!r} is {_show_vec(target)} bohr; "
                      f"a periodic grid spans its cell within 1e-6 relative "
                      f"(spec § 5.2)"))
    return errors


_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
_JPEG_SIGNATURE = b"\xff\xd8\xff"


def _preview_member_errors(sid: str, role: str, member: dict,
                           raw: bytes) -> list[tuple[str, str]]:
    """``preview`` (spec § 5.15): the bytes are what the ``media_type`` says."""
    media_type = member.get("media_type")
    where = f"section {sid!r} member {role!r}"
    if media_type == "image/png" and not raw.startswith(_PNG_SIGNATURE):
        return [("E_MEMBER_MALFORMED", f"{where}: image/png bytes do not begin "
                 f"with the PNG signature (spec § 5.15)")]
    if media_type == "image/jpeg" and not raw.startswith(_JPEG_SIGNATURE):
        return [("E_MEMBER_MALFORMED", f"{where}: image/jpeg bytes do not begin "
                 f"with FF D8 FF (spec § 5.15)")]
    if media_type == "image/svg+xml":
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            return [("E_UTF8", f"{where}: image/svg+xml bytes do not decode as "
                     f"UTF-8 ({exc}; spec § 5.15)")]
        if "<svg" not in text:
            return [("E_MEMBER_MALFORMED", f"{where}: image/svg+xml bytes "
                     f"contain no <svg element (spec § 5.15)")]
    return []


def _reaction_scheme_errors(sid: str, payload: Any,
                            kinds_by_id: dict[str, str]) -> list[tuple[str, str]]:
    """``reaction.scheme`` (spec § 5.14): participants name ``structure``
    sections, step ids are unique, and no structure repeats within a role."""
    out: list[tuple[str, str]] = []
    steps = payload.get("steps") if isinstance(payload, dict) else None
    if not isinstance(steps, list):
        return out
    seen: set[str] = set()
    for k, step in enumerate(steps):
        if not isinstance(step, dict):
            continue
        step_id = step.get("id")
        if isinstance(step_id, str):
            if step_id in seen:
                out.append(("E_MEMBER_MALFORMED",
                            f"section {sid!r}: step id {step_id!r} is used by "
                            f"more than one step (spec § 5.14)"))
            seen.add(step_id)
        for role in ("reactants", "products", "reagents"):
            items = step.get(role)
            if not isinstance(items, list):
                continue
            refs = []
            for item in items:
                ref = item.get("structure_ref") if isinstance(item, dict) else None
                if not isinstance(ref, str):
                    continue
                refs.append(ref)
                if kinds_by_id.get(ref) != "structure":
                    out.append(("E_REF_UNRESOLVED",
                                f"section {sid!r}: step {k} {role} structure_ref "
                                f"{ref!r} does not name a structure section "
                                f"(spec § 5.14)"))
            for ref in sorted({r for r in refs if refs.count(r) > 1}):
                out.append(("E_MEMBER_MALFORMED",
                            f"section {sid!r}: step {k} names {ref!r} more "
                            f"than once among its {role} (spec § 5.14)"))
        if isinstance(step.get("energies"), dict):
            out += _reaction_energies_errors(f"section {sid!r}: step {k}",
                                             step, kinds_by_id)
    return out


#: Reaction energy quantities and the term field holding each absolute value
#: (spec § 5.14, *Reaction energies*).
_REACTION_QUANTITIES = (("energy", "energy_eh"), ("enthalpy", "enthalpy_eh"),
                        ("gibbs_free_energy", "gibbs_free_energy_eh"))
_HARTREE_TO_KCAL_MOL = 627.509474          # CODATA 2018
_REACTION_DELTA_TOL_EH = 1e-8
_REACTION_KCAL_REL_TOL = 1e-6
_REACTION_COEFFICIENT_TOL = 1e-9


def _number(value: Any) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value))


def _reaction_energies_errors(where: str, step: dict,
                              kinds_by_id: Optional[dict[str, str]]
                              ) -> list[tuple[str, str]]:
    """A step's ``energies`` (spec § 5.14, *Reaction energies*): no energies
    on a resonance arrow; one term per reactant and product with the signed
    coefficient; each delta equals the sum of its terms and its kcal/mol
    value the conversion; one method; ``run_ref`` names a ``run.record``
    (``E_REF_UNRESOLVED``, only with ``kinds_by_id``). Shapes the schema
    rejects are skipped."""
    out: list[tuple[str, str]] = []
    bad = lambda msg: out.append(("E_MEMBER_MALFORMED",  # noqa: E731
                                  f"{where} energies: {msg} (spec § 5.14)"))
    energies = step["energies"]
    if step.get("arrow") == "resonance":
        bad("a resonance arrow has no reaction energy")
    if not isinstance(energies.get("terms"), list):
        return out
    expected: dict[str, float] = {}
    for role, sign in (("reactants", -1.0), ("products", 1.0)):
        for item in step.get(role) or ():
            if not isinstance(item, dict) or not isinstance(item.get("structure_ref"), str):
                continue
            c = item.get("coefficient", 1)
            if _number(c):
                ref = item["structure_ref"]
                expected[ref] = expected.get(ref, 0.0) + sign * c
    expected = {r: c for r, c in expected.items() if c != 0}
    terms = [t for t in energies["terms"] if isinstance(t, dict)]
    seen: set[str] = set()
    for t, term in enumerate(terms):
        ref, c = term.get("structure_ref"), term.get("coefficient")
        if not isinstance(ref, str):
            continue
        if ref in seen:
            bad(f"structure {ref!r} has more than one term")
            continue
        seen.add(ref)
        if ref not in expected:
            bad(f"term {t} names {ref!r}, which is not a reactant or product")
        elif _number(c) and abs(c - expected[ref]) > _REACTION_COEFFICIENT_TOL * max(1.0, abs(c)):
            bad(f"term {t} coefficient {c!r} for {ref!r}, but the step "
                f"gives {expected[ref]!r}")
        run_ref = term.get("run_ref")
        if (kinds_by_id is not None and isinstance(run_ref, str)
                and kinds_by_id.get(run_ref) != "run.record"):
            out.append(("E_REF_UNRESOLVED",
                        f"{where} energies: term {t} run_ref {run_ref!r} does "
                        f"not name a run.record section (spec § 5.14)"))
    for ref in sorted(set(expected) - seen):
        bad(f"no term for participant {ref!r}")
    methods = [term.get("method") for term in terms if isinstance(term.get("method"), dict)]
    if any(m != methods[0] for m in methods[1:]):
        bad("the terms have different methods")
    deltas = energies.get("deltas")
    if not isinstance(deltas, dict):
        return out
    for quantity, field in _REACTION_QUANTITIES:
        delta = deltas.get(quantity)
        if not isinstance(delta, dict) or not _number(delta.get("value_eh")):
            continue
        values = [(term.get("coefficient"), term.get(field)) for term in terms]
        missing = [t for t, (_, v) in enumerate(values) if not _number(v)]
        if missing:
            bad(f"{quantity} is given but term {missing[0]} has no {field}")
        elif all(_number(c) for c, _ in values):
            total = math.fsum(c * v for c, v in values)
            if abs(total - delta["value_eh"]) > _REACTION_DELTA_TOL_EH:
                bad(f"{quantity} value_eh {delta['value_eh']!r} is not the sum "
                    f"of its terms, {total!r}")
        kcal = delta.get("kcal_mol")
        if _number(kcal):
            converted = delta["value_eh"] * _HARTREE_TO_KCAL_MOL
            if abs(kcal - converted) > _REACTION_KCAL_REL_TOL * max(1.0, abs(converted)):
                bad(f"{quantity} kcal_mol {kcal!r} is not value_eh x "
                    f"627.509474 = {converted!r}")
    return out


def _identity_errors(sid: str, section: dict, sum_z: Optional[int],
                     provenance: Any, only_structure: bool
                     ) -> list[tuple[str, str, str]]:
    """The chemical identity of a ``structure`` section (spec § 5.1): a charge
    and multiplicity possible for its atoms, and, in an archive with one
    structure, equal to ``provenance``. Yields (code, field, message)."""
    out: list[tuple[str, str, str]] = []
    charge = section.get("charge")
    mult = section.get("multiplicity")
    charge = charge if _is_int(charge) else None
    mult = mult if _is_int(mult) and mult >= 1 else None
    if charge is not None and sum_z is not None:
        electrons = sum_z - charge
        if electrons < 0:
            out.append(("E_CHARGE_MULTIPLICITY", "charge",
                        f"structure {sid!r}: charge {charge} leaves {electrons} "
                        f"electrons for atoms whose atomic numbers sum to "
                        f"{sum_z} (spec § 5.1)"))
        elif mult is not None and (mult - 1 > electrons
                                   or (electrons - (mult - 1)) % 2):
            out.append(("E_CHARGE_MULTIPLICITY", "multiplicity",
                        f"structure {sid!r}: multiplicity {mult} is impossible "
                        f"with {electrons} electrons (charge {charge}): the "
                        f"unpaired electrons must not exceed the electron count "
                        f"and must have its parity (spec § 5.1)"))
    if only_structure and isinstance(provenance, dict):
        for field, value in (("charge", charge), ("multiplicity", mult)):
            stated = provenance.get(field)
            if value is not None and _is_int(stated) and stated != value:
                out.append(("E_ROOT_INCONSISTENT", field,
                            f"structure {sid!r}: {field} {value} differs from "
                            f"provenance.{field} {stated}, which describes the "
                            f"archive's only structure (spec § 5.1)"))
    return out


class Level(IntEnum):
    """Cumulative stages, matching the C++ validator. All is the default."""
    Container = 0
    Manifest = 1
    Payloads = 2
    Semantics = 3
    References = 4
    All = 5


def _level(value):
    if isinstance(value, str):
        try:
            return Level[value.capitalize()]
        except KeyError:
            raise ValueError(f"unknown validation level {value!r}") from None
    return Level(value)


class _Findings:
    """Collects findings; derives the legacy ``errors`` / ``warnings`` lists."""

    def __init__(self) -> None:
        from ._legacy_writer import TOOLKIT_VERSION
        self.references = True
        self.items: list[dict[str, str]] = []
        self.context = {"reader": "qvf-python", "reader_version": TOOLKIT_VERSION,
                        "archive_format_version": None, "archive_schema_uri": None,
                        "schema_toolkit_version": TOOLKIT_VERSION,
                        "schema_id": None, "schema_sha256": None,
                        "schema_checked": False, "supported_format_versions": [1],
                        "unsupported_format_version": None,
                        "scientific_validity": "not_assessed"}
        try:
            with open(_SCHEMA_PATH, "rb") as source:
                raw = source.read()
            self.context.update(schema_id=json.loads(raw).get("$id"),
                                schema_sha256=hashlib.sha256(raw).hexdigest())
        except (OSError, ValueError):
            pass  # The existing schema-skipped finding remains authoritative.

    def error(self, code: str, where: str, message: str) -> None:
        self._add("error", code, where, message)

    def warning(self, code: str, where: str, message: str) -> None:
        self._add("warning", code, where, message)

    def _add(self, severity: str, code: str, where: str, message: str) -> None:
        if code == "E_REF_UNRESOLVED" and not self.references:
            return
        assert FINDING_CODES.get(code) == severity, (code, severity)
        self.items.append({"severity": severity, "code": code,
                           "where": where, "message": message})

    def report(self, n_sections: int) -> dict[str, Any]:
        errors = [f"{f['code']}: {f['message']}" for f in self.items
                  if f["severity"] == "error"]
        warnings = [f"{f['code']}: {f['message']}" for f in self.items
                    if f["severity"] == "warning"]
        return {"ok": not errors, "errors": errors, "warnings": warnings,
                "n_sections": n_sections, "findings": list(self.items),
                "validation_context": dict(self.context)}


def _context_summary(report):
    context = report["validation_context"]
    return (f"archive format={context['archive_format_version']!r}; "
            f"reader={context['reader']} {context['reader_version']}; "
            f"bundled schema={context['schema_id']} "
            f"sha256={context['schema_sha256']}; "
            f"schema checked={context['schema_checked']}. "
            "Format version is not the toolkit release: compare the producer's "
            "validator pin before attributing a cross-version rejection to the producer. "
            "Numerical accuracy and electronic stability are not assessed.")


def validate_qvf(path: str, *, limits: Optional[Limits] = None, level=Level.All) -> dict[str, Any]:
    """Validate a ``.qvf`` archive (QVF spec § 6).

    ``level`` selects cumulative Container, Manifest, Payloads, Semantics,
    References or All stages. Container never decompresses a member; Manifest
    reads only manifest.json. Full validation is unchanged by default.

    Returns ``{ok, errors, warnings, n_sections, findings}``. ``findings`` is a
    list of ``{"severity", "code", "where", "message"}`` dicts whose ``code`` is
    a stable identifier from the spec's Appendix B; ``errors`` / ``warnings``
    are the same findings as ``"CODE: message"`` strings.

    Container checks: ZIP readability, unique entry names (several
    ``manifest.json`` entries are an archive updated in place, spec § 2.3: the
    last one is the manifest and each earlier one is the warning
    ``W_ORPHAN_MEMBER``), legal entry paths, per-member / total
    uncompressed-size and compression-ratio caps (``limits``; checked from the
    central directory before anything is decompressed), and orphan entries (a
    warning). Manifest checks: size/depth caps, JSON
    parseability, ``qvf_version`` / ``source``, section well-formedness, id
    uniqueness, reserved ids, kinds and extension-critical governance. Member
    checks: path rules, existence, a declared ``size`` against the central
    directory (before reading), sha256, JSON size/depth/parseability, dtype
    and binary byte-size, ``bands`` spin-axis shape agreement, the per-frame
    members of ``trajectory`` and ``reaction.path`` against ``coords`` (with
    the ``precision`` of ``coords``, the binary ``energies`` against
    ``metadata.energies`` and a trajectory ``lattice`` against
    ``metadata.pbc``), ``forces`` and ``hessian`` shapes against their
    structure, ``excited_states`` state indices, ``n_states`` and orbital
    rows, run.record UTF-8 and ``files`` index, reference resolution (a
    ``job.spec`` payload's ``structure_ref`` included), the ``structure``
    payload's periodicity
    contract, the shared array length of a 1-D ``spectrum`` payload, the
    ``meta`` entry of every ``atom_properties`` custom role, the atom
    indices of a standalone ``bonds`` section against the structure its
    ``structure_ref`` names, and the
    ``wavefunction.gto`` rules of spec § 5.3 (coefficient rank against
    ``complex``, ``symmetry_labels`` lengths, ``overlap`` and density-matrix
    shapes and pairing, ``basis.n_electrons`` against the occupations and
    ``provenance.n_electrons``, and the Appendix A.1 diagonal of ``overlap``,
    ``E_NORMALIZATION``), the data shapes of the vector, complex, series and
    Laplacian volumes, a series ``axis`` against its slabs, a periodic
    grid's span against the cell of its structure, and a grid
    ``units.data`` against the unit § 4.3 fixes for its volume kind
    (``E_UNIT_CONTRADICTED``) (spec § 5.2), and the
    ``properties.*`` rules (symmetric elastic matrices, ``static`` =
    ``electronic`` + ``ionic``, one Born tensor per atom, and a dipole that
    agrees with the root ``dipole_moment``), and the ``wavefunction.bloch``
    member shapes, k-point weights and periodic structure, the
    ``reaction.scheme`` participants, step ids and repeated structures, the
    bytes of ``preview`` members against their ``media_type``, a
    ``run.record`` ``job_ref``, and a ``structure`` section's ``charge`` and
    ``multiplicity`` against its atoms (``E_CHARGE_MULTIPLICITY``). Root
    checks: ``energies`` agrees with ``provenance`` and ``thermochemistry``,
    and the only ``structure`` section's charge and multiplicity with
    ``provenance``. The withdrawn ``qvf_version`` 2 is checked as 1 and
    reported as the warning ``W_VERSION_WITHDRAWN`` (spec § 3.1). Unit
    names: a declared unit outside the vocabulary of spec § 4.3 (a payload's
    ``units`` values, a 1-D spectrum's ``x_units`` / ``y_units`` /
    ``broadening.units``, the ``units`` of a ``{value, units}`` object in
    ``provenance``) is the warning ``W_UNIT_UNKNOWN``. Schema conformance — of
    the manifest (``E_SCHEMA``, which includes a ``byte_order`` other than
    ``"little"``) and of every JSON member with a payload definition
    (``E_PAYLOAD_SCHEMA``) — is checked when ``jsonschema`` is installed.
    """
    level = _level(level)
    lim = limits or Limits()
    out = _Findings()

    try:
        zf = _open_zip(path)
    except Exception as exc:  # noqa: BLE001 — any failure means unreadable
        reason = str(exc).removeprefix("E_ZIP_UNREADABLE: ").removeprefix("not a readable zip: ")
        out.error("E_ZIP_UNREADABLE", "archive", f"not a readable zip: {reason}")
        return out.report(0)
    with zf:
        n_sections = _validate_open_zip(zf, lim, out, level)
    return out.report(n_sections)


def _validate_open_zip(zf: zipfile.ZipFile, lim: Limits, out: _Findings,
                       level=Level.All) -> int:
    level = _level(level)
    payloads = level >= Level.Payloads
    semantics = level >= Level.Semantics
    out.references = level >= Level.References
    infos = zf.infolist()

    # --- container: central-directory checks, before decompressing anything ---
    counts = Counter(info.orig_filename for info in infos)
    for name, n in counts.items():
        if n > 1 and name == _MANIFEST:
            # An archive updated in place: the last manifest.json is the
            # manifest, every earlier one is superseded (spec § 2.3).
            for copy in range(1, n):
                out.warning("W_ORPHAN_MEMBER", f"entry:{name}",
                            f"ZIP entry {name!r} ({copy} of {n}) is superseded "
                            f"by the last {name!r} entry (archive updated in "
                            f"place; consumers ignore it)")
        elif n > 1:
            out.error("E_ZIP_DUPLICATE_ENTRY", f"entry:{name}",
                      f"ZIP entry name {name!r} occurs {n} times; readers "
                      f"disagree on which entry wins")
    for name in counts:
        why = _path_problem(name)
        if why:
            out.error("E_MEMBER_PATH", f"entry:{name}",
                      f"ZIP entry {name!r}: illegal path ({why})")

    unsafe: set[str] = set()   # entries that must not be decompressed
    total = 0
    for info in infos:
        name, usize, csize = info.orig_filename, info.file_size, info.compress_size
        total += usize
        if usize > lim.max_member_bytes:
            out.error("E_MEMBER_TOO_LARGE", f"entry:{name}",
                      f"ZIP entry {name!r} exceeds the per-member size cap "
                      f"({usize} > {lim.max_member_bytes} bytes uncompressed)")
            unsafe.add(name)
        if (lim.max_ratio > 0 and info.compress_type != zipfile.ZIP_STORED
                and usize >= lim.ratio_exempt_bytes
                and (csize == 0 or usize > lim.max_ratio * csize)):
            ratio = f"{usize / csize:.0f}" if csize else "inf"
            out.error("E_COMPRESSION_RATIO", f"entry:{name}",
                      f"ZIP entry {name!r} expands {usize} bytes from {csize} "
                      f"(ratio {ratio} > {lim.max_ratio:g})")
            unsafe.add(name)
    payloads_blocked = total > lim.max_total_bytes
    if payloads_blocked:
        out.error("E_TOTAL_TOO_LARGE", "archive",
                  f"archive expands to {total} bytes uncompressed, over the "
                  f"{lim.max_total_bytes}-byte total cap; member payloads were "
                  f"not read")

    # --- manifest ---
    if "manifest.json" not in counts:
        out.error("E_MANIFEST_MISSING", "manifest.json", "missing manifest.json")
        return 0
    # getinfo() resolves a duplicated name to its last central-directory
    # record: the current manifest of an archive updated in place (§ 2.3).
    if level == Level.Container:
        return 0
    minfo = zf.getinfo("manifest.json")
    if minfo.file_size > lim.max_json_bytes:
        out.error("E_JSON_TOO_LARGE", "manifest.json",
                  f"manifest.json is {minfo.file_size} bytes, over the "
                  f"{lim.max_json_bytes}-byte JSON cap")
        return 0
    if "manifest.json" in unsafe:
        return 0
    try:
        raw = _read_entry(zf, minfo)
    except Exception as exc:  # noqa: BLE001
        out.error("E_ZIP_UNREADABLE", "entry:manifest.json",
                  f"ZIP entry 'manifest.json' cannot be read: {_reason(exc)}")
        return 0
    if _json_depth(raw, lim.max_json_depth) > lim.max_json_depth:
        out.error("E_JSON_TOO_DEEP", "manifest.json",
                  f"manifest.json nests deeper than {lim.max_json_depth} levels")
        return 0
    try:
        manifest = _parse_json(raw)
    except Exception as exc:  # noqa: BLE001
        out.error("E_MANIFEST_JSON", "manifest.json",
                  f"manifest.json is not valid JSON: {exc}")
        return 0
    if not isinstance(manifest, dict):
        out.error("E_MANIFEST_JSON", "manifest.json",
                  "manifest.json must be a JSON object")
        return 0

    # Root sanity. The withdrawn version 2 is read as 1 (spec § 3.1).
    version = manifest.get("qvf_version")
    out.context.update(archive_format_version=version,
                       archive_schema_uri=manifest.get("schema_uri"),
                       unsupported_format_version=isinstance(version, bool) or version not in (1, 2))
    checked_manifest = manifest
    if version == 2 and not isinstance(version, bool):
        out.warning("W_VERSION_WITHDRAWN", "/qvf_version",
                    "qvf_version 2 was withdrawn: the archive is checked as "
                    "version 1, which producers must write (spec § 3.1)")
        checked_manifest = dict(manifest, qvf_version=1)
    elif isinstance(version, bool) or version != 1:
        out.error("E_VERSION", "/qvf_version",
                  f"unsupported archive format version {version!r}: "
                  f"qvf-python {out.context['reader_version']} selects the bundled "
                  "format-1 schema (withdrawn version 2 is read as 1); "
                  "use a reader supporting the producer's format, not a rewritten archive")
    src = manifest.get("source", {})
    if not isinstance(src, dict):
        out.error("E_SOURCE", "/source", "source must be a JSON object")
        src = {}
    for key in ("program", "version", "calculation"):
        if key not in src:
            out.error("E_SOURCE", "/source", f"source.{key} is missing")
    if not src.get("program") or not src.get("version"):
        out.error("E_SOURCE", "/source",
                  "source.program and source.version must be non-empty")
    if semantics:
        for where, msg in _root_consistency_errors(manifest):
            out.error("E_ROOT_INCONSISTENT", where, msg)

    # Optional full schema validation.
    validator = _load_schema_validator()
    out.context["schema_checked"] = validator is not None
    if validator is not None:
        for err in validator.iter_errors(checked_manifest):
            loc = "/".join(str(p) for p in err.absolute_path) or "(root)"
            out.error("E_SCHEMA", _pointer(*err.absolute_path) or "manifest.json",
                      f"{loc}: {err.message}")
    else:
        out.warning("W_SCHEMA_SKIPPED", "manifest.json",
                    "jsonschema not installed; structural schema check skipped")
    vocabulary = _unit_vocabulary()

    # Declared units of {value, units} objects in provenance (spec § 4.3).
    provenance = manifest.get("provenance")
    if semantics and vocabulary is not None and isinstance(provenance, dict):
        names, legacy = vocabulary
        for key, value in provenance.items():
            if not (isinstance(value, dict) and "value" in value
                    and "units" in value):
                continue
            unit = value["units"]
            if not (isinstance(unit, str) and (unit in names or unit in legacy)):
                out.warning("W_UNIT_UNKNOWN", _pointer("provenance", key, "units"),
                            f"provenance.{key}.units {unit!r} is not a unit "
                            f"name of spec § 4.3 or a legacy spelling")
    payload_validators: dict[str, Any] = {}

    def _payload_validator(kind: str, role: str):
        """jsonschema validator for the payload definition of (kind, role), or None."""
        defname = _PAYLOAD_DEFS.get((kind, role))
        if defname is None:
            return None
        if defname not in payload_validators:
            payload_validators[defname] = _load_payload_validator(defname)
        return payload_validators[defname]

    sections = manifest.get("sections", [])
    if not isinstance(sections, list):
        out.error("E_SECTION_MALFORMED", "/sections", "sections must be an array")
        sections = []
    section_ids = {s.get("id") for s in sections
                   if isinstance(s, dict) and isinstance(s.get("id"), str)}
    kinds_by_id: dict[str, str] = {}
    for s in sections:
        if (isinstance(s, dict) and isinstance(s.get("id"), str)
                and isinstance(s.get("kind"), str)):
            kinds_by_id.setdefault(s["id"], s["kind"])
    structure_ids = [s["id"] for s in sections
                     if isinstance(s, dict) and s.get("kind") == "structure"
                     and isinstance(s.get("id"), str)]
    only_structure = structure_ids[0] if len(structure_ids) == 1 else None
    n_atoms_of: dict[str, int] = {}    # structure id -> atom count
    sum_z_of: dict[str, int] = {}      # structure id -> sum of atomic numbers
    identity_sections: list[tuple] = []     # structure sections, checked last
    cells_of: dict[str, tuple[Any, Any]] = {}   # structure id -> (pbc, lattice)
    derivative_sections: list[tuple] = []   # forces / hessian, checked last
    bonds_sections: list[tuple] = []        # standalone bonds, checked last
    # Excited states (spec § 5.12): orbital rows of every wavefunction.gto,
    # state indices of every excited_states payload, and the sections whose
    # references are checked once every payload has been read.
    mo_rows_of: dict[str, dict[str, int]] = {}
    for s in sections:
        if (isinstance(s, dict) and s.get("kind") == "wavefunction.gto"
                and isinstance(s.get("id"), str) and s["id"] not in mo_rows_of):
            mo_rows_of[s["id"]] = _mo_rows(s.get("members"))
    state_indices_of: dict[str, Optional[set[int]]] = {}
    nto_sections: list[tuple] = []
    excited_sections: list[tuple] = []
    states_ref_spectra: list[tuple] = []
    transition_densities: list[tuple] = []
    grid_members: list[tuple] = []   # grid payloads, checked last (spec § 5.2)
    property_sections: list[tuple] = []     # properties.*, checked last
    periodic_of: dict[str, bool] = {}  # structure id -> any(pbc)
    bloch_sections: list[tuple] = []        # wavefunction.bloch, checked last
    # Section id -> kind, for reference resolution (spec § 4.1, § 5.5). With
    # duplicate ids (an error of its own) the first section with the id wins.
    section_kinds: dict[str, Any] = {}
    for s in sections:
        if isinstance(s, dict) and isinstance(s.get("id"), str):
            section_kinds.setdefault(s["id"], s.get("kind"))
    section_ids = set(section_kinds)
    declared_ext = manifest.get("extensions", {})
    if not isinstance(declared_ext, dict):
        declared_ext = {}
    declared_paths: set[str] = set()
    seen_ids: set[str] = set()

    for i, s in enumerate(sections):
        here = _pointer("sections", i)
        if not isinstance(s, dict):
            out.error("E_SECTION_MALFORMED", here,
                      f"section {i} must be a JSON object")
            continue
        sid = s.get("id")
        if not isinstance(sid, str):
            out.error("E_SECTION_MALFORMED", here,
                      f"section {i}: id must be a string, got {sid!r}")
            sid = f"<section {i}>"
        else:
            if sid in seen_ids:
                out.error("E_SECTION_ID_DUPLICATE", here + "/id",
                          f"duplicate section id {sid!r}")
            seen_ids.add(sid)
            if sid in _RESERVED_SECTION_IDS:
                out.error("E_RESERVED_ID", here + "/id",
                          f"section id {sid!r} is reserved (a viewer_defaults "
                          f"key, spec § 3.1, § 3.4)")

        kind = s.get("kind")
        if not isinstance(kind, str):
            out.error("E_SECTION_MALFORMED", here,
                      f"section {sid!r}: kind must be a string, got {kind!r}")
            kind = ""
        elif kind not in _CANONICAL_KINDS and not kind.startswith("x_"):
            out.error("E_KIND_UNKNOWN", here + "/kind",
                      f"section {sid!r}: unknown non-vendor kind {kind!r}")

        # critical governance (spec § 7.4).
        if s.get("critical") and kind.startswith("x_"):
            vendor = kind.split(".", 1)[0]
            if vendor not in declared_ext:
                out.error("E_CRITICAL_UNDECLARED", here + "/critical",
                          f"section {sid!r} is critical but its namespace "
                          f"{vendor!r} is not declared in root extensions")

        members = s.get("members")
        if not isinstance(members, dict):
            out.error("E_SECTION_MALFORMED", here,
                      f"section {sid!r}: members must be a JSON object")
            members = {}

        docs: dict[str, Any] = {}   # parsed JSON members of this section, by role
        overlap_raw: Optional[bytes] = None   # wavefunction.gto overlap bytes
        energies_raw: Optional[bytes] = None   # trajectory-family energies member
        for role, m in members.items():
            mwhere = _pointer("sections", i, "members", role)
            if (not isinstance(m, dict) or "path" not in m
                    or not isinstance(m["path"], str)):
                out.error("E_MEMBER_MALFORMED", mwhere,
                          f"section {sid!r} member {role!r}: malformed spec")
                continue
            mpath = m["path"]
            declared_paths.add(mpath)
            why = _path_problem(mpath)
            if why:
                out.error("E_MEMBER_PATH", mwhere,
                          f"member {mpath!r}: illegal path ({why})")
                continue
            if mpath not in counts:
                out.error("E_MEMBER_MISSING", mwhere,
                          f"member {mpath!r} declared but not in archive")
                continue
            info = zf.getinfo(mpath)
            # Declared size against the central directory, before reading.
            problem = _size_problem(m, info)
            if problem:
                out.error(problem[0], mwhere, problem[1])
            if not payloads or payloads_blocked or mpath in unsafe:
                continue   # the container finding already explains why
            fmt = m.get("format")
            if fmt == "json" and info.file_size > lim.max_json_bytes:
                out.error("E_JSON_TOO_LARGE", mwhere,
                          f"member {mpath!r}: {info.file_size} bytes, over the "
                          f"{lim.max_json_bytes}-byte JSON cap")
                continue
            try:
                raw = _read_entry(zf, info)
            except Exception as exc:  # noqa: BLE001
                out.error("E_ZIP_UNREADABLE", f"entry:{mpath}",
                          f"ZIP entry {mpath!r} cannot be read: {_reason(exc)}")
                continue
            if "sha256" in m and hashlib.sha256(raw).hexdigest() != m["sha256"]:
                out.error("E_SHA256_MISMATCH", mwhere,
                          f"member {mpath!r}: sha256 mismatch")
            if fmt == "json":
                if _json_depth(raw, lim.max_json_depth) > lim.max_json_depth:
                    out.error("E_JSON_TOO_DEEP", mwhere,
                              f"member {mpath!r}: JSON nests deeper than "
                              f"{lim.max_json_depth} levels")
                    continue
                try:
                    doc = _parse_json(raw)
                except Exception as exc:  # noqa: BLE001
                    out.error("E_JSON_INVALID", mwhere,
                              f"member {mpath!r}: not valid JSON ({exc})")
                else:
                    docs[role] = doc
                    if semantics and kind == "structure" and role == "structure":
                        for msg in _structure_payload_errors(sid, doc):
                            out.error("E_STRUCTURE_PERIODICITY", mwhere, msg)
                        if (isinstance(doc, dict)
                                and isinstance(doc.get("atoms"), list)):
                            n_atoms_of.setdefault(sid, len(doc["atoms"]))
                            zs = [a.get("atomic_number") if isinstance(a, dict)
                                  else None for a in doc["atoms"]]
                            if all(_is_int(z) and z >= 0 for z in zs):
                                sum_z_of.setdefault(sid, sum(zs))
                        if isinstance(doc, dict):
                            cells_of.setdefault(sid, (doc.get("pbc"),
                                                      doc.get("lattice_vectors")))
                            pbc = doc.get("pbc")
                            periodic_of.setdefault(sid, isinstance(pbc, list)
                                                   and any(x is True for x in pbc))
                    if (semantics and _PAYLOAD_DEFS.get((kind, role)) == "VolumeGridPayload"
                            and isinstance(doc, dict)):
                        grid_members.append((mwhere, sid, s, doc))
                    pv = _payload_validator(kind, role)
                    if pv is not None:
                        for err in pv.iter_errors(doc):
                            loc = "/".join(str(p) for p in err.absolute_path)
                            out.error("E_PAYLOAD_SCHEMA", mwhere,
                                      f"schema: {mpath}: {loc}: "
                                      f"{err.message}")
                    if semantics and vocabulary is not None and kind in _CANONICAL_KINDS:
                        for loc, unit in _declared_units(kind, doc, role):
                            if not (isinstance(unit, str)
                                    and unit in vocabulary[0]):
                                out.warning(
                                    "W_UNIT_UNKNOWN", mwhere,
                                    f"member {mpath!r}: declared unit "
                                    f"{unit!r} at {loc} is not a unit name "
                                    f"of spec § 4.3")
            elif fmt == "binary" and "dtype" in m:
                dtype = m["dtype"]
                itemsize = (_DTYPE_ITEMSIZE.get(dtype)
                            if isinstance(dtype, str) else None)
                shape = m.get("shape", [])
                if itemsize is None:
                    out.error("E_DTYPE", mwhere,
                              f"member {mpath!r}: bad dtype {dtype!r}")
                elif (not isinstance(shape, list) or not all(
                        isinstance(d, int) and not isinstance(d, bool) and d >= 0
                        for d in shape)):
                    out.error("E_MEMBER_MALFORMED", mwhere,
                              f"member {mpath!r}: shape must be a list of "
                              f"non-negative integers, got {shape!r}")
                else:
                    n = 1
                    for d in shape:
                        n *= d
                    if len(raw) != itemsize * n:
                        out.error("E_BINARY_SIZE", mwhere,
                                  f"member {mpath!r}: byte length {len(raw)} != "
                                  f"itemsize*prod(shape) {itemsize * n}")
                    elif kind == "wavefunction.gto" and role == "overlap":
                        overlap_raw = raw
            if (kind in _TRAJECTORY_KINDS and role == "energies"
                    and fmt == "binary"):
                energies_raw = raw
            if kind == "run.record" and role in ("input", "log"):
                try:
                    raw.decode("utf-8")
                except UnicodeDecodeError as exc:
                    out.error("E_UTF8", mwhere,
                              f"member {mpath!r}: run.record {role} must "
                              f"decode as UTF-8 ({exc})")
            if kind == "preview" and fmt == "binary":
                for code, msg in _preview_member_errors(sid, role, m, raw):
                    out.error(code, mwhere, msg)

        if not semantics:
            continue

        # 1-D spectra: one shared array length (spec § 5.4).
        if _PAYLOAD_DEFS.get((kind, "spectrum")) == "SpectrumPayload":
            for msg in _spectrum_payload_errors(sid, docs.get("spectrum")):
                out.error("E_PAYLOAD_SCHEMA",
                          _pointer("sections", i, "members", "spectrum"), msg)

        # Cross-references.
        if kind == "volume.difference":
            for op in ("operand_a", "operand_b"):
                ref = s.get(op)
                if ref is not None and (not isinstance(ref, str)
                                        or ref not in section_ids):
                    out.error("E_REF_UNRESOLVED", here + "/" + op,
                              f"section {sid!r}: {op} {ref!r} does not resolve")
        if kind == "reaction.waypoints":
            ref = s.get("trajectory_ref")
            if not isinstance(ref, str) or ref not in section_ids:
                out.error("E_REF_UNRESOLVED", here + "/trajectory_ref",
                          f"section {sid!r}: trajectory_ref {ref!r} does not "
                          f"resolve")
            elif section_kinds[ref] != "trajectory":
                out.error("E_REF_UNRESOLVED", here + "/trajectory_ref",
                          f"section {sid!r}: trajectory_ref {ref!r} names a "
                          f"section of kind {section_kinds[ref]!r}, not "
                          f"'trajectory' (spec § 5.5)")
        # Generic section references, on any kind (spec § 4.1).
        for field, target in _SECTION_REFS:
            ref = s.get(field)
            if ref is None:
                continue
            if not isinstance(ref, str) or ref not in section_ids:
                out.error("E_REF_UNRESOLVED", here + "/" + field,
                          f"section {sid!r}: {field} {ref!r} does not resolve")
            elif section_kinds[ref] != target:
                out.error("E_REF_UNRESOLVED", here + "/" + field,
                          f"section {sid!r}: {field} {ref!r} names a section "
                          f"of kind {section_kinds[ref]!r}, not {target!r} "
                          f"(spec § 4.1)")
        if kind == "bands":
            for role, msg in _bands_shape_errors(sid, members):
                out.error("E_MEMBER_MALFORMED",
                          _pointer("sections", i, "members", role), msg)
        if kind in _TRAJECTORY_KINDS:
            for role, msg in _trajectory_errors(sid, kind, members,
                                                docs.get("metadata"),
                                                energies_raw):
                out.error("E_MEMBER_MALFORMED",
                          _pointer("sections", i, "members", role), msg)
        if kind == "vibrations" and "hessian_ref" in s:
            ref = s["hessian_ref"]
            if not isinstance(ref, str) or kinds_by_id.get(ref) != "hessian":
                out.error("E_REF_UNRESOLVED", here + "/hessian_ref",
                          f"section {sid!r}: hessian_ref {ref!r} does not name "
                          f"a hessian section (spec § 5.5)")
        if kind == "job.spec" and isinstance(docs.get("spec"), dict) \
                and "structure_ref" in docs["spec"]:
            ref = docs["spec"]["structure_ref"]
            if not isinstance(ref, str) or kinds_by_id.get(ref) != "structure":
                out.error("E_REF_UNRESOLVED",
                          _pointer("sections", i, "members", "spec"),
                          f"section {sid!r}: spec.structure_ref {ref!r} does "
                          f"not name a structure section (spec § 5.9)")
        if kind in ("forces", "hessian"):
            derivative_sections.append((i, sid, kind, s, members,
                                        docs.get("meta")))
        if kind == "bonds" and "structure_ref" in s:
            bonds_sections.append((i, sid, s, docs.get("bonds")))
        if kind == "excited_states":
            excited_sections.append((i, sid, s, docs.get("states")))
            if sid not in state_indices_of:
                state_indices_of[sid] = _state_indices(docs.get("states"))
        if (_PAYLOAD_DEFS.get((kind, "spectrum")) == "SpectrumPayload"
                and isinstance(docs.get("spectrum"), dict)
                and "states_ref" in docs["spectrum"]):
            states_ref_spectra.append((i, sid, docs["spectrum"]))
        if kind == "volume.transition_density":
            transition_densities.append((i, sid, s))
        for code, role, msg in _volume_errors(sid, kind, members, docs,
                                              kinds_by_id):
            out.error(code, _pointer("sections", i, "members", role), msg)
        if kind in _PROPERTY_KINDS and "properties" in docs:
            property_sections.append((i, sid, kind, s, docs["properties"]))
        if kind == "wavefunction.bloch":
            bloch_sections.append((i, sid, s, members, docs))
        if kind == "atom_properties":
            for role, msg in _atom_properties_meta_errors(sid, members, docs):
                out.error("E_MEMBER_MALFORMED",
                          _pointer("sections", i, "members", role), msg)
        if kind == "wavefunction.gto":
            nto_sections.append((i, sid, docs.get("mo_metadata"), members))
            for msg in _unrestricted_nto_errors(docs.get("mo_metadata"), members):
                out.error("E_MEMBER_MALFORMED", _pointer("sections", i, "members", "mo_metadata"), msg)
            for code, role, msg in _wavefunction_errors(
                    sid, members, docs, manifest.get("provenance")):
                out.error(code, _pointer("sections", i, "members", role), msg)
            for msg in _normalization_errors(sid, docs.get("basis"),
                                             members.get("overlap"),
                                             overlap_raw):
                out.error("E_NORMALIZATION",
                          _pointer("sections", i, "members", "overlap"), msg)
        if kind == "reaction.scheme" and isinstance(docs.get("scheme"), dict):
            for code, msg in _reaction_scheme_errors(sid, docs["scheme"],
                                                     kinds_by_id):
                out.error(code, _pointer("sections", i, "members", "scheme"), msg)
        if kind == "run.record" and "job_ref" in s:
            ref = s["job_ref"]
            if not isinstance(ref, str) or kinds_by_id.get(ref) != "job.spec":
                out.error("E_REF_UNRESOLVED", here + "/job_ref",
                          f"section {sid!r}: job_ref {ref!r} does not name a "
                          f"job.spec section (spec § 5.8)")
        if kind == "structure":
            identity_sections.append((i, sid, s))
        if kind == "run.record" and isinstance(docs.get("files"), dict):
            for frole in docs["files"]:
                if frole not in members:
                    out.error("E_FILES_INDEX",
                              _pointer("sections", i, "members", "files"),
                              f"section {sid!r}: files index names role "
                              f"{frole!r} with no matching member")

    if semantics:
        # Chemical identity of the structures (spec § 5.1).
        for i, sid, s in identity_sections:
            for code, field, msg in _identity_errors(
                    sid, s, sum_z_of.get(sid), manifest.get("provenance"),
                    len(structure_ids) == 1):
                out.error(code, _pointer("sections", i, field), msg)

        # forces / hessian: shapes against the structure, which may come later.
        for i, sid, kind, s, members, meta in derivative_sections:
            for code, role, msg in _derivative_errors(
                    sid, kind, s, members, meta, n_atoms_of, kinds_by_id,
                    only_structure):
                out.error(code, _pointer("sections", i, "members", role), msg)

        # Standalone bonds: atom indices against the structure they name, which
        # may come later (spec § 5.1).
        for i, sid, s, payload in bonds_sections:
            for msg in _bonds_errors(sid, s, payload, kinds_by_id, n_atoms_of):
                out.error("E_MEMBER_MALFORMED",
                          _pointer("sections", i, "members", "bonds"), msg)

        for i, sid, meta, members in nto_sections:
            for code, msg in _nto_pair_errors(meta, members, nto_sections, kinds_by_id, state_indices_of):
                out.error(code, _pointer("sections", i, "members", "mo_metadata"), msg)

        # Excited states, the spectra that name them, and transition densities
        # (spec § 5.2, § 5.4, § 5.12): the sections they name may come later.
        for i, sid, s, payload in excited_sections:
            for code, field, msg in _excited_states_errors(
                    sid, s, payload, kinds_by_id, mo_rows_of):
                where = (_pointer("sections", i, field) if field
                         else _pointer("sections", i, "members", "states"))
                out.error(code, where, msg)
        for i, sid, payload in states_ref_spectra:
            for msg in _spectrum_states_errors(sid, payload, kinds_by_id,
                                               state_indices_of):
                out.error("E_REF_UNRESOLVED",
                          _pointer("sections", i, "members", "spectrum"), msg)
        for i, sid, s in transition_densities:
            for field, msg in _transition_density_errors(sid, s, kinds_by_id,
                                                         state_indices_of):
                out.error("E_REF_UNRESOLVED", _pointer("sections", i, field), msg)

        # Grids: references and periodic spans against the structure (spec § 5.2).
        for mwhere, sid, s, grid in grid_members:
            for code, msg in _grid_errors(sid, s, grid, kinds_by_id, cells_of,
                                          only_structure):
                out.error(code, mwhere, msg)
            msg = _grid_data_unit_error(sid, s.get("kind"), grid)
            if msg is not None:
                out.error("E_UNIT_CONTRADICTED", mwhere, msg)

        # properties.*: rules JSON Schema cannot express (spec § 5.13); the Born
        # charges are counted against a structure that may come later.
        for i, sid, kind, s, payload in property_sections:
            for code, msg in _property_errors(sid, kind, s, payload, manifest,
                                              n_atoms_of, only_structure):
                out.error(code, _pointer("sections", i, "members", "properties"),
                          msg)

        # wavefunction.bloch: shapes, and a periodic structure, which may come later.
        for i, sid, s, members, docs in bloch_sections:
            for code, where, msg in _bloch_errors(sid, s, members, docs,
                                                  kinds_by_id, periodic_of):
                out.error(code, _pointer("sections", i) + "/" + where, msg)

    # --- orphans: entries no manifest member references (spec § 2.1) ---
    for name in counts:
        if (name != "manifest.json" and not name.endswith("/")
                and name not in declared_paths):
            out.warning("W_ORPHAN_MEMBER", f"entry:{name}",
                        f"ZIP entry {name!r} is not referenced by any manifest "
                        f"member (orphan; consumers ignore it)")

    return len(sections)


def _main(argv: list[str]) -> int:
    args = argv[1:]
    as_json = "--json" in args
    files = [a for a in args if a != "--json"]
    if len(files) != 1 or any(a.startswith("-") for a in files):
        print("usage: python qvf_reader.py [--json] FILE.qvf", file=sys.stderr)
        return 2
    path = files[0]
    report = validate_qvf(path)
    if as_json:
        json.dump({"path": path, "ok": report["ok"],
                   "n_sections": report["n_sections"],
                   "findings": report["findings"],
                   "validation_context": report["validation_context"]}, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print(f"QVF validation: {path}")
        print(f"  {_context_summary(report)}")
        print(f"  sections: {report['n_sections']}")
        for f in report["findings"]:
            print(f"  {f['code']}: {f['message']}")
        print("  OK" if report["ok"] else "  FAILED")
    return 0 if report["ok"] else 1


def _main_entrypoint() -> int:
    """Console-script entry point (`qvf-validate [--json] FILE.qvf`)."""
    return _main(sys.argv)


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv))


def run(source, *, limits=None, level=Level.All):
    """Validate a path or open Reader with the shared finding codes."""
    if hasattr(source, "validate"):
        if limits is not None:
            raise TypeError("set Limits when opening Reader")
        return source.validate(level=level)
    return validate_qvf(source, limits=limits, level=level)
