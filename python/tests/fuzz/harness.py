"""The invariants every fuzz input is checked against.

`check(raw)` runs the validator (`qvf.validate.validate_qvf` on an in-memory
stream) and the readers (`qvf.reader.Reader` and the compatibility
`qvf.validate.QvfReader`) over one input and returns a list of problems,
empty when the input was handled correctly:

* the validator never raises, and every finding carries an Appendix B code of
  the matching severity; a rejected input has at least one error;
* a reader refuses an input only with a `ValueError` whose message starts with
  an Appendix B code, or with the API errors its documentation names
  (`KeyError` for an absent section, `TypeError` for a member of the wrong
  format, `zipfile.BadZipFile` for bytes that are not a ZIP);
* no call takes longer than `TIME_LIMIT` seconds;
* the peak Python allocation stays within what the input justifies: a fixed
  base plus a small multiple of the bytes its entries can really inflate to
  (the declared size, capped by the DEFLATE limit of 1032:1).
"""
from __future__ import annotations

import io
import re
import signal
import struct
import sys
import time
import traceback
import tracemalloc
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from qvf import validate as V  # noqa: E402
from qvf.reader import Reader  # noqa: E402

TIME_LIMIT = 10.0
MEMORY_BASE = 48 * 1024 ** 2
MEMORY_FACTOR = 4
CODES = dict(V.FINDING_CODES)
_CODE = re.compile(r"^(E_[A-Z0-9_]+|W_[A-Z0-9_]+)\b")

# Small caps keep each input cheap (as in the C++ fuzz target), and the
# default caps check the reference values against lying headers.
SMALL = V.Limits(max_member_bytes=1 << 24, max_total_bytes=1 << 26,
                 max_json_bytes=1 << 22)
DEFAULT = V.Limits()


class Timeout(BaseException):
    """Raised by the interval timer; a BaseException so `except Exception` cannot eat it."""


def _alarm(signum, frame):
    raise Timeout()


def _inflatable(raw: bytes) -> int:
    """Upper bound on the bytes the entries of `raw` can legitimately inflate to.

    Read from the central directory with a tolerant scan: for every record, the
    smaller of the declared size and 1032 times the compressed size (the DEFLATE
    expansion limit), plus the stored bytes themselves.
    """
    total = len(raw)
    pos = raw.find(b"PK\x01\x02")
    while pos != -1 and pos + 46 <= len(raw):
        csize, usize = struct.unpack_from("<II", raw, pos + 20)
        if usize == 0xFFFFFFFF or csize == 0xFFFFFFFF:
            usize = min(usize, 1 << 34)
        total += min(usize, 1032 * min(csize, len(raw)) + 65536)
        pos = raw.find(b"PK\x01\x02", pos + 4)
    return total


def _signature(exc: BaseException) -> str:
    tb = traceback.extract_tb(exc.__traceback__)
    frames = [f for f in tb if "/qvf/" in f.filename or "zipfile" in f.filename] or tb
    last = frames[-1] if frames else None
    where = f"{Path(last.filename).name}:{last.lineno}" if last else "?"
    return f"{type(exc).__name__} at {where}"


def check_validator(raw: bytes, limits) -> list[str]:
    problems = []
    try:
        report = V.validate_qvf(io.BytesIO(raw), limits=limits)
    except Timeout:
        raise
    except BaseException as exc:  # noqa: BLE001
        return [f"validator raised {_signature(exc)}: {exc!r:.200}"]
    errors = 0
    for f in report["findings"]:
        sev = CODES.get(f["code"])
        if sev is None:
            problems.append(f"validator: undocumented code {f['code']}")
        elif sev != f["severity"]:
            problems.append(f"validator: {f['code']} reported as {f['severity']}")
        errors += f["severity"] == "error"
    if report["ok"] != (errors == 0):
        problems.append("validator: ok disagrees with the findings")
    return problems


def _allowed(exc: BaseException) -> bool:
    if isinstance(exc, ValueError):
        if str(exc).startswith("tolist(): "):
            return True     # the documented cap of Tensor.tolist
        m = _CODE.match(str(exc))
        return bool(m) and m.group(1) in CODES
    if isinstance(exc, zipfile.BadZipFile):
        return True
    if isinstance(exc, KeyError):
        return True
    if isinstance(exc, TypeError):
        return bool(re.search(r"is not JSON|is not a tensor|expected .*got", str(exc)))
    return False


def check_reader(raw: bytes, limits) -> list[str]:
    problems = []

    def guard(label, fn):
        try:
            return fn()
        except Timeout:
            raise
        except BaseException as exc:  # noqa: BLE001
            if not _allowed(exc):
                problems.append(f"{label} raised {_signature(exc)}: {exc!r:.200}")
            return None

    reader = guard("Reader()", lambda: Reader.from_bytes(raw, limits=limits))
    if reader is not None:
        with reader:
            for section in reader.sections:
                for role, member in (section.members or {}).items():
                    if member.format == "json":
                        guard("Reader.json", lambda: reader.json(section, role))
                    elif member.format == "binary" and member.dtype is not None:
                        t = guard("Reader.tensor", lambda: reader.tensor(section, role))
                        if t is not None and len(t.data) <= 1 << 20:
                            guard("Tensor.tolist", t.tolist)
                    else:
                        guard("Reader.read_member", lambda: reader.read_member(member))
                if isinstance(section.id, str):
                    guard("Reader.get", lambda: reader.get(section.id))
            guard("Reader.validate", reader.validate)

    legacy = guard("QvfReader()", lambda: V.QvfReader(io.BytesIO(raw), limits=limits))
    if legacy is not None:
        with legacy:
            sections = legacy.sections if isinstance(legacy.manifest, dict) else []
            for s in sections if isinstance(sections, list) else []:
                members = s.get("members") if isinstance(s, dict) else None
                for m in (members.values() if isinstance(members, dict) else []):
                    if isinstance(m, dict) and isinstance(m.get("path"), str):
                        if m.get("format") == "json":
                            guard("QvfReader.read_json", lambda: legacy.read_json(m))
                        else:
                            guard("QvfReader.read_member", lambda: legacy.read_member(m))
    return problems


def check(raw: bytes, *, limits=SMALL, reader: bool = True,
          memory: bool = True) -> tuple[list[str], float, int]:
    """(problems, seconds, peak bytes) for one input."""
    problems = []
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, TIME_LIMIT)
    if memory:
        tracemalloc.start()
        tracemalloc.reset_peak()
    t0 = time.perf_counter()
    try:
        problems += check_validator(raw, limits)
        if reader:
            problems += check_reader(raw, limits)
    except Timeout:
        problems.append(f"timeout: over {TIME_LIMIT:g} s")
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)
        elapsed = time.perf_counter() - t0
        peak = 0
        if memory:
            peak = tracemalloc.get_traced_memory()[1]
            tracemalloc.stop()
    bound = MEMORY_BASE + MEMORY_FACTOR * _inflatable(raw)
    if memory and peak > bound:
        problems.append(f"memory: peak {peak} bytes over the bound {bound}")
    return problems, elapsed, peak
