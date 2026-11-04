#!/usr/bin/env python3
"""Turn .qvf archives into seed inputs for cpp/fuzz/manifest_fuzzer.

Each seed is one flags byte, the archive's manifest.json, then every other
entry as ``\\n=qvf=\\n<name>\\n<bytes>`` (the input format described in
cpp/fuzz/manifest_fuzzer.cpp). Archives Python's zipfile cannot open, or
whose entries do not decompress, are skipped: zip_reader_fuzzer covers those.

Usage::

    python scripts/fuzz_manifest_seeds.py OUT_DIR ARCHIVE_OR_DIR [...]
"""
import hashlib
import pathlib
import sys
import zipfile

SEP = b"\n=qvf=\n"
FLAGS = 1  # rewrite member sha256 values, STORE entries, schema check on


def seed_of(path: pathlib.Path) -> bytes | None:
    try:
        with zipfile.ZipFile(path) as z:
            names = z.namelist()
            if "manifest.json" not in names:
                return None
            out = bytes([FLAGS]) + z.read("manifest.json")
            for name in names:
                if name != "manifest.json" and not name.endswith("/"):
                    out += SEP + name.encode() + b"\n" + z.read(name)
            return out
    except Exception:  # noqa: BLE001 (any unreadable archive is skipped)
        return None


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    out_dir = pathlib.Path(argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for arg in argv[2:]:
        p = pathlib.Path(arg)
        for f in sorted(p.rglob("*.qvf")) if p.is_dir() else [p]:
            seed = seed_of(f)
            if seed is not None:
                (out_dir / hashlib.sha1(seed).hexdigest()).write_bytes(seed)
                n += 1
    print(f"{n} seeds written to {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
