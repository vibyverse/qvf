"""Drift guard and unit tests of the tables generated from the schema.

``scripts/gen_kinds.py`` reads ``spec/qvf_manifest.schema.json``,
``registry.json`` and spec Appendix B, and writes
``cpp/include/qvf/kinds_generated.hpp``, ``cpp/src/qvf_schema_table.hpp``,
``js/src/schemaTable.ts`` and ``python/qvf_kinds.py``. These tests regenerate
them and compare with the committed files (a stale table fails here; run
``python3 scripts/gen_kinds.py``), check the content against the schema, and
pin what the generator refuses.
"""

import copy
import functools
import importlib.util
import json
import os
import random
import re
import subprocess
import sys

import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
_SCRIPT = os.path.join(ROOT, "scripts", "gen_kinds.py")
sys.path.insert(0, os.path.join(HERE, ".."))

import qvf_kinds  # noqa: E402
import qvf_reader  # noqa: E402


@functools.lru_cache(maxsize=None)
def _load():
    spec = importlib.util.spec_from_file_location("gen_kinds", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _gen():
    if not os.path.exists(_SCRIPT):
        pytest.skip("scripts/gen_kinds.py not present (standalone checkout)")
    return _load()


def _inputs():
    with open(os.path.join(ROOT, "spec", "qvf_manifest.schema.json"), encoding="utf-8") as fh:
        schema = json.load(fh)
    with open(os.path.join(ROOT, "registry.json"), encoding="utf-8") as fh:
        registry = json.load(fh)
    with open(os.path.join(ROOT, "spec", "qvf-format-spec.md"), encoding="utf-8") as fh:
        spec = fh.read()
    return schema, registry, spec


# --- drift ---------------------------------------------------------------------

def test_generated_files_are_current():
    gen = _gen()
    outputs = gen.build()
    assert set(outputs) == {gen.OUT_CPP_KINDS, gen.OUT_CPP_TABLE, gen.OUT_TS_TABLE,
                            gen.OUT_PY_KINDS, gen.OUT_CPP_TYPED, gen.OUT_C_KINDS,
                            gen.OUT_KIND_DOCS, gen.OUT_RUST_TABLE, gen.OUT_RUST_KINDS}
    for rel, text in outputs.items():
        path = os.path.join(ROOT, rel)
        assert os.path.exists(path), f"{rel} missing; run python3 scripts/gen_kinds.py"
        with open(path, encoding="utf-8", newline="") as fh:
            assert fh.read() == text, f"{rel} is stale; run python3 scripts/gen_kinds.py"
        assert gen.GENERATED_NOTE in text.splitlines()[0], rel


def test_kind_catalogue_covers_every_kind():
    """docs/kinds.md has one entry per canonical kind, in registry order, each
    with a corpus example, and every in-page link names an anchor."""
    gen = _gen()
    _, registry, _ = _inputs()
    text = gen.build()[gen.OUT_KIND_DOCS]
    headings = re.findall(r"^## `([^`]+)`$", text, re.M)
    assert headings == registry["canonical_kinds"]
    entries = re.split(r"^## `", text, flags=re.M)[1:]
    for kind, entry in zip(headings, entries):
        assert "**Example** (from `conformance/corpus/" in entry, kind
    for defname in set(gen.PAYLOAD_DEFS.values()):
        assert f"### `{defname}`" in text, defname
    anchors = set(re.findall(r'<a id="([^"]+)"></a>', text))
    assert set(re.findall(r"\]\(#([^)]+)\)", text)) <= anchors
    for target in re.findall(r"\]\((\.\./[^)#]+)\)", text):
        assert os.path.exists(os.path.join(ROOT, "docs", target)), target


def test_generation_is_deterministic():
    gen = _gen()
    assert gen.build() == gen.build()


def test_check_option_reports_up_to_date():
    _gen()
    run = subprocess.run([sys.executable, _SCRIPT, "--check"], capture_output=True,
                         text=True, check=False)
    assert run.returncode == 0, run.stderr
    assert "up to date" in run.stdout


# --- content ----------------------------------------------------------------------

def test_kinds_follow_registry_and_schema():
    schema, registry, _ = _inputs()
    assert list(qvf_kinds.CANONICAL_KINDS) == registry["canonical_kinds"]
    assert qvf_kinds.SCHEMA_ID == schema["$id"] == registry["schema_id"]
    assert list(qvf_kinds.READABLE_VERSIONS) == registry["readable_versions"]
    defs = schema["$defs"]
    for kind in qvf_kinds.CANONICAL_KINDS:
        info = qvf_kinds.kind_info(kind)
        branch = next(defs[b["$ref"].rsplit("/", 1)[1]] for b in defs["Section"]["oneOf"]
                      if defs[b["$ref"].rsplit("/", 1)[1]].get("properties", {})
                      .get("kind", {}).get("const") == kind)
        members = branch["properties"]["members"]
        named = [r.role for r in info.roles if not r.pattern]
        assert named == list(members.get("properties", {})), kind
        assert {r.role for r in info.roles if r.required} == set(members.get("required", []))
        assert [r.role for r in info.roles if r.pattern] == \
            list(members.get("patternProperties", {}))
        assert info.allow_extra_members == (members.get("additionalProperties", True)
                                            is not False)
    assert qvf_kinds.kind_info("x_vendor.thing") is None


def test_role_details():
    wf = qvf_kinds.kind_info("wavefunction.gto")
    mo = wf.role("mo_coefficients")
    assert (mo.format, mo.dtypes, mo.min_rank, mo.max_rank) == ("binary", ("float64",), 2, 3)
    assert wf.role("basis").payload_def == "BasisPayload"
    data = qvf_kinds.kind_info("volume.density").role("data")
    assert data.dtypes == ("float32", "float64") and data.min_rank == data.max_rank == 3
    lattice = qvf_kinds.kind_info("trajectory").role("lattice")
    assert (lattice.min_rank, lattice.max_rank) == (2, 3)    # oneOf [3,3] | [n,3,3]
    energies = qvf_kinds.kind_info("dos.total").role("energies")
    assert energies.dtypes == qvf_kinds.DTYPES and energies.max_rank is None
    rr = qvf_kinds.kind_info("run.record")
    assert rr.role("log").format == "opaque" and rr.role("log").dtypes == ()
    assert [r.role for r in rr.roles if r.pattern] == [r"^attachment\.[A-Za-z0-9_.-]+$"]
    assert qvf_kinds.kind_info("citations").role("references").format == "opaque"


def test_payload_defs_cover_every_payload_definition():
    schema, _, _ = _inputs()
    defs = schema["$defs"]
    payloads = {n for n in defs if n.endswith("Payload")}
    assert set(qvf_kinds.PAYLOAD_DEFS.values()) == payloads
    assert dict(qvf_kinds.PAYLOAD_DEFS) == qvf_reader._PAYLOAD_DEFS
    for (kind, role), defname in qvf_kinds.PAYLOAD_DEFS.items():
        r = qvf_kinds.kind_info(kind).role(role)
        assert r.format == "json" and r.payload_def == defname, (kind, role)


def test_units_and_finding_codes():
    schema, _, spec = _inputs()
    defs = schema["$defs"]
    assert list(qvf_kinds.UNITS) == defs["UnitName"]["enum"]
    assert list(qvf_kinds.LEGACY_UNITS) == defs["LegacyUnitName"]["enum"]
    assert set(qvf_kinds.LEGACY_CANONICAL.values()) <= set(qvf_kinds.UNITS)
    assert list(qvf_kinds.SCHEMA_DEFS) == list(defs)
    assert qvf_kinds.FIXED_UNITS["VolumeGridPayload"] == {"origin": "bohr",
                                                          "voxel_vectors": "bohr"}
    appendix = spec.split("## Appendix B.", 1)[1]
    rows = dict(re.findall(r"^\| `([EW]_[A-Z0-9_]+)` \| (error|warning) \|", appendix, re.M))
    assert dict(qvf_kinds.FINDING_CODES) == rows == qvf_reader.FINDING_CODES


def test_cpp_unit_enumerators_are_stable():
    """Existing C++ enumerator names never change (qvf::Unit is public API)."""
    path = os.path.join(ROOT, "cpp", "include", "qvf", "kinds_generated.hpp")
    if not os.path.exists(path):
        pytest.skip("C++ sources not present")
    with open(path, encoding="utf-8") as fh:
        header = fh.read()
    enum = header[header.index("enum class Unit {"):]
    enum = enum[:enum.index("};")]
    names = dict((u, e) for e, u in re.findall(r"^\s+(\w+),\s*//\s*(\S+)\s*$", enum, re.M))
    assert names["ev"] == "ElectronVolt" and names["cm-1"] == "Wavenumber"
    assert names["kj_per_mol"] == "KjPerMol" and names["e"] == "ElementaryCharge"
    assert names["dimensionless"] == "Dimensionless" and names["a3"] == "Angstrom3"


# --- what the generator refuses ---------------------------------------------------------

def _model(schema, registry, spec):
    gen = _gen()
    gen.SchemaModel(schema, registry, spec)
    gen.SchemaCompiler(schema)


def test_refuses_unsupported_keywords_and_refs():
    gen = _gen()
    schema, registry, spec = _inputs()
    for patch in ({"unevaluatedProperties": False}, {"$ref": "other.json#/x"},
                  {"$dynamicRef": "#x"}, {"items": [{}]}):
        bad = copy.deepcopy(schema)
        bad["$defs"]["Vec3"] = {**bad["$defs"]["Vec3"], **patch}
        with pytest.raises(gen.GenerationError):
            _model(bad, registry, spec)


def test_refuses_an_unmapped_payload_definition():
    gen = _gen()
    schema, registry, spec = _inputs()
    bad = copy.deepcopy(schema)
    bad["$defs"]["NewThingPayload"] = {"type": "object"}
    with pytest.raises(gen.GenerationError, match="PAYLOAD_DEFS"):
        _model(bad, registry, spec)


def test_refuses_registry_drift_and_unmapped_legacy_units():
    gen = _gen()
    schema, registry, spec = _inputs()
    reg = copy.deepcopy(registry)
    reg["canonical_kinds"] = reg["canonical_kinds"][1:]
    with pytest.raises(gen.GenerationError, match="canonical_kinds"):
        _model(schema, reg, spec)
    bad = copy.deepcopy(schema)
    bad["$defs"]["LegacyUnitName"]["enum"].append("kcal")
    with pytest.raises(gen.GenerationError, match="LEGACY_CANONICAL"):
        _model(bad, registry, spec)


@pytest.mark.parametrize("pattern", [r"(?=a)", r"\1", r"\bword", r"(?<n>x)", r"a{2,1}",
                                     r"[z-a]", r"*a", r"\p{L}", "(a", "a)"])
def test_refuses_patterns_outside_the_subset(pattern):
    gen = _gen()
    with pytest.raises(gen.GenerationError):
        gen.compile_regex(pattern)


# --- the compiled patterns ---------------------------------------------------------------

def _pike(prog, text):
    """A reference run of a compiled program (what both interpreters do)."""
    cps = [ord(c) for c in text]
    matched = False

    def add(lst, pc, pos, seen):
        nonlocal matched
        stack = [pc]
        while stack:
            pc = stack.pop()
            if pc in seen:
                continue
            seen.add(pc)
            op, x, y = prog[pc]
            if op == "Jump":
                stack.append(x)
            elif op == "Split":
                stack.extend([y, x])
            elif op == "LineStart":
                if pos == 0:
                    stack.append(pc + 1)
            elif op == "LineEnd":
                if pos == len(cps):
                    stack.append(pc + 1)
            elif op == "Match":
                matched = True
            else:
                lst.append(pc)

    cur = []
    add(cur, 0, 0, set())
    for pos, c in enumerate(cps):
        if matched:
            break
        nxt, seen = [], set()
        for pc in cur:
            op, ranges, _ = prog[pc]
            hit = any(lo <= c <= hi for lo, hi in ranges)
            if hit == (op == "Class"):
                add(nxt, pc + 1, pos + 1, seen)
        add(nxt, 0, pos + 1, seen)
        cur = nxt
    return matched


def _schema_patterns():
    schema, _, _ = _inputs()
    out = set()

    def walk(node):
        if isinstance(node, dict):
            for k, v in node.items():
                if k == "pattern" and isinstance(v, str):
                    out.add(v)
                elif k == "patternProperties" and isinstance(v, dict):
                    out.update(v)
                if k not in ("enum", "const", "examples", "default"):
                    walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)

    walk(schema)
    return sorted(out)


def test_compiled_patterns_agree_with_re():
    """Every schema pattern, plus constructs the schema may grow into, gives
    the same answer as Python's re on strings built from its own alphabet
    (no trailing newline: there ECMA-262 and Python's `$` differ)."""
    gen = _gen()
    extra = [r"a|b|", r"^(ab)*c?$", r"x{2,}", r"[^\d\s]+", r"\w-\W", r"^.$", r"[\]a-]",
             r"é{1,3}", r"(?:ab|cd){2}$", r"\x41\t"]
    rnd = random.Random(41)
    alphabet = list("abcdx_-./0123456789:TZ\\ AéÅ\té") + ["custom.", "attachment.", "x_", ".."]
    for pattern in _schema_patterns() + extra:
        prog = gen.compile_regex(pattern)
        rx = re.compile(pattern, re.ASCII)   # ECMA-262 \d \w \s are ASCII
        samples = ["", "2026-01-02T03:04:05Z", "2026-06-30T23:59:60.25Z", "a" * 64, "x_a.b",
                   "custom.flux", "attachment.geom", "input", "a/../b", "/abs",
                   "3f2504e0-4f89-41d3-9a0c-0305e82c3301"]
        samples += ["".join(rnd.choice(alphabet) for _ in range(rnd.randint(0, 12)))
                    for _ in range(300)]
        for text in samples:
            if text.endswith("\n"):
                continue
            want = rx.search(text) is not None
            assert _pike(prog, text) == want, (pattern, text)


def test_c_entries_cover_every_canonical_kind():
    gen = _gen()
    _, registry, _ = _inputs()
    entries = re.findall(r'QVF_C_KIND\((\w+), "([^"]+)"\)', gen.build()[gen.OUT_C_KINDS])
    assert [kind for _, kind in entries] == registry["canonical_kinds"]
    assert len({name for name, _ in entries}) == len(entries)
    assert ("structure_payload", "structure") in entries
