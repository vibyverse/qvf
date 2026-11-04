"""Parity: the generated schema checker (C++ and JavaScript) ≡ jsonschema.

The C++ and JavaScript validators check ``manifest.json`` and every JSON member
payload with the table-driven checker generated from the schema
(``scripts/gen_kinds.py``); the Python reference uses ``jsonschema``. These
tests drive the built C++ tool (``qvf_schema_check``, built by CMake beside
``run_conformance``) and the JavaScript twin (``js/scripts/schema-check.mjs``
over ``js/dist``) and require:

* every positive corpus archive validates with the same findings as in Python
  (no errors);
* every negative and legacy archive gives the same finding codes as Python,
  the schema-only codes included, and every ``E_SCHEMA`` / ``E_PAYLOAD_SCHEMA``
  finding sits at the same location (JSON pointer);
* with the schema check disabled, the same codes as Python without
  ``jsonschema`` (``W_SCHEMA_SKIPPED``);
* on thousands of documents mutated from the corpus manifests and payloads,
  the same multiset of error locations as ``jsonschema``.

Where the binaries live: ``QVF_CPP_BUILD_DIR`` (the CMake build directory) and
``QVF_NODE`` (the node executable); otherwise a few conventional build
directories and ``node`` on PATH are tried, and a missing implementation is
skipped.
"""

import copy
import glob
import json
import os
import random
import shutil
import subprocess
import sys
import zipfile

import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, ".."))

import qvf_reader  # noqa: E402
from qvf_reader import Limits, validate_qvf  # noqa: E402

jsonschema = pytest.importorskip("jsonschema")

CORPUS = os.path.join(ROOT, "conformance", "corpus")
NEGATIVE = os.path.join(ROOT, "conformance", "negative")
LEGACY = os.path.join(ROOT, "conformance", "legacy")
SCHEMA_CODES = {"E_SCHEMA", "E_PAYLOAD_SCHEMA"}


# --- the implementations ---------------------------------------------------------

def _cpp_tool():
    names = ["qvf_schema_check", "qvf_schema_check.exe"]
    dirs = []
    if os.environ.get("QVF_CPP_BUILD_DIR"):
        dirs.append(os.environ["QVF_CPP_BUILD_DIR"])
    dirs += [os.path.join(ROOT, d) for d in ("build", "build-cpp", "cpp/build", "build/san")]
    for d in dirs:
        for sub in ("", "Release", "Debug"):
            for n in names:
                path = os.path.join(d, sub, n)
                if os.path.isfile(path) and os.access(path, os.X_OK):
                    return [path]
    return None


def _js_tool():
    node = os.environ.get("QVF_NODE") or shutil.which("node")
    script = os.path.join(ROOT, "js", "scripts", "schema-check.mjs")
    dist = os.path.join(ROOT, "js", "dist", "schemaCheck.js")
    if not node or not os.path.exists(script) or not os.path.exists(dist):
        return None
    return [node, script]


TOOLS = {"cpp": _cpp_tool, "js": _js_tool}


@pytest.fixture(params=["cpp", "js"])
def impl(request):
    cmd = TOOLS[request.param]()
    if cmd is None:
        pytest.skip(f"the {request.param} implementation is not built here "
                    f"(set QVF_CPP_BUILD_DIR / build js/dist)")
    return request.param, cmd


def _run(cmd, requests):
    """One process for many requests; one JSON reply per request."""
    text = "".join(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n" for r in requests)
    proc = subprocess.run(cmd, input=text.encode("utf-8"), capture_output=True,
                          check=False, timeout=600, cwd=ROOT)
    assert proc.returncode == 0, proc.stderr.decode("utf-8", "replace")
    replies = [json.loads(line) for line in proc.stdout.decode("utf-8").splitlines() if line]
    assert len(replies) == len(requests)
    for req, rep in zip(requests, replies):
        assert not (isinstance(rep, dict) and "error" in rep), (req, rep)
    return replies


# --- archives -----------------------------------------------------------------------

def _expectations(directory):
    out = []
    for ef in sorted(glob.glob(os.path.join(directory, "*.expected.json"))):
        with open(ef, encoding="utf-8") as fh:
            exp = json.load(fh)
        out.append((os.path.join(directory, exp["archive"]), exp.get("limits")))
    return out


def _python_findings(path, limits, schema=True, monkeypatch=None):
    if not schema:
        monkeypatch.setattr(qvf_reader, "_load_schema_validator", lambda: None)
        monkeypatch.setattr(qvf_reader, "_load_payload_validator", lambda _d: None)
    try:
        return validate_qvf(path, limits=Limits.from_dict(limits))["findings"]
    finally:
        if not schema:
            monkeypatch.undo()


def _compare_archives(impl, archives, schema=True, monkeypatch=None):
    name, cmd = impl
    requests = [{"archive": p, "limits": lim, "schema_check": schema} for p, lim in archives]
    problems = []
    for (path, limits), got in zip(archives, _run(cmd, requests)):
        archive = os.path.basename(path)
        got_codes = {f["code"] for f in got}
        want = _python_findings(path, limits, schema, monkeypatch)
        want_codes = {f["code"] for f in want}
        if got_codes != want_codes:
            problems.append(f"{archive}: codes {sorted(got_codes)} != Python {sorted(want_codes)}")
        got_at = {(f["code"], f["where"]) for f in got if f["code"] in SCHEMA_CODES}
        want_at = {(f["code"], f["where"]) for f in want if f["code"] in SCHEMA_CODES}
        if got_at != want_at:
            problems.append(f"{archive}: schema findings {sorted(got_at)} != Python "
                            f"{sorted(want_at)}")
        if any(f["severity"] == "error" for f in got) != any(f["severity"] == "error"
                                                             for f in want):
            problems.append(f"{archive}: verdict differs from Python")
    assert not problems, "\n".join(problems)


def test_positive_corpus_matches_python(impl):
    archives = _expectations(CORPUS)
    assert len(archives) >= 40
    _compare_archives(impl, archives)


def test_negative_and_legacy_corpora_match_python(impl):
    archives = _expectations(NEGATIVE) + _expectations(LEGACY)
    assert len(archives) >= 60
    _compare_archives(impl, archives)


def test_without_the_schema_check_matches_python_without_jsonschema(impl, monkeypatch):
    archives = _expectations(CORPUS) + _expectations(NEGATIVE) + _expectations(LEGACY)
    _compare_archives(impl, archives, schema=False, monkeypatch=monkeypatch)


# --- documents --------------------------------------------------------------------------

def _pointer(path):
    return "".join("/" + str(p).replace("~", "~0").replace("/", "~1") for p in path)


def _validator(defname):
    return (qvf_reader._load_schema_validator() if defname is None
            else qvf_reader._load_payload_validator(defname))


def _base_documents():
    """(def, document) pairs: every corpus manifest (def None) and every JSON
    member with a payload definition."""
    docs = []
    for path in sorted(glob.glob(os.path.join(CORPUS, "*.qvf"))):
        if os.path.basename(path) == "zip64.qvf":
            continue
        with zipfile.ZipFile(path) as zf:
            manifest = json.loads(zf.read("manifest.json"))
            docs.append((None, manifest))
            for section in manifest.get("sections", []):
                for role, member in section.get("members", {}).items():
                    defname = qvf_reader._PAYLOAD_DEFS.get((section["kind"], role))
                    if defname and member.get("format") == "json":
                        docs.append((defname, json.loads(zf.read(member["path"]))))
    return docs


_STRINGS = ["", "a b", "x_v.k", "x_", "custom.flux", "attachment.log", "../x", "a/../b",
            "/abs", "ev", "EV", "Eh", "bohr", "hartree", "float16", "float64", "json",
            "binary", "little", "big", "unrestricted", "restricted", "ÅÅÅÅ", "😀",
            "2026-13-01T00:00:00Z", "2026-01-02T03:04:05Z", "a" * 64, "A" * 64,
            "3f2504e0-4f89-41d3-9a0c-0305e82c3301", "structure", "bonds", "volume.density",
            "spectra.ir", "trajectory", "C:/x", "x" * 300]
_SCALARS = [None, True, False, 0, 1, -1, 2, 3, 4, 1.0, 2.5, -1e-9, 1e300, 230, 231]
_KEYS = ["x_extra", "bogus", "units", "kind", "spin", "pbc", "lattice_vectors",
         "dimensionality", "members", "id", "size", "dtype", "shape", "custom.q",
         "attachment.a", "operand_a", "metadata", "mass_weighted", "alpha"]


def _nodes(doc, path=()):
    yield path, doc
    if isinstance(doc, dict):
        for k, v in doc.items():
            yield from _nodes(v, path + (k,))
    elif isinstance(doc, list):
        for i, v in enumerate(doc):
            yield from _nodes(v, path + (i,))


def _set(doc, path, value):
    if not path:
        return value
    parent = doc
    for p in path[:-1]:
        parent = parent[p]
    parent[path[-1]] = value
    return doc


def _random_value(rnd, depth=0):
    roll = rnd.random()
    if roll < 0.45:
        return rnd.choice(_SCALARS)
    if roll < 0.8 or depth > 1:
        return rnd.choice(_STRINGS)
    if roll < 0.9:
        return [_random_value(rnd, depth + 1) for _ in range(rnd.randint(0, 4))]
    return {rnd.choice(_KEYS): _random_value(rnd, depth + 1) for _ in range(rnd.randint(0, 3))}


def _mutate(doc, rnd):
    nodes = list(_nodes(doc))
    path, node = rnd.choice(nodes)
    roll = rnd.random()
    if isinstance(node, dict) and node and roll < 0.3:
        del node[rnd.choice(list(node))]
    elif isinstance(node, dict) and roll < 0.55:
        node[rnd.choice(_KEYS + list(node))] = _random_value(rnd)
    elif isinstance(node, list) and node and roll < 0.3:
        del node[rnd.randrange(len(node))]
    elif isinstance(node, list) and roll < 0.55:
        node.append(copy.deepcopy(rnd.choice(node)) if node and rnd.random() < 0.6
                    else _random_value(rnd))
    elif isinstance(node, (int, float)) and not isinstance(node, bool) and roll < 0.8:
        doc = _set(doc, path, rnd.choice([-node, node + 1, float(node), 0, 2.5, -1]))
    else:
        doc = _set(doc, path, _random_value(rnd))
    return doc


def _mutated_documents(n_per_doc=6, seed=20261009):
    rnd = random.Random(seed)
    out = []
    for defname, base in _base_documents():
        out.append((defname, base))
        for _ in range(n_per_doc):
            doc = copy.deepcopy(base)
            for _ in range(rnd.randint(1, 3)):
                doc = _mutate(doc, rnd)
            out.append((defname, doc))
    return out


def test_mutated_documents_match_jsonschema(impl):
    _, cmd = impl
    docs = _mutated_documents()
    assert len(docs) > 1000
    replies = _run(cmd, [{"def": d, "doc": doc} for d, doc in docs])
    validators = {}
    mismatches = []
    n_errors = 0
    for (defname, doc), got in zip(docs, replies):
        if defname not in validators:
            validators[defname] = _validator(defname)
        want = sorted(_pointer(e.absolute_path) for e in validators[defname].iter_errors(doc))
        n_errors += len(want)
        if sorted(p for p, _ in got) != want:
            mismatches.append((defname, json.dumps(doc)[:300], sorted(p for p, _ in got), want))
    assert n_errors > 1000, "the mutations should break the documents often"
    assert not mismatches, "\n".join(map(str, mismatches[:10]))


def test_targeted_documents_match_jsonschema(impl):
    """Documents that exercise each keyword of the supported subset."""
    _, cmd = impl
    hexa = "a" * 64
    lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    docs = [
        ("Sha256Hex", hexa), ("Sha256Hex", hexa.upper()), ("Sha256Hex", 5),
        ("ZipPath", ""), ("ZipPath", "a/../b"), ("ZipPath", "a/..b"), ("ZipPath", "/x"),
        ("MemberSize", 2.0), ("MemberSize", 2.5), ("MemberSize", True), ("MemberSize", -3),
        ("ValueWithUnits", {"value": 1, "units": "Eh"}), ("ValueWithUnits", {"value": "1"}),
        ("StructureAtom", {"symbol": "😀😀😀", "position": [0, 0, 0], "atomic_number": 0}),
        ("StructureAtom", {"symbol": "😀😀😀😀", "position": [0, 0], "atomic_number": -1}),
        ("StructurePayload", {"atoms": [], "pbc": [True, False, False], "dimensionality": 2,
                              "lattice_vectors": lattice}),
        ("StructurePayload", {"atoms": [], "pbc": [True, True, True], "dimensionality": 0,
                              "lattice_vectors": None}),
        ("MoMetadataPayload", {"spin": "unrestricted", "energies": [1], "occupations": [1]}),
        ("MoMetadataPayload", {"spin": "restricted", "alpha": {}}),
        ("MoMetadataPayload", {"spin": "restricted", "energies": [], "occupations": [],
                               "coefficient_components": ["real", "real", "x"]}),
        ("HessianMetaPayload", {"mass_weighted": True}),
        ("AtomPropertiesMetaPayload", {"custom.q": {"units": "e"}, "other": 1}),
        ("RunRecordFilesPayload", {"input": {"filename": ""}, "attachment.a": {"x": 1},
                                   "log": 3}),
        ("Section", 42), ("Section", {"id": "a b", "kind": "x_v.k", "members": {}}),
        ("Section", {"id": "s", "kind": "atom_properties", "members": {"meta": {}}}),
        ("Section", {"id": "s", "kind": "volume.difference", "operand_a": "x",
                     "members": {}}),
        ("SectionAtomProperties", {"members": {"custom.q": {}}}),
        ("SectionAtomProperties", {"members": {"meta": {"path": "m", "format": "json",
                                                        "sha256": hexa}}}),
        ("Camera", {"position": [0, 0, 1], "focal_point": [0, 0, 0], "view_up": [0, 1, 0],
                    "view_angle": 10, "parallel_scale": 1}),
        ("Camera", {"position": [0, 0, 1], "focal_point": [0, 0, 0], "view_up": [0, 1, 0]}),
        ("TrajectoryLattice", {"path": "p", "format": "binary", "dtype": "float64",
                               "shape": [3, 3, 3], "sha256": hexa}),
        ("TrajectoryLattice", {"path": "p", "format": "binary", "dtype": "float64",
                               "shape": [0, 3, 3], "sha256": hexa}),
        ("Extensions", {"x_a": {"version": ""}, "y": {}}),
        ("EnergyComponents", {"x_a/b~c": {"value": "no", "units": "ev"}}),
        (None, {}), (None, []), (None, {"qvf_version": 1.0, "source": {"program": "p",
                                        "version": "1", "calculation": "c"},
                                        "sections": [], "byte_order": "big"}),
    ]
    replies = _run(cmd, [{"def": d, "doc": doc} for d, doc in docs])
    for (defname, doc), got in zip(docs, replies):
        want = sorted(_pointer(e.absolute_path) for e in _validator(defname).iter_errors(doc))
        assert sorted(p for p, _ in got) == want, (defname, doc, got)
