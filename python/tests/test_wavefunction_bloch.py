"""The provisional ``wavefunction.bloch`` kind (spec § 5.3.1; P15, WP-30).

Covers the writer (round trip and refusals), the payload definitions, the
validator rules that JSON Schema cannot express (member shapes against each
other, the k-point weights, a periodic structure), and the Bloch-sum
convention of the conformance archive ``bloch.qvf``: its orbitals are
orthonormal in the overlap S(k) built with the phase e^{+ik.T} from the
basis it carries.
"""

import hashlib
import json
import math
import os
import sys
import zipfile

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))

from qvf_reader import QvfReader, validate_qvf  # noqa: E402
from qvf_writer import BOHR_PER_ANGSTROM, QvfWriter, primitive_norm  # noqa: E402

_SCHEMA_CODES = {"E_SCHEMA", "E_PAYLOAD_SCHEMA", "W_SCHEMA_SKIPPED"}
_CORPUS = os.path.normpath(os.path.join(HERE, "..", "..", "conformance", "corpus"))

CHAIN = [{"symbol": "H", "position": [0.0, 0.0, 0.0], "atomic_number": 1},
         {"symbol": "H", "position": [0.8, 0.0, 0.0], "atomic_number": 1}]
CELL = [[2.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
SHELLS = [{"center": c, "l": 0, "exponents": [0.4166], "coefficients": [1.0]}
          for c in (0, 1)]
KPTS = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]


def _codes(report):
    return {f["code"] for f in report["findings"]} - _SCHEMA_CODES


def _where(report, code):
    return sorted(f["where"] for f in report["findings"] if f["code"] == code)


def _orbitals(n_k=2, n_spin=2, n_mo=2, n_ao=2):
    rng = np.random.default_rng(7)
    coeff = (rng.standard_normal((n_k, n_spin, n_mo, n_ao))
             + 1j * rng.standard_normal((n_k, n_spin, n_mo, n_ao)))
    energies = np.sort(rng.standard_normal((n_k, n_spin, n_mo)), axis=-1)
    occ = np.zeros((n_k, n_spin, n_mo))
    occ[:, :, 0] = 1.0 if n_spin == 2 else 2.0
    return coeff, energies, occ


def _writer(periodic=True):
    w = QvfWriter("t", "1", "bloch", auto_metadata=False)
    if periodic:
        w.add_structure(CHAIN, pbc=[True, False, False], lattice_vectors=CELL)
    else:
        w.add_structure(CHAIN)
    return w


# --- writer -------------------------------------------------------------------

def test_writer_round_trip(tmp_path):
    coeff, eps, occ = _orbitals()
    w = _writer()
    sid = w.add_wavefunction_bloch(SHELLS, KPTS, coeff, eps, occ,
                                   weights=[0.5, 0.5], spin="unrestricted",
                                   labels=["G", "X"], mesh=[2, 1, 1],
                                   shift=[0, 0, 0], time_reversal=False,
                                   fermi_energy_ev=-3.5, n_mo_per_k=[2, 2])
    assert sid == "wf_bloch"
    path = str(tmp_path / "b.qvf")
    w.write(path)
    report = validate_qvf(path)
    assert report["ok"] and not report["warnings"], report["findings"]
    with QvfReader(path) as r:
        sec = next(s for s in r.sections if s["id"] == sid)
        assert sec["kind"] == "wavefunction.bloch"
        m = sec["members"]
        c = r.read_array(m["mo_coefficients"])
        assert c.shape == (2, 2, 2, 2, 2)
        np.testing.assert_array_equal(c[..., 0] + 1j * c[..., 1], coeff)
        np.testing.assert_array_equal(r.read_array(m["energies"]), eps)
        np.testing.assert_array_equal(r.read_array(m["occupations"]), occ)
        kp = r.read_json(m["kpoints"])
        assert kp == {"kpoints": KPTS, "weights": [0.5, 0.5], "labels": ["G", "X"],
                      "mesh": [2, 1, 1], "shift": [0.0, 0.0, 0.0],
                      "time_reversal": False}
        meta = r.read_json(m["mo_metadata"])
        assert meta == {"spin": "unrestricted", "n_spin": 2, "n_mo": 2, "n_ao": 2,
                        "orbital_kind": "canonical", "complex": True,
                        "k_point_convention": "fractional", "n_mo_per_k": [2, 2],
                        "fermi_energy": {"value": -3.5, "units": "ev"}}
        basis = r.read_json(m["basis"])
        assert basis["structure_ref"] == "structure" and basis["n_ao"] == 2


def test_writer_accepts_split_real_coefficients_and_divides_libint_norms(tmp_path):
    coeff, eps, occ = _orbitals(n_spin=1)
    split = np.stack([coeff.real, coeff.imag], axis=-1)
    w = _writer()
    w.add_wavefunction_bloch(SHELLS, KPTS, split, eps, occ,
                             coeffs_are_libint_normalized=True)
    path = str(tmp_path / "b.qvf")
    w.write(path)
    assert validate_qvf(path)["ok"]
    with QvfReader(path) as r:
        m = next(s for s in r.sections if s["id"] == "wf_bloch")["members"]
        np.testing.assert_array_equal(r.read_array(m["mo_coefficients"]), split)
        assert r.read_json(m["mo_metadata"])["spin"] == "restricted"
        assert "weights" not in r.read_json(m["kpoints"])
        c0 = r.read_json(m["basis"])["shells"][0]["coefficients"][0]
        assert c0 == pytest.approx(1.0 / primitive_norm(0.4166, 0))


@pytest.mark.parametrize("change, message", [
    ({"spin": "generalized"}, "spin"),
    ({"spin": "restricted"}, "n_spin"),
    ({"kpoints": [[0, 0]]}, "kpoints"),
    ({"energies": np.zeros((2, 2, 3))}, "energies"),
    ({"occupations": np.zeros((2, 1, 2))}, "occupations"),
    ({"mo_coefficients": np.zeros((2, 2, 2, 2, 1))}, "mo_coefficients"),
    ({"mo_coefficients": np.zeros((2, 2, 2, 2), dtype=complex)[..., :1]}, "AOs"),
    ({"weights": [0.5, 0.4]}, "weights"),
    ({"weights": [1.0]}, "weights"),
    ({"labels": ["G"]}, "labels"),
    ({"mesh": [2, 0, 1]}, "mesh"),
    ({"shift": [0, 0]}, "shift"),
    ({"n_mo_per_k": [2, 3]}, "n_mo_per_k"),
    ({"n_electrons": 1000, "weights": [0.5, 0.5]}, "n_electrons"),
    ({"n_electrons": -1}, "n_electrons"),
    ({"ecp": [{"center": 0}]}, "ecp"),
])
def test_writer_refusals(change, message):
    coeff, eps, occ = _orbitals()
    args = {"shells": SHELLS, "kpoints": KPTS, "mo_coefficients": coeff,
            "energies": eps, "occupations": occ, "spin": "unrestricted"}
    args.update(change)
    w = _writer()
    with pytest.raises(ValueError, match=message):
        w.add_wavefunction_bloch(**args)
    assert [s["id"] for s in w.build_manifest()["sections"]] == ["structure"]


# --- raw archives, so the validator sees what the writer would refuse ----------

def _json_bytes(obj):
    return json.dumps(obj).encode("utf-8")


class _Raw:
    def __init__(self, structure_first=True, periodic=True):
        self.manifest = {"qvf_version": 1,
                         "source": {"program": "t", "version": "1",
                                    "calculation": "bloch"},
                         "sections": []}
        self.files = {}
        self.periodic = periodic
        if structure_first:
            self.structure()

    def _member(self, path, data, dtype=None, shape=None):
        self.files[path] = data
        m = {"path": path, "format": "json" if dtype is None else "binary",
             "sha256": hashlib.sha256(data).hexdigest()}
        if dtype is not None:
            m["dtype"], m["shape"] = dtype, list(shape)
        return m

    def structure(self, sid="structure", periodic=None):
        periodic = self.periodic if periodic is None else periodic
        payload = ({"atoms": CHAIN, "pbc": [True, False, False],
                    "lattice_vectors": CELL} if periodic
                   else {"atoms": CHAIN, "pbc": [False] * 3,
                         "lattice_vectors": None})
        self.manifest["sections"].append({"id": sid, "kind": "structure", "members": {
            "structure": self._member(f"{sid}/structure.json", _json_bytes(payload))}})
        return self

    def bloch(self, *, coeff=(2, 2, 2, 2, 2), energies=(2, 2, 2), occ=(2, 2, 2),
              kpoints=None, meta=None, basis=None, fields=None):
        kp = {"kpoints": KPTS, "weights": [0.5, 0.5]} if kpoints is None else kpoints
        mt = ({"spin": "unrestricted", "n_spin": 2, "n_mo": 2, "n_ao": 2,
               "complex": True, "k_point_convention": "fractional"}
              if meta is None else meta)
        bs = ({"structure_ref": "structure", "pure": True, "n_ao": 2,
               "shells": SHELLS} if basis is None else basis)
        members = {
            "basis": self._member("wf/basis.json", _json_bytes(bs)),
            "kpoints": self._member("wf/kpoints.json", _json_bytes(kp)),
            "mo_metadata": self._member("wf/mo_metadata.json", _json_bytes(mt)),
        }
        for role, shape in (("mo_coefficients", coeff), ("energies", energies),
                            ("occupations", occ)):
            members[role] = self._member(f"wf/{role}.bin",
                                         np.zeros(shape, "<f8").tobytes(),
                                         "float64", shape)
        self.manifest["sections"].append({"id": "wf", "kind": "wavefunction.bloch",
                                          **(fields or {}), "members": members})
        return self

    def validate(self, tmp_path):
        path = str(tmp_path / "raw.qvf")
        with zipfile.ZipFile(path, "w") as zf:
            for p, data in self.files.items():
                zf.writestr(p, data)
            zf.writestr("manifest.json", _json_bytes(self.manifest))
        return validate_qvf(path)


def test_valid_raw_archive_is_clean(tmp_path):
    report = _Raw().bloch().validate(tmp_path)
    assert report["ok"] and not report["findings"], report["findings"]


_META = {"spin": "unrestricted", "n_spin": 2, "n_mo": 2, "n_ao": 2,
         "complex": True, "k_point_convention": "fractional"}


@pytest.mark.parametrize("kwargs, where", [
    ({"coeff": (2, 2, 2, 2, 1)}, "mo_coefficients"),
    ({"coeff": (2, 3, 2, 2, 2), "energies": (2, 3, 2), "occ": (2, 3, 2)},
     "mo_coefficients"),
    ({"coeff": (2, 2, 2, 2)}, "mo_coefficients"),
    ({"energies": (2, 2, 3)}, "energies"),
    ({"energies": (2, 2)}, "energies"),
    ({"occ": (1, 2, 2)}, "occupations"),
    ({"kpoints": {"kpoints": KPTS + [[0.25, 0, 0]], "weights": [0.5, 0.25, 0.25]}},
     "kpoints"),
    ({"kpoints": {"kpoints": KPTS, "weights": [1.0]}}, "kpoints"),
    ({"kpoints": {"kpoints": KPTS, "weights": [0.5, 0.25]}}, "kpoints"),
    ({"kpoints": {"kpoints": KPTS, "weights": [0.5, 0.5], "labels": ["G"]}},
     "kpoints"),
    ({"meta": {**_META, "n_mo": 3}}, "mo_metadata"),
    ({"meta": {**_META, "n_ao": 3}}, "mo_metadata"),
    ({"meta": {**_META, "spin": "restricted"}}, "mo_metadata"),
    ({"meta": {**_META, "n_mo_per_k": [2]}}, "mo_metadata"),
    ({"meta": {**_META, "n_mo_per_k": [2, 3]}}, "mo_metadata"),
    ({"basis": {"structure_ref": "structure", "n_ao": 3, "shells": SHELLS}},
     "basis"),
])
def test_shape_rules(tmp_path, kwargs, where):
    report = _Raw().bloch(**kwargs).validate(tmp_path)
    assert _codes(report) == {"E_MEMBER_MALFORMED"}, report["findings"]
    assert f"/sections/1/members/{where}" in _where(report, "E_MEMBER_MALFORMED")


def test_weights_within_tolerance_and_absent_weights_are_fine(tmp_path):
    ok = _Raw().bloch(kpoints={"kpoints": KPTS, "weights": [0.5, 0.5 + 5e-7]})
    assert ok.validate(tmp_path)["ok"]
    assert _Raw().bloch(kpoints={"kpoints": KPTS}).validate(tmp_path)["ok"]


def test_molecular_structure_is_rejected(tmp_path):
    report = _Raw(periodic=False).bloch().validate(tmp_path)
    assert _codes(report) == {"E_MEMBER_MALFORMED"}
    assert _where(report, "E_MEMBER_MALFORMED") == ["/sections/1/members/basis"]


def test_structure_after_the_section_is_found(tmp_path):
    raw = _Raw(structure_first=False).bloch().structure()
    assert raw.validate(tmp_path)["ok"]
    raw = _Raw(structure_first=False, periodic=False).bloch().structure()
    assert _codes(raw.validate(tmp_path)) == {"E_MEMBER_MALFORMED"}


def test_basis_structure_ref_must_name_a_structure(tmp_path):
    basis = {"structure_ref": "nowhere", "n_ao": 2, "shells": SHELLS}
    report = _Raw().bloch(basis=basis).validate(tmp_path)
    assert _codes(report) == {"E_REF_UNRESOLVED"}
    assert _where(report, "E_REF_UNRESOLVED") == ["/sections/1/members/basis"]
    basis = {"structure_ref": "wf", "n_ao": 2, "shells": SHELLS}
    assert _codes(_Raw().bloch(basis=basis).validate(tmp_path)) == {"E_REF_UNRESOLVED"}


def test_section_level_structure_ref_takes_precedence(tmp_path):
    # The section names the periodic structure; the basis names a molecular one.
    raw = _Raw().structure("mol", periodic=False)
    basis = {"structure_ref": "mol", "n_ao": 2, "shells": SHELLS}
    raw.bloch(basis=basis, fields={"structure_ref": "structure"})
    assert raw.validate(tmp_path)["ok"]
    # The section names the molecular structure: reported at structure_ref.
    raw = _Raw().structure("mol", periodic=False)
    raw.bloch(fields={"structure_ref": "mol"})
    report = raw.validate(tmp_path)
    assert _codes(report) == {"E_MEMBER_MALFORMED"}
    assert _where(report, "E_MEMBER_MALFORMED") == ["/sections/2/structure_ref"]
    # An unresolved section-level reference is the generic finding only.
    raw = _Raw().bloch(fields={"structure_ref": "nowhere"})
    report = raw.validate(tmp_path)
    assert _codes(report) == {"E_REF_UNRESOLVED"}
    assert _where(report, "E_REF_UNRESOLVED") == ["/sections/1/structure_ref"]


# --- payload definitions --------------------------------------------------------

jsonschema = pytest.importorskip("jsonschema")


def _schema_codes(report):
    return {f["code"] for f in report["findings"]}


@pytest.mark.parametrize("kwargs", [
    {"kpoints": {"kpoints": KPTS, "extra": 1}},
    {"kpoints": {"weights": [0.5, 0.5]}},
    {"kpoints": {"kpoints": KPTS, "weights": [1.5, -0.5]}},
    {"kpoints": {"kpoints": KPTS, "units": {"kpoints": "bohr"}}},
    {"meta": {**_META, "complex": False}},
    {"meta": {**_META, "spin": "restricted", "n_spin": 2}},
    {"meta": {**_META, "k_point_convention": "cartesian"}},
    {"meta": {k: v for k, v in _META.items() if k != "n_ao"}},
    {"meta": {**_META, "energies": [0.0, 1.0]}},
])
def test_payload_definitions(tmp_path, kwargs):
    report = _Raw().bloch(**kwargs).validate(tmp_path)
    assert "E_PAYLOAD_SCHEMA" in _schema_codes(report), report["findings"]


def test_payloads_with_x_keys_and_units_are_accepted(tmp_path):
    kp = {"kpoints": KPTS, "weights": [0.5, 0.5], "x_vendor": True,
          "units": {"kpoints": "dimensionless", "weights": "dimensionless"}}
    meta = {**_META, "x_vendor_restart_kind": "closed_shell",
            "units": {"energies": "hartree"},
            "fermi_energy": {"value": -0.2, "units": "hartree"}}
    report = _Raw().bloch(kpoints=kp, meta=meta).validate(tmp_path)
    assert report["ok"] and not report["findings"], report["findings"]


# --- the Bloch-sum convention of the conformance archive ------------------------

def test_corpus_orbitals_are_orthonormal_in_the_stated_bloch_metric():
    """spec § 5.3.1: S_μν(k) = Σ_T e^{+ik·T} ⟨χ_μ(r)|χ_ν(r − T)⟩ with k·T =
    2π k_frac·n; the corpus orbitals satisfy C†S(k)C = 1 with S built from
    the archive's own structure and basis, and the bands section repeats the
    orbital energies in eV."""
    path = os.path.join(_CORPUS, "bloch.qvf")
    if not os.path.exists(path):
        pytest.skip("conformance corpus not present")
    with QvfReader(path) as r:
        by_id = {s["id"]: s for s in r.sections}
        st = r.read_json(by_id["structure"]["members"]["structure"])
        m = by_id["wf_bloch"]["members"]
        basis, kp = r.read_json(m["basis"]), r.read_json(m["kpoints"])
        c = r.read_array(m["mo_coefficients"])
        eps = r.read_array(m["energies"])
        eig = r.read_array(by_id["bands"]["members"]["eigenvalues"])
    a = np.array(st["lattice_vectors"][0]) * BOHR_PER_ANGSTROM
    centers = [np.array(st["atoms"][sh["center"]]["position"]) * BOHR_PER_ANGSTROM
               for sh in basis["shells"]]
    alphas = [sh["exponents"][0] for sh in basis["shells"]]
    assert math.isclose(sum(kp["weights"]), 1.0)
    for ik, k in enumerate(kp["kpoints"]):
        s_k = np.zeros((2, 2), dtype=complex)
        for mu in range(2):
            for nu in range(2):
                for n in range(-12, 13):
                    r2 = float(np.sum((centers[nu] + n * a - centers[mu]) ** 2))
                    p, q = alphas[mu], alphas[nu]
                    ovl = (2 * math.sqrt(p * q) / (p + q)) ** 1.5 \
                        * math.exp(-p * q * r2 / (p + q))
                    s_k[mu, nu] += np.exp(2j * math.pi * k[0] * n) * ovl
        for spin in range(2):
            cmat = (c[ik, spin, :, :, 0] + 1j * c[ik, spin, :, :, 1]).T
            np.testing.assert_allclose(cmat.conj().T @ s_k @ cmat, np.eye(2),
                                       atol=1e-10)
    assert np.abs(c[0, :, :, :, 1]).max() == 0.0          # real at Γ
    assert np.abs(c[1, :, :, :, 1]).max() > 0.1           # complex at X
    np.testing.assert_allclose(eig, np.transpose(eps, (1, 0, 2)) * 27.211386245988,
                               rtol=0, atol=1e-9)
