"""Energy decomposition analysis (spec § 5.16): the ``energy_decomposition`` kind, its
reference writer, the rules validators check, and the ``volume.difference``
fields that name a decomposition's orbital pairs (§ 5.2)."""
from __future__ import annotations

import copy
import hashlib
import json
import zipfile

import pytest

from qvf import validate_qvf
from qvf._kinds import PAYLOAD_DEFS
from qvf.writer import QvfWriter

np = pytest.importorskip("numpy")

WATER_DIMER = [{"symbol": s, "position": p} for s, p in [
    ("O", [0.0, 0.0, 0.0]), ("H", [0.0, 0.0, 0.96]), ("H", [0.0, 0.929, -0.24]),
    ("O", [0.0, 0.0, 2.95]), ("H", [0.0, 0.759, 3.538]), ("H", [0.0, -0.759, 3.538])]]

ZR = {
    "scheme": "ziegler-rauk",
    "variant": "eda-nocv",
    "method": {"method": "rks", "functional": "pbe", "basis": "tz2p", "dispersion": "d3bj"},
    "fragments": [{"label": "donor", "atoms": [0, 1, 2], "charge": 0, "multiplicity": 1,
                   "basis": "own"},
                  {"label": "acceptor", "atoms": [3, 4, 5], "basis": "own"}],
    "terms": [{"name": "interaction", "energy": -0.0080},
              {"name": "electrostatic", "energy": -0.0125},
              {"name": "pauli", "energy": 0.0132},
              {"name": "orbital", "energy": -0.0073},
              {"name": "dispersion", "energy": {"value": -0.0014 * 627.5094740631,
                                                "units": "kcal_per_mol"}}],
    "sum_rules": [{"total": 0, "parts": [1, 2, 3, 4]}],
    "orbital_pairs": {"type": "nocv", "eigenvalue_convention": "per_spin",
                      "pairs": [{"index": 1, "eigenvalue": 0.045, "energy": -0.0049,
                                 "density_ref": "nocv_1"},
                                {"index": 2, "eigenvalue": 0.017, "energy": -0.0011}]},
}


def _archive(tmp_path, payload=ZR, *, densities=True, name="energy_decomposition.qvf"):
    w = QvfWriter("test", "1.0", auto_metadata=False)
    w.add_structure(WATER_DIMER)
    if densities:
        w.add_volume("volume.difference", np.zeros((3, 3, 3)), [0.0, 0.0, 0.0],
                     np.eye(3) * 0.5, section_id="nocv_1", decomposition_ref="energy_decomposition",
                     decomposition_role="orbital_pair", pair_index=1)
    w.add_energy_decomposition(payload)
    path = tmp_path / name
    w.write(str(path))
    return path


def _codes(path) -> set[str]:
    return {f["code"] for f in validate_qvf(str(path))["findings"]}


def _rewrite(src, dst, *, payload_edit=None, section_edit=None, sid="energy_decomposition"):
    """Copy ``src`` to ``dst`` with the energy_decomposition payload or a section edited and the
    member's sha256 and size updated."""
    with zipfile.ZipFile(src) as zin:
        manifest = json.loads(zin.read("manifest.json"))
        replace = {}
        for sec in manifest["sections"]:
            if sec["id"] != sid:
                continue
            if section_edit is not None:
                section_edit(sec)
            if payload_edit is not None:
                member = sec["members"]["decomposition"]
                payload = json.loads(zin.read(member["path"]))
                payload_edit(payload)
                data = json.dumps(payload).encode()
                member["sha256"] = hashlib.sha256(data).hexdigest()
                member["size"] = len(data)
                replace[member["path"]] = data
        with zipfile.ZipFile(dst, "w") as zout:
            for info in zin.infolist():
                if info.filename != "manifest.json":
                    zout.writestr(info.filename, replace.get(info.filename,
                                                             zin.read(info.filename)))
            zout.writestr("manifest.json", json.dumps(manifest))
    return dst


def test_payload_definition():
    assert PAYLOAD_DEFS[("energy_decomposition", "decomposition")] == "EnergyDecompositionPayload"


def test_ziegler_rauk_nocv_round_trip(tmp_path):
    path = _archive(tmp_path)
    report = validate_qvf(str(path))
    assert report["ok"], report["findings"]
    with zipfile.ZipFile(path) as z:
        manifest = json.loads(z.read("manifest.json"))
        payload = json.loads(z.read("energy_decomposition/decomposition.json"))
    sec = {s["id"]: s for s in manifest["sections"]}
    assert sec["energy_decomposition"]["kind"] == "energy_decomposition"
    assert sec["nocv_1"]["decomposition_ref"] == "energy_decomposition" and sec["nocv_1"]["pair_index"] == 1
    # numbers are written in hartree, {value, units} objects as given
    assert payload["terms"][0]["energy"] == {"value": -0.008, "units": "hartree"}
    assert payload["terms"][4]["energy"]["units"] == "kcal_per_mol"


SAPT = {
    "scheme": "sapt", "variant": "sapt0",
    "fragments": [{"label": "A", "atoms": [0, 1, 2]}, {"label": "B", "atoms": [3, 4, 5]},
                  {"label": "A: O", "atoms": [0], "parent": 0},
                  {"label": "B: O", "atoms": [3], "parent": 1}],
    "terms": [{"name": "interaction", "energy": -0.0078},
              {"name": "electrostatic", "energy": -0.0134},
              {"name": "exchange", "energy": 0.0136},
              {"name": "induction", "energy": -0.0039},
              {"name": "dispersion", "energy": -0.0041},
              {"name": "x_ind20r", "parent": 3, "order": "(20),r", "energy": -0.0045},
              {"name": "exchange_induction", "parent": 3, "energy": 0.0024},
              {"name": "delta_hf", "parent": 3, "energy": -0.0018},
              {"name": "electrostatic", "fragments": [2, 3], "energy": -0.010}],
    "sum_rules": [{"total": 0, "parts": [1, 2, 3, 4]}, {"total": 3, "parts": [5, 6, 7]}],
}

LED = {
    "scheme": "led", "variant": "dlpno-ccsd(t)",
    "fragments": [{"label": f"W{k}", "atoms": [k]} for k in range(3)],
    "terms": [{"name": "interaction", "energy": -0.03},
              {"name": "preparation", "fragments": [0], "energy": 0.002},
              {"name": "electrostatic", "fragments": [0, 1], "energy": -0.010},
              {"name": "electrostatic", "fragments": [1, 2], "energy": -0.011},
              {"name": "dispersion", "fragments": [2, 0], "energy": -0.011}],
    "sum_rules": [{"total": 0, "parts": [1, 2, 3, 4]}],
}

ALMO = {
    "scheme": "almo-eda",
    "fragments": [{"label": "A", "atoms": [0, 1, 2]}, {"label": "B", "atoms": [3, 4, 5]}],
    "terms": [{"name": "interaction", "energy": -0.010},
              {"name": "frozen", "energy": -0.004},
              {"name": "polarization", "energy": -0.002},
              {"name": "charge_transfer", "energy": -0.004},
              {"name": "charge_transfer", "fragments": [0, 1], "directed": True, "energy": -0.003},
              {"name": "charge_transfer", "fragments": [1, 0], "directed": True, "energy": -0.001}],
    "sum_rules": [{"total": 0, "parts": [1, 2, 3]}, {"total": 3, "parts": [4, 5]}],
    "orbital_pairs": {"type": "covp", "pairs": [
        {"index": 1, "eigenvalue": 0.01, "energy": -0.003, "donor": 0, "acceptor": 1}]},
}


@pytest.mark.parametrize("payload", [SAPT, LED, ALMO], ids=["sapt", "led", "almo"])
def test_other_schemes_valid(tmp_path, payload):
    path = _archive(tmp_path, payload, densities=False)
    report = validate_qvf(str(path))
    assert report["ok"], report["findings"]


def test_periodic_per_k_pairs(tmp_path):
    w = QvfWriter("test", "1.0", auto_metadata=False)
    w.add_structure(WATER_DIMER, lattice_vectors=[[6.0, 0, 0], [0, 6.0, 0], [0, 0, 20.0]],
                    pbc=[True, True, False])
    payload = copy.deepcopy(ZR)
    payload["variant"] = "peda"
    payload["orbital_pairs"]["pairs"] = [
        {"index": 1, "eigenvalue": 0.04, "energy": -0.004,
         "kpoint_eigenvalues": [0.05, 0.03], "kpoint_energies": [-0.0025, -0.0015]}]
    payload["periodic"] = {"kpoints": {"kpoints": [[0, 0, 0], [0.5, 0, 0]],
                                       "weights": [0.5, 0.5]}, "mesh": [2, 1, 1]}
    w.add_energy_decomposition(payload)
    path = tmp_path / "p.qvf"
    w.write(str(path))
    report = validate_qvf(str(path))
    assert report["ok"], report["findings"]
    bad = _rewrite(path, tmp_path / "pk.qvf", payload_edit=lambda p: p["orbital_pairs"]
                   ["pairs"][0].__setitem__("kpoint_energies", [-0.004]))
    assert "E_MEMBER_MALFORMED" in _codes(bad)


def _bad(**changes):
    payload = copy.deepcopy(ZR)
    for key, value in changes.items():
        payload[key] = value
    return payload


@pytest.mark.parametrize("payload, match", [
    (_bad(fragments=[{"label": "A", "atoms": [0, 1, 2]}]), "at least two"),
    (_bad(fragments=[{"label": "A", "atoms": [0, 1, 2]}, {"label": "B", "atoms": [2, 3]}]),
     "disjoint"),
    (_bad(fragments=ZR["fragments"] + [{"label": "s", "atoms": [4], "parent": 0}]),
     "does not contain"),
    (_bad(terms=ZR["terms"] + [{"name": "interaction", "energy": -0.008}]), "repeats"),
    (_bad(terms=ZR["terms"][1:]), "exactly one"),
    (_bad(terms=ZR["terms"] + [{"name": "x_a", "parent": 9, "energy": 0.0}]),
     "not the index"),
    (_bad(terms=ZR["terms"] + [{"name": "pauli", "fragments": [5], "energy": 0.0}]),
     "fragment 5"),
    (_bad(sum_rules=[{"total": 0, "parts": [1, 2, 3]}]), "add up"),
    (_bad(sum_rules=[{"total": 0, "parts": [1, 7]}]), "names term"),
    (_bad(terms=ZR["terms"][:4] + [{"name": "dispersion",
                                    "energy": {"value": 1.0, "units": "angstrom"}}]),
     "not an energy unit"),
    (_bad(orbital_pairs={"type": "nocv", "pairs": []}), "eigenvalue_convention"),
])
def test_writer_refuses(payload, match):
    w = QvfWriter("test", "1.0", auto_metadata=False)
    w.add_structure(WATER_DIMER)
    with pytest.raises(ValueError, match=match):
        w.add_energy_decomposition(payload)


def test_sum_rule_tolerance(tmp_path):
    loose = _bad(sum_rules=[{"total": 0, "parts": [1, 2, 3, 4],
                             "tolerance": {"value": 0.5, "units": "kcal_per_mol"}}])
    loose["terms"] = copy.deepcopy(ZR["terms"])
    loose["terms"][0]["energy"] = -0.0083  # 0.19 kcal/mol off
    del loose["orbital_pairs"]["pairs"][0]["density_ref"]
    path = _archive(tmp_path, loose, densities=False)
    assert validate_qvf(str(path))["ok"]


@pytest.mark.parametrize("edit, code", [
    (lambda p: p["sum_rules"][0]["parts"].pop(), "E_SUM_RULE"),
    (lambda p: p["fragments"][1]["atoms"].append(9), "E_MEMBER_MALFORMED"),
    (lambda p: p["fragments"][1]["atoms"].append(0), "E_MEMBER_MALFORMED"),
    (lambda p: p["terms"][1].__setitem__("parent", 2) or p["terms"][2].__setitem__(
        "parent", 1), "E_MEMBER_MALFORMED"),
    (lambda p: p["orbital_pairs"]["pairs"][1].__setitem__("index", 1), "E_MEMBER_MALFORMED"),
    (lambda p: p["orbital_pairs"]["pairs"][0].__setitem__("density_ref", "missing"),
     "E_REF_UNRESOLVED"),
    (lambda p: p["terms"][1]["energy"].__setitem__("units", "bohr"), "E_UNIT_CONTRADICTED"),
    (lambda p: p.__setitem__("scheme", "unknown-scheme"), "E_PAYLOAD_SCHEMA"),
])
def test_validator_findings(tmp_path, edit, code):
    pytest.importorskip("jsonschema") if code == "E_PAYLOAD_SCHEMA" else None
    path = _rewrite(_archive(tmp_path), tmp_path / "bad.qvf", payload_edit=edit)
    assert code in _codes(path)


@pytest.mark.parametrize("edit", [
    lambda s: s.__setitem__("pair_index", 7),
    lambda s: s.__setitem__("pair_spin", "alpha"),
    lambda s: s.__setitem__("decomposition_ref", "structure"),
    lambda s: s.pop("pair_index"),
])
def test_density_references(tmp_path, edit):
    path = _rewrite(_archive(tmp_path), tmp_path / "bad.qvf", section_edit=edit, sid="nocv_1")
    assert "E_REF_UNRESOLVED" in _codes(path) or "E_SCHEMA" in _codes(path)


def test_wavefunction_ref_resolved(tmp_path):
    path = _rewrite(_archive(tmp_path), tmp_path / "bad.qvf",
                    section_edit=lambda s: s.__setitem__("wavefunction_ref", "nocv_1"))
    assert "E_REF_UNRESOLVED" in _codes(path)


def test_volume_fields_checked_by_writer():
    w = QvfWriter("test", "1.0", auto_metadata=False)
    grid = (np.zeros((2, 2, 2)), [0.0, 0.0, 0.0], np.eye(3))
    with pytest.raises(ValueError, match="volume.difference only"):
        w.add_volume("volume.density", *grid, decomposition_ref="energy_decomposition")
    with pytest.raises(ValueError, match="needs pair_index"):
        w.add_volume("volume.difference", *grid, decomposition_ref="energy_decomposition", decomposition_role="orbital_pair")
    with pytest.raises(ValueError, match="need"):
        w.add_volume("volume.difference", *grid, pair_index=1)
