"""X-ray absorption/emission and Mössbauer spectra (spec § 5.4, format
review P13): the reference writer's builders, the shared one-dimensional
rules on ``spectra.xas`` / ``spectra.xes``, and the strict payloads."""
from __future__ import annotations

import hashlib
import json
import zipfile

import pytest

from qvf import validate_qvf
from qvf._kinds import PAYLOAD_DEFS
from qvf.writer import QvfWriter

FECL4 = [{"symbol": "Fe", "position": [0.0, 0.0, 0.0]}] + [
    {"symbol": "Cl", "position": p} for p in
    ([1.32, 1.32, 1.32], [-1.32, -1.32, 1.32], [-1.32, 1.32, -1.32],
     [1.32, -1.32, -1.32])]


def _writer() -> QvfWriter:
    w = QvfWriter("test", "1.0", auto_metadata=False)
    w.add_structure(FECL4)
    return w


def _codes(path) -> set[str]:
    return {f["code"] for f in validate_qvf(str(path))["findings"]}


def _rewrite_spectrum(src, dst, sid, edit):
    """Copy ``src`` to ``dst`` with section ``sid``'s spectrum payload edited
    and the member's sha256 and size updated to the new bytes."""
    with zipfile.ZipFile(src) as zin:
        manifest = json.loads(zin.read("manifest.json"))
        sec = next(s for s in manifest["sections"] if s["id"] == sid)
        member = sec["members"]["spectrum"]
        payload = json.loads(zin.read(member["path"]))
        edit(payload)
        data_new = json.dumps(payload).encode()
        member["sha256"] = hashlib.sha256(data_new).hexdigest()
        if "size" in member:
            member["size"] = len(data_new)
        with zipfile.ZipFile(dst, "w") as zout:
            for info in zin.infolist():
                if info.filename == "manifest.json":
                    continue
                data = (data_new
                        if info.filename == member["path"]
                        else zin.read(info.filename))
                zout.writestr(info.filename, data)
            zout.writestr("manifest.json", json.dumps(manifest))


def test_payload_definitions():
    assert PAYLOAD_DEFS[("spectra.xas", "spectrum")] == "SpectrumPayload"
    assert PAYLOAD_DEFS[("spectra.xes", "spectrum")] == "SpectrumPayload"
    assert (PAYLOAD_DEFS[("spectra.mossbauer", "spectrum")]
            == "MossbauerSpectrumPayload")


def test_xas_and_xes_round_trip(tmp_path):
    w = _writer()
    w.add_excited_states([{"index": 1, "excitation_energy": 7112.0}],
                         section_id="es")
    w.add_xray_spectrum("spectra.xas", [7112.0], [2e-5], edge="K",
                        element="Fe", atom_indices=[0], core_orbital="1s",
                        intensity="oscillator_strength", states_ref="es",
                        peaks=[{"x": 7112.0, "assignment": {"state_index": 1,
                                                            "atom_index": 0}}])
    w.add_xray_spectrum("spectra.xes", [7105.0, 7108.0], [0.4, 1.0],
                        edge="K", line="Kβ2,5", intensity="normalized")
    path = tmp_path / "x.qvf"
    w.write(str(path))
    report = validate_qvf(str(path))
    assert report["ok"], report["findings"]
    with zipfile.ZipFile(path) as z:
        xas = json.loads(z.read("spectra/spectra_xas.json"))
    assert xas["xray"] == {"edge": "K", "element": "Fe", "atom_indices": [0],
                           "core_orbital": "1s",
                           "intensity": "oscillator_strength"}


@pytest.mark.parametrize("kwargs, message", [
    ({"edge": "K-edge"}, "IUPAC"),
    ({"edge": "L4"}, "IUPAC"),
    ({"edge": "K", "core_orbital": "1S"}, "n l notation"),
    ({"edge": "K", "atom_indices": [0, 0]}, "distinct"),
    ({"edge": "K", "intensity": "counts"}, "not one of"),
    ({"edge": "K", "intensity": "cross_section_mb",
      "y_units": "dimensionless"}, "megabarn"),
    ({"edge": "K", "intensity": "oscillator_strength",
      "y_units": "km_per_mol"}, "dimensionless"),
    ({"edge": "K", "element": "fe"}, "chemical symbol"),
])
def test_xray_builder_refuses(kwargs, message):
    with pytest.raises(ValueError, match=message):
        _writer().add_xray_spectrum("spectra.xas", [1.0], [1.0], **kwargs)


def test_xray_builder_refuses_other_kinds():
    with pytest.raises(ValueError, match="spectra.xas or spectra.xes"):
        _writer().add_xray_spectrum("spectra.uvvis", [1.0], [1.0], edge="K")


def test_xes_length_rule_applies(tmp_path):
    w = _writer()
    w.add_xray_spectrum("spectra.xes", [7105.0, 7108.0], [0.4, 1.0], edge="K")
    good, bad = tmp_path / "good.qvf", tmp_path / "bad.qvf"
    w.write(str(good))
    _rewrite_spectrum(good, bad, "spectra_xes",
                      lambda p: p.__setitem__("intensities", [1.0]))
    assert _codes(bad) == {"E_PAYLOAD_SCHEMA"}


@pytest.mark.parametrize("edit", [
    lambda p: p["xray"].__setitem__("edge", "L4"),
    lambda p: p["xray"].__setitem__("absorber", 0),
    lambda p: p.update({"y_units": "dimensionless",
                        "xray": dict(p["xray"], intensity="cross_section_mb")}),
])
def test_xray_block_schema(tmp_path, edit):
    w = _writer()
    w.add_xray_spectrum("spectra.xas", [7112.0], [1.0], edge="K",
                        intensity="normalized")
    good, bad = tmp_path / "good.qvf", tmp_path / "bad.qvf"
    w.write(str(good))
    assert validate_qvf(str(good))["ok"]
    _rewrite_spectrum(good, bad, "spectra_xas", edit)
    assert _codes(bad) == {"E_PAYLOAD_SCHEMA"}


def test_xray_x_keys_allowed(tmp_path):
    w = _writer()
    w.add_spectrum("spectra.xas", payload={
        "frequencies": [7112.0], "intensities": [1.0],
        "xray": {"edge": "L2,3", "x_orca_orbwin": [0, 0]}})
    path = tmp_path / "x.qvf"
    w.write(str(path))
    assert validate_qvf(str(path))["ok"]


def test_mossbauer_round_trip(tmp_path):
    w = _writer()
    w.add_mossbauer([{"atom_index": 0, "isomer_shift_mm_s": 0.31,
                      "quadrupole_splitting_mm_s": 0.0, "asymmetry": 0.0,
                      "hyperfine_field_t": -33.0, "x_site_note": "demo"}])
    path = tmp_path / "m.qvf"
    w.write(str(path))
    assert validate_qvf(str(path))["ok"]
    with zipfile.ZipFile(path) as z:
        payload = json.loads(z.read("spectra/spectra_mossbauer.json"))
    assert "isotope" not in payload            # absent means 57Fe
    assert payload["sites"][0]["hyperfine_field_t"] == -33.0


@pytest.mark.parametrize("sites, kwargs, message", [
    ([], {}, "non-empty"),
    ([{"isomer_shift_mm_s": 0.3}], {}, "atom_index"),
    ([{"atom_index": 0, "isomer_shift": 0.3}], {}, "unknown keys"),
    ([{"atom_index": 0, "asymmetry": 1.5}], {}, "at most 1"),
    ([{"atom_index": 0, "linewidth_mm_s": 0.0}], {}, "positive"),
    ([{"atom_index": 0, "efg_au": [[1, 0], [0, 1]]}], {}, "efg_au"),
    ([{"atom_index": 0}], {"isotope": "Fe57"}, "mass number first"),
    ([{"atom_index": 0}], {"calibration": {"alpha": -0.3}}, "alpha and c"),
])
def test_mossbauer_builder_refuses(sites, kwargs, message):
    with pytest.raises(ValueError, match=message):
        _writer().add_mossbauer(sites, **kwargs)


@pytest.mark.parametrize("edit", [
    lambda p: p["sites"][0].__setitem__("asymmetry", -0.1),
    lambda p: p["sites"][0].__setitem__("delta", 0.3),
    lambda p: p.__setitem__("sites", []),
    lambda p: p.__setitem__("units", {"isomer_shift_mm_s": "dimensionless"}),
])
def test_mossbauer_payload_is_strict(tmp_path, edit):
    w = _writer()
    w.add_mossbauer([{"atom_index": 0, "isomer_shift_mm_s": 0.31}])
    good, bad = tmp_path / "good.qvf", tmp_path / "bad.qvf"
    w.write(str(good))
    _rewrite_spectrum(good, bad, "spectra_mossbauer", edit)
    assert _codes(bad) == {"E_PAYLOAD_SCHEMA"}
