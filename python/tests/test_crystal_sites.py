"""Crystallographic atom sites: displacement parameters and disorder groups
(spec § 5.1, *Crystallographic atom sites*).

An atom of a ``structure`` payload may carry ``adp`` (CIF core
``_atom_site.ADP_type``, ``_atom_site_aniso.U_ij`` / ``B_ij``,
``_atom_site.U_iso_or_equiv``) and ``disorder_assembly`` /
``disorder_group`` (``_atom_site.disorder_assembly``,
``_atom_site.disorder_group``). These tests write the consumer's fixtures,
check every rule the validator adds, the spec's frame formulas, and the CIF
import and export of both.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from qvf import Reader, validate_qvf
from qvf._cif import _reciprocal_lengths, _rotate_crystal_u
from qvf.convert import _lattice_inverse, export_cif, import_cif
from qvf.writer import QvfWriter

BETA = math.radians(100.0)
MONOCLINIC = [[7.12, 0.0, 0.0], [0.0, 9.48, 0.0],
              [10.31 * math.cos(BETA), 0.0, 10.31 * math.sin(BETA)]]
U_CRYSTAL = [0.0381, 0.0455, 0.0402, 0.0023, 0.0081, -0.0017]


def _cart(frac, lattice=MONOCLINIC):
    return [sum(frac[c] * lattice[c][r] for c in range(3)) for r in range(3)]


def _archive(tmp_path, atoms, *, pbc=(True, True, True), lattice=MONOCLINIC,
             name="crystal.qvf"):
    w = QvfWriter("test", "0", "crystal sites")
    w.add_structure(atoms, pbc=list(pbc), lattice_vectors=lattice if any(pbc) else None)
    path = tmp_path / name
    w.write(str(path))
    return path


def _codes(path):
    return [f["code"] for f in validate_qvf(str(path))["findings"] if f["severity"] == "error"]


def _matrix(u):
    return np.array([[u[0], u[3], u[4]], [u[3], u[1], u[5]], [u[4], u[5], u[2]]])


def _crystal_to_cartesian(u, lattice):
    """U_cart = A N U N Aᵀ (spec § 5.1), A the lattice vectors as columns."""
    a = np.array(lattice).T
    n = np.diag(np.linalg.norm(np.linalg.inv(np.array(lattice)), axis=0))
    return a @ n @ _matrix(u) @ n @ a.T


# --- writer and validator -----------------------------------------------------

def _disordered_atoms():
    """The request's fixtures: one anisotropic atom, and a site disordered
    over two groups with occupancy below 1."""
    return [
        {"symbol": "Fe", "position": _cart([0.25, 0.25, 0.25]), "atom_name": "Fe1",
         "adp": {"type": "Uani", "frame": "crystal", "iso": 0.0412, "aniso": U_CRYSTAL}},
        {"symbol": "O", "position": _cart([0.6, 0.4, 0.1]), "atom_name": "O1A",
         "occupancy": 0.55, "disorder_assembly": "A", "disorder_group": "1",
         "adp": {"type": "Uiso", "iso": 0.031}},
        {"symbol": "O", "position": _cart([0.62, 0.43, 0.12]), "atom_name": "O1B",
         "occupancy": 0.45, "disorder_assembly": "A", "disorder_group": "2",
         "adp": {"type": "Biso", "iso": 2.4}},
    ]


def test_writer_round_trips_the_fields(tmp_path):
    atoms = _disordered_atoms()
    path = _archive(tmp_path, atoms)
    report = validate_qvf(str(path))
    assert report["ok"], report["errors"]
    with Reader(str(path)) as reader:
        back = reader.json("structure", "structure")["atoms"]
    for want, got in zip(atoms, back):
        for key in ("adp", "occupancy", "disorder_assembly", "disorder_group"):
            assert got.get(key) == want.get(key)


def test_cartesian_frame_is_valid_on_a_molecule(tmp_path):
    atom = {"symbol": "O", "position": [0.0, 0.0, 0.0],
            "adp": {"type": "Bani", "frame": "cartesian",
                    "aniso": [1.2, 1.5, 1.1, 0.1, 0.0, -0.05]},
            "disorder_group": "1"}
    assert validate_qvf(str(_archive(tmp_path, [atom], pbc=(False,) * 3)))["ok"]


@pytest.mark.parametrize("pbc", [(False, False, False), (True, True, False)])
def test_crystal_frame_needs_three_periodic_axes(tmp_path, pbc):
    atom = {"symbol": "Fe", "position": [1.0, 1.0, 1.0],
            "adp": {"type": "Uani", "frame": "crystal", "aniso": U_CRYSTAL}}
    assert _codes(_archive(tmp_path, [atom], pbc=pbc)) == ["E_MEMBER_MALFORMED"]


@pytest.mark.parametrize("aniso", [
    [0.02, 0.02, 0.02, 0.03, 0.0, 0.0],      # second minor negative
    [0.02, 0.02, 0.02, 0.02, 0.0, 0.0],      # singular: positive semidefinite only
    [-0.01, 0.02, 0.02, 0.0, 0.0, 0.0],      # a negative axis
    [0.02, 0.02, 0.001, 0.0, 0.019, 0.019],  # leading minors 1, 2 positive, det < 0
])
def test_tensor_must_be_positive_definite(tmp_path, aniso):
    atom = {"symbol": "Fe", "position": [1.0, 1.0, 1.0],
            "adp": {"type": "Uani", "frame": "cartesian", "aniso": aniso}}
    assert _codes(_archive(tmp_path, [atom])) == ["E_MEMBER_MALFORMED"]


@pytest.mark.parametrize("extra", [
    {"adp": {"type": "Uani", "aniso": U_CRYSTAL}},                       # no frame
    {"adp": {"type": "Uani", "frame": "crystal"}},                       # no components
    {"adp": {"type": "Uiso"}},                                           # no value
    {"adp": {"type": "Uiso", "iso": 0.0}},                               # not positive
    {"adp": {"type": "Biso", "iso": 1.0, "frame": "crystal"}},           # frame on iso
    {"adp": {"type": "Uovl", "iso": 0.02}},                              # not carried
    {"adp": {"type": "Uani", "frame": "crystal", "aniso": U_CRYSTAL[:5]}},
    {"adp": {"type": "Uiso", "iso": 0.02, "units": "angstrom^2"}},       # closed object
    {"disorder_assembly": "A"},                                          # no group
    {"disorder_group": "?"},
    {"disorder_group": "1 2"},
    {"disorder_group": 1},                                               # a code, not a number
])
def test_schema_refuses_malformed_fields(tmp_path, extra):
    atom = {"symbol": "Fe", "position": [1.0, 1.0, 1.0]}
    atom.update(extra)
    assert _codes(_archive(tmp_path, [atom])) == ["E_PAYLOAD_SCHEMA"]


# --- the spec's frame formulas ----------------------------------------------

def test_u_equiv_of_the_dictionary_is_a_third_of_the_cartesian_trace():
    """CIF core `_atom_site.U_iso_or_equiv`: U(equiv) = 1/3 Σ_ij U^ij a*_i
    a*_j a_i·a_j, which the spec's U_cart = A N U N Aᵀ must reproduce as
    trace(U_cart) / 3."""
    a = np.array(MONOCLINIC)
    nstar = np.linalg.norm(np.linalg.inv(a), axis=0)
    u = _matrix(U_CRYSTAL)
    dictionary = sum(u[i, j] * nstar[i] * nstar[j] * a[i] @ a[j]
                     for i in range(3) for j in range(3)) / 3
    assert np.trace(_crystal_to_cartesian(U_CRYSTAL, MONOCLINIC)) / 3 == pytest.approx(dictionary)


def test_symmetry_copy_tensor_in_both_frames():
    """The spec's copy rule: N⁻¹ R N U N Rᵀ N⁻¹ in the crystal frame is
    (A R A⁻¹) U_cart (A R A⁻¹)ᵀ in the Cartesian frame (here the 2-fold axis
    along b and the c-glide of P 21/c, in a monoclinic cell)."""
    nstar = _reciprocal_lengths(_lattice_inverse(MONOCLINIC))
    a = np.array(MONOCLINIC).T
    for rot in ([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], [[1, 0, 0], [0, -1, 0], [0, 0, 1]]):
        copy = _rotate_crystal_u(U_CRYSTAL, rot, nstar)
        q = a @ np.array(rot, dtype=float) @ np.linalg.inv(a)
        want = q @ _crystal_to_cartesian(U_CRYSTAL, MONOCLINIC) @ q.T
        assert np.allclose(_crystal_to_cartesian(copy, MONOCLINIC), want, atol=1e-9)
        assert np.allclose(np.linalg.eigvalsh(_matrix(copy) * np.outer(nstar, nstar)),
                           np.linalg.eigvalsh(_matrix(U_CRYSTAL) * np.outer(nstar, nstar)))


# --- CIF import and export ----------------------------------------------------

P21C = """data_disordered
_cell_length_a 7.12
_cell_length_b 9.48
_cell_length_c 10.31
_cell_angle_alpha 90
_cell_angle_beta 100
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P 1 21/c 1'
loop_
_space_group_symop_operation_xyz
'x, y, z'
'-x, y+1/2, -z+1/2'
'-x, -y, -z'
'x, -y+1/2, z+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_disorder_assembly
_atom_site_disorder_group
Fe1 Fe 0.2 0.3 0.4 1 0.0412(3) Uani . .
C7A C 0.4 0.1 0.3 0.62 0.051 Uiso A 1
C7B C 0.42 0.12 0.27 0.38 0.055 Uiso A 2
Bad1 Cl 0.7 0.8 0.9 1 0.04 Uani . .
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
Fe1 0.0381(4) 0.0455(5) 0.0402(4) 0.0023(3) 0.0081(3) -0.0017(3)
Bad1 0.02 0.02 0.02 0.03 0 0
"""


def _import(tmp_path, text=P21C):
    src = tmp_path / "in.cif"
    src.write_text(text)
    out = tmp_path / "in.qvf"
    report = import_cif(src, out)
    assert validate_qvf(str(out))["ok"]
    with Reader(str(out)) as reader:
        return report, reader.json("structure", "structure"), out


def test_cif_import_carries_displacements_and_disorder(tmp_path):
    report, geometry, _ = _import(tmp_path)
    atoms = geometry["atoms"]
    assert len(atoms) == 16 and report["anisotropic_sites"] == 1
    assert report["disordered_sites"] == 2
    assert report["displacement_parameters_dropped"] == 1      # Bad1 is not positive
    fe = [a for a in atoms if a["atom_name"] == "Fe1"]
    assert fe[0]["adp"] == {"type": "Uani", "frame": "crystal", "iso": 0.0412,
                            "aniso": U_CRYSTAL}
    # Each copy carries the tensor of its own orientation: the 2-fold copy
    # changes the sign of U12 and U23, the inversion leaves U unchanged.
    tensors = sorted(tuple(a["adp"]["aniso"]) for a in fe)
    assert tensors == sorted([tuple(U_CRYSTAL)] * 2 + [
        (0.0381, 0.0455, 0.0402, -0.0023, 0.0081, 0.0017)] * 2)
    assert all(a["adp"] == {"type": "Uiso", "iso": 0.04}
               for a in atoms if a["atom_name"] == "Bad1")
    c7a = [a for a in atoms if a["atom_name"] == "C7A"]
    assert {(a["disorder_assembly"], a["disorder_group"], a["occupancy"]) for a in c7a} \
        == {("A", "1", 0.62)} and len(c7a) == 4
    assert c7a[0]["adp"] == {"type": "Uiso", "iso": 0.051}
    assert all("disorder_group" not in a for a in fe)


def test_cif_export_round_trips_and_converts_cartesian_tensors(tmp_path):
    _, geometry, qvf_in = _import(tmp_path)
    cif = tmp_path / "out.cif"
    export_cif(qvf_in, cif)
    again = tmp_path / "again.qvf"
    import_cif(cif, again)
    assert validate_qvf(str(again))["ok"]
    with Reader(str(again)) as reader:
        back = reader.json("structure", "structure")["atoms"]
    assert [a.get("adp") for a in back] == [a.get("adp") for a in geometry["atoms"]]
    assert [(a.get("disorder_assembly"), a.get("disorder_group")) for a in back] \
        == [(a.get("disorder_assembly"), a.get("disorder_group")) for a in geometry["atoms"]]

    # A Cartesian-frame B tensor leaves as the crystal-frame U it equals.
    u_cart = _crystal_to_cartesian(U_CRYSTAL, MONOCLINIC)
    b_cart = 8 * math.pi ** 2 * u_cart
    atom = {"symbol": "Fe", "position": _cart([0.1, 0.2, 0.3]), "atom_name": "Fe1",
            "adp": {"type": "Bani", "frame": "cartesian",
                    "aniso": [b_cart[0, 0], b_cart[1, 1], b_cart[2, 2],
                              b_cart[0, 1], b_cart[0, 2], b_cart[1, 2]]}}
    src = _archive(tmp_path, [atom], name="cartesian.qvf")
    export_cif(src, tmp_path / "cartesian.cif")
    import_cif(tmp_path / "cartesian.cif", tmp_path / "cartesian_back.qvf")
    with Reader(str(tmp_path / "cartesian_back.qvf")) as reader:
        adp = reader.json("structure", "structure")["atoms"][0]["adp"]
    assert adp["type"] == "Uani" and adp["frame"] == "crystal"
    # The exported cell is re-oriented (a along x), so compare frame-free
    # quantities: the crystal-frame components themselves.
    assert adp["aniso"] == pytest.approx(U_CRYSTAL, abs=2e-6)
