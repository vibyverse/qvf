"""Periodic bond images and the symmetry site map (spec § 5.1, § 5.7).

A bond of a periodic structure may join atom i to an image of atom j in a
neighbouring cell (``pairs[].image``), and a ``structure.symmetry`` payload
may say which site, operation and lattice translation generated each atom
(``sites``, ``atom_sites``). These tests write the fixtures the consumer
request named, check that they validate, and check every rule the validator
enforces with an archive that breaks it.
"""
from __future__ import annotations

import math

import pytest

from qvf import validate_qvf
from qvf.writer import QvfWriter


# --- fixtures -------------------------------------------------------------------

A_GRAPHENE = 2.46
GRAPHENE_LATTICE = [[A_GRAPHENE, 0.0, 0.0],
                    [-A_GRAPHENE / 2, A_GRAPHENE * math.sqrt(3) / 2, 0.0],
                    [0.0, 0.0, 15.0]]
#: Each carbon bonds to the other through three images (1.42 Å each).
GRAPHENE_BONDS = [{"i": 0, "j": 1, "image": [0, 0, 0]},
                  {"i": 0, "j": 1, "image": [-1, 0, 0]},
                  {"i": 0, "j": 1, "image": [0, 1, 0]}]

A_CUBIC = 3.35
#: One atom bonded to its own six neighbours: three bonds, each written once.
CUBIC_BONDS = [{"i": 0, "j": 0, "image": [1, 0, 0]},
               {"i": 0, "j": 0, "image": [0, 1, 0]},
               {"i": 0, "j": 0, "image": [0, 0, 1]}]


def _cart(frac, lattice):
    return [sum(frac[c] * lattice[c][r] for c in range(3)) for r in range(3)]


def _graphene_atoms():
    return [{"symbol": "C", "position": _cart(f, GRAPHENE_LATTICE)}
            for f in ([1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5])]


def _parse_op(text):
    """An operation 'x+1/2,-y,z' as (R, t) on fractional coordinates."""
    rot, shift = [], []
    for part in text.split(","):
        row, t, sign, i = [0, 0, 0], 0.0, 1, 0
        while i < len(part):
            ch = part[i]
            if ch in "+-":
                sign = 1 if ch == "+" else -1
            elif ch in "xyz":
                row["xyz".index(ch)] = sign
                sign = 1
            elif ch.isdigit():
                j = i
                while j < len(part) and (part[j].isdigit() or part[j] == "/"):
                    j += 1
                num, den = part[i:j].split("/")
                t += sign * int(num) / int(den)
                sign, i = 1, j - 1
            i += 1
        rot.append(row)
        shift.append(t)
    return rot, shift


def _expand(op_texts, sites, lattice):
    """Expand ``sites`` (label, symbol, fractional position) with the
    operations into the reference cell, removing copies that coincide, and
    record each atom's site, operation and translation."""
    ops = [_parse_op(t) for t in op_texts]
    atoms, atom_sites = [], []
    for s, (_, symbol, pos) in enumerate(sites):
        seen = []
        for o, (rot, shift) in enumerate(ops):
            raw = [sum(rot[r][c] * pos[c] for c in range(3)) + shift[r] for r in range(3)]
            n = [-math.floor(x + 1e-9) for x in raw]
            frac = [raw[r] + n[r] for r in range(3)]
            if any(max(abs(frac[r] - f[r]) for r in range(3)) < 1e-6 for f in seen):
                continue
            seen.append(frac)
            atoms.append({"symbol": symbol, "position": _cart(frac, lattice)})
            atom_sites.append({"site": s, "operation": o, "translation": n})
    operations = [{"rotation": r, "translation": t} for r, t in ops]
    return atoms, operations, atom_sites


#: P 42/m n m (No. 136), the sixteen operations of the standard setting.
P42MNM = ["x,y,z", "-x,-y,z", "-y+1/2,x+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2",
          "-x+1/2,y+1/2,-z+1/2", "x+1/2,-y+1/2,-z+1/2", "y,x,-z", "-y,-x,-z",
          "-x,-y,-z", "x,y,-z", "y+1/2,-x+1/2,-z+1/2", "-y+1/2,x+1/2,-z+1/2",
          "x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,z+1/2", "-y,-x,z", "y,x,z"]
RUTILE_LATTICE = [[4.594, 0.0, 0.0], [0.0, 4.594, 0.0], [0.0, 0.0, 2.959]]
#: Ti on 2a, a special position; O on 4f, deposited outside the reference
#: cell so that its images carry nonzero translations.
RUTILE_SITES = [("Ti1", "Ti", [0.0, 0.0, 0.0]), ("O1", "O", [1.3049, 0.3049, 0.0])]

#: P 1 21/n 1, a non-default setting of P 21/c (No. 14), in a skew cell.
P21N = ["x,y,z", "-x+1/2,y+1/2,-z+1/2", "-x,-y,-z", "x+1/2,-y+1/2,z+1/2"]
BETA = math.radians(97.0)
P21N_LATTICE = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0],
                [7.0 * math.cos(BETA), 0.0, 7.0 * math.sin(BETA)]]
P21N_SITES = [("C1", "C", [0.11, 0.22, 0.33])]


def _symmetry_archive(tmp_path, op_texts, sites, lattice, *, edit=None,
                      atoms_edit=None, pbc=(True, True, True), name="site_map.qvf"):
    atoms, operations, atom_sites = _expand(op_texts, sites, lattice)
    data = {"source": "deposited", "operations": operations,
            "sites": [{"label": label, "position": pos} for label, _, pos in sites],
            "atom_sites": atom_sites}
    if edit is not None:
        edit(data)
    if atoms_edit is not None:
        atoms_edit(atoms)
    w = QvfWriter("test", "1.0", auto_metadata=False)
    w.add_structure(atoms, pbc=list(pbc), lattice_vectors=lattice)
    w.add_structure_symmetry(data)
    path = tmp_path / name
    w.write(str(path))
    return path


def _bonds_archive(tmp_path, atoms, pairs, *, pbc=None, lattice=None,
                   standalone=False, bond_orders=False, name="bonds.qvf"):
    w = QvfWriter("test", "1.0", auto_metadata=False)
    if standalone or bond_orders:
        w.add_structure(atoms, pbc=pbc, lattice_vectors=lattice)
        if bond_orders:
            w.add_bond_orders("mayer", [dict(p, order=p.get("order", 1.0)) for p in pairs])
        else:
            w.add_bonds(pairs, structure_ref="structure")
    else:
        w.add_structure(atoms, pbc=pbc, lattice_vectors=lattice, bonds=pairs)
    path = tmp_path / name
    w.write(str(path))
    return path


def _findings(path):
    return [(f["code"], f["where"]) for f in validate_qvf(str(path))["findings"]
            if f["severity"] == "error"]


def _codes(path):
    return {code for code, _ in _findings(path)}


# --- periodic bond images: archives that validate ----------------------------------

@pytest.mark.parametrize("standalone", [False, True])
def test_graphene_cell_bonds_through_three_images(tmp_path, standalone):
    path = _bonds_archive(tmp_path, _graphene_atoms(), GRAPHENE_BONDS,
                          pbc=[True, True, False], lattice=GRAPHENE_LATTICE,
                          standalone=standalone)
    report = validate_qvf(str(path))
    assert report["ok"], report["findings"]
    assert report["findings"] == []


def test_cubic_cell_bonds_to_its_own_images(tmp_path):
    path = _bonds_archive(tmp_path, [{"symbol": "Po", "position": [0.0, 0.0, 0.0]}],
                          CUBIC_BONDS, pbc=[True, True, True],
                          lattice=[[A_CUBIC, 0, 0], [0, A_CUBIC, 0], [0, 0, A_CUBIC]])
    assert validate_qvf(str(path))["ok"]


def test_molecule_bonds_without_images_keep_their_meaning(tmp_path):
    water = [{"symbol": "O", "position": [0.0, 0.0, 0.117]},
             {"symbol": "H", "position": [0.757, 0.0, -0.470]},
             {"symbol": "H", "position": [-0.757, 0.0, -0.470]}]
    # A repeated pair without images stays valid, as it was before images.
    pairs = [{"i": 0, "j": 1}, {"i": 0, "j": 2}, {"i": 1, "j": 0}]
    assert validate_qvf(str(_bonds_archive(tmp_path, water, pairs)))["ok"]


def test_bond_orders_carry_images(tmp_path):
    path = _bonds_archive(tmp_path, _graphene_atoms(), GRAPHENE_BONDS,
                          pbc=[True, True, False], lattice=GRAPHENE_LATTICE,
                          bond_orders=True)
    assert validate_qvf(str(path))["ok"]


def test_image_along_a_leading_aperiodic_axis_is_zero(tmp_path):
    """The periodic axes need not lead (spec § 5.1): a slab periodic along
    y and z bonds through a [0, 0, 1] image."""
    lattice = [[20.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]]
    atoms = [{"symbol": "H", "position": [0.0, 0.0, 0.0]}]
    path = _bonds_archive(tmp_path, atoms, [{"i": 0, "j": 0, "image": [0, 0, 1]}],
                          pbc=[False, True, True], lattice=lattice)
    assert validate_qvf(str(path))["ok"]


# --- periodic bond images: archives that must fail -------------------------------

def _graphene_fail(tmp_path, pairs, **kw):
    kw.setdefault("pbc", [True, True, False])
    kw.setdefault("lattice", GRAPHENE_LATTICE)
    return _bonds_archive(tmp_path, _graphene_atoms(), pairs, **kw)


def test_self_bond_with_a_zero_image_is_refused(tmp_path):
    path = _graphene_fail(tmp_path, [{"i": 1, "j": 1, "image": [0, 0, 0]}])
    assert _findings(path) == [("E_MEMBER_MALFORMED", "/sections/0/members/bonds")]


def test_image_that_is_not_an_integer_is_refused(tmp_path):
    path = _graphene_fail(tmp_path, [{"i": 0, "j": 1, "image": [0.5, 0, 0]}])
    assert _codes(path) == {"E_PAYLOAD_SCHEMA"}


def test_image_with_two_components_is_refused(tmp_path):
    path = _graphene_fail(tmp_path, [{"i": 0, "j": 1, "image": [1, 0]}])
    assert _codes(path) == {"E_PAYLOAD_SCHEMA"}


def test_image_on_a_molecule_is_refused(tmp_path):
    atoms = [{"symbol": "H", "position": [0.0, 0.0, 0.0]},
             {"symbol": "H", "position": [0.0, 0.0, 0.74]}]
    path = _bonds_archive(tmp_path, atoms, [{"i": 0, "j": 1, "image": [0, 0, 0]}],
                          standalone=True)
    assert _findings(path) == [("E_MEMBER_MALFORMED", "/sections/1/members/bonds")]


def test_image_along_an_aperiodic_axis_is_refused(tmp_path):
    path = _graphene_fail(tmp_path, [{"i": 0, "j": 1, "image": [0, 0, 1]}])
    assert _codes(path) == {"E_MEMBER_MALFORMED"}


@pytest.mark.parametrize("second", [
    {"i": 0, "j": 1, "image": [-1, 0, 0]},     # the same record twice
    {"i": 1, "j": 0, "image": [1, 0, 0]},      # its reversal
])
def test_one_bond_written_twice_is_refused(tmp_path, second):
    path = _graphene_fail(tmp_path, [GRAPHENE_BONDS[1], second])
    assert _codes(path) == {"E_MEMBER_MALFORMED"}


def test_self_bond_and_its_reversal_are_one_bond(tmp_path):
    path = _bonds_archive(tmp_path, [{"symbol": "Po", "position": [0.0, 0.0, 0.0]}],
                          [{"i": 0, "j": 0, "image": [1, 0, 0]},
                           {"i": 0, "j": 0, "image": [-1, 0, 0]}],
                          pbc=[True, True, True],
                          lattice=[[A_CUBIC, 0, 0], [0, A_CUBIC, 0], [0, 0, A_CUBIC]])
    assert _codes(path) == {"E_MEMBER_MALFORMED"}


def test_bond_orders_follow_the_same_rules(tmp_path):
    path = _graphene_fail(tmp_path, [GRAPHENE_BONDS[0], dict(GRAPHENE_BONDS[0], i=1, j=0)],
                          bond_orders=True)
    assert _findings(path) == [("E_MEMBER_MALFORMED", "/sections/1/members/bond_orders")]


# --- the site map: archives that validate ------------------------------------------

def test_rutile_site_map_on_a_special_position(tmp_path):
    path = _symmetry_archive(tmp_path, P42MNM, RUTILE_SITES, RUTILE_LATTICE)
    report = validate_qvf(str(path))
    assert report["ok"], report["findings"]
    atoms, _, atom_sites = _expand(P42MNM, RUTILE_SITES, RUTILE_LATTICE)
    assert len(atoms) == 6                       # Ti on 2a, O on 4f
    assert [e["site"] for e in atom_sites] == [0, 0, 1, 1, 1, 1]
    assert any(e["translation"] != [0, 0, 0] for e in atom_sites)


def test_non_default_setting_in_a_skew_cell(tmp_path):
    path = _symmetry_archive(tmp_path, P21N, P21N_SITES, P21N_LATTICE)
    assert validate_qvf(str(path))["ok"]


def test_molecule_site_map_in_cartesian_coordinates(tmp_path):
    data = {"point_group": "C2v", "source": "detected", "symprec": 1e-3,
            "operations": [{"rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]},
                           {"rotation": [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]}],
            "sites": [{"label": "O1", "position": [0.0, 0.0, 0.117]},
                      {"label": "H1", "position": [0.757, 0.0, -0.470]}],
            "atom_sites": [{"site": 0, "operation": 0}, {"site": 1, "operation": 0},
                           {"site": 1, "operation": 1, "translation": [0, 0, 0]}]}
    water = [{"symbol": "O", "position": [0.0, 0.0, 0.117]},
             {"symbol": "H", "position": [0.757, 0.0, -0.470]},
             {"symbol": "H", "position": [-0.757, 0.0, -0.470]}]
    w = QvfWriter("test", "1.0", auto_metadata=False)
    w.add_structure(water)
    w.add_structure_symmetry(data)
    path = tmp_path / "water.qvf"
    w.write(str(path))
    assert validate_qvf(str(path))["ok"]
    data["atom_sites"][2]["translation"] = [1, 0, 0]
    w = QvfWriter("test", "1.0", auto_metadata=False)
    w.add_structure(water)
    w.add_structure_symmetry(data)
    w.write(str(path))
    assert _codes(path) == {"E_MEMBER_MALFORMED"}


def test_symmetry_without_a_site_map_is_unchanged(tmp_path):
    path = _symmetry_archive(tmp_path, P42MNM, RUTILE_SITES, RUTILE_LATTICE,
                             edit=lambda d: [d.pop("atom_sites"), d.pop("sites")])
    assert validate_qvf(str(path))["ok"]


def test_site_tolerance_widens_the_position_check(tmp_path):
    def shift(atoms):
        atoms[3]["position"][0] += 0.05
    assert _codes(_symmetry_archive(tmp_path, P42MNM, RUTILE_SITES, RUTILE_LATTICE,
                                    atoms_edit=shift)) == {"E_MEMBER_MALFORMED"}
    path = _symmetry_archive(tmp_path, P42MNM, RUTILE_SITES, RUTILE_LATTICE,
                             atoms_edit=shift,
                             edit=lambda d: d.update(site_tolerance=0.1))
    assert validate_qvf(str(path))["ok"]


# --- the site map: archives that must fail --------------------------------------

def _bad_entry(**change):
    def edit(data):
        data["atom_sites"][2].update(change)
    return edit


@pytest.mark.parametrize("edit, code", [
    (_bad_entry(operation=16), "E_MEMBER_MALFORMED"),
    (_bad_entry(site=2), "E_MEMBER_MALFORMED"),
    (_bad_entry(operation=3), "E_MEMBER_MALFORMED"),          # a wrong position
    (_bad_entry(translation=[0, 0, 0.5]), "E_PAYLOAD_SCHEMA"),
    (_bad_entry(translation=[1, 0, 0]), "E_MEMBER_MALFORMED"),
    (lambda d: d["atom_sites"].pop(), "E_MEMBER_MALFORMED"),  # one entry short
    (lambda d: d.pop("sites"), "E_MEMBER_MALFORMED"),
    (lambda d: d.pop("operations"), "E_MEMBER_MALFORMED"),
    (lambda d: d.update(source="guessed"), "E_PAYLOAD_SCHEMA"),
    (lambda d: d["sites"][0].pop("position"), "E_PAYLOAD_SCHEMA"),
])
def test_site_map_rules(tmp_path, edit, code):
    path = _symmetry_archive(tmp_path, P42MNM, RUTILE_SITES, RUTILE_LATTICE, edit=edit)
    assert _codes(path) == {code}


def test_translation_along_an_aperiodic_axis_is_refused(tmp_path):
    def edit(data):
        data["atom_sites"][0]["translation"] = [0, 0, 1]
    def shift(atoms):
        atoms[0]["position"][2] += RUTILE_LATTICE[2][2]
    path = _symmetry_archive(tmp_path, P42MNM, RUTILE_SITES, RUTILE_LATTICE,
                             edit=edit, atoms_edit=shift, pbc=(True, True, False))
    assert _findings(path) == [("E_MEMBER_MALFORMED", "/sections/1/members/data")]


def test_an_edited_structure_disagrees_with_its_site_map(tmp_path):
    """A writer that moves an atom must rewrite or drop the map (spec § 5.7)."""
    def move(atoms):
        atoms[0]["position"] = [1.0, 1.0, 1.0]
    path = _symmetry_archive(tmp_path, P42MNM, RUTILE_SITES, RUTILE_LATTICE,
                             atoms_edit=move)
    assert _findings(path) == [("E_MEMBER_MALFORMED", "/sections/1/members/data")]
