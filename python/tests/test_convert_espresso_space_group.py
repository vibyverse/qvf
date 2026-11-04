"""pw.x space_group input with crystal_sg sites (qvf.convert.import_espresso_in, WP-48).

The expansion is checked against the CIF importer on the same structures,
against ase.spacegroup.crystal and against spglib's symmetry database and
space-group finder (the last two skip when the package is absent).
"""
import math

import pytest

from qvf import Reader
from qvf._spacegroups import (ITA_HALL, ITA_HALL_ORIGIN_2, ITA_HALL_RHOMBOHEDRAL,
                              ITA_HALL_UNIQUE_C, ita_hall, operators)
from qvf.convert import _lattice_inverse, import_cif, import_espresso_in
from qvf.validate import run

R_GROUPS = (146, 148, 155, 160, 161, 166, 167)


def pw_input(number, sites, *, celldm=None, extra="", flags=None):
    """A pw.x input of space group ``number`` with (label, x, y, z) sites."""
    celldm = celldm or {1: 9.0}
    labels = list(dict.fromkeys(s[0] for s in sites))
    cell = ", ".join(f"celldm({k}) = {v!r}" for k, v in celldm.items())
    lines = [f"&SYSTEM\n  space_group = {number}, {cell}, nat = {len(sites)}, "
             f"ntyp = {len(labels)}{extra}\n/", "ATOMIC_SPECIES"]
    lines += [f"{label} 1.0 {label}.UPF" for label in labels]
    lines.append("ATOMIC_POSITIONS crystal_sg")
    for k, (label, *x) in enumerate(sites):
        lines.append(f"{label} " + " ".join(str(v) for v in x)
                     + (f" {flags[k]}" if flags and flags[k] else ""))
    return "\n".join(lines) + "\n"


def imported(tmp_path, text, name="sg"):
    source, out = tmp_path / f"{name}.pwi", tmp_path / f"{name}.qvf"
    source.write_text(text)
    report = import_espresso_in(source, out)
    check = run(out)
    assert check["ok"] and check["findings"] == [] and report["findings"] == []
    with Reader(out) as reader:
        geometry = reader.json("structure", "structure")
    return report, geometry


def fractional(geometry):
    inverse = _lattice_inverse(geometry["lattice_vectors"])
    return [(a["symbol"], [sum(a["position"][k] * inverse[k][i] for k in range(3)) % 1.0
                           for i in range(3)]) for a in geometry["atoms"]]


def cellpar(lattice):
    lengths = [math.sqrt(sum(v * v for v in row)) for row in lattice]

    def angle(u, v):
        return math.degrees(math.acos(sum(x * y for x, y in zip(lattice[u], lattice[v]))
                                      / (lengths[u] * lengths[v])))
    return lengths + [angle(1, 2), angle(0, 2), angle(0, 1)]


def same_sites(ours, reference, tol=1e-6):
    """Both lists of (symbol, fractional) hold the same atoms modulo the lattice."""
    assert len(ours) == len(reference)
    left = list(reference)
    for symbol, f in ours:
        for k, (other, g) in enumerate(left):
            if other == symbol and all(abs((a - b + 0.5) % 1.0 - 0.5) < tol
                                       for a, b in zip(f, g)):
                del left[k]
                break
        else:
            raise AssertionError(f"{symbol} {f} has no counterpart")


# (number, settings, celldm, sites, atoms in the conventional cell)
STRUCTURES = [
    (225, "", {1: 10.66}, [("Na", 0, 0, 0), ("Cl", 0.5, 0.5, 0.5)], 8),
    (227, ", origin_choice = 2", {1: 10.26}, [("Si", 0.125, 0.125, 0.125)], 8),
    (227, ", origin_choice = 1", {1: 10.26}, [("Si", 0, 0, 0)], 8),
    (136, "", {1: 8.68, 3: 0.644}, [("Ti", 0, 0, 0), ("O", 0.3049, 0.3049, 0)], 6),
    (186, "", {1: 6.026, 3: 1.626}, [("Ga", 1 / 3, 2 / 3, 0), ("N", 1 / 3, 2 / 3, 0.377)], 4),
    (167, ", rhombohedral = .false.", {1: 8.99, 3: 2.73},
     [("Al", 0, 0, 0.3523), ("O", 0.3064, 0, 0.25)], 30),
    (167, ", rhombohedral = .true.", {1: 9.71, 4: 0.0857},
     [("Al", 0.3523, 0.3523, 0.3523), ("O", 0.5560, -0.0560, 0.25)], 10),
    (62, "", {1: 10.0, 2: 0.7, 3: 1.3}, [("Ca", 0.01, 0.25, 0.99), ("O", 0.2, 0.05, 0.3)], 12),
    (14, ", uniqueb = .true.", {1: 9.0, 2: 1.2, 3: 1.4, 5: -0.25},
     [("C", 0.1, 0.2, 0.3), ("O", 0, 0, 0.5)], 6),
    (14, "", {1: 9.0, 2: 1.2, 3: 1.4, 4: -0.25}, [("C", 0.1, 0.2, 0.3), ("O", 0, 0, 0.5)], 6),
    (15, "", {1: 9.0, 2: 1.2, 3: 1.4, 4: -0.25}, [("C", 0.1, 0.2, 0.3), ("O", 0, 0, 0)], 12),
    (2, "", {1: 9.0, 2: 1.1, 3: 1.2, 4: 0.1, 5: 0.2, 6: 0.3}, [("H", 0.1, 0.2, 0.3)], 2),
]


@pytest.mark.parametrize("number, settings, celldm, sites, count", STRUCTURES)
def test_expansion_matches_cif_import(tmp_path, number, settings, celldm, sites, count):
    report, geometry = imported(tmp_path, pw_input(number, sites, celldm=celldm, extra=settings))
    assert report["space_group"] == number and report["atoms"] == count
    assert report["positions_units"] == "crystal_sg" and report["sites"] == len(sites)
    assert len(geometry["atoms"]) == count
    a, b, c, alpha, beta, gamma = cellpar(geometry["lattice_vectors"])
    rows = [f"_cell_length_a {a:.10f}", f"_cell_length_b {b:.10f}", f"_cell_length_c {c:.10f}",
            f"_cell_angle_alpha {alpha:.10f}", f"_cell_angle_beta {beta:.10f}",
            f"_cell_angle_gamma {gamma:.10f}"]
    cif = ["data_sg", f"_space_group_name_Hall '{report['hall']}'", *rows, "loop_",
           "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
           "_atom_site_fract_y", "_atom_site_fract_z"]
    cif += [f"{s}{k} {s} {x!r} {y!r} {z!r}" for k, (s, x, y, z) in enumerate(sites)]
    (tmp_path / "sg.cif").write_text("\n".join(cif) + "\n")
    cif_report = import_cif(tmp_path / "sg.cif", tmp_path / "cif.qvf")
    assert cif_report["findings"] == [] and cif_report["atoms"] == count
    with Reader(tmp_path / "cif.qvf") as reader:
        reference = reader.json("structure", "structure")
    assert cellpar(reference["lattice_vectors"]) == pytest.approx([a, b, c, alpha, beta, gamma],
                                                                   rel=1e-6)
    same_sites(fractional(geometry), fractional(reference), tol=1e-5)


def test_diamond_origins_agree_after_shift(tmp_path):
    _, one = imported(tmp_path, pw_input(227, [("Si", 0, 0, 0)], celldm={1: 10.26},
                                         extra=", origin_choice = 1"), "o1")
    _, two = imported(tmp_path, pw_input(227, [("Si", 0.125, 0.125, 0.125)],
                                         celldm={1: 10.26}, extra=", origin_choice = 2"), "o2")
    shifted = [(s, [(v - 0.125) % 1.0 for v in f]) for s, f in fractional(one)]
    same_sites(shifted, fractional(two))


def test_cell_from_a_b_c_and_frozen_sites(tmp_path):
    text = pw_input(221, [("Cs", 0, 0, 0), ("Cl", 0.5, 0.5, 0.5)], celldm={1: 1.0},
                    flags=["0 0 0", None]).replace("celldm(1) = 1.0", "A = 4.12")
    report, geometry = imported(tmp_path, text)
    assert report["cell_source"] == "A" and report["atoms"] == 2
    assert geometry["lattice_vectors"][0] == pytest.approx([4.12, 0, 0])
    with Reader(tmp_path / "sg.qvf") as reader:
        assert reader.raw_manifest["constraints"]["frozen_atoms"] == [0]


def test_frozen_site_freezes_every_image(tmp_path):
    text = pw_input(225, [("Na", 0, 0, 0), ("Cl", 0.5, 0.5, 0.5)], celldm={1: 10.66},
                    flags=["0 0 0", "1 1 1"])
    report, _ = imported(tmp_path, text)
    assert report["frozen_atoms"] == 4
    with Reader(tmp_path / "sg.qvf") as reader:
        assert reader.raw_manifest["constraints"]["frozen_atoms"] == [0, 1, 2, 3]


def test_settings_that_do_not_apply_are_reported(tmp_path):
    text = pw_input(225, [("Na", 0, 0, 0)], celldm={1: 10.66},
                    extra=", ibrav = 2, uniqueb = .true., origin_choice = 2, "
                          "rhombohedral = .false.")
    report, geometry = imported(tmp_path, text)
    assert report["settings_ignored"] == ["ibrav", "origin_choice", "rhombohedral", "uniqueb"]
    assert len(geometry["atoms"]) == 4


@pytest.mark.parametrize("text, message", [
    (pw_input(225, [("Na", "4a")]), "Wyckoff label '4a'"),
    (pw_input(225, [("Na", "8c", 0.25)]), "Wyckoff label '8c'"),
    (pw_input(231, [("Na", 0, 0, 0)]), "1 to 230"),
    (pw_input(225, [("Na", 0, 0, 0)], extra=", origin_choice = 3"), "origin_choice"),
    (pw_input(225, [("Na", 0, 0, 0)], extra=", uniqueb = 1"), "logicals"),
    (pw_input(225, [("Na", 0, 0, 0)]).replace("crystal_sg", "crystal"), "require each other"),
    (pw_input(225, [("Na", 0, 0, 0)]).replace("space_group = 225,", "ibrav = 2,"),
     "require each other"),
    (pw_input(225, [("Na", 0, 0, 0)]) + "CELL_PARAMETERS bohr\n1 0 0\n0 1 0\n0 0 1\n",
     "not from CELL_PARAMETERS"),
    (pw_input(225, [("Na", 0, 0, 0)], celldm={2: 1.0}), "celldm(1) or A"),
    (pw_input(225, [("Na", 0, 0, 0), ("Cl", 0.5, 0.5, 0)]), "overlapping atoms"),
])
def test_refusals(tmp_path, text, message):
    (tmp_path / "bad.pwi").write_text(text)
    with pytest.raises(ValueError, match=message.replace("(", r"\(").replace(")", r"\)")):
        import_espresso_in(tmp_path / "bad.pwi", tmp_path / "bad.qvf")
    assert not (tmp_path / "bad.qvf").exists()


# -- independent checks ------------------------------------------------------

def generic_celldm(number, uniqueb=False, rhombohedral=True):
    if number <= 2:
        return {1: 9.0, 2: 1.1, 3: 1.2, 4: 0.1, 5: 0.2, 6: 0.3}
    if number <= 15:
        return {1: 9.0, 2: 1.2, 3: 1.4, (5 if uniqueb else 4): -0.25}
    if number <= 74:
        return {1: 9.0, 2: 1.2, 3: 1.4}
    if number <= 194:
        if number in R_GROUPS and rhombohedral:
            return {1: 9.0, 4: 0.3}
        return {1: 9.0, 3: 1.4}
    return {1: 9.0}


def settings():
    """(number, pw.x settings, keyword arguments of ita_hall) of every setting."""
    for n in range(1, 231):
        if 3 <= n <= 15:
            yield n, ", uniqueb = .true.", dict(unique_axis="b")
            yield n, "", dict(unique_axis="c")
        elif n in ITA_HALL_ORIGIN_2:
            yield n, ", origin_choice = 1", {}
            yield n, ", origin_choice = 2", dict(origin_choice=2)
        elif n in R_GROUPS:
            yield n, ", rhombohedral = .false.", {}
            yield n, "", dict(rhombohedral=True)
        else:
            yield n, "", {}


SETTINGS = list(settings())
# Three species on general positions: one orbit alone can have more
# symmetry than its group (one atom in P 1 is P -1).
GENERAL = [("C", 0.1234, 0.3417, 0.2791), ("N", 0.2917, 0.0731, 0.4489),
           ("O", 0.6123, 0.8341, 0.1577)]


def test_hall_table_shape():
    assert len(ITA_HALL) == 230 and len(SETTINGS) == 230 + 13 + 24 + 7
    assert set(ITA_HALL_UNIQUE_C) == set(range(3, 16)) and set(ITA_HALL_RHOMBOHEDRAL) == set(
        R_GROUPS)
    # The PDB table and the ITA table agree on the chiral groups both hold.
    assert ita_hall(19) == "P 2ac 2ab" and ita_hall(146, rhombohedral=True) == "P 3*"
    assert ita_hall(5, unique_axis="c") == "B 2"


def test_hall_symbols_against_spglib_database():
    spglib = pytest.importorskip("spglib")
    from fractions import Fraction
    by_symbol = {spglib.get_spacegroup_type(h).hall_symbol: h for h in range(1, 531)}
    for number, _, kwargs in SETTINGS:
        hall = ita_hall(number, **kwargs)
        h = by_symbol[hall]
        assert spglib.get_spacegroup_type(h).number == number
        data = spglib.get_symmetry_from_database(h)
        reference = {(tuple(map(tuple, r.astype(int).tolist())),
                      tuple(Fraction(round(v * 24), 24) % 1 for v in t))
                     for r, t in zip(data["rotations"], data["translations"])}
        assert set(operators(hall)) == reference, (number, hall)


def test_every_setting_expands_to_its_space_group(tmp_path):
    """General sites expanded in every pw.x setting have full orbits and
    spglib finds the space group back from the atoms alone."""
    spglib = pytest.importorskip("spglib")
    for k, (number, extra, kwargs) in enumerate(SETTINGS):
        celldm = generic_celldm(number, "uniqueb = .true." in extra, "false" not in extra)
        report, geometry = imported(tmp_path, pw_input(number, GENERAL, celldm=celldm,
                                                       extra=extra), f"g{k}")
        assert report["atoms"] == 3 * len(operators(ita_hall(number, **kwargs)))
        sites = fractional(geometry)
        dataset = spglib.get_symmetry_dataset(
            (geometry["lattice_vectors"], [f for _, f in sites],
             ["CNO".index(s) for s, _ in sites]), symprec=1e-4)
        assert dataset.number == number, (number, extra)


@pytest.mark.parametrize("number, settings, celldm, sites, count", STRUCTURES)
def test_against_ase_spacegroup(tmp_path, number, settings, celldm, sites, count):
    spacegroup = pytest.importorskip("ase.spacegroup")
    if 3 <= number <= 15 and "uniqueb" not in settings:
        pytest.skip("ASE tabulates monoclinic groups with unique axis b only")
    setting = 2 if ("origin_choice = 2" in settings or "rhombohedral = .true." in settings) \
        else 1
    report, geometry = imported(tmp_path, pw_input(number, sites, celldm=celldm, extra=settings))
    par = cellpar(geometry["lattice_vectors"])
    atoms = spacegroup.crystal([s for s, *_ in sites], [x for _, *x in sites],
                               spacegroup=number, setting=setting, cellpar=par)
    assert len(atoms) == report["atoms"] == count
    assert atoms.cell.cellpar().tolist() == pytest.approx(par, rel=1e-9)
    reference = [(s, [v % 1.0 for v in f]) for s, f in
                 zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions().tolist())]
    same_sites(fractional(geometry), reference)
