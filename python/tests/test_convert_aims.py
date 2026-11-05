"""FHI-aims geometry.in import/export (qvf.convert, WP-48)."""
import warnings

import pytest

from qvf import Reader, Writer
from qvf.convert import export_aims, import_aims
from qvf.validate import run

REL = 1e-10

# Synthetic periodic geometry with fractional and Cartesian atoms.
GEOMETRY = """# synthetic GaN cell
lattice_vector   3.1890000000    0.0000000000    0.0000000000
lattice_vector  -1.5945000000    2.7617646000    0.0000000000
lattice_vector   0.0000000000    0.0000000000    5.1850000000
atom_frac  0.3333333333  0.6666666667  0.0000000000  Ga
    constrain_relaxation .true.
    initial_moment 0.5
atom_frac  0.6666666667  0.3333333333  0.5000000000  Ga
    initial_charge 0.1
atom       1.5945000000  0.9205882000  1.9547450000  N_tight   # a species label
    constrain_relaxation x
    constrain_relaxation y
atom_frac  0.6666666667  0.3333333333  0.8770000000  N
    constrain_relaxation x
    constrain_relaxation y
    constrain_relaxation z
    velocity 0.1 0.0 0.0
"""

MOLECULE = """atom 0.0 0.0 0.1173 O
atom 0.0 0.7572 -0.4692 H
atom 0.0 -0.7572 -0.4692 H
empty 0.0 0.0 2.0 H
"""


# ASE reads only all-fractional or all-Cartesian files.
ASE_GEOMETRY = GEOMETRY.replace(
    "atom       1.5945000000  0.9205882000  1.9547450000  N_tight   # a species label",
    "atom_frac  0.3333333333  0.6666666667  0.3770000000  N").replace(
    "    velocity 0.1 0.0 0.0\n", "")


def write(tmp_path, text, name="geometry.in"):
    path = tmp_path / name
    path.write_text(text)
    return path


def imported(tmp_path, text):
    out = tmp_path / f'out{len(list(tmp_path.glob("out*.qvf")))}.qvf'
    report = import_aims(write(tmp_path, text), out)
    check = run(out)
    assert check["ok"] and check["findings"] == [] and report["findings"] == []
    with Reader(out) as reader:
        geometry = reader.json("structure", "structure")
        raw = reader.raw_manifest
    return report, geometry, raw, out


def flat(rows):
    return [v for row in rows for v in row]


LATTICE = [[3.189, 0, 0], [-1.5945, 2.7617646, 0], [0, 0, 5.185]]


def frac(f):
    return [sum(f[i] * LATTICE[i][j] for i in range(3)) for j in range(3)]


def test_periodic_fractional_and_constraints(tmp_path):
    report, geometry, raw, _ = imported(tmp_path, GEOMETRY)
    assert report["periodic"] is True and report["fractional_atoms"] == 3
    assert report["species_labels"] == {"Ga": "Ga", "N_tight": "N", "N": "N"}
    assert report["frozen_atoms"] == 2 and report["partially_fixed_atoms_omitted"] == 1
    assert report["keywords_omitted"] == {"initial_charge": 1, "initial_moment": 1,
                                          "velocity": 1}
    assert raw["constraints"] == {"frozen_atoms": [0, 3]}
    assert geometry["pbc"] == [True] * 3
    assert flat(geometry["lattice_vectors"]) == pytest.approx(flat(LATTICE), rel=REL)
    assert report["cell_volume_angstrom3"] == pytest.approx(3.189 * 2.7617646 * 5.185, rel=REL)
    assert [a["symbol"] for a in geometry["atoms"]] == ["Ga", "Ga", "N", "N"]
    assert geometry["atoms"][1]["position"] == pytest.approx(
        frac([0.6666666667, 0.3333333333, 0.5]), rel=REL)
    assert geometry["atoms"][2]["position"] == pytest.approx([1.5945, 0.9205882, 1.954745],
                                                             rel=REL)


def test_molecule_and_empty_sites(tmp_path):
    report, geometry, raw, _ = imported(tmp_path, MOLECULE)
    assert report["periodic"] is False and report["cell_volume_angstrom3"] is None
    assert report["keywords_omitted"] == {"empty": 1}
    assert geometry["pbc"] == [False] * 3 and "lattice_vectors" not in geometry
    assert len(geometry["atoms"]) == 3 and "constraints" not in raw


@pytest.mark.parametrize("text, message", [
    (GEOMETRY.replace("lattice_vector  -1.5945000000    2.7617646000    0.0000000000\n", ""),
     "no lattice_vector"),
    (GEOMETRY.replace("2.7617646000", "0.0"), "linearly dependent"),
    (MOLECULE + "atom_frac 0.5 0.5 0.5 H\n", "atom_frac requires three lattice_vector"),
    (MOLECULE.replace("O\n", "Qq\n"), "does not start with an element"),
    (MOLECULE + "pseudocore 0 0 0 Na\n", "keyword 'pseudocore' is not supported"),
    ("constrain_relaxation .true.\n" + MOLECULE, "must follow an atom"),
    (MOLECULE.replace("O\n", "O\n  constrain_relaxation maybe\n"), ".true., .false."),
    (MOLECULE.replace("0.1173 O", "0.1173"), "three coordinates and a species"),
    (MOLECULE.replace("0.1173", "x"), "three numbers"),
    ("# only a comment\n", "at least one atom"),
])
def test_refusals(tmp_path, text, message):
    with pytest.raises(ValueError, match=message):
        import_aims(write(tmp_path, text), tmp_path / "bad.qvf")
    assert not (tmp_path / "bad.qvf").exists()


@pytest.mark.parametrize("text", [GEOMETRY, MOLECULE])
def test_export_round_trip(tmp_path, text):
    report, geometry, raw, qvf = imported(tmp_path, text)
    out = tmp_path / "geometry.out"
    exported = export_aims(qvf, out)
    assert exported["frozen_atoms"] == report["frozen_atoms"]
    assert exported["periodic"] is report["periodic"]
    again, back, raw_back, _ = imported(tmp_path, out.read_text())
    assert raw_back.get("constraints") == raw.get("constraints")
    assert again["fractional_atoms"] == 0
    if report["periodic"]:
        assert flat(back["lattice_vectors"]) == pytest.approx(flat(geometry["lattice_vectors"]),
                                                              rel=REL, abs=1e-14)
        assert again["cell_volume_angstrom3"] == pytest.approx(report["cell_volume_angstrom3"],
                                                               rel=REL)
    for a, b in zip(back["atoms"], geometry["atoms"]):
        assert a["symbol"] == b["symbol"]
        assert a["position"] == pytest.approx(b["position"], rel=REL, abs=1e-14)


def test_export_refuses_slab(tmp_path):
    path = tmp_path / "slab.qvf"
    with Writer(path, "qvf-test", "1") as writer:
        writer.add_section({"id": "structure", "kind": "structure"})
        writer.add_json("structure", "structure", {
            "atoms": [{"symbol": "H", "atomic_number": 1, "position": [0, 0, 0]}],
            "pbc": [True, True, False], "lattice_vectors": [[3, 0, 0], [0, 3, 0], [0, 0, 3]]})
    with pytest.raises(ValueError, match="periodic on all three axes or a molecule"):
        export_aims(path, tmp_path / "geometry.in")


def ase_read(path):
    aims = pytest.importorskip("ase.io.aims")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with open(path) as stream:
            return aims.read_aims(stream)


@pytest.mark.parametrize("text", [ASE_GEOMETRY, MOLECULE.split("empty")[0]])
def test_independent_ase_reader(tmp_path, text):
    report, geometry, _, qvf = imported(tmp_path, text)
    for path in (write(tmp_path, text, "ase.in"), tmp_path / "export.in"):
        if path.name == "export.in":
            export_aims(qvf, path)
        atoms = ase_read(path)
        assert atoms.get_chemical_symbols() == [a["symbol"] for a in geometry["atoms"]]
        assert atoms.positions.ravel().tolist() == pytest.approx(
            [v for a in geometry["atoms"] for v in a["position"]], rel=REL, abs=1e-14)
        if report["periodic"]:
            assert atoms.pbc.tolist() == [True] * 3
            assert atoms.cell.array.ravel().tolist() == pytest.approx(
                flat(geometry["lattice_vectors"]), rel=REL)
            assert atoms.cell.volume == pytest.approx(report["cell_volume_angstrom3"], rel=REL)
            assert sorted(set(atoms.constraints[0].get_indices().tolist())) == [0, 3]
        else:
            assert not atoms.pbc.any()


def test_independent_ase_writer(tmp_path):
    aims = pytest.importorskip("ase.io.aims")
    report, geometry, _, _ = imported(tmp_path, ASE_GEOMETRY)
    atoms = ase_read(write(tmp_path, ASE_GEOMETRY, "ase.in"))
    for scaled in (False, True):
        path = tmp_path / f"ase-{scaled}.in"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with open(path, "w") as stream:
                aims.write_aims(stream, atoms, scaled=scaled)
        again, back, raw, _ = imported(tmp_path, path.read_text())
        assert raw["constraints"] == {"frozen_atoms": [0, 3]}
        assert again["fractional_atoms"] == (4 if scaled else 0)
        assert flat(back["lattice_vectors"]) == pytest.approx(flat(geometry["lattice_vectors"]),
                                                              rel=REL)
        for a, b in zip(back["atoms"], geometry["atoms"]):
            assert a["position"] == pytest.approx(b["position"], rel=1e-9, abs=1e-9)


def test_cli_aims(tmp_path):
    from qvf.cli import main
    src = write(tmp_path, GEOMETRY)
    qvf = tmp_path / "cli.qvf"
    assert main(["convert", str(src), "-o", str(qvf), "--from", "aims", "--to", "qvf"]) == 0
    assert main(["convert", str(qvf), "-o", str(tmp_path / "out.in"), "--from", "qvf",
                 "--to", "aims"]) == 0
    assert (tmp_path / "out.in").read_text().startswith("# structure\nlattice_vector")
