"""pp.x filplot import (qvf.convert.import_espresso_filplot, WP-48).

The fixtures in fixtures/espresso are written by
scripts/gen_espresso_fixtures.py: each field as a filplot file and as the
pp.x cube of the same data, so the filplot import is checked against the
espresso-cube import and ASE's cube reader.
"""
from pathlib import Path
import struct
import subprocess
import sys

import pytest

from qvf import Reader
from qvf._legacy_writer import BOHR_PER_ANGSTROM
from qvf.convert import import_espresso_cube, import_espresso_filplot
from qvf.validate import run

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).resolve().parent / "fixtures" / "espresso"

# name -> (kind, factor, declared unit, atoms, shape)
CASES = {
    "si_rho": ("volume.density", 1.0, "e_per_bohr3", 2, [6, 6, 6]),
    "gan_v": ("volume.potential", -0.5, "hartree_per_e", 4, [5, 5, 8]),
    "fe_spin": ("volume.spin", 1.0, "e_per_bohr3", 1, [4, 4, 4]),
}


def tensor(reader, sid):
    return [v for (v,) in struct.iter_unpack("<d", reader.tensor(sid, "data").data)]


def flat(rows):
    return [v for row in rows for v in row]


def contents(path):
    with Reader(path) as reader:
        section = next(s for s in reader.raw_manifest["sections"] if s["id"] == "volume")
        return (section, reader.json("volume", "grid"), tensor(reader, "volume"),
                reader.json("structure", "structure"))


def imported(tmp_path, name, convert=import_espresso_filplot, suffix="filplot", **kwargs):
    out = tmp_path / f"{name}-{suffix}.qvf"
    report = convert(FIXTURES / f"{name}.{suffix}", out, **kwargs)
    check = run(out)
    assert check["ok"] and check["findings"] == [] and report["findings"] == []
    return report, out


def test_fixtures_are_fresh():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "gen_espresso_fixtures.py"),
                    "--check"], check=True)


@pytest.mark.parametrize("name", CASES)
def test_filplot_matches_the_pp_cube_of_the_same_data(tmp_path, name):
    kind, factor, unit, atoms, shape = CASES[name]
    report, out = imported(tmp_path, name)
    assert (report["kind"], report["factor"], report["units"]) == (kind, factor, unit)
    assert report["atoms"] == atoms and report["shape"] == shape
    section, grid, data, structure = contents(out)
    assert section["kind"] == kind and section["label"].endswith(f"plot_num {report['plot_num']})")
    assert grid["periodic"] == [True] * 3 and grid["includes_endpoint"] is False
    assert grid["units"] == {"data": unit} and grid["shape"] == shape
    assert flat([[v * n for v in row] for n, row in zip(shape, grid["voxel_vectors"])]) \
        == pytest.approx(flat([[v * BOHR_PER_ANGSTROM for v in row]
                               for row in structure["lattice_vectors"]]), rel=1e-12, abs=1e-12)
    # The cube holds the values as pp.x wrote them, to e13.5 and 6 decimals.
    _, cube = imported(tmp_path, name, import_espresso_cube, "cube", kind="volume.generic")
    _, cube_grid, cube_data, cube_structure = contents(cube)
    assert cube_grid["shape"] == shape and cube_grid["includes_endpoint"] is False
    assert flat(grid["voxel_vectors"]) == pytest.approx(flat(cube_grid["voxel_vectors"]),
                                                        abs=1e-6)
    assert flat(structure["lattice_vectors"]) == pytest.approx(
        flat(cube_structure["lattice_vectors"]), abs=1e-5)
    assert [a["symbol"] for a in structure["atoms"]] == \
        [a["symbol"] for a in cube_structure["atoms"]]
    assert flat(a["position"] for a in structure["atoms"]) == pytest.approx(
        flat(a["position"] for a in cube_structure["atoms"]), abs=1e-5)
    assert data == pytest.approx([v * factor for v in cube_data], rel=1e-5, abs=1e-9)


def test_density_integral_and_header_report(tmp_path):
    report, _ = imported(tmp_path, "si_rho")
    assert report["integral"] == pytest.approx(8.0, rel=1e-8)
    assert report["ibrav"] == 2 and report["alat_angstrom"] == pytest.approx(10.26 / 1.8897261246)
    assert report["species_labels"] == {"Si": "Si"}
    assert report["valence_charges_omitted"] == {"Si": 4.0}
    assert (report["gcutm"], report["dual"], report["ecutwfc_ry"]) == (60.0, 4.0, 30.0)
    spin, _ = imported(tmp_path, "fe_spin")
    assert spin["integral"] == pytest.approx(0.0, abs=1e-12)
    assert spin["species_labels"] == {"Fe": "Fe"}


def test_padding_and_three_digit_exponent(tmp_path):
    report, out = imported(tmp_path, "gan_v")
    assert report["fft_dimensions"] == [6, 7, 8] and "integral" not in report
    _, grid, data, structure = contents(out)
    nx, ny, nz = grid["shape"]
    # The generator writes 1.25e-120 Ry (as 1.250000000-120) at x, y, z = 2, 3, 5.
    assert data[(2 * ny + 3) * nz + 5] == pytest.approx(-0.5 * 1.25e-120, rel=1e-12)
    with Reader(out) as reader:
        labels = {x["id"]: x.get("label") for x in reader.raw_manifest["sections"]}
    assert labels["structure"] == "synthetic GaN potential"
    assert structure["lattice_vectors"][1] == pytest.approx(
        [-0.5 * 6.026 / BOHR_PER_ANGSTROM, 3 ** 0.5 / 2 * 6.026 / BOHR_PER_ANGSTROM, 0.0])


def test_cube_values_against_ase(tmp_path):
    cube = pytest.importorskip("ase.io.cube")
    _, out = imported(tmp_path, "si_rho")
    _, _, data, _ = contents(out)
    with open(FIXTURES / "si_rho.cube") as stream:
        reference = cube.read_cube(stream)
    assert data == pytest.approx(reference["data"].ravel().tolist(), rel=1e-5, abs=1e-9)


def edited(tmp_path, old, new, name="si_rho"):
    text = (FIXTURES / f"{name}.filplot").read_text()
    assert old in text
    path = tmp_path / "edited.filplot"
    path.write_text(text.replace(old, new, 1))
    return path


@pytest.mark.parametrize("plot_num, kind, factor, unit", [
    (1, "volume.generic", 0.5, "hartree"),
    (8, "volume.elf", 1.0, "dimensionless"),
    (9, "volume.difference", 1.0, "e_per_bohr3"),
    (19, "volume.rdg", 1.0, "dimensionless"),
    (7, "volume.generic", 1.0, None),
])
def test_plot_num_table(tmp_path, plot_num, kind, factor, unit):
    path = edited(tmp_path, "30.0000000000     0", f"30.0000000000{plot_num:6d}")
    report = import_espresso_filplot(path, tmp_path / "out.qvf")
    assert report["findings"] == [] and run(tmp_path / "out.qvf")["findings"] == []
    assert (report["kind"], report["factor"], report["units"]) == (kind, factor, unit)
    with Reader(tmp_path / "out.qvf") as reader:
        grid = reader.json("volume", "grid")
    assert grid.get("units") == ({"data": unit} if unit else None)


@pytest.mark.parametrize("old, new, message", [
    ("       6       6       6       6       6       6       2       1",
     "       6       6       6       6       6       6       2", "plot_io layout"),
    ("       6       6       6       6       6       6       2       1",
     "       6       6       6       7       6       6       2       1", "nr <= nrx"),
    ("     2 10.26000000", "     2  0.00000000", "celldm(1) must be positive"),
    ("     2 10.26000000", "    99 10.26000000", "ibrav = 99"),
    ("   1   Si    4.00", "   2   Si    4.00", "species indices"),
    ("   1   Si    4.00", "   1   Xx    4.00", "element symbol"),
    ("  4.962837931E-02", "", "found 215"),
    ("  4.962837931E-02", "  4.96283793QE-02", "not a number"),
])
def test_refusals(tmp_path, old, new, message):
    path = edited(tmp_path, old, new)
    with pytest.raises(ValueError, match=message.replace("(", r"\(").replace(")", r"\)")):
        import_espresso_filplot(path, tmp_path / "bad.qvf")
    assert not (tmp_path / "bad.qvf").exists()


def test_cli_filplot(tmp_path, capsys):
    from qvf.cli import main
    out = tmp_path / "cli.qvf"
    assert main(["convert", str(FIXTURES / "gan_v.filplot"), "-o", str(out), "--from",
                 "espresso-filplot", "--to", "qvf"]) == 0
    assert run(out)["findings"] == []
    capsys.readouterr()
    assert main(["convert", str(out), "-o", str(tmp_path / "x.filplot"), "--from", "qvf",
                 "--to", "espresso-filplot"]) == 2
    assert "import-only" in capsys.readouterr().err
