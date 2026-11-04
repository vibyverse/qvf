"""Quantum ESPRESSO pw.x input, XML data file and pp.x cube import/export (qvf.convert, WP-48)."""
import math
import struct
import warnings

import pytest

from qvf import Reader, Writer
from qvf._legacy_writer import BOHR_PER_ANGSTROM
from qvf.convert import (export_cube, export_espresso_in, import_espresso_cube,
                         import_espresso_in, import_espresso_xml)
from qvf.validate import run

BOHR = 1 / BOHR_PER_ANGSTROM
REL = 1e-10

# Synthetic inputs written for these tests.
SI_FCC = """&CONTROL
  calculation = 'scf', pseudo_dir = './pseudo/'   ! a slash inside a string
  title = 'synthetic Si'
/
&SYSTEM
  ibrav = 2, celldm(1) = 10.26, nat = 2, ntyp = 1,
  ecutwfc = 30.0d0
/
&ELECTRONS
  conv_thr = 1.0d-8
/
ATOMIC_SPECIES
  Si1  28.086  Si.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
  Si1  0.00 0.00 0.00   0 0 0
  Si1  1/4  0.25 (1+1)/8
K_POINTS automatic
  4 4 4 0 0 0
"""

GAN_ANGSTROM = """&control
  title = 'synthetic GaN'
/
&system
  ibrav = 0, nat = 4, ntyp = 2, ecutwfc = 40
/
&electrons
/
ATOMIC_SPECIES
Ga 69.723 Ga.UPF
N  14.007 N.UPF
CELL_PARAMETERS angstrom
   3.1890000000    0.0000000000    0.0000000000
  -1.5945000000    2.7617646000    0.0000000000
   0.0000000000    0.0000000000    5.1850000000
ATOMIC_POSITIONS crystal
Ga  0.3333333333  0.6666666667  0.0000000000  0 0 0
Ga  0.6666666667  0.3333333333  0.5000000000
N   0.3333333333  0.6666666667  0.3770000000  0 0 1
N   0.6666666667  0.3333333333  0.8770000000
K_POINTS gamma
"""


def cubic(cell, positions, system="", species="Na 22.99 Na.UPF\nCl 35.45 Cl.UPF"):
    return (f"&SYSTEM\n ibrav = 0, nat = 2, ntyp = 2{system}\n/\nATOMIC_SPECIES\n{species}\n"
            f"CELL_PARAMETERS {cell[0]}\n" + "\n".join(cell[1:]) +
            f"\nATOMIC_POSITIONS {positions[0]}\n" + "\n".join(positions[1:]) + "\n")


def write(tmp_path, text, name="pw.in"):
    path = tmp_path / name
    path.write_text(text)
    return path


def imported(tmp_path, text, convert=import_espresso_in, name="pw.in"):
    out = tmp_path / f'out{len(list(tmp_path.glob("out*.qvf")))}.qvf'
    report = convert(write(tmp_path, text, name), out)
    check = run(out)
    assert check["ok"] and check["findings"] == [] and report["findings"] == []
    with Reader(out) as reader:
        geometry = reader.json("structure", "structure")
        raw = reader.raw_manifest
    return report, geometry, raw, out


def flat(rows):
    return [v for row in rows for v in row]


def det(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def test_ibrav2_crystal_expressions_and_if_pos(tmp_path):
    report, geometry, raw, _ = imported(tmp_path, SI_FCC)
    a = 10.26 * BOHR
    assert report["ibrav"] == 2 and report["cell_source"] == "celldm"
    assert report["alat_angstrom"] == pytest.approx(a, rel=REL)
    assert report["species_labels"] == {"Si1": "Si"}
    assert report["cards_omitted"] == ["K_POINTS"]
    assert report["frozen_atoms"] == 1 and raw["constraints"] == {"frozen_atoms": [0]}
    assert raw["sections"][0]["label"] == "synthetic Si"
    assert geometry["pbc"] == [True] * 3
    assert flat(geometry["lattice_vectors"]) == pytest.approx(
        flat([[-a / 2, 0, a / 2], [0, a / 2, a / 2], [-a / 2, a / 2, 0]]), rel=REL)
    assert report["cell_volume_angstrom3"] == pytest.approx(a ** 3 / 4, rel=REL)
    assert geometry["atoms"][1]["position"] == pytest.approx([-a / 4, a / 4, a / 4], rel=REL)


@pytest.mark.parametrize("cell, positions, system, factor", [
    (["alat", "1 0 0", "0 1 0", "0 0 1"], ["alat", "Na 0 0 0", "Cl 0.5 0.5 0.5"],
     ", celldm(1) = 10.0", 10.0 * BOHR),
    (["alat", "1 0 0", "0 1 0", "0 0 1"], ["angstrom", "Na 0 0 0", "Cl 2.645886 2.645886 2.645886"],
     ", A = 5.291772", 5.291772),
    (["bohr", "10 0 0", "0 10 0", "0 0 10"], ["bohr", "Na 0 0 0", "Cl 5 5 5"], "", 10 * BOHR),
    (["angstrom", "5.0 0 0", "0 5.0 0", "0 0 5.0"], ["crystal", "Na 0 0 0", "Cl 0.5 0.5 0.5"],
     "", 5.0),
    # alat of positions is the length of the first vector without celldm(1).
    (["bohr", "10 0 0", "0 10 0", "0 0 10"], ["alat", "Na 0 0 0", "Cl 0.5 0.5 0.5"], "",
     10 * BOHR),
    # Deprecated forms: CELL_PARAMETERS without unit (alat here) and positions without unit.
    (["", "1 0 0", "0 1 0", "0 0 1"], ["", "Na 0 0 0", "Cl 0.5 0.5 0.5"], ", celldm(1)=10.0",
     10.0 * BOHR),
])
def test_ibrav0_units(tmp_path, cell, positions, system, factor):
    report, geometry, _, _ = imported(tmp_path, cubic(cell, positions, system))
    assert flat(geometry["lattice_vectors"]) == pytest.approx(
        flat([[factor, 0, 0], [0, factor, 0], [0, 0, factor]]), rel=REL)
    assert geometry["atoms"][1]["position"] == pytest.approx([factor / 2] * 3, rel=1e-6)
    assert report["cell_volume_angstrom3"] == pytest.approx(factor ** 3, rel=REL)
    assert report["cell_units_defaulted"] is (cell[0] == "")
    assert report["positions_units_defaulted"] is (positions[0] == "")


# Conventional cells in the QE orientation: (multiplicity, a, b, c, alpha, beta, gamma).
A, B, C = 4.1, 5.3, 6.7
CASES = {
    1: ("", 1, (A, A, A, 90, 90, 90)),
    2: ("", 4, (A, A, A, 90, 90, 90)),
    3: ("", 2, (A, A, A, 90, 90, 90)),
    -3: ("", 2, (A, A, A, 90, 90, 90)),
    4: (f", C = {C}", 1, (A, A, C, 90, 90, 120)),
    6: (f", C = {C}", 1, (A, A, C, 90, 90, 90)),
    7: (f", C = {C}", 2, (A, A, C, 90, 90, 90)),
    8: (f", B = {B}, C = {C}", 1, (A, B, C, 90, 90, 90)),
    9: (f", B = {B}, C = {C}", 2, (A, B, C, 90, 90, 90)),
    -9: (f", B = {B}, C = {C}", 2, (A, B, C, 90, 90, 90)),
    91: (f", B = {B}, C = {C}", 2, (A, B, C, 90, 90, 90)),
    10: (f", B = {B}, C = {C}", 4, (A, B, C, 90, 90, 90)),
    11: (f", B = {B}, C = {C}", 2, (A, B, C, 90, 90, 90)),
    12: (f", B = {B}, C = {C}, cosAB = 0.3", 1, (A, B, C, 90, 90, math.degrees(math.acos(0.3)))),
    -12: (f", B = {B}, C = {C}, cosAC = -0.2", 1, (A, B, C, 90, math.degrees(math.acos(-0.2)), 90)),
    13: (f", B = {B}, C = {C}, cosAB = 0.3", 2, (A, B, C, 90, 90, math.degrees(math.acos(0.3)))),
    -13: (f", B = {B}, C = {C}, cosAC = -0.2", 2,
          (A, B, C, 90, math.degrees(math.acos(-0.2)), 90)),
    14: (f", B = {B}, C = {C}, cosBC = 0.1, cosAC = -0.2, cosAB = 0.3", 1,
         (A, B, C, math.degrees(math.acos(0.1)), math.degrees(math.acos(-0.2)),
          math.degrees(math.acos(0.3)))),
}


def conventional(a, b, c, alpha, beta, gamma):
    """Rows of a cell with a along x and b in the xy plane (crystallographic standard)."""
    ca, cb, cg = (math.cos(math.radians(v)) for v in (alpha, beta, gamma))
    sg = math.sin(math.radians(gamma))
    cy = (ca - cb * cg) / sg
    return [[a, 0, 0], [b * cg, b * sg, 0], [c * cb, c * cy, c * math.sqrt(1 - cb * cb - cy * cy)]]


def ibrav_input(ibrav, extra):
    return (f"&SYSTEM\n ibrav = {ibrav}, A = {A}{extra}, nat = 1, ntyp = 1\n/\n"
            "ATOMIC_SPECIES\n X1 1.0 x.UPF\nATOMIC_POSITIONS crystal\n X1 0.1 0.2 0.3\n"
            ).replace("X1", "Fe1")


@pytest.mark.parametrize("ibrav", sorted(CASES))
def test_ibrav_lattices_against_conventional_cells(tmp_path, ibrav):
    extra, multiplicity, cellpar = CASES[ibrav]
    _, geometry, _, _ = imported(tmp_path, ibrav_input(ibrav, extra))
    primitive = geometry["lattice_vectors"]
    conv = conventional(*cellpar)
    # Independent check: the conventional vectors are integer combinations of
    # the primitive ones with the centring multiplicity as determinant.
    inverse = [[0.0] * 3 for _ in range(3)]
    d = det(primitive)
    for i in range(3):
        for j in range(3):
            inverse[i][j] = (primitive[(j + 1) % 3][(i + 1) % 3] * primitive[(j + 2) % 3][(i + 2) % 3]
                             - primitive[(j + 1) % 3][(i + 2) % 3]
                             * primitive[(j + 2) % 3][(i + 1) % 3]) / d
    m = [[sum(conv[r][k] * inverse[k][s] for k in range(3)) for s in range(3)] for r in range(3)]
    assert all(abs(v - round(v)) < 1e-9 for row in m for v in row), m
    assert abs(round(det([[round(v) for v in row] for row in m]))) == multiplicity
    assert abs(d) == pytest.approx(abs(det(conv)) / multiplicity, rel=REL)
    if multiplicity == 1:
        assert flat(primitive) == pytest.approx(flat(conv), rel=REL, abs=1e-12)
    position = geometry["atoms"][0]["position"]
    assert position == pytest.approx(
        [sum(f * primitive[i][j] for i, f in enumerate((0.1, 0.2, 0.3))) for j in range(3)],
        rel=REL)


@pytest.mark.parametrize("ibrav", [5, -5])
def test_ibrav_trigonal(tmp_path, ibrav):
    for extra in (", cosAB = 0.25", None):
        text = ibrav_input(ibrav, extra or "")
        if extra is None:  # the same cell from celldm
            text = text.replace(f"A = {A}", f"celldm(1) = {A * BOHR_PER_ANGSTROM!r}, "
                                            "celldm(4) = 0.25")
        _, geometry, _, _ = imported(tmp_path, text)
        rows = geometry["lattice_vectors"]
        for i in range(3):
            assert math.dist(rows[i], [0, 0, 0]) == pytest.approx(A, rel=REL)
            j = (i + 1) % 3
            assert sum(x * y for x, y in zip(rows[i], rows[j])) / A ** 2 == pytest.approx(
                0.25, rel=REL)
        axis = [sum(r[k] for r in rows) for k in range(3)]
        expected = [0, 0, 1] if ibrav == 5 else [1 / math.sqrt(3)] * 3
        norm = math.dist(axis, [0, 0, 0])
        assert [v / norm for v in axis] == pytest.approx(expected, abs=1e-12)
        assert abs(det(rows)) == pytest.approx(
            A ** 3 * math.sqrt(1 - 3 * 0.25 ** 2 + 2 * 0.25 ** 3), rel=REL)


@pytest.mark.parametrize("text, message", [
    (SI_FCC.replace("&SYSTEM", "&SYS"), "requires a &SYSTEM"),
    (SI_FCC.replace("ibrav = 2,", ""), "integer ibrav"),
    (SI_FCC.replace("ibrav = 2", "ibrav = 15"), "not a Bravais lattice index"),
    (SI_FCC.replace("celldm(1) = 10.26,", ""), "requires celldm\\(1\\) or A"),
    (SI_FCC.replace("celldm(1) = 10.26,", "celldm(1) = 10.26, A = 5.4,"), "both celldm and A"),
    (SI_FCC.replace("nat = 2", "nat = 3"), "nat position lines"),
    (SI_FCC.replace("ntyp = 1", "ntyp = 2"), "ntyp species lines"),
    (SI_FCC.replace("{crystal}", "{crystal_sg}"), "crystal_sg"),
    (SI_FCC.replace("{crystal}", "{furlongs}"), "not alat, bohr, angstrom or crystal"),
    (SI_FCC.replace("1/4  0.25", "1/0  0.25"), "arithmetic expression"),
    (SI_FCC.replace("1/4  0.25", "x  0.25"), "number or arithmetic"),
    (SI_FCC.replace("0 0 0\n", "0 2 0\n"), "0 or 1"),
    (SI_FCC.replace("  Si1  1/4", "  Xx1  1/4"), "not in ATOMIC_SPECIES"),
    (SI_FCC.replace("Si1  28.086", "Qq1  28.086").replace("  Si1", "  Qq1"),
     "does not start with an element"),
    (SI_FCC.replace("conv_thr = 1.0d-8\n/", "conv_thr = 1.0d-8\n"), "closing"),
    (SI_FCC.replace("ecutwfc = 30.0d0", "ecutwfc = thirty"), "not a number"),
    (SI_FCC.replace("ibrav = 2,", "ibrav 2,"), "key = value"),
    (SI_FCC.replace("K_POINTS", "HELLO\nK_POINTS").replace("ATOMIC_SPECIES", "stray\nATOMIC_SPECIES"),
     "neither a namelist nor a card"),
    (cubic(["bohr", "10 0 0", "0 10 0", "0 0 10"], ["bohr", "Na 0 0 0", "Cl 5 5 5"],
           ", celldm(1) = 10"), "given twice"),
    (cubic(["alat", "1 0 0", "0 1 0", "0 0 1"], ["bohr", "Na 0 0 0", "Cl 5 5 5"]),
     "alat requires celldm"),
    (cubic(["bohr", "10 0 0", "20 0 0", "0 0 10"], ["bohr", "Na 0 0 0", "Cl 5 5 5"]),
     "linearly dependent"),
    (cubic(["bohr", "10 0 0", "0 10 0"], ["bohr", "Na 0 0 0", "Cl 5 5 5"]), "three vectors"),
    (SI_FCC.replace("ibrav = 2", "ibrav = 0"), "requires a CELL_PARAMETERS"),
    (GAN_ANGSTROM.replace("ibrav = 0", "ibrav = 4, A = 3.189"), "requires ibrav = 0"),
    (ibrav_input(14, ", B = 1, C = 1, cosBC = 0.9, cosAC = -0.9, cosAB = 0.9"),
     "do not describe a cell"),
    (ibrav_input(8, ", C = 2"), "positive b/a"),
])
def test_input_refusals(tmp_path, text, message):
    with pytest.raises(ValueError, match=message):
        import_espresso_in(write(tmp_path, text), tmp_path / "bad.qvf")
    assert not (tmp_path / "bad.qvf").exists()


def test_export_round_trip(tmp_path):
    for text in (SI_FCC, GAN_ANGSTROM):
        report, geometry, raw, qvf = imported(tmp_path, text)
        out = tmp_path / "export.in"
        out.unlink(missing_ok=True)
        exported = export_espresso_in(qvf, out)
        assert exported["frozen_atoms"] == report["frozen_atoms"]
        again, back, raw_back, _ = imported(tmp_path, out.read_text())
        assert again["ibrav"] == 0 and again["cell_units"] == "angstrom"
        assert raw_back.get("constraints") == raw.get("constraints")
        assert flat(back["lattice_vectors"]) == pytest.approx(flat(geometry["lattice_vectors"]),
                                                              rel=REL, abs=1e-14)
        assert again["cell_volume_angstrom3"] == pytest.approx(report["cell_volume_angstrom3"],
                                                               rel=REL)
        for a, b in zip(back["atoms"], geometry["atoms"]):
            assert a["symbol"] == b["symbol"]
            assert a["position"] == pytest.approx(b["position"], rel=REL, abs=1e-14)
    assert GAN_ANGSTROM and exported["species"] == ["Ga", "N"]


def test_export_refuses_slab_and_reports_constraints(tmp_path):
    path = tmp_path / "slab.qvf"
    with Writer(path, "qvf-test", "1") as writer:
        writer.add_section({"id": "structure", "kind": "structure"})
        writer.add_json("structure", "structure", {
            "atoms": [{"symbol": "H", "atomic_number": 1, "position": [0, 0, 0]}],
            "pbc": [True, True, False], "lattice_vectors": [[3, 0, 0], [0, 3, 0], [0, 0, 3]]})
    with pytest.raises(ValueError, match="periodic on all three axes"):
        export_espresso_in(path, tmp_path / "slab.in")
    path = tmp_path / "dist.qvf"
    with Writer(path, "qvf-test", "1") as writer:
        writer.add_section({"id": "structure", "kind": "structure"})
        writer.add_json("structure", "structure", {
            "atoms": [{"symbol": "H", "atomic_number": 1, "position": [0, 0, 0]},
                      {"symbol": "H", "atomic_number": 1, "position": [0.74, 0, 0]}],
            "pbc": [True] * 3, "lattice_vectors": [[5, 0, 0], [0, 5, 0], [0, 0, 5]]})
        writer.set_constraints(distance_constraints=[{"atoms": [0, 1], "target_angstrom": 0.74}])
    report = export_espresso_in(path, tmp_path / "dist.in")
    assert report["constraints_not_written"] == ["distance_constraints"]


def ase_bohr_ratio():
    from ase.io.espresso import units
    return BOHR / units["Bohr"]


@pytest.mark.parametrize("text", [
    GAN_ANGSTROM,
    cubic(["alat", "1 0 0", "0 1 0.1", "0 0 1"], ["alat", "Na 0 0 0", "Cl 0.5 0.5 0.5"],
          ", celldm(1) = 10.0"),
    cubic(["bohr", "10 0 0", "0 10 0", "1 0 10"], ["bohr", "Na 0 0 0", "Cl 5 5 5"]),
    cubic(["angstrom", "5 0 0", "0 5 0", "0 0 5"], ["angstrom", "Na 0 0 0", "Cl 2.5 2.5 2.5"]),
    cubic(["angstrom", "5 0 0", "0 5 0", "0 0 5"], ["crystal", "Na 0 0 0", "Cl 1/2 0.5 0.5"]),
])
def test_independent_ase_reader(tmp_path, text):
    espresso = pytest.importorskip("ase.io.espresso")
    report, geometry, _, _ = imported(tmp_path, text)
    with open(write(tmp_path, text, "ase.in")) as stream:
        atoms = espresso.read_espresso_in(stream)
    # ASE converts bohr with the CODATA 2006 value of QE; QVF uses its own
    # constant, so values in bohr or alat differ by that ratio only.
    ratio = 1.0 if "angstrom" in text.split("CELL_PARAMETERS")[1].split("\n")[0] else \
        ase_bohr_ratio()
    assert flat(geometry["lattice_vectors"]) == pytest.approx(
        [v * ratio for v in atoms.cell.array.ravel()], rel=REL, abs=1e-14)
    assert report["cell_volume_angstrom3"] == pytest.approx(atoms.cell.volume * ratio ** 3,
                                                            rel=REL)
    units = text.split("ATOMIC_POSITIONS")[1].split("\n")[0]
    pratio = 1.0 if "angstrom" in units or "crystal" in units and ratio == 1.0 else ratio
    assert [v for a in geometry["atoms"] for v in a["position"]] == pytest.approx(
        [v * pratio for v in atoms.positions.ravel()], rel=REL, abs=1e-14)
    assert [a["symbol"] for a in geometry["atoms"]] == atoms.get_chemical_symbols()


def test_independent_ase_reads_export_and_writes_input(tmp_path):
    espresso = pytest.importorskip("ase.io.espresso")
    _, geometry, _, qvf = imported(tmp_path, GAN_ANGSTROM)
    out = tmp_path / "export.in"
    export_espresso_in(qvf, out)
    with open(out) as stream:
        atoms = espresso.read_espresso_in(stream)
    assert atoms.get_chemical_symbols() == ["Ga", "Ga", "N", "N"]
    assert atoms.cell.array.ravel().tolist() == pytest.approx(
        flat(geometry["lattice_vectors"]), rel=REL)
    assert atoms.positions.ravel().tolist() == pytest.approx(
        [v for a in geometry["atoms"] for v in a["position"]], rel=REL, abs=1e-14)
    assert [c.get_indices().tolist() for c in atoms.constraints] == [[0]]
    written = tmp_path / "ase.pwi"
    with open(written, "w") as stream:
        espresso.write_espresso_in(stream, atoms, input_data={"system": {"ecutwfc": 40}},
                                   pseudopotentials={"Ga": "Ga.UPF", "N": "N.UPF"},
                                   kpts=(2, 2, 2))
    report, back, raw, _ = imported(tmp_path, written.read_text())
    assert report["frozen_atoms"] == 1 and raw["constraints"] == {"frozen_atoms": [0]}
    assert "K_POINTS" in report["cards_omitted"]
    assert flat(back["lattice_vectors"]) == pytest.approx(flat(geometry["lattice_vectors"]),
                                                          rel=REL)
    # ASE writes positions with 10 decimals.
    for a, b in zip(back["atoms"], geometry["atoms"]):
        assert a["position"] == pytest.approx(b["position"], abs=1e-9)


# -- XML data file ----------------------------------------------------------

def xml(scf="true", positions=None, extra_output="", stress=True, free="1 1 1 0 0 0"):
    positions = positions or ('<atomic_positions>\n'
                              '        <atom name="Si" index="1">0.0 0.0 0.0</atom>\n'
                              '        <atom name="Si" index="2">2.565 2.565 2.565</atom>\n'
                              '      </atomic_positions>')
    stress_text = ('<stress rank="2" dims="3 3" order="F">\n 1.0e-4 2.0e-6 0.0\n'
                   ' 2.0e-6 1.0e-4 0.0\n 0.0 0.0 1.0e-4\n </stress>') if stress else ""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0" Units="Hartree atomic units">
  <general_info><creator NAME="PWSCF" VERSION="synthetic">XML file generated by PWSCF</creator></general_info>
  <input>
    <free_positions rank="2" dims="3 2" order="F">{free}</free_positions>
  </input>
  <output>
    <convergence_info>
      <scf_conv><convergence_achieved>{scf}</convergence_achieved><n_scf_steps>9</n_scf_steps></scf_conv>
    </convergence_info>
    <atomic_species ntyp="1"><species name="Si"><mass>28.086</mass><pseudo_file>Si.UPF</pseudo_file></species></atomic_species>
    <atomic_structure nat="2" alat="10.26" bravais_index="2">
      {positions}
      <cell>
        <a1>-5.13 0.0 5.13</a1>
        <a2>0.0 5.13 5.13</a2>
        <a3>-5.13 5.13 0.0</a3>
      </cell>
    </atomic_structure>
    <dft><functional>PBE</functional></dft>
    <total_energy><etot>-1.5850000000000e1</etot><eband>1.0</eband></total_energy>
    <band_structure><lsda>false</lsda></band_structure>
    <forces rank="2" dims="3 2" order="F">
 1.0e-3 -2.0e-3 3.0e-3
 -1.0e-3 2.0e-3 -3.0e-3
    </forces>
    {stress_text}{extra_output}
  </output>
  <exit_status>0</exit_status>
  <closed DATE="synthetic" TIME="synthetic"></closed>
</qes:espresso>
"""


def xml_import(tmp_path, text):
    out = tmp_path / f'xml{len(list(tmp_path.glob("xml*.qvf")))}.qvf'
    report = import_espresso_xml(write(tmp_path, text, "data-file-schema.xml"), out)
    check = run(out)
    assert check["ok"] and check["findings"] == [] and report["findings"] == []
    return report, out


def tensor(reader, sid, role):
    return [v for (v,) in struct.iter_unpack("<d", reader.tensor(sid, role).data)]


def test_xml_structure_energy_forces_stress(tmp_path):
    report, out = xml_import(tmp_path, xml())
    assert report["scf_converged"] is True and report["exit_status"] == "0"
    assert report["functional"] == "PBE" and report["bravais_index"] == "2"
    assert report["energy_hartree"] == -15.85 and report["forces"] and report["stress"]
    assert report["frozen_atoms"] == 1 and report["output_elements_omitted"] == [
        "atomic_species", "band_structure", "dft"]
    a = 10.26 * BOHR
    assert report["cell_volume_angstrom3"] == pytest.approx(a ** 3 / 4, rel=REL)
    with Reader(out) as reader:
        geometry = reader.json("structure", "structure")
        raw = reader.raw_manifest
        assert raw["energies"]["total"] == {"value": -15.85, "units": "hartree"}
        assert raw["constraints"] == {"frozen_atoms": [1]}
        assert tensor(reader, "forces", "gradient") == [-1e-3, 2e-3, -3e-3, 1e-3, -2e-3, 3e-3]
        stress = tensor(reader, "forces", "stress")
    assert stress == [-1e-4, -2e-6, 0, -2e-6, -1e-4, 0, 0, 0, -1e-4]
    assert flat(geometry["lattice_vectors"]) == pytest.approx(
        flat([[-a / 2, 0, a / 2], [0, a / 2, a / 2], [-a / 2, a / 2, 0]]), rel=REL)
    assert geometry["atoms"][1]["position"] == pytest.approx([2.565 * BOHR] * 3, rel=REL)


def test_xml_crystal_positions_and_unconverged(tmp_path):
    crystal = ('<crystal_positions>\n <atom name="Si1" index="1">0 0 0</atom>\n'
               ' <atom name="Si1" index="2">0.25 0.25 0.25</atom>\n</crystal_positions>')
    report, out = xml_import(tmp_path, xml(positions=crystal, scf="false"))
    assert report["results_omitted_unconverged"] is True
    assert report["energy_hartree"] is None and report["forces"] is False
    with Reader(out) as reader:
        assert "energies" not in reader.raw_manifest
        assert [s.kind for s in reader.sections] == ["structure"]
        geometry = reader.json("structure", "structure")
    a = 10.26 * BOHR
    assert geometry["atoms"][1]["position"] == pytest.approx([-a / 4, a / 4, a / 4], rel=REL)


@pytest.mark.parametrize("text, message", [
    ('<?xml version="1.0"?><!DOCTYPE x [<!ENTITY a "b">]><qes:espresso/>', "DOCTYPE"),
    ("<espresso><output>", "not well-formed"),
    ("<other/>", "qes:espresso root"),
    (xml().replace('nat="2"', 'nat="3"'), "one <atom> per atom"),
    (xml().replace("<a3>-5.13 5.13 0.0</a3>", "<a3>-5.13 5.13</a3>"), "a3"),
    (xml().replace("<a3>-5.13 5.13 0.0</a3>", "<a3>-5.13 0.0 5.13</a3>"), "linearly dependent"),
    (xml().replace("-1.0e-3 2.0e-3 -3.0e-3", ""), "forces"),
    (xml().replace(" 2.0e-6 1.0e-4", " 9.0e-6 1.0e-4"), "not symmetric"),
    (xml(free="1 1 1 0 2 0"), "free_positions"),
    (xml(scf="maybe"), "true or false"),
    (xml().replace('name="Si" index="2"', 'name="Qq" index="2"'), "element symbol"),
])
def test_xml_refusals(tmp_path, text, message):
    with pytest.raises(ValueError, match=message):
        import_espresso_xml(write(tmp_path, text, "bad.xml"), tmp_path / "bad.qvf")
    assert not (tmp_path / "bad.qvf").exists()


# -- pp.x cube --------------------------------------------------------------

def pp_cube(shape=(4, 5, 6)):
    """A cube laid out as pp.x writes it: FFT voxels at/nr (no endpoint plane)."""
    lattice = [[-5.13, 0.0, 5.13], [0.0, 5.13, 5.13], [-5.13, 5.13, 0.0]]
    lines = [" Cubefile created from PWScf calculation",
             " Contains the selected quantity on a FFT grid",
             "    2    0.000000    0.000000    0.000000"]
    for n, row in zip(shape, lattice):
        lines.append(f"{n:5d}" + "".join(f"{v / n:12.6f}" for v in row))
    lines += ["   14    4.000000    0.000000    0.000000    0.000000",
              "   14    4.000000   -2.565000    2.565000    2.565000"]
    values = [0.01 * (1 + i % 7) for i in range(math.prod(shape))]
    for i in range(0, len(values), shape[2]):
        lines.append("".join(f"{v:13.5E}" for v in values[i:i + shape[2]]))
    return "\n".join(lines) + "\n", lattice, values


def test_pp_cube_is_periodic_without_endpoint(tmp_path):
    text, lattice, values = pp_cube()
    out = tmp_path / "rho.qvf"
    report = import_espresso_cube(write(tmp_path, text, "rho.cube"), out)
    check = run(out)
    assert check["ok"] and check["findings"] == [] and report["findings"] == []
    assert report["kind"] == "volume.density" and report["includes_endpoint"] is False
    assert report["nonzero_atom_charges_omitted"] == 2
    volume_bohr = abs(det(lattice))
    assert report["cell_volume_angstrom3"] == pytest.approx(volume_bohr * BOHR ** 3, rel=1e-6)
    assert report["integral"] == pytest.approx(sum(values) / len(values) * volume_bohr, rel=1e-5)
    with Reader(out) as reader:
        grid = reader.json("volume", "grid")
        geometry = reader.json("structure", "structure")
        stored = tensor(reader, "volume", "data")
    assert grid["periodic"] == [True] * 3 and grid["includes_endpoint"] is False
    assert flat([[v * n for v in row] for n, row in zip(grid["shape"], grid["voxel_vectors"])]) \
        == pytest.approx(flat(lattice), abs=1e-5)
    assert flat(geometry["lattice_vectors"]) == pytest.approx(
        flat([[v * BOHR_PER_ANGSTROM ** -1 for v in row] for row in lattice]), rel=1e-5)
    assert stored == pytest.approx(values, rel=1e-5)
    # export_cube writes the same voxels; importing that again keeps the grid.
    export_cube(out, tmp_path / "again.cube")
    second = import_espresso_cube(tmp_path / "again.cube", tmp_path / "again.qvf",
                                  kind="volume.spin")
    assert second["kind"] == "volume.spin"
    with Reader(tmp_path / "again.qvf") as reader:
        assert reader.json("volume", "grid") == grid
        assert tensor(reader, "volume", "data") == stored


def test_pp_cube_against_ase(tmp_path):
    cube = pytest.importorskip("ase.io.cube")
    text, lattice, values = pp_cube()
    path = write(tmp_path, text, "rho.cube")
    import_espresso_cube(path, tmp_path / "rho.qvf")
    with open(path) as stream:
        data = cube.read_cube(stream)
    from ase.units import Bohr
    with Reader(tmp_path / "rho.qvf") as reader:
        geometry = reader.json("structure", "structure")
        stored = tensor(reader, "volume", "data")
    ratio = BOHR / Bohr
    assert flat(geometry["lattice_vectors"]) == pytest.approx(
        [v * ratio for v in data["atoms"].cell.array.ravel()], rel=REL, abs=1e-14)
    assert stored == pytest.approx(data["data"].ravel().tolist(), rel=REL)


def test_cli_espresso(tmp_path, capsys):
    from qvf.cli import main
    src = write(tmp_path, GAN_ANGSTROM)
    qvf = tmp_path / "cli.qvf"
    assert main(["convert", str(src), "-o", str(qvf), "--from", "espresso-in", "--to", "qvf"]) == 0
    assert main(["convert", str(qvf), "-o", str(tmp_path / "out.in"), "--from", "qvf",
                 "--to", "espresso-in"]) == 0
    xml_path = write(tmp_path, xml(), "data-file-schema.xml")
    assert main(["convert", str(xml_path), "-o", str(tmp_path / "x.qvf"), "--from",
                 "espresso-xml", "--to", "qvf"]) == 0
    cube_path = write(tmp_path, pp_cube()[0], "rho.cube")
    assert main(["convert", str(cube_path), "-o", str(tmp_path / "c.qvf"), "--from",
                 "espresso-cube", "--to", "qvf", "--kind", "volume.density"]) == 0
    capsys.readouterr()
    assert main(["convert", str(qvf), "-o", str(tmp_path / "x.xml"), "--from", "qvf",
                 "--to", "espresso-xml"]) == 2
    assert "import-only" in capsys.readouterr().err
    assert main(["convert", str(cube_path), "-o", str(tmp_path / "d.qvf"), "--from",
                 "espresso-cube", "--to", "qvf", "--periodic", "true", "true", "true"]) == 2
    assert "always periodic" in capsys.readouterr().err
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert main(["convert", str(src), "-o", str(tmp_path / "e.qvf"), "--from", "espresso-in",
                     "--to", "qvf", "--species", "Ga"]) == 2
