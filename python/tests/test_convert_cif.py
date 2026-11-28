"""Core crystallographic CIF import/export (qvf.convert, issue #34)."""
import math

import pytest

from qvf import Reader, Writer
from qvf.convert import export_cif, import_cif
from qvf.validate import run

# Rutile TiO2, P 42/m n m, with all 16 operations listed (synthetic values).
RUTILE = """# synthetic rutile example
data_rutile
_cell_length_a    4.5937(3)
_cell_length_b    4.5937(3)
_cell_length_c    2.9587(2)
_cell_angle_alpha 90
_cell_angle_beta  90.0
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P 42/m n m'
_space_group_IT_number 136
_journal_name_full
;
A multi-line text field
with 'quotes' that are not values
;
loop_
_space_group_symop_id
_space_group_symop_operation_xyz
1 'x, y, z'
2 '-x, -y, z'
3 '-y+1/2, x+1/2, z+1/2'
4 'y+1/2, -x+1/2, z+1/2'
5 '-x+1/2, y+1/2, -z+1/2'
6 'x+1/2, -y+1/2, -z+1/2'
7 'y, x, -z'
8 '-y, -x, -z'
9 '-x, -y, -z'
10 'x, y, -z'
11 'y+1/2, -x+1/2, -z+1/2'
12 '-y+1/2, x+1/2, -z+1/2'
13 'x+1/2, -y+1/2, z+1/2'
14 '-x+1/2, y+1/2, z+1/2'
15 '-y, -x, z'
16 'y, x, z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Ti1 Ti4+ 0 0 0 1 0.0050(2)
O1 O2- 0.30478(6) 0.30478(6) 0 1.0 ?
"""

# Rock salt from its Hall symbol only (192 operations).
NACL = """data_nacl
_cell_length_a 5.64
_cell_length_b 5.64
_cell_length_c 5.64
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_Hall '-F 4 2 3'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na Na 0 0 0
Cl Cl 0.5 0.5 0.5
"""


def imported(tmp_path, text, name='in.cif', **kw):
    src = tmp_path/name
    src.write_text(text)
    out = tmp_path/(name + '.qvf')
    report = import_cif(src, out, **kw)
    assert run(out)['ok']
    with Reader(out) as reader:
        geometry = reader.json('structure', 'structure')
        symmetry = (reader.json('symmetry', 'data')
                    if any(s.id == 'symmetry' for s in reader.sections) else None)
    return report, geometry, symmetry, out


def fractional(geometry):
    rows = geometry['lattice_vectors']
    return sorted((a['symbol'], *(round(p / rows[i][i], 5) % 1 for i, p in enumerate(a['position'])))
                  for a in geometry['atoms'])


def test_listed_operations_build_the_unit_cell(tmp_path):
    report, geometry, symmetry, _ = imported(tmp_path, RUTILE)
    assert report['block'] == 'rutile' and report['operator_source'] == 'listed'
    assert (report['sites'], report['operators'], report['atoms']) == (2, 16, 6)
    assert report['special_position_atoms_merged'] == 16 * 2 - 6
    assert report['charges_omitted'] == 2 and report['partial_occupancy_sites'] == 0
    assert report['unused_tags'] == 2  # the symop id and the journal name
    assert symmetry == {'space_group_symbol': 'P 42/m n m', 'space_group_number': 136}
    assert geometry['pbc'] == [True] * 3
    a, c = 4.5937, 2.9587
    lattice = geometry['lattice_vectors']
    assert lattice[0] == [a, 0, 0] and abs(lattice[2][2] - c) < 1e-12
    ti = [x for x in geometry['atoms'] if x['symbol'] == 'Ti']
    assert ti[0] == {'symbol': 'Ti', 'atomic_number': 22, 'position': [0, 0, 0],
                     'atom_name': 'Ti1', 'occupancy': 1,
                     'b_factor': pytest.approx(8 * math.pi ** 2 * 0.005),
                     'adp': {'type': 'Uiso', 'iso': 0.005}}
    assert fractional(geometry) == [
        ('O', 0.19522, 0.80478, 0.5), ('O', 0.30478, 0.30478, 0.0),
        ('O', 0.69522, 0.69522, 0.0), ('O', 0.80478, 0.19522, 0.5),
        ('Ti', 0.0, 0.0, 0.0), ('Ti', 0.5, 0.5, 0.5)]
    assert all(0 <= v < 1 for x in fractional(geometry) for v in x[1:])


def test_hall_symbol_builds_the_unit_cell(tmp_path):
    report, geometry, symmetry, _ = imported(tmp_path, NACL)
    assert report['operator_source'] == 'hall' and report['operators'] == 192
    assert report['atoms'] == 8 and report['special_position_atoms_merged'] == 2 * 192 - 8
    assert symmetry == {'hall': '-F 4 2 3'}
    assert fractional(geometry)[:4] == [('Cl', 0.0, 0.0, 0.5), ('Cl', 0.0, 0.5, 0.0),
                                        ('Cl', 0.5, 0.0, 0.0), ('Cl', 0.5, 0.5, 0.5)]


def test_p1_needs_no_operations_and_blocks_are_selected(tmp_path):
    p1 = ("data_first\n_cell_length_a 3\n_cell_length_b 4\n_cell_length_c 5\n"
          "_cell_angle_alpha 90\n_cell_angle_beta 100\n_cell_angle_gamma 90\n"
          "_symmetry_space_group_name_H-M 'P1'\nloop_\n_atom_site_label\n"
          "_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n"
          "_atom_site_fract_z\n_atom_site_occupancy\nC1 C 0.1 0.2 0.3 0.5\n"
          "C2 C 0.1 0.2 0.3 0.5\n")
    text = p1 + p1.replace('data_first', 'data_second') + 'data_meta\n_audit_creation_date 2026\n'
    with pytest.raises(ValueError, match='select a CIF data block'):
        import_cif(_write(tmp_path, text), tmp_path/'x.qvf')
    report, geometry, symmetry, _ = imported(tmp_path, text, block='SECOND')
    assert report['block'] == 'second' and report['blocks'] == ['first', 'second', 'meta']
    assert report['operator_source'] == 'p1' and report['atoms'] == 2
    # Two sites at one place are different sites (disorder), not merged.
    assert report['special_position_atoms_merged'] == 0
    assert report['partial_occupancy_sites'] == 2
    assert symmetry == {'space_group_symbol': 'P1'}
    with pytest.raises(ValueError, match='has no atom sites'):
        import_cif(_write(tmp_path, text), tmp_path/'y.qvf', block='meta')
    with pytest.raises(ValueError, match="no data block 'other'"):
        import_cif(_write(tmp_path, text), tmp_path/'z.qvf', block='other')


def _write(tmp_path, text, name='case.cif'):
    path = tmp_path/name
    path.write_text(text)
    return path


SITES = ("loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n"
         "_atom_site_fract_y\n_atom_site_fract_z\n")
CELL = ("data_x\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n_cell_angle_alpha 90\n"
        "_cell_angle_beta 90\n_cell_angle_gamma 90\n")


@pytest.mark.parametrize('text, message', [
    (CELL + "_space_group_IT_number 2\n" + SITES + "C1 C 0 0 0\n", 'not P 1'),
    (CELL + "loop_\n_symmetry_equiv_pos_as_xyz\n-x,-y,-z\n" + SITES + "C1 C 0 0 0\n",
     'include x,y,z'),
    (CELL + "loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\n-x,y,z+1/4\n" + SITES + "C1 C 0 0 0\n",
     'do not form a group'),
    (CELL + "loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\n2x,y,z\n" + SITES + "C1 C 0 0 0\n",
     'malformed'),
    (CELL + "loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\nx,x,z\n" + SITES + "C1 C 0 0 0\n",
     'not a lattice symmetry'),
    (CELL + "_symmetry_space_group_name_Hall 'Q 9'\n" + SITES + "C1 C 0 0 0\n",
     'Hall symbol'),
    (CELL + "_space_group_IT_number 1\n" + SITES + "C1 Xx 0 0 0\n", 'type_symbol'),
    (CELL + "_space_group_IT_number 1\n" + SITES + "C1 ? 0 0 0\n", 'type_symbol'),
    (CELL + "_space_group_IT_number 1\nloop_\n_atom_site_label\n_atom_site_fract_x\n"
     "_atom_site_fract_y\n_atom_site_fract_z\nC1 0 0 0\n", 'type_symbol'),
    (CELL + "_space_group_IT_number 1\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
     "_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\nC1 C 0 0 0\n",
     'Cartesian sites'),
    (CELL + "_space_group_IT_number 1\n" + SITES + "C1 C 0 0\n", 'whole rows'),
    (CELL + "_space_group_IT_number 1\n" + SITES + "C1 C 0 0 1.2.3\n", 'requires a number'),
    (CELL.replace('_cell_angle_gamma 90', '_cell_angle_gamma 190') + SITES + "C1 C 0 0 0\n",
     'angles in'),
    (CELL + "_cell_length_a 5\n" + SITES + "C1 C 0 0 0\n", 'repeats a tag'),
    (CELL + "_space_group_IT_number 1\n_title 'unclosed\n" + SITES + "C1 C 0 0 0\n",
     'not closed'),
    (CELL + "_title\n;\nopen text\n", 'not closed'),
    ("_cell_length_a 4\n" + CELL, 'before its first data_'),
    (CELL + "save_frame\n", 'not supported'),
    (CELL + "_space_group_IT_number 1\n" + SITES.replace('loop_', 'loop_\n_atom_site_occupancy')
     + "2 C1 C 0 0 0\n", r'occupancy must lie in \(0, 1\]'),
])
def test_refusals(tmp_path, text, message):
    with pytest.raises(ValueError, match=message):
        import_cif(_write(tmp_path, text), tmp_path/'bad.qvf')
    assert not (tmp_path/'bad.qvf').exists()


def test_export_writes_a_p1_cif_that_imports_back(tmp_path):
    _, geometry, _, qvf = imported(tmp_path, RUTILE)
    out = tmp_path/'p1.cif'
    report = export_cif(qvf, out)
    assert report['format'] == 'cif' and report['atoms'] == 6
    assert report['labels'] == 'generated'  # atom names repeat across copies
    assert report['max_rounding_error_angstrom'] < 1e-7
    assert report['omitted_sections'] == ['symmetry']
    text = out.read_text()
    assert "_space_group_symop_operation_xyz\nx,y,z\n" in text
    assert "\nTi1 Ti 0.00000000 0.00000000 0.00000000 1.0000 0.3948 0.005000 Uiso\n" in text
    again_report, again, symmetry, _ = imported(tmp_path, text, name='again.cif')
    assert again_report['operator_source'] == 'listed' and again_report['atoms'] == 6
    assert symmetry == {'space_group_symbol': 'P 1', 'space_group_number': 1}
    for x, y in zip(again['atoms'], geometry['atoms']):
        assert x['position'] == pytest.approx(y['position'], abs=1e-6)
        assert x['symbol'] == y['symbol'] and x['occupancy'] == y['occupancy']


def test_export_converts_any_lattice_orientation(tmp_path):
    rotated = [[0, 3, 0], [-4, 0, 0], [0, 0, 5]]
    path = tmp_path/'rotated.qvf'
    with Writer(path, 'qvf-test', '1') as writer:
        writer.add_section({'id': 'structure', 'kind': 'structure'})
        writer.add_json('structure', 'structure', {
            'atoms': [{'symbol': 'Si', 'atomic_number': 14, 'position': [-1, 0.75, 2.5],
                       'atom_name': 'Si_a'}],
            'pbc': [True] * 3, 'lattice_vectors': rotated})
    out = tmp_path/'rotated.cif'
    report = export_cif(path, out)
    assert report['labels'] == 'atom_name'
    text = out.read_text()
    assert "_cell_length_a      3.000000" in text and "_cell_angle_gamma   90.000000" in text
    assert "\nSi_a Si 0.25000000 0.25000000 0.50000000\n" in text


def test_export_refuses_nonperiodic_structures(tmp_path):
    path = tmp_path/'molecule.qvf'
    with Writer(path, 'qvf-test', '1') as writer:
        writer.add_section({'id': 'structure', 'kind': 'structure'})
        writer.add_json('structure', 'structure', {
            'atoms': [{'symbol': 'H', 'atomic_number': 1, 'position': [0, 0, 0]}],
            'pbc': [False] * 3, 'lattice_vectors': None})
    with pytest.raises(ValueError, match='periodic on all three axes'):
        export_cif(path, tmp_path/'molecule.cif')


@pytest.mark.parametrize('text', [RUTILE, NACL])
def test_unit_cell_matches_gemmi(tmp_path, text):
    gemmi = pytest.importorskip('gemmi')
    src = _write(tmp_path, text)
    _, geometry, _, _ = imported(tmp_path, text, name='ref.cif')
    small = gemmi.read_small_structure(str(src))
    theirs = sorted((s.element.name, *(round(v, 5) % 1 for v in s.fract.tolist()))
                    for s in small.get_all_unit_cell_sites())
    assert fractional(geometry) == theirs


def test_gemmi_reads_the_export(tmp_path):
    gemmi = pytest.importorskip('gemmi')
    _, geometry, _, qvf = imported(tmp_path, RUTILE)
    out = tmp_path/'p1.cif'; export_cif(qvf, out)
    small = gemmi.read_small_structure(str(out))
    assert small.cell.a == pytest.approx(4.5937) and small.cell.c == pytest.approx(2.9587)
    assert len(small.get_all_unit_cell_sites()) == 6
    assert sorted(s.element.name for s in small.sites) == ['O'] * 4 + ['Ti'] * 2


def test_cli_cif_round_trip(tmp_path, capsys):
    from qvf.cli import main
    src = _write(tmp_path, RUTILE)
    qvf = tmp_path/'cli.qvf'
    assert main(['convert', str(src), '-o', str(qvf), '--from', 'cif', '--to', 'qvf',
                 '--block', 'rutile']) == 0
    assert main(['convert', str(qvf), '-o', str(tmp_path/'cli.cif'), '--from', 'qvf',
                 '--to', 'cif']) == 0
    capsys.readouterr()
    assert main(['convert', str(src), '-o', str(tmp_path/'x.qvf'), '--from', 'xyz', '--to',
                 'qvf', '--block', 'rutile']) == 2
    assert '--block applies only to CIF/mmCIF import' in capsys.readouterr().err


def test_hall_parser_matches_gemmi_for_every_setting():
    gemmi = pytest.importorskip('gemmi')
    from fractions import Fraction
    from qvf._spacegroups import operators
    for sg in gemmi.spacegroup_table():
        ours = operators(sg.hall.strip())
        theirs = {(tuple(tuple(v // op.DEN for v in row) for row in op.rot),
                   tuple(Fraction(v, op.DEN) % 1 for v in op.tran)) for op in sg.operations()}
        assert len(ours) == len(theirs) and set(ours) == theirs, sg.hall
