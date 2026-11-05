"""Gaussian fchk: AO order and normalization, refusals, and independent checks.

Synthetic fixtures come from fchk_fixtures/generate.py (written from the
documented layout, independent of qvf._fchk). PySCF calculations are exported
with qvf's exporter and read back by cclib and IOData when installed, and by
the importer, whose archives are evaluated with an Appendix A evaluator.
"""
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from qvf import Reader, Writer
from qvf.convert import _binary, export_fchk, import_fchk
from qvf._fchk import _layout, _real
from qvf.validate import run
from fchk_fixtures import generate
from test_convert_molden import qvf_orbitals

FIXTURES = HERE / 'fchk_fixtures'
NAMES = list(generate.FIXTURES)


def clean(path):
    report = run(path)
    assert report['ok'] and not report['findings'], report


def fixture(name):
    return FIXTURES / f'{name}.fchk'


def test_fixtures_match_their_generation_script(tmp_path):
    generate.main(tmp_path)
    for name in NAMES:
        assert (tmp_path / f'{name}.fchk').read_bytes() == fixture(name).read_bytes(), name


def test_documented_layout_tables():
    # QVF Cartesian d is xx, xy, xz, yy, yz, zz; Gaussian xx, yy, zz, xy, xz, yz.
    indices, factors = _layout(2, False)
    assert indices == [0, 3, 5, 1, 2, 4]
    assert factors == pytest.approx([1, 1, 1] + [1 / math.sqrt(3)] * 3)
    # Gaussian f: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz.
    indices, factors = _layout(3, False)
    assert indices == [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]
    assert factors == pytest.approx([1, 1, 1] + [1 / math.sqrt(5)] * 6 + [1 / math.sqrt(15)])
    # g and higher: the reverse of the QVF (libint) order, zzzz first.
    assert _layout(4, False)[0] == list(range(14, -1, -1))
    assert _layout(4, False)[1][0] == 1 and _layout(4, False)[1][-1] == 1
    # Pure: m = 0, +1, -1, +2, -2, ...; QVF m = -l .. +l.
    assert _layout(2, True) == ([2, 3, 1, 4, 0], [1.0] * 5)
    assert _layout(3, True)[0] == [3, 4, 2, 5, 1, 6, 0]
    # p is x, y, z in fchk; QVF pure p is y, z, x.
    assert _layout(1, True)[0] == [2, 0, 1] and _layout(1, False)[0] == [0, 1, 2]


def test_fortran_number_forms():
    assert _real('1.23456789-100') == pytest.approx(1.23456789e-100, rel=1e-12)
    assert _real('-4.5D+01') == -45.0 and _real('2.0E-03') == 0.002


def dfact(n):
    return math.prod(range(n, 0, -2)) if n > 0 else 1


GAUSSIAN_CARTESIAN = {
    1: ['x', 'y', 'z'],
    2: ['xx', 'yy', 'zz', 'xy', 'xz', 'yz'],
    3: ['xxx', 'yyy', 'zzz', 'xyy', 'xxy', 'xxz', 'xzz', 'yzz', 'yyz', 'xyz'],
    4: ['zzzz', 'yzzz', 'yyzz', 'yyyz', 'yyyy', 'xzzz', 'xyzz', 'xyyz', 'xyyy', 'xxzz',
        'xxyz', 'xxyy', 'xxxz', 'xxxy', 'xxxx'],
}


def solid_harmonic(np, l, m, x, y, z):
    """Real solid harmonic without Condon-Shortley phase (Appendix A.2)."""
    r = np.sqrt(x * x + y * y + z * z); t = z / r; phi = np.arctan2(y, x); k = abs(m)
    legendre = np.polynomial.legendre.Legendre.basis(l)
    p = (1 - t * t) ** (k / 2) * legendre.deriv(k)(t) if k else legendre(t)
    scale = math.sqrt(2 * math.factorial(l - k) / math.factorial(l + k)) if k else 1
    trig = 1 if m == 0 else (np.cos(k * phi) if m > 0 else np.sin(k * phi))
    return r ** l * scale * p * trig


def gaussian_orbitals(name, points):
    """Evaluate the fixture's orbitals with Gaussian's documented AOs.

    Cartesian AOs are unit-normalized per component in the order above, pure
    AOs are m = 0, +1, -1, ... real solid harmonics on axially normalized
    primitives. Uses the fixture definition, not the importer.
    """
    np = pytest.importorskip('numpy')
    shells, channels = generate.orbitals(name)
    atoms = generate.FIXTURES[name]['atoms']; columns = []
    for shell_type, atom, exponents, first, second in shells:
        delta = points - np.array(atoms[atom - 1][1]); x, y, z = delta.T
        r2 = (delta * delta).sum(axis=1)
        parts = [(0, first), (1, second)] if shell_type == -1 else [(abs(shell_type), first)]
        for l, coefficients in parts:
            if shell_type < -1:
                radial = sum(c * (2 * a / math.pi) ** .75 * (4 * a) ** (l / 2) / math.sqrt(dfact(2 * l - 1))
                             * np.exp(-a * r2) for a, c in zip(exponents, coefficients))
                for m in [0] + [s * k for k in range(1, l + 1) for s in (1, -1)]:
                    columns.append(radial * solid_harmonic(np, l, m, x, y, z))
                continue
            for label in GAUSSIAN_CARTESIAN.get(l, ['']):
                lx, ly, lz = (label.count(a) for a in 'xyz')
                radial = sum(c * (2 * a / math.pi) ** .75 * (4 * a) ** (l / 2)
                             / math.sqrt(dfact(2 * lx - 1) * dfact(2 * ly - 1) * dfact(2 * lz - 1))
                             * np.exp(-a * r2) for a, c in zip(exponents, coefficients))
                columns.append(radial * x ** lx * y ** ly * z ** lz)
    ao = np.array(columns).T
    return [ao @ np.array(rows).T for rows in channels.values()]


def random_points(np, centres, n=24, seed=7):
    rng = np.random.default_rng(seed)
    centres = np.asarray(centres, dtype=float)
    return centres[rng.integers(len(centres), size=n)] + rng.normal(scale=.9, size=(n, 3))


@pytest.mark.parametrize('name', NAMES)
def test_fixture_orbital_values_against_documented_aos(tmp_path, name):
    np = pytest.importorskip('numpy')
    path = tmp_path / 'import.qvf'; report = import_fchk(fixture(name), path); clean(path)
    spec = generate.FIXTURES[name]
    points = random_points(np, [xyz for _, xyz in spec['atoms']])
    for actual, expected in zip(qvf_orbitals(path, points), gaussian_orbitals(name, points)):
        np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-9)
    assert report['contraction_norms'] == pytest.approx([1] * len(report['contraction_norms']), abs=1e-7)
    assert report['omitted_text_fields'] == ['Route'] and 'Info1-9' in report['omitted_fields']
    with Reader(path) as reader:
        metadata = reader.json('wavefunction', 'mo_metadata')
        basis = reader.json('wavefunction', 'basis')
        provenance = reader.raw_manifest['provenance']
        assert provenance['charge'] == spec['charge'] and provenance['multiplicity'] == spec['multiplicity']
        assert basis['n_electrons'] == spec['n_alpha'] + spec['n_beta'] == provenance['n_electrons']
        energies = reader.raw_manifest.get('energies')
        if spec['energy'] is None:
            assert energies is None
        else:
            assert energies['total']['value'] == spec['energy']
            assert energies['method'] == spec['method'].lower()
        if spec['mulliken'] is None:
            assert 'atom_properties' not in [s.id for s in reader.sections]
        else:
            assert reader.tensor('atom_properties', 'mulliken_charge').tolist() == spec['mulliken']
        # Densities in the QVF basis agree with the derived occupations.
        if spec['unrestricted']:
            assert metadata['spin'] == 'unrestricted'
            pairs = [('density_matrix_alpha', 'alpha'), ('density_matrix_beta', 'beta')]
            assert metadata['alpha']['occupations'][:spec['n_alpha'] + 1] == [1] * spec['n_alpha'] + [0]
            assert sum(metadata['beta']['occupations']) == spec['n_beta']
        else:
            assert metadata['spin'] == ('restricted' if spec['n_alpha'] == spec['n_beta'] else 'restricted_open')
            pairs = [('density_matrix', None)]
            assert sum(metadata['occupations']) == spec['n_alpha'] + spec['n_beta']
        for role, channel in pairs:
            occ = np.array((metadata[channel] if channel else metadata)['occupations'])
            c = np.array(reader.tensor('wavefunction', 'mo_coefficients' + ('_' + channel if channel else '')).tolist())
            p = np.array(reader.tensor('wavefunction', role).tolist())
            np.testing.assert_allclose(p, (c.T * occ) @ c, atol=1e-7)


@pytest.mark.parametrize('name', NAMES)
def test_fixture_mo_overlap_matches_iodata(tmp_path, name):
    """IOData's own fchk conventions and overlap code give the same MO metric."""
    np = pytest.importorskip('numpy'); pytest.importorskip('iodata')
    from iodata import load_one
    from iodata.overlap import compute_overlap
    data = load_one(str(fixture(name)))
    s_file = compute_overlap(data.obasis, data.atcoords)
    path = tmp_path / 'import.qvf'; import_fchk(fixture(name), path)
    with Reader(path) as reader:
        s = np.array(reader.tensor('wavefunction', 'overlap').tolist())
        unrestricted = data.mo.kind == 'unrestricted'
        for channel, coeffs in [('alpha', data.mo.coeffsa)] + ([('beta', data.mo.coeffsb)] if unrestricted else []):
            c = np.array(reader.tensor('wavefunction', 'mo_coefficients' + ('_' + channel if unrestricted else '')).tolist())
            np.testing.assert_allclose(c @ s @ c.T, coeffs.T @ s_file @ coeffs, atol=5e-8)


@pytest.mark.parametrize('name', NAMES)
def test_fixture_round_trip(tmp_path, name):
    np = pytest.importorskip('numpy')
    first = tmp_path / 'first.qvf'; import_fchk(fixture(name), first)
    output = tmp_path / 'export.fchk'; report = export_fchk(first, output)
    spec = generate.FIXTURES[name]
    assert report['charge'] == spec['charge'] and report['multiplicity'] == spec['multiplicity']
    assert report['density_source'] == 'archive' and report['omitted_members'] == ['overlap']
    second = tmp_path / 'second.qvf'; again = import_fchk(output, second); clean(second)
    assert again['spin'] == ('unrestricted' if spec['unrestricted'] else
                             'restricted' if spec['n_alpha'] == spec['n_beta'] else 'restricted_open')
    points = random_points(np, [xyz for _, xyz in spec['atoms']], seed=3)
    for a, b in zip(qvf_orbitals(first, points), qvf_orbitals(second, points)):
        np.testing.assert_allclose(a, b, rtol=1e-7, atol=1e-9)
    text = output.read_text()
    # A pure d/f basis keeps its shell-type signs and flags.
    assert ('Pure/Cartesian d shells' in text) and ('Total SCF Density' in text)
    assert ('Spin SCF Density' in text) == (spec['n_alpha'] != spec['n_beta'])


# -- PySCF calculations through the exporter, read by cclib, IOData and qvf --

CASES = [
    # (label, method, cartesian, basis, spin 2S, independent cclib check)
    ('rhf-sph-f', 'RHF', False, 'cc-pvtz', 0, True),
    ('uhf-cart-d', 'UHF', True, '6-31g**', 1, True),
    ('rohf-sph-d', 'ROHF', False, 'cc-pvdz', 2, True),
    ('rhf-cart-g', 'RHF', True, {'F': 'cc-pvqz', 'H': 'sto-3g'}, 0, False),
    ('uhf-sph-g', 'UHF', False, {'F': 'cc-pvqz', 'H': 'sto-3g'}, 2, False),
]


@pytest.mark.parametrize('label, method, cart, basis, spin, with_cclib', CASES, ids=[c[0] for c in CASES])
def test_pyscf_export_import_independent(tmp_path, label, method, cart, basis, spin, with_cclib):
    np = pytest.importorskip('numpy'); pytest.importorskip('pyscf')
    from pyscf import gto, scf, lib
    from qvf.adapters.pyscf import write_scf
    lib.num_threads(1)
    atom = ('O 0 0 0.22; H 0 1.43 -0.89; H 0 -1.43 -0.89' if 'g' not in label[-1]
            else 'F 0 0 0; H 0.3 0.2 1.73')
    charge = 1 if label == 'uhf-cart-d' else 0
    mol = gto.M(atom=atom, unit='Bohr', basis=basis, cart=cart, spin=spin, charge=charge, verbose=0)
    mf = getattr(scf, method)(mol); mf.conv_tol = 1e-11; mf.run()
    assert mf.converged
    source = tmp_path / 'pyscf.qvf'; write_scf(mf, source)
    output = tmp_path / 'pyscf.fchk'; report = export_fchk(source, output)
    assert report['charge'] == mol.charge and report['multiplicity'] == mol.spin + 1
    assert report['mulliken_charges'] and report['total_energy'] and report['omitted_sections'] == ['job_spec']
    unrestricted = method == 'UHF'
    mo = list(mf.mo_coeff) if unrestricted else [mf.mo_coeff]
    energies = list(mf.mo_energy) if unrestricted else [mf.mo_energy]
    with Reader(source) as reader:
        mulliken = reader.tensor('atom_properties', 'mulliken_charge').tolist()
        reference_overlap = np.array(reader.tensor('wf', 'overlap').tolist())
    if with_cclib:
        cclib = pytest.importorskip('cclib')
        data = cclib.io.ccread(str(output))
        assert data.charge == mol.charge and data.mult == mol.spin + 1
        assert list(data.atomnos) == list(mol.atom_charges())
        np.testing.assert_allclose(data.atomcoords[0] * (1 / 0.529177210903), mol.atom_coords(), atol=1e-7)
        assert data.nbasis == mol.nao_nr() and data.nmo == mo[0].shape[1]
        assert list(data.homos) == ([mol.nelec[0] - 1, mol.nelec[1] - 1] if unrestricted
                                    else [mol.nelec[0] - 1])
        # cclib reports eV with its own constant; compare in hartree.
        from cclib.parser.utils import convertor
        for parsed, reference in zip(data.moenergies, energies):
            np.testing.assert_allclose(convertor(np.asarray(parsed), 'eV', 'hartree'), reference, rtol=2e-8, atol=1e-8)
        np.testing.assert_allclose(data.atomcharges['mulliken'], mulliken, atol=1e-7)
        assert convertor(data.scfenergies[0], 'eV', 'hartree') == pytest.approx(mf.e_tot, abs=1e-9)
        assert len(data.mocoeffs) == len(mo) and data.mocoeffs[0].shape == (mo[0].shape[1], mol.nao_nr())
    try:
        from iodata import load_one
        from iodata.overlap import compute_overlap
    except ImportError:
        load_one = None
    if load_one is not None:
        # IOData reads the exported file with its own conventions: the SCF
        # orbitals must be orthonormal in its overlap of that basis.
        data = load_one(str(output))
        s_file = compute_overlap(data.obasis, data.atcoords)
        for coeffs in ([data.mo.coeffsa, data.mo.coeffsb] if unrestricted else [data.mo.coeffsa]):
            np.testing.assert_allclose(coeffs.T @ s_file @ coeffs, np.eye(coeffs.shape[1]), atol=2e-7)
        assert data.energy == pytest.approx(mf.e_tot, abs=1e-10)
    back = tmp_path / 'back.qvf'; again = import_fchk(output, back); clean(back)
    assert again['spin'] == {'RHF': 'restricted', 'UHF': 'unrestricted', 'ROHF': 'restricted_open'}[method]
    points = random_points(np, mol.atom_coords(), seed=11)
    ao = mol.eval_gto('GTOval_cart' if cart else 'GTOval_sph', points)
    for actual, coeff in zip(qvf_orbitals(back, points), mo):  # PySCF's own AOs
        np.testing.assert_allclose(actual, ao @ coeff, rtol=1e-6, atol=1e-7)
    with Reader(back) as reader:
        s = np.array(reader.tensor('wavefunction', 'overlap').tolist())
        np.testing.assert_allclose(s, reference_overlap, atol=1e-8)
        roles = ['mo_coefficients_alpha', 'mo_coefficients_beta'] if unrestricted else ['mo_coefficients']
        for role in roles:
            c = np.array(reader.tensor('wavefunction', role).tolist())
            np.testing.assert_allclose(c @ s @ c.T, np.eye(len(c)), atol=2e-7)
        densities = ['density_matrix_alpha', 'density_matrix_beta'] if unrestricted else ['density_matrix']
        count = sum(np.trace(np.array(reader.tensor('wavefunction', d).tolist()) @ s) for d in densities)
        assert count == pytest.approx(mol.nelectron, abs=1e-6)
        assert reader.raw_manifest['energies']['total']['value'] == pytest.approx(mf.e_tot, abs=1e-12)
        np.testing.assert_allclose(reader.tensor('atom_properties', 'mulliken_charge').tolist(), mulliken, atol=1e-7)


# -- refusals -----------------------------------------------------------------

def minimal(**changes):
    """A one-atom, one-s-shell restricted fchk built with the fixture writer."""
    fields = dict(charge=0, multiplicity=1, electrons=2, alpha=1, beta=1, nbasis=1, z=2,
                  exponent=1.2, coefficient=1.0, coefficients=[1.0], shell_type=0)
    fields.update(changes)
    lines = ['minimal', 'SP        RHF                           Gen']
    for key, value in (('Number of atoms', 1), ('Charge', fields['charge']),
                       ('Multiplicity', fields['multiplicity']), ('Number of electrons', fields['electrons']),
                       ('Number of alpha electrons', fields['alpha']), ('Number of beta electrons', fields['beta']),
                       ('Number of basis functions', fields['nbasis'])):
        lines.append(generate.header(key, 'I', value=value))
    lines += generate.integers('Atomic numbers', [fields['z']])
    lines += generate.reals('Current cartesian coordinates', [0.0, 0.0, 0.0])
    lines += generate.integers('Shell types', [fields['shell_type']]) + generate.integers('Number of primitives per shell', [1])
    lines += generate.integers('Shell to atom map', [1])
    lines += generate.reals('Primitive exponents', [fields['exponent']])
    lines += generate.reals('Contraction coefficients', [fields['coefficient']])
    lines += generate.reals('Alpha Orbital Energies', [-0.5])
    lines += generate.reals('Alpha MO coefficients', fields['coefficients'])
    return '\n'.join(lines) + '\n'


def write(tmp_path, text):
    path = tmp_path / 'input.fchk'; path.write_text(text); return path


def test_minimal_imports(tmp_path):
    path = tmp_path / 'min.qvf'; report = import_fchk(write(tmp_path, minimal()), path); clean(path)
    assert report['spin'] == 'restricted' and report['density_members'] == []
    with Reader(path) as reader:
        assert reader.json('wavefunction', 'mo_metadata')['occupations'] == [2]
        assert reader.tensor('wavefunction', 'overlap').tolist() == [[pytest.approx(1)]]


MIN = minimal()
DENSITY = generate.reals('Total SCF Density', [1.0])
D_SHELL = minimal(nbasis=5, coefficients=[0.1] * 5, shell_type=-2)


@pytest.mark.parametrize('text, message', [
    ('', 'title line'),
    (MIN.replace('Atomic numbers', 'Atomic Numbers'), "'Atomic numbers'"),
    (MIN.replace('           2\nCurrent', '           0\nCurrent'), 'ghost'),
    (MIN + '\n'.join(generate.reals('Nuclear charges', [1.0])) + '\n', 'effective core'),
    (MIN + '\n'.join(generate.integers('ECP-LenNCF', [0])) + '\n', 'effective core'),
    (minimal(exponent=-1.2), 'positive'),
    (minimal(coefficient=0.5), 'squared norm'),
    (minimal(nbasis=2), 'Number of basis functions'),
    (minimal(multiplicity=3), 'disagree'),
    (minimal(charge=1), 'disagree'),
    (minimal(alpha=0, beta=2), 'disagree'),
    (MIN.rsplit('\n', 2)[0] + '\n', 'truncated'),
    (MIN + generate.header('Charge', 'I', value=0) + '\n', 'repeated'),
    (MIN.replace('Charge  ', 'charge = 0\nCharge  '), 'field header'),
    (MIN.replace('Shell to atom map                          I   N=           1\n           1',
                 'Shell to atom map                          I   N=           1\n           2'), 'no atom'),
    (MIN.replace('Shell types                                I   N=           1\n           0',
                 'Shell types                                I   N=           1\n           9'), 'l <= 7'),
    (MIN.replace('Shell types                                I   N=           1\n           0',
                 'Shell types                                I   N=           1\n          -1'), 'P\\(S=P\\)'),
    (minimal(coefficients=[1.0, 0.0]), 'expected 1'),
    (MIN + '\n'.join(generate.reals('Coordinates of each shell', [0.0, 0.0, 0.5])) + '\n', 'not centred'),
    (MIN + '\n'.join(generate.reals('Beta Orbital Energies', [-0.5])) + '\n', 'Beta Orbital Energies'),
    (MIN + '\n'.join(generate.reals('Total SCF Density', [1.0])) + '\n', 'do not reproduce'),
    (MIN + '\n'.join(generate.reals('Total SCF Density', [1.0, 2.0])) + '\n', 'packed lower triangle'),
    (MIN + '\n'.join(generate.reals('Mulliken Charges', [0.0, 0.0])) + '\n', 'one value per atom'),
    (D_SHELL + generate.header('Pure/Cartesian d shells', 'I', value=1) + '\n', 'disagree with'),
    (MIN.replace('1.20000000E+00', 'nan'), 'finite'),
])
def test_import_refusal_is_atomic(tmp_path, text, message):
    source = write(tmp_path, text)
    with pytest.raises(ValueError, match=message):
        import_fchk(source, tmp_path / 'invalid.qvf')
    assert [p.name for p in tmp_path.iterdir()] == ['input.fchk']


def archive(tmp_path, *, metadata, shells=None, rows=None, geometry=None, basis_extra=None, name='source.qvf'):
    """A synthetic QVF wavefunction written without the fchk importer."""
    shells = shells or [dict(center=0, l=0, pure=True, exponents=[1.2], coefficients=[1.0])]
    nao = sum(2 * s['l'] + 1 if s.get('pure', True) else (s['l'] + 1) * (s['l'] + 2) // 2 for s in shells)
    path = tmp_path / name
    with Writer(path, 'synthetic', '1') as writer:
        writer.add_section({'id': 'structure', 'kind': 'structure'})
        writer.add_json('structure', 'structure', geometry or {
            'atoms': [{'symbol': 'He', 'atomic_number': 2, 'position': [0, 0, 0]}], 'pbc': [False] * 3})
        writer.add_section({'id': 'wavefunction', 'kind': 'wavefunction.gto', 'structure_ref': 'structure'})
        writer.add_json('wavefunction', 'basis', dict({'structure_ref': 'structure', 'shells': shells, 'n_ao': nao},
                                            **(basis_extra or {})))
        writer.add_json('wavefunction', 'mo_metadata', metadata)
        channels = ['alpha', 'beta'] if metadata['spin'] == 'unrestricted' else [None]
        for channel in channels:
            payload = metadata[channel] if channel else metadata
            n = len(payload['energies'])
            data = rows or [[1.0 if i == j else 0.0 for j in range(nao)] for i in range(n)]
            _binary(writer, 'wavefunction', 'mo_coefficients' + ('_' + channel if channel else ''),
                    (v for row in data[:n] for v in row), [n, nao])
    return path


def test_export_derives_charge_multiplicity_and_density(tmp_path):
    np = pytest.importorskip('numpy')
    shells = [dict(center=0, l=0, pure=True, exponents=[1.2], coefficients=[0.5]),
              dict(center=0, l=1, pure=True, exponents=[0.9], coefficients=[2.0])]
    metadata = dict(spin='restricted_open', energies=[-1, -.5, -.4, .1], occupations=[2, 1, 0, 0])
    path = archive(tmp_path, metadata=metadata, shells=shells); clean(path)
    output = tmp_path / 'out.fchk'; report = export_fchk(path, output)
    assert (report['charge'], report['multiplicity']) == (-1, 2)
    assert report['density_source'] == 'orbitals' and report['contraction_norms'] == [.5, 2.0]
    back = tmp_path / 'back.qvf'; again = import_fchk(output, back); clean(back)
    assert again['spin'] == 'restricted_open' and 'Spin SCF Density' in again['omitted_fields']
    with Reader(back) as reader:
        c = np.array(reader.tensor('wavefunction', 'mo_coefficients').tolist())
        # Contractions are renormalized (norms 0.5 and 2), and without a pure
        # d shell in the file the p shell returns Cartesian: y, z, x -> x, y, z.
        np.testing.assert_allclose(c, [[.5, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2], [0, 2, 0, 0]], atol=1e-8)
        assert reader.json('wavefunction', 'basis')['shells'][1]['pure'] is False
    points = random_points(np, [[0, 0, 0]])
    np.testing.assert_allclose(qvf_orbitals(back, points)[0], qvf_orbitals(path, points)[0], rtol=1e-8, atol=1e-12)


@pytest.mark.parametrize('kwargs, message', [
    (dict(metadata=dict(spin='restricted', energies=[-.5], occupations=[1.5])), 'aufbau'),
    (dict(metadata=dict(spin='restricted', energies=[-.5, .1], occupations=[0, 2]),
          shells=[dict(center=0, l=0, exponents=[1.2], coefficients=[1.0]),
                  dict(center=0, l=0, exponents=[.3], coefficients=[1.0])]), 'aufbau'),
    (dict(metadata=dict(spin='restricted', energies=[-.5], occupations=[1])), 'restricted_open'),
    (dict(metadata=dict(spin='unrestricted', alpha=dict(energies=[-.5, .1], occupations=[1, 0]),
                        beta=dict(energies=[-.4], occupations=[1])),
          shells=[dict(center=0, l=0, exponents=[1.2], coefficients=[1.0]),
                  dict(center=0, l=0, exponents=[.3], coefficients=[1.0])]), 'same number'),
    (dict(metadata=dict(spin='restricted', energies=[-.5], occupations=[2]),
          geometry={'atoms': [{'symbol': 'He', 'atomic_number': 2, 'position': [0, 0, 0]}],
                    'pbc': [True] * 3, 'lattice_vectors': [[3, 0, 0], [0, 3, 0], [0, 0, 3]]}), 'molecular'),
    (dict(metadata=dict(spin='restricted', energies=[-.5, -.4, .1, .2, .3, .4, .5, .6, .7, .8, .9],
                        occupations=[2] + [0] * 10),
          shells=[dict(center=0, l=2, pure=True, exponents=[1.0], coefficients=[1.0]),
                  dict(center=0, l=2, pure=False, exponents=[.5], coefficients=[1.0])]), 'mixes'),
    (dict(metadata=dict(spin='restricted', energies=[-.5], occupations=[2]),
          basis_extra={'ecp': [{'center': 0, 'n_core_electrons': 0}]}), 'ECP'),
])
def test_export_refusals(tmp_path, kwargs, message):
    path = archive(tmp_path, **kwargs)
    with pytest.raises(ValueError, match=message):
        export_fchk(path, tmp_path / 'invalid.fchk')
    assert not (tmp_path / 'invalid.fchk').exists()


def test_cli_and_standard_library_only(tmp_path, capsys):
    from qvf.cli import main
    path = tmp_path / 'cli.qvf'
    assert main(['convert', str(fixture('rhf_pure')), '--from', 'fchk', '--to', 'qvf', '-o', str(path)]) == 0
    assert json.loads(capsys.readouterr().out)['n_ao'] == 30
    with pytest.raises(FileExistsError):
        import_fchk(fixture('rhf_pure'), path)
    output = tmp_path / 'cli.fchk'
    assert main(['convert', str(path), '--from', 'qvf', '--to', 'fchk', '-o', str(output)]) == 0
    assert json.loads(capsys.readouterr().out)['spin'] == 'restricted'
    assert main(['convert', str(path), '--from', 'qvf', '--to', 'fchk', '-o', str(output)]) == 2
    result = subprocess.run([sys.executable, '-S', '-m', 'qvf', 'convert', str(fixture('uhf_cart')),
                             '--from', 'fchk', '--to', 'qvf', '-o', str(tmp_path / 'stdlib.qvf')],
                            env=dict(os.environ, PYTHONPATH=str(HERE.parent)), text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    clean(tmp_path / 'stdlib.qvf')
