"""Molden producer dialects checked against PySCF AO values.

Each file is written from PySCF orbitals: PySCF's own Molden writer, or a
hand-built file that applies a dialect's documented transformation to the
same orbitals. The transformations here are written from the dialect table in
docs/converters.md (IOData 1.0.1 documents the ORCA, Turbomole, CFOUR and
Psi4 deviations), with their own textbook normalization, not with the
importer's helpers. Orbital values after import are compared with PySCF's
AO evaluation at off-axis points on two centers.
"""
import json
import math
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qvf import Reader
from qvf.convert import import_molden
from qvf._legacy_writer import BOHR_PER_ANGSTROM
from test_convert_molden import SIMPLE, clean, qvf_orbitals, source

np = pytest.importorskip('numpy')

# Molden's Cartesian component orders, from the Molden format description.
CARTESIAN = {2: ['xx', 'yy', 'zz', 'xy', 'xz', 'yz'],
             3: ['xxx', 'yyy', 'zzz', 'xyy', 'xxy', 'xxz', 'xzz', 'yzz', 'yyz', 'xyz'],
             4: ['xxxx', 'yyyy', 'zzzz', 'xxxy', 'xxxz', 'yyyx', 'yyyz', 'zzzx', 'zzzy',
                 'xxyy', 'xxzz', 'yyzz', 'xxyz', 'yyxz', 'zzxy']}
POINTS = [[.2, .3, .4], [-.7, .2, .5], [1.1, -.3, -.8], [-.6, -.5, .9]]
ORCA_TITLE = 'Molden file created by orca_2mkl for BaseName=synthetic'


def dfact(n):
    return math.prod(range(n, 0, -2)) if n > 0 else 1


def cart_norm(alpha, powers):
    """Unit norm of x^a y^b z^c exp(-alpha r^2)."""
    l = sum(powers)
    return ((2 * alpha / math.pi) ** .75 * (4 * alpha) ** (l / 2)
            / math.sqrt(math.prod(dfact(2 * n - 1) for n in powers)))


def molecule(cart):
    pytest.importorskip('pyscf')
    from pyscf import gto, lib
    lib.num_threads(1)
    return gto.M(atom='He 0 .2 -.1; He .7 -.3 1.1', unit='Bohr', cart=cart, verbose=0,
                 basis={'He': [[0, [2.3, .4], [.7, .6]], [1, [1.9, .6], [.5, .4]],
                               [2, [1.3, .7], [.3, .4]], [3, [1.1, .8], [.4, .3]], [4, [.9, 1.]]]})


def molden_ao_values(pure):
    """PySCF AO values in Molden order and normalization, shell by shell."""
    pytest.importorskip('pyscf')
    from pyscf.tools import molden
    blocks = {}
    for cart in (False, True):
        mol = molecule(cart)
        ao = mol.eval_gto('GTOval_cart' if cart else 'GTOval_sph', np.array(POINTS))
        if cart:  # Molden Cartesian AOs are unit-normalized (PySCF's writer agrees).
            ao = ao / np.sqrt(mol.intor('int1e_ovlp').diagonal())
        ao = ao[:, molden.order_ao_index(mol)]
        offset = 0; blocks[cart] = []
        for ib in range(mol.nbas):
            l = mol.bas_angular(ib); size = (l + 1) * (l + 2) // 2 if cart else 2 * l + 1
            blocks[cart].append(ao[:, offset:offset + size]); offset += size
    mol = molecule(False)
    shells = [(mol.bas_atom(ib), mol.bas_angular(ib), mol.bas_exp(ib), mol.bas_ctr_coeff(ib)[:, 0])
              for ib in range(mol.nbas)]
    values = [blocks[not pure.get(l, l < 2)][i] for i, (_, l, _, _) in enumerate(shells)]
    return mol, shells, np.hstack(values)


def primitive_written(dialect, l, pure, alpha, c):
    """The contraction coefficient a producer writes for a normalized one."""
    if dialect == 'orca':  # norm of the (0,0,0), (1,0,0), (1,1,0), (1,1,1), (2,1,1) primitive
        return c * cart_norm(alpha, [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)][l])
    if dialect == 'xtb':
        return c * cart_norm(alpha, (l, 0, 0))
    if dialect == 'turbomole' and not pure and l >= 2:
        return c / math.sqrt(dfact(2 * l - 1))
    if dialect == 'renormalize':
        return c * 1.37
    return c


def component_written(dialect, l, pure):
    """Written MO coefficient / Molden-standard coefficient per component."""
    if pure or l < 2:
        size = 2 * l + 1 if l >= 2 else (1 if l == 0 else 3)
        if dialect == 'orca' and l >= 3:  # m = +-3 (f, g) and +-4 (g) negated
            return [1.] * 5 + [-1.] * (size - 5)
        return [(-1.) ** l if dialect == 'xtb' else 1.] * size
    counts = [[name.count(a) for a in 'xyz'] for name in CARTESIAN[l]]
    if dialect == 'cfour':
        return [1 / math.sqrt(math.prod(dfact(2 * n - 1) for n in p)) for p in counts]
    if dialect == 'psi4-1.3':
        return [math.sqrt(dfact(2 * l - 1) / math.prod(dfact(2 * n - 1) for n in p)) for p in counts]
    return [(-1.) ** l if dialect == 'xtb' else 1.] * len(counts)


def write(path, *, dialect='standard', pure=None, flags=None, unit='AU', spin='restricted',
          title=None, labels=False, fortran=False, compact=False, header='[Molden Format]',
          order=('Sym', 'Ene', 'Spin', 'Occup'), syms=None, occupations=None):
    """Write one dialect file; return the PySCF reference orbital values."""
    pure = pure or {}
    mol, shells, ao = molden_ao_values(pure)
    nao = ao.shape[1]; rng = np.random.default_rng(7)
    sets = [('Alpha', rng.normal(size=(4, nao)))]
    if spin == 'unrestricted':
        sets.append(('Beta', rng.normal(size=(3, nao))))
    num = (lambda v: format(v, '.17E').replace('E', 'D')) if fortran else (lambda v: format(v, '.17g'))
    lines = [header] + (['[Title]', title] if title else [])
    lines.append('[Atoms] ' + unit)
    for ia in range(mol.natm):
        xyz = mol.atom_coord(ia) / (1 if unit in ('AU', '(AU)') else BOHR_PER_ANGSTROM)
        label = f'HE{ia + 1}' if labels else 'He'
        lines.append(f'{label} {ia + 1} 2 ' + ' '.join(num(v) for v in xyz))
    lines.append('[GTO]'); center = None; signs = []
    for atom, l, exponents, coefficients in shells:
        if atom != center:
            if center is not None and not compact:
                lines.append('')
            center = atom; lines.append(f'{atom + 1}' if compact else f'{atom + 1} 0')
        lines.append(f'{"spdfg"[l]} {len(exponents)} 1.00')
        is_pure = pure.get(l, l < 2)
        lines.extend(f'{num(a)} {num(primitive_written(dialect, l, is_pure, a, c))}'
                     for a, c in zip(exponents, coefficients))
        signs.extend(component_written(dialect, l, is_pure))
    lines.append('')
    lines.extend(flags or [])
    lines.append('[MO]')
    signs = np.array(signs)
    for spin_label, matrix in sets:
        for i, row in enumerate(matrix):
            occupation = (occupations or [2, 1.5, .25, 0])[i] / (2 if spin == 'unrestricted' else 1)
            values = {'Sym': (syms or ['1a1', "2a'", '1.1', 'A'])[i], 'Ene': num(-.5 + .2 * i),
                      'Spin': spin_label, 'Occup': num(occupation)}
            lines.extend(f' {key}= {values[key]}' for key in order)
            lines.extend(f' {k} {num(v)}' for k, v in enumerate(row * signs, 1))
    path.write_text('\n'.join(lines) + '\n')
    return [ao @ matrix.T for _, matrix in sets]


def check(tmp_path, expected, *, spin='restricted', dialect='standard', file='dialect.molden'):
    path = tmp_path / 'dialect.qvf'
    report = import_molden(tmp_path / file, path, spin=spin, dialect=dialect)
    clean(path)
    assert report['dialect'] == dialect and report['findings'] == []
    for actual, reference in zip(qvf_orbitals(path, np.array(POINTS)), expected):
        np.testing.assert_allclose(actual, reference, rtol=3e-12, atol=3e-13)
    return path, report


ALL_PURE = {2: True, 3: True, 4: True}
ALL_CART = {2: False, 3: False, 4: False}
FLAG_CASES = {
    'absent': (ALL_CART, []),
    '6D10F15G': (ALL_CART, ['[6D]', '[10F]', '[15G]']),
    '5D': ({2: True, 3: True, 4: False}, ['[5D]']),
    '5D7F': ({2: True, 3: True, 4: False}, ['[5D7F]']),
    '5D7F 9G': (ALL_PURE, ['[5D7F]', '[9G]']),
    '5d 7f 9g lines': (ALL_PURE, ['[5d]', '[7f]', '[9g]']),
    '5D10F': ({2: True, 3: False, 4: False}, ['[5D10F]']),
    '5D10F 9G': ({2: True, 3: False, 4: True}, ['[5D10F]', '[9G]']),
    '7F': ({2: False, 3: True, 4: False}, ['[7F]']),
    '7F 9G': ({2: False, 3: True, 4: True}, ['[7F]', '[9G]']),
    '9G': ({2: False, 3: False, 4: True}, ['[9G]']),
}


@pytest.mark.parametrize('case', FLAG_CASES)
def test_angular_flags_including_mixed_conventions(tmp_path, case):
    pure, flags = FLAG_CASES[case]
    expected = write(tmp_path / 'dialect.molden', pure=pure, flags=flags)
    path, _ = check(tmp_path, expected)
    with Reader(path) as reader:
        shells = reader.json('wavefunction', 'basis')['shells']
    assert {s['l']: s['pure'] for s in shells if s['l'] >= 2} == pure


@pytest.mark.parametrize('variant', [
    dict(unit='Angs'), dict(unit='(AU)'), dict(unit='(Angs)'), dict(header='[MOLDEN FORMAT]'),
    dict(header='[Molden format]'), dict(fortran=True), dict(compact=True), dict(labels=True),
    dict(order=('Ene', 'Spin', 'Occup', 'Sym')), dict(order=('Spin', 'Sym', 'Occup', 'Ene')),
    dict(syms=['1', 'b1u', 'Ag', '12.3']), dict(occupations=[1.9987, 1.0, .0013, 0]),
    dict(title='synthetic Molpro/Psi4/Multiwfn style file'),
])
def test_text_variants_of_common_producers(tmp_path, variant):
    expected = write(tmp_path / 'dialect.molden', pure=ALL_PURE, flags=['[5D]', '[9G]'], **variant)
    path, _ = check(tmp_path, expected)
    with Reader(path) as reader:
        metadata = reader.json('wavefunction', 'mo_metadata')
    assert metadata['symmetry_labels'] == variant.get('syms', ['1a1', "2a'", '1.1', 'A'])
    assert metadata['occupations'] == pytest.approx(variant.get('occupations', [2, 1.5, .25, 0]))


@pytest.mark.parametrize('cart', [False, True])
def test_pyscf_writer_output_is_standard(tmp_path, cart):
    pytest.importorskip('pyscf')
    from pyscf.tools import molden
    mol = molecule(cart); nao = mol.nao_nr()
    coeff = np.random.default_rng(3).normal(size=(nao, 3))
    with (tmp_path / 'dialect.molden').open('w') as stream:
        molden.header(mol, stream)
        molden.orbital_coeff(mol, stream, coeff, ene=[-1, 0, 1], occ=[2, 2, 0])
    ao = mol.eval_gto('GTOval_cart' if cart else 'GTOval_sph', np.array(POINTS))
    check(tmp_path, [ao @ coeff])


@pytest.mark.parametrize('dialect, pure, flags', [
    ('orca', ALL_PURE, ['[5D]', '[7F]', '[9G]']),
    ('turbomole', ALL_CART, []),
    ('turbomole', {2: True, 3: False, 4: False}, ['[5D10F]']),
    ('cfour', ALL_CART, []),
    ('cfour', {2: False, 3: True, 4: False}, ['[7F]']),
    ('psi4-1.3', ALL_CART, []),
    ('xtb', ALL_CART, []),
    ('xtb', ALL_PURE, ['[5D7F]', '[9G]']),
    ('renormalize', ALL_PURE, ['[5D]', '[9G]']),
    ('renormalize', ALL_CART, []),
])
@pytest.mark.parametrize('spin', ['restricted', 'unrestricted'])
def test_producer_dialects(tmp_path, dialect, pure, flags, spin):
    title = ORCA_TITLE if dialect == 'orca' else None
    expected = write(tmp_path / 'dialect.molden', dialect=dialect, pure=pure, flags=flags,
                     spin=spin, title=title)
    _, report = check(tmp_path, expected, spin=spin, dialect=dialect)
    if dialect == 'renormalize':
        assert report['contraction_norms'] == pytest.approx([1.37] * len(report['contraction_norms']))
    else:
        assert report['contraction_norms'] == pytest.approx([1] * len(report['contraction_norms']))


@pytest.mark.parametrize('dialect, pure, flags', [
    ('orca', ALL_PURE, ['[5D]', '[7F]', '[9G]']), ('turbomole', ALL_CART, []),
    ('cfour', ALL_CART, []), ('psi4-1.3', ALL_CART, []),
])
def test_standard_reading_of_a_dialect_file_is_wrong_or_refused(tmp_path, dialect, pure, flags):
    """The dialects matter: read as standard, each file is refused or differs."""
    expected = write(tmp_path / 'dialect.molden', dialect=dialect, pure=pure, flags=flags,
                     title=ORCA_TITLE if dialect == 'orca' else None)
    if dialect in ('orca', 'turbomole'):
        with pytest.raises(ValueError, match="squared norm.*" + ("dialect 'orca'" if dialect == 'orca' else '')):
            import_molden(tmp_path / 'dialect.molden', tmp_path / 'wrong.qvf', spin='restricted')
        assert not (tmp_path / 'wrong.qvf').exists()
    else:
        import_molden(tmp_path / 'dialect.molden', tmp_path / 'wrong.qvf', spin='restricted')
        actual = qvf_orbitals(tmp_path / 'wrong.qvf', np.array(POINTS))[0]
        assert not np.allclose(actual, expected[0], rtol=1e-6, atol=1e-8)


def test_orca_sign_convention_is_needed_beyond_normalization(tmp_path):
    """Renormalizing an ORCA file fixes the norms but not the f/g phases."""
    expected = write(tmp_path / 'dialect.molden', dialect='orca', pure=ALL_PURE,
                     flags=['[5D]', '[7F]', '[9G]'], title=ORCA_TITLE)
    import_molden(tmp_path / 'dialect.molden', tmp_path / 'r.qvf', spin='restricted', dialect='renormalize')
    actual = qvf_orbitals(tmp_path / 'r.qvf', np.array(POINTS))[0]
    assert not np.allclose(actual, expected[0], rtol=1e-6, atol=1e-8)


GTO_D = SIMPLE.replace('s 1 1.00', 'd 1 1.00').replace('1 1.0\n', ''.join(f'{i} .1\n' for i in range(1, 7)))


@pytest.mark.parametrize('text, dialect, message', [
    (SIMPLE.replace('Occup= 2', 'Occup= -0.00001'), 'standard', 'negative'),
    (SIMPLE.replace('Occup= 2', 'Occup= 2.0001'), 'standard', 'outside 0 to 2'),
    (SIMPLE.replace('He 1 2', 'He 1 1'), 'standard', 'ECP effective core charge'),
    (SIMPLE.replace('He 1 2', 'C1 1 2'), 'standard', 'nuclear charge'),
    (SIMPLE.replace('He 1 2', 'Bq 1 2'), 'standard', 'element symbol'),
    (SIMPLE + '[core]\n1 : 2\n', 'standard', 'ECP'),
    (SIMPLE + '[Pseudo]\nHe 1 0\n', 'standard', 'ECP'),
    (SIMPLE + '[STO]\n', 'standard', 'unsupported Molden section'),
    (SIMPLE.replace('s 1 1.00', 'h 1 1.00'), 'standard', 's, p, sp, d, f and g'),
    (SIMPLE.replace('Spin= Alpha\n', ''), 'standard', 'Sym/Ene/Spin/Occup'),
    (SIMPLE.replace('1 1.0\n', '1 1.0 0.0\n'), 'standard', 'one real value'),
    (SIMPLE.replace('1.2 1.0', '1.2 0'), 'renormalize', 'nonzero contraction norm'),
    (SIMPLE.replace('1.2 1.0', '1.2 .4'), 'standard', "'renormalize'"),
    (SIMPLE, 'gaussian', 'unknown Molden dialect'),
    (GTO_D, 'orca', 'Cartesian shell'),
    (GTO_D.replace('3 .1\n', ''), 'standard', 'each AO coefficient'),
])
def test_unrepresentable_input_is_refused_with_a_reason(tmp_path, text, dialect, message):
    with pytest.raises(ValueError, match=message):
        import_molden(source(tmp_path, text), tmp_path / 'invalid.qvf', spin='restricted', dialect=dialect)
    assert [p.name for p in tmp_path.iterdir()] == ['input.molden']


def test_cli_dialect_option(tmp_path, capsys):
    from qvf.cli import main
    write(tmp_path / 'dialect.molden', dialect='orca', pure=ALL_PURE, flags=['[5D]', '[7F]', '[9G]'],
          title=ORCA_TITLE)
    base = ['convert', str(tmp_path / 'dialect.molden'), '--from', 'molden', '--to', 'qvf', '--spin', 'restricted']
    assert main(base + ['-o', str(tmp_path / 'a.qvf')]) == 2 and not (tmp_path / 'a.qvf').exists()
    capsys.readouterr()
    assert main(base + ['--molden-dialect', 'orca', '-o', str(tmp_path / 'b.qvf')]) == 0
    assert json.loads(capsys.readouterr().out)['dialect'] == 'orca'
    clean(tmp_path / 'b.qvf')
    xyz = tmp_path / 'x.xyz'; xyz.write_text('1\n\nHe 0 0 0\n')
    assert main(['convert', str(xyz), '--from', 'xyz', '--to', 'qvf', '--molden-dialect', 'orca',
                 '-o', str(tmp_path / 'c.qvf')]) == 2


def test_integral_float_nuclear_charge_and_lowercase_symbol(tmp_path):
    path = tmp_path / 'charge.qvf'
    import_molden(source(tmp_path, SIMPLE.replace('He 1 2', 'he 1 2.0')), path, spin='restricted')
    clean(path)
    with Reader(path) as reader:
        assert reader.json('structure', 'structure')['atoms'][0]['symbol'] == 'He'
