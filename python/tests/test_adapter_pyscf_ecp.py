"""ECP projectors from the PySCF adapter (WP-62), checked against PySCF itself.

* The radial functions in the archive equal those evaluated from PySCF's
  basis-set library (``gto.basis.load_ecp``), channel by channel.
* Radial quadrature of the archive's potential reproduces PySCF's
  one-centre ``ECPscalar`` integrals for s, p, d and f primitives, which
  fixes the radial-power convention (``r^(n-2)``) and the meaning of the
  semi-local channels (they act in addition to the local part).
* A molecule rebuilt from the archive's basis and ECP alone reproduces the
  SCF energy, for def2 (Stuttgart-Köln), LANL2DZ and CRENBL potentials.
"""
from __future__ import annotations

import pytest

pytest.importorskip("pyscf")
np = pytest.importorskip("numpy")

from pyscf import gto, scf  # noqa: E402

from qvf import Reader, validate_qvf  # noqa: E402
from qvf.adapters.pyscf import add_scf, write_scf  # noqa: E402
from qvf.units import BOHR_PER_ANGSTROM  # noqa: E402

HI = "H 0 0 0; I 0 0 1.609"


def _mol(atom, basis, ecp, **kw):
    return gto.M(atom=atom, basis=basis, ecp=ecp, unit="Angstrom", verbose=0, **kw)


def _write(tmp_path, mf, name, **options):
    path = tmp_path / f"{name}.qvf"
    write_scf(mf, str(path), **options)
    assert validate_qvf(str(path))["findings"] == []
    with Reader(str(path)) as reader:
        wf = next(s for s in reader.sections if s.kind == "wavefunction.gto")
        basis = reader.json(wf, "basis")
        atoms = reader.json("structure", "structure")["atoms"]
    return path, basis, atoms


def _u_archive(projector, r):
    """U(r) = Σ c r^(n-2) exp(-α r²) of spec § 5.3."""
    return sum(c * r ** (n - 2) * np.exp(-a * r * r) for a, c, n in
               zip(projector["exponents"], projector["coefficients"], projector["powers"]))


def _u_library(channel, r):
    """A channel of PySCF's library form: one ``[[α, c], ...]`` list per r order."""
    total = np.zeros_like(r)
    for rorder, terms in enumerate(channel):
        for term in terms:
            total = total + term[1] * r ** (rorder - 2) * np.exp(-term[0] * r * r)
    return total


def _channels(entry):
    local = [p for p in entry["projectors"] if p.get("local")]
    semi = {p["l"]: p for p in entry["projectors"] if not p.get("local")}
    return (local[0] if local else None), semi


@pytest.mark.parametrize("ecp,symbol", [("def2-svp", "I"), ("lanl2dz", "I"),
                                        ("crenbl", "I"), ("stuttgart_rsc", "Ag"),
                                        ("stuttgart", "I")])
def test_radial_functions_match_library(ecp, symbol):
    from qvf.adapters.pyscf import _ecp_entries
    mol = _mol(f"{symbol} 0 0 0", "def2-svp" if ecp == "stuttgart_rsc" else ecp, ecp,
               charge=-1 if symbol == "I" else 0, spin=0 if symbol == "I" else 1)
    (entry,) = _ecp_entries(mol)
    n_core, library = gto.basis.load_ecp(ecp, symbol)
    assert entry["n_core_electrons"] == n_core
    local, semi = _channels(entry)
    r = np.linspace(0.05, 6.0, 400)
    by_l = {l: channel for l, channel in library}
    scalar = {l: [[t[:2] for t in terms] for terms in channel] for l, channel in by_l.items()}
    assert set(semi) == {l for l in scalar if l >= 0}
    if not any(scalar.get(-1, [])):  # no local part: U_L = 0, nothing written
        assert local is None
    else:
        assert local["l"] == max(semi) + 1
        np.testing.assert_allclose(_u_archive(local, r), _u_library(scalar[-1], r),
                                   rtol=1e-12, atol=1e-12)
    for l, projector in semi.items():
        np.testing.assert_allclose(_u_archive(projector, r), _u_library(scalar[l], r),
                                   rtol=1e-12, atol=1e-12)
    has_so = any(len(t) == 3 for channel in by_l.values() for terms in channel for t in terms)
    assert entry.get("x_pyscf_spin_orbit_omitted", False) is has_so


def _radial_integral(f, a, l):
    """∫ r^(2l+2) exp(-2 a r²) f(r) dr by adaptive quadrature (f may carry r⁻²)."""
    from scipy.integrate import quad
    return quad(lambda r: r ** (2 * l + 2) * np.exp(-2 * a * r * r) * f(r), 0.0, np.inf,
                epsabs=1e-14, epsrel=1e-13, limit=500)[0]


@pytest.mark.parametrize("ecp,symbol,charge", [("def2-svp", "I", -1), ("lanl2dz", "I", -1),
                                               ("stuttgart_rsc", "Ag", 1)])
def test_quadrature_reproduces_ecpscalar(ecp, symbol, charge):
    """One-centre ⟨g|V|g⟩/⟨g|g⟩ from the archive's projectors equals PySCF's
    ECPscalar for primitives of l = 0 … 3 (for def2 and LANL2DZ, l = 3 sees
    only the local part; Stuttgart RSC Ag has an f channel and no local part)."""
    from qvf.adapters.pyscf import _ecp_entries
    exps = [3.0, 0.7, 0.2]
    basis = {symbol: [[l, [a, 1.0]] for l in range(4) for a in exps]}
    mol = gto.M(atom=f"{symbol} 0 0 0", basis=basis, ecp=ecp, charge=charge, verbose=0)
    (entry,) = _ecp_entries(mol)
    local, semi = _channels(entry)
    v, s = mol.intor("ECPscalar"), mol.intor("int1e_ovlp")
    ao = 0
    for l in range(4):
        for a in exps:
            def u(r, l=l):
                total = _u_archive(local, r) if local else 0.0
                return total + (_u_archive(semi[l], r) if l in semi else 0.0)
            expected = _radial_integral(u, a, l) / _radial_integral(lambda r: 1.0, a, l)
            got = v[ao, ao] / s[ao, ao]
            assert got == pytest.approx(expected, rel=1e-10, abs=1e-12), (l, a)
            ao += 2 * l + 1


def _rebuild(basis, atoms, charge=0, spin=0):
    """A PySCF molecule from the archive's basis and ECP payload alone."""
    labels = [f"{a['symbol']}{i + 1}" for i, a in enumerate(atoms)]
    atom = [[label, [x * BOHR_PER_ANGSTROM for x in a["position"]]]
            for label, a in zip(labels, atoms)]
    shells = {label: [] for label in labels}
    for shell in basis["shells"]:
        shells[labels[shell["center"]]].append(
            [shell["l"]] + [[e, c] for e, c in zip(shell["exponents"], shell["coefficients"])])
    ecp = {}
    for entry in basis.get("ecp", []):
        channels = []
        for p in entry["projectors"]:
            by_power = [[] for _ in range(max(p["powers"]) + 1)]
            for a, c, n in zip(p["exponents"], p["coefficients"], p["powers"]):
                by_power[n].append([a, c])
            channels.append([-1 if p.get("local") else p["l"], by_power])
        ecp[labels[entry["center"]]] = [entry["n_core_electrons"], channels]
    return gto.M(atom=atom, unit="Bohr", basis=shells, ecp=ecp, charge=charge,
                 spin=spin, cart=not basis.get("pure", True), verbose=0)


@pytest.mark.parametrize("name,atom,basis,ecp,spin", [
    ("hi-def2svp", HI, "def2-svp", "def2-svp", 0),
    ("i2-lanl2dz", "I 0 0 0; I 0 0 2.666", "lanl2dz", "lanl2dz", 0),
    ("agh-stuttgart-rsc", "Ag 0 0 0; H 0 0 1.618", "def2-svp", "stuttgart_rsc", 0),
    ("hi-crenbl", HI, "crenbl", "crenbl", 0),
])
def test_rebuilt_molecule_reproduces_scf(tmp_path, name, atom, basis, ecp, spin):
    mol = _mol(atom, basis, ecp, spin=spin)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    mf.run()
    _, payload, atoms = _write(tmp_path, mf, name)
    assert payload["n_electrons"] == mol.nelectron
    assert all("projectors" in e for e in payload["ecp"])
    rebuilt = _rebuild(payload, atoms, spin=spin)
    assert rebuilt.nelectron == mol.nelectron
    np.testing.assert_allclose(rebuilt.intor("ECPscalar"), mol.intor("ECPscalar"),
                               atol=1e-10)
    again = scf.RHF(rebuilt)
    again.conv_tol = 1e-11
    again.kernel(dm0=mf.make_rdm1())
    assert again.e_tot == pytest.approx(mf.e_tot, abs=1e-8)


def test_projectors_disabled_keeps_label_and_core(tmp_path):
    mf = scf.RHF(_mol(HI, "def2-svp", "def2-svp")).run()
    _, payload, _ = _write(tmp_path, mf, "hi-noproj", ecp_projectors=False)
    assert payload["ecp"] == [{"center": 1, "n_core_electrons": 28, "label": "def2-svp"}]


def test_add_scf_option_and_td_ntos(tmp_path):
    from qvf.adapters.pyscf import write_td
    mf = scf.RHF(_mol(HI, "def2-svp", "def2-svp")).run(conv_tol=1e-10)
    td = mf.TDA().set(nstates=2).run()
    path = tmp_path / "hi-tda.qvf"
    ids = write_td(td, str(path))
    assert validate_qvf(str(path))["findings"] == []
    with Reader(str(path)) as reader:
        for sid in [ids["wavefunction"], *ids["nto_holes"], *ids["nto_particles"]]:
            (entry,) = reader.json(sid, "basis")["ecp"]
            assert entry["projectors"]
    path = tmp_path / "hi-tda-noproj.qvf"
    ids = write_td(td, str(path), ecp_projectors=False)
    with Reader(str(path)) as reader:
        assert "projectors" not in reader.json(ids["nto_holes"][0], "basis")["ecp"][0]

    import qvf
    path = tmp_path / "hi-add.qvf"
    with qvf.Writer(str(path), program="PySCF", version="test", calculation="rhf") as writer:
        add_scf(writer, mf, ecp_projectors=False)
    with Reader(str(path)) as reader:
        assert "projectors" not in reader.json("wf", "basis")["ecp"][0]
