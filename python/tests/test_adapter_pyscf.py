"""WP-62 tranche 1: the PySCF adapter for molecular SCF and DFT (#23).

The archives are generated here from synthetic molecules. Every check reads
them back through the ``qvf`` reader and compares with PySCF:

* the reference validator reports nothing, ``E_NORMALIZATION`` included;
* ``C S Cᵀ = 1`` for every orbital set, with the stored ``overlap``;
* the populations add up to the electron count (and the spin populations to
  ``n_α − n_β``);
* the stored densities equal ``Σ n C Cᵀ`` and PySCF's density in the QVF AO
  basis;
* an evaluator written from spec Appendix A alone (``N_i`` primitives, the
  A.2 order, real spherical harmonics from SciPy without the Condon–Shortley
  phase) reproduces every PySCF orbital and the density on random points, so
  AO order, phases and normalization are checked end to end.

Tests that need PySCF are skipped when it is not importable.
"""
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

from qvf import Reader, Writer, validate_qvf

HAVE_PYSCF = importlib.util.find_spec("pyscf") is not None
requires_pyscf = pytest.mark.skipif(not HAVE_PYSCF, reason="PySCF is not installed")
PYTHON_DIR = Path(__file__).resolve().parents[1]

WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
OH = "O 0 0 0; H 0.12 0.05 0.96"


def test_import_without_pyscf_and_clear_error():
    """The module imports without PySCF; using it names the extra to install."""
    code = (
        "import sys\n"
        "sys.modules['pyscf'] = None\n"
        "import qvf.adapters.pyscf as adapter\n"
        "assert adapter.real_solid_harmonics(1)\n"
        "try:\n"
        "    adapter.write_scf(object(), 'never.qvf')\n"
        "except ImportError as exc:\n"
        "    assert \"qvf-writer[pyscf]\" in str(exc), exc\n"
        "else:\n"
        "    raise SystemExit('no ImportError')\n")
    result = subprocess.run([sys.executable, "-c", code], cwd=PYTHON_DIR,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# -- an evaluator written from the specification alone --------------------------

def _spec_norm(alpha, l):
    """Appendix A.1 N_i, written out independently of the toolkit."""
    double_factorial = math.prod(range(1, 2 * l, 2)) if l else 1
    return ((2 * alpha / math.pi) ** 0.75 * (4 * alpha) ** (l / 2)
            / math.sqrt(double_factorial))


def _real_spherical(l, m, theta, phi):
    """Orthonormal real spherical harmonic: cos(mφ) for m > 0, sin(|m|φ) for
    m < 0, no Condon–Shortley phase (SciPy's complex harmonic includes it)."""
    import numpy as np
    from scipy.special import sph_harm_y
    y = sph_harm_y(l, abs(m), theta, phi)
    if m == 0:
        return y.real
    sign = -1.0 if abs(m) % 2 else 1.0
    return math.sqrt(2) * sign * (y.real if m > 0 else y.imag) * np.ones_like(theta)


def _spec_ao_values(basis, centers, points):
    """QVF AO values ``[n_points, n_ao]`` from a ``basis`` payload (bohr)."""
    import numpy as np
    columns = []
    pure_default = basis.get("pure", True)
    for shell in basis["shells"]:
        l = shell["l"]
        d = points - centers[shell["center"]]
        r2 = (d ** 2).sum(axis=1)
        radial = sum(c * _spec_norm(a, l) * np.exp(-a * r2)
                     for a, c in zip(shell["exponents"], shell["coefficients"]))
        if shell.get("pure", pure_default):
            r = np.sqrt(r2)
            theta = np.arccos(np.clip(d[:, 2] / r, -1.0, 1.0))
            phi = np.arctan2(d[:, 1], d[:, 0])
            scale = math.sqrt(4 * math.pi / (2 * l + 1))  # the norm of x^l
            for m in range(-l, l + 1):
                columns.append(radial * r ** l * scale * _real_spherical(l, m, theta, phi))
        else:
            for lx in range(l, -1, -1):
                for ly in range(l - lx, -1, -1):
                    lz = l - lx - ly
                    columns.append(radial * d[:, 0] ** lx * d[:, 1] ** ly * d[:, 2] ** lz)
    return np.stack(columns, axis=1)


# -- archives -------------------------------------------------------------------

def _make(name):
    from pyscf import dft, gto, scf
    water = dict(atom=WATER, verbose=0)
    oh = dict(atom=OH, spin=1, verbose=0)
    builders = {
        "rhf_sph": lambda: scf.RHF(gto.M(basis="def2-svp", **water)),
        "rhf_cart": lambda: scf.RHF(gto.M(basis="def2-svp", cart=True, **water)),
        "rks_sph": lambda: dft.RKS(gto.M(basis="cc-pvdz", **water), xc="pbe0"),
        "rks_cart": lambda: dft.RKS(gto.M(basis="cc-pvdz", cart=True, **water), xc="b3lyp"),
        "uhf_sph": lambda: scf.UHF(gto.M(basis="def2-svp", **oh)),
        "uks_cart": lambda: dft.UKS(gto.M(basis="def2-svp", cart=True, **oh), xc="pbe"),
        "rohf_sph": lambda: scf.ROHF(gto.M(basis="cc-pvdz", **oh)),
    }
    mf = builders[name]()
    if hasattr(mf, "grids"):
        mf.grids.level = 2
    mf.conv_tol = 1e-10
    mf.kernel()
    assert mf.converged
    return mf


CASES = ["rhf_sph", "rhf_cart", "rks_sph", "rks_cart", "uhf_sph", "uks_cart", "rohf_sph"]


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    """Each case run once and written once: ``{case: (mf, path)}``."""
    if not HAVE_PYSCF:
        pytest.skip("PySCF is not installed")
    from qvf.adapters.pyscf import write_scf
    out = {}
    root = tmp_path_factory.mktemp("pyscf")
    for name in CASES:
        mf = _make(name)
        path = root / f"{name}.qvf"
        write_scf(mf, path)
        out[name] = (mf, path)
    return out


def _load(path):
    """Read an archive back through the qvf reader."""
    with Reader(path) as reader:
        manifest = reader.raw_manifest
        wf = reader.sections_of("wavefunction.gto")[0]
        data = {"manifest": manifest,
                "structure": reader.json("structure", "structure"),
                "basis": reader.json(wf, "basis"),
                "meta": reader.json(wf, "mo_metadata"),
                "overlap": reader.tensor(wf, "overlap").numpy(),
                "members": dict(wf.members),
                "spec": reader.json("job_spec", "spec"),
                "props": {role: reader.tensor("atom_properties", role).numpy()
                          for role in reader.section("atom_properties").members}}
        for role in data["members"]:
            if role.startswith(("mo_coefficients", "density_matrix")):
                data[role] = reader.tensor(wf, role).numpy()
    return data


def _orbital_sets(data):
    meta = data["meta"]
    if meta["spin"] == "unrestricted":
        return [(data[f"mo_coefficients_{ch}"], meta[ch]["occupations"])
                for ch in ("alpha", "beta")]
    return [(data["mo_coefficients"], meta["occupations"])]


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_reference_validator_reports_nothing(written, case):
    report = validate_qvf(written[case][1])
    assert report["ok"], report["findings"]
    assert report["findings"] == []


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_orbitals_orthonormal_with_stored_overlap(written, case):
    import numpy as np
    data = _load(written[case][1])
    for c, _ in _orbital_sets(data):
        assert np.abs(c @ data["overlap"] @ c.T - np.eye(c.shape[0])).max() < 1e-9


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_overlap_diagonal_is_appendix_a1(written, case):
    """Unit diagonal for spherical and axial AOs; 1/3 for d_xy (Appendix A.1)."""
    import numpy as np
    data = _load(written[case][1])
    diagonal = np.diag(data["overlap"])
    if data["basis"]["pure"]:
        assert np.abs(diagonal - 1).max() < 1e-12
    else:
        ao = 0
        for shell in data["basis"]["shells"]:
            l = shell["l"]
            for lx in range(l, -1, -1):
                for ly in range(l - lx, -1, -1):
                    lz = l - lx - ly
                    want = (math.prod(range(1, 2 * lx, 2)) * math.prod(range(1, 2 * ly, 2))
                            * math.prod(range(1, 2 * lz, 2)) / math.prod(range(1, 2 * l, 2)))
                    assert diagonal[ao] == pytest.approx(want, rel=1e-12)
                    ao += 1
        assert any(s["l"] == 2 for s in data["basis"]["shells"])


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_density_round_trip(written, case):
    """Stored densities = Σ n C Cᵀ = PySCF's density mapped into QVF AOs."""
    import numpy as np
    from qvf.adapters.pyscf import ao_map
    mf, path = written[case]
    data = _load(path)
    amap = ao_map(mf.mol)
    dm = np.asarray(mf.make_rdm1())
    if data["meta"]["spin"] == "unrestricted":
        pairs = [("density_matrix_alpha", dm[0]), ("density_matrix_beta", dm[1])]
    else:
        pairs = [("density_matrix", dm if dm.ndim == 2 else dm[0] + dm[1])]
    assert len(pairs) == len(_orbital_sets(data))
    for (role, reference), (c, occ) in zip(pairs, _orbital_sets(data)):
        from_orbitals = (c.T * np.asarray(occ)) @ c
        assert np.abs(data[role] - from_orbitals).max() < 1e-10
        assert np.abs(data[role] - amap.density(reference)).max() < 1e-10
    total = sum(data[role] for role, _ in pairs)
    assert np.trace(total @ data["overlap"]) == pytest.approx(mf.mol.nelectron, abs=1e-8)


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_spec_evaluator_reproduces_pyscf_orbitals(written, case):
    """Every orbital and the density, on points, from Appendix A alone."""
    import numpy as np
    from qvf.units import BOHR_PER_ANGSTROM
    mf, path = written[case]
    mol = mf.mol
    data = _load(path)
    centers = np.array([a["position"] for a in data["structure"]["atoms"]]) * BOHR_PER_ANGSTROM
    assert np.abs(centers - mol.atom_coords()).max() < 1e-9
    rng = np.random.default_rng(7)
    points = rng.uniform(-2.5, 2.5, size=(300, 3)) + centers.mean(axis=0)
    ours = _spec_ao_values(data["basis"], centers, points)
    theirs = mol.eval_gto("GTOval_cart" if mol.cart else "GTOval_sph", points)
    mo = mf.mo_coeff if data["meta"]["spin"] == "unrestricted" else [mf.mo_coeff]
    for (c, _), reference in zip(_orbital_sets(data), mo):
        psi_qvf = ours @ c.T
        psi_pyscf = theirs @ reference
        assert np.abs(psi_qvf - psi_pyscf).max() < 1e-9 * max(1.0, np.abs(psi_pyscf).max())
    dm = np.asarray(mf.make_rdm1())
    total = dm if dm.ndim == 2 else dm[0] + dm[1]
    stored = sum(data[r] for r in data["members"] if r.startswith("density_matrix"))
    rho_qvf = np.einsum("pi,ij,pj->p", ours, stored, ours)
    rho_pyscf = np.einsum("pi,ij,pj->p", theirs, total, theirs)
    assert np.abs(rho_qvf - rho_pyscf).max() < 1e-9 * max(1.0, rho_pyscf.max())


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_populations(written, case):
    import numpy as np
    from scipy.linalg import sqrtm
    mf, path = written[case]
    mol = mf.mol
    data = _load(path)
    props = data["props"]
    z_eff = mol.atom_charges()
    for role in ("mulliken_charge", "loewdin_charge"):
        assert (z_eff - props[role]).sum() == pytest.approx(mol.nelectron, abs=1e-8)
    assert props["mulliken_charge"] == pytest.approx(mf.mulliken_pop(verbose=0)[1], abs=1e-9)
    open_shell = data["meta"]["spin"] != "restricted"
    assert ("spin_population" in props) == open_shell
    if open_shell:
        assert props["spin_population"].sum() == pytest.approx(mol.spin, abs=1e-8)
    if not mol.cart:  # PySCF's own AOs: Löwdin populations computed independently
        s = mol.intor_symmetric("int1e_ovlp")
        half = sqrtm(s).real
        dm = np.asarray(mf.make_rdm1())
        total = dm if dm.ndim == 2 else dm[0] + dm[1]
        q = np.diag(half @ total @ half)
        pops = np.array([q[a0:a1].sum() for _, _, a0, a1 in mol.aoslice_by_atom()])
        assert props["loewdin_charge"] == pytest.approx(z_eff - pops, abs=1e-9)


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_manifest_job_spec_and_energies(written, case):
    import pyscf
    mf, path = written[case]
    data = _load(path)
    manifest, spec = data["manifest"], data["spec"]
    method = {"rhf": "rhf", "rks": "rks", "uhf": "uhf", "uks": "uks", "rohf": "rohf"}[case.split("_")[0]]
    assert manifest["source"]["program"] == "PySCF"
    assert manifest["source"]["version"] == pyscf.__version__
    prov = manifest["provenance"]
    assert prov["method"] == spec["method"] == method
    assert prov["scf_energy"] == {"value": mf.e_tot, "units": "hartree"}
    assert prov["run_status"] == "converged" and prov["scf_converged"] is True
    assert prov["n_electrons"] == data["basis"]["n_electrons"] == mf.mol.nelectron
    assert prov["software"][0] == {"name": "PySCF", "version": pyscf.__version__}
    energies = manifest["energies"]
    assert energies["method"] == method and energies["total"]["value"] == mf.e_tot
    components = energies["components"]
    assert components["nuclear_repulsion"]["value"] == pytest.approx(mf.energy_nuc(), abs=1e-12)
    assert (components["electronic"]["value"] + components["nuclear_repulsion"]["value"]
            == pytest.approx(mf.e_tot, abs=1e-10))
    assert spec["job_type"] == "molecular" and spec["structure_ref"] == "structure"
    assert spec["basis"] == mf.mol.basis
    assert spec["charge"] == 0 and spec["multiplicity"] == mf.mol.spin + 1
    assert "functional" in spec  # null and absent differ (spec § 5.9)
    if method.endswith("ks"):
        assert spec["functional"] == prov["functional"] == mf.xc.lower()
    else:
        assert spec["functional"] is None and "functional" not in prov
    section = {s["id"]: s for s in manifest["sections"]}
    assert section["wf"]["structure_ref"] == data["basis"]["structure_ref"] == "structure"
    assert section["atom_properties"]["structure_ref"] == "structure"
    assert section["structure"]["multiplicity"] == mf.mol.spin + 1
    spin = {"r": "restricted", "u": "unrestricted"}.get(method[0])
    assert data["meta"]["spin"] == ("restricted_open" if method.startswith("ro") else spin)


# -- the AO map ------------------------------------------------------------------

def test_real_solid_harmonics_match_the_spherical_harmonics():
    """The adapter's table is the evaluator's definition, l = 0..5."""
    pytest.importorskip("scipy")
    import numpy as np
    from qvf.adapters.pyscf import cartesian_powers, real_solid_harmonics
    rng = np.random.default_rng(3)
    d = rng.normal(size=(64, 3))
    r = np.linalg.norm(d, axis=1)
    theta, phi = np.arccos(d[:, 2] / r), np.arctan2(d[:, 1], d[:, 0])
    for l in range(6):
        table = np.asarray(real_solid_harmonics(l))
        monomials = np.stack([d[:, 0] ** a * d[:, 1] ** b * d[:, 2] ** c
                              for a, b, c in cartesian_powers(l)], axis=1)
        scale = math.sqrt(4 * math.pi / (2 * l + 1))
        for m in range(-l, l + 1):
            want = r ** l * scale * _real_spherical(l, m, theta, phi)
            assert np.allclose(monomials @ table[:, m + l], want, atol=1e-12, rtol=1e-12)
    # A.2 anchors: p is (y, z, x) and d_{-2} is sqrt(3) xy.
    assert real_solid_harmonics(1) == [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    assert real_solid_harmonics(2)[1][0] == pytest.approx(math.sqrt(3))


@requires_pyscf
def test_ao_map_is_the_documented_permutation():
    """p: (x, y, z) -> m = -1, 0, +1; l >= 2 spherical: identity; Cartesian:
    identity order scaled by the l >= 2 axial norm sqrt(4π/(2l+1))."""
    import numpy as np
    from pyscf import gto
    from qvf.adapters.pyscf import ao_map
    basis = {"Ne": [[0, [1.5, 1.0]], [1, [1.1, 1.0]], [2, [0.9, 1.0]], [3, [0.7, 1.0]],
                    [4, [0.5, 1.0]]]}
    sph = ao_map(gto.M(atom="Ne 0 0 0", basis=basis, verbose=0))
    t = sph.matrix
    assert np.array_equal(t[1:4, 1:4], np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]))
    assert np.array_equal(t[4:, 4:], np.eye(t.shape[0] - 4))
    cart = ao_map(gto.M(atom="Ne 0 0 0", basis=basis, cart=True, verbose=0))
    diagonal = np.diag(cart.matrix)
    assert np.allclose(diagonal[:4], 1.0)
    start = 4
    for l in (2, 3, 4):
        n = (l + 1) * (l + 2) // 2
        assert np.allclose(diagonal[start:start + n], math.sqrt(4 * math.pi / (2 * l + 1)))
        start += n
    assert np.count_nonzero(cart.matrix - np.diag(diagonal)) == 0


@requires_pyscf
def test_general_contractions_split_in_pyscf_ao_order():
    from pyscf import gto
    from qvf.adapters.pyscf import ao_map
    mol = gto.M(atom=WATER, basis="cc-pvdz", verbose=0)
    assert any(mol.bas_nctr(i) > 1 for i in range(mol.nbas))
    amap = ao_map(mol)
    assert len(amap.shells) == sum(mol.bas_nctr(i) for i in range(mol.nbas))
    assert amap.n_ao == mol.nao_nr()
    assert all(0.0 not in s["coefficients"] for s in amap.shells)
    labels = mol.ao_labels(fmt=False)
    assert [labels[i][0] for i in range(mol.nao_nr())] == amap.ao_atoms()


# -- further behaviour -------------------------------------------------------------

@requires_pyscf
def test_symmetry_frame_and_labels(tmp_path):
    from pyscf import gto, scf
    from qvf.adapters.pyscf import write_scf
    mol = gto.M(atom=WATER, basis="def2-svp", symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    path = tmp_path / "sym.qvf"
    write_scf(mf, path)
    assert validate_qvf(path)["findings"] == []
    data = _load(path)
    labels = data["meta"]["symmetry_labels"]
    assert labels[:3] == ["A1", "A1", "B2"] and set(labels) <= {"A1", "A2", "B1", "B2"}


@requires_pyscf
def test_ecp_core_count_and_electron_count(tmp_path):
    from pyscf import gto, scf
    from qvf.adapters.pyscf import write_scf
    mol = gto.M(atom="H 0 0 0; I 0 0 1.61", basis="def2-svp", ecp="def2-svp", verbose=0)
    mf = scf.RHF(mol).run()
    path = tmp_path / "hi.qvf"
    write_scf(mf, path)
    report = validate_qvf(path)
    assert report["findings"] == [], report["findings"]
    data = _load(path)
    (entry,) = data["basis"]["ecp"]
    assert {k: entry[k] for k in ("center", "n_core_electrons", "label")} == {
        "center": 1, "n_core_electrons": 28, "label": "def2-svp"}
    assert entry["projectors"]  # checked in test_adapter_pyscf_ecp.py
    assert data["basis"]["n_electrons"] == mol.nelectron == 26
    assert (mol.atom_charges() - data["props"]["mulliken_charge"]).sum() == pytest.approx(26)


@requires_pyscf
def test_add_scf_into_an_open_writer(tmp_path):
    import pyscf
    from pyscf import gto, scf
    from qvf.adapters.pyscf import add_scf
    mf = scf.RHF(gto.M(atom=WATER, basis="def2-svp", verbose=0)).run()
    path = tmp_path / "custom.qvf"
    with Writer(path, program="PySCF", version=pyscf.__version__) as writer:
        ids = add_scf(writer, mf, structure_id="mol", wavefunction_id="orbitals",
                      populations=False, job_spec=False, density=False)
    assert ids == {"structure": "mol", "wavefunction": "orbitals"}
    with Reader(path) as reader:
        kinds = [(s.id, s.kind) for s in reader.sections]
        members = set(reader.section("orbitals").members)
        assert reader.json("orbitals", "basis")["structure_ref"] == "mol"
    assert kinds == [("mol", "structure"), ("orbitals", "wavefunction.gto")]
    assert "overlap" in members and not any(m.startswith("density") for m in members)
    assert validate_qvf(path)["findings"] == []


@requires_pyscf
def test_unsupported_inputs_are_refused(tmp_path):
    from pyscf import gto, scf
    from qvf.adapters.pyscf import write_scf
    mol = gto.M(atom=WATER, basis="sto-3g", verbose=0)
    with pytest.raises(ValueError, match="run it first"):
        write_scf(scf.RHF(mol), tmp_path / "unrun.qvf")
    with pytest.raises(NotImplementedError, match="later tranches"):
        write_scf(scf.GHF(mol).run(), tmp_path / "ghf.qvf")
    ghost = gto.M(atom="O 0 0 0; ghost-H 0 0 1", basis="sto-3g", spin=0, charge=0,
                  verbose=0)
    with pytest.raises(ValueError, match="ghost"):
        write_scf(scf.RHF(ghost).run(), tmp_path / "ghost.qvf")
    assert not any(tmp_path.iterdir())


@requires_pyscf
def test_written_archive_json_is_plain(written):
    """Payloads hold plain JSON numbers (no NumPy scalars leaked)."""
    with Reader(written["uks_cart"][1]) as reader:
        meta = reader.json("wf", "mo_metadata")
    json.dumps(meta)
    assert isinstance(meta["alpha"]["energies"][0], float)
