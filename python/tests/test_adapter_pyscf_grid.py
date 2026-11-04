"""WP-62: densities, spin densities and orbitals of the PySCF adapter on a grid.

Every archive is written from a small synthetic molecule on a coarse grid
and read back through the ``qvf`` reader. The checks are independent of the
adapter's own evaluation:

* the reference validator reports nothing;
* the density integrates to the electron count, the spin density to
  ``n_α − n_β`` and every orbital's square to one (grid quadrature, with
  tolerances measured for these grids; see ``docs/adapters/pyscf.md``);
* on sampled grid points the stored values equal ρ and ψ evaluated from the
  archive's own ``basis`` and ``density_matrix`` members through the spec
  Appendix A evaluator of ``test_adapter_pyscf.py``;
* the stored values equal PySCF's ``tools.cubegen`` output on the same grid,
  and the cube's grid header, imported through ``qvf.convert.import_cube``,
  is the archive's grid.

Tests that need PySCF are skipped when it is not importable.
"""
import importlib.util
import io

import pytest

from qvf import Reader, Writer, validate_qvf

HAVE_PYSCF = importlib.util.find_spec("pyscf") is not None
requires_pyscf = pytest.mark.skipif(not HAVE_PYSCF, reason="PySCF is not installed")

WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
OH = "O 0 0 0; H 0.12 0.05 0.96"
SPACING = 0.15  # bohr; with the default 4 bohr margin about 1.5e5 points for water


def _make(name):
    from pyscf import dft, gto, scf
    builders = {
        "rhf_sph": lambda: scf.RHF(gto.M(atom=WATER, basis="def2-svp", verbose=0)),
        "rks_cart": lambda: dft.RKS(gto.M(atom=WATER, basis="def2-svp", cart=True, verbose=0),
                                    xc="pbe"),
        "uhf_sph": lambda: scf.UHF(gto.M(atom=OH, basis="def2-svp", spin=1, verbose=0)),
        "rohf_sph": lambda: scf.ROHF(gto.M(atom=OH, basis="def2-svp", spin=1, verbose=0)),
    }
    mf = builders[name]()
    if hasattr(mf, "grids"):
        mf.grids.level = 2
    mf.conv_tol = 1e-10
    mf.kernel()
    assert mf.converged
    return mf


CASES = ["rhf_sph", "rks_cart", "uhf_sph", "rohf_sph"]


def _write(mf, path, **density_options):
    import pyscf
    from qvf.adapters.pyscf import add_density_grid, add_orbital_grid, add_scf
    with Writer(path, program="PySCF", version=pyscf.__version__) as writer:
        add_scf(writer, mf, populations=False, job_spec=False)
        ids = add_density_grid(writer, mf, **density_options)
        grid_options = {k: v for k, v in density_options.items()
                        if k in ("grid", "spacing", "shape", "margin")}
        ids["orbitals"] = add_orbital_grid(writer, mf, ["homo", "lumo", 0], **grid_options)
        if isinstance(mf.mo_coeff, tuple) or getattr(mf.mo_coeff, "ndim", 2) == 3:
            ids["orbitals"] += add_orbital_grid(writer, mf, ["homo"], spin="beta",
                                                **grid_options)
    return ids


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    """Each case run and written once: ``{case: (mf, path, ids)}``."""
    if not HAVE_PYSCF:
        pytest.skip("PySCF is not installed")
    root = tmp_path_factory.mktemp("pyscf_grid")
    out = {}
    for name in CASES:
        mf = _make(name)
        path = root / f"{name}.qvf"
        out[name] = (mf, path, _write(mf, path, spacing=SPACING))
    return out


def _volume(reader, sid):
    import numpy as np
    grid = reader.json(sid, "grid")
    data = reader.tensor(sid, "data").numpy().astype(np.float64)
    return grid, data


def _points(grid, flat_index):
    import numpy as np
    ijk = np.stack(np.unravel_index(flat_index, grid["shape"]), axis=1)
    return np.asarray(grid["origin"]) + ijk @ np.asarray(grid["voxel_vectors"])


def _voxel_volume(grid):
    import numpy as np
    return abs(np.linalg.det(np.asarray(grid["voxel_vectors"])))


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_reference_validator_reports_nothing(written, case):
    report = validate_qvf(written[case][1])
    assert report["ok"], report["findings"]
    assert report["findings"] == []


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_sections_kinds_and_links(written, case):
    mf, path, ids = written[case]
    open_shell = case != "rhf_sph" and case != "rks_cart"
    assert ("spin" in ids) == open_shell
    with Reader(path) as reader:
        sections = {s["id"]: s for s in reader.raw_manifest["sections"]}
        for sid in [ids["density"], *([ids["spin"]] if open_shell else []), *ids["orbitals"]]:
            section = sections[sid]
            assert section["structure_ref"] == "structure" and section["group"] == "wf"
            grid = reader.json(sid, "grid")
            assert grid["structure_ref"] == "structure"
            assert set(grid) == {"origin", "voxel_vectors", "shape", "structure_ref"}
            assert reader.tensor(sid, "data").dtype == "float32"
    kinds = {sid: sections[sid]["kind"] for sid in sections}
    assert kinds[ids["density"]] == "volume.density"
    assert all(kinds[sid] == "volume.orbital" for sid in ids["orbitals"])
    if open_shell:
        assert kinds[ids["spin"]] == "volume.spin"
    homo = int(mf.mol.nelectron // 2) - 1 if not open_shell else None
    if case == "rhf_sph":
        assert ids["orbitals"] == [f"mo_{homo}", f"mo_{homo + 1}", "mo_0"]
        assert sections[f"mo_{homo}"]["label"] == f"MO {homo} (HOMO)"
    if case == "uhf_sph":
        n_alpha, n_beta = mf.nelec
        assert ids["orbitals"] == [f"mo_alpha_{n_alpha - 1}", f"mo_alpha_{n_alpha}",
                                   "mo_alpha_0", f"mo_beta_{n_beta - 1}"]


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_grid_integrals(written, case):
    """∫ρ = N, ∫(ρ_α − ρ_β) = 2S and ∫ψ² = 1 by the grid sum Σ f ΔV.

    On a uniform grid the sum aliases the tight O 1s primitives (α up to
    2266 bohr⁻² in def2-SVP): at 0.15 bohr and a 4 bohr margin the measured
    |ΔN| is 0.012 to 0.015 e, twice the 1s orbital's error, and it changes
    sign with the spacing. Valence orbitals and the spin density are
    smooth (measured errors below 1e-4); the diffuse LUMO loses about 2 %
    outside the box. The bounds leave a factor of 2 to 5.
    """
    mf, path, ids = written[case]
    with Reader(path) as reader:
        grid, rho = _volume(reader, ids["density"])
        dv = _voxel_volume(grid)
        assert rho.sum() * dv == pytest.approx(mf.mol.nelectron, abs=0.03)
        assert rho.min() > -1e-6
        if "spin" in ids:
            _, spin = _volume(reader, ids["spin"])
            assert spin.sum() * dv == pytest.approx(mf.mol.spin, abs=1e-3)
        homo, lumo, core = ids["orbitals"][:3]
        bounds = {homo: 5e-4, core: 0.02, lumo: 0.05}
        for sid in ids["orbitals"][3:]:
            bounds[sid] = 5e-4  # the beta HOMO
        for sid, bound in bounds.items():
            _, psi = _volume(reader, sid)
            assert (psi ** 2).sum() * dv == pytest.approx(1.0, abs=bound), sid


@requires_pyscf
@pytest.mark.parametrize("atoms,spin", [("H 0 0 0; H 0 0 0.74", 0),
                                        ("H 0 0 0; H 0 0 0.9; H 0 0.8 0.45", 1)])
def test_grid_integral_without_a_tight_core(atoms, spin):
    """Hydrogen only (α ≤ 13 bohr⁻²): the grid sum is converged in the spacing
    and limited by the box. Measured at 0.25 bohr and a 6 bohr margin:
    |ΔN| and |Δ(2S)| below 5e-5 (H2 RHF, H3 UHF); 1e-3 to 1e-2 at 4 bohr."""
    import numpy as np
    from pyscf import gto, scf
    from qvf.adapters.pyscf import add_density_grid, add_scf
    mol = gto.M(atom=atoms, basis="def2-svp", spin=spin, verbose=0)
    mf = (scf.UHF if spin else scf.RHF)(mol).run(conv_tol=1e-10)
    buf = io.BytesIO()
    with Writer(buf, program="PySCF", version="test") as writer:
        add_scf(writer, mf, populations=False, job_spec=False)
        ids = add_density_grid(writer, mf, spacing=0.25, margin=6.0, dtype="float64")
    buf.seek(0)
    with Reader(buf) as reader:
        grid, rho = _volume(reader, ids["density"])
        dv = _voxel_volume(grid)
        assert abs(rho.sum() * dv - mol.nelectron) < 2e-4
        if spin:
            _, s = _volume(reader, ids["spin"])
            assert abs(s.sum() * dv - spin) < 2e-4
        assert np.isfinite(rho).all()


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_values_match_the_archive_basis_through_appendix_a(written, case):
    """Sampled points: stored ρ, spin and ψ equal the Appendix A evaluation of
    the archive's own basis, density matrices and coefficients."""
    import numpy as np
    from qvf.units import BOHR_PER_ANGSTROM
    from test_adapter_pyscf import _spec_ao_values
    mf, path, ids = written[case]
    with Reader(path) as reader:
        basis = reader.json("wf", "basis")
        meta = reader.json("wf", "mo_metadata")
        members = set(reader.section("wf").members)
        atoms = reader.json("structure", "structure")["atoms"]
        centers = np.array([a["position"] for a in atoms]) * BOHR_PER_ANGSTROM
        matrices = {r: reader.tensor("wf", r).numpy() for r in members
                    if r.startswith(("density_matrix", "mo_coefficients"))}
        grid, rho = _volume(reader, ids["density"])
        rng = np.random.default_rng(11)
        flat = rng.choice(rho.size, size=400, replace=False)
        # Include the densest points, near the nuclei, where errors would show.
        flat = np.unique(np.concatenate([flat, np.argsort(rho, axis=None)[-20:]]))
        ao = _spec_ao_values(basis, centers, _points(grid, flat))
        if meta["spin"] == "unrestricted":
            p_alpha, p_beta = matrices["density_matrix_alpha"], matrices["density_matrix_beta"]
        else:
            p_alpha = p_beta = None
        total = (p_alpha + p_beta) if p_alpha is not None else matrices["density_matrix"]
        want = np.einsum("pi,ij,pj->p", ao, total, ao)
        # float32 storage: relative 1e-6 of each value, plus a small floor.
        assert np.abs(rho.ravel()[flat] - want).max() <= 2e-6 * np.abs(want).max()
        assert np.allclose(rho.ravel()[flat], want, rtol=2e-6, atol=1e-7)
        if "spin" in ids:
            _, spin = _volume(reader, ids["spin"])
            if p_alpha is None:
                # ROHF: the archive holds the total density only; the spin
                # density is the singly occupied orbitals' density.
                c = matrices["mo_coefficients"]
                occ = np.asarray(meta["occupations"])
                singly = c[occ == 1]
                spin_dm = singly.T @ singly
            else:
                spin_dm = p_alpha - p_beta
            want_spin = np.einsum("pi,ij,pj->p", ao, spin_dm, ao)
            assert np.allclose(spin.ravel()[flat], want_spin, rtol=2e-6,
                               atol=2e-6 * np.abs(want_spin).max())
        sections = {s["id"]: s for s in reader.raw_manifest["sections"]}
        for sid in ids["orbitals"]:
            _, psi = _volume(reader, sid)
            row = int(sid.rsplit("_", 1)[1])
            beta = "_beta_" in sid
            c = matrices["mo_coefficients_beta" if beta else
                         ("mo_coefficients_alpha" if meta["spin"] == "unrestricted"
                          else "mo_coefficients")]
            want_psi = ao @ c[row]
            assert np.allclose(psi.ravel()[flat], want_psi, rtol=2e-6,
                               atol=2e-6 * np.abs(want_psi).max()), sections[sid]["label"]


@requires_pyscf
@pytest.mark.parametrize("case", ["rhf_sph", "uhf_sph"])
def test_same_grid_as_cubegen(tmp_path, written, case):
    """``tools.cubegen`` on ``shape`` and ``margin`` gives the archive's grid
    and values; the cube header, imported by ``qvf.convert``, is that grid."""
    import numpy as np
    from pyscf.tools import cubegen
    from qvf.convert import import_cube
    mf = written[case][0]
    shape, margin = (31, 35, 33), 3.0
    path = tmp_path / "grid.qvf"
    ids = _write(mf, path, shape=shape, margin=margin)
    dm = np.asarray(mf.make_rdm1())
    total = dm if dm.ndim == 2 else dm[0] + dm[1]
    cube = tmp_path / "density.cube"
    want = cubegen.density(mf.mol, str(cube), total, nx=shape[0], ny=shape[1], nz=shape[2],
                           margin=margin)
    imported = tmp_path / "density_cube.qvf"
    import_cube(cube, imported, kind="volume.density")
    with Reader(path) as reader, Reader(imported) as from_cube:
        grid, rho = _volume(reader, ids["density"])
        cube_grid, cube_rho = _volume(from_cube, "volume")
        assert grid["shape"] == list(shape) == cube_grid["shape"]
        # The cube header keeps 6 decimals; cubegen's own points are exact.
        assert np.allclose(grid["origin"], cube_grid["origin"], rtol=0, atol=1e-6)
        assert np.allclose(grid["voxel_vectors"], cube_grid["voxel_vectors"], rtol=0, atol=1e-6)
        points = _points(grid, np.arange(rho.size))
        reference = cubegen.Cube(mf.mol, *shape, margin=margin).get_coords()
        assert np.abs(points - reference).max() < 1e-12
        assert np.allclose(rho, want, rtol=2e-6, atol=1e-7 * want.max())
        # The cube text keeps 5 significant digits (PySCF writes %.5E).
        assert np.allclose(rho, cube_rho, rtol=1e-5, atol=1e-6 * want.max())
        homo_id = ids["orbitals"][0]
        _, psi = _volume(reader, homo_id)
    coeff = mf.mo_coeff[0] if isinstance(mf.mo_coeff, tuple) or np.ndim(mf.mo_coeff) == 3 \
        else mf.mo_coeff
    homo = int(homo_id.rsplit("_", 1)[1])
    want_psi = cubegen.orbital(mf.mol, str(tmp_path / "homo.cube"), np.asarray(coeff)[:, homo],
                               nx=shape[0], ny=shape[1], nz=shape[2], margin=margin)
    assert np.allclose(psi, want_psi, rtol=2e-6, atol=2e-6 * np.abs(want_psi).max())
    assert validate_qvf(path)["findings"] == []


@requires_pyscf
def test_default_grid():
    """The default box: 4 bohr margin, 0.2 bohr steps, centred on the nuclei.
    Water/def2-SVP on it: measured |ΔN| = 0.061 e (O 1s aliasing)."""
    import numpy as np
    from qvf.adapters.pyscf import (DEFAULT_GRID_MARGIN_BOHR, DEFAULT_GRID_SPACING_BOHR,
                                    add_density_grid, add_scf, molecular_grid)
    mf = _make("rhf_sph")
    mol = mf.mol
    grid = molecular_grid(mol)
    assert DEFAULT_GRID_MARGIN_BOHR == 4.0 and DEFAULT_GRID_SPACING_BOHR == 0.2
    coords = mol.atom_coords()
    v = np.asarray(grid["voxel_vectors"])
    assert np.array_equal(v, 0.2 * np.eye(3))
    far = np.asarray(grid["origin"]) + (np.asarray(grid["shape"]) - 1) * 0.2
    assert (np.asarray(grid["origin"]) <= coords.min(axis=0) - 4.0 + 1e-9).all()
    assert (far >= coords.max(axis=0) + 4.0 - 1e-9).all()
    assert np.allclose((np.asarray(grid["origin"]) + far) / 2,
                       (coords.min(axis=0) + coords.max(axis=0)) / 2)
    buf = io.BytesIO()
    with Writer(buf, program="PySCF", version="test") as writer:
        add_scf(writer, mf, populations=False, job_spec=False)
        ids = add_density_grid(writer, mf, dtype="float64")
    buf.seek(0)
    with Reader(buf) as reader:
        g, rho = _volume(reader, ids["density"])
        assert g["shape"] == grid["shape"] and g["origin"] == grid["origin"]
        assert reader.tensor(ids["density"], "data").dtype == "float64"
    assert abs(rho.sum() * _voxel_volume(g) - mol.nelectron) < 0.12


@requires_pyscf
def test_explicit_skewed_grid_and_chunking():
    """A non-orthogonal grid payload is honoured, and the block size does
    not change the stored values."""
    import numpy as np
    from qvf.adapters.pyscf import add_density_grid, add_orbital_grid, add_scf
    from test_adapter_pyscf import _spec_ao_values
    mf = _make("rhf_sph")
    grid = {"origin": [-3.0, -3.5, -3.2], "voxel_vectors": [[0.3, 0.0, 0.0], [0.1, 0.3, 0.0],
                                                           [0.0, 0.05, 0.3]],
            "shape": [20, 22, 21]}
    stored = []
    for chunk in (None, 97):
        buf = io.BytesIO()
        with Writer(buf, program="PySCF", version="test") as writer:
            structure = add_scf(writer, mf, populations=False, job_spec=False)
            ids = add_density_grid(writer, mf, grid=grid, chunk=chunk, dtype="float64")
            add_orbital_grid(writer, mf, [2], grid=grid, chunk=chunk, dtype="float64")
        buf.seek(0)
        with Reader(buf) as reader:
            g, rho = _volume(reader, ids["density"])
            _, psi = _volume(reader, "mo_2")
            assert g["voxel_vectors"] == grid["voxel_vectors"] and g["origin"] == grid["origin"]
            basis = reader.json(structure["wavefunction"], "basis")
            dm = reader.tensor("wf", "density_matrix").numpy()
            c = reader.tensor("wf", "mo_coefficients").numpy()
        stored.append((rho, psi))
    # Equal to the last bits that BLAS block sizes may change.
    for a, b in zip(*stored):
        assert np.allclose(a, b, rtol=1e-13, atol=1e-15 * np.abs(a).max())
    flat = np.arange(0, rho.size, 251)
    ao = _spec_ao_values(basis, mf.mol.atom_coords(), _points(grid, flat))
    assert np.allclose(rho.ravel()[flat], np.einsum("pi,ij,pj->p", ao, dm, ao),
                       rtol=1e-10, atol=1e-12)
    assert np.allclose(psi.ravel()[flat], ao @ c[2], rtol=1e-10, atol=1e-12)


@requires_pyscf
def test_refusals(tmp_path):
    import pyscf
    from qvf.adapters.pyscf import add_density_grid, add_orbital_grid, add_scf
    mf = _make("rhf_sph")
    with Writer(tmp_path / "x.qvf", program="PySCF", version=pyscf.__version__) as writer:
        with pytest.raises(ValueError, match="no structure section"):
            add_density_grid(writer, mf)
        add_scf(writer, mf, populations=False, job_spec=False)
        with pytest.raises(ValueError, match="closed-shell"):
            add_density_grid(writer, mf, spin=True)
        with pytest.raises(ValueError, match="wavefunction.gto"):
            add_density_grid(writer, mf, wavefunction_ref="nope")
        with pytest.raises(ValueError, match="outside"):
            add_orbital_grid(writer, mf, [mf.mo_coeff.shape[1]])
        with pytest.raises(TypeError):
            add_orbital_grid(writer, mf, [1.5])
        with pytest.raises(ValueError, match="one set"):
            add_orbital_grid(writer, mf, [0], spin="beta")
        with pytest.raises(ValueError, match="twice"):
            add_orbital_grid(writer, mf, [4, "homo"])
        with pytest.raises(ValueError, match="linearly dependent"):
            add_density_grid(writer, mf, grid={"origin": [0, 0, 0], "shape": [2, 2, 2],
                                               "voxel_vectors": [[1, 0, 0], [2, 0, 0],
                                                                 [0, 0, 1]]})
        with pytest.raises(ValueError, match="spacing"):
            add_density_grid(writer, mf, spacing=0.0)
        # wavefunction_ref=None writes no group, and the archive still validates.
        ids = add_density_grid(writer, mf, spacing=0.5, wavefunction_ref=None)
    with Reader(tmp_path / "x.qvf") as reader:
        section = next(s for s in reader.raw_manifest["sections"] if s["id"] == ids["density"])
    assert "group" not in section
    assert validate_qvf(tmp_path / "x.qvf")["findings"] == []
