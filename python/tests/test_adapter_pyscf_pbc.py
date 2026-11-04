"""WP-62: the PySCF adapter for periodic Γ-point SCF and DFT.

Small synthetic cells (diamond in its two-atom primitive cell, a hydrogen
atom in a cubic box) with GTH pseudopotentials and a coarse FFT mesh. Every
check reads the archive back through the ``qvf`` reader:

* the reference validator reports nothing, ``E_NORMALIZATION`` included, so
  the stored ``overlap`` is the AO overlap of spec § 5.3;
* the orbitals are orthonormal and ``tr(P S)`` is the electron count in the
  lattice-summed Γ overlap, computed twice: by numerical integration over
  the cell with an AO evaluator written from spec Appendix A alone, and by
  PySCF from a cell rebuilt from the archive;
* the lattice vectors and fractional coordinates round-trip;
* ``energies.components.nuclear_repulsion`` equals an Ewald sum written
  here, from the archive's lattice, positions and effective charges, and
  ``total = electronic + nuclear_repulsion``;
* the cell rebuilt from the archive alone reproduces ``e_tot``, from the
  stored density and by a fresh SCF.

Tests that need PySCF are skipped when it is not importable.
"""
import importlib.util
import itertools
import math

import pytest

from qvf import Reader, Writer, validate_qvf

HAVE_PYSCF = importlib.util.find_spec("pyscf") is not None
requires_pyscf = pytest.mark.skipif(not HAVE_PYSCF, reason="PySCF is not installed")

DIAMOND = dict(atom="C 0 0 0; C 0.8917 0.8917 0.8917",
               a="0 1.7834 1.7834; 1.7834 0 1.7834; 1.7834 1.7834 0")
H_BOX = dict(atom="H 0.1 0.2 0.3", a="4 0 0; 0 4 0; 0 0 4", spin=1)
MESH = [15, 15, 15]


def _cell(geometry):
    from pyscf.pbc import gto
    return gto.M(basis="gth-szv", pseudo="gth-pade", mesh=MESH, verbose=0, **geometry)


def _make(name):
    from pyscf.pbc import dft, scf
    builders = {
        "diamond_rhf": lambda: scf.RHF(_cell(DIAMOND)),
        "diamond_rks": lambda: dft.RKS(_cell(DIAMOND), xc="pbe"),
        "h_uhf": lambda: scf.UHF(_cell(H_BOX)),
        "h_uks": lambda: dft.UKS(_cell(H_BOX), xc="lda,vwn"),
    }
    mf = builders[name]()
    mf.conv_tol = 1e-10
    mf.kernel()
    assert mf.converged
    return mf


CASES = ["diamond_rhf", "diamond_rks", "h_uhf", "h_uks"]


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    """Each case run once and written once: ``{case: (mf, path)}``."""
    if not HAVE_PYSCF:
        pytest.skip("PySCF is not installed")
    from qvf.adapters.pyscf import write_scf
    root = tmp_path_factory.mktemp("pyscf_pbc")
    out = {}
    for name in CASES:
        mf = _make(name)
        path = root / f"{name}.qvf"
        write_scf(mf, path)
        out[name] = (mf, path)
    return out


def _load(path):
    import numpy as np
    with Reader(path) as reader:
        wf = reader.sections_of("wavefunction.gto")[0]
        data = {"manifest": reader.raw_manifest,
                "structure": reader.json("structure", "structure"),
                "basis": reader.json(wf, "basis"),
                "meta": reader.json(wf, "mo_metadata"),
                "overlap": reader.tensor(wf, "overlap").numpy(),
                "spec": reader.json("job_spec", "spec"),
                "props": {role: reader.tensor("atom_properties", role).numpy()
                          for role in reader.section("atom_properties").members}}
        for role in wf.members:
            if role.startswith(("mo_coefficients", "density_matrix")):
                data[role] = reader.tensor(wf, role).numpy()
    if data["meta"]["spin"] == "unrestricted":
        data["sets"] = [(data[f"mo_coefficients_{ch}"], data["meta"][ch]["occupations"])
                        for ch in ("alpha", "beta")]
        data["density"] = data["density_matrix_alpha"] + data["density_matrix_beta"]
    else:
        data["sets"] = [(data["mo_coefficients"], data["meta"]["occupations"])]
        data["density"] = data["density_matrix"]
    data["sets"] = [(c, np.asarray(n)) for c, n in data["sets"]]
    return data


def _geometry_bohr(data):
    """Lattice rows and positions in bohr, through the format's factor."""
    import numpy as np
    from qvf.units import BOHR_PER_ANGSTROM
    lattice = np.asarray(data["structure"]["lattice_vectors"]) * BOHR_PER_ANGSTROM
    positions = np.asarray([a["position"] for a in data["structure"]["atoms"]]) * BOHR_PER_ANGSTROM
    return lattice, positions


def _effective_charges(data):
    import numpy as np
    z = np.asarray([a["atomic_number"] for a in data["structure"]["atoms"]], dtype=float)
    for entry in data["basis"].get("ecp", []):
        z[entry["center"]] -= entry["n_core_electrons"]
    return z


# -- independent references -------------------------------------------------------

def _gamma_overlap_numerical(data, n_grid=30, cutoff=13.0):
    """``Σ_T ⟨χ_μ | χ_ν(r − T)⟩`` by quadrature over one cell.

    The periodic Γ-point AOs ``Σ_T χ_μ(r − T)`` are evaluated with the
    Appendix A evaluator of ``test_adapter_pyscf`` on a uniform, shifted
    grid. Their integral over the cell is the lattice-summed overlap. Images
    farther than ``cutoff`` bohr from every grid point are left out; the
    most diffuse GTH-SZV exponent (0.12 bohr⁻²) is below 1e-8 there.
    """
    import numpy as np
    from test_adapter_pyscf import _spec_ao_values
    lattice, positions = _geometry_bohr(data)
    axis = (np.arange(n_grid) + 0.37) / n_grid  # no point on a nucleus
    frac = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), -1).reshape(-1, 3)
    points = frac @ lattice
    middle = points.mean(axis=0)
    radius = np.linalg.norm(points - middle, axis=1).max()
    heights = 1.0 / np.linalg.norm(np.linalg.inv(lattice), axis=0)  # plane spacings
    reach = [int(math.ceil((cutoff + radius) / h)) + 1 for h in heights]
    phi = 0.0
    for n in itertools.product(*(range(-k, k + 1) for k in reach)):
        centers = positions + np.asarray(n, dtype=float) @ lattice
        if np.linalg.norm(centers - middle, axis=1).min() > cutoff + radius:
            continue
        phi = phi + _spec_ao_values(data["basis"], centers, points)
    weight = abs(np.linalg.det(lattice)) / len(points)
    return weight * phi.T @ phi


def _ewald(lattice, positions, charges, eta=1.0, cut=10.0):
    """Ewald energy per cell of point charges in a neutralizing background."""
    import numpy as np
    from scipy.special import erfc
    volume = abs(np.linalg.det(lattice))
    real = 0.0
    reach = [int(math.ceil(cut / eta / np.linalg.norm(v))) + 1 for v in lattice]
    for n in itertools.product(*(range(-k, k + 1) for k in reach)):
        shift = np.asarray(n, dtype=float) @ lattice
        d = positions[:, None, :] - positions[None, :, :] + shift
        r = np.linalg.norm(d, axis=-1)
        mask = r > 1e-12
        real += 0.5 * (np.outer(charges, charges)[mask] * erfc(eta * r[mask]) / r[mask]).sum()
    recip = 2 * math.pi * np.linalg.inv(lattice).T
    reach = [int(math.ceil(2 * cut * eta / np.linalg.norm(b))) + 1 for b in recip]
    reciprocal = 0.0
    for m in itertools.product(*(range(-k, k + 1) for k in reach)):
        if not any(m):
            continue
        g = np.asarray(m, dtype=float) @ recip
        g2 = g @ g
        structure_factor = (charges * np.exp(1j * positions @ g)).sum()
        reciprocal += (2 * math.pi / volume * math.exp(-g2 / (4 * eta ** 2)) / g2
                       * abs(structure_factor) ** 2)
    self_term = -eta / math.sqrt(math.pi) * (charges ** 2).sum()
    background = -math.pi / (2 * volume * eta ** 2) * charges.sum() ** 2
    return real + reciprocal + self_term + background


def _rebuild(data):
    """A PySCF cell from the archive alone: structure, basis shells, labels, mesh."""
    from pyscf.pbc import gto
    lattice, positions = _geometry_bohr(data)
    atoms = data["structure"]["atoms"]
    labels = [f"{a['symbol']}{i + 1}" for i, a in enumerate(atoms)]
    basis = {label: [] for label in labels}
    for shell in data["basis"]["shells"]:
        basis[labels[shell["center"]]].append(
            [shell["l"]] + [[e, c] for e, c in zip(shell["exponents"], shell["coefficients"])])
    pseudo = {labels[e["center"]]: e["label"] for e in data["basis"].get("ecp", [])}
    options = data["spec"]["options"]
    assert options["pseudo"] == "gth-pade"
    return gto.M(atom=[[label, tuple(p)] for label, p in zip(labels, positions)],
                 a=lattice, unit="Bohr", basis=basis, pseudo=pseudo,
                 mesh=options["mesh"], spin=data["spec"]["multiplicity"] - 1, charge=data["spec"]["charge"], cart=not data["basis"].get("pure", True),
                 verbose=0)


def _rebuild_scf(data, cell):
    from pyscf.pbc import dft, scf
    spec = data["spec"]
    method, functional = spec["method"], spec["functional"]
    cls = {"rhf": scf.RHF, "uhf": scf.UHF, "rks": dft.RKS, "uks": dft.UKS}[method]
    mf = cls(cell) if functional is None else cls(cell, xc=functional)
    mf.exxdiv = spec["options"]["exxdiv"]
    assert spec["options"]["integrals"] == "fftdf"
    mf.conv_tol = 1e-10
    return mf


# -- checks ----------------------------------------------------------------------

@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_reference_validator_reports_nothing(written, case):
    assert validate_qvf(written[case][1])["findings"] == []


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_periodic_structure_and_round_trip(written, case):
    """Lattice rows and fractional coordinates of the cell are recovered."""
    import numpy as np
    mf, path = written[case]
    data = _load(path)
    structure = data["structure"]
    assert structure["pbc"] == [True, True, True]
    lattice, positions = _geometry_bohr(data)
    cell = mf.mol
    assert np.abs(lattice - cell.lattice_vectors()).max() < 1e-12
    frac = positions @ np.linalg.inv(lattice)
    frac_ref = cell.atom_coords() @ np.linalg.inv(cell.lattice_vectors())
    assert np.abs(frac - frac_ref).max() < 1e-12
    rebuilt = _rebuild(data)
    assert np.abs(rebuilt.lattice_vectors() - cell.lattice_vectors()).max() < 1e-12
    assert np.abs(rebuilt.atom_coords() - cell.atom_coords()).max() < 1e-12


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_gamma_metadata_and_pseudopotentials(written, case):
    mf, path = written[case]
    data = _load(path)
    meta, basis = data["meta"], data["basis"]
    assert meta["k_point"] == [0.0, 0.0, 0.0]
    assert "symmetry_labels" not in meta
    natm = mf.mol.natm
    cores = {"C": 2, "H": 0}
    assert basis["ecp"] == [
        {"center": i, "n_core_electrons": cores[mf.mol.atom_symbol(i)], "label": "gth-pade"}
        for i in range(natm)]
    assert "projectors" not in basis["ecp"][0]
    z = sum(a["atomic_number"] for a in data["structure"]["atoms"])
    assert basis["n_electrons"] == z - sum(e["n_core_electrons"] for e in basis["ecp"]) \
        == mf.mol.nelectron
    spec = data["spec"]
    assert spec["job_type"] == "periodic" and spec["kpoints"] == [1, 1, 1]
    assert spec["options"]["mesh"] == MESH


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_stored_overlap_is_the_ao_overlap(written, case):
    """``overlap`` is spec § 5.3's ⟨χ_μ|χ_ν⟩, not the lattice sum."""
    import numpy as np
    mf, path = written[case]
    data = _load(path)
    from qvf.adapters.pyscf import ao_map
    expected = ao_map(mf.mol).overlap(mf.mol.intor_symmetric("int1e_ovlp"))
    assert np.abs(data["overlap"] - expected).max() < 1e-12
    assert np.allclose(np.diag(data["overlap"]), 1.0, atol=1e-10)


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_gamma_metric_orthonormality_and_electron_count(written, case):
    """``Cᵀ S_Γ C = 1`` and ``tr(P S_Γ) = N`` with two independent ``S_Γ``."""
    import numpy as np
    from qvf.adapters.pyscf import ao_map
    mf, path = written[case]
    data = _load(path)
    numerical = _gamma_overlap_numerical(data)
    rebuilt = _rebuild(data)
    analytic = ao_map(rebuilt).overlap(rebuilt.pbc_intor("int1e_ovlp", hermi=1))
    assert np.abs(numerical - analytic).max() < 1e-6
    for s_gamma, tol in ((analytic, 1e-8), (numerical, 1e-5)):
        for c, occ in data["sets"]:
            assert np.abs(c @ s_gamma @ c.T - np.eye(len(c))).max() < tol
            density = (c.T * occ) @ c
            assert abs(np.trace(density @ s_gamma) - occ.sum()) < tol
        assert abs(np.trace(data["density"] @ s_gamma) - mf.mol.nelectron) < tol
    # The stored overlap alone does not give the electron count of a cell.
    assert abs(np.trace(data["density"] @ data["overlap"]) - mf.mol.nelectron) > 1e-2


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_density_and_populations(written, case):
    import numpy as np
    from qvf.adapters.pyscf import ao_map
    mf, path = written[case]
    data = _load(path)
    amap = ao_map(mf.mol)
    dm = np.asarray(mf.make_rdm1())
    total = dm[0] + dm[1] if dm.ndim == 3 else dm
    assert np.abs(data["density"] - amap.density(total)).max() < 1e-12
    z_eff = _effective_charges(data)
    charges = data["props"]["mulliken_charge"]
    assert abs(charges.sum() - (z_eff.sum() - mf.mol.nelectron)) < 1e-8
    assert abs(data["props"]["loewdin_charge"].sum() - charges.sum()) < 1e-8
    # PySCF's Mulliken analysis with the lattice-summed overlap.
    from pyscf.scf import hf, uhf
    population = uhf.mulliken_pop if dm.ndim == 3 else hf.mulliken_pop
    reference = population(mf.mol, dm, s=mf.get_ovlp(), verbose=0)[1]
    assert np.abs(charges - reference).max() < 1e-8
    if data["meta"]["spin"] == "unrestricted":
        spin = data["props"]["spin_population"]
        n_alpha, n_beta = (data["meta"][ch]["occupations"] for ch in ("alpha", "beta"))
        assert abs(spin.sum() - (sum(n_alpha) - sum(n_beta))) < 1e-8


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_energies_and_ewald_nuclear_energy(written, case):
    mf, path = written[case]
    data = _load(path)
    energies = data["manifest"]["energies"]
    total = energies["total"]["value"]
    nuclear = energies["components"]["nuclear_repulsion"]["value"]
    electronic = energies["components"]["electronic"]["value"]
    assert total == mf.e_tot == data["manifest"]["provenance"]["scf_energy"]["value"]
    assert abs(total - (electronic + nuclear)) < 1e-12
    lattice, positions = _geometry_bohr(data)
    reference = _ewald(lattice, positions, _effective_charges(data))
    assert abs(nuclear - reference) < 1e-9
    assert abs(nuclear - mf.mol.energy_nuc()) < 1e-12


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_rebuilt_cell_reproduces_the_energy(written, case):
    """The archive alone suffices to repeat the Γ-point calculation."""
    import numpy as np
    from qvf.adapters.pyscf import ao_map
    mf, path = written[case]
    data = _load(path)
    rebuilt = _rebuild(data)
    again = _rebuild_scf(data, rebuilt)
    amap = ao_map(rebuilt)
    inverse = np.linalg.inv(amap.matrix)
    back = (lambda p: inverse @ p @ inverse.T)
    if data["meta"]["spin"] == "unrestricted":
        dm = np.stack([back(data["density_matrix_alpha"]), back(data["density_matrix_beta"])])
    else:
        dm = back(data["density_matrix"])
    total = data["manifest"]["energies"]["total"]["value"]
    assert abs(again.energy_tot(dm=dm) - total) < 1e-8
    again.kernel()
    assert again.converged
    assert abs(again.e_tot - total) < 1e-7


@requires_pyscf
def test_low_dimensional_structure(tmp_path):
    """PySCF's periodic axes are the first ``dimension`` lattice vectors."""
    from pyscf.pbc import gto
    from qvf.adapters.pyscf import add_structure
    slab = gto.M(atom="H 0 0 0; H 0 0 0.74", a="3 0 0; 0 3 0; 0 0 10", dimension=2,
                 basis="gth-szv", pseudo="gth-pade", mesh=[11, 11, 41], verbose=0)
    path = tmp_path / "slab.qvf"
    with Writer(path, program="PySCF", version="test") as writer:
        add_structure(writer, slab)
    with Reader(path) as reader:
        structure = reader.json("structure", "structure")
    assert structure["pbc"] == [True, True, False]
    assert validate_qvf(path)["findings"] == []


@requires_pyscf
def test_unsupported_periodic_inputs_are_refused(tmp_path):
    import numpy as np
    from pyscf.pbc import scf
    from qvf.adapters.pyscf import write_scf, write_td
    cell = _cell(H_BOX)
    with pytest.raises(NotImplementedError, match="k-point sampling"):
        write_scf(scf.KUHF(cell, kpts=cell.make_kpts([1, 1, 1])).run(), tmp_path / "k.qvf")
    kpt = cell.make_kpts([2, 1, 1])[1]
    shifted = scf.UHF(cell, kpt=kpt).run()
    with pytest.raises(NotImplementedError, match="not the Γ point"):
        write_scf(shifted, tmp_path / "kpt.qvf")
    with pytest.raises(NotImplementedError, match="RHF, UHF, RKS and UKS"):
        write_scf(scf.ROHF(cell).run(), tmp_path / "rohf.qvf")
    complex_mf = scf.UHF(cell).run()
    complex_mf.mo_coeff = np.asarray(complex_mf.mo_coeff, dtype=complex)
    with pytest.raises(NotImplementedError, match="complex"):
        write_scf(complex_mf, tmp_path / "complex.qvf")
    closed = scf.RHF(_cell(DIAMOND)).run()
    with pytest.raises(NotImplementedError, match="periodic response"):
        write_td(closed.TDA().run(nstates=1), tmp_path / "td.qvf")
    assert not any(tmp_path.iterdir())
