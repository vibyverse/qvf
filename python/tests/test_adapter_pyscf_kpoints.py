"""WP-62: the PySCF adapter for k-point sampled periodic SCF and DFT.

Small cells (diamond in its two-atom primitive cell, H₂ and a hydrogen atom
in a cubic box) with GTH-SZV, GTH-Pade and a coarse FFT mesh, on k meshes of
at most eight points. Every check reads the archive back through the ``qvf``
reader:

* the reference validator reports nothing;
* S(k) built from the archive alone with the Bloch sum of spec § 5.3.1,
  ``Σ_T e^{+ik·T} ⟨χ_μ(r)|χ_ν(r − T)⟩``, equals PySCF's
  ``cell.pbc_intor('int1e_ovlp', kpts=...)``, and the opposite phase and a
  phase that includes the atomic positions do not. The lattice sum is built
  from molecular overlaps between the cell and its images, and checked once
  more by quadrature over the cell with the Appendix A AO evaluator of
  ``test_adapter_pyscf``;
* ``C(k)ᴴ S(k) C(k) = 1`` at every k and spin, and
  ``Σ_k w_k tr(P(k) S(k))`` is the electron count of a cell;
* the k-points round-trip between PySCF's absolute k (bohr⁻¹) and the spec's
  fractional coordinates, and ``mesh`` and ``shift`` describe them;
* a cell rebuilt from the archive, with densities built from the stored
  coefficients, reproduces ``e_tot``.

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
H2_BOX = dict(atom="H 0 0 0; H 0.74 0.1 0.05", a="3 0 0; 0 3 0; 0 0 3")
H_BOX = dict(atom="H 0.1 0.2 0.3", a="4 0 0; 0 4 0; 0 0 4")
MESH = [15, 15, 15]
EXPLICIT = [[0.0, 0.0, 0.0], [0.2, 0.1, 0.0], [-0.3, 0.0, 0.4]]  # fractional, no mesh


def _cell(geometry, **extra):
    from pyscf.pbc import gto
    return gto.M(basis="gth-szv", pseudo="gth-pade", mesh=MESH, verbose=0,
                 **geometry, **extra)


def _make(name):
    from pyscf.pbc import dft, scf
    if name == "diamond_krhf_222":
        cell = _cell(DIAMOND)
        mf = scf.KRHF(cell, kpts=cell.make_kpts([2, 2, 2]))
    elif name == "diamond_krks_311":
        cell = _cell(DIAMOND)
        mf = dft.KRKS(cell, kpts=cell.make_kpts([3, 1, 1]), xc="pbe")
    elif name == "h2_kuhf_explicit":
        cell = _cell(H2_BOX)
        mf = scf.KUHF(cell, kpts=cell.get_abs_kpts(EXPLICIT))
    elif name == "h_kuks_twisted":
        # PySCF's KUHF spin is that of the supercell: one unpaired electron per cell.
        cell = _cell(H_BOX, spin=2)
        mf = dft.KUKS(cell, kpts=cell.make_kpts([2, 1, 1], scaled_center=[0.1, 0, 0]),
                      xc="lda,vwn")
    else:
        raise KeyError(name)
    mf.conv_tol = 1e-10
    mf.kernel()
    assert mf.converged
    return mf


CASES = ["diamond_krhf_222", "diamond_krks_311", "h2_kuhf_explicit", "h_kuks_twisted"]


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    """Each case run once and written once: ``{case: (mf, path)}``."""
    if not HAVE_PYSCF:
        pytest.skip("PySCF is not installed")
    from qvf.adapters.pyscf import write_scf
    root = tmp_path_factory.mktemp("pyscf_kpoints")
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
        (wf,) = reader.sections_of("wavefunction.bloch")
        assert not reader.sections_of("wavefunction.gto")
        data = {"manifest": reader.raw_manifest,
                "structure": reader.json("structure", "structure"),
                "basis": reader.json(wf, "basis"),
                "kpoints": reader.json(wf, "kpoints"),
                "meta": reader.json(wf, "mo_metadata"),
                "spec": reader.json("job_spec", "spec")}
        raw = reader.tensor(wf, "mo_coefficients").numpy()
        data["coefficients"] = raw[..., 0] + 1j * raw[..., 1]  # [n_k, n_spin, n_mo, n_ao]
        data["energies"] = reader.tensor(wf, "energies").numpy()
        data["occupations"] = reader.tensor(wf, "occupations").numpy()
    data["frac"] = np.asarray(data["kpoints"]["kpoints"])
    data["weights"] = np.asarray(data["kpoints"]["weights"])
    return data


def _geometry_bohr(data):
    """Lattice rows and positions in bohr, through the format's factor."""
    import numpy as np
    from qvf.units import BOHR_PER_ANGSTROM
    lattice = np.asarray(data["structure"]["lattice_vectors"]) * BOHR_PER_ANGSTROM
    positions = np.asarray([a["position"] for a in data["structure"]["atoms"]]) * BOHR_PER_ANGSTROM
    return lattice, positions


def _basis_by_label(data):
    atoms = data["structure"]["atoms"]
    labels = [f"{a['symbol']}{i + 1}" for i, a in enumerate(atoms)]
    basis = {label: [] for label in labels}
    for shell in data["basis"]["shells"]:
        basis[labels[shell["center"]]].append(
            [shell["l"]] + [[e, c] for e, c in zip(shell["exponents"], shell["coefficients"])])
    return labels, basis


def _complex_overlap(amap, s):
    """A complex PySCF AO matrix in the QVF AO basis (the map is real)."""
    return amap.overlap(s.real) + 1j * amap.overlap(s.imag)


# -- S(k) from the archive alone ---------------------------------------------------

def _translation_blocks(data, cutoff=20.0):
    """``{n: ⟨χ_μ(r)|χ_ν(r − T_n)⟩}`` in the QVF AO basis, molecular integrals only.

    A molecule of the cell's atoms with the archive's shells, and one of all
    images within ``cutoff`` bohr; the cross overlap is split per image.
    Overlaps between primitives farther apart than ``cutoff`` are below 1e-10
    for GTH-SZV (``exp(−αβR²/(α + β))``, smallest exponent 0.12 bohr⁻²).
    """
    import numpy as np
    from pyscf import gto
    from qvf.adapters.pyscf import ao_map
    lattice, positions = _geometry_bohr(data)
    labels, basis = _basis_by_label(data)
    pbc = data["structure"]["pbc"]
    heights = 1.0 / np.linalg.norm(np.linalg.inv(lattice), axis=0)
    span = np.linalg.norm(positions[:, None] - positions[None], axis=-1).max()
    reach = [int(math.ceil((cutoff + span) / h)) + 1 if p else 0
             for h, p in zip(heights, pbc)]

    def molecule(atoms):
        from pyscf.data.elements import charge
        mol = gto.Mole()
        mol.atom, mol.basis, mol.unit, mol.verbose = atoms, basis, "Bohr", 0
        mol.spin = sum(charge(a[0].rstrip("0123456789")) for a in atoms) % 2
        mol.build()
        return mol

    home = molecule([[label, tuple(p)] for label, p in zip(labels, positions)])
    amap = ao_map(home)
    images, atoms = [], []
    for n in itertools.product(*(range(-k, k + 1) for k in reach)):
        shifted = positions + np.asarray(n, dtype=float) @ lattice
        if np.linalg.norm(shifted[:, None] - positions[None], axis=-1).min() > cutoff:
            continue
        images.append(n)
        atoms += [[label, tuple(p)] for label, p in zip(labels, shifted)]
    many = molecule(atoms)
    cross = gto.intor_cross("int1e_ovlp", home, many)
    n_ao = home.nao_nr()
    return {n: amap.overlap(cross[:, i * n_ao:(i + 1) * n_ao]) for i, n in enumerate(images)}


def _bloch_overlap(blocks, frac, sign=+1, positions_phase=None):
    """``S(k) = Σ_T e^{sign·2πi k·n} S_T`` at fractional ``frac``."""
    import numpy as np
    out = sum(np.exp(sign * 2j * math.pi * (np.asarray(n) @ frac)) * s
              for n, s in blocks.items())
    if positions_phase is not None:  # the other convention: e^{ik·(T + R_ν − R_μ)}
        out = positions_phase.conj()[:, None] * out * positions_phase[None, :]
    return out


def _quadrature_overlap(data, frac, n_grid=24, cutoff=13.0):
    """``∫_cell φ̄_μk φ_νk`` with ``φ_μk = Σ_T e^{ik·T} χ_μ(r − T)`` by quadrature.

    The AOs are evaluated by the Appendix A evaluator of
    ``test_adapter_pyscf`` on a uniform shifted grid of one cell.
    """
    import numpy as np
    from test_adapter_pyscf import _spec_ao_values
    lattice, positions = _geometry_bohr(data)
    axis = (np.arange(n_grid) + 0.37) / n_grid
    points = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), -1).reshape(-1, 3) @ lattice
    middle = points.mean(axis=0)
    radius = np.linalg.norm(points - middle, axis=1).max()
    heights = 1.0 / np.linalg.norm(np.linalg.inv(lattice), axis=0)
    reach = [int(math.ceil((cutoff + radius) / h)) + 1 for h in heights]
    phi = [0.0 for _ in frac]
    for n in itertools.product(*(range(-k, k + 1) for k in reach)):
        centers = positions + np.asarray(n, dtype=float) @ lattice
        if np.linalg.norm(centers - middle, axis=1).min() > cutoff + radius:
            continue
        values = _spec_ao_values(data["basis"], centers, points)
        for i, k in enumerate(frac):
            phi[i] = phi[i] + np.exp(2j * math.pi * (np.asarray(n) @ k)) * values
    weight = abs(np.linalg.det(lattice)) / len(points)
    return [weight * p.conj().T @ p for p in phi]


@pytest.fixture(scope="module")
def overlaps(written):
    """``{case: [S(k)]}`` built from each archive alone."""
    out = {}
    for case, (_, path) in written.items():
        data = _load(path)
        blocks = _translation_blocks(data)
        out[case] = [_bloch_overlap(blocks, k) for k in data["frac"]]
    return out


# -- checks ------------------------------------------------------------------------

@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_reference_validator_reports_nothing(written, case):
    assert validate_qvf(written[case][1])["findings"] == []


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_kpoints_round_trip(written, case):
    """Spec fractional k ↔ PySCF absolute k, with the archive's own lattice."""
    import numpy as np
    from qvf.adapters.pyscf import absolute_kpoints, fractional_kpoints
    mf, path = written[case]
    data = _load(path)
    lattice, _ = _geometry_bohr(data)
    reciprocal = 2 * math.pi * np.linalg.inv(lattice).T  # rows b_j, a_i · b_j = 2π δ_ij
    assert np.abs(lattice @ reciprocal.T - 2 * math.pi * np.eye(3)).max() < 1e-12
    absolute = data["frac"] @ reciprocal
    assert np.abs(absolute - np.asarray(mf.kpts)).max() < 1e-12
    assert np.abs(absolute_kpoints(mf.mol, data["frac"]) - mf.kpts).max() < 1e-12
    assert np.abs(fractional_kpoints(mf.mol, mf.kpts) - data["frac"]).max() == 0.0
    assert np.abs(mf.mol.get_scaled_kpts(mf.kpts) - data["frac"]).max() < 1e-12
    assert np.allclose(data["weights"], 1.0 / len(mf.kpts)) and abs(data["weights"].sum() - 1) < 1e-12


@requires_pyscf
def test_mesh_and_shift(written):
    import numpy as np
    kp = {case: _load(path)["kpoints"] for case, (_, path) in written.items()}
    assert kp["diamond_krhf_222"]["mesh"] == [2, 2, 2] and "shift" not in kp["diamond_krhf_222"]
    assert kp["diamond_krks_311"]["mesh"] == [3, 1, 1] and "shift" not in kp["diamond_krks_311"]
    assert kp["diamond_krks_311"]["kpoints"][1] == [1 / 3, 0.0, 0.0]
    assert "mesh" not in kp["h2_kuhf_explicit"] and "shift" not in kp["h2_kuhf_explicit"]
    assert np.abs(np.asarray(kp["h2_kuhf_explicit"]["kpoints"]) - EXPLICIT).max() < 1e-12
    assert kp["h_kuks_twisted"]["mesh"] == [2, 1, 1]
    assert kp["h_kuks_twisted"]["shift"] == [0.1, 0.0, 0.0]
    specs = {case: _load(path)["spec"] for case, (_, path) in written.items()}
    assert specs["diamond_krhf_222"]["kpoints"] == [2, 2, 2]
    assert "kpoint_shift" not in specs["diamond_krhf_222"]["options"]
    assert specs["h_kuks_twisted"]["options"]["kpoint_shift"] == [0.1, 0.0, 0.0]
    assert "kpoints" not in specs["h2_kuhf_explicit"]
    assert np.abs(np.asarray(specs["h2_kuhf_explicit"]["options"]["kpoints_fractional"])
                  - EXPLICIT).max() < 1e-12
    assert all(s["job_type"] == "periodic" for s in specs.values())


@requires_pyscf
def test_kpoint_mesh_detection():
    from qvf.adapters.pyscf import kpoint_mesh
    assert kpoint_mesh([[0, 0, 0]]) == ([1, 1, 1], None)
    assert kpoint_mesh([[-0.25, 0, 0], [0.25, 0, 0]]) == ([2, 1, 1], [0.25, 0.0, 0.0])
    assert kpoint_mesh([[0, 0, 0], [0.5, 0, 0], [1.0, 0, 0]]) == (None, None)  # repeats Γ
    assert kpoint_mesh([[0, 0, 0], [0.3, 0, 0]]) == (None, None)
    assert kpoint_mesh([[0, 0, 0], [0.5, 0.5, 0]]) == (None, None)  # 2 of 2 × 2


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_bloch_overlap_convention_matches_pyscf(written, overlaps, case):
    """The spec's S(k) from the archive alone is PySCF's ``pbc_intor`` S(k)."""
    import numpy as np
    from qvf.adapters.pyscf import ao_map
    mf, path = written[case]
    data = _load(path)
    cell = mf.mol
    amap = ao_map(cell)
    pyscf_sk = cell.pbc_intor("int1e_ovlp", hermi=1, kpts=mf.kpts)
    blocks = _translation_blocks(data)
    _, positions = _geometry_bohr(data)
    ao_positions = positions[amap.ao_atoms()]
    complex_somewhere = False
    for k, frac in enumerate(data["frac"]):
        reference = _complex_overlap(amap, pyscf_sk[k])
        assert np.abs(overlaps[case][k] - reference).max() < 1e-9
        if np.abs(reference.imag).max() > 1e-3:
            complex_somewhere = True
            # Neither the opposite phase nor one with the atomic positions fits.
            assert np.abs(_bloch_overlap(blocks, frac, sign=-1) - reference).max() > 1e-2
            kvec = frac @ (2 * math.pi * np.linalg.inv(_geometry_bohr(data)[0]).T)
            with_positions = _bloch_overlap(blocks, frac,
                                            positions_phase=np.exp(1j * ao_positions @ kvec))
            assert np.abs(with_positions - reference).max() > 1e-2
    # A 2 × 2 × 2 Γ-centred mesh has real phases only, and a single s AO has a
    # real S(k) since S(T) = S(−T). The other two cases test the sign.
    assert complex_somewhere == (case in ("diamond_krks_311", "h2_kuhf_explicit"))


@requires_pyscf
def test_bloch_overlap_by_quadrature(written, overlaps):
    """The lattice sum once more, by quadrature with the Appendix A AO evaluator."""
    import numpy as np
    data = _load(written["diamond_krks_311"][1])
    numerical = _quadrature_overlap(data, data["frac"])
    for k in range(len(numerical)):
        assert np.abs(numerical[k] - overlaps["diamond_krks_311"][k]).max() < 1e-4


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_orthonormality_and_electron_count(written, overlaps, case):
    """``C(k)ᴴ S(k) C(k) = 1`` and ``Σ_k w_k tr(P(k) S(k)) = N`` per cell."""
    import numpy as np
    mf, path = written[case]
    data = _load(path)
    n_k, n_spin, n_mo, n_ao = data["coefficients"].shape
    assert n_spin == (1 if data["meta"]["spin"] == "restricted" else 2)
    assert data["meta"]["n_mo"] == n_mo and data["meta"]["n_ao"] == n_ao == data["basis"]["n_ao"]
    total = 0.0
    for k in range(n_k):
        s = overlaps[case][k]
        for spin in range(n_spin):
            c = data["coefficients"][k, spin]  # rows are orbitals
            assert np.abs(c.conj() @ s @ c.T - np.eye(n_mo)).max() < 1e-8
            density = c.T @ np.diag(data["occupations"][k, spin]) @ c.conj()
            assert np.abs(density - density.conj().T).max() < 1e-12
            total += data["weights"][k] * np.trace(density @ s)
    assert abs(total.imag) < 1e-10
    assert abs(total.real - mf.mol.nelectron) < 1e-8
    assert data["basis"]["n_electrons"] == mf.mol.nelectron


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_energies_occupations_and_metadata(written, case):
    import numpy as np
    mf, path = written[case]
    data = _load(path)
    meta = data["meta"]
    assert meta["complex"] is True and meta["k_point_convention"] == "fractional"
    assert "n_mo_per_k" not in meta
    restricted = meta["spin"] == "restricted"
    energies = [mf.mo_energy] if restricted else list(mf.mo_energy)
    occupations = [mf.mo_occ] if restricted else list(mf.mo_occ)
    for spin in range(len(energies)):
        for k in range(len(mf.kpts)):
            assert np.array_equal(data["energies"][k, spin], energies[spin][k])
            assert np.array_equal(data["occupations"][k, spin], occupations[spin][k])
    entries = data["basis"]["ecp"]
    cores = {"C": 2, "H": 0}
    assert entries == [{"center": i, "n_core_electrons": cores[mf.mol.atom_symbol(i)],
                        "label": "gth-pade"} for i in range(mf.mol.natm)]
    manifest = data["manifest"]
    assert manifest["energies"]["total"]["value"] == mf.e_tot
    nuclear = manifest["energies"]["components"]["nuclear_repulsion"]["value"]
    assert abs(nuclear - mf.energy_nuc()) < 1e-12
    multiplicity = {"h_kuks_twisted": 2}.get(case, 1)
    (structure,) = [s for s in manifest["sections"] if s["id"] == "structure"]
    assert structure["multiplicity"] == multiplicity
    assert data["spec"]["multiplicity"] == manifest["provenance"]["multiplicity"] == multiplicity
    assert data["structure"]["pbc"] == [True, True, True]


@requires_pyscf
@pytest.mark.parametrize("case", CASES)
def test_rebuilt_cell_reproduces_the_energy(written, case):
    """The stored orbitals in a cell rebuilt from the archive give ``e_tot``."""
    import numpy as np
    from pyscf.pbc import dft, gto, scf
    from qvf.adapters.pyscf import ao_map
    mf, path = written[case]
    data = _load(path)
    lattice, positions = _geometry_bohr(data)
    labels, basis = _basis_by_label(data)
    pseudo = {labels[e["center"]]: e["label"] for e in data["basis"]["ecp"]}
    spec = data["spec"]
    n_k = len(data["frac"])
    spin_per_cell = spec["multiplicity"] - 1
    cell = gto.M(atom=[[label, tuple(p)] for label, p in zip(labels, positions)],
                 a=lattice, unit="Bohr", basis=basis, pseudo=pseudo,
                 mesh=spec["options"]["mesh"], spin=spin_per_cell * n_k, verbose=0)
    kpts = data["frac"] @ (2 * math.pi * np.linalg.inv(lattice).T)
    cls = {"rhf": scf.KRHF, "uhf": scf.KUHF, "rks": dft.KRKS, "uks": dft.KUKS}[spec["method"]]
    again = cls(cell, kpts=kpts)
    if spec["functional"] is not None:
        again.xc = spec["functional"]
    again.exxdiv = spec["options"]["exxdiv"]
    inverse = np.linalg.inv(ao_map(cell).matrix)
    dms = []
    for spin in range(data["coefficients"].shape[1]):
        per_k = []
        for k in range(n_k):
            c = inverse @ data["coefficients"][k, spin].T  # PySCF [n_ao, n_mo]
            per_k.append((c * data["occupations"][k, spin]) @ c.conj().T)
        dms.append(np.asarray(per_k))
    dm = dms[0] if len(dms) == 1 else np.asarray(dms)
    assert abs(again.energy_tot(dm=dm) - mf.e_tot) < 1e-8


@requires_pyscf
def test_gamma_only_kpoint_object(tmp_path):
    """A one-point KRHF at Γ is written as wavefunction.bloch, not wavefunction.gto."""
    from pyscf.pbc import scf
    from qvf.adapters.pyscf import write_scf
    cell = _cell(H2_BOX)
    path = tmp_path / "gamma.qvf"
    write_scf(scf.KRHF(cell, kpts=cell.make_kpts([1, 1, 1])).run(), path)
    assert validate_qvf(path)["findings"] == []
    data = _load(path)
    assert data["kpoints"]["kpoints"] == [[0.0, 0.0, 0.0]] and data["kpoints"]["mesh"] == [1, 1, 1]


@requires_pyscf
def test_add_wavefunction_dispatches(tmp_path):
    from pyscf.pbc import scf
    from qvf.adapters.pyscf import add_bloch_wavefunction, add_structure, add_wavefunction
    cell = _cell(H2_BOX)
    mf = scf.KRHF(cell, kpts=cell.make_kpts([2, 1, 1])).run()
    path = tmp_path / "parts.qvf"
    with Writer(path, program="PySCF", version="test") as writer:
        add_structure(writer, cell)
        assert add_wavefunction(writer, mf, section_id="a") == "a"
        assert add_bloch_wavefunction(writer, mf, section_id="b") == "b"
        with pytest.raises(NotImplementedError, match="density"):
            add_wavefunction(writer, mf, section_id="c", density=True)
        with pytest.raises(TypeError, match="not a k-point sampled"):
            add_bloch_wavefunction(writer, scf.RHF(cell).run(), section_id="d")
    assert validate_qvf(path)["findings"] == []
    with Reader(path) as reader:
        assert [s.id for s in reader.sections_of("wavefunction.bloch")] == ["a", "b"]
        assert reader.json("a", "kpoints") == reader.json("b", "kpoints")


@requires_pyscf
def test_unsupported_kpoint_inputs_are_refused(tmp_path):
    import numpy as np
    from pyscf.pbc import gto, scf
    from qvf.adapters.pyscf import add_populations, write_scf
    target = tmp_path / "refused.qvf"
    h = _cell(H_BOX, spin=1)
    with pytest.raises(NotImplementedError, match="KRHF, KUHF, KRKS and KUKS"):
        write_scf(scf.KROHF(h, kpts=h.make_kpts([1, 1, 1])).run(), target)
    with pytest.raises(NotImplementedError, match="KRHF, KUHF, KRKS and KUKS"):
        write_scf(scf.KGHF(h, kpts=h.make_kpts([1, 1, 1])).run(), target)
    # Two unpaired electrons over four k-points is half a spin per cell.
    half = _cell(H_BOX, spin=2)
    with pytest.raises(NotImplementedError, match="integer spin per cell"):
        write_scf(scf.KUHF(half, kpts=half.make_kpts([4, 1, 1])).run(), target)
    sym_cell = _cell(H2_BOX, space_group_symmetry=True)
    symmetric = scf.KRHF(sym_cell, kpts=sym_cell.make_kpts([2, 1, 1], space_group_symmetry=True))
    symmetric.mo_coeff = symmetric.mo_occ = np.zeros(1)  # refused before it is read
    with pytest.raises(NotImplementedError, match="k-point symmetry"):
        write_scf(symmetric, target)
    slab = gto.M(atom="H 0 0 0; H 0 0 0.74", a="3 0 0; 0 3 0; 0 0 10", dimension=2,
                 basis="gth-szv", pseudo="gth-pade", mesh=[11, 11, 41], verbose=0)
    lifted = scf.KRHF(slab, kpts=slab.get_abs_kpts([[0, 0, 0], [0, 0, 0.25]]))
    lifted.mo_coeff = lifted.mo_occ = np.zeros(1)
    with pytest.raises(NotImplementedError, match="non-periodic axis 2"):
        write_scf(lifted, target)
    cell = _cell(H2_BOX)
    mf = scf.KRHF(cell, kpts=cell.make_kpts([2, 1, 1])).run()
    with pytest.raises(NotImplementedError, match="populations"):
        write_scf(mf, target, populations=True)
    with pytest.raises(NotImplementedError, match="density"):
        write_scf(mf, target, density=True)
    with Writer(tmp_path / "pop.qvf", program="PySCF", version="test") as writer:
        with pytest.raises(NotImplementedError, match="populations"):
            add_populations(writer, mf)
    (tmp_path / "pop.qvf").unlink()
    assert not target.exists()
