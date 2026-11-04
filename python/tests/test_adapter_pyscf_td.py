"""Convert source TD results without changing their meaning or amplitudes."""
import copy

import pytest

np = pytest.importorskip("numpy")
pyscf = pytest.importorskip("pyscf")
from pyscf import dft, gto, lib, scf
from pyscf.tdscf.rhf import TDA, TDHF

from qvf import Reader, Writer, validate_qvf
from qvf.adapters.pyscf import add_td, method_name, write_td
from qvf.units import BOHR_PER_ANGSTROM, HARTREE_IN_EV


@pytest.fixture(scope="module")
def references():
    old_threads = lib.num_threads()
    lib.num_threads(1)
    mol = gto.M(atom="O 0 0 0; H 0 0.76 0.59; H 0 -0.76 0.59",
                basis="def2-svp", verbose=0, max_memory=256)
    rhf = scf.RHF(mol).run(conv_tol=1e-11)
    # A small response space avoids making conversion tests depend on a
    # large iterative TD eigensolver's convergence. Water still tests d AOs,
    # rectangular NTO spaces and non-contiguous frozen orbital indices.
    response = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="6-31g", verbose=0, max_memory=256)
    response_rhf = scf.RHF(response).run(conv_tol=1e-11)
    hybrid = dft.RKS(response, xc="pbe0")
    hybrid.grids.level = 0
    hybrid.kernel()
    pure = dft.RKS(response, xc="pbe")
    pure.grids.level = 0
    pure.kernel()
    cart = mol.copy()
    cart.cart = True
    cart.build()
    cartesian = scf.RHF(cart).run(conv_tol=1e-11)
    small = gto.M(atom=mol.atom, basis="sto-3g", verbose=0, max_memory=256)
    minimal = scf.RHF(small).run(conv_tol=1e-11)
    yield {"rhf": rhf, "response": response_rhf, "hybrid": hybrid, "pure": pure,
           "cartesian": cartesian, "minimal": minimal}
    lib.num_threads(old_threads)


@pytest.fixture(scope="module")
def results(references):
    variants = {
        "rhf-tda": references["rhf"].TDA(),
        "rhf-full": references["response"].TDHF(),
        "rks-tda": references["hybrid"].TDA(),
        "rks-full": references["hybrid"].TDDFT(),
        "casida": references["pure"].TDDFT(),
        "cartesian": references["cartesian"].TDA(),
        "minimal": references["minimal"].TDA(),
        "minimal-full": references["minimal"].TDHF(),
        "frozen-full": references["minimal"].TDHF(),
        "frozen": references["rhf"].TDA(),
        "triplet": references["rhf"].TDA(),
        "triplet-full": references["response"].TDHF(),
        "triplet-rks": references["hybrid"].TDDFT(),
    }
    variants["frozen"].frozen = [0, 7]  # includes a virtual MO: no contiguous-index shortcut
    variants["frozen-full"].frozen = [0, 5]
    for name in ("triplet", "triplet-full", "triplet-rks"):
        variants[name].singlet = False
    for name, td in variants.items():
        td.nstates = 2
        td.conv_tol = 1e-8
        td.kernel()
        assert np.asarray(td.converged).all(), name
    return variants


@pytest.mark.parametrize("case", ["rhf-tda", "rhf-full", "rks-tda", "rks-full",
                                  "casida", "cartesian", "minimal", "frozen", "triplet",
                                  "triplet-full", "triplet-rks", "minimal-full", "frozen-full"])
def test_roundtrip_preserves_source_and_linked_results(results, tmp_path, case):
    td = results[case]
    mf = td._scf
    before = [mf.mo_coeff.copy(), mf.mo_occ.copy(), mf.mo_energy.copy(), td.e.copy()]
    xy_before = [(x.copy(), np.array(y, copy=True)) for x, y in td.xy]
    geometry = mf.mol.atom_coords().copy()
    path = tmp_path / f"{case}.qvf"
    ids = write_td(td, path)
    report = validate_qvf(path)
    assert report["ok"] and not report["findings"], report
    for actual, expected in zip([mf.mo_coeff, mf.mo_occ, mf.mo_energy, td.e], before):
        np.testing.assert_array_equal(actual, expected)
    for (x, y), (old_x, old_y) in zip(td.xy, xy_before):
        np.testing.assert_array_equal(x, old_x)
        np.testing.assert_array_equal(y, old_y)
    np.testing.assert_array_equal(mf.mol.atom_coords(), geometry)

    with Reader(path) as reader:
        manifest = reader.raw_manifest
        assert manifest["source"]["program"] == "PySCF"
        assert manifest["source"]["version"] == pyscf.__version__
        assert manifest["provenance"]["method"] == method_name(mf)
        assert manifest["provenance"]["scf_energy"]["value"] == mf.e_tot
        assert manifest["energies"]["total"]["value"] == mf.e_tot
        assert manifest["energies"]["method"] == method_name(mf)
        positions = [a["position"] for a in reader.json(ids["structure"], "structure")["atoms"]]
        np.testing.assert_allclose(np.array(positions) * BOHR_PER_ANGSTROM, geometry, atol=1e-14)
        payload = reader.json(ids["excited_states"], "states")
        spectrum = reader.json(ids["spectrum"], "spectrum")
        tda = isinstance(td, TDA) and not isinstance(td, TDHF)
        assert payload["meta"]["tda"] is tda
        assert payload["meta"]["method"] == (
            "tddft" if isinstance(mf, dft.rks.KohnShamDFT) else "cis" if tda else "tdhf")
        assert payload["meta"]["n_states"] == len(td.e)
        assert payload["meta"]["reference_state"] == "SCF reference"
        assert spectrum["states_ref"] == ids["excited_states"]
        assert spectrum["x_units"] == "ev" and spectrum["y_units"] == "dimensionless"
        by_id = {s["id"]: s for s in manifest["sections"]}
        assert by_id[ids["excited_states"]]["wavefunction_ref"] == ids["wavefunction"]
        for section in by_id.values():
            if section["kind"] in ("wavefunction.gto", "excited_states", "spectra.uvvis"):
                assert section["structure_ref"] == ids["structure"]
        mask = td.get_frozen_mask()
        occupied = np.flatnonzero(mask & (mf.mo_occ == 2))
        virtual = np.flatnonzero(mask & (mf.mo_occ == 0))
        for k, state in enumerate(payload["states"]):
            assert state["index"] == k + 1
            assert state["excitation_energy"] == {"value": td.e[k] * HARTREE_IN_EV, "units": "ev"}
            assert state["total_energy"] == {"value": mf.e_tot + td.e[k], "units": "hartree"}
            assert state["multiplicity"] == (1 if td.singlet else 3)
            np.testing.assert_allclose(state["transition_dipole"], td.transition_dipole()[k], atol=1e-14)
            assert state["oscillator_strength"] == pytest.approx(td.oscillator_strength()[k], abs=1e-14)
            assert state["oscillator_strength"] == pytest.approx(
                2 / 3 * td.e[k] * np.dot(state["transition_dipole"], state["transition_dipole"]), abs=1e-14)
            peak = spectrum["peaks"][k]
            assert peak["assignment"]["state_index"] == k + 1
            assert peak["x"] == spectrum["frequencies"][k] == state["excitation_energy"]["value"]
            assert peak["y"] == spectrum["intensities"][k] == state["oscillator_strength"]
            x, y = td.xy[k]
            expected = x * x - np.asarray(y) ** 2
            expected /= expected.sum()
            contributions = state["contributions"]
            assert [(c["from_mo"], c["to_mo"]) for c in contributions] == [
                (i, a) for i in occupied for a in virtual]
            np.testing.assert_allclose([c["weight"] for c in contributions], expected.ravel(), atol=1e-14)
            assert sum(c["weight"] for c in contributions) == pytest.approx(1)

        assert len(ids["nto_holes"]) == len(td.e) == len(ids["nto_particles"])
        c = reader.tensor(ids["wavefunction"], "mo_coefficients").numpy()
        s = reader.tensor(ids["wavefunction"], "overlap").numpy()
        r = mf.mol.intor_symmetric("int1e_r", comp=3)
        for k, (hole_id, particle_id) in enumerate(zip(ids["nto_holes"], ids["nto_particles"])):
            assert payload["states"][k]["nto_ref"] == hole_id
            hm = reader.json(hole_id, "mo_metadata")
            pm = reader.json(particle_id, "mo_metadata")
            x, y = td.xy[k]
            # The singlet-path matrix: normalized X for TDA, sqrt(2)(X+Y) otherwise.
            singlet_path = x / np.linalg.norm(x) if tda else np.sqrt(2) * (x + y)
            pairs = _nto_channels(reader, hole_id, particle_id)
            assert [ch for ch, *_ in pairs] == (["alpha", "beta"] if not td.singlet else [None])
            total = 0.0
            moments = []
            for channel, h, p, weights, energies in pairs:
                assert len(weights) == min(len(occupied), len(virtual))
                assert np.all(weights[:-1] >= weights[1:])
                assert energies == [0.0] * len(weights)
                total += weights.sum()
                np.testing.assert_allclose(h @ s @ h.T, np.eye(len(weights)), atol=1e-10)
                np.testing.assert_allclose(p @ s @ p.T, np.eye(len(weights)), atol=1e-10)
                np.testing.assert_allclose(h @ s @ p.T, 0, atol=1e-10)
                # Recover the rotations from archive orbitals/overlap, without
                # the adapter's AO transform or its SVD implementation.
                u = c[occupied] @ s @ h.T
                v = c[virtual] @ s @ p.T
                reconstructed = (u * np.sqrt(weights)) @ v.T
                expected = (singlet_path if td.singlet else
                            singlet_path / np.sqrt(2) * (1 if channel == "alpha" else -1))
                np.testing.assert_allclose(reconstructed, expected, atol=1e-10)
                # Contract the archived pairs through source AO integrals.
                # The origin shift cancels because hole and particle spaces
                # are orthogonal.
                hole_ao = mf.mo_coeff[:, occupied] @ u
                particle_ao = mf.mo_coeff[:, virtual] @ v
                moments.append(np.einsum("pk,xpq,qk,k->x", hole_ao, r,
                                         particle_ao, np.sqrt(weights)))
            if tda:
                assert total == pytest.approx(1)
            else:
                assert total == pytest.approx(2 * np.sum((x + y) ** 2))
            for meta, role in [(hm, "hole"), (pm, "particle")]:
                assert meta["orbital_kind"] == "nto" and meta["nto_role"] == role
                assert meta["occupation_semantics"] == "transition_weight"
                assert meta["state_ref"] == ids["excited_states"] and meta["state_index"] == k + 1
                assert meta["spin"] == ("restricted" if td.singlet else "unrestricted")
            if td.singlet:
                # PySCF's singlet moment includes sqrt(2) for the two spins.
                np.testing.assert_allclose(np.sqrt(2) * moments[0], td.transition_dipole()[k], atol=1e-9)
            else:
                # The triplet charge transition (alpha + beta) vanishes, as
                # PySCF's zero triplet dipole says; alpha - beta does not.
                np.testing.assert_allclose(moments[0] + moments[1], 0, atol=1e-12)
                np.testing.assert_allclose(td.transition_dipole()[k], 0, atol=1e-14)


def _nto_channels(reader, hole_id, particle_id):
    """(channel, hole rows, particle rows, weights, energies) per spin channel."""
    hm = reader.json(hole_id, "mo_metadata")
    pm = reader.json(particle_id, "mo_metadata")
    if hm["spin"] == "restricted":
        assert hm["occupations"] == pm["occupations"]
        return [(None, reader.tensor(hole_id, "mo_coefficients").numpy(),
                 reader.tensor(particle_id, "mo_coefficients").numpy(),
                 np.array(hm["occupations"]), hm["energies"])]
    out = []
    for ch in ("alpha", "beta"):
        assert hm[ch]["occupations"] == pm[ch]["occupations"]
        out.append((ch, reader.tensor(hole_id, f"mo_coefficients_{ch}").numpy(),
                    reader.tensor(particle_id, f"mo_coefficients_{ch}").numpy(),
                    np.array(hm[ch]["occupations"]), hm[ch]["energies"]))
    return out


def test_open_writer_ids_and_optional_reference_sections(results, tmp_path):
    path = tmp_path / "custom.qvf"
    with Writer(path, program="example", version="1") as writer:
        ids = add_td(writer, results["rhf-tda"], structure_id="geom", wavefunction_id="reference",
                     states_id="states", spectrum_id="sticks", nto_prefix="pair",
                     populations=False, job_spec=False, density=False, provenance=False)
    with Reader(path) as reader:
        assert reader.validate()["findings"] == []
        assert "provenance" not in reader.raw_manifest
        assert "job_spec" not in ids and "atom_properties" not in ids
        assert ids["nto_holes"] == ["pair_1_hole", "pair_2_hole"]
        assert reader.json("states", "states")["states"][0]["nto_ref"] == "pair_1_hole"


@pytest.mark.parametrize("case", ["rhf-tda", "rhf-full", "casida", "triplet"])
def test_ntos_can_be_omitted(results, tmp_path, case):
    path = tmp_path / "no-ntos.qvf"
    ids = write_td(results[case], path, ntos=False)
    assert ids["nto_holes"] == ids["nto_particles"] == []
    assert validate_qvf(path)["findings"] == []


@pytest.mark.parametrize("kind", ["TDA", "TDHF"])
def test_triplet_ntos_match_unrestricted_response(references, tmp_path, kind):
    # Independent route: PySCF's spin-orbital UHF response of the same closed
    # shell contains the M_S = 0 triplets with their own alpha/beta amplitudes.
    mf = references["minimal"]
    td = getattr(mf, kind)()
    td.singlet = False
    td.nstates = 2
    td.conv_tol = 1e-12
    td.kernel()
    utd = getattr(mf.to_uhf(), kind)()
    utd.nstates = 10
    utd.conv_tol = 1e-12
    utd.kernel()
    path = tmp_path / "triplet.qvf"
    ids = write_td(td, path, ntos=True)
    assert validate_qvf(path)["findings"] == []
    occupied, virtual = mf.mo_occ == 2, mf.mo_occ == 0
    with Reader(path) as reader:
        c = reader.tensor(ids["wavefunction"], "mo_coefficients").numpy()
        s = reader.tensor(ids["wavefunction"], "overlap").numpy()
        for k, energy in enumerate(td.e):
            j = int(np.argmin(abs(utd.e - energy)))
            assert utd.e[j] == pytest.approx(energy, abs=1e-8)
            (xa, xb), (ya, yb) = utd.xy[j]
            ya, yb = (np.zeros_like(xa), np.zeros_like(xb)) if np.ndim(ya) == 0 else (ya, yb)
            channels = {}
            for ch, h, p, weights, _ in _nto_channels(
                    reader, ids["nto_holes"][k], ids["nto_particles"][k]):
                u, v = c[occupied] @ s @ h.T, c[virtual] @ s @ p.T
                channels[ch] = (u * np.sqrt(weights)) @ v.T
            # Eigenvectors converge like sqrt(conv_tol); a sign or sqrt(2)
            # error in the mapping would differ by about 0.1.
            phase = np.sign(np.vdot(channels["alpha"], xa + ya))
            np.testing.assert_allclose(phase * channels["alpha"], xa + ya, atol=1e-5)
            np.testing.assert_allclose(phase * channels["beta"], xb + yb, atol=1e-5)


def _copy_td(td):
    result = copy.copy(td)
    result._scf = copy.copy(td._scf)
    result.e = td.e.copy()
    result.xy = [(x.copy(), np.array(y, copy=True)) for x, y in td.xy]
    result.converged = np.array(td.converged, copy=True)
    return result


@pytest.mark.parametrize("direction", ["parallel", "opposite", "off-diagonal"])
def test_full_response_retains_y_and_nonunit_weight_scale(results, tmp_path, direction):
    # Controlled amplitudes on the actual source orbitals discriminate X-only,
    # X-Y and independently normalized mappings without a new solver claim.
    td = _copy_td(results["minimal-full"])
    x = np.zeros_like(td.xy[0][0])
    y = np.zeros_like(x)
    x[0, 0] = np.sqrt(0.5 + 0.3 ** 2)
    y[1 if direction == "off-diagonal" else 0,
      1 if direction == "off-diagonal" else 0] = -0.3 if direction == "opposite" else 0.3
    td.xy = [(x, y)]
    td.e = td.e[:1]
    td.converged = td.converged[:1]
    path = tmp_path / "controlled.qvf"
    ids = write_td(td, path, ntos=True)
    assert validate_qvf(path)["findings"] == []
    with Reader(path) as reader:
        c = reader.tensor(ids["wavefunction"], "mo_coefficients").numpy()
        s = reader.tensor(ids["wavefunction"], "overlap").numpy()
        h = reader.tensor(ids["nto_holes"][0], "mo_coefficients").numpy()
        p = reader.tensor(ids["nto_particles"][0], "mo_coefficients").numpy()
        weights = np.array(reader.json(ids["nto_holes"][0], "mo_metadata")["occupations"])
        u = c[td._scf.mo_occ == 2] @ s @ h.T
        v = c[td._scf.mo_occ == 0] @ s @ p.T
        transition = (u * np.sqrt(weights)) @ v.T
        np.testing.assert_allclose(transition, np.sqrt(2) * (x + y), atol=1e-10)
        assert weights.sum() == pytest.approx(2 * np.sum((x + y) ** 2))
        assert abs(weights.sum() - 1) > 0.1
        assert not np.allclose(transition, np.sqrt(2) * x)
        assert not np.allclose(transition, np.sqrt(2) * (x - y))
        if direction == "opposite":
            assert weights.sum() < 1
        else:
            assert weights.sum() > 1
    np.testing.assert_array_equal(td.xy[0][0], x)
    np.testing.assert_array_equal(td.xy[0][1], y)


@pytest.mark.parametrize("bad_input", ["unconverged", "scf_unconverged", "empty", "nan_energy",
                                      "negative_energy", "missing_pair", "complex", "bad_shape",
                                      "bad_y", "nonzero_y", "renormalized", "fractional"])
def test_invalid_source_is_rejected_before_writing(results, tmp_path, bad_input):
    td = _copy_td(results["rhf-tda"])
    if bad_input == "unconverged": td.converged[0] = False
    elif bad_input == "scf_unconverged": td._scf.converged = False
    elif bad_input == "empty": td.e = np.array([])
    elif bad_input == "nan_energy": td.e[0] = np.nan
    elif bad_input == "negative_energy": td.e[0] = -1
    elif bad_input == "missing_pair": td.xy.pop()
    elif bad_input == "complex": td.xy[0] = (td.xy[0][0].astype(complex), 0)
    elif bad_input == "bad_shape": td.xy[0] = (td.xy[0][0][:-1], 0)
    elif bad_input == "bad_y": td.xy[0] = (td.xy[0][0], np.zeros((1, 1)))
    elif bad_input == "nonzero_y": td.xy[0] = (td.xy[0][0], td.xy[0][0] * 0.1)
    elif bad_input == "renormalized": td.get_nto(state=1, verbose=0)
    elif bad_input == "fractional":
        td._scf.mo_occ = td._scf.mo_occ.copy()
        td._scf.mo_occ[0] = 1.5
    path = tmp_path / "retained.qvf"
    path.write_bytes(b"original")
    with pytest.raises((ValueError, NotImplementedError)):
        write_td(td, path)
    assert path.read_bytes() == b"original"

