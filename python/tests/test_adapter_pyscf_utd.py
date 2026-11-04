"""Unrestricted (UHF/UKS) response: states, spin contributions and NTO channels."""
import copy

import pytest

np = pytest.importorskip("numpy")
pyscf = pytest.importorskip("pyscf")
from pyscf import dft, gto, lib, scf
from pyscf.tdscf import uhf as tduhf

from qvf import Reader, validate_qvf
from qvf.adapters.pyscf import write_td
from qvf.units import HARTREE_IN_EV

CHANNELS = ("alpha", "beta")


@pytest.fixture(scope="module")
def references():
    old_threads = lib.num_threads()
    lib.num_threads(1)
    # Bent NH2 (2B1) has no zero-frequency rotation root, unlike linear OH (2Pi).
    nh2 = gto.M(atom="N 0 0 0.15; H 0 0.80 -0.40; H 0 -0.80 -0.40",
                basis="6-31g", spin=1, verbose=0, max_memory=256)
    uhf = scf.UHF(nh2).run(conv_tol=1e-11)
    hybrid = dft.UKS(nh2, xc="pbe0")
    hybrid.grids.level = 0
    hybrid.conv_tol = 1e-11
    hybrid.kernel()
    pure = dft.UKS(nh2, xc="pbe")
    pure.grids.level = 0
    pure.conv_tol = 1e-11
    pure.kernel()
    cart = nh2.copy()
    cart.cart = True
    cart.build()
    cartesian = scf.UHF(cart).run(conv_tol=1e-11)
    # One alpha electron, no beta electron: an empty beta occupied space.
    hydrogen = scf.UHF(gto.M(atom="H 0 0 0", basis="6-31++g", spin=1, verbose=0)).run(conv_tol=1e-11)
    yield {"uhf": uhf, "hybrid": hybrid, "pure": pure, "cartesian": cartesian,
           "hydrogen": hydrogen}
    lib.num_threads(old_threads)


@pytest.fixture(scope="module")
def results(references):
    variants = {
        "uhf-tda": references["uhf"].TDA(),
        "uhf-full": references["uhf"].TDHF(),
        "uks-tda": references["hybrid"].TDA(),
        "uks-full": references["hybrid"].TDDFT(),
        "casida": references["pure"].TDDFT(),
        "cartesian": references["cartesian"].TDA(),
        "frozen": references["uhf"].TDHF(),
        "hydrogen": references["hydrogen"].TDA(),
    }
    variants["frozen"].frozen = [0]
    for name, td in variants.items():
        td.nstates = 2 if name != "hydrogen" else 1
        td.conv_tol = 1e-9
        td.kernel()
        assert np.asarray(td.converged).all(), name
    return variants


def _channels(reader, hole_id, particle_id):
    hm = reader.json(hole_id, "mo_metadata")
    pm = reader.json(particle_id, "mo_metadata")
    assert hm["spin"] == pm["spin"] == "unrestricted"
    out = {}
    for ch in CHANNELS:
        assert hm[ch]["occupations"] == pm[ch]["occupations"]
        assert hm[ch]["energies"] == [0.0] * len(hm[ch]["occupations"])
        out[ch] = (reader.tensor(hole_id, f"mo_coefficients_{ch}").numpy(),
                   reader.tensor(particle_id, f"mo_coefficients_{ch}").numpy(),
                   np.array(hm[ch]["occupations"]))
    return out


def _amplitudes(xy):
    (xa, xb), (ya, yb) = xy
    ys = [np.zeros_like(x) if np.ndim(y) == 0 else y for x, y in ((xa, ya), (xb, yb))]
    return [xa, xb], ys


@pytest.mark.parametrize("case", ["uhf-tda", "uhf-full", "uks-tda", "uks-full", "casida",
                                  "cartesian", "frozen", "hydrogen"])
def test_roundtrip_preserves_source_and_channels(results, tmp_path, case):
    td = results[case]
    mf = td._scf
    before = [np.array(a, copy=True) for a in (mf.mo_coeff, mf.mo_occ, mf.mo_energy, td.e)]
    xy_before = copy.deepcopy(td.xy)
    path = tmp_path / f"{case}.qvf"
    ids = write_td(td, path)
    report = validate_qvf(path)
    assert report["ok"] and not report["findings"], report
    for actual, expected in zip((mf.mo_coeff, mf.mo_occ, mf.mo_energy, td.e), before):
        np.testing.assert_array_equal(actual, expected)
    for (x, y), (old_x, old_y) in zip(td.xy, xy_before):
        for a, b in zip(x, old_x):
            np.testing.assert_array_equal(a, b)
        for a, b in zip(y, old_y):
            np.testing.assert_array_equal(a, b)

    tda = isinstance(td, tduhf.TDA) and not isinstance(td, tduhf.TDHF)
    masks = td.get_frozen_mask()
    spaces = {ch: (np.flatnonzero(m & (o == 1)), np.flatnonzero(m & (o == 0)))
              for ch, m, o in zip(CHANNELS, masks, mf.mo_occ)}
    r = mf.mol.intor_symmetric("int1e_r", comp=3)
    with Reader(path) as reader:
        payload = reader.json(ids["excited_states"], "states")
        assert payload["meta"]["tda"] is tda and payload["meta"]["spin_flip"] is False
        assert payload["meta"]["method"] == (
            "tddft" if isinstance(mf, dft.rks.KohnShamDFT) else "cis" if tda else "tdhf")
        c = {ch: reader.tensor(ids["wavefunction"], f"mo_coefficients_{ch}").numpy() for ch in CHANNELS}
        s = reader.tensor(ids["wavefunction"], "overlap").numpy()
        assert len(ids["nto_holes"]) == len(td.e) == len(ids["nto_particles"])
        for k, state in enumerate(payload["states"]):
            assert state["index"] == k + 1
            assert "multiplicity" not in state  # not a spin eigenstate
            assert state["excitation_energy"] == {"value": td.e[k] * HARTREE_IN_EV, "units": "ev"}
            assert state["total_energy"] == {"value": mf.e_tot + td.e[k], "units": "hartree"}
            np.testing.assert_allclose(state["transition_dipole"], td.transition_dipole()[k], atol=1e-14)
            assert state["oscillator_strength"] == pytest.approx(td.oscillator_strength()[k], abs=1e-14)
            xs, ys = _amplitudes(td.xy[k])
            expected = []
            for ch, x, y in zip(CHANNELS, xs, ys):
                occupied, virtual = spaces[ch]
                expected += [(i, a, ch, w) for (i, a), w in
                             zip([(i, a) for i in occupied for a in virtual], (x * x - y * y).ravel())]
            got = [(cn["from_mo"], cn["to_mo"], cn["spin"], cn["weight"]) for cn in state["contributions"]]
            assert [g[:3] for g in got] == [e[:3] for e in expected]
            np.testing.assert_allclose([g[3] for g in got], [e[3] for e in expected], atol=1e-12)
            assert sum(g[3] for g in got) == pytest.approx(1)

            assert state["nto_ref"] == ids["nto_holes"][k]
            pairs = _channels(reader, ids["nto_holes"][k], ids["nto_particles"][k])
            total, moment = 0.0, np.zeros(3)
            for ch, x, y in zip(CHANNELS, xs, ys):
                h, p, weights = pairs[ch]
                occupied, virtual = spaces[ch]
                assert len(weights) == min(len(occupied), len(virtual)) == h.shape[0] == p.shape[0]
                assert np.all(weights >= 0) and np.all(weights[:-1] >= weights[1:])
                total += weights.sum()
                if not len(weights):
                    continue
                np.testing.assert_allclose(h @ s @ h.T, np.eye(len(weights)), atol=1e-10)
                np.testing.assert_allclose(p @ s @ p.T, np.eye(len(weights)), atol=1e-10)
                # Rotations from archive orbitals/overlap only.
                u = c[ch][occupied] @ s @ h.T
                v = c[ch][virtual] @ s @ p.T
                transition = (u * np.sqrt(weights)) @ v.T
                np.testing.assert_allclose(transition, x if tda else x + y, atol=1e-10)
                mo = mf.mo_coeff[CHANNELS.index(ch)]
                moment += np.einsum("pk,xpq,qk,k->x", mo[:, occupied] @ u, r,
                                    mo[:, virtual] @ v, np.sqrt(weights))
            if tda:
                assert total == pytest.approx(1)
            else:
                assert total == pytest.approx(sum(np.sum((x + y) ** 2) for x, y in zip(xs, ys)))
            # The alpha + beta pair contraction is PySCF's own transition dipole.
            np.testing.assert_allclose(moment, td.transition_dipole()[k], atol=1e-9)


def test_empty_beta_channel(results, tmp_path):
    path = tmp_path / "h.qvf"
    ids = write_td(results["hydrogen"], path)
    with Reader(path) as reader:
        pairs = _channels(reader, ids["nto_holes"][0], ids["nto_particles"][0])
        assert len(pairs["beta"][2]) == 0 and len(pairs["alpha"][2]) == 1
        assert pairs["alpha"][2].sum() == pytest.approx(1)


def test_closed_shell_uhf_matches_restricted_singlet(tmp_path):
    # Independent route: a closed-shell UHF TDA contains the restricted
    # singlets; its alpha/beta channels must equal the restricted singlet
    # transition split equally between the spins.
    mol = gto.M(atom="O 0 0 0; H 0 0.76 0.59; H 0 -0.76 0.59", basis="sto-3g", verbose=0)
    rhf = scf.RHF(mol).run(conv_tol=1e-12)
    rtd = rhf.TDA()
    rtd.nstates = 2
    rtd.conv_tol = 1e-10
    rtd.kernel()
    utd = rhf.to_uhf().TDA()
    utd.nstates = 20  # the whole space: partial Davidson runs can skip roots here
    utd.conv_tol = 1e-10
    utd.max_cycle = 200
    utd.kernel()
    assert np.asarray(utd.converged).all()
    restricted = write_td(rtd, tmp_path / "r.qvf")
    unrestricted = write_td(utd, tmp_path / "u.qvf")
    assert validate_qvf(tmp_path / "u.qvf")["findings"] == []
    with Reader(tmp_path / "r.qvf") as rr, Reader(tmp_path / "u.qvf") as ur:
        for k, energy in enumerate(rtd.e):
            j = int(np.argmin(abs(utd.e - energy)))
            assert utd.e[j] == pytest.approx(energy, abs=1e-8)
            rw = np.array(rr.json(restricted["nto_holes"][k], "mo_metadata")["occupations"])
            um = ur.json(unrestricted["nto_holes"][j], "mo_metadata")
            for ch in CHANNELS:
                np.testing.assert_allclose(um[ch]["occupations"], rw / 2, atol=1e-6)


def test_restricted_open_reference_is_refused(tmp_path):
    mol = gto.M(atom="N 0 0 0.15; H 0 0.80 -0.40; H 0 -0.80 -0.40", basis="sto-3g",
                spin=1, verbose=0)
    rohf = scf.ROHF(mol).run(conv_tol=1e-10)
    path = tmp_path / "retained.qvf"
    path.write_bytes(b"original")
    with pytest.raises(NotImplementedError, match="restricted-open"):
        write_td(tduhf.TDA(rohf), path)
    assert path.read_bytes() == b"original"


def _copy_td(td):
    result = copy.copy(td)
    result._scf = copy.copy(td._scf)
    result.e = td.e.copy()
    result.xy = copy.deepcopy(td.xy)
    result.converged = np.array(td.converged, copy=True)
    return result


@pytest.mark.parametrize("bad_input", ["unconverged", "renormalized", "bad_beta_shape",
                                      "nonzero_tda_y", "flat_pair", "complex"])
def test_invalid_source_is_rejected_before_writing(results, tmp_path, bad_input):
    td = _copy_td(results["uhf-tda"])
    (xa, xb), _ = td.xy[0]
    if bad_input == "unconverged":
        td.converged[0] = False
    elif bad_input == "renormalized":
        td.xy[0] = ((xa * 2, xb * 2), (0, 0))
    elif bad_input == "bad_beta_shape":
        td.xy[0] = ((xa, xb[:-1]), (0, 0))
    elif bad_input == "nonzero_tda_y":
        td.xy[0] = ((xa, xb), (xa * 0.1, 0))
    elif bad_input == "flat_pair":
        td.xy[0] = (xa, 0)
    elif bad_input == "complex":
        td.xy[0] = ((xa.astype(complex), xb), (0, 0))
    path = tmp_path / "retained.qvf"
    path.write_bytes(b"original")
    with pytest.raises((ValueError, NotImplementedError)):
        write_td(td, path)
    assert path.read_bytes() == b"original"
