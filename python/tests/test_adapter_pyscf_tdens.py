"""Transition densities on a grid, checked against the archive's own NTOs."""
import pytest

np = pytest.importorskip("numpy")
pyscf = pytest.importorskip("pyscf")
pytest.importorskip("scipy")
from pyscf import gto, lib, scf

from qvf import Reader, validate_qvf
from qvf.adapters.pyscf import write_td
from qvf.units import BOHR_PER_ANGSTROM

from test_adapter_pyscf import _spec_ao_values

GRID = {"spacing": 0.3, "margin": 4.0, "dtype": "float64"}


@pytest.fixture(scope="module")
def results():
    old_threads = lib.num_threads()
    lib.num_threads(1)
    water = gto.M(atom="O 0 0 0; H 0 0.76 0.59; H 0 -0.76 0.59", basis="6-31g",
                  verbose=0, max_memory=256)
    rhf = scf.RHF(water).run(conv_tol=1e-11)
    nh2 = gto.M(atom="N 0 0 0.15; H 0 0.80 -0.40; H 0 -0.80 -0.40", basis="6-31g",
                spin=1, verbose=0, max_memory=256)
    uhf = scf.UHF(nh2).run(conv_tol=1e-11)
    out = {"rhf-tda": rhf.TDA(), "rhf-full": rhf.TDHF(), "uhf-tda": uhf.TDA(),
           "uhf-full": uhf.TDHF(), "triplet": rhf.TDA()}
    out["triplet"].singlet = False
    for td in out.values():
        td.nstates = 2
        td.conv_tol = 1e-9
        td.kernel()
        assert np.asarray(td.converged).all()
    yield out
    lib.num_threads(old_threads)


def _nto_pairs(reader, hole_id, particle_id):
    meta = reader.json(hole_id, "mo_metadata")
    if meta["spin"] == "restricted":
        return math_sqrt2(), [(reader.tensor(hole_id, "mo_coefficients").numpy(),
                               reader.tensor(particle_id, "mo_coefficients").numpy(),
                               np.array(meta["occupations"]))]
    return 1.0, [(reader.tensor(hole_id, f"mo_coefficients_{ch}").numpy(),
                  reader.tensor(particle_id, f"mo_coefficients_{ch}").numpy(),
                  np.array(meta[ch]["occupations"])) for ch in ("alpha", "beta")]


def math_sqrt2():
    return float(np.sqrt(2.0))


@pytest.mark.parametrize("case", ["rhf-tda", "rhf-full", "uhf-tda", "uhf-full"])
def test_transition_densities(results, tmp_path, case):
    td = results[case]
    path = tmp_path / f"{case}.qvf"
    ids = write_td(td, path, transition_densities=GRID)
    report = validate_qvf(path)
    assert report["ok"] and not report["findings"], report
    assert len(ids["transition_densities"]) == len(td.e)
    with Reader(path) as reader:
        sections = {s["id"]: s for s in reader.raw_manifest["sections"]}
        states = reader.json(ids["excited_states"], "states")["states"]
        basis = reader.json(ids["wavefunction"], "basis")
        atoms = reader.json(ids["structure"], "structure")["atoms"]
        centers = np.array([a["position"] for a in atoms]) * BOHR_PER_ANGSTROM
        for k, (state, sid) in enumerate(zip(states, ids["transition_densities"])):
            section = sections[sid]
            assert section["kind"] == "volume.transition_density"
            assert section["excited_states_ref"] == ids["excited_states"]
            assert section["state_index"] == state["index"] == k + 1
            assert state["density_ref"] == sid
            grid = reader.json(sid, "grid")
            rho = reader.tensor(sid, "data").numpy()
            dv = abs(np.linalg.det(np.asarray(grid["voxel_vectors"])))
            ijk = np.stack(np.unravel_index(np.arange(rho.size), grid["shape"]), axis=1)
            points = np.asarray(grid["origin"]) + ijk @ np.asarray(grid["voxel_vectors"])
            values = rho.ravel()
            # Integrates to zero; the first moment is PySCF's transition dipole.
            scale = np.abs(values).sum() * dv
            # Measured: at most 0.55 % of the integral of |rho| and 0.0013 au
            # on this mesh, from the tight O 1s functions. The pointwise
            # check below is the exact one.
            assert abs(values.sum() * dv) <= 1e-2 * scale
            moment = (values[:, None] * points).sum(axis=0) * dv
            dipole = td.transition_dipole()[k]
            assert np.abs(moment - dipole).max() <= 0.01 * np.linalg.norm(dipole) + 1e-3
            # Pointwise: the archive's own NTO pairs through an Appendix A
            # evaluator rebuild the stored values.
            factor, channels = _nto_pairs(reader, state["nto_ref"],
                                          state["nto_ref"].replace("_hole", "_particle"))
            rng = np.random.default_rng(k + 5)
            flat = np.unique(np.concatenate([rng.choice(values.size, 300, replace=False),
                                             np.argsort(np.abs(values))[-20:]]))
            ao = _spec_ao_values(basis, centers, points[flat])
            want = np.zeros(len(flat))
            for h, p, w in channels:
                if len(w):
                    want += np.einsum("pk,pk,k->p", ao @ h.T, ao @ p.T, np.sqrt(w))
            want *= factor
            assert np.allclose(values[flat], want, rtol=1e-8, atol=1e-10 * np.abs(want).max())


def test_default_writes_no_densities(results, tmp_path):
    ids = write_td(results["rhf-tda"], tmp_path / "plain.qvf")
    assert ids["transition_densities"] == []
    with Reader(tmp_path / "plain.qvf") as reader:
        states = reader.json(ids["excited_states"], "states")["states"]
        assert all("density_ref" not in s for s in states)


@pytest.mark.parametrize("option,error", [
    (True, ValueError),  # applied to the triplet result below
    ({"spacing": 0.3, "colour": "red"}, TypeError),
    ("yes", TypeError),
])
def test_refusals_before_opening_target(results, tmp_path, option, error):
    td = results["triplet"] if option is True else results["rhf-tda"]
    path = tmp_path / "retained.qvf"
    path.write_bytes(b"original")
    with pytest.raises(error):
        write_td(td, path, transition_densities=option)
    assert path.read_bytes() == b"original"
