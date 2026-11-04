"""Tests for the standalone QVF reference writer + reader.

These use only pytest + numpy. Two cross-checks run opportunistically when the
optional dependency is importable:

* ``jsonschema`` — full structural schema validation of the manifest.
* ``vibeqc`` — validates the produced archive with vibe-qc's own
  ``validate_qvf`` (the reference producer/consumer), proving the standalone
  writer is contract-compatible with the canonical implementation.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))

from qvf_writer import (  # noqa: E402
    BOHR_PER_ANGSTROM, Lattice, QvfWriter, primitive_norm)
from qvf_reader import QvfReader, validate_qvf  # noqa: E402


def _full_writer() -> QvfWriter:
    """A writer exercising every canonical section kind + root metadata."""
    w = QvfWriter(program="test-code", version="9.9", calculation="everything")
    origin = [-2.0, -2.0, -2.0]
    vox = [[0.2, 0, 0], [0, 0.2, 0], [0, 0, 0.2]]
    dens = np.abs(np.random.default_rng(1).standard_normal((8, 8, 8))).astype(np.float32)

    w.add_structure(
        [{"symbol": "O", "position": [0, 0, 0.117], "atomic_number": 8,
          "atom_name": "O", "residue_name": "HOH", "residue_seq": 1,
          "chain_id": "A", "b_factor": 11.25},
         {"symbol": "H", "position": [0, 0.757, -0.469], "atomic_number": 1},
         {"symbol": "H", "position": [0, -0.757, -0.469], "atomic_number": 1}],
        bonds=[{"i": 0, "j": 1, "order": 1.0}],
        biomolecule={
            "chains": ["A"],
            "residues": [{"name": "HOH", "seq": 1, "chain": "A",
                          "atom_indices": [0, 1, 2]}],
            "secondary_structure": [{"type": "coil", "chain": "A",
                                     "start_seq": 1, "end_seq": 1}],
            "b_factors": [11.25, 14.0, 13.5],
        },
    )
    w.add_bonds([{"i": 0, "j": 2, "order": 1.0}])
    w.add_bond_orders("mayer", [{"i": 0, "j": 1, "order": 0.98, "distance_ang": 0.96}])
    a = w.add_density(dens, origin, vox, label="rho")
    b = w.add_orbital(dens, origin, vox, label="homo", component="real")
    w.add_volume("volume.spin", dens, origin, vox)
    w.add_volume("volume.elf", dens, origin, vox)
    w.add_volume("volume.potential", dens.astype(np.float64), origin, vox, dtype="float64")
    w.add_volume("volume.rdg", dens, origin, vox)
    w.add_volume("volume.generic", dens, origin, vox, label="custom")
    w.add_volume("volume.difference", dens, origin, vox, operand_a=a, operand_b=b,
                 description="rho - homo")
    w.add_basis_ao({"atom_index": 0, "atom_symbol": "O", "shell_index": 0,
                    "primitive_index": 0, "angular_momentum": [0, 0],
                    "shell_type": "s", "exponent": 130.7, "coefficient": 0.15,
                    "is_primitive": False, "is_contracted": True, "ao_index": 0},
                   dens, origin, vox)

    shells = [{"center": 0, "l": 0, "exponents": [1.0, 0.5],
               "coefficients": [0.6, 0.4]},
              {"center": 0, "l": 1, "exponents": [1.0], "coefficients": [1.0]}]
    C = np.random.default_rng(2).standard_normal((4, 4))
    w.add_wavefunction_gto(shells, mo_coefficients=C, energies=[-1, -0.5, 0.5, 1],
                           occupations=[2, 2, 0, 0])

    w.add_atom_properties(mulliken_charge=[-0.7, 0.35, 0.35],
                          spin_population=[0.0, 0.0, 0.0])
    coords = np.random.default_rng(3).standard_normal((5, 3, 3))
    tid = w.add_trajectory({"symbols": ["O", "H", "H"],
                            "energies": [-74.9, -74.95, -74.96, -74.962, -74.963]},
                           coords)
    w.add_reaction_waypoints(tid, {"waypoints": [{"frame_index": 0, "label": "R",
                                                  "kind": "reactant"}]})
    w.add_reaction_path({"symbols": ["O", "H", "H"],
                         "waypoints": [{"frame_index": 0, "label": "R",
                                        "kind": "reactant"}]}, coords)
    w.add_scan_surface({"axis_a_label": "r1", "axis_b_label": "r2"},
                       [0.0, 1.0], [0.0, 1.0], np.zeros((2, 2)))
    disp = np.random.default_rng(4).standard_normal((3, 3, 3))
    w.add_vibrations({"symbols": ["O", "H", "H"], "frequencies": [1595, 3657, 3756]},
                     disp)
    for kind in ("spectra.ir", "spectra.raman", "spectra.uvvis", "spectra.ecd",
                 "spectra.vcd", "spectra.generic"):
        w.add_spectrum(kind, frequencies=[100.0, 200.0], intensities=[1.0, 2.0],
                       x_units="cm-1", y_units="km_per_mol",
                       peaks=[{"x": 100.0, "y": 1.0, "label": "mode 1",
                               "assignment": {"mode_index": 0}}],
                       broadening={"type": "lorentzian", "fwhm": 8.0,
                                   "units": "cm-1"})
    w.add_spectrum("spectra.nmr", payload={"isotope": "1H", "reference": "TMS",
                   "chemical_shifts": [{"atom_index": 1, "symbol": "H",
                                        "isotropic_shift_ppm": 4.6}]})
    w.add_spectrum("spectra.epr", payload={
        "g_tensor": {"principal": [2.0023, 2.0021, 2.0089], "isotropic": 2.0044},
        "hyperfine": [{"atom_index": 0, "symbol": "N", "a_iso_mhz": 45.2}],
        "zero_field_splitting": {"d_mhz": 1200.0, "e_mhz": 30.0}})
    eig = np.random.default_rng(5).standard_normal((1, 10, 4))
    w.add_bands({"n_spin": 1, "n_kpoints": 10, "n_bands": 4,
                 "segments": [{"label_start": "G", "label_end": "X",
                               "k_start": [0, 0, 0], "k_end": [0.5, 0, 0],
                               "n_points": 10}]}, eig)
    w.add_dos_total([-2, -1, 0, 1, 2], [0.1, 0.5, 1.0, 0.5, 0.1], smearing=0.05,
                    fermi_energy_ev=0.0, n_electrons=10.0)
    w.add_dos_projected([-2, -1, 0, 1, 2], np.zeros((2, 5)),
                        channels=[{"atom_index": 0, "symbol": "O", "l": 0,
                                   "label": "O-2s"},
                                  {"atom_index": 0, "symbol": "O", "l": 1,
                                   "label": "O-2p"}])
    w.add_dos_coop_cohp("dos.coop", [-1, 0, 1], np.zeros((1, 3)), np.zeros((1, 3)),
                        {"pairs": [[0, 1]]})
    w.add_dos_coop_cohp("dos.cohp", [-1, 0, 1], np.zeros((1, 3)), np.zeros((1, 3)),
                        {"pairs": [[0, 1]]})
    w.add_structure_symmetry({"space_group": "C2v", "point_group": "C2v"})
    w.add_scf_history([{"iter": 1, "energy_eh": -74.9, "delta_e": 1.0,
                        "diis_error": 0.1}])
    w.add_citations("@article{x2026, title={T}, author={A}, year={2026}}")
    w.add_run_record(program="test-code", program_version="9.9",
                     command="test-code water.inp", exit_status=0,
                     started_utc="2026-07-24T12:00:00Z",
                     finished_utc="2026-07-24T12:00:41Z", sequence=0,
                     input_text="! RHF STO-3G\n* xyz 0 1\nO 0 0 0\n*\n",
                     log_text="SCF CONVERGED\nFINAL ENERGY -74.96\n",
                     files={"input": {"filename": "water.inp"},
                            "log": {"filename": "water.out"},
                            "attachment.geom": {"filename": "water.xyz"}},
                     attachments={"geom": b"1\nwater\nO 0 0 0\n"})
    fs = np.random.default_rng(6).standard_normal((4, 4, 4, 2))
    w.add_fermi_surface({"nk1": 4, "nk2": 4, "nk3": 4, "n_spin": 1,
                         "fermi_energy_ev": 0.0, "band_indices": [3, 4],
                         "lattice_vectors": [[4.2, 0, 0], [0, 4.2, 0], [0, 0, 4.2]]},
                        fs)
    w.add_phonon_bands({"n_atoms": 2, "n_modes": 6, "has_eigenvectors": False,
                        "segments": [{"label_start": "G", "label_end": "X",
                                      "k_start": [0, 0, 0], "k_end": [0.5, 0, 0],
                                      "n_points": 10}]},
                       np.abs(np.random.default_rng(7).standard_normal((10, 6))))
    w.add_phonon_dos({"n_atoms": 2, "n_modes": 6}, [0, 100, 200],
                     [0.1, 1.0, 0.1])
    w.add_equation_of_state([90, 100, 110], [-5.1, -5.2, -5.15],
                            {"model": "birch_murnaghan", "V0": 100.0, "E0": -5.2,
                             "B0": 74.0, "B0_prime": 4.0, "energy_unit": "eV",
                             "volume_unit": "angstrom^3", "pressure_unit": "GPa"})
    w.add_topology_qtaim({"points": [{"type": "bcp", "position": [0.7, 0, 0],
                                      "rho": 0.26, "laplacian": -0.54,
                                      "atom_pair": [0, 1]}]},
                         units={"position": "angstrom"})
    w.set_extensions({"x_orca": {"version": "1.0", "critical": True}})
    w.add_vendor_section("x_orca.epr",
                         json_members={"epr": {"g_tensor": [[2.0, 0, 0], [0, 2.0, 0],
                                                            [0, 0, 2.0]]}},
                         critical=True)
    w.add_forces(np.zeros((3, 3)), meta={"structure_ref": "structure"})
    w.add_hessian(np.eye(9), meta={"atom_masses_amu": [15.999, 1.008, 1.008]})
    w.set_thermochemistry(zpve_eh=0.021, gibbs_free_energy_eh=-74.94,
                          temperature_k=298.15, pressure_atm=1.0)
    w.set_energies(-74.963, method="rhf", components={"zpve": 0.021})
    w.set_constraints(frozen_atoms=[0])
    w.set_viewer_defaults({"auto_open": [a]})
    return w


def test_all_kinds_validate(tmp_path):
    w = _full_writer()
    path = str(tmp_path / "full.qvf")
    w.write(path)
    report = validate_qvf(path)
    assert report["ok"], report["errors"]
    # Every canonical kind + one vendor section present.
    kinds = {s["kind"] for s in w.build_manifest()["sections"]}
    assert "wavefunction.gto" in kinds
    assert "x_orca.epr" in kinds
    assert len(w.build_manifest()["sections"]) >= 35


def test_roundtrip_binary_values(tmp_path):
    w = QvfWriter(program="t", version="1")
    coords = np.arange(30, dtype=np.float64).reshape(5, 2, 3)
    w.add_structure([{"symbol": "H", "position": [0, 0, 0], "atomic_number": 1},
                     {"symbol": "H", "position": [0, 0, 0.74], "atomic_number": 1}])
    w.add_trajectory({"symbols": ["H", "H"]}, coords)
    path = str(tmp_path / "t.qvf")
    w.write(path)
    with QvfReader(path) as r:
        traj = next(s for s in r.sections if s["kind"] == "trajectory")
        got = r.read_array(traj["members"]["coords"])
        np.testing.assert_array_equal(got, coords)


def test_wavefunction_primitive_norm_division(tmp_path):
    # A libint-normalized coefficient of 1.0 on an s primitive (alpha=1) must be
    # divided by primitive_norm(1, 0) on the way out.
    shells = [{"center": 0, "l": 0, "exponents": [1.0], "coefficients": [1.0]}]
    w = QvfWriter(program="t", version="1")
    w.add_structure([{"symbol": "H", "position": [0, 0, 0], "atomic_number": 1}])
    w.add_wavefunction_gto(shells, mo_coefficients=np.eye(1),
                           coeffs_are_libint_normalized=True)
    path = str(tmp_path / "wf.qvf")
    w.write(path)
    with QvfReader(path) as r:
        wf = next(s for s in r.sections if s["kind"] == "wavefunction.gto")
        basis = r.read_json(wf["members"]["basis"])
        expected = 1.0 / primitive_norm(1.0, 0)
        assert abs(basis["shells"][0]["coefficients"][0] - expected) < 1e-12


class TestPeriodicReactionPathLattice:
    """`reaction.path` carries an optional `lattice` member under qvf_version 1.

    Governance ruling 2026-07-10: a periodic reaction path is signalled by the
    presence of `lattice`, never by a manifest version bump.
    """

    @staticmethod
    def _writer():
        w = QvfWriter(program="t", version="1", calculation="neb")
        w.add_structure([{"symbol": "H", "position": [0, 0, 0], "atomic_number": 1}])
        return w, np.zeros((5, 1, 3))

    def test_fixed_cell_lattice_round_trips(self, tmp_path):
        w, coords = self._writer()
        lat = np.eye(3) * 5.0
        w.add_reaction_path({"symbols": ["H"]}, coords, lattice=lat)
        path = str(tmp_path / "rp.qvf")
        w.write(path)
        assert validate_qvf(path)["ok"]
        # Archive stays v1 — the lattice member, not a version bump, marks it periodic.
        with QvfReader(path) as r:
            assert r.manifest["qvf_version"] == 1
            sec = next(s for s in r.sections if s["kind"] == "reaction.path")
            assert "lattice" in sec["members"]
            np.testing.assert_allclose(r.read_array(sec["members"]["lattice"]), lat)

    def test_per_frame_lattice_round_trips(self, tmp_path):
        w, coords = self._writer()
        lat = np.stack([np.eye(3) * (1.0 + i) for i in range(5)])
        w.add_reaction_path({"symbols": ["H"]}, coords, lattice=lat)
        path = str(tmp_path / "rp2.qvf")
        w.write(path)
        assert validate_qvf(path)["ok"]

    def test_molecular_path_has_no_lattice(self, tmp_path):
        w, coords = self._writer()
        w.add_reaction_path({"symbols": ["H"]}, coords)
        path = str(tmp_path / "rp3.qvf")
        w.write(path)
        assert validate_qvf(path)["ok"]
        with QvfReader(path) as r:
            sec = next(s for s in r.sections if s["kind"] == "reaction.path")
            assert "lattice" not in sec["members"]

    def test_bad_lattice_shape_rejected(self):
        w, coords = self._writer()
        with pytest.raises(ValueError):
            w.add_reaction_path({}, coords, lattice=np.zeros((2, 2)))

    def test_per_frame_lattice_frame_count_must_match(self):
        w, coords = self._writer()
        with pytest.raises(ValueError):
            w.add_reaction_path({}, coords, lattice=np.zeros((3, 3, 3)))


class TestLattice:
    """The `Lattice` helper and its two on-disk conventions (spec § 4.3).

    The reference numbers are the "Lattice conventions" worked example: a
    hexagonal cell with a = 2.46 A and a 15.875 A normal, as
    `structure.lattice_vectors` (rows, A) and as the equivalent
    `reaction.path.lattice` (columns, bohr, 4 decimals).
    """

    ROWS_ANGSTROM = [[2.46, 0.0, 0.0], [1.23, 2.1304, 0.0], [0.0, 0.0, 15.875]]
    COLUMNS_BOHR_SPEC = [[4.6487, 2.3244, 0.0],
                         [0.0, 4.0259, 0.0],
                         [0.0, 0.0, 29.9994]]
    SPEC_ATOL = 5e-5   # the spec shows 4 decimals

    def test_constant_is_the_rounded_codata_2018_reciprocal(self):
        assert BOHR_PER_ANGSTROM == 1.8897261246
        assert abs(BOHR_PER_ANGSTROM - 1.0 / 0.529177210903) < 1e-10

    def test_rows_to_columns_matches_spec_example(self):
        cell = Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        np.testing.assert_allclose(cell.columns_bohr(), self.COLUMNS_BOHR_SPEC,
                                   rtol=0, atol=self.SPEC_ATOL)
        assert cell.rows_angstrom() == self.ROWS_ANGSTROM

    def test_columns_are_the_scaled_rows(self):
        cell = Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        rows = np.array(cell.rows_angstrom())
        cols = np.array(cell.columns_bohr())
        for j in range(3):
            np.testing.assert_allclose(cols[:, j], BOHR_PER_ANGSTROM * rows[j],
                                       rtol=0, atol=1e-12)
        # The bohr matrix is NOT its own transpose here: reading it by rows
        # would be a different cell (a = (4.6487, 2.3244, 0)).
        assert not np.allclose(cols, cols.T)

    def test_columns_to_rows_matches_spec_example(self):
        cell = Lattice.from_columns_bohr(self.COLUMNS_BOHR_SPEC)
        np.testing.assert_allclose(
            cell.rows_angstrom(), self.ROWS_ANGSTROM, rtol=0, atol=1e-4)

    def test_round_trip_rows_columns_rows(self):
        cell = Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        back = Lattice.from_columns_bohr(cell.columns_bohr())
        np.testing.assert_allclose(back.rows_angstrom(), self.ROWS_ANGSTROM,
                                   rtol=0, atol=1e-9)

    def test_round_trip_columns_rows_columns(self):
        cols = np.random.default_rng(15).uniform(-9.0, 9.0, (3, 3))
        cell = Lattice.from_columns_bohr(cols)
        np.testing.assert_allclose(cell.columns_bohr(), cols, rtol=0, atol=1e-9)
        np.testing.assert_allclose(
            Lattice.from_rows_angstrom(cell.rows_angstrom()).columns_bohr(),
            cols, rtol=0, atol=1e-9)

    def test_nested_lists_of_floats(self):
        cell = Lattice.from_rows_angstrom(np.array(self.ROWS_ANGSTROM))
        for out in (cell.rows_angstrom(), cell.columns_bohr()):
            assert isinstance(out, list) and len(out) == 3
            for row in out:
                assert isinstance(row, list) and len(row) == 3
                assert all(type(x) is float for x in row)
        # The result is a copy: editing it cannot change the lattice.
        cell.rows_angstrom()[0][0] = 99.0
        assert cell.rows_angstrom() == self.ROWS_ANGSTROM

    def test_volume(self):
        cell = Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        expected = 2.46 * 2.1304 * 15.875   # triangular cell: product of the diagonal
        assert cell.volume_angstrom3() == pytest.approx(expected, rel=1e-12)
        assert cell.volume_bohr3() == pytest.approx(
            expected * BOHR_PER_ANGSTROM ** 3, rel=1e-12)
        assert cell.volume_bohr3() == pytest.approx(
            abs(np.linalg.det(np.array(cell.columns_bohr()))), rel=1e-12)

    def test_volume_is_positive_for_a_left_handed_cell(self):
        left = Lattice.from_rows_angstrom([[0, 2, 0], [1, 0, 0], [0, 0, 3]])
        assert left.volume_angstrom3() == pytest.approx(6.0)

    @pytest.mark.parametrize("bad", [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],                  # 2 x 3
        [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]],                # 3 x 2
        [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],       # flat
        np.zeros((2, 3, 3)),                                 # stack, not one cell
        [[1.0, 0.0, 0.0], [0.0, 1.0], [0.0, 0.0, 1.0]],      # ragged
        [["a", "b", "c"]] * 3,                               # not numbers
        None,
    ])
    def test_shape_errors(self, bad):
        with pytest.raises(ValueError):
            Lattice.from_rows_angstrom(bad)
        with pytest.raises(ValueError):
            Lattice.from_columns_bohr(bad)

    @pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
    def test_non_finite_rejected(self, value):
        bad = [[1.0, 0.0, 0.0], [0.0, value, 0.0], [0.0, 0.0, 1.0]]
        with pytest.raises(ValueError):
            Lattice.from_rows_angstrom(bad)
        with pytest.raises(ValueError):
            Lattice.from_columns_bohr(bad)

    def test_equality_and_repr(self):
        a = Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        assert a == Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        assert a != Lattice.from_rows_angstrom(np.eye(3))
        assert a != self.ROWS_ANGSTROM
        assert "2.46" in repr(a)

    def test_exported(self):
        import qvf_writer
        assert "Lattice" in qvf_writer.__all__
        assert "BOHR_PER_ANGSTROM" in qvf_writer.__all__

    # -- writer integration ------------------------------------------------

    @staticmethod
    def _write(tmp_path, name, build):
        w = QvfWriter(program="t", version="1", calculation="lattice")
        build(w)
        path = str(tmp_path / name)
        w.write(path)
        assert validate_qvf(path)["ok"]
        return path

    @staticmethod
    def _slab(w, lattice_vectors):
        w.add_structure(
            [{"symbol": "C", "position": [0, 0, 0], "atomic_number": 6}],
            pbc=[True, True, False], dimensionality=2,
            lattice_vectors=lattice_vectors)

    @staticmethod
    def _path(w, lattice, n_frames=3):
        w.add_reaction_path({"symbols": ["C"]}, np.zeros((n_frames, 1, 3)),
                            lattice=lattice)

    def test_add_structure_lattice_matches_nested_list(self, tmp_path):
        p_list = self._write(tmp_path, "list.qvf", lambda w: self._slab(
            w, self.ROWS_ANGSTROM))
        p_cell = self._write(tmp_path, "cell.qvf", lambda w: self._slab(
            w, Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)))
        with QvfReader(p_list) as a, QvfReader(p_cell) as b:
            ma = next(s for s in a.sections if s["kind"] == "structure")["members"]
            mb = next(s for s in b.sections if s["kind"] == "structure")["members"]
            # Same payload, byte for byte (the sha256 covers the stored bytes).
            assert ma["structure"]["sha256"] == mb["structure"]["sha256"]
            payload = b.read_json(mb["structure"])
        assert payload["lattice_vectors"] == self.ROWS_ANGSTROM
        assert payload["pbc"] == [True, True, False]

    def test_add_structure_lattice_from_columns_is_written_as_angstrom_rows(
            self, tmp_path):
        cell = Lattice.from_columns_bohr(self.COLUMNS_BOHR_SPEC)
        path = self._write(tmp_path, "cols.qvf", lambda w: self._slab(w, cell))
        with QvfReader(path) as r:
            sec = next(s for s in r.sections if s["kind"] == "structure")
            payload = r.read_json(sec["members"]["structure"])
        np.testing.assert_allclose(payload["lattice_vectors"],
                                   self.ROWS_ANGSTROM, rtol=0, atol=1e-4)

    def test_add_structure_still_requires_a_lattice_for_pbc(self):
        w = QvfWriter(program="t", version="1")
        with pytest.raises(ValueError):
            w.add_structure([{"symbol": "C", "position": [0, 0, 0]}],
                            pbc=[True, True, False])

    def test_add_reaction_path_lattice_writes_bohr_columns(self, tmp_path):
        cell = Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        path = self._write(tmp_path, "rp.qvf", lambda w: self._path(w, cell))
        with QvfReader(path) as r:
            sec = next(s for s in r.sections if s["kind"] == "reaction.path")
            member = sec["members"]["lattice"]
            assert member["dtype"] == "float64" and member["shape"] == [3, 3]
            lat = r.read_array(member)
        # Hand-computed values from the spec example (4 decimals) ...
        np.testing.assert_allclose(lat, self.COLUMNS_BOHR_SPEC,
                                   rtol=0, atol=self.SPEC_ATOL)
        # ... and exactly the transpose of the A rows times the constant.
        np.testing.assert_allclose(
            lat, np.array(self.ROWS_ANGSTROM).T * BOHR_PER_ANGSTROM,
            rtol=0, atol=1e-12)

    def test_add_reaction_path_per_frame_lattices(self, tmp_path):
        cells = [Lattice.from_rows_angstrom(np.array(self.ROWS_ANGSTROM) * s)
                 for s in (1.0, 1.01, 1.02)]
        for seq in (cells, tuple(cells)):
            path = self._write(tmp_path, "rpf.qvf",
                               lambda w, seq=seq: self._path(w, seq))
            with QvfReader(path) as r:
                sec = next(s for s in r.sections if s["kind"] == "reaction.path")
                member = sec["members"]["lattice"]
                assert member["shape"] == [3, 3, 3]
                lat = r.read_array(member)
            for i, cell in enumerate(cells):
                np.testing.assert_allclose(lat[i], cell.columns_bohr(),
                                           rtol=0, atol=1e-12)

    def test_per_frame_lattice_sequence_errors(self):
        cell = Lattice.from_rows_angstrom(self.ROWS_ANGSTROM)
        w = QvfWriter(program="t", version="1")
        with pytest.raises(ValueError):            # 2 cells for 3 frames
            self._path(w, [cell, cell], n_frames=3)
        with pytest.raises(TypeError):             # Lattice mixed with arrays
            self._path(w, [cell, np.eye(3), cell], n_frames=3)
        # A rejected lattice leaves the writer untouched: the same section id
        # is still free.
        self._path(w, cell)

    def test_plain_arrays_are_written_unchanged(self, tmp_path):
        # No Lattice involved: the array is already columns in bohr, as before.
        arr = np.array(self.COLUMNS_BOHR_SPEC)
        path = self._write(tmp_path, "plain.qvf", lambda w: self._path(w, arr))
        with QvfReader(path) as r:
            sec = next(s for s in r.sections if s["kind"] == "reaction.path")
            np.testing.assert_array_equal(
                r.read_array(sec["members"]["lattice"]), arr)


def test_duplicate_id_rejected():
    w = QvfWriter(program="t", version="1")
    w.add_structure([{"symbol": "H", "position": [0, 0, 0], "atomic_number": 1}],
                    section_id="s")
    with pytest.raises(ValueError):
        w.add_bonds([{"i": 0, "j": 0, "order": 1.0}], section_id="s")


def test_bad_dtype_rejected():
    w = QvfWriter(program="t", version="1")
    with pytest.raises(ValueError):
        w.add_volume("volume.density", np.zeros((2, 2, 2)), [0, 0, 0],
                     [[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype="float16")


def test_manifest_schema_conformance(tmp_path):
    jsonschema = pytest.importorskip("jsonschema")
    import json
    schema_path = os.path.join(HERE, "..", "..", "spec", "qvf_manifest.schema.json")
    with open(schema_path, encoding="utf-8") as fh:
        schema = json.load(fh)
    validator = jsonschema.Draft202012Validator(schema)
    manifest = _full_writer().build_manifest()
    errs = sorted(validator.iter_errors(manifest), key=lambda e: list(e.path))
    assert not errs, "\n".join(f"{list(e.path)}: {e.message}" for e in errs[:20])


def test_vibeqc_validate_qvf_accepts_our_archive(tmp_path):
    try:
        import vibeqc.output.formats.qvf as vibeqc_qvf
    except ImportError as exc:
        pytest.skip(f"vibeqc not importable ({exc}); cross-check deferred to M6")
    path = str(tmp_path / "full.qvf")
    _full_writer().write(path)
    report = vibeqc_qvf.validate_qvf(path)
    # vibe-qc's validate_qvf returns a report dict; accept either {"valid": ...}
    # or {"ok": ...} conventions.
    ok = report.get("valid", report.get("ok"))
    assert ok, report


def test_bundled_schema_matches_canonical():
    """Guard: the bundled schema copy is byte-identical to vibe-qc's canonical
    schema when the repo is available (skipped in a standalone checkout)."""
    canonical = os.path.join(HERE, "..", "..", "..", "python", "vibeqc", "output",
                             "formats", "qvf_manifest.schema.json")
    if not os.path.exists(canonical):
        pytest.skip("canonical schema not present (standalone checkout)")
    bundled = os.path.join(HERE, "..", "..", "spec", "qvf_manifest.schema.json")
    with open(canonical, "rb") as a, open(bundled, "rb") as b:
        assert a.read() == b.read(), "bundled schema has drifted from canonical"


class TestBiomoleculeMetadata:
    """Spec § 5.1 biomolecule fields: section peer keys + per-atom carriers."""

    def _structure_section(self, tmp_path, **kwargs):
        w = QvfWriter(program="t", version="0", calculation="bio")
        w.add_structure(
            [{"symbol": "C", "position": [0, 0, 0], "atomic_number": 6,
              "atom_name": "CA", "residue_name": "ALA", "residue_seq": 4,
              "chain_id": "B", "b_factor": 21.5}],
            **kwargs)
        path = os.path.join(str(tmp_path), "bio.qvf")
        w.write(path)
        validate_qvf(path)
        r = QvfReader(path)
        section = next(s for s in r.sections if s["kind"] == "structure")
        return section, r.read_json(section["members"]["structure"])

    def test_all_four_peer_fields_round_trip(self, tmp_path):
        section, _ = self._structure_section(
            tmp_path,
            biomolecule={
                "chains": ["B"],
                "residues": [{"name": "ALA", "seq": 4, "chain": "B",
                              "atom_indices": [0]}],
                "secondary_structure": [{"type": "helix", "chain": "B",
                                         "start_seq": 4, "end_seq": 4}],
                "b_factors": [21.5],
            })
        assert section["chains"] == ["B"]
        assert section["residues"][0]["atom_indices"] == [0]
        assert section["secondary_structure"][0]["type"] == "helix"
        assert section["b_factors"] == [21.5]

    def test_per_atom_fields_pass_through_to_payload(self, tmp_path):
        _, payload = self._structure_section(tmp_path)
        atom = payload["atoms"][0]
        assert atom["atom_name"] == "CA"
        assert atom["residue_name"] == "ALA"
        assert atom["residue_seq"] == 4
        assert atom["chain_id"] == "B"
        assert atom["b_factor"] == 21.5

    def test_plain_structure_carries_none_of_them(self, tmp_path):
        w = QvfWriter(program="t", version="0", calculation="plain")
        w.add_structure([{"symbol": "C", "position": [0, 0, 0],
                          "atomic_number": 6}])
        path = os.path.join(str(tmp_path), "plain.qvf")
        w.write(path)
        validate_qvf(path)
        r = QvfReader(path)
        section = next(s for s in r.sections if s["kind"] == "structure")
        for key in ("chains", "residues", "secondary_structure", "b_factors"):
            assert key not in section
        atom = r.read_json(section["members"]["structure"])["atoms"][0]
        assert set(atom) == {"symbol", "position", "atomic_number"}

    def _residue_pair(self, gly_code):
        """SER 10 with OG alternates A/B, then GLY 10A (spec § 5.1, #18)."""
        atoms = [
            {"symbol": "C", "position": [0, 0, 0], "atom_name": "CA",
             "residue_seq": 10, "chain_id": "A", "hetero": False},
            {"symbol": "O", "position": [1, 0, 0], "atom_name": "OG",
             "residue_seq": 10, "chain_id": "A", "alt_loc": "A",
             "occupancy": 0.6},
            {"symbol": "O", "position": [0, 1, 0], "atom_name": "OG",
             "residue_seq": 10, "chain_id": "A", "alt_loc": "B",
             "occupancy": 0.4},
            {"symbol": "C", "position": [3, 0, 0], "atom_name": "CA",
             "residue_seq": 10, "chain_id": "A", "insertion_code": "A"},
        ]
        gly = {"name": "GLY", "seq": 10, "chain": "A", "atom_indices": [3]}
        if gly_code is not None:
            gly["insertion_code"] = gly_code
        residues = [{"name": "SER", "seq": 10, "chain": "A",
                     "atom_indices": [0, 1, 2]}, gly]
        return atoms, {"residues": residues}

    def test_atom_site_identity_round_trips(self, tmp_path):
        atoms, bio = self._residue_pair("A")
        w = QvfWriter(program="t", version="0", calculation="altloc")
        w.add_structure(atoms, biomolecule=bio)
        path = os.path.join(str(tmp_path), "altloc.qvf")
        w.write(path)
        validate_qvf(path)
        r = QvfReader(path)
        section = next(s for s in r.sections if s["kind"] == "structure")
        payload = r.read_json(section["members"]["structure"])["atoms"]
        assert [a.get("alt_loc") for a in payload] == [None, "A", "B", None]
        assert [a.get("occupancy") for a in payload] == [None, 0.6, 0.4, None]
        assert payload[0]["hetero"] is False
        assert payload[3]["insertion_code"] == "A"
        assert section["residues"][1]["insertion_code"] == "A"

    @pytest.mark.parametrize("gly_code", [None, "B"])
    def test_residue_insertion_code_must_match_atoms(self, gly_code):
        atoms, bio = self._residue_pair(gly_code)
        w = QvfWriter(program="t", version="0", calculation="mismatch")
        with pytest.raises(ValueError, match="insertion_code"):
            w.add_structure(atoms, biomolecule=bio)


class TestRunRecord:
    def test_requires_program_and_payload(self):
        w = QvfWriter(program="t", version="1")
        with pytest.raises(ValueError):
            w.add_run_record(program="", input_text="x")
        with pytest.raises(ValueError):
            w.add_run_record(program="t")

    def test_roundtrip_text_and_metadata(self, tmp_path):
        w = QvfWriter(program="packer", version="0.1", calculation="repack")
        w.add_run_record(program="orca", program_version="6.0.1",
                         command="orca water.inp", exit_status=0,
                         input_text="! RHF STO-3G\n",
                         log_text="FINAL SINGLE POINT ENERGY -74.96\n")
        path = str(tmp_path / "r.qvf")
        w.write(path)
        report = validate_qvf(path)
        assert report["ok"], report["errors"]
        with QvfReader(path) as r:
            s = next(x for x in r.sections if x["kind"] == "run.record")
            # section program identity is independent of root source.program
            assert s["program"] == "orca"
            assert r.manifest["source"]["program"] == "packer"
            assert r.read_member(s["members"]["input"]).decode("utf-8") == \
                "! RHF STO-3G\n"
            assert "ENERGY" in r.read_member(s["members"]["log"]).decode("utf-8")

    def test_validator_flags_non_utf8_log(self, tmp_path):
        w = QvfWriter(program="t", version="1")
        w.add_run_record(program="t", log_text=b"\xff\xfe\x00 not utf-8")
        path = str(tmp_path / "bad.qvf")
        w.write(path)
        report = validate_qvf(path)
        assert not report["ok"]
        assert any("UTF-8" in e for e in report["errors"]), report["errors"]

    def test_validator_flags_dangling_files_role(self, tmp_path):
        w = QvfWriter(program="t", version="1")
        w.add_run_record(program="t", input_text="x",
                         files={"log": {"filename": "a.out"}})
        path = str(tmp_path / "dangling.qvf")
        w.write(path)
        report = validate_qvf(path)
        assert not report["ok"]
        assert any("files index" in e for e in report["errors"]), \
            report["errors"]

    def test_attachments_may_be_arbitrary_bytes(self, tmp_path):
        w = QvfWriter(program="t", version="1")
        w.add_run_record(program="t", input_text="x",
                         attachments={"restart": b"\x00\x01\xfe\xff"})
        path = str(tmp_path / "att.qvf")
        w.write(path)
        report = validate_qvf(path)
        assert report["ok"], report["errors"]


class TestViewerDefaultsHints:
    """Spec § 3.4: `render` and `transfer_function` are per-section hints.

    Recorded in the spec only (governance 2026-10-07): `$defs/ViewerHint` stays
    open, so both keys validate, and a malformed one is never an error, because
    a consumer must not refuse a file over a hint.
    """

    SPIN_HINT = {
        "render": "volume",
        "transfer_function": {
            "range": [-0.02, 0.02],
            "color": [[-0.02, 0.23, 0.30, 0.75], [0.0, 0.95, 0.95, 0.95],
                      [0.02, 0.71, 0.02, 0.15]],
            "opacity": [[-0.02, 0.5], [-0.002, 0.0], [0.002, 0.0], [0.02, 0.5]],
            "colormap": "custom",
        },
    }

    @staticmethod
    def _write(tmp_path, hint):
        w = QvfWriter(program="t", version="1", calculation="spin")
        w.add_structure([{"symbol": "O", "position": [0, 0, 0], "atomic_number": 8}])
        data = np.linspace(-0.03, 0.03, 64, dtype=np.float32).reshape(4, 4, 4)
        sid = w.add_volume("volume.spin", data, [-1.0, -1.0, -1.0],
                           [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
                           section_id="spin")
        w.set_viewer_defaults({sid: hint})
        path = str(tmp_path / "spin.qvf")
        w.write(path)
        return w, path

    def test_spin_example_round_trips_and_validates(self, tmp_path):
        pytest.importorskip("jsonschema")
        w, path = self._write(tmp_path, self.SPIN_HINT)
        expected = {"spin": self.SPIN_HINT}
        assert w.build_manifest()["viewer_defaults"] == expected
        with QvfReader(path) as r:
            assert r.manifest["viewer_defaults"] == expected
        report = validate_qvf(path)
        assert report["ok"], report["errors"]
        assert not any("jsonschema" in m for m in report["warnings"])

    @pytest.mark.parametrize("hint", [
        {"render": "VOLUME "},
        {"render": 3},
        {"transfer_function": "volume"},
        {"transfer_function": {"range": [1, 0]}},
        {"transfer_function": {"range": [0, 1], "color": 5, "opacity": 0.5}},
    ])
    def test_malformed_hint_is_never_an_error(self, tmp_path, hint):
        pytest.importorskip("jsonschema")
        _, path = self._write(tmp_path, hint)
        report = validate_qvf(path)
        assert report["ok"], report["errors"]


class TestProvenanceHostname:
    """Spec § 3.2: ``provenance.hostname`` is deprecated; the writer drops it
    by default and keeps it only on explicit opt-in. Both archives stay valid."""

    @staticmethod
    def _schema_errors(manifest):
        jsonschema = pytest.importorskip("jsonschema")
        import json
        schema_path = os.path.join(HERE, "..", "..", "spec",
                                   "qvf_manifest.schema.json")
        with open(schema_path, encoding="utf-8") as fh:
            schema = json.load(fh)
        validator = jsonschema.Draft202012Validator(schema)
        return [f"{list(e.path)}: {e.message}"
                for e in validator.iter_errors(manifest)]

    @staticmethod
    def _writer():
        w = QvfWriter(program="t", version="1", calculation="prov")
        w.add_structure([{"symbol": "H", "position": [0, 0, 0],
                          "atomic_number": 1}])
        return w

    def test_hostname_dropped_by_default(self, tmp_path):
        w = self._writer()
        w.set_provenance(method="RHF", hostname="example-host")
        manifest = w.build_manifest()
        assert "hostname" not in manifest["provenance"]
        assert manifest["provenance"]["method"] == "RHF"
        path = str(tmp_path / "no_host.qvf")
        w.write(path)
        with QvfReader(path) as r:
            assert "hostname" not in r.manifest["provenance"]
        report = validate_qvf(path)
        assert report["ok"], report["errors"]
        assert not self._schema_errors(manifest)

    def test_hostname_kept_with_allow_hostname(self, tmp_path):
        w = self._writer()
        w.set_provenance(method="RHF", hostname="example-host",
                         allow_hostname=True)
        manifest = w.build_manifest()
        assert manifest["provenance"]["hostname"] == "example-host"
        assert "allow_hostname" not in manifest["provenance"]
        path = str(tmp_path / "host.qvf")
        w.write(path)
        report = validate_qvf(path)
        assert report["ok"], report["errors"]
        assert not self._schema_errors(manifest)

    def test_allow_hostname_is_keyword_only(self):
        with pytest.raises(TypeError):
            self._writer().set_provenance(True)


# --- writer/schema agreement: every optional argument (D02, D03, D08) ---------
import functools  # noqa: E402
import inspect  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402

from qvf_writer import _ATOM_PROPERTY_ROLES, _ATOM_VECTOR_ROLES  # noqa: E402

_ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
_SCHEMA_PATH = os.path.join(_ROOT, "spec", "qvf_manifest.schema.json")
_SCALAR_ROLES = [
    "mulliken_charge", "loewdin_charge", "hirshfeld_charge", "iao_charge",
    "spin_population",
    "cm5_charge", "chelpg_charge", "resp_charge", "esp_charge", "nbo_charge",
    "mbis_charge", "bader_charge", "ddec_charge",
    "valence", "bonded_valence", "fukui_plus", "fukui_minus", "fukui_zero",
    "dual_descriptor", "atomic_energy_eh", "atomic_volume_bohr3",
    "nmr_shielding_iso_ppm"]
_VECTOR_ROLES = ["magnetic_moment", "force_eh_bohr"]
_ROLES = _SCALAR_ROLES + _VECTOR_ROLES


def _schema_validator():
    jsonschema = pytest.importorskip("jsonschema")
    with open(_SCHEMA_PATH, encoding="utf-8") as fh:
        return jsonschema.Draft202012Validator(json.load(fh))


def _schema_errors(manifest):
    return [f"{list(e.path)}: {e.message}"
            for e in _schema_validator().iter_errors(manifest)]


def _assert_valid_archive(writer, path):
    """The archive passes the manifest schema and the full reference validator,
    with jsonschema active (no W_SCHEMA_SKIPPED)."""
    errors = _schema_errors(writer.build_manifest())
    assert not errors, "\n".join(errors[:20])
    writer.write(path)
    report = validate_qvf(path)
    assert report["ok"], report["errors"]
    assert not any(f["code"] == "W_SCHEMA_SKIPPED" for f in report["findings"])


def _axial_self_overlap(l, exponents, coefficients):
    """Self-overlap of the axial component of a contracted shell on N_i-
    normalized primitives (spec Appendix A.1, "Verifying with overlap")."""
    return sum(ca * cb * (2.0 * np.sqrt(a * b) / (a + b)) ** (l + 1.5)
               for a, ca in zip(exponents, coefficients)
               for b, cb in zip(exponents, coefficients))


def _tiny_png(width: int = 2, height: int = 2) -> bytes:
    """A valid grey PNG image, for ``preview`` members."""
    import struct
    import zlib

    def chunk(tag: bytes, data: bytes) -> bytes:
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))
    raw = b"".join(b"\x00" + b"\x80" * width for _ in range(height))
    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


def _every_option_writer() -> QvfWriter:
    """A sibling of ``_full_writer`` that passes every optional argument of
    every ``add_*`` method at least once; ``test_every_add_option_is_exercised``
    keeps it that way."""
    w = QvfWriter(program="test-code", version="9.9", calculation="every option")
    origin = [-1.0, -1.0, -1.0]
    vox = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]
    dens = np.abs(np.random.default_rng(11).standard_normal((4, 4, 4))).astype(np.float32)
    water = [{"symbol": "O", "position": [0, 0, 0.117], "atomic_number": 8},
             {"symbol": "H", "position": [0, 0.757, -0.469], "atomic_number": 1},
             {"symbol": "H", "position": [0, -0.757, -0.469], "atomic_number": 1}]
    cell = Lattice.from_rows_angstrom([[6.0, 0, 0], [0, 6.0, 0], [0, 0, 15.0]])
    sid = w.add_structure(water, pbc=[True, True, False], lattice_vectors=cell,
                          dimensionality=2, bonds=[{"i": 0, "j": 1, "order": 1.0}],
                          biomolecule={"chains": ["A"], "b_factors": [10.0, 12.0, 12.0]},
                          charge=0, multiplicity=1, name="water", smiles="O",
                          inchi="InChI=1S/H2O/h1H2",
                          inchikey="XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                          label="structure", section_id="slab")
    w.add_bonds([{"i": 0, "j": 2, "order": 1.0}], structure_ref=sid, label="bonds",
                section_id="bonds_extra")
    w.add_bond_orders("mayer", [{"i": 0, "j": 1, "order": 0.98, "distance_ang": 0.96}],
                      label="bond_orders", section_id="mayer")
    a = w.add_density(dens, origin, vox, label="rho", section_id="rho")
    b = w.add_orbital(dens, origin, vox, label="homo", component="real",
                      section_id="homo")
    w.add_volume("volume.potential", dens.astype(np.float64), origin, vox,
                 dtype="float64", units={"origin": "bohr", "voxel_vectors": "bohr"},
                 section_id="esp")
    w.add_volume("volume.difference", dens, origin, vox, label="rho - homo",
                 operand_a=a, operand_b=b, description="rho - homo",
                 section_id="diff")
    # A cube-style periodic grid over the slab's two periodic axes (spec § 5.2):
    # voxel_vectors[a] x (shape[a] - 1) is lattice row a in bohr.
    step = 6.0 * BOHR_PER_ANGSTROM / 4
    w.add_volume("volume.density", np.ones((5, 5, 3), np.float32), [0, 0, 0],
                 [[step, 0, 0], [0, step, 0], [0, 0, 0.5]],
                 periodic=[True, True, False], includes_endpoint=True,
                 structure_ref=sid, section_id="rho_periodic")
    w.add_volume("volume.vector", np.zeros((4, 4, 4, 3)), origin, vox,
                 units={"data": "dimensionless"}, section_id="unit_field")
    w.add_volume("volume.complex", dens + 0.5j * dens, origin, vox,
                 section_id="psi_complex")
    w.add_volume("volume.laplacian", dens, origin, vox, section_id="lap")
    w.add_basis_ao({"atom_index": 0, "atom_symbol": "O", "shell_index": 0,
                    "primitive_index": 0, "angular_momentum": [0, 0],
                    "shell_type": "s", "exponent": 130.7, "coefficient": 0.15,
                    "is_primitive": False, "is_contracted": True, "ao_index": 0},
                   dens, origin, vox, label="O 1s", dtype="float64", section_id="ao0")

    shells = [{"center": 0, "l": 0, "exponents": [1.0, 0.5],
               "coefficients": [0.6, 0.4]},
              {"center": 0, "l": 1, "exponents": [1.0], "coefficients": [1.0]}]
    C = np.random.default_rng(12).standard_normal((4, 4))
    # The AOs as written, after the N_i division: one centre, so s and p are
    # orthogonal and the overlap is diagonal (spec Appendix A.1).
    written = [[c / primitive_norm(a, sh["l"])
                for a, c in zip(sh["exponents"], sh["coefficients"])]
               for sh in shells]
    s_s, s_p = (_axial_self_overlap(sh["l"], sh["exponents"], c)
                for sh, c in zip(shells, written))
    w.add_wavefunction_gto(shells, mo_coefficients=C, energies=[-1, -0.5, 0.5, 1],
                           occupations=[2, 2, 0, 0], structure_ref=sid,
                           orbital_kind="canonical", pure=True, k_point=[0, 0, 0],
                           coeffs_are_libint_normalized=True,
                           ecp=[{"center": 0, "n_core_electrons": 2,
                                 "label": "toy-ECP",
                                 "projectors": [
                                     {"l": 1, "exponents": [3.0, 1.2],
                                      "coefficients": [-1.5, -0.4],
                                      "powers": [2, 2], "local": True},
                                     {"l": 0, "exponents": [2.0],
                                      "coefficients": [4.0], "powers": [0]}]}],
                           n_electrons=4, symmetry_labels=["a1", "a1", "b2", "b1"],
                           overlap=np.diag([s_s, s_p, s_p, s_p]),
                           density_matrix=2.0 * C[:2].T @ C[:2],
                           localization_reference_basis="minao",
                           label="wavefunction_gto", section_id="wf_rhf")
    w.add_wavefunction_gto(shells, mo_coefficients_alpha=C, mo_coefficients_beta=C,
                           energies_alpha=[-1, -0.5, 0.5, 1],
                           occupations_alpha=[1, 1, 0, 0],
                           energies_beta=[-1, -0.5, 0.5, 1],
                           occupations_beta=[1, 0, 0, 0],
                           symmetry_labels={"alpha": ["a1", "a1", "b2", "b1"],
                                            "beta": ["a1", "a1", "b2", "b1"]},
                           density_matrix_alpha=C[:2].T @ C[:2],
                           density_matrix_beta=C[:1].T @ C[:1],
                           structure_ref=sid, section_id="wf_uhf")
    # Excited states: an NTO pair, a transition density, the table and the
    # spectrum that names it (spec § 5.2, § 5.3, § 5.4, § 5.12).
    for role, rows in (("hole", C[:2]), ("particle", C[2:])):
        w.add_wavefunction_gto(shells, mo_coefficients=rows, occupations=[0.9, 0.1],
                               structure_ref=sid, orbital_kind="nto",
                               state_ref="tddft", state_index=1, nto_role=role,
                               section_id=f"nto_{role}")
    w.add_volume("volume.transition_density", dens, origin, vox,
                 excited_states_ref="tddft", state_index=1, section_id="tden")
    w.add_excited_states(
        [{"index": 1, "excitation_energy": 7.4, "total_energy": -76.1,
          "multiplicity": 2, "symmetry": "B1", "oscillator_strength": 0.05,
          "rotatory_strength_cgs": 0.0, "transition_dipole": [0.5, 0, 0],
          "transition_magnetic_dipole": [0, 0.2, 0],
          "contributions": [{"from_mo": 1, "to_mo": 2, "spin": "alpha", "weight": 0.6},
                            {"from_mo": 0, "to_mo": 3, "spin": "beta", "weight": 0.3}],
          "nto_ref": "nto_hole", "density_ref": "tden", "x_note": "demo"},
         {"index": 2, "excitation_energy": {"value": 9.2, "units": "ev"},
          "density_ref": "diff"}],
        meta={"method": "tddft", "reference_state": "D0", "n_states": 2,
              "tda": True, "spin_flip": False,
              "units": {"excitation_energy": "ev", "total_energy": "hartree",
                        "oscillator_strength": "dimensionless"}},
        wavefunction_ref="wf_uhf", label="excited_states", section_id="tddft")
    w.add_spectrum("spectra.uvvis", payload={
        "frequencies": [7.4, 9.2], "intensities": [0.05, 0.0],
        "peaks": [{"x": 7.4, "assignment": {"state_index": 1}}],
        "states_ref": "tddft"}, section_id="uvvis")
    G = (np.random.default_rng(20).standard_normal((8, 8))
         + 1j * np.random.default_rng(21).standard_normal((8, 8)))
    w.add_wavefunction_gto(shells, mo_coefficients=G, spin="generalized",
                           complex_coefficients=True,
                           occupations=[1, 1, 1, 1, 0, 0, 0, 0],
                           structure_ref=sid, section_id="wf_ghf")
    Ck = np.random.default_rng(20).standard_normal((2, 2, 4, 4, 2))
    Ck[1, :, 3] = 0.0                     # padding beyond n_mo_per_k at X
    occ_k = np.zeros((2, 2, 4))
    occ_k[:, :, 0] = 1.0
    w.add_wavefunction_bloch(shells, [[0, 0, 0], [0.5, 0, 0]], Ck,
                             np.sort(np.random.default_rng(21).standard_normal((2, 2, 4))),
                             occ_k, weights=[0.5, 0.5], spin="unrestricted",
                             orbital_kind="canonical", structure_ref=sid,
                             coeffs_are_libint_normalized=True,
                             fermi_energy_ev=-4.2, pure=True, labels=["G", "X"],
                             mesh=[2, 1, 1], shift=[0, 0, 0], time_reversal=False,
                             n_mo_per_k=[4, 3], label="wavefunction_bloch", section_id="wf_bloch")

    w.add_atom_properties(**{role: [0.1, -0.2, 0.1] for role in _SCALAR_ROLES},
                          magnetic_moment=[[0, 0, 1.0], [0, 0, -0.5], [0, 0, 0]],
                          force_eh_bohr=[[0, 0, 0.01], [0, 0.02, -0.01],
                                         [0, -0.02, -0.01]],
                          custom={"net_flux": ([0.1, -0.2, 0.1],
                                               {"units": "e", "kind": "scalar",
                                                "description": "Net flux."}),
                                  "dipole_like": np.zeros((3, 3))},
                          meta={"custom.dipole_like": {"units": "debye",
                                                       "description": "Dipole.",
                                                       "kind": "vector"}},
                          label="atom_properties", section_id="charges")
    coords = np.random.default_rng(13).standard_normal((2, 3, 3))
    tid = w.add_trajectory({"symbols": ["O", "H", "H"], "energies": [-74.9, -74.95]},
                           coords, units={"energies": "hartree", "coords": "angstrom"},
                           label="trajectory", section_id="opt")
    w.add_reaction_waypoints(tid, {"waypoints": [{"frame_index": 0, "label": "R",
                                                  "kind": "reactant"}]},
                             label="reaction_waypoints", section_id="opt_marks")
    # A periodic NPT run: one cell per frame in rows, Å (spec § 5.5).
    rng = np.random.default_rng(20)
    md_cells = [Lattice.from_rows_angstrom([[6.0 + 0.1 * k, 0, 0], [0.2, 6.0, 0],
                                            [0, 0, 15.0]]) for k in range(2)]
    w.add_trajectory({"symbols": ["O", "H", "H"], "kind": "md", "ensemble": "npt",
                      "timestep_fs": 0.5, "pbc": [True, True, False],
                      "energies": [-76.30, -76.31]},
                     coords, lattice=md_cells,
                     velocities=rng.standard_normal((2, 3, 3)) * 1e-3,
                     gradient=rng.standard_normal((2, 3, 3)) * 1e-2,
                     times=[0.0, 0.5], energies=[-76.30, -76.31],
                     cell_volumes=[c.volume_angstrom3() for c in md_cells],
                     temperatures=[300.0, 301.5], precision="float32",
                     section_id="md")
    w.add_volume("volume.series", np.stack([dens, 2 * dens]), origin, vox,
                 axis={"kind": "time", "values": [0.0, 0.5], "units": "fs",
                       "labels": ["t = 0", "t = 0.5 fs"], "trajectory_ref": tid},
                 section_id="rho_series")
    # An orbital volume that names its wavefunction, row and channel (spec § 5.2).
    w.add_volume("volume.orbital", dens, origin, vox, wavefunction_ref="wf_uhf",
                 mo_index=1, spin="beta", section_id="mo_beta_1")
    w.add_reaction_path({"symbols": ["O", "H", "H"],
                         "waypoints": [{"frame_index": 0, "label": "R",
                                        "kind": "reactant"}]},
                        coords, lattice=cell,
                        velocities=np.zeros((2, 3, 3)),
                        gradient=rng.standard_normal((2, 3, 3)) * 1e-2,
                        times=[0.0, 1.0], energies=[-76.2, -76.25],
                        label="reaction_path", section_id="neb")
    w.add_scan_surface({"axis_a_label": "r1", "axis_b_label": "r2"},
                       [0.0, 1.0], [0.0, 1.0], np.zeros((2, 2)),
                       geometries=np.zeros((2, 2, 3, 3)),
                       units={"axis_a": "bohr", "axis_b": "radian",
                              "energies": "hartree"}, label="scan_surface", section_id="scan")
    hess = np.random.default_rng(18).standard_normal((9, 9))
    w.add_hessian(hess + hess.T,
                  meta={"structure_ref": sid, "mass_weighted": True,
                        "projected": True, "atom_masses_amu": [15.999, 1.008, 1.008],
                        "units": {"energy": "hartree", "length": "bohr"}},
                  label="hessian", section_id="hess")
    w.add_vibrations({"symbols": ["O", "H", "H"], "frequencies": [1595, 3657, 3756]},
                     np.random.default_rng(14).standard_normal((3, 3, 3)),
                     hessian_ref="hess", label="vibrations", section_id="modes")
    w.add_forces(np.random.default_rng(19).standard_normal((3, 3)),
                 stress=np.diag([-1e-4, -1e-4, -2e-4]),
                 meta={"structure_ref": sid, "method": "rks", "convention": "gradient",
                       "units": {"energy": "hartree", "length": "bohr"}},
                 label="forces", section_id="grad")
    # Property kinds (spec § 5.13); the electric dipole restates the root one.
    w.set_dipole_moment(total_debye=1.0, vector_debye=[0.0, 0.6, 0.8],
                        origin="center_of_mass")
    w.add_electric_properties(
        {"origin": [0.0, 0.0, 0.0],
         "dipole": {"vector": [0.0, 0.6, 0.8], "total": 1.0, "units": "debye"},
         "quadrupole": {"matrix": np.diag([1.0, -0.5, -0.5]), "traceless": True,
                        "convention": "buckingham"},
         "octupole": {"tensor": np.zeros((3, 3, 3)), "convention": "cartesian"},
         "polarizability": {"matrix": np.diag([9.0, 9.0, 9.0]),
                            "frequency_ev": 1.5, "isotropic": 9.0},
         "hyperpolarizability": {"tensor": np.zeros(27),
                                 "frequencies_ev": [1.5, 1.5]},
         "units": {"dipole": "debye", "origin": "bohr"}},
        label="electric_properties", section_id="multipoles")
    stiffness = np.diag([200.0, 200.0, 200.0, 80.0, 80.0, 80.0])
    stiffness[0, 1] = stiffness[1, 0] = 50.0
    w.add_elastic_properties(
        stiffness, compliance=np.linalg.inv(stiffness),
        moduli={"bulk_modulus_gpa": 100.0, "shear_modulus_gpa": 60.0,
                "youngs_modulus_gpa": 150.0, "poisson_ratio": 0.25},
        average="hill",
        units={"stiffness": "gpa", "poisson_ratio": "dimensionless"},
        label="elastic_properties", section_id="elastic_constants")
    w.add_dielectric_properties(
        np.diag([3.0, 3.0, 3.0]), ionic=np.diag([1.0, 1.0, 1.0]),
        static=np.diag([4.0, 4.0, 4.0]), born_charges=np.zeros((3, 3, 3)),
        piezoelectric=np.zeros((3, 6)),
        units={"electronic": "dimensionless", "born_charges": "e"},
        label="dielectric_properties", section_id="permittivity")
    w.set_provenance(method="rks", functional="pbe",
                     scf_energy={"value": -76.38, "units": "hartree"})
    w.set_energies(-76.38, method="rks", units="hartree",
                   components={"nuclear_repulsion": 9.2, "electronic": -85.58,
                               "exchange_correlation": -9.1, "dispersion": -0.001},
                   by_method=[{"method": "rks", "total": -76.38,
                               "components": {"dispersion": -0.001}}],
                   units_note="all energies in hartree")
    w.add_spectrum("spectra.ir", frequencies=[1595.0, 3657.0], intensities=[67.0, 5.0],
                   x_units="cm-1", y_units="km_per_mol",
                   peaks=[{"x": 1595.0, "y": 67.0, "label": "bend",
                           "assignment": {"mode_index": 0}}],
                   broadening={"type": "lorentzian", "fwhm": 8.0, "units": "cm-1"},
                   label="IR", section_id="ir")
    w.add_spectrum("spectra.nmr", payload={"isotope": "1H", "reference": "TMS",
                   "chemical_shifts": [{"atom_index": 1, "symbol": "H",
                                        "isotropic_shift_ppm": 4.6}]},
                   section_id="nmr")
    w.add_bands({"n_spin": 1, "n_kpoints": 4, "n_bands": 2,
                 "segments": [{"label_start": "G", "label_end": "X",
                               "k_start": [0, 0, 0], "k_end": [0.5, 0, 0],
                               "n_points": 4}],
                 "channels": [{"atom": 0, "l": 0, "label": "O-2s"},
                              {"atom": 0, "l": 1, "label": "O-2p"}]},
                np.random.default_rng(15).standard_normal((1, 4, 2)),
                occupations=np.ones((1, 4, 2)),
                projections=np.zeros((4, 2, 2)), label="bands", section_id="fat_bands")
    energies = [-2, -1, 0, 1, 2]
    w.add_dos_total(energies, [0.1, 0.5, 1.0, 0.5, 0.1], smearing=0.05,
                    smearing_type="gaussian", fermi_energy_ev=0.0, n_electrons=10.0,
                    n_spin=1, label="dos_total", section_id="dos")
    w.add_dos_projected(energies, np.zeros((2, 5)),
                        channels=[{"atom": 0, "l": 0, "m": 0, "spin": 0,
                                   "label": "O-2s"},
                                  {"atom": 0, "l": 1, "m": -1, "spin": 0,
                                   "label": "O-2p"}],
                        n_spin=1, fermi_energy_ev=0.0, label="dos_projected", section_id="pdos")
    w.add_dos_coop_cohp("dos.coop", [-1, 0, 1], np.zeros((1, 3)), np.zeros((1, 3)),
                        {"pairs": [[0, 1]]}, label="dos_coop_cohp", section_id="coop")
    w.add_structure_symmetry({"space_group": "C2v", "point_group": "C2v"},
                             label="structure_symmetry", section_id="sym")
    w.add_scf_history([{"iter": 1, "energy_eh": -74.9, "delta_e": 1.0,
                        "diis_error": 0.1}], label="scf_history", section_id="scf")
    w.add_citations("@article{x2026, title={T}, author={A}, year={2026}}",
                    label="citations", section_id="refs")
    w.add_run_record(program="test-code", program_version="9.9",
                     command="test-code water.inp", exit_status=0,
                     started_utc="2026-07-24T12:00:00Z",
                     finished_utc="2026-07-24T12:00:41Z", sequence=0,
                     input_text="! RHF STO-3G\n", log_text="SCF CONVERGED\n",
                     files={"input": {"filename": "water.inp"},
                            "attachment.geom": {"filename": "water.xyz"}},
                     attachments={"geom": b"3\nwater\n"}, label="run_record", section_id="run")
    w.add_job_spec(job_type="periodic", method="rks", basis="pob-tzvp-rev2",
                   functional="pbe", charge=0, multiplicity=1, kpoints=[4, 4, 1],
                   tasks=["single_point"], options={"cutoff_ha": 300.0},
                   structure_ref=sid,
                   solvation={"model": "cpcm", "solvent": "water", "epsilon": 78.4},
                   dispersion="d3bj", relativistic="x2c",
                   scf={"convergence": 1e-8, "max_iterations": 100, "guess": "sad",
                        "damping": 0.3},
                   geometry={"convergence": 4.5e-4, "max_steps": 50,
                             "coordinates": "redundant"},
                   resources={"n_cores": 16, "memory_mb": 32000,
                              "walltime_s": 86400},
                   label="job_spec", section_id="job")
    # The run of that job (a job set's bookkeeping, spec § 5.8, § 5.9).
    w.add_run_record(program="test-code", log_text="done\n", job_ref="job",
                     run_status="converged", section_id="run_job")
    w.add_reaction_scheme(
        [{"id": "r1", "label": "a", "reactants": [sid],
          "reagents": [{"structure_ref": sid, "coefficient": 2, "x_role": "catalyst"}],
          "products": [{"structure_ref": sid, "coefficient": 1}],
          "arrow": "equilibrium", "conditions": ["THF, -78 °C", "2 h"],
          "yield_percent": 87, "yield_text": "87%", "x_note": "every option"}],
        label="scheme 1", section_id="scheme")
    w.add_preview(_tiny_png(), thumbnail_media_type="image/png",
                  image="<svg xmlns='http://www.w3.org/2000/svg' width='4' height='4'/>",
                  image_media_type="image/svg+xml", label="picture",
                  section_id="picture")
    # The typed provenance fields (spec § 3.2); repeated calls merge.
    w.set_provenance(software=[{"name": "test-code", "version": "9.9"},
                               {"name": "integral-lib", "version": "2.9.0",
                                "url": "https://example.org/integral-lib"}],
                     git_commit="9f2b6c1d4e7a", n_threads=8, n_mpi_ranks=4,
                     n_gpus=0, started_utc="2026-07-24T12:00:00Z",
                     finished_utc="2026-07-24T12:00:41Z")
    w.add_fermi_surface({"nk1": 2, "nk2": 2, "nk3": 2, "n_spin": 1,
                         "fermi_energy_ev": 0.0, "band_indices": [3, 4]},
                        np.zeros((2, 2, 2, 2)), label="fermi_surface", section_id="fs")
    eig = np.random.default_rng(16).standard_normal((4, 9, 3, 3))
    w.add_phonon_bands({"n_atoms": 3, "n_modes": 9, "has_eigenvectors": True,
                        "segments": [{"label_start": "G", "label_end": "X",
                                      "k_start": [0, 0, 0], "k_end": [0.5, 0, 0],
                                      "n_points": 4}]},
                       np.abs(np.random.default_rng(17).standard_normal((4, 9))),
                       eigenvectors=eig + 0.5j * eig, label="phonon_bands", section_id="phonons")
    w.add_phonon_dos({"n_atoms": 3, "n_modes": 9,
                      "channels": [{"atom": 0, "label": "O"},
                                   {"atom": 1, "label": "H"}]},
                     [0, 100, 200], [0.1, 1.0, 0.1],
                     projected=[[0.05, 0.6, 0.05], [0.05, 0.4, 0.05]],
                     label="phonon_dos", section_id="phonon_dos_proj")
    w.add_equation_of_state([90, 100, 110], [-5.1, -5.2, -5.15],
                            {"model": "birch_murnaghan", "V0": 100.0, "E0": -5.2,
                             "B0": 74.0, "B0_prime": 4.0}, label="equation_of_state", section_id="eos_fit")
    w.add_topology_qtaim({"points": [{"type": "bcp", "position": [0.7, 0, 0],
                                      "rho": 0.26, "laplacian": -0.54}]},
                         units={"position": "angstrom"}, label="topology_qtaim", section_id="aim")
    w.set_extensions({"x_demo": {"version": "1.0", "critical": True}})
    w.add_vendor_section("x_demo.custom", json_members={"payload": {"answer": 42}},
                         binary_members={"table": (np.arange(3), "int32")},
                         opaque_members={"drawing": (b"<svg/>", "image/svg+xml"),
                                         "blob": b"\x00\x01"},
                         critical=True, schema_uri="https://example.org/x_demo.json",
                         label="demo", section_id="demo")
    # Section references and the root metadata block (spec § 4.1, § 3.5).
    w.annotate_section("rho", structure_ref=sid, run_ref="run", group="densities",
                       label="electron density")
    w.annotate_section("wf_rhf", structure_ref=sid, run_ref="run", group="orbitals")
    w.annotate_section("demo", structure_ref=sid, group="vendor")
    w.set_metadata(title="every option", description="writer coverage",
                   modified_utc="2026-10-07T12:00:00Z",
                   authors=[{"name": "A. Researcher", "orcid": "0000-0002-1825-0097"}],
                   license="Apache-2.0", doi="10.5281/zenodo.1234567",
                   keywords=["coverage"], derived_from=["10.5281/zenodo.1234567"],
                   tags=["test"])
    return w


def _optional_parameters(fn):
    return {name for name, p in inspect.signature(fn).parameters.items()
            if p.default is not inspect.Parameter.empty}


def test_every_add_option_is_exercised(monkeypatch):
    """``_every_option_writer`` passes every optional argument of every
    ``add_*`` method, so a writer option the schema does not accept cannot
    hide (the D03 failure mode)."""
    seen: dict[str, set[str]] = {}
    methods = {name: fn for name, fn in vars(QvfWriter).items()
               if name.startswith("add_") and inspect.isfunction(fn)}

    def recording(name, fn):
        sig = inspect.signature(fn)

        @functools.wraps(fn)
        def wrapper(self, *args, **kwargs):
            seen.setdefault(name, set()).update(
                sig.bind(self, *args, **kwargs).arguments)
            return fn(self, *args, **kwargs)
        return wrapper

    for name, fn in methods.items():
        monkeypatch.setattr(QvfWriter, name, recording(name, fn))
    _every_option_writer()
    assert set(seen) == set(methods), f"add_* never called: {set(methods) - set(seen)}"
    missing = {name: sorted(_optional_parameters(fn) - seen[name])
               for name, fn in methods.items()}
    assert not {k: v for k, v in missing.items() if v}, missing


def test_every_option_writer_validates(tmp_path):
    """Regression for D02/D03/D08: an archive using every optional writer
    argument validates against the schema, with jsonschema active."""
    _assert_valid_archive(_every_option_writer(), str(tmp_path / "every.qvf"))


class TestPhononOptionalMembers:
    """Spec § 5.6: ``phonon_bands.eigenvectors`` and ``phonon_dos.projected``."""

    QPATH = {"n_atoms": 1, "n_modes": 3, "has_eigenvectors": True,
             "segments": [{"label_start": "G", "label_end": "X",
                           "k_start": [0, 0, 0], "k_end": [0.5, 0, 0],
                           "n_points": 4}]}

    @staticmethod
    def _writer():
        w = QvfWriter(program="t", version="1", calculation="phonons")
        w.add_structure([{"symbol": "Si", "position": [0, 0, 0], "atomic_number": 14}],
                        pbc=[True, True, True],
                        lattice_vectors=[[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])
        return w

    def test_d03_reproduction_now_validates(self, tmp_path):
        """docs/format_review.md D03: both optional members used to fail the
        schema; the archive must now validate."""
        w = self._writer()
        w.add_phonon_bands(self.QPATH, np.zeros((4, 3)),
                           eigenvectors=np.zeros((4, 3, 1, 3)))
        w.add_phonon_dos({"channels": [{"atom": 0, "label": "A"}]}, [0, 1, 2], [0, 1, 0],
                         projected=np.zeros((1, 3)))
        _assert_valid_archive(w, str(tmp_path / "d03.qvf"))

    def test_eigenvector_forms_round_trip(self, tmp_path):
        rng = np.random.default_rng(21)
        real = rng.standard_normal((4, 3, 1, 3))
        cplx = real + 1j * rng.standard_normal((4, 3, 1, 3))
        stacked = np.stack([cplx.real, cplx.imag], axis=-1)
        w = self._writer()
        w.add_phonon_bands(self.QPATH, np.ones((4, 3)), eigenvectors=real,
                           section_id="real")
        w.add_phonon_bands(self.QPATH, np.ones((4, 3)), eigenvectors=cplx,
                           section_id="complex")
        w.add_phonon_bands(self.QPATH, np.ones((4, 3)), eigenvectors=stacked,
                           section_id="stacked")
        path = str(tmp_path / "eig.qvf")
        _assert_valid_archive(w, path)
        with QvfReader(path) as r:
            by_id = {s["id"]: s for s in r.sections}
            got_real = r.read_array(by_id["real"]["members"]["eigenvectors"])
            got_cplx = r.read_array(by_id["complex"]["members"]["eigenvectors"])
            got_stacked = r.read_array(by_id["stacked"]["members"]["eigenvectors"])
        assert by_id["real"]["members"]["eigenvectors"]["shape"] == [4, 3, 1, 3]
        assert by_id["complex"]["members"]["eigenvectors"]["shape"] == [4, 3, 1, 3, 2]
        assert by_id["complex"]["members"]["eigenvectors"]["dtype"] == "float64"
        np.testing.assert_array_equal(got_real, real)
        np.testing.assert_array_equal(got_cplx, stacked)
        np.testing.assert_array_equal(got_stacked, stacked)

    @pytest.mark.parametrize("eigenvectors, frequencies", [
        (np.zeros((4, 3, 1)), np.zeros((4, 3))),             # rank 3
        (np.zeros((4, 3, 1, 2)), np.zeros((4, 3))),          # trailing axis 2
        (np.zeros((4, 3, 1, 3, 3)), np.zeros((4, 3))),       # complex axis 3
        (np.zeros((4, 3, 1, 3, 2, 1)), np.zeros((4, 3))),    # rank 6
        (np.zeros((5, 3, 1, 3)), np.zeros((4, 3))),          # n_q mismatch
        (np.zeros((4, 2, 1, 3)), np.zeros((4, 3))),          # n_modes mismatch
        (np.zeros((4, 3, 1, 3)), np.zeros(12)),              # 1-D frequencies
        (np.full((4, 3, 1, 3), "x"), np.zeros((4, 3))),      # not numeric
    ])
    def test_bad_eigenvectors_rejected_before_writing(self, eigenvectors, frequencies):
        w = self._writer()
        before = json.loads(json.dumps(w.build_manifest()))
        with pytest.raises(ValueError, match="phonon_bands"):
            w.add_phonon_bands(self.QPATH, frequencies, eigenvectors=eigenvectors)
        assert w.build_manifest() == before
        # Nothing was registered: the default id and paths are still free.
        w.add_phonon_bands(self.QPATH, np.zeros((4, 3)))

    def test_projected_dos_round_trips(self, tmp_path):
        projected = np.array([[0.0, 0.5, 1.0], [1.0, 0.5, 0.0]])
        meta = {"channels": [{"atom": 0, "label": "Si"}, {"label": "rest"}]}
        w = self._writer()
        w.add_phonon_dos(meta, [0, 100, 200], [0, 1, 0], projected=projected)
        path = str(tmp_path / "pdos.qvf")
        _assert_valid_archive(w, path)
        with QvfReader(path) as r:
            sec = next(s for s in r.sections if s["kind"] == "phonon_dos")
            np.testing.assert_array_equal(
                r.read_array(sec["members"]["projected"]), projected)
            assert r.read_json(sec["members"]["meta"]) == meta
        assert sec["members"]["projected"]["shape"] == [2, 3]

    @pytest.mark.parametrize("projected, meta", [
        (np.zeros(3), {}),                                    # rank 1
        (np.zeros((2, 3, 1)), {}),                            # rank 3
        (np.zeros((2, 4)), {}),                               # n_points mismatch
        (np.zeros((2, 3)), {"channels": [{"atom": 0}]}),      # channel count
        (np.zeros((2, 3)) + 1j, {}),                          # complex
        (np.zeros((2, 3), dtype=bool), {}),                   # not a number
    ])
    def test_bad_projected_rejected_before_writing(self, projected, meta):
        w = self._writer()
        before = json.loads(json.dumps(w.build_manifest()))
        with pytest.raises(ValueError, match="phonon_dos"):
            w.add_phonon_dos(meta, [0, 100, 200], [0, 1, 0], projected=projected)
        assert w.build_manifest() == before
        w.add_phonon_dos(meta, [0, 100, 200], [0, 1, 0])


class TestAtomPropertyRoles:
    """D02, P10: one role list in the schema, spec § 5 / § 5.7, both writers and
    the C ABI."""

    def test_one_role_list_everywhere(self):
        with open(_SCHEMA_PATH, encoding="utf-8") as fh:
            schema = json.load(fh)
        defs = schema["$defs"]
        members = defs["SectionAtomProperties"]["properties"]["members"]
        # The roles, then the `meta` member that describes custom roles.
        assert list(members["properties"]) == _ROLES + ["meta"]
        for role in _SCALAR_ROLES:
            assert members["properties"][role]["$ref"] == "#/$defs/AtomChargeBinary"
        for role in _VECTOR_ROLES:
            assert members["properties"][role]["$ref"] == "#/$defs/AtomVectorBinary"
        assert list(_ATOM_PROPERTY_ROLES) == _ROLES
        assert list(_ATOM_VECTOR_ROLES) == _VECTOR_ROLES
        params = list(inspect.signature(QvfWriter.add_atom_properties).parameters)
        assert params == ["self", *_ROLES, "custom", "meta", "label", "section_id"]

        sources = {
            "spec": os.path.join(_ROOT, "spec", "qvf-format-spec.md"),
            "cpp": os.path.join(_ROOT, "cpp", "src", "qvf.cpp"),
            "c_abi": os.path.join(_ROOT, "cpp", "include", "qvf", "qvf_c.h"),
            "guide": os.path.join(_ROOT, "docs", "library_guide.md"),
        }
        if not all(os.path.exists(p) for p in sources.values()):
            pytest.skip("spec / C++ sources not present (standalone checkout)")
        text = {}
        for key, path in sources.items():
            with open(path, encoding="utf-8") as fh:
                text[key] = fh.read()

        def ticked(s):
            return re.findall(r"`([a-z_0-9]+)`", s)

        row = next(line for line in text["spec"].splitlines()
                   if line.startswith("| `atom_properties` |"))
        assert ticked(row.split("|")[2]) == _ROLES + ["meta"], "spec § 5 table row"
        # The § 5.7 role table: header "| Role | Shape |", then one row per role.
        lines = text["spec"].splitlines()
        head = next(i for i, line in enumerate(lines)
                    if line.startswith("| Role | Shape |"))
        table = []
        for line in lines[head + 2:]:
            if not line.startswith("|"):
                break
            table.append(line.split("|")[1:3])
        assert [ticked(c[0])[0] for c in table] == _ROLES, "spec § 5.7 table"
        for cells in table:
            role = ticked(cells[0])[0]
            want = "`[n_atoms, 3]`" if role in _VECTOR_ROLES else "`[n_atoms]`"
            assert cells[1].strip() == want, f"spec § 5.7 shape of {role}"
        cpp = re.search(r"kAtomPropertyRoles\[\]\s*=\s*\{([^}]*)\}", text["cpp"])
        assert cpp and re.findall(r'"([a-z_0-9]+)"', cpp.group(1)) == _ROLES, \
            "C++ writer"
        cpp_vec = re.search(r"kAtomVectorRoles\[\]\s*=\s*\{([^}]*)\}", text["cpp"])
        assert cpp_vec and re.findall(r'"([a-z_0-9]+)"', cpp_vec.group(1)) == \
            _VECTOR_ROLES, "C++ writer vector roles"
        assert all(f'"{role}"' in text["c_abi"] for role in _ROLES), "C ABI docs"
        cheat = next(line for line in text["guide"].splitlines()
                     if "| `add_atom_properties`" in line)
        assert ticked(cheat.split("|")[1]) == _ROLES + ["meta"], \
            "library guide cheat sheet"

    def test_all_roles_round_trip(self, tmp_path):
        values = {role: [-0.1 * (i + 1), 0.1 * (i + 1)]
                  for i, role in enumerate(_SCALAR_ROLES)}
        values["magnetic_moment"] = [[0.0, 0.0, 1.1], [0.0, 0.0, -0.1]]
        values["force_eh_bohr"] = [[0.0, 0.0, 0.02], [0.0, 0.0, -0.02]]
        w = QvfWriter(program="t", version="1", calculation="charges")
        w.add_structure([{"symbol": "Li", "position": [0, 0, 0], "atomic_number": 3},
                         {"symbol": "H", "position": [0, 0, 1.6], "atomic_number": 1}])
        w.add_atom_properties(**values)
        path = str(tmp_path / "charges.qvf")
        _assert_valid_archive(w, path)
        with QvfReader(path) as r:
            sec = next(s for s in r.sections if s["kind"] == "atom_properties")
            assert list(sec["members"]) == _ROLES
            for role in _ROLES:
                member = sec["members"][role]
                assert member["dtype"] == "float64"
                assert member["shape"] == ([2, 3] if role in _VECTOR_ROLES else [2])
                np.testing.assert_array_equal(r.read_array(member), values[role])

    @pytest.mark.parametrize("kwargs, match", [
        ({}, "at least one"),
        ({"mulliken_charge": [0.1, -0.1], "iao_charge": [0.1, -0.1, 0.0]},
         "differ in length"),
        ({"mulliken_charge": [0.1, -0.1],
          "magnetic_moment": [[0, 0, 1.0], [0, 0, 0], [0, 0, 0]]},
         "differ in length"),
        ({"hirshfeld_charge": [[0.1, -0.1]]}, r"\[n_atoms\]"),
        ({"iao_charge": []}, "no values"),
        ({"loewdin_charge": [0.1 + 1j, -0.1]}, "real numbers"),
        ({"magnetic_moment": [0.0, 1.0]}, r"\[n_atoms, 3\]"),
        ({"force_eh_bohr": [[0.0, 1.0], [0.0, 1.0]]}, r"\[n_atoms, 3\]"),
        ({"force_eh_bohr": np.zeros((2, 3, 1))}, r"\[n_atoms, 3\]"),
        ({"magnetic_moment": np.zeros((0, 3))}, "no values"),
        ({"cm5_charge": [[0.1, 0.2, 0.3]]}, r"\[n_atoms\]"),
    ])
    def test_bad_roles_rejected_before_writing(self, kwargs, match):
        w = QvfWriter(program="t", version="1")
        with pytest.raises(ValueError, match=match):
            w.add_atom_properties(**kwargs)
        assert w.build_manifest()["sections"] == []
        w.add_atom_properties(mulliken_charge=[0.1, -0.1], iao_charge=[0.2, -0.2])

    def test_unknown_role_is_not_a_keyword(self):
        with pytest.raises(TypeError):
            QvfWriter(program="t", version="1").add_atom_properties(
                bogus_charge=[0.1, -0.1])


class TestAtomPropertyCustomRoles:
    """P10: ``custom.<name>`` roles with a required ``meta`` member."""

    ENTRY = {"units": "e", "description": "Born effective charge, mean."}

    @staticmethod
    def _writer():
        w = QvfWriter(program="t", version="1", calculation="custom")
        w.add_structure([{"symbol": "Li", "position": [0, 0, 0], "atomic_number": 3},
                         {"symbol": "H", "position": [0, 0, 1.6], "atomic_number": 1}])
        return w

    def test_custom_roles_round_trip(self, tmp_path):
        w = self._writer()
        w.add_atom_properties(
            mulliken_charge=[0.3, -0.3],
            custom={"born_mean": ([0.9, -0.9], self.ENTRY),
                    "custom.polar_axis": ([[0, 0, 1.0], [0, 0, -1.0]],
                                          {"units": "debye", "description": "d",
                                           "kind": "vector", "x_note": 1})})
        path = str(tmp_path / "custom.qvf")
        _assert_valid_archive(w, path)
        with QvfReader(path) as r:
            sec = next(s for s in r.sections if s["kind"] == "atom_properties")
            assert list(sec["members"]) == [
                "mulliken_charge", "custom.born_mean", "custom.polar_axis", "meta"]
            assert sec["members"]["custom.polar_axis"]["shape"] == [2, 3]
            meta = r.read_json(sec["members"]["meta"])
            assert meta == {
                "custom.born_mean": {**self.ENTRY, "kind": "scalar"},
                "custom.polar_axis": {"units": "debye", "description": "d",
                                      "kind": "vector", "x_note": 1}}
            np.testing.assert_array_equal(
                r.read_array(sec["members"]["custom.born_mean"]), [0.9, -0.9])

    def test_meta_keyword_supplies_the_entries(self, tmp_path):
        w = self._writer()
        w.add_atom_properties(custom={"born_mean": [0.9, -0.9]},
                              meta={"born_mean": self.ENTRY})
        _assert_valid_archive(w, str(tmp_path / "meta.qvf"))

    def test_custom_roles_alone_make_a_section(self, tmp_path):
        w = self._writer()
        w.add_atom_properties(custom={"x": ([1.0, 2.0], self.ENTRY)})
        _assert_valid_archive(w, str(tmp_path / "alone.qvf"))

    @pytest.mark.parametrize("kwargs, match", [
        ({"custom": {"a": [1.0, 2.0]}}, "needs a meta entry"),
        ({"custom": {"a": ([1.0, 2.0], {"units": "e"})}}, "description"),
        ({"custom": {"a": ([1.0, 2.0], {"description": "d"})}}, "units"),
        ({"custom": {"a": ([1.0, 2.0], {"units": "", "description": "d"})}},
         "units"),
        ({"custom": {"a": ([1.0, 2.0], {"units": "e", "description": "d",
                                        "extra": 1})}}, "unknown key"),
        ({"custom": {"a": ([1.0, 2.0], {"units": "e", "description": "d",
                                        "kind": "vector"})}}, "disagrees"),
        ({"custom": {"a": ([[0, 0, 1.0], [0, 0, 0]], {"units": "e",
                                                      "description": "d",
                                                      "kind": "scalar"})}},
         "disagrees"),
        ({"custom": {"a b": ([1.0, 2.0], {"units": "e", "description": "d"})}},
         "must match"),
        ({"custom": {"a.b": ([1.0, 2.0], {"units": "e", "description": "d"})}},
         "must match"),
        ({"custom": {"": ([1.0, 2.0], {"units": "e", "description": "d"})}},
         "must match"),
        ({"custom": {"a": ([1.0, 2.0], {"units": "e", "description": "d"})},
          "meta": {"a": {"units": "e", "description": "d"}}}, "both"),
        ({"custom": {"a": ([1.0, 2.0], {"units": "e", "description": "d"})},
          "meta": {"b": {"units": "e", "description": "d"}}}, "does not carry"),
        ({"custom": {"a": ([1.0, 2.0, 3.0], {"units": "e", "description": "d"})},
          "mulliken_charge": [0.1, 0.2]}, "differ in length"),
        ({"custom": {"a": ([[1.0, 2.0]], {"units": "e", "description": "d"})}},
         r"\[n_atoms, 3\]"),
        ({"custom": {"a": ((1.0, 2.0), 3, 4)}}, "pair"),
        ({"custom": {"a": ([1.0, 2.0], "e")}}, "mapping"),
        ({"meta": {"a": {"units": "e", "description": "d"}}}, "does not carry"),
    ])
    def test_bad_custom_roles_rejected_before_writing(self, kwargs, match):
        w = self._writer()
        before = json.loads(json.dumps(w.build_manifest()))
        with pytest.raises(ValueError, match=match):
            w.add_atom_properties(**kwargs)
        assert w.build_manifest() == before


def _atom_properties_archive(tmp_path, members, files, name="custom_meta.qvf"):
    """A one-structure archive with an atom_properties section built from
    explicit member specs (the writer refuses to produce the broken ones)."""
    import hashlib
    import zipfile

    struct_doc = json.dumps({"atoms": [
        {"symbol": "Li", "position": [0, 0, 0], "atomic_number": 3},
        {"symbol": "H", "position": [0, 0, 1.6], "atomic_number": 1}],
        "pbc": [False, False, False], "lattice_vectors": None}).encode()
    manifest = {
        "qvf_version": 1,
        "source": {"program": "t", "version": "1", "calculation": "c"},
        "sections": [
            {"id": "structure", "kind": "structure", "members": {"structure": {
                "path": "s.json", "format": "json",
                "sha256": hashlib.sha256(struct_doc).hexdigest()}}},
            {"id": "ap", "kind": "atom_properties", "members": members}]}
    path = str(tmp_path / name)
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("s.json", struct_doc)
        for fname, data in files.items():
            zf.writestr(fname, data)
        zf.writestr("manifest.json", json.dumps(manifest))
    return path


def _binary_spec(path, values, shape):
    import hashlib
    import struct
    raw = struct.pack(f"<{len(values)}d", *values)
    return ({"path": path, "format": "binary", "dtype": "float64", "shape": shape,
             "sha256": hashlib.sha256(raw).hexdigest()}, raw)


def _json_spec(path, doc):
    import hashlib
    raw = json.dumps(doc).encode()
    return ({"path": path, "format": "json",
             "sha256": hashlib.sha256(raw).hexdigest()}, raw)


class TestAtomPropertyCustomValidation:
    """Spec § 5.7, § 6: a custom role without a ``meta`` entry is
    ``E_MEMBER_MALFORMED`` (with or without a JSON Schema engine); the schema
    pins the shapes and the closed role list."""

    ENTRY = {"units": "e", "description": "d", "kind": "scalar"}

    def _validate(self, tmp_path, roles, meta_doc="absent", **kw):
        members, files = {}, {}
        for role, (values, shape) in roles.items():
            spec, raw = _binary_spec(f"ap/{role}.bin", values, shape)
            members[role], files[f"ap/{role}.bin"] = spec, raw
        if meta_doc != "absent":
            spec, raw = _json_spec("ap/meta.json", meta_doc)
            members["meta"], files["ap/meta.json"] = spec, raw
        path = _atom_properties_archive(tmp_path, members, files, **kw)
        return validate_qvf(path)

    @staticmethod
    def _codes(report):
        return sorted({(f["code"], f["where"]) for f in report["findings"]
                       if f["severity"] == "error"})

    def test_valid_custom_role_passes(self, tmp_path):
        r = self._validate(tmp_path, {"custom.a": ([1.0, 2.0], [2])},
                           {"custom.a": self.ENTRY})
        assert r["ok"], r["errors"]

    def test_no_meta_member(self, tmp_path):
        r = self._validate(tmp_path, {"custom.a": ([1.0, 2.0], [2]),
                                      "custom.b": ([1.0, 2.0], [2])})
        assert not r["ok"]
        assert ("E_MEMBER_MALFORMED", "/sections/1/members/custom.a") in self._codes(r)
        assert ("E_MEMBER_MALFORMED", "/sections/1/members/custom.b") in self._codes(r)

    def test_meta_without_the_entry(self, tmp_path):
        r = self._validate(tmp_path, {"custom.a": ([1.0, 2.0], [2]),
                                      "custom.b": ([1.0, 2.0], [2])},
                           {"custom.a": self.ENTRY})
        assert [c for c in self._codes(r) if c[0] == "E_MEMBER_MALFORMED"] == [
            ("E_MEMBER_MALFORMED", "/sections/1/members/custom.b")]

    @pytest.mark.parametrize("meta_doc", [[], "text", None,
                                          {"custom.a": None},
                                          {"custom.a": "scalar"},
                                          {"custom.A": ENTRY}])
    def test_meta_that_does_not_describe_the_role(self, tmp_path, meta_doc):
        r = self._validate(tmp_path, {"custom.a": ([1.0, 2.0], [2])}, meta_doc)
        assert ("E_MEMBER_MALFORMED", "/sections/1/members/custom.a") in self._codes(r)

    def test_check_runs_without_the_schema_engine(self, tmp_path, monkeypatch):
        import qvf_reader
        monkeypatch.setattr(qvf_reader, "_load_schema_validator", lambda: None)
        monkeypatch.setattr(qvf_reader, "_load_payload_validator", lambda _d: None)
        r = self._validate(tmp_path, {"custom.a": ([1.0, 2.0], [2])})
        assert ("E_MEMBER_MALFORMED", "/sections/1/members/custom.a") in self._codes(r)
        r = self._validate(tmp_path, {"custom.a": ([1.0, 2.0], [2])},
                           {"custom.a": self.ENTRY}, name="ok.qvf")
        assert not [f for f in r["findings"] if f["severity"] == "error"]

    def test_a_listed_role_needs_no_meta(self, tmp_path):
        r = self._validate(tmp_path, {"mulliken_charge": ([0.1, -0.1], [2])})
        assert r["ok"], r["errors"]

    def test_unreadable_meta_is_not_judged_twice(self, tmp_path):
        members, files = {}, {}
        members["custom.a"], files["ap/custom.a.bin"] = _binary_spec(
            "ap/custom.a.bin", [1.0, 2.0], [2])
        members["meta"], files["ap/meta.json"] = _json_spec("ap/meta.json", {})
        files["ap/meta.json"] = b"{not json"
        path = _atom_properties_archive(tmp_path, members, files)
        codes = {f["code"] for f in validate_qvf(path)["findings"]}
        assert "E_JSON_INVALID" in codes
        assert "E_MEMBER_MALFORMED" not in codes

    # -- the schema --------------------------------------------------------

    @pytest.mark.parametrize("roles, meta_doc", [
        ({"magnetic_moment": ([0.0, 1.0], [2])}, "absent"),               # [n]
        ({"force_eh_bohr": ([0.0] * 4, [2, 2])}, "absent"),               # [n, 2]
        ({"cm5_charge": ([0.0] * 6, [2, 3])}, "absent"),                  # vector shape
        ({"bogus_charge": ([0.1, -0.1], [2])}, "absent"),                 # unknown
        ({"custom.a": ([1.0, 2.0], [2])}, "absent"),                      # no meta
        ({"custom.a-b": ([1.0, 2.0], [2])}, {"custom.a-b": ENTRY}),       # name
        ({"custom.a": ([1.0] * 8, [2, 4])}, {"custom.a": ENTRY}),         # shape
        ({"custom.a": ([1.0, 2.0], [2])}, {"custom.a": {"units": "e"}}),  # entry
        ({"custom.a": ([1.0, 2.0], [2])},
         {"custom.a": {**ENTRY, "kind": "tensor"}}),
        ({"custom.a": ([1.0, 2.0], [2])}, {"custom.a": ENTRY, "stray": 1}),
    ])
    def test_schema_rejects(self, tmp_path, roles, meta_doc):
        pytest.importorskip("jsonschema")
        r = self._validate(tmp_path, roles, meta_doc)
        assert not r["ok"]
        assert {f["code"] for f in r["findings"] if f["severity"] == "error"} & \
            {"E_SCHEMA", "E_PAYLOAD_SCHEMA"}, r["findings"]

    @pytest.mark.parametrize("doc, valid", [
        ({}, True),
        ({"custom.a": ENTRY}, True),
        ({"custom.a": ENTRY, "custom.b_2": {**ENTRY, "kind": "vector",
                                            "x_note": [1]}}, True),
        ({"x_producer": "demo", "custom.a": ENTRY}, True),
        ([], False),
        ({"a": ENTRY}, False),                                  # not custom.<name>
        ({"custom.": ENTRY}, False),
        ({"custom.a-b": ENTRY}, False),
        ({"custom.a": {"units": "e", "description": "d"}}, False),     # no kind
        ({"custom.a": {"units": "e", "kind": "scalar"}}, False),       # no description
        ({"custom.a": {"description": "d", "kind": "scalar"}}, False), # no units
        ({"custom.a": {**ENTRY, "kind": "tensor"}}, False),
        ({"custom.a": {**ENTRY, "unit": "e"}}, False),                 # unknown key
        ({"custom.a": {**ENTRY, "units": 1}}, False),
        ({"custom.a": {**ENTRY, "description": ""}}, False),
        ({"custom.a": "scalar"}, False),
    ])
    def test_meta_payload_definition(self, doc, valid):
        pytest.importorskip("jsonschema")
        import qvf_reader
        validator = qvf_reader._load_payload_validator("AtomPropertiesMetaPayload")
        assert validator is not None
        assert (not list(validator.iter_errors(doc))) is valid, doc

    def test_meta_only_section_is_rejected(self, tmp_path):
        pytest.importorskip("jsonschema")
        r = self._validate(tmp_path, {}, {"custom.a": self.ENTRY})
        assert any(f["code"] == "E_SCHEMA" for f in r["findings"])

    def test_custom_roles_need_meta_in_the_schema_too(self, tmp_path):
        pytest.importorskip("jsonschema")
        r = self._validate(tmp_path, {"custom.a": ([1.0, 2.0], [2])})
        assert any(f["code"] == "E_SCHEMA" for f in r["findings"])


class TestSchemaNamesNewFields:
    """The new members and the DOS section fields are typed in the schema, and
    ``Section`` stays open (maintainer decision 3 is pending)."""

    @staticmethod
    def _manifest():
        w = QvfWriter(program="t", version="1", calculation="schema")
        w.add_structure([{"symbol": "Si", "position": [0, 0, 0], "atomic_number": 14}],
                        pbc=[True, True, True],
                        lattice_vectors=[[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])
        w.add_phonon_bands({"n_atoms": 1, "n_modes": 3}, np.zeros((2, 3)),
                           eigenvectors=np.zeros((2, 3, 1, 3)))
        w.add_phonon_dos({"channels": [{"atom": 0}]}, [0, 1, 2], [0, 1, 0],
                         projected=np.zeros((1, 3)))
        w.add_dos_total([-1, 0, 1], [0, 1, 0], smearing=0.1, smearing_type="gaussian",
                        fermi_energy_ev=0.0, n_electrons=8.0)
        w.add_dos_projected([-1, 0, 1], np.zeros((1, 3)),
                            channels=[{"atom": 0, "l": 1, "m": 0, "spin": 0,
                                       "label": "Si-3p"}])
        manifest = json.loads(json.dumps(w.build_manifest()))
        return manifest, {s["id"]: s for s in manifest["sections"]}

    def test_written_fields_validate(self):
        manifest, _ = self._manifest()
        assert not _schema_errors(manifest)

    @pytest.mark.parametrize("section, key, value", [
        ("phonon_bands", ("members", "eigenvectors", "shape"), [2, 3, 1]),
        ("phonon_bands", ("members", "eigenvectors", "shape"), [2, 3, 1, 2]),
        ("phonon_bands", ("members", "eigenvectors", "shape"), [2, 3, 1, 3, 3]),
        ("phonon_bands", ("members", "eigenvectors", "dtype"), "float32"),
        ("phonon_dos", ("members", "projected", "shape"), [3]),
        ("phonon_dos", ("members", "projected", "dtype"), "int32"),
        ("dos_total", ("n_spin",), 3),
        ("dos_total", ("n_spin",), 1.5),
        ("dos_total", ("smearing",), -0.1),
        ("dos_total", ("smearing_type",), 5),
        ("dos_total", ("fermi_energy_ev",), "0.0"),
        ("dos_total", ("n_electrons",), -1),
        ("dos_projected", ("n_spin",), 0),
        ("dos_projected", ("channels",), {"atom": 0}),
        ("dos_projected", ("channels", 0, "atom"), -1),
        ("dos_projected", ("channels", 0, "l"), -1),
        ("dos_projected", ("channels", 0, "m"), 0.5),
        ("dos_projected", ("channels", 0, "spin"), 2),
        ("dos_projected", ("channels", 0, "label"), 3),
    ])
    def test_malformed_value_fails(self, section, key, value):
        manifest, by_id = self._manifest()
        target = by_id[section]
        for step in key[:-1]:
            target = target[step]
        target[key[-1]] = value
        assert _schema_errors(manifest)

    def test_section_object_stays_open(self):
        manifest, by_id = self._manifest()
        by_id["dos_total"]["x_note"] = "unnamed fields still validate"
        by_id["dos_projected"]["channels"][0]["symbol"] = "Si"
        assert not _schema_errors(manifest)


class TestVolumesOfAWavefunction:
    """Spec § 5.2, Volumes of a wavefunction: the writer refuses a reference,
    row or spin the added wavefunction does not support."""

    @staticmethod
    def _writer():
        w = QvfWriter(program="t", version="1", calculation="volumes")
        w.add_structure([{"symbol": "H", "atomic_number": 1, "position": [0, 0, 0]},
                         {"symbol": "H", "atomic_number": 1, "position": [0, 0, 0.74]}])
        shells = [{"center": a, "l": 0, "exponents": [1.0], "coefficients": [1.0]}
                  for a in range(2)]
        w.add_wavefunction_gto(shells, mo_coefficients=np.eye(2), energies=[-0.5, 0.5],
                               occupations=[2, 0], section_id="wf")
        w.add_wavefunction_gto(shells, mo_coefficients_alpha=np.eye(2),
                               mo_coefficients_beta=np.eye(2)[:1],
                               energies_alpha=[-0.5, 0.5], occupations_alpha=[1, 0],
                               energies_beta=[-0.4], occupations_beta=[1],
                               section_id="wf_u")
        return w

    def _add(self, w, kind="volume.orbital", **fields):
        return w.add_volume(kind, np.zeros((2, 2, 2)), [0, 0, 0], np.eye(3), **fields)

    def test_fields_are_written_and_validate(self, tmp_path):
        w = self._writer()
        self._add(w, "volume.density", wavefunction_ref="wf", section_id="rho")
        self._add(w, wavefunction_ref="wf", mo_index=1, section_id="mo_1")
        self._add(w, wavefunction_ref="wf_u", mo_index=0, spin="beta", section_id="mo_b0")
        sections = {s["id"]: s for s in w.build_manifest()["sections"]}
        assert sections["rho"]["wavefunction_ref"] == "wf"
        assert "mo_index" not in sections["rho"]
        assert sections["mo_1"]["mo_index"] == 1 and "spin" not in sections["mo_1"]
        assert sections["mo_b0"]["spin"] == "beta"
        _assert_valid_archive(w, str(tmp_path / "linked.qvf"))

    @pytest.mark.parametrize("kind, fields, match", [
        ("volume.elf", {"wavefunction_ref": "wf"}, "applies to"),
        ("volume.density", {"wavefunction_ref": "wf", "mo_index": 0}, "volume.orbital only"),
        ("volume.orbital", {"mo_index": 0}, "needs wavefunction_ref"),
        ("volume.orbital", {"wavefunction_ref": "wf", "spin": "alpha"}, "needs mo_index"),
        ("volume.orbital", {"wavefunction_ref": "structure"}, "wavefunction.gto"),
        ("volume.orbital", {"wavefunction_ref": "wf", "mo_index": 2}, "outside"),
        ("volume.orbital", {"wavefunction_ref": "wf", "mo_index": -1}, "integer >= 0"),
        ("volume.orbital", {"wavefunction_ref": "wf", "mo_index": 0, "spin": "alpha"},
         "one set"),
        ("volume.orbital", {"wavefunction_ref": "wf_u", "mo_index": 0}, "unrestricted"),
        ("volume.orbital", {"wavefunction_ref": "wf_u", "mo_index": 1, "spin": "beta"},
         "outside"),
    ])
    def test_refusals(self, kind, fields, match):
        w = self._writer()
        before = len(w.build_manifest()["sections"])
        with pytest.raises(ValueError, match=match):
            self._add(w, kind, **fields)
        assert len(w.build_manifest()["sections"]) == before
