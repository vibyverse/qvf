//! The trajectory, bands, DOS and spectra builders: every archive must have
//! no findings from the Rust and the Python validators and read back
//! unchanged.

mod common;

use qvf::reader::DType;
use qvf::{
    Atom, BandProjections, Bands, BondingAnalysis, Broadening, DosBonding, DosChannel, DosInfo, DosPair, DosProjected, DosTotal,
    ExcitedState, ExcitedStates, Kpath, KpathSegment, Peak, Precision, Reader, Spectrum, Structure, Trajectory, TrajectoryCell,
    Writer,
};
use serde_json::{json, Map, Value};

fn source() -> qvf::Source {
    qvf::Source::new("qvf-rust-test", "0.1.0").with_calculation("builders")
}

fn ramp(n: usize, scale: f64, offset: f64) -> Vec<f64> {
    (0..n).map(|i| offset + scale * ((i % 37) as f64 - 18.0) / 37.0).collect()
}

fn water_frames(n_frames: usize) -> Vec<f64> {
    let base = [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692];
    (0..n_frames).flat_map(|f| base.iter().enumerate().map(move |(k, x)| x + 0.001 * (f * (k + 1)) as f64)).collect()
}

fn units(pairs: &[(&str, &str)]) -> Map<String, Value> {
    pairs.iter().map(|(k, v)| ((*k).to_owned(), Value::from(*v))).collect()
}

/// An optimization, a float32 run, a periodic NPT run and a slab with a fixed cell.
fn add_trajectories<W: std::io::Write>(w: &mut Writer<W>) {
    let symbols: Vec<String> = ["O", "H", "H"].iter().map(|s| s.to_string()).collect();
    let mut opt = Trajectory::new("opt", symbols.clone(), water_frames(4));
    opt.kind = Some("optimization".into());
    opt.energies = Some(vec![-76.01, -76.02, -76.025, -76.026]);
    opt.energies_in_metadata = true;
    opt.units = Some(units(&[("energies", "hartree"), ("coords", "angstrom")]));
    w.add_trajectory(&opt).unwrap();

    let mut f32run = Trajectory::md("md_f32", symbols, water_frames(5), 0.5);
    f32run.precision = Precision::Float32;
    f32run.ensemble = Some("nve".into());
    f32run.label = Some("NVE, float32 frames".into());
    w.add_trajectory(&f32run).unwrap();

    // Periodic MD: two Si atoms in a breathing cubic cell (NPT).
    let n = 6;
    let a0 = 5.43;
    let mut coords = Vec::new();
    let mut cells = Vec::new();
    for f in 0..n {
        let a = a0 * (1.0 + 0.002 * f as f64);
        coords.extend([0.0, 0.0, 0.0, a / 4.0, a / 4.0, a / 4.0 + 0.001 * f as f64]);
        cells.push([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]]);
    }
    let volumes: Vec<f64> = cells.iter().map(|c| c[0][0] * c[1][1] * c[2][2]).collect();
    let mut npt = Trajectory::md("md_npt", vec!["Si".into(), "Si".into()], coords, 1.0);
    npt.ensemble = Some("npt".into());
    npt.pbc = Some([true, true, true]);
    npt.cell = Some(TrajectoryCell::PerFrame(cells));
    npt.velocities = Some(ramp(n * 6, 0.01, 0.0));
    npt.gradient = Some(ramp(n * 6, 0.02, 0.0));
    npt.times = Some((0..n).map(|f| f as f64).collect());
    npt.energies = Some(ramp(n, 0.001, -578.0));
    npt.cell_volumes = Some(volumes);
    npt.temperatures = Some(ramp(n, 10.0, 300.0));
    npt.units = Some(units(&[
        ("lattice", "angstrom"),
        ("velocities", "angstrom_per_fs"),
        ("gradient", "hartree_per_bohr"),
        ("times", "fs"),
        ("cell_volumes", "a3"),
        ("temperatures", "kelvin"),
    ]));
    w.add_trajectory(&npt).unwrap();

    let mut slab = Trajectory::new(
        "slab",
        vec!["C".into(), "C".into()],
        vec![0.0, 0.0, 0.0, 1.42, 0.0, 0.0, 0.0, 0.0, 0.01, 1.42, 0.0, -0.01],
    );
    slab.kind = Some("md".into());
    slab.pbc = Some([true, true, false]);
    slab.cell = Some(TrajectoryCell::Fixed([[2.46, 0.0, 0.0], [1.23, 2.1304, 0.0], [0.0, 0.0, 15.875]]));
    w.add_trajectory(&slab).unwrap();
}

fn add_bands_and_dos<W: std::io::Write>(w: &mut Writer<W>) {
    let (n_spin, n_k, n_bands, n_ch) = (2, 5, 3, 2);
    let mut kpath = Kpath {
        kpoints: (0..n_k).map(|k| [k as f64 * 0.125, 0.0, 0.0]).collect(),
        labels: vec!["G".into(), "".into(), "".into(), "".into(), "X".into()],
        segments: vec![KpathSegment::between("G", [0.0, 0.0, 0.0], "X", [0.5, 0.0, 0.0], n_k)],
        fermi_ev: Some(5.1),
        channels: vec![DosChannel::shell("Si-s", 0, 0), DosChannel::shell("Si-p", 0, 1)],
        title: Some("Si along G-X".into()),
        units: Some(units(&[("fermi", "ev"), ("eigenvalues", "ev")])),
        ..Kpath::default()
    };
    kpath.extra.insert("x_rusttest_note".into(), json!("vendor key"));
    let mut bands = Bands::new("bands", kpath, ramp(n_spin * n_k * n_bands, 8.0, 5.0), [n_spin, n_k, n_bands]);
    bands.projections =
        Some(BandProjections { data: ramp(n_spin * n_k * n_bands * n_ch, 0.5, 0.5), n_channels: n_ch, spin_axis: true });
    bands.occupations = Some((0..n_spin * n_k * n_bands).map(|i| if i % n_bands == 0 { 1.0 } else { 0.0 }).collect());
    bands.label = Some("spin-polarized bands".into());
    w.add_bands(&bands).unwrap();

    let mut plain = Bands::new("bands_rank3", Kpath::default(), ramp(n_k * 4, 6.0, 0.0), [1, n_k, 4]);
    plain.kpath.channels = vec![DosChannel::new("total")];
    plain.projections = Some(BandProjections { data: vec![1.0; n_k * 4], n_channels: 1, spin_axis: false });
    w.add_bands(&plain).unwrap();

    let energies: Vec<f64> = (0..11).map(|i| -5.0 + i as f64).collect();
    let mut total = DosTotal::new("dos", energies.clone(), ramp(2 * 11, 1.0, 1.0));
    total.n_spin = Some(2);
    total.info = DosInfo {
        smearing: Some(0.1),
        smearing_type: Some("gaussian".into()),
        fermi_energy_ev: Some(5.1),
        n_electrons: Some(8.0),
    };
    w.add_dos_total(&total).unwrap();
    w.add_dos_total(&DosTotal::new("dos_rank1", energies.clone(), ramp(11, 1.0, 1.0))).unwrap();

    let channels = vec![DosChannel::shell("Si1-s", 0, 0), DosChannel::shell("Si1-p", 0, 1), DosChannel::new("all-d")];
    let mut pdos = DosProjected::new("pdos", energies.clone(), ramp(2 * 3 * 11, 0.3, 0.4), channels.clone());
    pdos.n_spin = Some(2);
    pdos.info.fermi_energy_ev = Some(5.1);
    w.add_dos_projected(&pdos).unwrap();
    w.add_dos_projected(&DosProjected::new("pdos_rank2", energies.clone(), ramp(3 * 11, 0.3, 0.4), channels)).unwrap();

    let pairs = vec![
        DosPair { i: 0, j: 1, label: Some("Si1-Si2".into()), distance_ang: Some(2.35) },
        DosPair { i: 1, j: 0, ..DosPair::default() },
    ];
    let mut coop = DosBonding::new("coop", BondingAnalysis::Coop, energies.clone(), pairs.clone());
    coop.projections = ramp(2 * 11, 0.1, 0.0);
    coop.integrated = vec![0.3, 0.2];
    coop.info.fermi_energy_ev = Some(5.1);
    w.add_dos_bonding(&coop).unwrap();
    let mut cohp = DosBonding::new("cohp", BondingAnalysis::Cohp, energies, pairs);
    cohp.n_spin = Some(2);
    cohp.projections = ramp(2 * 2 * 11, -0.1, 0.0);
    cohp.integrated = vec![-1.0, -0.9, -1.1, -0.8];
    cohp.info.smearing = Some(0.05);
    w.add_dos_bonding(&cohp).unwrap();
}

fn add_spectra<W: std::io::Write>(w: &mut Writer<W>) {
    let mut es = ExcitedStates { id: "es".into(), ..ExcitedStates::default() };
    es.states = vec![ExcitedState::new(1, 4.2), ExcitedState::new(2, 5.0), ExcitedState::new(3, 6.1)];
    w.add_excited_states(&es).unwrap();

    let mut ir = Spectrum::new("spectra.ir", vec![1595.0, 3657.0, 3756.0], vec![70.0, 4.0, 45.0]);
    ir.peaks = vec![Peak {
        x: 1595.0,
        y: Some(70.0),
        width: Some(10.0),
        label: Some("bend".into()),
        mode_index: Some(0),
        ..Peak::default()
    }];
    ir.broadening = Some(Broadening { kind: "lorentzian".into(), fwhm: Some(10.0), units: Some("cm-1".into()) });
    ir.temperature_k = Some(298.15);
    w.add_spectrum(&ir).unwrap();

    let mut raman = Spectrum::new("spectra.raman", vec![520.0, 960.0], vec![1.0, 0.1]);
    raman.x_units = Some("cm-1".into());
    raman.y_units = Some("dimensionless".into());
    raman.raman_laser_nm = Some(532.0);
    raman.raman_depolarization_ratios = Some(vec![0.1, 0.75]);
    w.add_spectrum(&raman).unwrap();

    let mut uv = Spectrum::new("spectra.uvvis", vec![4.2, 5.0, 6.1], vec![0.01, 0.3, 0.05]);
    uv.states_ref = Some("es".into());
    uv.peaks = (1..=3).map(|i| Peak { x: [4.2, 5.0, 6.1][i - 1], state_index: Some(i as u32), ..Peak::default() }).collect();
    uv.uvvis_x_axis = Some("wavelength_nm".into());
    uv.broadening = Some(Broadening { kind: "gaussian".into(), fwhm: Some(0.3), units: Some("ev".into()) });
    w.add_spectrum(&uv).unwrap();

    let mut ecd = Spectrum::new("spectra.ecd", vec![4.2, 5.0], vec![12.0, -8.0]);
    ecd.x_units = Some("ev".into());
    w.add_spectrum(&ecd).unwrap();
    let mut vcd = Spectrum::new("spectra.vcd", vec![1100.0, 1250.0], vec![-3.0, 2.5]);
    vcd.x_units = Some("cm-1".into());
    w.add_spectrum(&vcd).unwrap();
    let mut generic = Spectrum::new("spectra.generic", vec![1.0, 2.0, 3.0], vec![0.5, 1.0, 0.25]);
    generic.id = "xas".into();
    generic.write_x_alias = true;
    generic.x_units = Some("ev".into());
    generic.y_units = Some("dimensionless".into());
    generic.label = Some("X-ray absorption".into());
    generic.extra.insert("x_rusttest_edge".into(), json!("K"));
    w.add_spectrum(&generic).unwrap();

    let nmr = json!({"isotope": "1H", "reference": "TMS",
        "chemical_shifts": [{"atom_index": 1, "symbol": "H", "isotropic_shift_ppm": 31.2}]});
    w.add_spectrum_payload("nmr", "spectra.nmr", &nmr, Some("1H shieldings")).unwrap();
    let epr = json!({"g_tensor": {"isotropic": 2.0023, "principal": [2.0021, 2.0023, 2.0025]}});
    w.add_spectrum_payload("epr", "spectra.epr", &epr, None).unwrap();
}

#[test]
fn every_builder_writes_a_valid_archive() {
    let path = common::scratch("rust_builders.qvf");
    let mut w = Writer::create(&path, source()).unwrap();
    w.add_structure(&Structure::periodic(
        "si",
        vec![Atom::new("Si", [0.0, 0.0, 0.0]), Atom::new("Si", [1.3575, 1.3575, 1.3575])],
        [[5.43, 0.0, 0.0], [0.0, 5.43, 0.0], [0.0, 0.0, 5.43]],
        [true, true, true],
    ))
    .unwrap();
    add_trajectories(&mut w);
    add_bands_and_dos(&mut w);
    add_spectra(&mut w);
    w.finish().unwrap();
    common::assert_clean(&path);

    let r = Reader::open(&path).unwrap();
    let kinds: Vec<&str> = r.sections().iter().map(|s| s.kind).collect();
    for k in qvf::sections::SPECTRUM_1D_KINDS.iter().chain(&[
        "spectra.nmr",
        "spectra.epr",
        "trajectory",
        "bands",
        "dos.total",
        "dos.projected",
        "dos.coop",
        "dos.cohp",
    ]) {
        assert!(kinds.contains(k), "{k} missing");
    }
    let coords = r.read_tensor("md_f32", "coords").unwrap();
    assert_eq!((coords.dtype, coords.shape.clone()), (DType::Float32, vec![5, 3, 3]));
    assert_eq!(r.read_json("md_f32", "metadata").unwrap()["precision"], "float32");
    let lattice = r.read_tensor("md_npt", "lattice").unwrap();
    assert_eq!(lattice.shape, vec![6, 3, 3]);
    assert!((lattice.at(&[5, 2, 2]).unwrap() - 5.43 * 1.01).abs() < 1e-12);
    assert_eq!(r.read_tensor("slab", "lattice").unwrap().shape, vec![3, 3]);
    assert_eq!(r.read_tensor("md_npt", "velocities").unwrap().shape, vec![6, 2, 3]);
    assert_eq!(r.read_tensor("bands", "projections").unwrap().shape, vec![2, 5, 3, 2]);
    assert_eq!(r.read_tensor("bands_rank3", "projections").unwrap().shape, vec![5, 4, 1]);
    assert_eq!(r.read_tensor("cohp", "integrated").unwrap().shape, vec![2, 2]);
    assert_eq!(r.read_tensor("dos_rank1", "dos").unwrap().shape, vec![11]);
    assert_eq!(r.section("pdos").unwrap().raw["n_spin"], 2);
    let kp = r.read_json("bands", "kpath").unwrap();
    assert_eq!((kp["n_spin"].as_u64(), kp["n_kpoints"].as_u64(), kp["fermi"].as_f64()), (Some(2), Some(5), Some(5.1)));
    assert_eq!(r.read_json("xas", "spectrum").unwrap()["x"], json!([1.0, 2.0, 3.0]));
    assert_eq!(r.member_spec("spectra_uvvis", "spectrum").unwrap()["path"], "spectra/spectra_uvvis.json");

    // The typed models read every new section.
    let t: qvf::kinds::Trajectory = r.get("md_npt").unwrap();
    assert!(t.lattice.is_some() && t.temperatures.is_some());
    let _: qvf::kinds::Bands = r.get("bands").unwrap();
    let _: qvf::kinds::SpectraUVVis = r.get("spectra_uvvis").unwrap();
    let _: qvf::kinds::DosCohp = r.get("cohp").unwrap();
}

#[test]
fn builder_misuse_is_refused_before_writing() {
    let mut w = Writer::new(Vec::new(), source()).unwrap();
    let sym = |n: usize| vec!["H".to_string(); n];
    // Coords that are not whole frames, a short per-frame member, a cell without pbc, pbc without a cell.
    assert!(w.add_trajectory(&Trajectory::new("t", sym(2), vec![0.0; 7])).is_err());
    let mut t = Trajectory::new("t", sym(1), vec![0.0; 6]);
    t.times = Some(vec![0.0]);
    assert!(w.add_trajectory(&t).is_err());
    let mut t = Trajectory::new("t", sym(1), vec![0.0; 6]);
    t.cell = Some(TrajectoryCell::Fixed([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]));
    assert!(w.add_trajectory(&t).is_err());
    let mut t = Trajectory::new("t", sym(1), vec![0.0; 6]);
    t.pbc = Some([true, false, false]);
    assert!(w.add_trajectory(&t).is_err());
    let mut t = Trajectory::new("t", sym(1), vec![0.0; 6]);
    t.units = Some(units(&[("coords", "bohr")]));
    assert!(w.add_trajectory(&t).is_err(), "a declaration must restate the fixed unit");
    let mut t = Trajectory::new("t", sym(1), vec![0.0; 6]);
    t.cell = Some(TrajectoryCell::PerFrame(vec![[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]]));
    t.pbc = Some([true, true, true]);
    assert!(w.add_trajectory(&t).is_err(), "one cell for two frames");

    assert!(w.add_bands(&Bands::new("b", Kpath::default(), vec![0.0; 5], [3, 1, 5])).is_err(), "n_spin 3");
    let mut b = Bands::new("b", Kpath::default(), vec![0.0; 4], [1, 2, 2]);
    b.kpath.labels = vec!["G".into()];
    assert!(w.add_bands(&b).is_err(), "one label for two k-points");
    let mut b = Bands::new("b", Kpath::default(), vec![0.0; 4], [1, 2, 2]);
    b.occupations = Some(vec![0.0; 3]);
    assert!(w.add_bands(&b).is_err());

    let mut d = DosTotal::new("d", vec![0.0, 1.0], vec![0.0; 4]);
    assert!(w.add_dos_total(&d).is_err(), "two spin rows without the spin axis");
    d.n_spin = Some(3);
    assert!(w.add_dos_total(&d).is_err());
    let p = DosProjected::new("p", vec![0.0, 1.0], vec![0.0; 2], vec![DosChannel::new("a"), DosChannel::new("b")]);
    assert!(w.add_dos_projected(&p).is_err());

    assert!(w.add_spectrum(&Spectrum::new("spectra.ir", vec![1.0, 2.0], vec![1.0])).is_err(), "one length");
    assert!(w.add_spectrum(&Spectrum::new("spectra.nmr", vec![1.0], vec![1.0])).is_err(), "not 1-D");
    let mut s = Spectrum::new("spectra.ir", vec![1.0], vec![1.0]);
    s.x_units = Some("wavenumbers".into());
    assert!(w.add_spectrum(&s).is_err(), "unknown unit name");
    let mut s = Spectrum::new("spectra.ir", vec![1.0], vec![1.0]);
    s.states_ref = Some("nothing".into());
    assert!(w.add_spectrum(&s).is_err(), "states_ref to a section not added");
    let mut s = Spectrum::new("spectra.ir", vec![1.0], vec![1.0]);
    s.uvvis_x_axis = Some("energy_ev".into());
    assert!(w.add_spectrum(&s).is_err(), "a uvvis block on an IR spectrum");
    assert!(w.add_spectrum_payload("n", "spectra.nmr", &json!({"chemical_shifts": [{"atom_index": -1}]}), None).is_err());

    // Nothing was written: the archive is empty but valid.
    assert_eq!(w.manifest()["sections"], json!([]));
    let bytes = w.finish().unwrap();
    assert!(common::rust_findings(&bytes).is_empty());
}
