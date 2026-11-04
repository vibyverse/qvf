//! Archives written by the Rust writer: valid for the Rust validator and the
//! Python reference validator (`python -m qvf validate`, skipped when Python
//! or the `qvf` package is unavailable), and read back unchanged by the Rust
//! reader.

use std::path::{Path, PathBuf};
use std::process::Command;

use qvf::reader::DType;
use qvf::sections::SCALAR_VOLUME_KINDS;
use qvf::{
    validate_bytes, Atom, Basis, Bond, Channel, Compression, Contribution, ExcitedState, ExcitedStates, NtoLink, NtoRole,
    Orbitals, Reader, SeriesAxis, Shell, SpinOrbitals, Structure, ValidateOptions, Volume, WavefunctionGto, Writer,
    WriterOptions,
};
use serde_json::{json, Map, Value};

const BOHR_PER_ANGSTROM: f64 = 1.8897261246;
const H_EXPONENTS: [f64; 3] = [3.42525091, 0.62391373, 0.1688554];
const H_COEFFICIENTS: [f64; 3] = [0.15432897, 0.53532814, 0.44463454];
const H2_BOND: f64 = 0.74;

/// Overlap of two contracted s functions on `N_i`-normalized primitives, a
/// distance `r` (bohr) apart.
fn s_overlap(r: f64) -> f64 {
    let mut s = 0.0;
    for (a, ca) in H_EXPONENTS.iter().zip(H_COEFFICIENTS) {
        for (b, cb) in H_EXPONENTS.iter().zip(H_COEFFICIENTS) {
            let p = a + b;
            s += ca * cb * (2.0 * (a * b).sqrt() / p).powf(1.5) * (-a * b / p * r * r).exp();
        }
    }
    s
}

struct H2 {
    overlap: Vec<f64>,
    /// Rows σg and σu, orthonormal in `overlap`.
    mo: Vec<f64>,
}

fn h2() -> H2 {
    let (s_aa, s_ab) = (s_overlap(0.0), s_overlap(H2_BOND * BOHR_PER_ANGSTROM));
    let g = 1.0 / (2.0 * (s_aa + s_ab)).sqrt();
    let u = 1.0 / (2.0 * (s_aa - s_ab)).sqrt();
    H2 { overlap: vec![s_aa, s_ab, s_ab, s_aa], mo: vec![g, g, u, -u] }
}

fn density(mo: &[f64], occupations: &[f64], n_ao: usize) -> Vec<f64> {
    let mut p = vec![0.0; n_ao * n_ao];
    for (i, n) in occupations.iter().enumerate() {
        for mu in 0..n_ao {
            for nu in 0..n_ao {
                p[mu * n_ao + nu] += n * mo[i * n_ao + mu] * mo[i * n_ao + nu];
            }
        }
    }
    p
}

fn h2_structure(id: &str) -> Structure {
    let mut s = Structure::molecule(id, vec![Atom::new("H", [0.0, 0.0, 0.0]), Atom::new("H", [0.0, 0.0, H2_BOND])]);
    s.bonds = vec![Bond { i: 0, j: 1, order: Some(1.0) }];
    s.charge = Some(0);
    s
}

fn h2_basis(structure: &str) -> Basis {
    let shell = |c| Shell::new(c, 0, H_EXPONENTS.to_vec(), H_COEFFICIENTS.to_vec());
    let mut b = Basis::new(structure, true, vec![shell(0), shell(1)]);
    b.n_electrons = Some(2);
    b
}

fn grid(n: usize) -> ([f64; 3], [[f64; 3]; 3], [usize; 3]) {
    ([-2.0, -2.0, -2.0], [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]], [n, n + 1, n + 2])
}

fn ramp(n: usize, scale: f64) -> Vec<f64> {
    (0..n).map(|i| scale * ((i % 97) as f64 - 48.0) / 97.0).collect()
}

/// What an archive should hold, checked through the reader.
#[derive(Default)]
struct Expect {
    json: Vec<(String, String, Value)>,
    f64s: Vec<(String, String, Vec<usize>, Vec<f64>)>,
    f32s: Vec<(String, String, Vec<usize>, Vec<f32>)>,
}

/// An archive with every builder: a restricted and an unrestricted
/// wavefunction, NTOs, all volume kinds and excited states.
fn write_full(options: WriterOptions) -> (Vec<u8>, Expect) {
    let mut ex = Expect::default();
    let source = qvf::Source::new("qvf-rust-test", "0.1.0").with_calculation("h2-full");
    let mut w = Writer::with_options(Vec::new(), source, options).unwrap();
    w.set_root("provenance", json!({"method": "rhf", "basis": "sto-3g"})).unwrap();

    // structures
    let mol = h2_structure("mol");
    w.add_structure(&mol).unwrap();
    ex.json.push(("mol".into(), "structure".into(), mol.payload()));
    let mut triplet = h2_structure("mol_triplet");
    triplet.multiplicity = Some(3);
    triplet.name = Some("dihydrogen, triplet".into());
    triplet.label = Some("H2 (T)".into());
    w.add_structure(&triplet).unwrap();
    let a = 3.0;
    let mut cell = Structure::periodic(
        "cell",
        vec![Atom::new("H", [0.0, 0.0, 0.0]), Atom::new("H", [0.0, 0.0, H2_BOND])],
        [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
        [true, true, true],
    );
    cell.charge = Some(0);
    w.add_structure(&cell).unwrap();
    ex.json.push(("cell".into(), "structure".into(), cell.payload()));

    // restricted wavefunction
    let h = h2();
    let occ = vec![2.0, 0.0];
    let mut rhf = WavefunctionGto::new(
        "wf",
        h2_basis("mol"),
        SpinOrbitals::Restricted(Orbitals::new(h.mo.clone(), vec![-0.578, 0.670], occ.clone())),
    );
    rhf.overlap = Some(h.overlap.clone());
    rhf.density_matrix = Some(density(&h.mo, &occ, 2));
    w.add_wavefunction_gto(&rhf).unwrap();
    ex.json.push(("wf".into(), "basis".into(), rhf.basis.payload()));
    ex.f64s.push(("wf".into(), "mo_coefficients".into(), vec![2, 2], h.mo.clone()));
    ex.f64s.push(("wf".into(), "overlap".into(), vec![2, 2], h.overlap.clone()));

    // unrestricted wavefunction (triplet: both electrons α)
    let (occ_a, occ_b) = (vec![1.0, 1.0], vec![0.0, 0.0]);
    let mut uhf = WavefunctionGto::new(
        "wf_triplet",
        h2_basis("mol_triplet"),
        SpinOrbitals::Unrestricted {
            alpha: Orbitals::new(h.mo.clone(), vec![-0.40, 0.10], occ_a.clone()),
            beta: Orbitals::new(h.mo.clone(), vec![0.20, 0.90], occ_b.clone()),
        },
    );
    uhf.overlap = Some(h.overlap.clone());
    uhf.density_matrices_ab = Some((density(&h.mo, &occ_a, 2), density(&h.mo, &occ_b, 2)));
    uhf.label = Some("UHF triplet".into());
    w.add_wavefunction_gto(&uhf).unwrap();
    ex.f64s.push(("wf_triplet".into(), "mo_coefficients_beta".into(), vec![2, 2], h.mo.clone()));

    // natural transition orbitals of state 1 (σg -> σu)
    for (id, role, row) in [("nto_s1_hole", NtoRole::Hole, 0), ("nto_s1_particle", NtoRole::Particle, 1)] {
        let coeffs = h.mo[row * 2..row * 2 + 2].to_vec();
        let mut nto =
            WavefunctionGto::new(id, h2_basis("mol"), SpinOrbitals::Restricted(Orbitals::new(coeffs, vec![], vec![1.0])));
        nto.basis.n_electrons = None;
        nto.orbital_kind = "nto".into();
        nto.nto = Some(NtoLink { state_ref: "excited_states".into(), state_index: 1, role });
        w.add_wavefunction_gto(&nto).unwrap();
    }

    // volumes: every scalar kind, then vector, complex and series
    let (origin, voxels, shape) = grid(4);
    let npts = shape.iter().product::<usize>();
    for kind in SCALAR_VOLUME_KINDS {
        let id = kind.replace('.', "_");
        let mut v = Volume::new(&id, kind, origin, voxels, shape);
        v.structure_ref = Some("mol".into());
        match *kind {
            "volume.difference" => v.operands = Some(("volume_density".into(), "volume_spin".into())),
            "volume.transition_density" => v.excited_state = Some(("excited_states".into(), Some(1))),
            "volume.orbital" => v.component = Some("sigma_g".into()),
            _ => {}
        }
        let data: Vec<f32> = ramp(npts, 0.1).into_iter().map(|x| x as f32).collect();
        w.add_volume(&v, &data).unwrap();
        ex.f32s.push((id, "data".into(), shape.to_vec(), data));
    }
    let vector = Volume::new("field", "volume.vector", origin, voxels, shape);
    let data = ramp(npts * 3, 1.0);
    w.add_volume(&vector, &data).unwrap();
    ex.f64s.push(("field".into(), "data".into(), vec![shape[0], shape[1], shape[2], 3], data));
    let complex = Volume::new("bloch_orbital", "volume.complex", origin, voxels, shape);
    w.add_volume(&complex, &ramp(npts * 2, 0.5)).unwrap();
    let mut series = Volume::new("density_series", "volume.series", origin, voxels, shape);
    series.axis = Some(SeriesAxis {
        kind: "time".into(),
        values: vec![0.0, 0.5, 1.0],
        units: Some("fs".into()),
        labels: Some(vec!["t0".into(), "t1".into(), "t2".into()]),
        trajectory_ref: None,
    });
    // streamed in slabs of one time step each
    let mut sink = w.begin_volume(&series, DType::Float32).unwrap();
    let mut all = Vec::new();
    for t in 0..3 {
        let slab: Vec<f32> = ramp(npts, 1.0 + t as f64).into_iter().map(|x| x as f32).collect();
        sink.write(&slab).unwrap();
        all.extend(slab);
    }
    sink.finish().unwrap();
    ex.f32s.push(("density_series".into(), "data".into(), vec![3, shape[0], shape[1], shape[2]], all));

    // a periodic volume: CHGCAR style, spanning the lattice exactly
    let n = 6;
    let step = a * BOHR_PER_ANGSTROM / n as f64;
    let mut periodic =
        Volume::new("cell_density", "volume.density", [0.0; 3], [[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, step]], [n; 3]);
    periodic.periodic = Some([true; 3]);
    periodic.includes_endpoint = Some(false);
    periodic.structure_ref = Some("cell".into());
    w.add_volume(&periodic, &vec![0.25f64; n * n * n]).unwrap();
    ex.json.push(("cell_density".into(), "grid".into(), periodic.grid_payload()));

    // excited states of the RHF reference
    let mut s1 = ExcitedState::new(1, 15.2);
    s1.multiplicity = Some(1);
    s1.symmetry = Some("B1u".into());
    s1.oscillator_strength = Some(0.48);
    s1.transition_dipole = Some([0.0, 0.0, 1.05]);
    s1.contributions = vec![Contribution { from_mo: 0, to_mo: 1, spin: None, weight: 1.0 }];
    s1.nto_ref = Some("nto_s1_hole".into());
    s1.density_ref = Some("volume_transition_density".into());
    let mut t1 = ExcitedState::new(2, 10.6);
    t1.multiplicity = Some(3);
    t1.oscillator_strength = Some(0.0);
    t1.contributions = vec![Contribution { from_mo: 0, to_mo: 1, spin: None, weight: 1.0 }];
    let mut meta = Map::new();
    meta.insert("method".into(), Value::from("tdhf"));
    meta.insert("n_states".into(), Value::from(2));
    meta.insert("tda".into(), Value::from(false));
    let es = ExcitedStates {
        id: "excited_states".into(),
        states: vec![s1, t1],
        meta: Some(meta),
        wavefunction_ref: Some("wf".into()),
        label: Some("TDHF".into()),
    };
    w.add_excited_states(&es).unwrap();

    (w.finish().unwrap(), ex)
}

fn archive_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

fn rust_findings(bytes: &[u8]) -> Vec<String> {
    let report = validate_bytes(bytes.to_vec(), &ValidateOptions::default());
    report.findings.iter().map(|f| format!("{} {} {}", f.code, f.location, f.message)).collect()
}

/// The findings of the Python reference validator, or `None` when Python or
/// the `qvf` package cannot be run.
fn python_findings(path: &Path) -> Option<Vec<String>> {
    let python_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("python");
    let interpreter = std::env::var("PYTHON").unwrap_or_else(|_| "python3".into());
    let probe = Command::new(&interpreter).env("PYTHONPATH", &python_dir).args(["-c", "import qvf"]).output().ok()?;
    if !probe.status.success() {
        return None;
    }
    let out = Command::new(&interpreter)
        .env("PYTHONPATH", &python_dir)
        .args(["-m", "qvf", "validate", "--json"])
        .arg(path)
        .output()
        .expect("python runs");
    let report: Value = serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("python validator printed no JSON ({e}): {}", String::from_utf8_lossy(&out.stderr)));
    let findings = report["findings"].as_array().expect("findings").iter();
    Some(findings.map(|f| f.to_string()).collect())
}

fn check_archive(name: &str, bytes: &[u8], ex: &Expect) {
    assert_eq!(rust_findings(bytes), Vec::<String>::new(), "Rust validator findings for {name}");

    let path = archive_path(name);
    std::fs::write(&path, bytes).unwrap();
    match python_findings(&path) {
        Some(f) => {
            // Without jsonschema the Python validator can only skip the schema check.
            let f: Vec<String> = f.into_iter().filter(|s| !s.contains("W_SCHEMA_SKIPPED")).collect();
            assert_eq!(f, Vec::<String>::new(), "Python validator findings for {name}");
        }
        None => eprintln!("SKIP Python validation of {name}: python3 with the qvf package is not available"),
    }

    let reader = Reader::open(&path).unwrap();
    for (sid, role, want) in &ex.json {
        assert_eq!(&reader.read_json(sid, role).unwrap(), want, "{sid}/{role}");
    }
    for (sid, role, shape, want) in &ex.f64s {
        let t = reader.read_tensor(sid, role).unwrap();
        assert_eq!((&t.shape, t.dtype), (shape, DType::Float64), "{sid}/{role}");
        assert_eq!(&t.to_f64_vec(), want, "{sid}/{role}");
    }
    for (sid, role, shape, want) in &ex.f32s {
        let t = reader.read_tensor(sid, role).unwrap();
        assert_eq!((&t.shape, t.dtype), (shape, DType::Float32), "{sid}/{role}");
        let want: Vec<f64> = want.iter().map(|&x| f64::from(x)).collect();
        assert_eq!(t.to_f64_vec(), want, "{sid}/{role}");
    }
}

#[test]
fn full_archive_is_valid_and_round_trips() {
    let (bytes, ex) = write_full(WriterOptions::default());
    let reader = Reader::from_bytes(bytes.clone(), qvf::Limits::default()).unwrap();
    assert_eq!(reader.manifest()["source"]["calculation"], "h2-full");
    assert!(reader.manifest()["metadata"]["created_utc"].as_str().is_some_and(|s| s.ends_with('Z')));
    assert_eq!(reader.zip().entries().last().unwrap().name, "manifest.json");
    assert_eq!(reader.zip().entries().last().unwrap().method, qvf::zip::METHOD_STORE);
    let kinds: Vec<&str> = reader.sections().iter().map(|s| s.kind).collect();
    for kind in ["structure", "wavefunction.gto", "volume.series", "volume.complex", "volume.vector", "excited_states"] {
        assert!(kinds.contains(&kind), "{kind} missing");
    }
    check_archive("rust_full.qvf", &bytes, &ex);
}

#[test]
fn zip64_archive_is_valid_and_round_trips() {
    let options = WriterOptions { force_zip64: true, ..WriterOptions::default() };
    let (bytes, ex) = write_full(options);
    assert!(bytes.windows(4).any(|w| w == [0x50, 0x4b, 0x06, 0x06]), "no ZIP64 end record");
    check_archive("rust_zip64.qvf", &bytes, &ex);
}

#[test]
fn stored_archive_is_valid_and_round_trips() {
    let options = WriterOptions { compression: Compression::Stored, created_utc: false, ..WriterOptions::default() };
    let (bytes, ex) = write_full(options);
    let reader = Reader::from_bytes(bytes.clone(), qvf::Limits::default()).unwrap();
    assert!(reader.zip().entries().iter().all(|e| e.method == qvf::zip::METHOD_STORE));
    assert!(reader.manifest().get("metadata").is_none());
    check_archive("rust_stored.qvf", &bytes, &ex);
}

#[test]
fn every_dtype_round_trips() {
    let mut w = Writer::new(Vec::new(), qvf::Source::new("qvf-rust-test", "0.1.0")).unwrap();
    w.add_section("t", "x_rusttest.tensors", Map::new()).unwrap();
    w.add_tensor("t", "i8", &[-1i8, 2, -3, 4], &[2, 2]).unwrap();
    w.add_tensor("t", "i16", &[-300i16, 2], &[2]).unwrap();
    w.add_tensor("t", "i32", &[-70_000i32], &[1]).unwrap();
    w.add_tensor("t", "i64", &[-(1i64 << 40)], &[1, 1]).unwrap();
    w.add_tensor("t", "u8", &[255u8], &[1]).unwrap();
    w.add_tensor("t", "u16", &[65_535u16], &[1]).unwrap();
    w.add_tensor("t", "u32", &[4_000_000_000u32], &[1]).unwrap();
    w.add_tensor("t", "u64", &[1u64 << 50], &[1]).unwrap();
    w.add_tensor("t", "f32", &[1.5f32, -2.25], &[2, 1]).unwrap();
    w.add_tensor("t", "f64", &[f64::MIN_POSITIVE], &[1]).unwrap();
    w.add_tensor::<f64>("t", "empty", &[], &[0, 3]).unwrap();
    w.add_opaque("t", "png", b"\x89PNG\r\n\x1a\n", Some("image/png")).unwrap();
    let bytes = w.finish().unwrap();
    assert_eq!(rust_findings(&bytes), Vec::<String>::new());
    let r = Reader::from_bytes(bytes, qvf::Limits::default()).unwrap();
    assert_eq!(r.read_tensor("t", "i8").unwrap().at(&[1, 0]), Some(-3.0));
    assert_eq!(r.read_tensor("t", "i16").unwrap().at(&[0]), Some(-300.0));
    assert_eq!(r.read_tensor("t", "i32").unwrap().at(&[0]), Some(-70_000.0));
    assert_eq!(r.read_tensor("t", "i64").unwrap().at(&[0, 0]), Some(-(2f64.powi(40))));
    assert_eq!(r.read_tensor("t", "u8").unwrap().at(&[0]), Some(255.0));
    assert_eq!(r.read_tensor("t", "u16").unwrap().at(&[0]), Some(65_535.0));
    assert_eq!(r.read_tensor("t", "u32").unwrap().at(&[0]), Some(4e9));
    assert_eq!(r.read_tensor("t", "u64").unwrap().at(&[0]), Some(2f64.powi(50)));
    assert_eq!(r.read_tensor("t", "f32").unwrap().at(&[1, 0]), Some(-2.25));
    assert_eq!(r.read_tensor("t", "f64").unwrap().at(&[0]), Some(f64::MIN_POSITIVE));
    assert_eq!(r.read_tensor("t", "empty").unwrap().shape, vec![0, 3]);
    assert_eq!(r.read_bytes("t", "png").unwrap(), b"\x89PNG\r\n\x1a\n");
    assert_eq!(r.zip().find("t/png.bin").unwrap().method, qvf::zip::METHOD_STORE);
}

#[test]
fn misuse_is_refused_before_writing() {
    let mut w = Writer::new(Vec::new(), qvf::Source::new("qvf-rust-test", "0.1.0")).unwrap();
    w.add_structure(&h2_structure("mol")).unwrap();
    let err = |r: Result<String, qvf::WriteError>| r.unwrap_err().to_string();
    assert!(err(w.add_structure(&h2_structure("mol"))).contains("duplicate section id"));
    assert!(err(w.add_structure(&h2_structure("camera"))).contains("reserved"));
    assert!(err(w.add_structure(&h2_structure("a b"))).contains("must match"));
    let mut charged = h2_structure("ion");
    charged.multiplicity = Some(2);
    assert!(err(w.add_structure(&charged)).contains("impossible"));
    let mut slab = h2_structure("slab");
    slab.pbc = [true, true, false];
    assert!(err(w.add_structure(&slab)).contains("lattice_vectors"));

    let h = h2();
    let orbitals = || SpinOrbitals::Restricted(Orbitals::new(h.mo.clone(), vec![], vec![2.0, 0.0]));
    let wf = WavefunctionGto::new("wf", h2_basis("nowhere"), orbitals());
    assert!(err(w.add_wavefunction_gto(&wf)).contains("names no structure"));
    let mut wf = WavefunctionGto::new("wf", h2_basis("mol"), orbitals());
    wf.basis.shells[1].center = 5;
    assert!(err(w.add_wavefunction_gto(&wf)).contains("centre 5"));
    let mut wf = WavefunctionGto::new("wf", h2_basis("mol"), orbitals());
    wf.basis.n_electrons = Some(3);
    assert!(err(w.add_wavefunction_gto(&wf)).contains("n_electrons"));
    let mut wf = WavefunctionGto::new("wf", h2_basis("mol"), orbitals());
    wf.overlap = Some(vec![1.0; 3]);
    assert!(err(w.add_wavefunction_gto(&wf)).contains("overlap"));
    w.add_wavefunction_gto(&WavefunctionGto::new("wf", h2_basis("mol"), orbitals())).unwrap();

    let (origin, voxels, shape) = grid(2);
    let v = Volume::new("rho", "volume.density", origin, voxels, shape);
    assert!(err(w.add_volume(&v, &[0.0f64; 3])).contains("do not fill"));
    let series = Volume::new("s", "volume.series", origin, voxels, shape);
    assert!(err(w.add_volume(&series, &[0.0f64; 3])).contains("needs an axis"));
    assert!(w.begin_volume(&v, DType::Int32).is_err());

    let mut s = ExcitedState::new(1, 5.0);
    s.contributions = vec![Contribution { from_mo: 0, to_mo: 7, spin: None, weight: 1.0 }];
    let es = ExcitedStates { id: "es".into(), states: vec![s], wavefunction_ref: Some("wf".into()), ..ExcitedStates::default() };
    assert!(err(w.add_excited_states(&es)).contains("outside"));
    let mut s = ExcitedState::new(1, 5.0);
    s.contributions = vec![Contribution { from_mo: 0, to_mo: 1, spin: Some(Channel::Beta), weight: 1.0 }];
    let es = ExcitedStates { id: "es".into(), states: vec![s.clone(), s], ..ExcitedStates::default() };
    assert!(err(w.add_excited_states(&es)).contains("above the previous"));

    // None of the refused calls left anything behind: the archive is valid.
    let bytes = w.finish().unwrap();
    assert_eq!(rust_findings(&bytes), Vec::<String>::new());
}

#[test]
fn short_tensor_poisons_the_writer() {
    let mut w = Writer::new(Vec::new(), qvf::Source::new("qvf-rust-test", "0.1.0")).unwrap();
    w.add_section("t", "x_rusttest.tensors", Map::new()).unwrap();
    let mut sink = w.begin_tensor("t", "a", DType::Float64, &[4]).unwrap();
    sink.write(&[1.0f64, 2.0]).unwrap();
    assert!(sink.write(&[1.0f64, 2.0, 3.0]).is_err());
    assert!(sink.finish().unwrap_err().to_string().contains("E_BINARY_SIZE"));
    assert!(w.finish().is_err());
}
