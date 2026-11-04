//! The editor (spec § 2.3) on a corpus archive: append, replace and remove,
//! in place and atomically, then repack. The Rust and the Python validators
//! must agree on every edited archive, report only the orphans
//! (`W_ORPHAN_MEMBER`) before the repack and nothing after it, and the
//! untouched members must keep their bytes, `sha256` and `size`.

mod common;

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use qvf::{
    validate_bytes, Code, Compression, Editor, EditorOptions, Limits, Peak, Reader, Spectrum, Trajectory, ValidateOptions,
};
use serde_json::{json, Value};

fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/corpus").join(name)
}

fn copy_of(name: &str, as_name: &str) -> PathBuf {
    let path = common::scratch(as_name);
    std::fs::copy(corpus(name), &path).unwrap();
    path
}

/// The Rust findings as sorted (severity, code, location).
fn rust_set(path: &Path, limits: Limits) -> Vec<(String, String, String)> {
    let report = validate_bytes(std::fs::read(path).unwrap(), &ValidateOptions { limits, ..ValidateOptions::default() });
    let mut v: Vec<_> = report
        .findings
        .iter()
        .map(|f| (f.severity.as_str().to_owned(), f.code.as_str().to_owned(), f.location.clone()))
        .collect();
    v.sort();
    v
}

/// Both validators agree, and every finding is an orphan warning; returns
/// the number of orphans.
fn orphans_only(path: &Path) -> usize {
    let rust = rust_set(path, Limits::default());
    for (severity, code, location) in &rust {
        assert_eq!((severity.as_str(), code.as_str()), ("warning", "W_ORPHAN_MEMBER"), "{location} in {}", path.display());
    }
    match common::python_findings(path) {
        Some(mut py) => {
            py.sort();
            assert_eq!(py, rust, "Python and Rust differ on {}", path.display());
        }
        None => eprintln!("SKIP Python validation of {}", path.display()),
    }
    rust.len()
}

fn member<'a>(r: &'a Reader, id: &str, role: &str) -> &'a serde_json::Map<String, Value> {
    r.member_spec(id, role).unwrap()
}

#[test]
fn edits_follow_the_update_rules_and_repack_drops_orphans() {
    let path = copy_of("tddft.qvf", "edited.qvf");
    let before = Reader::open(&path).unwrap();
    let wf_spec = member(&before, "wf", "mo_coefficients").clone();
    let wf_bytes = before.read_bytes("wf", "mo_coefficients").unwrap();
    let structure = before.read_json("structure", "structure").unwrap();
    let density = before.read_tensor("tden_s1", "data").unwrap();
    let created = before.manifest()["metadata"]["created_utc"].clone();
    let original_len = std::fs::metadata(&path).unwrap().len();
    let original_names: Vec<String> = before.zip().entries().iter().map(|e| e.name.clone()).collect();
    drop(before);

    // 1: in place. Append sections that reference the archive's own
    // sections, and replace a JSON and a tensor member.
    let mut ed = Editor::open(&path).unwrap();
    let mut uv = Spectrum::new("spectra.uvvis", vec![7.1, 9.3], vec![0.05, 0.2]);
    uv.id = "uv_new".into();
    uv.states_ref = Some("excited_states".into());
    uv.peaks = vec![Peak { x: 7.1, state_index: Some(1), ..Peak::default() }];
    ed.writer().add_spectrum(&uv).unwrap();
    let mut bad = uv.clone();
    bad.id = "bad".into();
    bad.peaks[0].state_index = Some(99);
    assert!(ed.writer().add_spectrum(&bad).is_err(), "the archive's state indices are known to the builders");
    let mut ir = Spectrum::new("spectra.ir", vec![1595.0], vec![70.0]);
    ir.id = "ir".into();
    ed.writer().add_spectrum(&ir).unwrap();
    let mut opt = Trajectory::new("opt", vec!["O".into(), "H".into(), "H".into()], vec![0.0; 18]);
    opt.energies = Some(vec![-76.0, -76.01]);
    ed.writer().add_trajectory(&opt).unwrap();
    ed.replace_json("structure", "structure", &structure).unwrap();
    let scaled: Vec<f32> = density.to_f64_vec().iter().map(|&x| (2.0 * x) as f32).collect();
    ed.replace_tensor("tden_s1", "data", &scaled, &density.shape).unwrap();
    assert!(ed.writer().add_json("structure", "structure", &json!({})).is_err(), "the role is taken again");
    ed.commit().unwrap();

    let r = Reader::open(&path).unwrap();
    assert_eq!(member(&r, "structure", "structure")["path"], "structure/structure.2.json");
    assert_eq!(member(&r, "tden_s1", "data")["path"], "volumes/tden_s1.2.dat");
    assert_eq!(member(&r, "wf", "mo_coefficients"), &wf_spec, "an untouched member keeps its spec");
    assert_eq!(r.read_bytes("wf", "mo_coefficients").unwrap(), wf_bytes);
    assert_eq!(r.read_json("structure", "structure").unwrap(), structure);
    assert_eq!(r.read_tensor("tden_s1", "data").unwrap().to_f64_vec()[5], f64::from(scaled[5]));
    assert_eq!(r.manifest()["metadata"]["created_utc"], created);
    assert!(r.manifest()["metadata"]["modified_utc"].as_str().is_some_and(|s| s.ends_with('Z')));
    let names: Vec<&str> = r.zip().entries().iter().map(|e| e.name.as_str()).collect();
    assert_eq!(names.iter().filter(|n| **n == "manifest.json").count(), 2);
    assert_eq!(names.last(), Some(&"manifest.json"), "the new manifest is the last entry");
    let first_new = names.iter().position(|n| *n == "spectra/uv_new.json").unwrap();
    assert_eq!(names[..first_new], original_names, "old entries first, in place, then the appended ones");
    assert!(std::fs::metadata(&path).unwrap().len() > original_len);
    drop(r);
    // the old manifest, the old structure and the old density
    assert_eq!(orphans_only(&path), 3);

    // 2: atomic. Replace the structure again and remove a section.
    let options = EditorOptions { atomic: true, compression: Compression::Stored, ..EditorOptions::default() };
    let mut ed = Editor::open_with(&path, options).unwrap();
    ed.replace_json("structure", "structure", &structure).unwrap();
    let removed = ed.remove_section("ir").unwrap();
    assert_eq!(removed["kind"], "spectra.ir");
    assert!(ed.remove_section("ir").is_err());
    ed.set_run_status("converged", None).unwrap();
    assert!(ed.set_run_status("pending", Some("wf")).is_err());
    ed.commit().unwrap();
    let r = Reader::open(&path).unwrap();
    assert_eq!(member(&r, "structure", "structure")["path"], "structure/structure.3.json");
    assert!(r.section("ir").is_none());
    assert_eq!(r.manifest()["provenance"]["run_status"], "converged");
    drop(r);
    // + the second manifest, structure.2.json and the spectrum of "ir"
    assert_eq!(orphans_only(&path), 6);

    // The orphans count towards the total-size cap (spec § 2.3, § 6.1).
    let repacked = common::scratch("repacked.qvf");
    Editor::repack(&path, &repacked, None).unwrap();
    let live = Reader::open(&repacked).unwrap().zip().total_size();
    let edited = Reader::open(&path).unwrap().zip().total_size();
    assert!(edited > live);
    let cap = Limits { max_total_bytes: live as f64, ..Limits::default() };
    assert!(rust_set(&path, cap).iter().any(|(_, c, _)| c == Code::E_TOTAL_TOO_LARGE.as_str()));
    assert!(!rust_set(&repacked, cap).iter().any(|(_, c, _)| c == Code::E_TOTAL_TOO_LARGE.as_str()));

    // The repack holds the referenced entries and the manifest, nothing else.
    common::assert_clean(&repacked);
    let a = Reader::open(&path).unwrap();
    let b = Reader::open(&repacked).unwrap();
    assert_eq!(a.manifest(), b.manifest(), "a repack changes no member spec");
    let referenced: BTreeSet<String> = a
        .sections()
        .iter()
        .flat_map(|s| s.members().values().filter_map(|m| m["path"].as_str().map(str::to_owned)).collect::<Vec<_>>())
        .collect();
    let mut entries: BTreeSet<String> = b.zip().entries().iter().map(|e| e.name.clone()).collect();
    assert!(entries.remove("manifest.json"));
    assert_eq!(entries, referenced);
    for s in a.sections() {
        for role in s.members().keys() {
            assert_eq!(a.read_bytes(s.id, role).unwrap(), b.read_bytes(s.id, role).unwrap(), "{}/{role}", s.id);
        }
    }
    // STORE everywhere on request, and in place over the source.
    Editor::repack(&repacked, &repacked, Some(Compression::Stored)).unwrap();
    let c = Reader::open(&repacked).unwrap();
    assert!(c.zip().entries().iter().all(|e| e.method == qvf::zip::METHOD_STORE));
    assert_eq!(c.manifest(), a.manifest());

    // The Python editor reads and repacks the Rust-edited archive the same way.
    if let Some(mut py) = common::python() {
        let out = common::scratch("repacked_by_python.qvf");
        let status = py
            .args(["-c", "import sys; from qvf import Editor; Editor.repack(sys.argv[1], sys.argv[2])"])
            .arg(&path)
            .arg(&out)
            .status()
            .unwrap();
        assert!(status.success());
        let p = Reader::open(&out).unwrap();
        let names = |r: &Reader| r.zip().entries().iter().map(|e| e.name.clone()).collect::<Vec<_>>();
        assert_eq!(names(&p), names(&Reader::open(common::scratch("repacked.qvf")).unwrap()));
        assert_eq!(p.manifest(), a.manifest());
    }
}

#[test]
fn an_aborted_or_dropped_update_leaves_the_archive_as_it_was() {
    let path = copy_of("tddft.qvf", "aborted.qvf");
    let original = std::fs::read(&path).unwrap();
    let mut ed = Editor::open(&path).unwrap();
    ed.writer().add_spectrum(&Spectrum::new("spectra.ir", vec![1.0; 50_000], vec![1.0; 50_000])).unwrap();
    ed.remove_section("wf").unwrap();
    ed.abort();
    assert_eq!(std::fs::read(&path).unwrap(), original);
    {
        let mut ed = Editor::open(&path).unwrap();
        ed.replace_json("structure", "structure", &json!({"atoms": []})).unwrap();
    }
    assert_eq!(std::fs::read(&path).unwrap(), original, "a dropped editor rolls back");

    let options = EditorOptions { atomic: true, ..EditorOptions::default() };
    let mut ed = Editor::open_with(&path, options).unwrap();
    ed.writer().add_spectrum(&Spectrum::new("spectra.ir", vec![1.0], vec![1.0])).unwrap();
    ed.abort();
    assert_eq!(std::fs::read(&path).unwrap(), original);
    let dir = path.parent().unwrap();
    let leftovers: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().starts_with(".qvf-"))
        .collect();
    assert!(leftovers.is_empty(), "temporary files left behind: {leftovers:?}");
}

#[test]
fn updates_keep_job_specs_and_refuse_ambiguous_archives() {
    let path = copy_of("job_spec.qvf", "job.qvf");
    let r = Reader::open(&path).unwrap();
    let job = r.sections_of("job.spec")[0].id.to_owned();
    let role = r.section(&job).unwrap().members().keys().next().cloned();
    drop(r);
    let mut ed = Editor::open(&path).unwrap();
    assert!(ed.remove_section(&job).is_err(), "job.spec is kept verbatim (spec § 5.9)");
    if let Some(role) = role {
        assert!(ed.remove_member(&job, &role).is_err());
    }
    ed.set_run_status("running", None).unwrap();
    ed.commit().unwrap();
    assert_eq!(orphans_only(&path), 1);

    let dup = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/negative/zip_duplicate_entry.qvf");
    if dup.exists() {
        let copy = common::scratch("dup.qvf");
        std::fs::copy(&dup, &copy).unwrap();
        assert!(Editor::open(&copy).is_err());
        assert!(Editor::repack(&copy, common::scratch("dup_out.qvf"), None).is_err());
    }
}

#[test]
fn a_withdrawn_version_2_archive_is_written_back_as_version_1() {
    let legacy = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/legacy/reaction_path_v2.qvf");
    let path = common::scratch("v2.qvf");
    std::fs::copy(&legacy, &path).unwrap();
    let mut ed = Editor::open(&path).unwrap();
    ed.writer().add_spectrum(&Spectrum::new("spectra.ir", vec![1.0], vec![1.0])).unwrap();
    ed.commit().unwrap();
    let r = Reader::open(&path).unwrap();
    assert_eq!(r.manifest()["qvf_version"], 1);
    assert_eq!(r.manifest()["schema_uri"], qvf::generated::SCHEMA_ID);
    drop(r);
    let findings = rust_set(&path, Limits::default());
    assert!(findings.iter().all(|(s, _, _)| s == "warning"), "{findings:?}");
}
