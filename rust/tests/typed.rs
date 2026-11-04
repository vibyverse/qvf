//! Typed kinds (`qvf::kinds`, generated from the schema `$defs`): every
//! section of every archive of the positive and legacy corpora decodes with
//! the model of its kind, writes back through the model, and the rewritten
//! archive holds the same section fields, payloads, tensors and opaque bytes,
//! decodes to an equal model and validates without errors.

use std::path::{Path, PathBuf};

use qvf::kinds::{AnyKind, ExcitedStates, SpectraIR, SpectraNMR, Structure, Trajectory, VolumeDensity, WavefunctionGTO};
use qvf::{validate_bytes, Code, Reader, Section, Severity, Source, TypedKind, ValidateOptions, Writer, WriterOptions};
use serde_json::{Map, Value};

fn corpus(dir: &str) -> Vec<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance").join(dir);
    let mut out: Vec<PathBuf> = std::fs::read_dir(&root)
        .unwrap_or_else(|e| panic!("{}: {e}", root.display()))
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "qvf"))
        .collect();
    out.sort();
    out
}

fn archive(name: &str) -> Reader {
    Reader::open(Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/corpus").join(name)).unwrap()
}

/// JSON equality with numbers compared as values: a typed `f64` field writes
/// the integer `2` back as `2.0`.
fn same_json(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Number(x), Value::Number(y)) => x.as_f64() == y.as_f64(),
        (Value::Array(x), Value::Array(y)) => x.len() == y.len() && x.iter().zip(y).all(|(p, q)| same_json(p, q)),
        (Value::Object(x), Value::Object(y)) => {
            x.len() == y.len() && x.iter().all(|(k, v)| y.get(k).is_some_and(|w| same_json(v, w)))
        }
        _ => a == b,
    }
}

fn fields(section: Section<'_>) -> Value {
    let mut map = section.raw.clone();
    map.remove("members");
    Value::Object(map)
}

/// Copy a vendor section member by member (it has no typed model).
fn copy_raw(reader: &Reader, section: Section<'_>, w: &mut Writer<Vec<u8>>) {
    let mut map = section.raw.clone();
    for key in ["id", "kind", "members"] {
        map.remove(key);
    }
    w.add_section(section.id, section.kind, map).unwrap();
    for (role, spec) in section.members() {
        let spec = spec.as_object().unwrap();
        match spec.get("format").and_then(Value::as_str) {
            Some("json") => drop(w.add_json(section.id, role, &reader.read_json(section.id, role).unwrap()).unwrap()),
            _ if spec.contains_key("dtype") => {
                let t = reader.read_tensor(section.id, role).unwrap();
                qvf::typed::put_tensor(w, section.id, role, &t).unwrap();
            }
            _ => {
                let mt = spec.get("media_type").and_then(Value::as_str);
                drop(w.add_opaque(section.id, role, &reader.read_bytes(section.id, role).unwrap(), mt).unwrap())
            }
        }
    }
}

/// Decode every section, write it back through its model and return the
/// rewritten archive with the decoded models.
fn rewrite(reader: &Reader) -> (Vec<u8>, Vec<AnyKind>) {
    let m = reader.manifest().as_object().unwrap();
    let src = m["source"].as_object().unwrap();
    let text = |k: &str| src.get(k).and_then(Value::as_str).unwrap_or("").to_owned();
    let source = Source::new(text("program"), text("version")).with_calculation(text("calculation"));
    let options = WriterOptions { created_utc: false, ..WriterOptions::default() };
    let mut w = Writer::with_options(Vec::new(), source, options).unwrap();
    for (key, value) in m {
        if !["qvf_version", "schema_uri", "byte_order", "source", "sections"].contains(&key.as_str()) {
            w.set_root(key, value.clone()).unwrap();
        }
    }
    let mut models = Vec::new();
    for section in reader.sections() {
        match AnyKind::from_section(reader, section) {
            Ok(Some(model)) => {
                assert_eq!(model.kind(), section.kind);
                assert_eq!(model.id(), section.id);
                model.to_writer(&mut w).unwrap_or_else(|e| panic!("{}: write: {e}", section.id));
                models.push(model);
            }
            Ok(None) => copy_raw(reader, section, &mut w),
            Err(e) => panic!("section {:?} ({}): {e}", section.id, section.kind),
        }
    }
    (w.finish().unwrap(), models)
}

fn check_same(a: &Reader, b: &Reader) -> Result<(), String> {
    let (sa, sb) = (a.sections(), b.sections());
    if sa.len() != sb.len() {
        return Err(format!("{} sections became {}", sa.len(), sb.len()));
    }
    for (x, y) in sa.iter().zip(&sb) {
        if !same_json(&fields(*x), &fields(*y)) {
            return Err(format!("section {:?}: fields {} became {}", x.id, fields(*x), fields(*y)));
        }
        let roles: Vec<&String> = x.members().keys().collect();
        if roles != y.members().keys().collect::<Vec<_>>() {
            return Err(format!("section {:?}: member roles changed", x.id));
        }
        for role in roles {
            let spec: &Map<String, Value> = x.members()[role].as_object().unwrap();
            let same = match spec.get("format").and_then(Value::as_str) {
                Some("json") => same_json(&a.read_json(x.id, role).unwrap(), &b.read_json(y.id, role).unwrap()),
                _ if spec.contains_key("dtype") => a.read_tensor(x.id, role).unwrap() == b.read_tensor(y.id, role).unwrap(),
                _ => {
                    a.read_bytes(x.id, role).unwrap() == b.read_bytes(y.id, role).unwrap()
                        && spec.get("media_type") == y.members()[role].get("media_type")
                }
            };
            if !same {
                return Err(format!("section {:?} member {role:?} changed", x.id));
            }
        }
    }
    Ok(())
}

#[test]
fn corpus_round_trips_through_typed_models() {
    let mut archives = corpus("corpus");
    archives.extend(corpus("legacy"));
    let (mut sections, mut failures) = (0, Vec::new());
    for path in &archives {
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        let reader = Reader::open(path).unwrap_or_else(|e| panic!("{name}: {e}"));
        let (bytes, models) = rewrite(&reader);
        sections += models.len();
        let report = validate_bytes(bytes.clone(), &ValidateOptions::default());
        let errors: Vec<_> = report.findings.iter().filter(|f| f.severity == Severity::Error).map(|f| f.code.as_str()).collect();
        if !errors.is_empty() {
            failures.push(format!("{name}: rewritten archive has errors {errors:?}"));
            continue;
        }
        let again = Reader::from_bytes(bytes, qvf::Limits::default()).unwrap();
        if let Err(e) = check_same(&reader, &again) {
            failures.push(format!("{name}: {e}"));
            continue;
        }
        for model in &models {
            let decoded = again.get_any(model.id()).unwrap().unwrap();
            if &decoded != model {
                failures.push(format!("{name}: section {:?} decodes differently after the rewrite", model.id()));
            }
        }
    }
    println!("typed round trip: {} archives, {sections} typed sections", archives.len());
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn structure_and_restricted_wavefunction() {
    let r = archive("tddft.qvf");
    let mol = r.sections_of("structure")[0].id.to_owned();
    let s: Structure = r.get(&mol).unwrap();
    assert!(!s.structure.atoms.is_empty());
    assert_eq!(s.structure.atoms[0].position.len(), 3);
    let wf_id = r.sections_of("wavefunction.gto")[0].id.to_owned();
    let wf: WavefunctionGTO = r.get(&wf_id).unwrap();
    let c = wf.mo_coefficients.as_ref().expect("restricted coefficients");
    assert_eq!(c.shape.len(), 2);
    assert!(wf.mo_coefficients_alpha.is_none());
    assert!(!wf.basis.shells.is_empty());
    assert_eq!(wf.mo_metadata.energies.as_ref().map(Vec::len), Some(c.shape[0]));
    let states: ExcitedStates = r.get(r.sections_of("excited_states")[0].id).unwrap();
    assert!(!states.states.states.is_empty());
}

#[test]
fn unrestricted_wavefunction() {
    let r = archive("nto_unrestricted.qvf");
    let wf = r
        .sections_of("wavefunction.gto")
        .into_iter()
        .map(|s| r.get::<WavefunctionGTO>(s.id).unwrap())
        .find(|w| w.mo_coefficients_alpha.is_some())
        .expect("an unrestricted wavefunction");
    assert!(wf.mo_coefficients.is_none());
    let (a, b) = (wf.mo_coefficients_alpha.unwrap(), wf.mo_coefficients_beta.unwrap());
    assert_eq!(a.shape[1], b.shape[1]);
}

#[test]
fn volumes_spectra_and_trajectory() {
    let r = archive("volume_density.qvf");
    let id = r.sections_of("volume.density")[0].id.to_owned();
    let v: VolumeDensity = r.get(&id).unwrap();
    let n: usize = v.grid.shape.iter().map(|&x| x as usize).product();
    assert_eq!(v.data.len(), n);
    let any = r.get_volume(&id).unwrap();
    assert_eq!((any.kind.as_str(), &any.grid, &any.data), ("volume.density", &v.grid, &v.data));
    assert_eq!(r.get::<Structure>(&id).unwrap_err().code, Code::E_KIND_UNKNOWN);
    assert_eq!(r.get::<Structure>("no-such-id").unwrap_err().code, Code::E_REF_UNRESOLVED);

    let r = archive("spectra_ir.qvf");
    let ir: SpectraIR = r.get(r.sections_of("spectra.ir")[0].id).unwrap();
    assert_eq!(SpectraIR::KIND, "spectra.ir");
    assert!(ir.spectrum.extra.is_empty() || !ir.spectrum.extra.contains_key("x"));
    let r = archive("nmr.qvf");
    let nmr: SpectraNMR = r.get(r.sections_of("spectra.nmr")[0].id).unwrap();
    assert!(!nmr.fields.id.is_empty());

    let r = archive("trajectory.qvf");
    let t: Trajectory = r.get(r.sections_of("trajectory")[0].id).unwrap();
    assert_eq!(t.coords.shape.len(), 3);
    assert_eq!(t.coords.shape[2], 3);
}
