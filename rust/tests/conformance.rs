//! The conformance runner: the Rust port of `python/qvf/conformance.py`.
//!
//! * `conformance/corpus`: every archive validates without errors, has the
//!   expected section kinds and source program, and reproduces every check of
//!   its `expected.json` through the reader.
//! * `conformance/negative`: every archive is judged as its `expected.json`
//!   says: rejected with at least the expected finding codes, or, for a
//!   warning-only archive, valid with the warning. `limits` overrides apply.
//! * `conformance/legacy`: the same rules as `negative/`, run separately.
//!
//! An archive the Rust implementation cannot judge yet would be listed in
//! `SKIPPED` with a reason and reported as skipped, never as passing. Each
//! test prints `passed / failed / skipped` counts (`cargo test -- --nocapture`).

use std::path::{Path, PathBuf};

use qvf::json::json_equal;
use qvf::{validate_file, Code, Limits, Reader, Severity, Tensor, TensorData, ValidateOptions};
use serde_json::Value;

/// Archives the Rust validator does not cover yet: (corpus, archive, reason).
/// Empty: every archive of every corpus is run.
const SKIPPED: &[(&str, &str, &str)] = &[];

fn conformance_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("conformance")
}

fn expected_files(dir: &Path) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("cannot list {}: {e}", dir.display()))
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.ends_with(".expected.json")))
        .collect();
    files.sort();
    files
}

fn read_json(path: &Path) -> Value {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{} is not JSON: {e}", path.display()))
}

/// The archive an expectation names, which must stay inside its directory.
fn archive_path(dir: &Path, name: &str) -> Result<PathBuf, String> {
    if name.contains('/') || name.contains('\\') || name == ".." {
        return Err(format!("corpus archive {name:?} must stay within the corpus directory"));
    }
    Ok(dir.join(name))
}

fn skip_reason(corpus: &str, archive: &str) -> Option<&'static str> {
    SKIPPED.iter().find(|(c, a, _)| *c == corpus && *a == archive).map(|(_, _, r)| *r)
}

/// Walk a value by keys and indices (negative indices count from the end, as in Python).
fn navigate<'a>(mut cur: &'a Value, pointer: &Value) -> Result<&'a Value, String> {
    for step in pointer.as_array().ok_or("pointer is not an array")? {
        cur = match (cur, step) {
            (Value::Object(o), Value::String(k)) => o.get(k).ok_or_else(|| format!("no key {k:?}"))?,
            (Value::Array(a), Value::Number(n)) => {
                let i = n.as_i64().ok_or("index is not an integer")?;
                let i = if i < 0 { a.len() as i64 + i } else { i };
                a.get(usize::try_from(i).map_err(|_| "index out of range")?).ok_or("index out of range")?
            }
            _ => return Err(format!("cannot step {step} into {}", if cur.is_array() { "an array" } else { "a non-container" })),
        };
    }
    Ok(cur)
}

fn element(t: &Tensor, k: usize) -> Value {
    match &t.data {
        TensorData::Int8(v) => Value::from(v[k]),
        TensorData::Int16(v) => Value::from(v[k]),
        TensorData::Int32(v) => Value::from(v[k]),
        TensorData::Int64(v) => Value::from(v[k]),
        TensorData::UInt8(v) => Value::from(v[k]),
        TensorData::UInt16(v) => Value::from(v[k]),
        TensorData::UInt32(v) => Value::from(v[k]),
        TensorData::UInt64(v) => Value::from(v[k]),
        TensorData::Float32(v) => Value::from(f64::from(v[k])),
        TensorData::Float64(v) => Value::from(v[k]),
    }
}

/// The tensor as nested JSON lists, like NumPy's `tolist()`.
fn to_list(t: &Tensor) -> Value {
    fn build(t: &Tensor, axis: usize, offset: usize) -> Value {
        if axis == t.shape.len() {
            return element(t, offset);
        }
        let stride: usize = t.shape[axis + 1..].iter().product();
        Value::Array((0..t.shape[axis]).map(|i| build(t, axis + 1, offset + i * stride)).collect())
    }
    build(t, 0, 0)
}

fn equal(got: &Value, want: &Value, tol: Option<f64>) -> bool {
    match tol {
        Some(tol) => match (got.as_f64(), want.as_f64()) {
            (Some(g), Some(w)) => (g - w).abs() <= tol,
            _ => false,
        },
        None => json_equal(got, want),
    }
}

fn run_check(reader: &Reader, chk: &Value) -> Result<(), String> {
    let tol = chk.get("tol").and_then(Value::as_f64);
    let want = chk.get("equals").cloned().unwrap_or(Value::Null);
    let got: Value = if let Some(p) = chk.get("manifest_pointer") {
        navigate(reader.manifest(), p)?.clone()
    } else {
        let section = chk.get("section").and_then(Value::as_str).ok_or("check without section")?;
        let member = chk.get("member").and_then(Value::as_str).ok_or("check without member")?;
        if let Some(p) = chk.get("json_pointer") {
            let doc = reader.read_json(section, member).map_err(|e| e.to_string())?;
            navigate(&doc, p)?.clone()
        } else if let Some(p) = chk.get("binary_index") {
            let t = reader.read_tensor(section, member).map_err(|e| e.to_string())?;
            navigate(&to_list(&t), p)?.clone()
        } else if chk.get("decode").and_then(Value::as_str) == Some("utf8") {
            let bytes = reader.read_bytes(section, member).map_err(|e| e.to_string())?;
            let text = String::from_utf8(bytes).map_err(|e| format!("not UTF-8: {e}"))?;
            if let Some(needle) = chk.get("contains").and_then(Value::as_str) {
                return if text.contains(needle) { Ok(()) } else { Err(format!("text does not contain {needle:?}")) };
            }
            Value::String(text)
        } else {
            return Err(format!("malformed check (no pointer): {chk}"));
        }
    };
    if equal(&got, &want, tol) {
        Ok(())
    } else {
        Err(format!("got {got}, want {want}{}", tol.map(|t| format!(" (tol {t})")).unwrap_or_default()))
    }
}

fn run_positive_one(path: &Path, expected: &Value) -> Vec<String> {
    let report = validate_file(path, &ValidateOptions::default());
    if !report.ok {
        return report
            .findings
            .iter()
            .filter(|f| f.severity == Severity::Error)
            .map(|f| format!("validation failed: {} {} {}", f.code, f.location, f.message))
            .collect();
    }
    let reader = match Reader::open(path) {
        Ok(r) => r,
        Err(e) => return vec![format!("reader refused the archive: {e}")],
    };
    let mut fails = Vec::new();
    let mut got_kinds: Vec<String> = reader.sections().iter().map(|s| s.kind.to_string()).collect();
    got_kinds.sort();
    got_kinds.dedup();
    if let Some(want) = expected.get("kinds").and_then(Value::as_array) {
        let mut want: Vec<String> = want.iter().filter_map(|k| k.as_str().map(String::from)).collect();
        want.sort();
        want.dedup();
        if got_kinds != want {
            fails.push(format!("kinds mismatch: got {got_kinds:?}, want {want:?}"));
        }
    }
    if let Some(want) = expected.get("source_program") {
        let got = reader.manifest().get("source").and_then(|s| s.get("program")).cloned().unwrap_or(Value::Null);
        if &got != want {
            fails.push(format!("source.program: got {got}, want {want}"));
        }
    }
    for chk in expected.get("checks").and_then(Value::as_array).into_iter().flatten() {
        if let Err(e) = run_check(&reader, chk) {
            let at = format!(
                "{}/{}",
                chk.get("section").and_then(Value::as_str).unwrap_or("manifest"),
                chk.get("member").and_then(Value::as_str).unwrap_or("")
            );
            fails.push(format!("check failed ({at}): {e}"));
        }
    }
    fails
}

/// The judgement of one archive of the negative format (`negative/`, `legacy/`).
/// This validator always runs the JSON Schema check, so the schema-only codes
/// bind it too.
fn run_negative_one(path: &Path, expected: &Value) -> Vec<String> {
    let limits = match Limits::from_json(expected.get("limits")) {
        Ok(l) => l,
        Err(e) => return vec![format!("bad limits in expected.json: {e}")],
    };
    let options = ValidateOptions { limits, ..ValidateOptions::default() };
    let report = validate_file(path, &options);
    let got: Vec<(&str, &str)> = report.findings.iter().map(|f| (f.code.as_str(), f.severity.as_str())).collect();
    let mut fails = Vec::new();
    let mut expects_error = false;
    for want in expected.get("expect").and_then(Value::as_array).into_iter().flatten() {
        let code = want.get("code").and_then(Value::as_str).unwrap_or("");
        let severity = want.get("severity").and_then(Value::as_str).unwrap_or("");
        if Code::parse(code).is_none() {
            fails.push(format!("expected.json names {code:?}, which is not an Appendix B code"));
        }
        expects_error |= severity == "error";
        if !got.contains(&(code, severity)) {
            let mut codes: Vec<&str> = got.iter().map(|(c, _)| *c).collect();
            codes.sort_unstable();
            fails.push(format!("expected {severity} {code}, got {codes:?}"));
        }
    }
    let errors: Vec<&str> = report.findings.iter().filter(|f| f.severity == Severity::Error).map(|f| f.code.as_str()).collect();
    if expects_error && errors.is_empty() {
        fails.push("validator accepted a must-fail archive".into());
    }
    if !expects_error && !errors.is_empty() {
        fails.push(format!("archive must validate ok here; errors: {errors:?}"));
    }
    fails
}

struct Tally {
    label: &'static str,
    passed: usize,
    failed: Vec<String>,
    skipped: Vec<String>,
}

impl Tally {
    fn report(&self) {
        let total = self.passed + self.failed.len() + self.skipped.len();
        println!(
            "{} conformance (Rust): {} passed, {} failed, {} skipped of {total}",
            self.label,
            self.passed,
            self.failed.len(),
            self.skipped.len()
        );
        for s in &self.skipped {
            println!("SKIP {s}");
        }
        for f in &self.failed {
            println!("FAIL {f}");
        }
        assert!(self.failed.is_empty(), "{} conformance: {} archive(s) failed", self.label, self.failed.len());
    }
}

fn run_corpus(label: &'static str, positive: bool) -> Tally {
    let dir = conformance_dir().join(label);
    let mut tally = Tally { label, passed: 0, failed: Vec::new(), skipped: Vec::new() };
    let files = expected_files(&dir);
    assert!(!files.is_empty(), "no *.expected.json in {}", dir.display());
    for ef in files {
        let expected = read_json(&ef);
        let name = expected.get("archive").and_then(Value::as_str).unwrap_or("").to_string();
        if let Some(reason) = skip_reason(label, &name) {
            tally.skipped.push(format!("{name}: {reason}"));
            continue;
        }
        let fails = match archive_path(&dir, &name) {
            Err(e) => vec![e],
            Ok(p) if !p.exists() => vec!["archive missing".into()],
            Ok(p) => {
                if positive {
                    run_positive_one(&p, &expected)
                } else {
                    run_negative_one(&p, &expected)
                }
            }
        };
        if fails.is_empty() {
            tally.passed += 1;
        } else {
            tally.failed.push(format!("{name}\n     - {}", fails.join("\n     - ")));
        }
    }
    tally
}

#[test]
fn positive_corpus() {
    run_corpus("corpus", true).report();
}

#[test]
fn negative_corpus() {
    run_corpus("negative", false).report();
}

#[test]
fn legacy_corpus() {
    run_corpus("legacy", false).report();
}
