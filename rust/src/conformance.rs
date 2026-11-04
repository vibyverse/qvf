//! The conformance runner, a port of `python/qvf/conformance.py`: certify
//! this implementation against the shared corpus (`conformance/`).
//!
//! * Positive corpus: every archive validates without errors, has the
//!   expected section kinds and source program, and reproduces every check
//!   of its `expected.json` through the reader.
//! * Negative corpus (`negative/`): every archive is judged as its
//!   `expected.json` says, rejected with at least the expected finding
//!   codes, or, for a warning-only archive, valid with the warning. An
//!   archive that names `limits` is validated with those caps.
//! * Legacy corpus (`legacy/`): the negative rules, reported separately.
//!
//! [`run`] prints the report of the Python runner (`PASS`/`FAIL` lines and
//! a summary per corpus) and returns its exit status. This validator always
//! runs the JSON Schema check, so the schema-only codes bind it.

use std::io::{self, Write};
use std::path::{Path, PathBuf};

use serde_json::Value;

use crate::error::Code;
use crate::error::Severity;
use crate::json::json_equal;
use crate::limits::Limits;
use crate::reader::{Reader, Tensor, TensorData};
use crate::validate::{validate_file, ValidateOptions};

/// The judgement of one archive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Outcome {
    /// The archive named by the `expected.json`.
    pub archive: String,
    /// Empty when the archive passed.
    pub failures: Vec<String>,
    /// What a `PASS` line reports: `N checks`, or the expected codes.
    pub note: String,
}

/// The `*.expected.json` files of a directory, sorted.
pub fn expected_files(dir: &Path) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|it| {
            it.filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.ends_with(".expected.json")))
                .collect()
        })
        .unwrap_or_default();
    files.sort();
    files
}

fn read_expected(path: &Path) -> Result<Value, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    serde_json::from_str(&text).map_err(|e| format!("{} is not JSON: {e}", path.display()))
}

/// The archive an expectation names, which must stay inside its directory.
fn archive_path(dir: &Path, name: &str) -> Result<PathBuf, String> {
    let bad = name.is_empty()
        || name.contains('\\')
        || Path::new(name).components().any(|c| !matches!(c, std::path::Component::Normal(_)));
    if bad {
        return Err("corpus archive must stay within the corpus directory".into());
    }
    Ok(dir.join(name))
}

/// Walk a value by keys and indices (negative indices count from the end).
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

fn equal(got: &Value, want: &Value, tol: Option<f64>) -> Result<bool, String> {
    match tol {
        Some(tol) => match (got.as_f64(), want.as_f64()) {
            (Some(g), Some(w)) => Ok((g - w).abs() <= tol),
            _ => Err(format!("cannot compare {got} and {want} within a tolerance")),
        },
        None => Ok(json_equal(got, want)),
    }
}

/// One check of a positive `expected.json`: `Err` holds the failure line.
fn run_check(reader: &Reader, chk: &Value) -> Result<(), String> {
    let tol = chk.get("tol").and_then(Value::as_f64);
    let want = chk.get("equals").cloned().unwrap_or(Value::Null);
    let at = format!(
        "{}/{}",
        chk.get("section").and_then(Value::as_str).unwrap_or("manifest"),
        chk.get("member").and_then(Value::as_str).unwrap_or("")
    );
    let errored = |e: String| format!("check errored {chk}: {e}");
    let got: Value = if let Some(p) = chk.get("manifest_pointer") {
        navigate(reader.manifest(), p).map_err(errored)?.clone()
    } else {
        let section = chk.get("section").and_then(Value::as_str).ok_or_else(|| errored("no section".into()))?;
        let member = chk.get("member").and_then(Value::as_str).ok_or_else(|| errored("no member".into()))?;
        if let Some(p) = chk.get("json_pointer") {
            let doc = reader.read_json(section, member).map_err(|e| errored(e.to_string()))?;
            navigate(&doc, p).map_err(errored)?.clone()
        } else if let Some(p) = chk.get("binary_index") {
            let t = reader.read_tensor(section, member).map_err(|e| errored(e.to_string()))?;
            navigate(&to_list(&t), p).map_err(errored)?.clone()
        } else if chk.get("decode").and_then(Value::as_str) == Some("utf8") {
            let bytes = reader.read_bytes(section, member).map_err(|e| errored(e.to_string()))?;
            let text = String::from_utf8(bytes).map_err(|e| errored(format!("not UTF-8: {e}")))?;
            if let Some(needle) = chk.get("contains").and_then(Value::as_str) {
                return if text.contains(needle) {
                    Ok(())
                } else {
                    Err(format!("check failed ({at}): text does not contain {needle:?}"))
                };
            }
            Value::String(text)
        } else {
            return Err(format!("malformed check (no pointer): {chk}"));
        }
    };
    if equal(&got, &want, tol).map_err(errored)? {
        Ok(())
    } else {
        let tol = tol.filter(|t| *t != 0.0).map(|t| format!(" (tol {t})")).unwrap_or_default();
        Err(format!("check failed ({at}): got {got}, want {want}{tol}"))
    }
}

/// The failures of one positive archive (empty: it passed).
pub fn check_positive(path: &Path, expected: &Value) -> Vec<String> {
    let report = validate_file(path, &ValidateOptions::default());
    if !report.ok {
        let errors: Vec<String> = report
            .findings
            .iter()
            .filter(|f| f.severity == Severity::Error)
            .map(|f| format!("validation failed: {}: {}", f.code, f.message))
            .collect();
        return if errors.is_empty() { vec!["validation failed".into()] } else { errors };
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
            fails.push(e);
        }
    }
    fails
}

/// The failures of one archive of the negative format (`negative/`, `legacy/`).
pub fn check_negative(path: &Path, expected: &Value) -> Vec<String> {
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

/// Judge every archive of a corpus directory: positive rules, or the
/// negative ones. `Err` when an `expected.json` cannot be read or names an
/// archive outside the directory.
pub fn run_dir(dir: &Path, positive: bool) -> Result<Vec<Outcome>, String> {
    let mut out = Vec::new();
    for ef in expected_files(dir) {
        let expected = read_expected(&ef)?;
        let name = expected.get("archive").and_then(Value::as_str).ok_or_else(|| format!("{}: no archive", ef.display()))?;
        let path = archive_path(dir, name)?;
        let note = if positive {
            format!("{} checks", expected.get("checks").and_then(Value::as_array).map_or(0, Vec::len))
        } else {
            let codes: Vec<&str> = expected
                .get("expect")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .filter_map(|w| w.get("code").and_then(Value::as_str))
                .collect();
            codes.join(", ")
        };
        let failures = if !path.exists() {
            vec!["archive missing".into()]
        } else if positive {
            check_positive(&path, &expected)
        } else {
            check_negative(&path, &expected)
        };
        out.push(Outcome { archive: name.to_owned(), failures, note });
    }
    Ok(out)
}

/// `name/` beside the corpus directory, when it holds `expected.json` files.
fn sibling(corpus: &Path, name: &str) -> Option<PathBuf> {
    let absolute = std::fs::canonicalize(corpus).ok().unwrap_or_else(|| corpus.to_path_buf());
    let candidate = absolute.parent()?.join(name);
    (!expected_files(&candidate).is_empty()).then_some(candidate)
}

/// Where [`run`] looks for the corpora.
#[derive(Debug, Clone, Default)]
pub struct RunOptions {
    /// The negative corpus (default: a `negative/` sibling of the corpus).
    pub negative_dir: Option<PathBuf>,
    /// Skip the negative corpus (and, unless `legacy` says otherwise, the legacy one).
    pub no_negative: bool,
    /// The legacy corpus (default: a `legacy/` sibling of the corpus).
    pub legacy_dir: Option<PathBuf>,
    /// Skip the legacy corpus.
    pub no_legacy: bool,
}

fn report_dir(out: &mut dyn Write, outcomes: &[Outcome], positive: bool) -> io::Result<()> {
    for o in outcomes {
        if o.failures.is_empty() {
            writeln!(out, "PASS {}  ({})", o.archive, o.note)?;
        } else if o.failures == ["archive missing"] && positive {
            writeln!(out, "FAIL {}: archive missing", o.archive)?;
        } else {
            writeln!(out, "FAIL {}", o.archive)?;
            for f in &o.failures {
                writeln!(out, "     - {f}")?;
            }
        }
    }
    Ok(())
}

/// Run the positive corpus, then the negative and the legacy ones, and
/// print the report of the Python runner to `out`. Returns the exit status:
/// 0 when every archive is judged as expected, 1 when one is not, 2 when a
/// corpus has no `expected.json`. `Err` for an unreadable expectation.
pub fn run(corpus: &Path, options: &RunOptions, out: &mut dyn Write) -> Result<u8, String> {
    let io = |e: io::Error| e.to_string();
    let mut rc = 0u8;
    let corpus_run = |dir: &Path, label: &str, positive: bool, out: &mut dyn Write| -> Result<u8, String> {
        if expected_files(dir).is_empty() {
            eprintln!("no *.expected.json found in {}", dir.display());
            return Ok(2);
        }
        let outcomes = run_dir(dir, positive)?;
        report_dir(out, &outcomes, positive).map_err(io)?;
        let passed = outcomes.iter().filter(|o| o.failures.is_empty()).count();
        let summary = match label {
            "corpus" => format!("conformance: {passed}/{} archives passed", outcomes.len()),
            "negative" => format!("negative conformance: {passed}/{} archives rejected as expected", outcomes.len()),
            _ => format!("{label} conformance: {passed}/{} archives judged as expected", outcomes.len()),
        };
        writeln!(out, "\n{summary}").map_err(io)?;
        Ok(u8::from(passed != outcomes.len()))
    };
    rc = rc.max(corpus_run(corpus, "corpus", true, out)?);
    let combine = |rc: u8, next: u8| if rc != 0 { rc } else { next };
    if !options.no_negative {
        if let Some(dir) = options.negative_dir.clone().or_else(|| sibling(corpus, "negative")) {
            writeln!(out).map_err(io)?;
            let r = corpus_run(&dir, "negative", false, out)?;
            rc = combine(rc, r);
        }
    }
    let legacy = if options.no_legacy {
        false
    } else if options.legacy_dir.is_some() {
        true
    } else {
        !options.no_negative
    };
    if legacy {
        if let Some(dir) = options.legacy_dir.clone().or_else(|| sibling(corpus, "legacy")) {
            writeln!(out).map_err(io)?;
            let r = corpus_run(&dir, "legacy", false, out)?;
            rc = combine(rc, r);
        }
    }
    Ok(rc)
}
