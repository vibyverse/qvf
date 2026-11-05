//! Helpers shared by the builder and editor tests: the Rust and Python
//! validators on an archive file.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::Command;

use qvf::{validate_bytes, ValidateOptions};
use serde_json::Value;

/// A path in the test's scratch directory (under `target/`).
pub fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("qvf-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// `code location message` of every Rust finding.
pub fn rust_findings(bytes: &[u8]) -> Vec<String> {
    let report = validate_bytes(bytes.to_vec(), &ValidateOptions::default());
    report.findings.iter().map(|f| format!("{} {} {}", f.code, f.location, f.message)).collect()
}

/// The interpreter that runs the Python package (`$PYTHON` or `python3`
/// with `PYTHONPATH=../python`), or `None` when it cannot import `qvf`.
pub fn python() -> Option<Command> {
    let python_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("python");
    let interpreter = std::env::var("PYTHON").unwrap_or_else(|_| "python3".into());
    let probe = Command::new(&interpreter).env("PYTHONPATH", &python_dir).args(["-c", "import qvf"]).output().ok()?;
    if !probe.status.success() {
        return None;
    }
    let mut c = Command::new(interpreter);
    c.env("PYTHONPATH", python_dir);
    Some(c)
}

/// The Python validator's findings as `(severity, code, where)`, ignoring
/// `W_SCHEMA_SKIPPED` (no `jsonschema`), or `None` without Python.
pub fn python_findings(path: &Path) -> Option<Vec<(String, String, String)>> {
    let out = python()?.args(["-m", "qvf", "validate", "--json"]).arg(path).output().expect("python runs");
    let report: Value = serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("python validator printed no JSON ({e}): {}", String::from_utf8_lossy(&out.stderr)));
    let s = |f: &Value, k: &str| f[k].as_str().unwrap_or_default().to_owned();
    Some(
        report["findings"]
            .as_array()
            .expect("findings")
            .iter()
            .map(|f| (s(f, "severity"), s(f, "code"), s(f, "where")))
            .filter(|(_, c, _)| c != "W_SCHEMA_SKIPPED")
            .collect(),
    )
}

/// Assert that both validators report nothing for the archive at `path`.
pub fn assert_clean(path: &Path) {
    let bytes = std::fs::read(path).unwrap();
    assert_eq!(rust_findings(&bytes), Vec::<String>::new(), "Rust validator findings for {}", path.display());
    match python_findings(path) {
        Some(f) => assert_eq!(f, Vec::new(), "Python validator findings for {}", path.display()),
        None => eprintln!("SKIP Python validation of {}: python3 with the qvf package is not available", path.display()),
    }
}
