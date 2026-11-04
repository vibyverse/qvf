//! The `qvf` binary against the Python reference CLI (`python -m qvf`) on
//! every archive of the positive, negative and legacy corpora: `validate`
//! (at every level) must give the same exit code, verdict, finding codes and
//! (code, severity, location) findings, and the same validation context, and
//! the human output the same verdict line and finding lines;
//! `info --json` and `cat` the same exit code and JSON. Skipped with a `SKIP`
//! line when Python or the `qvf` package is unavailable.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

const LEVELS: [&str; 6] = ["container", "manifest", "payloads", "semantics", "references", "all"];

/// Runs `qvf.cli.main` in one interpreter for every command line given as a
/// JSON list on stdin; prints one `{"code", "stdout"}` line per command.
const DRIVER: &str = r#"
import contextlib, io, json, sys
from qvf import cli
for line in sys.stdin:
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = cli.main(json.loads(line))
        except SystemExit as exc:
            code = exc.code
    print(json.dumps({"code": code, "stdout": out.getvalue()}), flush=True)
"#;

fn archives() -> Vec<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance");
    let mut out = Vec::new();
    for dir in ["corpus", "negative", "legacy"] {
        let mut v: Vec<PathBuf> = std::fs::read_dir(root.join(dir))
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.extension().is_some_and(|x| x == "qvf"))
            .collect();
        v.sort();
        out.extend(v);
    }
    out
}

/// Exit code and stdout of each command under the Python CLI, or `None`
/// when Python or the package is unavailable.
fn python(commands: &[Vec<String>]) -> Option<Vec<(i64, String)>> {
    let python_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../python");
    let interpreter = std::env::var("PYTHON").unwrap_or_else(|_| "python3".into());
    let probe = Command::new(&interpreter).env("PYTHONPATH", &python_dir).args(["-c", "import qvf.cli"]).output().ok()?;
    if !probe.status.success() {
        return None;
    }
    let dir = std::env::temp_dir().join(format!("qvf-rust-cli-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("commands.jsonl");
    let lines: Vec<String> = commands.iter().map(|c| serde_json::to_string(c).unwrap()).collect();
    std::fs::write(&input, lines.join("\n") + "\n").unwrap();
    let out = Command::new(&interpreter)
        .env("PYTHONPATH", &python_dir)
        .args(["-c", DRIVER])
        .stdin(std::fs::File::open(&input).unwrap())
        .output()
        .expect("python runs");
    let _ = std::fs::remove_dir_all(&dir);
    assert!(out.status.success(), "python driver failed: {}", String::from_utf8_lossy(&out.stderr));
    let results: Vec<(i64, String)> = String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|l| {
            let v: Value = serde_json::from_str(l).unwrap();
            (v["code"].as_i64().unwrap_or(2), v["stdout"].as_str().unwrap().to_owned())
        })
        .collect();
    assert_eq!(results.len(), commands.len());
    Some(results)
}

fn rust(command: &[String]) -> (i64, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_qvf")).args(command).output().expect("qvf runs");
    (out.status.code().unwrap_or(-1) as i64, String::from_utf8(out.stdout).unwrap())
}

fn json(text: &str) -> Value {
    serde_json::from_str(text).unwrap_or_else(|e| panic!("not JSON ({e}): {text}"))
}

/// The findings as (code, severity, where), and the codes alone.
fn findings(report: &Value) -> (BTreeSet<(String, String, String)>, BTreeSet<String>) {
    let list = report["findings"].as_array().unwrap();
    let s = |f: &Value, k: &str| f[k].as_str().unwrap().to_owned();
    (list.iter().map(|f| (s(f, "code"), s(f, "severity"), s(f, "where"))).collect(), list.iter().map(|f| s(f, "code")).collect())
}

/// `info --json` with each section's members in role order (Python keeps
/// manifest order, serde_json sorts keys).
fn info(mut v: Value) -> Value {
    for s in v["sections"].as_array_mut().into_iter().flatten() {
        if let Some(m) = s["members"].as_array_mut() {
            m.sort_by(|a, b| a["role"].as_str().cmp(&b["role"].as_str()));
        }
    }
    v
}

#[test]
fn cli_matches_python_on_the_corpus() {
    let archives = archives();
    let mut commands: Vec<Vec<String>> = Vec::new();
    for path in &archives {
        let p = path.to_string_lossy().into_owned();
        for level in LEVELS {
            commands.push(vec!["validate".into(), "--json".into(), "--level".into(), level.into(), p.clone()]);
        }
        commands.push(vec!["validate".into(), p.clone()]);
        commands.push(vec!["info".into(), "--json".into(), p.clone()]);
        commands.push(vec!["cat".into(), p.clone()]);
    }
    let Some(expected) = python(&commands) else {
        println!("SKIP: Python or the qvf package is unavailable");
        return;
    };
    let mut failures = Vec::new();
    for (command, (py_code, py_out)) in commands.iter().zip(&expected) {
        let name = Path::new(command.last().unwrap()).file_name().unwrap().to_string_lossy().into_owned();
        let what = format!("{name}: qvf {}", command[..command.len() - 1].join(" "));
        let (code, out) = rust(command);
        if code != *py_code {
            failures.push(format!("{what}: exit {code}, Python {py_code}"));
            continue;
        }
        match command[0].as_str() {
            "validate" if command.len() == 2 => {
                // Human output: the verdict line exactly, then one
                // "severity code where: message" line per finding.
                let lines = |t: &str| -> (Option<String>, BTreeSet<String>) {
                    let mut it = t.lines();
                    let head = it.next().map(str::to_owned);
                    let rest = it.skip(1).map(|l| l.split_once(": ").map_or(l, |(a, _)| a).to_owned()).collect();
                    (head, rest)
                };
                if lines(&out) != lines(py_out) {
                    failures.push(format!("{what}: human output {:?}, Python {:?}", lines(&out), lines(py_out)));
                }
            }
            "validate" => {
                let (r, p) = (json(&out), json(py_out));
                let ((rf, rc), (pf, pc)) = (findings(&r), findings(&p));
                if r["ok"] != p["ok"] || rc != pc {
                    failures.push(format!("{what}: codes {rc:?}, Python {pc:?}"));
                } else if rf != pf {
                    failures.push(format!("{what}: findings {rf:?}, Python {pf:?}"));
                }
                for key in
                    ["archive_format_version", "archive_schema_uri", "schema_id", "schema_sha256", "unsupported_format_version"]
                {
                    let (a, b) = (&r["validation_context"][key], &p["validation_context"][key]);
                    if a != b {
                        failures.push(format!("{what}: validation_context.{key} {a}, Python {b}"));
                    }
                }
            }
            "info" if code == 0 && info(json(&out)) != info(json(py_out)) => {
                failures.push(format!("{what}: output differs"));
            }
            "cat" if code == 0 && json(&out) != json(py_out) => {
                failures.push(format!("{what}: output differs"));
            }
            _ => {}
        }
    }
    println!("qvf CLI against Python: {} archives, {} commands, {} differences", archives.len(), commands.len(), failures.len());
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn extract_verifies_and_never_overwrites() {
    let archive = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/corpus/tddft.qvf");
    let dir = std::env::temp_dir().join(format!("qvf-rust-extract-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("mo.dat");
    let run = |args: &[&str]| Command::new(env!("CARGO_BIN_EXE_qvf")).args(args).output().unwrap();
    let a = archive.to_str().unwrap();
    let o = out.to_str().unwrap();
    let first = run(&["extract", a, "wf", "mo_coefficients", "-o", o]);
    assert_eq!(first.status.code(), Some(0), "{}", String::from_utf8_lossy(&first.stderr));
    let reader = qvf::Reader::open(&archive).unwrap();
    assert_eq!(std::fs::read(&out).unwrap(), reader.read_bytes("wf", "mo_coefficients").unwrap());
    let again = run(&["extract", a, "wf", "mo_coefficients", "-o", o]);
    assert_eq!(again.status.code(), Some(2), "an existing output is never overwritten");
    let missing = run(&["extract", a, "wf", "no_such_member", "-o", dir.join("x").to_str().unwrap()]);
    assert_eq!(missing.status.code(), Some(2));
    assert!(!dir.join("x").exists());
    let tampered = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/negative/sha256_mismatch.qvf");
    if tampered.exists() {
        let r = qvf::Reader::open(&tampered).unwrap();
        let s = r.sections()[0];
        let role = s.members().keys().next().unwrap().clone();
        let t = run(&["extract", tampered.to_str().unwrap(), s.id, &role, "-o", dir.join("t").to_str().unwrap()]);
        assert_eq!(t.status.code(), Some(2));
        assert!(!dir.join("t").exists());
    }
    assert_eq!(run(&["validate"]).status.code(), Some(2));
    assert_eq!(run(&["validate", a, "--level", "bogus"]).status.code(), Some(2));
    assert_eq!(run(&["frobnicate"]).status.code(), Some(2));
    let _ = std::fs::remove_dir_all(&dir);
}
