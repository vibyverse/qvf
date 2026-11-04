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

use serde_json::{json, Value};

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
        except Exception as exc:  # a defect of the Python CLI: reported as code -1
            code = -1
            print(f"uncaught {type(exc).__name__}: {exc}", file=sys.stderr)
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
    // One directory per call: the tests run in parallel threads.
    static CALLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let call = CALLS.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    let dir = std::env::temp_dir().join(format!("qvf-rust-cli-{}-{call}", std::process::id()));
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

/// The members of an archive as (role path, declared spec without the
/// path), and its manifest without the time stamp and UUID scrub renews.
fn content(path: &Path) -> (Vec<String>, Value) {
    let r = qvf::Reader::open(path).unwrap();
    let mut names: Vec<String> = r.zip().entries().iter().map(|e| e.name.clone()).collect();
    names.sort();
    let mut m = r.manifest().clone();
    if let Some(meta) = m.get_mut("metadata").and_then(Value::as_object_mut) {
        for key in ["uuid", "created_utc"] {
            meta.remove(key);
        }
    }
    (names, m)
}

#[test]
fn maintenance_commands_match_python_on_the_corpus() {
    let archives = archives();
    let base = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("qvf-cli-maint-{}", std::process::id()));
    let (py_dir, rs_dir) = (base.join("python"), base.join("rust"));
    std::fs::create_dir_all(&py_dir).unwrap();
    std::fs::create_dir_all(&rs_dir).unwrap();
    let s = |p: &Path| p.to_string_lossy().into_owned();
    // (command, output file name or "")
    let mut commands: Vec<(Vec<String>, String)> = vec![(vec!["schema".into()], String::new())];
    for (k, path) in archives.iter().enumerate() {
        let p = s(path);
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        let next = s(&archives[(k + 1) % archives.len()]);
        commands.push((vec!["diff".into(), "--json".into(), p.clone(), p.clone()], String::new()));
        commands.push((
            vec!["diff".into(), "--json".into(), "--max-differences".into(), "1000".into(), p.clone(), next.clone()],
            String::new(),
        ));
        commands.push((
            vec!["diff".into(), "--ignore-metadata".into(), "--atol".into(), "1e-6".into(), p.clone(), next],
            String::new(),
        ));
        commands.push((vec!["repack".into(), p.clone(), "-o".into(), format!("@/repack_{name}")], format!("repack_{name}")));
        commands.push((
            vec!["repack".into(), "--compression".into(), "store".into(), p.clone(), "-o".into(), format!("@/store_{name}")],
            format!("store_{name}"),
        ));
        commands.push((vec!["scrub".into(), "--dry-run".into(), p.clone()], String::new()));
        commands.push((vec!["scrub".into(), p.clone(), "-o".into(), format!("@/scrub_{name}")], format!("scrub_{name}")));
    }
    // A tolerance difference: the same archive against a copy with perturbed floats is covered by
    // tests/editor.rs; here the corpus pairs exercise the manifest and payload walk.
    let with_dir = |c: &[String], dir: &Path| -> Vec<String> { c.iter().map(|a| a.replace('@', &s(dir))).collect() };
    let py_commands: Vec<Vec<String>> = commands.iter().map(|(c, _)| with_dir(c, &py_dir)).collect();
    let Some(expected) = python(&py_commands) else {
        println!("SKIP: Python or the qvf package is unavailable");
        return;
    };
    let mut failures = Vec::new();
    let mut counts: std::collections::BTreeMap<String, usize> = Default::default();
    for (((command, output), py_command), (py_code, py_out)) in commands.iter().zip(&py_commands).zip(&expected) {
        let rs_command = with_dir(command, &rs_dir);
        let what = format!("qvf {}", py_command.join(" "));
        let (code, out) = rust(&rs_command);
        *counts.entry(command[0].clone()).or_default() += 1;
        if *py_code == -1 {
            // The Python CLI raised an uncaught exception (a defect of the
            // reference CLI on malformed manifests): there is nothing to match.
            *counts.entry("python_uncaught".into()).or_default() += 1;
            continue;
        }
        if code != *py_code {
            failures.push(format!("{what}: exit {code}, Python {py_code}"));
            continue;
        }
        match command[0].as_str() {
            "schema" if out != *py_out => failures.push(format!("{what}: schema text differs")),
            "diff" if code != 2 && command[1] == "--json" && json(&out) != json(py_out) => {
                failures.push(format!("{what}: {out} Python {py_out}"))
            }
            "diff" if command[1] != "--json" && out != *py_out => failures.push(format!("{what}: {out:?} Python {py_out:?}")),
            "scrub" if code == 0 && json(&out) != json(py_out) => failures.push(format!("{what}: {out} Python {py_out}")),
            _ => {}
        }
        if !output.is_empty() && code == 0 {
            let (rs, py) = (rs_dir.join(output), py_dir.join(output));
            if command[0] == "repack" && out.replace(&s(&rs_dir), "@") != py_out.replace(&s(&py_dir), "@") {
                failures.push(format!("{what}: stdout {out:?}, Python {py_out:?}"));
            }
            let (rc, pc) = (content(&rs), content(&py));
            // The Python model leaves a `sections` array that holds a
            // non-object untyped, so its repack copies no member at all and
            // writes an archive with E_MEMBER_MISSING. The Rust repack keeps
            // the referenced members.
            let untyped = qvf::Reader::open(&rs).unwrap().manifest()["sections"]
                .as_array()
                .is_some_and(|l| l.iter().any(|s| !s.is_object()));
            if untyped {
                *counts.entry("python_untyped_sections".into()).or_default() += 1;
                continue;
            }
            if rc != pc {
                failures.push(format!("{what}: archive content differs: {rc:?} Python {pc:?}"));
            }
            let (rv, pv) =
                (rust(&["validate".into(), "--json".into(), s(&rs)]), rust(&["validate".into(), "--json".into(), s(&py)]));
            let (rf, pf) = (findings(&json(&rv.1)).0, findings(&json(&pv.1)).0);
            if rf != pf {
                failures.push(format!("{what}: findings of the output {rf:?}, of Python's {pf:?}"));
            }
            if command[0] == "scrub" {
                let (r, p) = (qvf::Reader::open(&rs).unwrap(), qvf::Reader::open(&py).unwrap());
                let (rm, pm) = (&r.manifest()["metadata"], &p.manifest()["metadata"]);
                let uuid = rm["uuid"].as_str().unwrap_or_default();
                if uuid.len() != 36 || &uuid[14..15] != "4" || rm["uuid"] == pm["uuid"] {
                    failures.push(format!("{what}: uuid {} (Python {})", rm["uuid"], pm["uuid"]));
                }
                let created = rm["created_utc"].as_str().unwrap_or_default();
                if created.len() != pm["created_utc"].as_str().unwrap_or_default().len() || !created.ends_with('Z') {
                    failures.push(format!("{what}: created_utc {created} (Python {})", pm["created_utc"]));
                }
            }
            if command[1] == "--compression" {
                let r = qvf::Reader::open(&rs).unwrap();
                assert!(r.zip().entries().iter().all(|e| e.method == qvf::zip::METHOD_STORE), "{what}");
            }
            if command[0] == "repack" || command[0] == "scrub" {
                let again = rust(&rs_command);
                if again.0 != 2 {
                    failures.push(format!("{what}: an existing output was not refused"));
                }
            }
        }
    }
    let _ = std::fs::remove_dir_all(&base);
    println!("qvf maintenance commands against Python: {counts:?}, {} differences", failures.len());
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn the_schema_copy_is_the_normative_schema() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let copy = std::fs::read(root.join("qvf_manifest.schema.json")).unwrap();
    assert_eq!(
        copy,
        std::fs::read(root.join("../spec/qvf_manifest.schema.json")).unwrap(),
        "rust/qvf_manifest.schema.json is stale"
    );
    assert_eq!(qvf::reader::sha256_hex(&copy), qvf::generated::SCHEMA_SHA256);
}

#[test]
fn diff_tolerances_and_truncation_match_python() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("qvf-cli-diff-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let a = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance/corpus/tddft.qvf");
    let b = dir.join("perturbed.qvf");
    std::fs::copy(&a, &b).unwrap();
    // Perturb floats by about 1e-9 absolute (tensor) and 1e-8 relative (JSON).
    let r = qvf::Reader::open(&a).unwrap();
    let mo = r.read_tensor("wf", "mo_coefficients").unwrap();
    let mut structure = r.read_json("structure", "structure").unwrap();
    drop(r);
    let c: Vec<f64> = mo.to_f64_vec().iter().enumerate().map(|(i, x)| if i % 5 == 0 { x + 1e-9 } else { *x }).collect();
    for atom in structure["atoms"].as_array_mut().unwrap() {
        let x = atom["position"][0].as_f64().unwrap();
        atom["position"][0] = json!(x * (1.0 + 1e-8) + 1e-12);
    }
    let mut ed = qvf::Editor::open(&b).unwrap();
    ed.replace_tensor("wf", "mo_coefficients", &c, &mo.shape).unwrap();
    ed.replace_json("structure", "structure", &structure).unwrap();
    ed.commit().unwrap();
    let (a, b) = (a.to_string_lossy().into_owned(), b.to_string_lossy().into_owned());
    let mut commands: Vec<Vec<String>> = Vec::new();
    for opts in [
        vec!["--json"],
        vec!["--json", "--ignore-metadata"],
        vec!["--json", "--ignore-metadata", "--atol", "1e-6"],
        vec!["--json", "--ignore-metadata", "--rtol", "1e-7"],
        vec!["--json", "--ignore-metadata", "--rtol", "1e-9"],
        vec!["--json", "--ignore-metadata", "--max-differences", "3"],
        vec!["--ignore-metadata", "--max-differences", "2"],
        vec!["--ignore-metadata", "--atol", "1e-6"],
        vec!["--atol", "-1"],
        vec!["--atol", "nan"],
        vec!["--max-differences", "0"],
        vec!["--atol", "x"],
    ] {
        let mut c: Vec<String> = vec!["diff".into()];
        c.extend(opts.iter().map(|s| s.to_string()));
        c.push(a.clone());
        c.push(b.clone());
        commands.push(c);
    }
    let Some(expected) = python(&commands) else {
        println!("SKIP: Python or the qvf package is unavailable");
        return;
    };
    let mut failures = Vec::new();
    for (command, (py_code, py_out)) in commands.iter().zip(&expected) {
        let (code, out) = rust(command);
        let what = format!("qvf {}", command[..command.len() - 2].join(" "));
        if code != *py_code {
            failures.push(format!("{what}: exit {code}, Python {py_code}"));
        } else if code != 2 && command[1] == "--json" && json(&out) != json(py_out) {
            failures.push(format!("{what}: {out} Python {py_out}"));
        } else if code != 2 && command[1] != "--json" && out != *py_out {
            failures.push(format!("{what}: {out:?} Python {py_out:?}"));
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    let (strict, loose) = (json(&expected[1].1), json(&expected[2].1));
    assert!(strict["difference_count"].as_u64().unwrap() > 10, "{strict}");
    // Within tolerance; the update writes the default byte_order explicitly, as the Python editor does.
    assert_eq!(loose["differences"], json!([{"where": "manifest/byte_order", "message": "missing from left"}]), "{loose}");
}

#[test]
fn conformance_matches_python() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../conformance");
    let dir = |d: &str| root.join(d).to_string_lossy().into_owned();
    let empty = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("qvf-cli-empty-{}", std::process::id()));
    std::fs::create_dir_all(&empty).unwrap();
    let commands: Vec<Vec<String>> = vec![
        vec!["conformance".into(), dir("corpus")],
        vec!["conformance".into(), "--no-negative".into(), dir("corpus")],
        vec!["conformance".into(), "--no-legacy".into(), dir("corpus")],
        vec!["conformance".into(), dir("corpus"), "--negative".into(), dir("legacy"), "--no-legacy".into()],
        vec!["conformance".into(), "--no-negative".into(), "--legacy".into(), dir("legacy"), dir("corpus")],
        // The negative corpus run with the positive rules: every archive fails.
        vec!["conformance".into(), "--no-negative".into(), dir("negative")],
        vec!["conformance".into(), empty.to_string_lossy().into_owned()],
    ];
    let Some(expected) = python(&commands) else {
        println!("SKIP: Python or the qvf package is unavailable");
        return;
    };
    // The PASS/FAIL lines and the summaries; the failure details are free text.
    let lines = |t: &str| -> Vec<String> { t.lines().filter(|l| !l.starts_with("     - ")).map(str::to_owned).collect() };
    let mut failures = Vec::new();
    for (command, (py_code, py_out)) in commands.iter().zip(&expected) {
        let (code, out) = rust(command);
        let what = format!("qvf {}", command.join(" "));
        if code != *py_code {
            failures.push(format!("{what}: exit {code}, Python {py_code}"));
        } else if lines(&out) != lines(py_out) {
            failures.push(format!("{what}:\n{out}\nPython:\n{py_out}"));
        }
    }
    let _ = std::fs::remove_dir_all(&empty);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(expected[0].0, 0);
    assert_eq!(expected[5].0, 1);
}
