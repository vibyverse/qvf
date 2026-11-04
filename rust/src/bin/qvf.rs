//! `qvf`: validate and inspect QVF archives, the `validate`, `info`, `cat`
//! and `extract` commands of the Python `qvf` CLI with the same finding
//! codes, report layout and exit codes: 0 success, 1 validation failed,
//! 2 usage or I/O error.
//!
//! ```sh
//! qvf validate ARCHIVE [--level LEVEL] [--json]
//! qvf info ARCHIVE [--json]
//! qvf cat ARCHIVE [SECTION MEMBER]
//! qvf extract ARCHIVE SECTION MEMBER -o OUTPUT
//! ```

use std::fs::OpenOptions;
use std::io::Write;
use std::process::ExitCode;

use qvf::zip::ZipArchive;
use qvf::{generated, validate_bytes, Level, Limits, Reader, Report, ValidateOptions};
use serde_json::{json, Map, Value};

const USAGE: &str = "usage: qvf [--version] {validate,info,cat,extract} ...

  validate ARCHIVE [--level {container,manifest,payloads,semantics,references,all}] [--json]
                                validate cumulatively through a selected stage
  info ARCHIVE [--json]         show sections, kinds and member sizes
  cat ARCHIVE [SECTION MEMBER]  print the manifest or one JSON member
  extract ARCHIVE SECTION MEMBER -o OUTPUT
                                verify and extract one declared member (never overwrites)

exit status: 0 success, 1 validation failed, 2 usage or I/O error";

/// A failure that ends the command with exit status 2.
struct Fail(String);

impl<E: std::fmt::Display> From<E> for Fail {
    fn from(e: E) -> Self {
        Fail(e.to_string())
    }
}

fn usage(message: &str) -> Fail {
    Fail(format!("{message}\n{USAGE}"))
}

/// Positional arguments and flags of one command.
struct Args {
    positional: Vec<String>,
    json: bool,
    level: Option<String>,
    output: Option<String>,
}

fn parse(args: &[String], allowed: &[&str]) -> Result<Args, Fail> {
    let mut out = Args { positional: Vec::new(), json: false, level: None, output: None };
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let (flag, inline) = match a.split_once('=') {
            Some((f, v)) if a.starts_with("--") => (f, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let mut value = |name: &str| {
            inline.clone().or_else(|| it.next().cloned()).ok_or_else(|| usage(&format!("argument {name}: expected one argument")))
        };
        match flag {
            "--json" if allowed.contains(&"--json") => out.json = true,
            "--level" if allowed.contains(&"--level") => out.level = Some(value("--level")?),
            "-o" | "--output" if allowed.contains(&"-o") => out.output = Some(value("-o/--output")?),
            "-h" | "--help" => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            _ if a.starts_with('-') && a.len() > 1 => return Err(usage(&format!("unrecognized arguments: {a}"))),
            _ => out.positional.push(a.clone()),
        }
    }
    Ok(out)
}

fn level(name: &str) -> Option<Level> {
    Some(match name {
        "container" => Level::Container,
        "manifest" => Level::Manifest,
        "payloads" => Level::Payloads,
        "semantics" => Level::Semantics,
        "references" => Level::References,
        "all" => Level::All,
        _ => return None,
    })
}

fn print_json(value: &Value) {
    println!("{}", serde_json::to_string_pretty(value).expect("a Value serializes"));
}

/// A JSON value as Python's `repr` prints the parsed value, for the human
/// output the Python CLI formats with f-strings.
fn py_repr(v: &Value) -> String {
    match v {
        Value::Null => "None".into(),
        Value::Bool(b) => if *b { "True" } else { "False" }.into(),
        Value::Number(n) => match (n.as_i64(), n.as_u64(), n.as_f64()) {
            (Some(i), _, _) => i.to_string(),
            (_, Some(u), _) => u.to_string(),
            (_, _, Some(f)) if f.fract() == 0.0 && f.abs() < 1e16 => format!("{f:.1}"),
            (_, _, Some(f)) => f.to_string(),
            _ => n.to_string(),
        },
        Value::String(s) => {
            let quote = if s.contains('\'') && !s.contains('"') { '"' } else { '\'' };
            let mut out = String::from(quote);
            for c in s.chars() {
                match c {
                    '\\' => out.push_str("\\\\"),
                    '\n' => out.push_str("\\n"),
                    '\r' => out.push_str("\\r"),
                    '\t' => out.push_str("\\t"),
                    c if c == quote => {
                        out.push('\\');
                        out.push(c);
                    }
                    c if (c as u32) < 0x20 || c as u32 == 0x7f => out.push_str(&format!("\\x{:02x}", c as u32)),
                    c => out.push(c),
                }
            }
            out.push(quote);
            out
        }
        Value::Array(a) => format!("[{}]", a.iter().map(py_repr).collect::<Vec<_>>().join(", ")),
        Value::Object(o) => {
            // Python keeps file order; this keeps the `source` order of
            // spec § 3.2 and sorts the rest.
            let first = ["program", "version", "calculation"];
            let keys = first.iter().filter(|k| o.contains_key(**k)).map(|k| k.to_string());
            let keys: Vec<String> = keys.chain(o.keys().filter(|k| !first.contains(&k.as_str())).cloned()).collect();
            let items: Vec<String> =
                keys.iter().map(|k| format!("{}: {}", py_repr(&Value::from(k.as_str())), py_repr(&o[k]))).collect();
            format!("{{{}}}", items.join(", "))
        }
    }
}

/// The `validation_context` block of the Python report, filled from the
/// archive's manifest when the level reads it.
fn context(bytes: &[u8], level: Level, limits: &Limits) -> Value {
    let manifest = (level >= Level::Manifest)
        .then(|| {
            let zip = ZipArchive::open(bytes.to_vec()).ok()?;
            let entry = zip.find("manifest.json")?;
            // The caps under which the validators parse it at all.
            let too_large = entry.size as f64 > limits.max_json_bytes.min(limits.max_member_bytes);
            if too_large || qvf::zip::ratio_exceeded(entry, limits) {
                return None;
            }
            let raw = zip.read(entry, Some(limits.max_member_bytes)).ok()?;
            if qvf::json::json_depth(&raw, limits.max_json_depth) > limits.max_json_depth {
                return None;
            }
            serde_json::from_slice::<Value>(&raw).ok().filter(Value::is_object)
        })
        .flatten();
    let version = manifest.as_ref().and_then(|m| m.get("qvf_version")).cloned().unwrap_or(Value::Null);
    let unsupported = manifest
        .as_ref()
        .map(|_| Value::Bool(!matches!(version.as_f64(), Some(v) if !version.is_boolean() && (v == 1.0 || v == 2.0))));
    json!({
        "reader": "qvf-rust",
        "reader_version": env!("CARGO_PKG_VERSION"),
        "archive_format_version": version,
        "archive_schema_uri": manifest.as_ref().and_then(|m| m.get("schema_uri")).cloned().unwrap_or(Value::Null),
        "schema_toolkit_version": env!("CARGO_PKG_VERSION"),
        "schema_id": generated::SCHEMA_ID,
        "schema_sha256": generated::SCHEMA_SHA256,
        "schema_checked": manifest.is_some(),
        "supported_format_versions": generated::READABLE_VERSIONS,
        "unsupported_format_version": unsupported.unwrap_or(Value::Null),
        "scientific_validity": "not_assessed",
    })
}

fn report_json(report: &Report, context: Value) -> Value {
    let lines = |severity: qvf::Severity| -> Vec<String> {
        report.findings.iter().filter(|f| f.severity == severity).map(|f| format!("{}: {}", f.code.as_str(), f.message)).collect()
    };
    let findings: Vec<Value> = report
        .findings
        .iter()
        .map(|f| json!({"severity": f.severity.as_str(), "code": f.code.as_str(), "where": f.location, "message": f.message}))
        .collect();
    json!({
        "ok": report.ok,
        "errors": lines(qvf::Severity::Error),
        "warnings": lines(qvf::Severity::Warning),
        "n_sections": report.n_sections,
        "findings": findings,
        "validation_context": context,
    })
}

fn validate(args: &[String]) -> Result<u8, Fail> {
    let a = parse(args, &["--json", "--level"])?;
    let [archive] = a.positional.as_slice() else {
        return Err(usage("validate requires exactly one ARCHIVE"));
    };
    let level_name = a.level.unwrap_or_else(|| "all".into());
    let level = level(&level_name).ok_or_else(|| usage(&format!("argument --level: invalid choice: {level_name:?}")))?;
    let options = ValidateOptions { level, ..ValidateOptions::default() };
    let (report, context) = match std::fs::read(archive) {
        Ok(bytes) => {
            let context = context(&bytes, level, &options.limits);
            (validate_bytes(bytes, &options), context)
        }
        // validate_file reports the unreadable file as E_ZIP_UNREADABLE.
        Err(_) => (qvf::validate_file(archive, &options), context(&[], Level::Container, &options.limits)),
    };
    let value = report_json(&report, context);
    if a.json {
        print_json(&value);
    } else {
        let c = &value["validation_context"];
        println!("{level_name}: {} ({} sections)", if report.ok { "OK" } else { "FAILED" }, report.n_sections);
        println!(
            "archive format={}; reader={} {}; bundled schema={} sha256={}; schema checked={}. Format version is not the \
             toolkit release: compare the producer's validator pin before attributing a cross-version rejection to the \
             producer. Numerical accuracy and electronic stability are not assessed.",
            py_repr(&c["archive_format_version"]),
            c["reader"].as_str().unwrap_or_default(),
            c["reader_version"].as_str().unwrap_or_default(),
            c["schema_id"].as_str().unwrap_or_default(),
            c["schema_sha256"].as_str().unwrap_or_default(),
            py_repr(&c["schema_checked"]),
        );
        for f in &report.findings {
            println!("{} {} {}: {}", f.severity.as_str(), f.code.as_str(), f.location, f.message);
        }
    }
    Ok(if report.ok { 0 } else { 1 })
}

/// The member object as the Python model writes it, plus its ZIP sizes.
fn member_info(reader: &Reader, role: &str, spec: &Map<String, Value>) -> Result<Value, Fail> {
    let mut out = Map::new();
    for key in ["path", "format", "dtype", "shape", "sha256", "size", "media_type"] {
        if let Some(v) = spec.get(key) {
            let v = if key == "format" { Value::from(if v == "json" { "json" } else { "binary" }) } else { v.clone() };
            out.insert(key.into(), v);
        }
    }
    let path = spec.get("path").and_then(Value::as_str).ok_or_else(|| Fail(format!("member {role:?} has no path")))?;
    let entry = reader.zip().find(path).ok_or_else(|| Fail(format!("\"There is no item named '{path}' in the archive\"")))?;
    out.insert("role".into(), Value::from(role));
    out.insert("size".into(), Value::from(entry.size));
    out.insert("compressed_size".into(), Value::from(entry.compressed_size));
    Ok(Value::Object(out))
}

/// Open an archive with the checks the Python `Reader` makes at open: on
/// top of the Rust [`Reader`]'s, every entry path is legal, every entry is
/// within the size and ratio caps, and `qvf_version` is 1 (or the withdrawn 2).
fn open_checked(path: &str) -> Result<Reader, Fail> {
    let reader = Reader::open(path)?;
    let limits = reader.limits();
    for e in reader.zip().entries() {
        if let Some(why) = qvf::zip::path_problem(&e.name) {
            return Err(Fail(format!("E_MEMBER_PATH: {:?}: {why}", e.name)));
        }
        if e.size as f64 > limits.max_member_bytes {
            return Err(Fail(format!("E_MEMBER_TOO_LARGE: {:?}", e.name)));
        }
        if qvf::zip::ratio_exceeded(e, limits) {
            return Err(Fail(format!("E_COMPRESSION_RATIO: {:?}", e.name)));
        }
    }
    if reader.zip().total_size() as f64 > limits.max_total_bytes {
        return Err(Fail("E_TOTAL_TOO_LARGE: archive exceeds total expansion cap".into()));
    }
    let version = reader.manifest().get("qvf_version").cloned().unwrap_or(Value::Null);
    if version.is_boolean() || !matches!(version.as_f64(), Some(v) if v == 1.0 || v == 2.0) {
        return Err(Fail(format!("E_VERSION: unsupported qvf_version {}", py_repr(&version))));
    }
    Ok(reader)
}

/// The sections `info` lists: none when `sections` holds a non-object, as
/// the Python model then keeps the list untyped.
fn listed_sections(reader: &Reader) -> Vec<qvf::Section<'_>> {
    let all_objects = reader.manifest().get("sections").and_then(Value::as_array).is_some_and(|l| l.iter().all(Value::is_object));
    if all_objects {
        reader.sections()
    } else {
        Vec::new()
    }
}

fn info(args: &[String]) -> Result<u8, Fail> {
    let a = parse(args, &["--json"])?;
    let [archive] = a.positional.as_slice() else {
        return Err(usage("info requires exactly one ARCHIVE"));
    };
    let reader = open_checked(archive)?;
    let mut sections = Vec::new();
    for s in listed_sections(&reader) {
        let members = s
            .members()
            .iter()
            .map(|(role, spec)| member_info(&reader, role, spec.as_object().unwrap_or(&Map::new())))
            .collect::<Result<Vec<_>, _>>()?;
        sections.push(json!({"id": s.id, "kind": s.kind, "members": members}));
    }
    let m = reader.manifest();
    let value = json!({
        "qvf_version": m.get("qvf_version").cloned().unwrap_or(Value::Null),
        "source": m.get("source").cloned().unwrap_or(Value::Null),
        "sections": sections,
    });
    if a.json {
        print_json(&value);
        return Ok(0);
    }
    println!("QVF {} | {} sections | {}", py_repr(&value["qvf_version"]), sections.len(), py_repr(&value["source"]));
    for s in &sections {
        println!("  {} [{}]", s["id"].as_str().unwrap_or_default(), s["kind"].as_str().unwrap_or_default());
        for m in s["members"].as_array().into_iter().flatten() {
            let tensor = match m.get("dtype") {
                Some(d) => format!(
                    " {}{}",
                    d.as_str().map(str::to_owned).unwrap_or_else(|| py_repr(d)),
                    py_repr(m.get("shape").unwrap_or(&json!([])))
                ),
                None => String::new(),
            };
            let text = |k: &str| m[k].as_str().map(str::to_owned).unwrap_or_else(|| py_repr(&m[k]));
            println!("    {}: {} ({}{tensor}, {} bytes)", text("role"), text("path"), text("format"), m["size"]);
        }
    }
    Ok(0)
}

fn cat(args: &[String]) -> Result<u8, Fail> {
    let a = parse(args, &[])?;
    let reader;
    match a.positional.as_slice() {
        [archive] => {
            reader = open_checked(archive)?;
            print_json(reader.manifest());
        }
        [archive, section, member] => {
            reader = open_checked(archive)?;
            print_json(&reader.read_json(section, member)?);
        }
        [_, _] => return Err(Fail("cat requires both SECTION and MEMBER, or neither for the manifest".into())),
        _ => return Err(usage("cat requires ARCHIVE [SECTION MEMBER]")),
    }
    Ok(0)
}

fn extract(args: &[String]) -> Result<u8, Fail> {
    let a = parse(args, &["-o"])?;
    let [archive, section, member] = a.positional.as_slice() else {
        return Err(usage("extract requires ARCHIVE SECTION MEMBER"));
    };
    let output = a.output.ok_or_else(|| usage("the following arguments are required: -o/--output"))?;
    let reader = open_checked(archive)?;
    // Verified (declared size and sha256) before anything is written.
    let bytes = reader.read_bytes(section, member)?;
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)
        .map_err(|e| Fail(format!("refusing to overwrite or create {output}: {e}")))?;
    if let Err(e) = file.write_all(&bytes).and_then(|_| file.sync_all()) {
        drop(file);
        let _ = std::fs::remove_file(&output);
        return Err(e.into());
    }
    println!("Extracted {section}/{member} to {output}");
    Ok(0)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("validate") => validate(&args[1..]),
        Some("info") => info(&args[1..]),
        Some("cat") => cat(&args[1..]),
        Some("extract") => extract(&args[1..]),
        Some("--version") => {
            println!("{}", env!("CARGO_PKG_VERSION"));
            Ok(0)
        }
        Some("-h" | "--help") => {
            println!("{USAGE}");
            Ok(0)
        }
        Some(other) => Err(usage(&format!("invalid choice: {other:?}"))),
        None => Err(usage("the following arguments are required: command")),
    };
    match result {
        Ok(code) => ExitCode::from(code),
        Err(Fail(message)) => {
            eprintln!("qvf: {message}");
            ExitCode::from(2)
        }
    }
}
