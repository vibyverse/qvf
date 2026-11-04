//! Archive operations of the `qvf` command line, ported from the Python
//! `qvf.operations`: [`diff`] compares decoded content independent of the
//! ZIP layout, and [`scrub`] removes selected personal content into a new,
//! repacked archive. Repacking is [`crate::Editor::repack`].

use std::collections::HashMap;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::Path;

use serde_json::{json, Map, Value};

use crate::editor::{copy_member, referenced_members};
use crate::error::{Code, QvfError};
use crate::reader::Reader;
use crate::writer::{invalid, UtcTime, WriteError};
use crate::zip::METHOD_STORE;
use crate::zip_writer::{Method, ZipWriter};

/// The chunk size the Python CLI compares binary members in.
const CHUNK: usize = 1 << 20;

/// Options of [`diff`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiffOptions {
    /// Absolute tolerance of floating-point numbers (≥ 0).
    pub atol: f64,
    /// Relative tolerance of floating-point numbers (≥ 0).
    pub rtol: f64,
    /// Leave out the root `metadata` (UUID, time stamps, authors).
    pub ignore_metadata: bool,
    /// How many differences to list (≥ 1); all are counted.
    pub max_differences: usize,
}

impl Default for DiffOptions {
    fn default() -> Self {
        DiffOptions { atol: 0.0, rtol: 0.0, ignore_metadata: false, max_differences: 20 }
    }
}

/// One difference: a JSON-pointer-like location and what differs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Difference {
    pub location: String,
    pub message: String,
}

/// The result of [`diff`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiffReport {
    pub equal: bool,
    pub difference_count: usize,
    /// The first `max_differences` differences.
    pub differences: Vec<Difference>,
    pub truncated: bool,
}

impl DiffReport {
    /// The report as the Python CLI prints it with `--json`.
    pub fn to_json(&self) -> Value {
        let list: Vec<Value> = self.differences.iter().map(|d| json!({"where": d.location, "message": d.message})).collect();
        json!({"equal": self.equal, "difference_count": self.difference_count, "differences": list, "truncated": self.truncated})
    }
}

/// `math.isclose(a, b, rel_tol=rtol, abs_tol=atol)` of Python.
fn isclose(a: f64, b: f64, atol: f64, rtol: f64) -> bool {
    if a == b {
        return true;
    }
    if a.is_infinite() || b.is_infinite() {
        return false;
    }
    let d = (a - b).abs();
    d <= (rtol * b.abs()) || d <= (rtol * a.abs()) || d <= atol
}

struct Differ<'o> {
    options: &'o DiffOptions,
    count: usize,
    list: Vec<Difference>,
}

fn is_int(n: &serde_json::Number) -> bool {
    n.is_i64() || n.is_u64()
}

impl Differ<'_> {
    fn record(&mut self, location: &str, message: &str) {
        self.count += 1;
        if self.list.len() < self.options.max_differences {
            self.list.push(Difference { location: location.to_owned(), message: message.to_owned() });
        }
    }

    fn compare(&mut self, a: &Value, b: &Value, location: &str) {
        match (a, b) {
            (Value::Object(x), Value::Object(y)) => {
                let mut keys: Vec<&String> = x.keys().chain(y.keys().filter(|k| !x.contains_key(*k))).collect();
                keys.sort();
                for key in keys {
                    let at = format!("{location}/{}", key.replace('~', "~0").replace('/', "~1"));
                    match (x.get(key), y.get(key)) {
                        (Some(p), Some(q)) => self.compare(p, q, &at),
                        (None, _) => self.record(&at, "missing from left"),
                        (_, None) => self.record(&at, "missing from right"),
                    }
                }
            }
            (Value::Array(x), Value::Array(y)) => {
                if x.len() != y.len() {
                    self.record(location, "array lengths differ");
                }
                for (i, (p, q)) in x.iter().zip(y).enumerate() {
                    self.compare(p, q, &format!("{location}/{i}"));
                }
            }
            (Value::Number(x), Value::Number(y)) => {
                let equal = if is_int(x) && is_int(y) {
                    x.as_i128() == y.as_i128()
                } else {
                    let (p, q) = (x.as_f64().unwrap_or(f64::NAN), y.as_f64().unwrap_or(f64::NAN));
                    isclose(p, q, self.options.atol, self.options.rtol)
                };
                if !equal {
                    self.record(location, "numeric values differ");
                }
            }
            _ => {
                if a != b {
                    self.record(location, "values differ");
                }
            }
        }
    }
}

/// The manifest without the transport details of its members (`path`,
/// `size`, `sha256`), and without `metadata` when asked.
fn logical(manifest: &Value, ignore_metadata: bool) -> Value {
    let mut m = manifest.clone();
    if let Some(obj) = m.as_object_mut() {
        if ignore_metadata {
            obj.remove("metadata");
        }
        for s in obj.get_mut("sections").and_then(Value::as_array_mut).into_iter().flatten() {
            for member in s.get_mut("members").and_then(Value::as_object_mut).into_iter().flat_map(|m| m.values_mut()) {
                if let Some(spec) = member.as_object_mut() {
                    for key in ["path", "size", "sha256"] {
                        spec.remove(key);
                    }
                }
            }
        }
    }
    m
}

fn member_bytes(reader: &Reader, spec: &Map<String, Value>) -> Result<Vec<u8>, QvfError> {
    let path = spec.get("path").and_then(Value::as_str).unwrap_or_default();
    let mut out = Vec::new();
    copy_member(reader, path, spec, |chunk| {
        out.extend_from_slice(chunk);
        Ok(())
    })
    .map_err(|e| match e {
        WriteError::Read(e) => e,
        other => QvfError::new(Code::E_ZIP_UNREADABLE, other.to_string()),
    })?;
    Ok(out)
}

fn float_at(bytes: &[u8], dtype: &str, k: usize) -> f64 {
    match dtype {
        "float32" => f64::from(f32::from_le_bytes(bytes[4 * k..4 * k + 4].try_into().expect("4 bytes"))),
        _ => f64::from_le_bytes(bytes[8 * k..8 * k + 8].try_into().expect("8 bytes")),
    }
}

/// Compare two archives by their logical manifests and decoded payloads,
/// ignoring the ZIP layout. Member paths, sizes and hashes are transport
/// details: contents are compared, and checksums are always verified.
/// Tolerances apply to floating-point numbers (JSON numbers that are not
/// both integers, and float32/float64 tensors); integers and opaque bytes
/// compare exactly, and NaN differs from everything.
pub fn diff(left: &Reader, right: &Reader, options: &DiffOptions) -> Result<DiffReport, QvfError> {
    if !(options.atol.is_finite() && options.atol >= 0.0 && options.rtol.is_finite() && options.rtol >= 0.0) {
        return Err(QvfError::new(Code::E_SCHEMA, "tolerances must be finite and non-negative"));
    }
    if options.max_differences == 0 {
        return Err(QvfError::new(Code::E_SCHEMA, "max_differences must be positive"));
    }
    let mut d = Differ { options, count: 0, list: Vec::new() };
    d.compare(
        &logical(left.manifest(), options.ignore_metadata),
        &logical(right.manifest(), options.ignore_metadata),
        "manifest",
    );
    for (section, roles) in left.sections_with_roles() {
        let Some(other) = right.section(section.id) else { continue };
        for role in &roles {
            let member = &section.members()[role];
            let Some(counterpart) = other.members().get(role) else { continue };
            let location = format!("{}/{role}", section.id);
            let (ma, mb) =
                (member.as_object().cloned().unwrap_or_default(), counterpart.as_object().cloned().unwrap_or_default());
            let format = |m: &Map<String, Value>| m.get("format").and_then(Value::as_str).map(str::to_owned);
            if format(&ma).as_deref() == Some("json") && format(&mb).as_deref() == Some("json") {
                let (a, b) = (left.read_json(section.id, role)?, right.read_json(other.id, role)?);
                d.compare(&a, &b, &location);
                continue;
            }
            let dtype = ma.get("dtype").and_then(Value::as_str);
            let float = format(&ma).as_deref() == Some("binary")
                && format(&mb).as_deref() == Some("binary")
                && dtype.is_some()
                && dtype == mb.get("dtype").and_then(Value::as_str)
                && ma.get("shape") == mb.get("shape")
                && matches!(dtype, Some("float32" | "float64"));
            let (a, b) = (member_bytes(left, &ma)?, member_bytes(right, &mb)?);
            let (ca, cb): (Vec<&[u8]>, Vec<&[u8]>) = (a.chunks(CHUNK).collect(), b.chunks(CHUNK).collect());
            let mut index = 0usize;
            for k in 0..ca.len().max(cb.len()) {
                let (x, y) = (ca.get(k).copied().unwrap_or(&[]), cb.get(k).copied().unwrap_or(&[]));
                if float && x.len() == y.len() {
                    let dtype = dtype.unwrap_or("float64");
                    let size = if dtype == "float32" { 4 } else { 8 };
                    for e in 0..x.len() / size {
                        if !isclose(float_at(x, dtype, e), float_at(y, dtype, e), options.atol, options.rtol) {
                            d.record(&format!("{location}[{index}]"), "numeric values differ");
                        }
                        index += 1;
                    }
                } else if x != y {
                    d.record(&location, "raw payload chunks differ");
                }
            }
        }
    }
    let truncated = d.count > d.list.len();
    Ok(DiffReport { equal: d.count == 0, difference_count: d.count, differences: d.list, truncated })
}

/// A random version-4 UUID. The randomness comes from the keys the standard
/// library seeds `RandomState` with from the operating system; it is unique,
/// not cryptographically strong.
pub fn uuid4() -> String {
    use std::hash::{BuildHasher, Hasher};
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_nanos());
    let mut words = [0u64; 2];
    for (k, w) in words.iter_mut().enumerate() {
        let mut h = std::collections::hash_map::RandomState::new().build_hasher();
        h.write_u128(nanos);
        h.write_usize(k);
        h.write_u32(std::process::id());
        *w = h.finish();
    }
    let mut b = [0u8; 16];
    b[..8].copy_from_slice(&words[0].to_le_bytes());
    b[8..].copy_from_slice(&words[1].to_le_bytes());
    b[6] = (b[6] & 0x0f) | 0x40;
    b[8] = (b[8] & 0x3f) | 0x80;
    let h: String = b.iter().map(|x| format!("{x:02x}")).collect();
    format!("{}-{}-{}-{}-{}", &h[..8], &h[8..12], &h[12..16], &h[16..20], &h[20..])
}

/// Options of [`scrub`].
#[derive(Debug, Clone)]
pub struct ScrubOptions {
    /// Report what would change and write nothing.
    pub dry_run: bool,
    /// Patterns whose matches are redacted in the `run.record` command and
    /// in its UTF-8 `input`, `log`, `attachment.*` and `files` members.
    pub patterns: Vec<crate::regex::Regex>,
    /// The literal text every match is replaced by.
    pub replacement: String,
}

impl Default for ScrubOptions {
    fn default() -> Self {
        ScrubOptions { dry_run: false, patterns: Vec::new(), replacement: "[redacted]".into() }
    }
}

/// A JSON value that keeps the order of object keys, for the members scrub
/// rewrites as the Python `json` module writes them.
#[derive(Debug, Clone, PartialEq)]
enum OValue {
    Null,
    Bool(bool),
    Number(serde_json::Number),
    String(String),
    Array(Vec<OValue>),
    Object(Vec<(String, OValue)>),
}

impl<'de> serde::Deserialize<'de> for OValue {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> serde::de::Visitor<'de> for V {
            type Value = OValue;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("JSON")
            }
            fn visit_bool<E>(self, b: bool) -> Result<OValue, E> {
                Ok(OValue::Bool(b))
            }
            fn visit_i64<E>(self, n: i64) -> Result<OValue, E> {
                Ok(OValue::Number(n.into()))
            }
            fn visit_u64<E>(self, n: u64) -> Result<OValue, E> {
                Ok(OValue::Number(n.into()))
            }
            fn visit_f64<E: serde::de::Error>(self, x: f64) -> Result<OValue, E> {
                serde_json::Number::from_f64(x).map(OValue::Number).ok_or_else(|| E::custom("non-finite number"))
            }
            fn visit_str<E>(self, v: &str) -> Result<OValue, E> {
                Ok(OValue::String(v.to_owned()))
            }
            fn visit_string<E>(self, v: String) -> Result<OValue, E> {
                Ok(OValue::String(v))
            }
            fn visit_unit<E>(self) -> Result<OValue, E> {
                Ok(OValue::Null)
            }
            fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut seq: A) -> Result<OValue, A::Error> {
                let mut out = Vec::new();
                while let Some(v) = seq.next_element()? {
                    out.push(v);
                }
                Ok(OValue::Array(out))
            }
            fn visit_map<A: serde::de::MapAccess<'de>>(self, mut map: A) -> Result<OValue, A::Error> {
                let mut out: Vec<(String, OValue)> = Vec::new();
                while let Some((k, v)) = map.next_entry::<String, OValue>()? {
                    match out.iter_mut().find(|(key, _)| *key == k) {
                        Some(slot) => slot.1 = v,
                        None => out.push((k, v)),
                    }
                }
                Ok(OValue::Object(out))
            }
        }
        d.deserialize_any(V)
    }
}

/// `repr` of a Python float.
fn py_float(x: f64) -> String {
    let e = format!("{x:e}");
    let (mantissa, exp) = e.split_once('e').unwrap_or((&e, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let (sign, mantissa) = mantissa.strip_prefix('-').map_or(("", mantissa), |m| ("-", m));
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if (-4..16).contains(&exp) {
        let point = exp + 1;
        let body = if point <= 0 {
            format!("0.{}{digits}", "0".repeat((-point) as usize))
        } else if point as usize >= digits.len() {
            format!("{digits}{}.0", "0".repeat(point as usize - digits.len()))
        } else {
            format!("{}.{}", &digits[..point as usize], &digits[point as usize..])
        };
        format!("{sign}{body}")
    } else {
        let m = if digits.len() > 1 { format!("{}.{}", &digits[..1], &digits[1..]) } else { digits };
        format!("{sign}{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    }
}

fn py_string(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\x08' => out.push_str("\\b"),
            '\x0c' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

impl OValue {
    /// `json.dumps(value, ensure_ascii=False, separators=(",", ":"))`.
    fn dumps(&self, out: &mut String) {
        match self {
            OValue::Null => out.push_str("null"),
            OValue::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            OValue::Number(n) => match n.as_f64() {
                Some(x) if !(n.is_i64() || n.is_u64()) => out.push_str(&py_float(x)),
                _ => out.push_str(&n.to_string()),
            },
            OValue::String(s) => py_string(s, out),
            OValue::Array(a) => {
                out.push('[');
                for (i, v) in a.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    v.dumps(out);
                }
                out.push(']');
            }
            OValue::Object(o) => {
                out.push('{');
                for (i, (k, v)) in o.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    py_string(k, out);
                    out.push(':');
                    v.dumps(out);
                }
                out.push('}');
            }
        }
    }
}

struct Redactor<'o> {
    options: &'o ScrubOptions,
    changes: Vec<Value>,
}

impl Redactor<'_> {
    fn text(&mut self, s: &str, location: &str) -> Option<String> {
        let mut total = 0;
        let mut value = s.to_owned();
        for p in &self.options.patterns {
            let (next, n) = p.replace_all(&value, &self.options.replacement);
            value = next;
            total += n;
        }
        if total > 0 {
            self.changes.push(json!({"where": location, "action": "redact", "matches": total}));
            Some(value)
        } else {
            None
        }
    }

    fn value(&mut self, v: &mut Value, location: &str) {
        match v {
            Value::Object(o) => {
                for (k, item) in o.iter_mut() {
                    self.value(item, &format!("{location}/{k}"));
                }
            }
            Value::Array(a) => {
                for (i, item) in a.iter_mut().enumerate() {
                    self.value(item, &format!("{location}/{i}"));
                }
            }
            Value::String(s) => {
                if let Some(t) = self.text(s, location) {
                    *s = t;
                }
            }
            _ => {}
        }
    }

    fn ordered(&mut self, v: &mut OValue, location: &str) {
        match v {
            OValue::Object(o) => {
                for (k, item) in o.iter_mut() {
                    self.ordered(item, &format!("{location}/{k}"));
                }
            }
            OValue::Array(a) => {
                for (i, item) in a.iter_mut().enumerate() {
                    self.ordered(item, &format!("{location}/{i}"));
                }
            }
            OValue::String(s) => {
                if let Some(t) = self.text(s, location) {
                    *s = t;
                }
            }
            _ => {}
        }
    }
}

/// Remove `provenance.hostname` and the `email` and `affiliation` of every
/// `metadata.authors` entry, redact the matches of `options.patterns` in each
/// `run.record` command and in its UTF-8 `input`, `log`, `attachment.*` and
/// JSON `files` members (other members are reported as skipped, never
/// certified), and write the result to `destination` (which must not exist)
/// as a new archive: a fresh `metadata.uuid` (the old one is appended to
/// `metadata.derived_from`), a new `created_utc`, no `modified_utc`, and only
/// the members the manifest references. Returns the report the Python CLI
/// prints. Patterns are the Python `re` subset of [`crate::regex`].
pub fn scrub(reader: &Reader, destination: Option<&Path>, options: &ScrubOptions) -> Result<Value, WriteError> {
    if !options.dry_run && destination.is_none() {
        return invalid("scrub requires an explicit output or dry_run=True");
    }
    let mut manifest = reader.manifest().clone();
    let mut changes: Vec<Value> = Vec::new();
    let mut skipped: Vec<Value> = Vec::new();
    let Some(root) = manifest.as_object_mut() else { return invalid("manifest.json is not an object") };
    if let Some(p) = root.get_mut("provenance").and_then(Value::as_object_mut) {
        if p.remove("hostname").is_some() {
            changes.push(json!({"where": "/provenance/hostname", "action": "remove"}));
        }
    }
    let metadata = root.entry("metadata").or_insert_with(|| Value::Object(Map::new()));
    let Some(metadata) = metadata.as_object_mut() else { return invalid("metadata is not an object") };
    if let Some(authors) = metadata.get_mut("authors").and_then(Value::as_array_mut) {
        for (i, author) in authors.iter_mut().enumerate() {
            let Some(a) = author.as_object_mut() else { continue };
            for field in ["email", "affiliation"] {
                if a.remove(field).is_some() {
                    changes.push(json!({"where": format!("/metadata/authors/{i}/{field}"), "action": "remove"}));
                }
            }
        }
    }
    // As in Python, every member must be an object with a path.
    for s in manifest.get("sections").and_then(Value::as_array).into_iter().flatten() {
        for m in s.get("members").and_then(Value::as_object).into_iter().flat_map(|m| m.values()) {
            if !m.get("path").is_some_and(Value::is_string) {
                return invalid("scrub: a member spec has no path");
            }
        }
    }
    let members = referenced_members(&manifest)?;
    let mut uses: HashMap<String, usize> = HashMap::new();
    for s in manifest.get("sections").and_then(Value::as_array).into_iter().flatten() {
        for m in s.get("members").and_then(Value::as_object).into_iter().flat_map(|m| m.values()) {
            *uses.entry(m["path"].as_str().unwrap_or_default().to_owned()).or_default() += 1;
        }
    }
    let order = reader
        .zip()
        .find("manifest.json")
        .and_then(|e| reader.zip().read(e, None).ok())
        .and_then(|raw| crate::json::member_roles_in_file_order(&raw))
        .unwrap_or_default();
    let mut redactor = Redactor { options, changes };
    let mut replacements: HashMap<String, Vec<u8>> = HashMap::new();
    let limits = *reader.limits();
    let n_sections = manifest.get("sections").and_then(Value::as_array).map_or(0, Vec::len);
    for i in 0..n_sections {
        let s = &mut manifest["sections"][i];
        if s.get("kind").and_then(Value::as_str) != Some("run.record") {
            continue;
        }
        let sid = s.get("id").and_then(Value::as_str).unwrap_or_default().to_owned();
        if let Some(command) = s.get_mut("command") {
            redactor.value(command, &format!("{sid}/command"));
        }
        let mut roles: Vec<String> = order.get(i).cloned().unwrap_or_default();
        for r in s.get("members").and_then(Value::as_object).into_iter().flat_map(|m| m.keys()) {
            if !roles.contains(r) {
                roles.push(r.clone());
            }
        }
        for role in roles {
            let at = format!("{sid}/{role}");
            if options.patterns.is_empty() {
                skipped.push(json!({"where": at, "reason": "no patterns selected"}));
                continue;
            }
            let spec = s["members"][&role].as_object().cloned().unwrap_or_default();
            let path = spec.get("path").and_then(Value::as_str).unwrap_or_default().to_owned();
            let data: Vec<u8>;
            let changed: bool;
            if role == "files" && spec.get("format").and_then(Value::as_str) == Some("json") {
                let original = reader.read_json(&sid, &role)?;
                let raw = reader.read_bytes(&sid, &role)?;
                let mut doc: OValue = serde_json::from_slice(&raw).map_err(|e| WriteError::Invalid(e.to_string()))?;
                let before = doc.clone();
                redactor.ordered(&mut doc, &at);
                let mut text = String::new();
                doc.dumps(&mut text);
                text.push('\n');
                changed = doc != before;
                let _ = original;
                data = text.into_bytes();
            } else if role == "input" || role == "log" || role.starts_with("attachment.") {
                let size = reader.zip().find(&path).map_or(0, |e| e.size);
                if size as f64 > limits.max_json_bytes {
                    return invalid(format!("text member exceeds scrub cap: {path:?}"));
                }
                let raw = reader.read_bytes(&sid, &role)?;
                let Ok(text) = std::str::from_utf8(&raw) else {
                    skipped.push(json!({"where": at, "reason": "non-UTF-8 attachment"}));
                    if role == "input" || role == "log" {
                        return invalid(format!("E_UTF8: {path:?}"));
                    }
                    continue;
                };
                match redactor.text(text, &at) {
                    Some(t) => {
                        changed = t.as_bytes() != &raw[..];
                        data = t.into_bytes();
                    }
                    None => {
                        changed = false;
                        data = raw;
                    }
                }
            } else {
                skipped.push(json!({"where": at, "reason": "outside text scrub scope"}));
                continue;
            }
            if changed {
                if uses.get(&path).copied().unwrap_or(0) > 1 {
                    return invalid("cannot scrub a shared member; separate its declarations first");
                }
                let member = &mut s["members"][&role];
                member["size"] = Value::from(data.len() as u64);
                member["sha256"] = Value::from(crate::reader::sha256_hex(&data));
                replacements.insert(path, data);
            }
        }
    }
    let mut changes = redactor.changes;
    changes.push(json!({"where": "/metadata", "action": "new UUID and derivation on write"}));
    let report = json!({
        "dry_run": options.dry_run,
        "changes": changes,
        "skipped": skipped,
        "scope": "hostname, author contacts, selected run.record text; other content unchanged",
    });
    if options.dry_run {
        return Ok(report);
    }
    let destination = destination.expect("checked above");
    let metadata = manifest.get_mut("metadata").and_then(Value::as_object_mut).expect("inserted above");
    if let Some(previous) = metadata.get("uuid").filter(|u| u.as_str().is_some_and(|s| !s.is_empty())).cloned() {
        let derived = metadata.entry("derived_from").or_insert_with(|| Value::Array(Vec::new()));
        match derived.as_array_mut() {
            Some(list) => list.push(previous),
            None => return invalid("metadata.derived_from is not an array"),
        }
    }
    metadata.insert("uuid".into(), Value::from(uuid4()));
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap_or_default();
    let stamp = UtcTime::from_unix(now.as_secs()).iso8601();
    metadata.insert("created_utc".into(), Value::from(format!("{}.{:06}Z", &stamp[..19], now.subsec_micros())));
    metadata.remove("modified_utc");

    let file = OpenOptions::new().write(true).create_new(true).open(destination)?;
    let written = (|| -> Result<(), WriteError> {
        let mut out = ZipWriter::new(BufWriter::new(file), false);
        let zip = reader.zip();
        let mut by_path: HashMap<String, Map<String, Value>> = members.into_iter().collect();
        // In the order the manifest first references each member, as Python writes them.
        let order: Vec<String> = referenced_members(&manifest)?.into_iter().map(|(p, _)| p).collect();
        for path in order {
            let spec = by_path.remove(&path).unwrap_or_default();
            if let Some(data) = replacements.get(&path) {
                out.start_entry(&path, Method::Deflate, Some(data.len() as u64))?;
                out.write_data(data)?;
                out.finish_entry()?;
                continue;
            }
            let entry = zip.find(&path).ok_or_else(|| {
                WriteError::Read(QvfError::new(Code::E_MEMBER_MISSING, format!("member {path:?} declared but not in archive")))
            })?;
            let method = if entry.method == METHOD_STORE { Method::Store } else { Method::Deflate };
            out.start_entry(&path, method, Some(entry.size))?;
            copy_member(reader, &path, &spec, |chunk| out.write_data(chunk))?;
            out.finish_entry()?;
        }
        let raw = serde_json::to_vec(&manifest).expect("a Value serializes");
        out.start_entry("manifest.json", Method::Store, Some(raw.len() as u64))?;
        out.write_data(&raw)?;
        out.finish_entry()?;
        let mut sink = out.finish()?;
        sink.flush()?;
        sink.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        Ok(())
    })();
    if let Err(e) = written {
        let _ = std::fs::remove_file(destination);
        return Err(e);
    }
    Ok(report)
}
