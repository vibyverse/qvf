//! The archive writer: members are streamed into the ZIP as they are added,
//! and `manifest.json` (STORED) is written last (spec § 2.1).
//!
//! [`Writer`] writes typed members into sections: JSON documents, tensors of
//! any manifest `dtype` written little-endian with their `shape` (spec § 2.2,
//! § 4.2), and opaque bytes. Every member records its `sha256` and `size`.
//! A tensor can be streamed in slabs through a [`TensorSink`], so no member
//! has to be held in memory. The section builders for `structure`,
//! `wavefunction.gto`, `volume.*` and `excited_states` live in
//! [`crate::sections`].
//!
//! ```no_run
//! use qvf::writer::{Source, Writer};
//! use serde_json::json;
//!
//! let mut w = Writer::create("out.qvf", Source::new("my-code", "1.0"))?;
//! w.add_section("note", "x_mycode.note", serde_json::Map::new())?;
//! w.add_json("note", "text", &json!({"hello": "world"}))?;
//! w.finish()?;
//! # Ok::<(), qvf::writer::WriteError>(())
//! ```

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use serde_json::{Map, Value};
use sha2::{Digest, Sha256};

use crate::generated;
use crate::reader::DType;
use crate::zip::path_problem;
use crate::zip_writer::{Method, ZipWriter};

/// Why a write failed.
#[derive(Debug)]
pub enum WriteError {
    /// The sink failed.
    Io(io::Error),
    /// The caller asked for something the format does not allow, such as a
    /// duplicate section id, an illegal path or a tensor of the wrong size.
    /// Nothing was written for that call.
    Invalid(String),
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WriteError::Io(e) => write!(f, "write failed: {e}"),
            WriteError::Invalid(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for WriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            WriteError::Io(e) => Some(e),
            WriteError::Invalid(_) => None,
        }
    }
}

impl From<io::Error> for WriteError {
    fn from(e: io::Error) -> Self {
        WriteError::Io(e)
    }
}

pub(crate) fn invalid<T>(message: impl Into<String>) -> Result<T, WriteError> {
    Err(WriteError::Invalid(message.into()))
}

/// The producing program, the manifest `source` (spec § 3.2).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Source {
    pub program: String,
    pub version: String,
    /// A free label of the calculation (required, may be empty).
    pub calculation: String,
}

impl Source {
    pub fn new(program: impl Into<String>, version: impl Into<String>) -> Self {
        Source { program: program.into(), version: version.into(), calculation: String::new() }
    }

    /// The same source with a calculation label.
    pub fn with_calculation(mut self, calculation: impl Into<String>) -> Self {
        self.calculation = calculation.into();
        self
    }
}

/// How members other than `manifest.json` are compressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Compression {
    /// DEFLATE, which spec § 2.1 recommends.
    #[default]
    Deflate,
    /// STORE.
    Stored,
}

/// Options of a [`Writer`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WriterOptions {
    /// Compression of the members. `manifest.json` is always STORED.
    pub compression: Compression,
    /// Write ZIP64 fields and end records even when nothing needs them.
    pub force_zip64: bool,
    /// Write `metadata.created_utc` (the current UTC time) unless the
    /// caller's metadata sets it.
    pub created_utc: bool,
}

impl Default for WriterOptions {
    fn default() -> Self {
        WriterOptions { compression: Compression::Deflate, force_zip64: false, created_utc: true }
    }
}

/// A UTC calendar time, for the `created_utc` stamp and the ZIP entry times.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UtcTime {
    pub year: i64,
    pub month: u16,
    pub day: u16,
    pub hour: u16,
    pub minute: u16,
    pub second: u16,
}

impl UtcTime {
    /// The current time (1970-01-01 if the clock is before the epoch).
    pub fn now() -> Self {
        let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs());
        Self::from_unix(secs)
    }

    /// The UTC calendar time of a Unix time stamp.
    pub fn from_unix(secs: u64) -> Self {
        let days = (secs / 86_400) as i64;
        let rem = secs % 86_400;
        // Civil date from days since the epoch (proleptic Gregorian calendar).
        let z = days + 719_468;
        let era = z.div_euclid(146_097);
        let doe = z.rem_euclid(146_097);
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let day = (doy - (153 * mp + 2) / 5 + 1) as u16;
        let month = if mp < 10 { mp + 3 } else { mp - 9 } as u16;
        let year = yoe + era * 400 + i64::from(month <= 2);
        UtcTime { year, month, day, hour: (rem / 3600) as u16, minute: (rem % 3600 / 60) as u16, second: (rem % 60) as u16 }
    }

    /// `YYYY-MM-DDTHH:MM:SSZ`.
    pub fn iso8601(&self) -> String {
        format!("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z", self.year, self.month, self.day, self.hour, self.minute, self.second)
    }
}

/// A Rust number type that is one of the manifest dtypes.
pub trait Element: Copy {
    /// The dtype of this type.
    const DTYPE: DType;
    /// Append the little-endian bytes of `self`.
    fn put_le(self, out: &mut Vec<u8>);
}

macro_rules! element {
    ($t:ty, $d:ident) => {
        impl Element for $t {
            const DTYPE: DType = DType::$d;
            fn put_le(self, out: &mut Vec<u8>) {
                out.extend_from_slice(&self.to_le_bytes());
            }
        }
    };
}

element!(i8, Int8);
element!(i16, Int16);
element!(i32, Int32);
element!(i64, Int64);
element!(u8, UInt8);
element!(u16, UInt16);
element!(u32, UInt32);
element!(u64, UInt64);
element!(f32, Float32);
element!(f64, Float64);

/// The manifest name of a dtype.
pub fn dtype_name(d: DType) -> &'static str {
    match d {
        DType::Int8 => "int8",
        DType::Int16 => "int16",
        DType::Int32 => "int32",
        DType::Int64 => "int64",
        DType::UInt8 => "uint8",
        DType::UInt16 => "uint16",
        DType::UInt32 => "uint32",
        DType::UInt64 => "uint64",
        DType::Float32 => "float32",
        DType::Float64 => "float64",
    }
}

/// The byte length of a tensor of `dtype` and `shape`, if it fits in `u64`.
pub fn tensor_bytes(dtype: DType, shape: &[usize]) -> Option<u64> {
    shape.iter().try_fold(dtype.itemsize() as u64, |acc, &n| acc.checked_mul(n as u64))
}

/// The JSON bytes of a member or the manifest: UTF-8, two-space indent, as
/// the Python writer formats them. Non-finite numbers cannot occur, since
/// `serde_json::Value` cannot hold them.
pub fn json_bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec_pretty(value).expect("a serde_json::Value always serializes")
}

/// What the writer knows about an added section, for the builders' checks.
#[derive(Debug, Clone, Default)]
pub(crate) struct SectionInfo {
    pub kind: String,
    /// Atom count of a `structure`.
    pub atoms: Option<usize>,
    /// Orbital row count per spin channel ("" for a single set) of a
    /// `wavefunction.gto`.
    pub mo_rows: HashMap<String, usize>,
}

/// Root keys the writer itself owns.
const OWNED_ROOT_KEYS: &[&str] = &["qvf_version", "schema_uri", "byte_order", "source", "sections"];
/// `viewer_defaults` keys that are not section ids (spec § 3.1, § 3.4).
const RESERVED_IDS: &[&str] = &["auto_open", "bookmarks", "camera", "crossfade"];

/// Whether `id` is a legal section id: `[A-Za-z0-9_.-]+`.
pub fn is_section_id(id: &str) -> bool {
    !id.is_empty() && id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'.' || b == b'-')
}

/// Whether `kind` is a canonical kind or a vendor kind (`x_<ns>[.<name>]`).
pub fn is_known_kind(kind: &str) -> bool {
    if generated::CANONICAL_KINDS.contains(&kind) {
        return true;
    }
    let Some(rest) = kind.strip_prefix("x_") else { return false };
    let (ns, name) = match rest.split_once('.') {
        Some((ns, name)) => (ns, Some(name)),
        None => (rest, None),
    };
    !ns.is_empty()
        && ns.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
        && name.map_or(true, |n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b)))
}

/// A streaming QVF writer.
///
/// Sections are added with [`Writer::add_section`] (or a builder of
/// [`crate::sections`]) and their members with [`Writer::add_json`],
/// [`Writer::add_tensor`], [`Writer::begin_tensor`] and
/// [`Writer::add_opaque`]. Each member goes to the sink at once. Call
/// [`Writer::finish`] to write the manifest. A writer dropped without
/// `finish` leaves an archive without a manifest, which is not a QVF file.
pub struct Writer<W: Write> {
    zip: ZipWriter<W>,
    options: WriterOptions,
    source: Source,
    root: Map<String, Value>,
    sections: Vec<Map<String, Value>>,
    pub(crate) info: HashMap<String, SectionInfo>,
    paths: HashSet<String>,
    failed: bool,
}

impl Writer<BufWriter<File>> {
    /// Create (or truncate) the file at `path` and write an archive into it.
    pub fn create(path: impl AsRef<Path>, source: Source) -> Result<Self, WriteError> {
        Self::create_with(path, source, WriterOptions::default())
    }

    /// [`Writer::create`] with options.
    pub fn create_with(path: impl AsRef<Path>, source: Source, options: WriterOptions) -> Result<Self, WriteError> {
        let file = File::create(path)?;
        Writer::with_options(BufWriter::new(file), source, options)
    }
}

impl<W: Write> Writer<W> {
    /// A writer on any sink, which need not be seekable.
    pub fn new(sink: W, source: Source) -> Result<Self, WriteError> {
        Self::with_options(sink, source, WriterOptions::default())
    }

    /// A writer with options.
    pub fn with_options(sink: W, source: Source, options: WriterOptions) -> Result<Self, WriteError> {
        if source.program.is_empty() || source.version.is_empty() {
            return invalid("source.program and source.version must not be empty (spec § 3.2)");
        }
        Ok(Writer {
            zip: ZipWriter::new(sink, options.force_zip64),
            options,
            source,
            root: Map::new(),
            sections: Vec::new(),
            info: HashMap::new(),
            paths: HashSet::new(),
            failed: false,
        })
    }

    fn ensure_ok(&self) -> Result<(), WriteError> {
        if self.failed {
            return invalid("the writer failed earlier and cannot write a valid archive");
        }
        Ok(())
    }

    /// Set a root block of the manifest, such as `metadata`, `provenance`,
    /// `energies`, `extensions` or `viewer_defaults` (spec § 3). The keys the
    /// writer owns (`qvf_version`, `schema_uri`, `byte_order`, `source`,
    /// `sections`) are refused.
    pub fn set_root(&mut self, key: &str, value: Value) -> Result<(), WriteError> {
        if OWNED_ROOT_KEYS.contains(&key) {
            return invalid(format!("the writer owns the manifest key {key:?}"));
        }
        self.root.insert(key.to_owned(), value);
        Ok(())
    }

    /// Whether a section of this id was added, and its kind.
    pub fn section_kind(&self, id: &str) -> Option<&str> {
        self.info.get(id).map(|i| i.kind.as_str())
    }

    pub(crate) fn require_kind(&self, what: &str, id: &str, kinds: &[&str]) -> Result<&SectionInfo, WriteError> {
        match self.info.get(id) {
            Some(i) if kinds.contains(&i.kind.as_str()) => Ok(i),
            Some(i) => invalid(format!("{what} {id:?} must name a {} section, not a {:?} section", kinds.join(" or "), i.kind)),
            None => invalid(format!("{what} {id:?} names no {} section added so far", kinds.join(" or "))),
        }
    }

    /// Add a section without members. `fields` are its other keys (`label`,
    /// `structure_ref`, kind-specific fields). Members are attached with the
    /// `add_*` member methods.
    pub fn add_section(&mut self, id: &str, kind: &str, fields: Map<String, Value>) -> Result<(), WriteError> {
        self.ensure_ok()?;
        if !is_section_id(id) {
            return invalid(format!("section id {id:?} must match [A-Za-z0-9_.-]+ (spec § 4.1)"));
        }
        if RESERVED_IDS.contains(&id) {
            return invalid(format!("section id {id:?} is reserved (spec § 3.1, § 3.4)"));
        }
        if self.info.contains_key(id) {
            return invalid(format!("duplicate section id {id:?}"));
        }
        if !is_known_kind(kind) {
            return invalid(format!("unknown section kind {kind:?}"));
        }
        for key in ["id", "kind", "members"] {
            if fields.contains_key(key) {
                return invalid(format!("section field {key:?} is set by the writer"));
            }
        }
        let mut section = Map::new();
        section.insert("id".into(), Value::from(id));
        section.insert("kind".into(), Value::from(kind));
        section.extend(fields);
        section.insert("members".into(), Value::Object(Map::new()));
        self.sections.push(section);
        self.info.insert(id.to_owned(), SectionInfo { kind: kind.to_owned(), ..SectionInfo::default() });
        Ok(())
    }

    fn members_of(&mut self, id: &str) -> Result<&mut Map<String, Value>, WriteError> {
        let section = self.sections.iter_mut().find(|s| s.get("id").and_then(Value::as_str) == Some(id));
        match section.and_then(|s| s.get_mut("members")).and_then(Value::as_object_mut) {
            Some(m) => Ok(m),
            None => invalid(format!("no section {id:?}")),
        }
    }

    fn check_target(&mut self, id: &str, role: &str) -> Result<(), WriteError> {
        self.ensure_ok()?;
        if role.is_empty() {
            return invalid("member role must not be empty");
        }
        if self.members_of(id)?.contains_key(role) {
            return invalid(format!("section {id:?} already has a member {role:?}"));
        }
        Ok(())
    }

    fn claim_path(&mut self, path: &str) -> Result<(), WriteError> {
        if path.is_empty() || path == "manifest.json" || path.ends_with('/') {
            return invalid(format!("illegal member path {path:?}"));
        }
        if let Some(why) = path_problem(path) {
            return invalid(format!("illegal member path {path:?}: {why} (spec § 2.1)"));
        }
        if !self.paths.insert(path.to_owned()) {
            return invalid(format!("duplicate member path {path:?}"));
        }
        Ok(())
    }

    fn method(&self) -> Method {
        match self.options.compression {
            Compression::Deflate => Method::Deflate,
            Compression::Stored => Method::Store,
        }
    }

    /// Write one whole entry; a failure poisons the writer.
    fn write_entry(&mut self, path: &str, method: Method, data: &[u8]) -> Result<(), WriteError> {
        let result = (|| {
            self.zip.start_entry(path, method, Some(data.len() as u64))?;
            self.zip.write_data(data)?;
            self.zip.finish_entry()
        })();
        if let Err(e) = result {
            self.failed = true;
            return Err(e.into());
        }
        Ok(())
    }

    fn attach(&mut self, id: &str, role: &str, member: Map<String, Value>) -> Result<Value, WriteError> {
        let value = Value::Object(member);
        self.members_of(id)?.insert(role.to_owned(), value.clone());
        Ok(value)
    }

    /// Add a JSON member (path `<id>/<role>.json`). Returns the member object
    /// written into the manifest.
    pub fn add_json(&mut self, id: &str, role: &str, payload: &Value) -> Result<Value, WriteError> {
        self.add_json_at(id, role, payload, &format!("{id}/{role}.json"))
    }

    /// [`Writer::add_json`] at an explicit path.
    pub fn add_json_at(&mut self, id: &str, role: &str, payload: &Value, path: &str) -> Result<Value, WriteError> {
        self.check_target(id, role)?;
        self.claim_path(path)?;
        let data = json_bytes(payload);
        let method = self.method();
        self.write_entry(path, method, &data)?;
        let mut m = Map::new();
        m.insert("path".into(), Value::from(path));
        m.insert("format".into(), Value::from("json"));
        m.insert("sha256".into(), Value::from(hex(&Sha256::digest(&data))));
        m.insert("size".into(), Value::from(data.len() as u64));
        self.attach(id, role, m)
    }

    /// Add an opaque binary member (no dtype or shape, spec § 4.2) at
    /// `<id>/<role>.bin`, with an optional media type. Images (PNG, JPEG,
    /// SVG) are STORED.
    pub fn add_opaque(&mut self, id: &str, role: &str, data: &[u8], media_type: Option<&str>) -> Result<Value, WriteError> {
        self.check_target(id, role)?;
        if let Some(mt) = media_type {
            let ok = mt
                .split_once('/')
                .is_some_and(|(a, b)| !a.is_empty() && !b.is_empty() && !mt.contains(char::is_whitespace) && !b.contains('/'));
            if !ok {
                return invalid(format!("media type {mt:?} is not type/subtype (spec § 4.2)"));
            }
        }
        let path = format!("{id}/{role}.bin");
        self.claim_path(&path)?;
        let method = match media_type {
            Some("image/png" | "image/jpeg" | "image/svg+xml") => Method::Store,
            _ => self.method(),
        };
        self.write_entry(&path, method, data)?;
        let mut m = Map::new();
        m.insert("path".into(), Value::from(path));
        m.insert("format".into(), Value::from("binary"));
        m.insert("sha256".into(), Value::from(hex(&Sha256::digest(data))));
        m.insert("size".into(), Value::from(data.len() as u64));
        if let Some(mt) = media_type {
            m.insert("media_type".into(), Value::from(mt));
        }
        self.attach(id, role, m)
    }

    /// Add a tensor member (path `<id>/<role>.dat`): `data` in C order with
    /// `shape`, written little-endian whatever the host byte order.
    pub fn add_tensor<T: Element>(&mut self, id: &str, role: &str, data: &[T], shape: &[usize]) -> Result<Value, WriteError> {
        let expected = tensor_bytes(T::DTYPE, shape);
        if expected != Some((data.len() * T::DTYPE.itemsize()) as u64) {
            return invalid(format!("E_BINARY_SIZE: {} elements do not fill shape {shape:?} of member {id}/{role}", data.len()));
        }
        let mut sink = self.begin_tensor(id, role, T::DTYPE, shape)?;
        sink.write(data)?;
        sink.finish()
    }

    /// Begin a tensor member written in slabs. The sink must receive exactly
    /// the bytes of `shape` and be finished before anything else is written.
    pub fn begin_tensor(&mut self, id: &str, role: &str, dtype: DType, shape: &[usize]) -> Result<TensorSink<'_, W>, WriteError> {
        self.begin_tensor_at(id, role, dtype, shape, &format!("{id}/{role}.dat"))
    }

    /// [`Writer::begin_tensor`] at an explicit path.
    pub fn begin_tensor_at(
        &mut self,
        id: &str,
        role: &str,
        dtype: DType,
        shape: &[usize],
        path: &str,
    ) -> Result<TensorSink<'_, W>, WriteError> {
        self.check_target(id, role)?;
        let Some(expected) = tensor_bytes(dtype, shape) else {
            return invalid(format!("shape {shape:?} overflows"));
        };
        self.claim_path(path)?;
        let method = self.method();
        if let Err(e) = self.zip.start_entry(path, method, Some(expected)) {
            self.failed = true;
            return Err(e.into());
        }
        Ok(TensorSink {
            writer: self,
            id: id.to_owned(),
            role: role.to_owned(),
            path: path.to_owned(),
            dtype,
            shape: shape.to_vec(),
            expected,
            written: 0,
            hash: Sha256::new(),
            buf: Vec::new(),
            done: false,
        })
    }

    /// The manifest as it would be written now.
    pub fn manifest(&self) -> Value {
        let mut m = Map::new();
        m.insert("qvf_version".into(), Value::from(generated::QVF_VERSION));
        m.insert("schema_uri".into(), Value::from(generated::SCHEMA_ID));
        m.insert("byte_order".into(), Value::from("little"));
        let mut source = Map::new();
        source.insert("program".into(), Value::from(self.source.program.as_str()));
        source.insert("version".into(), Value::from(self.source.version.as_str()));
        source.insert("calculation".into(), Value::from(self.source.calculation.as_str()));
        m.insert("source".into(), Value::Object(source));
        for (k, v) in &self.root {
            m.insert(k.clone(), v.clone());
        }
        m.insert("sections".into(), Value::Array(self.sections.iter().cloned().map(Value::Object).collect()));
        Value::Object(m)
    }

    /// Write `manifest.json` (STORED, last) and the ZIP central directory,
    /// and return the sink.
    pub fn finish(mut self) -> Result<W, WriteError> {
        self.ensure_ok()?;
        if self.options.created_utc {
            let meta = self.root.entry("metadata").or_insert_with(|| Value::Object(Map::new()));
            if let Some(obj) = meta.as_object_mut() {
                obj.entry("created_utc").or_insert_with(|| Value::from(UtcTime::now().iso8601()));
            }
        }
        let raw = json_bytes(&self.manifest());
        self.write_entry("manifest.json", Method::Store, &raw)?;
        Ok(self.zip.finish()?)
    }
}

/// An open tensor member. Write slabs with [`TensorSink::write`] or
/// [`TensorSink::write_le_bytes`] and close it with [`TensorSink::finish`].
/// A sink dropped unfinished poisons its writer.
pub struct TensorSink<'a, W: Write> {
    writer: &'a mut Writer<W>,
    id: String,
    role: String,
    path: String,
    dtype: DType,
    shape: Vec<usize>,
    expected: u64,
    written: u64,
    hash: Sha256,
    buf: Vec<u8>,
    done: bool,
}

impl<W: Write> TensorSink<'_, W> {
    /// Write elements of the member's dtype, in C order.
    pub fn write<T: Element>(&mut self, data: &[T]) -> Result<(), WriteError> {
        if T::DTYPE != self.dtype {
            return invalid(format!("member {} has dtype {}, not {}", self.path, dtype_name(self.dtype), dtype_name(T::DTYPE)));
        }
        let mut buf = std::mem::take(&mut self.buf);
        for chunk in data.chunks(8192) {
            buf.clear();
            for &x in chunk {
                x.put_le(&mut buf);
            }
            let r = self.write_le_bytes(&buf);
            if r.is_err() {
                self.buf = buf;
                return r;
            }
        }
        self.buf = buf;
        Ok(())
    }

    /// Write raw bytes that are already little-endian. Elements may be split
    /// across calls.
    pub fn write_le_bytes(&mut self, bytes: &[u8]) -> Result<(), WriteError> {
        if self.written + bytes.len() as u64 > self.expected {
            self.writer.failed = true;
            return invalid(format!("E_BINARY_SIZE: member {} exceeds its declared shape {:?}", self.path, self.shape));
        }
        if let Err(e) = self.writer.zip.write_data(bytes) {
            self.writer.failed = true;
            return Err(e.into());
        }
        self.hash.update(bytes);
        self.written += bytes.len() as u64;
        Ok(())
    }

    /// Close the member and record it in its section. Fails, and poisons the
    /// writer, when fewer bytes than the shape needs were written.
    pub fn finish(mut self) -> Result<Value, WriteError> {
        self.done = true;
        if self.written != self.expected {
            self.writer.failed = true;
            return invalid(format!(
                "E_BINARY_SIZE: member {} holds {} bytes, its shape {:?} needs {}",
                self.path, self.written, self.shape, self.expected
            ));
        }
        if let Err(e) = self.writer.zip.finish_entry() {
            self.writer.failed = true;
            return Err(e.into());
        }
        let mut m = Map::new();
        m.insert("path".into(), Value::from(self.path.as_str()));
        m.insert("format".into(), Value::from("binary"));
        m.insert("dtype".into(), Value::from(dtype_name(self.dtype)));
        m.insert("shape".into(), Value::from(self.shape.clone()));
        m.insert("sha256".into(), Value::from(hex(&std::mem::take(&mut self.hash).finalize())));
        m.insert("size".into(), Value::from(self.written));
        let (id, role) = (std::mem::take(&mut self.id), std::mem::take(&mut self.role));
        self.writer.attach(&id, &role, m)
    }
}

impl<W: Write> Drop for TensorSink<'_, W> {
    fn drop(&mut self) {
        if !self.done {
            self.writer.zip.abandon_entry();
            self.writer.failed = true;
        }
    }
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn civil_dates() {
        assert_eq!(UtcTime::from_unix(0).iso8601(), "1970-01-01T00:00:00Z");
        assert_eq!(UtcTime::from_unix(951_782_400).iso8601(), "2000-02-29T00:00:00Z");
        assert_eq!(UtcTime::from_unix(1_798_761_599).iso8601(), "2026-12-31T23:59:59Z");
    }

    #[test]
    fn vendor_kinds() {
        assert!(is_known_kind("structure"));
        assert!(is_known_kind("x_orca.gbw"));
        assert!(is_known_kind("x_orca"));
        assert!(!is_known_kind("x_"));
        assert!(!is_known_kind("structures"));
    }
}
