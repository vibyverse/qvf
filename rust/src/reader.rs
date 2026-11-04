//! The reader: `manifest.json`, member access and little-endian tensors
//! (spec § 2, § 4.2). A port of `js/src/reader.ts`.
//!
//! The reader applies the resource caps of spec § 6.1 and checks what it
//! reads (sha256, dtype and byte length), but it is not the validator: run
//! [`crate::validate`] for the findings of spec § 6.

use std::path::Path;

use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::error::{Code, QvfError};
use crate::json::{json_depth, parse_json_bytes};
use crate::limits::Limits;
use crate::zip::{path_problem, ratio_exceeded, ZipArchive, ZipEntry};

/// Lowercase hexadecimal SHA-256 of `data`.
pub fn sha256_hex(data: &[u8]) -> String {
    Sha256::digest(data).iter().map(|b| format!("{b:02x}")).collect()
}

/// A dtype of spec § 4.2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DType {
    Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    Float32,
    Float64,
}

impl DType {
    /// The dtype for its spec § 4.2 name.
    pub fn parse(name: &str) -> Option<DType> {
        Some(match name {
            "int8" => DType::Int8,
            "int16" => DType::Int16,
            "int32" => DType::Int32,
            "int64" => DType::Int64,
            "uint8" => DType::UInt8,
            "uint16" => DType::UInt16,
            "uint32" => DType::UInt32,
            "uint64" => DType::UInt64,
            "float32" => DType::Float32,
            "float64" => DType::Float64,
            _ => return None,
        })
    }

    /// Bytes per element.
    pub fn itemsize(self) -> usize {
        match self {
            DType::Int8 | DType::UInt8 => 1,
            DType::Int16 | DType::UInt16 => 2,
            DType::Int32 | DType::UInt32 | DType::Float32 => 4,
            DType::Int64 | DType::UInt64 | DType::Float64 => 8,
        }
    }
}

/// The elements of a binary member, decoded from little-endian bytes.
#[derive(Debug, Clone, PartialEq)]
pub enum TensorData {
    Int8(Vec<i8>),
    Int16(Vec<i16>),
    Int32(Vec<i32>),
    Int64(Vec<i64>),
    UInt8(Vec<u8>),
    UInt16(Vec<u16>),
    UInt32(Vec<u32>),
    UInt64(Vec<u64>),
    Float32(Vec<f32>),
    Float64(Vec<f64>),
}

/// A binary member with its declared dtype and shape (row-major, C order).
#[derive(Debug, Clone, PartialEq)]
pub struct Tensor {
    pub dtype: DType,
    pub shape: Vec<usize>,
    pub data: TensorData,
}

macro_rules! decode_le {
    ($bytes:expr, $t:ty, $n:expr) => {
        $bytes
            .chunks_exact($n)
            .map(|c| {
                let mut b = [0u8; $n];
                b.copy_from_slice(c);
                <$t>::from_le_bytes(b)
            })
            .collect()
    };
}

impl Tensor {
    /// Decode little-endian `bytes` of `dtype` and `shape`. The byte length
    /// must be `itemsize * prod(shape)` (an empty shape is a scalar).
    pub fn decode(dtype: DType, shape: Vec<usize>, bytes: &[u8]) -> Result<Tensor, QvfError> {
        let want = shape
            .iter()
            .try_fold(dtype.itemsize(), |n, d| n.checked_mul(*d))
            .ok_or_else(|| QvfError::new(Code::E_BINARY_SIZE, "shape overflows the address space"))?;
        if bytes.len() != want {
            return Err(QvfError::new(
                Code::E_BINARY_SIZE,
                format!("byte length {} != itemsize*prod(shape) {want}", bytes.len()),
            ));
        }
        let data = match dtype {
            DType::Int8 => TensorData::Int8(bytes.iter().map(|&b| b as i8).collect()),
            DType::UInt8 => TensorData::UInt8(bytes.to_vec()),
            DType::Int16 => TensorData::Int16(decode_le!(bytes, i16, 2)),
            DType::UInt16 => TensorData::UInt16(decode_le!(bytes, u16, 2)),
            DType::Int32 => TensorData::Int32(decode_le!(bytes, i32, 4)),
            DType::UInt32 => TensorData::UInt32(decode_le!(bytes, u32, 4)),
            DType::Float32 => TensorData::Float32(decode_le!(bytes, f32, 4)),
            DType::Int64 => TensorData::Int64(decode_le!(bytes, i64, 8)),
            DType::UInt64 => TensorData::UInt64(decode_le!(bytes, u64, 8)),
            DType::Float64 => TensorData::Float64(decode_le!(bytes, f64, 8)),
        };
        Ok(Tensor { dtype, shape, data })
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.shape.iter().product()
    }

    /// Whether the tensor has no elements.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Element `flat` (row-major) as `f64`; 64-bit integers above 2^53 round.
    pub fn get_f64(&self, flat: usize) -> Option<f64> {
        Some(match &self.data {
            TensorData::Int8(v) => f64::from(*v.get(flat)?),
            TensorData::Int16(v) => f64::from(*v.get(flat)?),
            TensorData::Int32(v) => f64::from(*v.get(flat)?),
            TensorData::Int64(v) => *v.get(flat)? as f64,
            TensorData::UInt8(v) => f64::from(*v.get(flat)?),
            TensorData::UInt16(v) => f64::from(*v.get(flat)?),
            TensorData::UInt32(v) => f64::from(*v.get(flat)?),
            TensorData::UInt64(v) => *v.get(flat)? as f64,
            TensorData::Float32(v) => f64::from(*v.get(flat)?),
            TensorData::Float64(v) => *v.get(flat)?,
        })
    }

    /// The element at a full multi-index, or `None` when it is out of range.
    pub fn at(&self, index: &[usize]) -> Option<f64> {
        if index.len() != self.shape.len() {
            return None;
        }
        let mut flat = 0usize;
        for (i, d) in index.iter().zip(&self.shape) {
            if i >= d {
                return None;
            }
            flat = flat * d + i;
        }
        self.get_f64(flat)
    }

    /// Every element as `f64`, row-major.
    pub fn to_f64_vec(&self) -> Vec<f64> {
        (0..self.len()).filter_map(|k| self.get_f64(k)).collect()
    }
}

/// A section of the manifest with a string `id`, a string `kind` and an
/// object `members`.
#[derive(Debug, Clone, Copy)]
pub struct Section<'a> {
    pub id: &'a str,
    pub kind: &'a str,
    pub raw: &'a serde_json::Map<String, Value>,
}

impl<'a> Section<'a> {
    /// The member specs of the section.
    pub fn members(&self) -> &'a serde_json::Map<String, Value> {
        self.raw.get("members").and_then(Value::as_object).expect("checked at construction")
    }
}

/// Spec § 4.2: a declared `size` must equal the uncompressed size the ZIP
/// entry records.
pub(crate) fn declared_size_problem(
    spec: &serde_json::Map<String, Value>,
    entry_size: u64,
    name: &str,
) -> Option<(Code, String)> {
    let size = spec.get("size")?;
    match size.as_f64() {
        Some(x) if x.fract() == 0.0 && x >= 0.0 => {
            if x != entry_size as f64 {
                Some((
                    Code::E_BINARY_SIZE,
                    format!(
                        "member '{name}': declared size {} != {entry_size} bytes recorded in the ZIP central directory",
                        crate::json::num(x)
                    ),
                ))
            } else {
                None
            }
        }
        _ => Some((Code::E_MEMBER_MALFORMED, format!("member '{name}': size must be a non-negative integer, got {size}"))),
    }
}

fn check_entry_caps(entry: &ZipEntry, limits: &Limits) -> Result<(), QvfError> {
    if entry.size as f64 > limits.max_member_bytes {
        return Err(QvfError::new(
            Code::E_MEMBER_TOO_LARGE,
            format!(
                "ZIP entry {:?} exceeds the per-member size cap ({} > {} bytes uncompressed)",
                entry.name, entry.size, limits.max_member_bytes
            ),
        ));
    }
    if ratio_exceeded(entry, limits) {
        return Err(QvfError::new(
            Code::E_COMPRESSION_RATIO,
            format!(
                "ZIP entry {:?} expands {} bytes from {}, over the ratio cap {}",
                entry.name, entry.size, entry.compressed_size, limits.max_ratio
            ),
        ));
    }
    Ok(())
}

/// An opened `.qvf` archive.
#[derive(Debug)]
pub struct Reader {
    zip: ZipArchive,
    limits: Limits,
    manifest: Value,
}

impl Reader {
    /// Open an archive file with the default limits.
    pub fn open(path: impl AsRef<Path>) -> Result<Reader, QvfError> {
        let bytes = std::fs::read(path.as_ref())
            .map_err(|e| QvfError::new(Code::E_ZIP_UNREADABLE, format!("cannot read {}: {e}", path.as_ref().display())))?;
        Reader::from_bytes(bytes, Limits::default())
    }

    /// Open an archive held in memory.
    ///
    /// Refuses an archive with a duplicated entry name other than
    /// `manifest.json` (several `manifest.json` entries are an archive updated
    /// in place, and the last one is the manifest, spec § 2.3), a missing,
    /// oversized, too deep or invalid manifest, and a `byte_order` other than
    /// `"little"` (spec § 2.2).
    pub fn from_bytes(bytes: Vec<u8>, limits: Limits) -> Result<Reader, QvfError> {
        let zip = ZipArchive::open(bytes)?;
        let duplicates: Vec<&String> = zip.duplicates().iter().filter(|n| *n != "manifest.json").collect();
        if !duplicates.is_empty() {
            return Err(QvfError::new(
                Code::E_ZIP_DUPLICATE_ENTRY,
                format!("duplicate ZIP entry name(s) {duplicates:?}; refusing an ambiguous archive"),
            ));
        }
        let entry = zip.find("manifest.json").ok_or_else(|| QvfError::new(Code::E_MANIFEST_MISSING, "missing manifest.json"))?;
        if entry.size as f64 > limits.max_json_bytes {
            return Err(QvfError::new(
                Code::E_JSON_TOO_LARGE,
                format!("manifest.json is {} bytes, over the {}-byte JSON cap", entry.size, limits.max_json_bytes),
            ));
        }
        check_entry_caps(entry, &limits)?;
        let raw = zip.read(entry, Some(limits.max_member_bytes))?;
        if json_depth(&raw, limits.max_json_depth) > limits.max_json_depth {
            return Err(QvfError::new(
                Code::E_JSON_TOO_DEEP,
                format!("manifest.json nests deeper than {} levels", limits.max_json_depth),
            ));
        }
        let manifest = parse_json_bytes(&raw)
            .map_err(|e| QvfError::new(Code::E_MANIFEST_JSON, format!("manifest.json is not valid JSON: {e}")))?;
        let Some(obj) = manifest.as_object() else {
            return Err(QvfError::new(Code::E_MANIFEST_JSON, "manifest.json must be a JSON object"));
        };
        if let Some(order) = obj.get("byte_order") {
            if order != "little" {
                return Err(QvfError::new(
                    Code::E_SCHEMA,
                    format!("byte_order {order} is not \"little\", the only byte order of QVF v1"),
                ));
            }
        }
        Ok(Reader { zip, limits, manifest })
    }

    /// The parsed `manifest.json`.
    pub fn manifest(&self) -> &Value {
        &self.manifest
    }

    /// The limits the reader applies.
    pub fn limits(&self) -> &Limits {
        &self.limits
    }

    /// The underlying ZIP archive.
    pub fn zip(&self) -> &ZipArchive {
        &self.zip
    }

    /// Whether `qvf_version` is one the registry says a consumer opens.
    pub fn version_readable(&self) -> bool {
        self.manifest
            .get("qvf_version")
            .and_then(Value::as_f64)
            .is_some_and(|v| crate::generated::READABLE_VERSIONS.iter().any(|r| *r as f64 == v))
    }

    /// The well-formed sections (string `id`, string `kind`, object
    /// `members`), in manifest order.
    pub fn sections(&self) -> Vec<Section<'_>> {
        let Some(list) = self.manifest.get("sections").and_then(Value::as_array) else { return Vec::new() };
        list.iter()
            .filter_map(|s| {
                let raw = s.as_object()?;
                let id = raw.get("id")?.as_str()?;
                let kind = raw.get("kind")?.as_str()?;
                raw.get("members")?.as_object()?;
                Some(Section { id, kind, raw })
            })
            .collect()
    }

    /// The first well-formed section with that id.
    pub fn section(&self, id: &str) -> Option<Section<'_>> {
        self.sections().into_iter().find(|s| s.id == id)
    }

    /// The well-formed sections of a kind.
    pub fn sections_of(&self, kind: &str) -> Vec<Section<'_>> {
        self.sections().into_iter().filter(|s| s.kind == kind).collect()
    }

    /// The member spec of `(section, role)`.
    pub fn member_spec(&self, section_id: &str, role: &str) -> Result<&serde_json::Map<String, Value>, QvfError> {
        let section = self
            .section(section_id)
            .ok_or_else(|| QvfError::new(Code::E_REF_UNRESOLVED, format!("section id not found: {section_id}")))?;
        section
            .members()
            .get(role)
            .ok_or_else(|| QvfError::new(Code::E_MEMBER_MISSING, format!("section {section_id:?} has no member {role:?}")))?
            .as_object()
            .ok_or_else(|| {
                QvfError::new(Code::E_MEMBER_MALFORMED, format!("section {section_id:?} member {role:?}: malformed spec"))
            })
    }

    fn load(&self, section_id: &str, role: &str, verify: bool, json: bool) -> Result<Vec<u8>, QvfError> {
        let spec = self.member_spec(section_id, role)?;
        let Some(path) = spec.get("path").and_then(Value::as_str) else {
            return Err(QvfError::new(
                Code::E_MEMBER_MALFORMED,
                format!("section {section_id:?} member {role:?}: malformed spec"),
            ));
        };
        if let Some(problem) = path_problem(path) {
            return Err(QvfError::new(Code::E_MEMBER_PATH, format!("member {path:?}: illegal path ({problem})")));
        }
        let entry = self
            .zip
            .find(path)
            .ok_or_else(|| QvfError::new(Code::E_MEMBER_MISSING, format!("member {path:?} declared but not in archive")))?;
        if self.zip.total_size() as f64 > self.limits.max_total_bytes {
            return Err(QvfError::new(
                Code::E_TOTAL_TOO_LARGE,
                format!(
                    "archive expands to {} bytes, over the {}-byte total cap",
                    self.zip.total_size(),
                    self.limits.max_total_bytes
                ),
            ));
        }
        if json && entry.size as f64 > self.limits.max_json_bytes {
            return Err(QvfError::new(
                Code::E_JSON_TOO_LARGE,
                format!("member {path:?} is {} bytes, over the {}-byte JSON cap", entry.size, self.limits.max_json_bytes),
            ));
        }
        check_entry_caps(entry, &self.limits)?;
        if let Some((code, msg)) = declared_size_problem(spec, entry.size, path) {
            return Err(QvfError::new(code, msg));
        }
        let bytes = self.zip.read(entry, Some(self.limits.max_member_bytes))?;
        if verify {
            if let Some(declared) = spec.get("sha256") {
                if declared.as_str() != Some(sha256_hex(&bytes).as_str()) {
                    return Err(QvfError::new(Code::E_SHA256_MISMATCH, format!("sha256 mismatch for {path}")));
                }
            }
        }
        Ok(bytes)
    }

    /// The raw bytes of a member, its sha256 verified when declared.
    pub fn read_bytes(&self, section_id: &str, role: &str) -> Result<Vec<u8>, QvfError> {
        self.load(section_id, role, true, false)
    }

    /// A JSON member, parsed.
    pub fn read_json(&self, section_id: &str, role: &str) -> Result<Value, QvfError> {
        let bytes = self.load(section_id, role, true, true)?;
        if json_depth(&bytes, self.limits.max_json_depth) > self.limits.max_json_depth {
            return Err(QvfError::new(
                Code::E_JSON_TOO_DEEP,
                format!("member {section_id}/{role} nests deeper than {} levels", self.limits.max_json_depth),
            ));
        }
        parse_json_bytes(&bytes)
            .map_err(|e| QvfError::new(Code::E_JSON_INVALID, format!("member {section_id}/{role} is not valid JSON ({e})")))
    }

    /// A binary member decoded with its manifest `dtype` and `shape` (an
    /// absent shape is a scalar).
    pub fn read_tensor(&self, section_id: &str, role: &str) -> Result<Tensor, QvfError> {
        let spec = self.member_spec(section_id, role)?;
        let dtype_value = spec.get("dtype");
        let dtype = dtype_value.and_then(Value::as_str).and_then(DType::parse).ok_or_else(|| {
            QvfError::new(Code::E_DTYPE, format!("member {section_id}/{role}: bad dtype {}", crate::json::show(dtype_value)))
        })?;
        let shape: Vec<usize> = match spec.get("shape") {
            None => Vec::new(),
            Some(Value::Array(dims)) => {
                dims.iter().map(|d| d.as_f64().filter(|x| x.fract() == 0.0 && *x >= 0.0 && *x <= usize::MAX as f64).map(|x| x as usize)).collect::<Option<Vec<_>>>().ok_or_else(
                    || {
                        QvfError::new(
                            Code::E_MEMBER_MALFORMED,
                            format!("member {section_id}/{role}: shape must be a list of non-negative integers"),
                        )
                    },
                )?
            }
            Some(_) => {
                return Err(QvfError::new(
                    Code::E_MEMBER_MALFORMED,
                    format!("member {section_id}/{role}: shape must be a list of non-negative integers"),
                ))
            }
        };
        let bytes = self.load(section_id, role, true, false)?;
        Tensor::decode(dtype, shape, &bytes)
            .map_err(|e| QvfError::new(e.code, format!("member {section_id}/{role}: {}", e.message)))
    }
}
