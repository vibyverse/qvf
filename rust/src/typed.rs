//! Typed read access to sections, the protocol of the Python `qvf.reader`
//! (`Reader.get`) and the C++ `Reader::get<K>`.
//!
//! The models themselves are generated from the schema `$defs` into
//! [`crate::kinds`]. Each canonical kind is a struct that implements
//! [`TypedKind`]: its section fields (typed, open fields in `extra`), each
//! JSON member as its payload model, each tensor member as a [`Tensor`] and
//! each opaque member as [`OpaqueData`]. Decoding checks the kind and the
//! payload shapes the models carry. It is no substitute for the validator:
//! value ranges, patterns and cross-references stay with [`crate::validate`].
//!
//! ```no_run
//! use qvf::kinds::{Structure, WavefunctionGTO};
//! let reader = qvf::Reader::open("h2.qvf")?;
//! let mol: Structure = reader.get("mol")?;
//! println!("{} atoms", mol.structure.atoms.len());
//! let wf: WavefunctionGTO = reader.get("wf")?;
//! println!("{} shells", wf.basis.shells.len());
//! # Ok::<(), qvf::QvfError>(())
//! ```

use std::collections::BTreeMap;
use std::io::Write;

use serde::de::DeserializeOwned;
use serde::Serialize;
use serde_json::{Map, Value};

use crate::error::{Code, QvfError};
use crate::kinds::{AnyKind, VolumeGridPayload};
use crate::reader::{Reader, Section, Tensor, TensorData};
use crate::writer::{WriteError, Writer};

/// A canonical section kind with a typed model.
pub trait TypedKind: Sized {
    /// The kind string of the registry.
    const KIND: &'static str;

    /// Decode a section of this kind, reading and verifying its members.
    fn from_section(reader: &Reader, section: Section<'_>) -> Result<Self, QvfError>;

    /// Write the section and its members (at the writer's default paths);
    /// returns the section id.
    fn to_writer<W: Write>(&self, writer: &mut Writer<W>) -> Result<String, WriteError>;

    /// The section id.
    fn id(&self) -> &str;
}

/// An opaque member: its bytes and declared media type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OpaqueData {
    pub bytes: Vec<u8>,
    pub media_type: Option<String>,
}

impl Reader {
    /// The section `id` decoded as the kind `K`. Fails with
    /// `E_REF_UNRESOLVED` for an unknown id and `E_KIND_UNKNOWN` when the
    /// section is of another kind.
    pub fn get<K: TypedKind>(&self, id: &str) -> Result<K, QvfError> {
        K::from_section(self, self.section_or_err(id)?)
    }

    /// The section `id` decoded with the model of its kind, or `None` for a
    /// vendor kind.
    pub fn get_any(&self, id: &str) -> Result<Option<AnyKind>, QvfError> {
        AnyKind::from_section(self, self.section_or_err(id)?)
    }

    /// The section `id` as a [`Volume`], for any `volume.*` kind.
    pub fn get_volume(&self, id: &str) -> Result<Volume, QvfError> {
        Volume::from_section(self, self.section_or_err(id)?)
    }

    fn section_or_err(&self, id: &str) -> Result<Section<'_>, QvfError> {
        self.section(id).ok_or_else(|| QvfError::new(Code::E_REF_UNRESOLVED, format!("section not found: {id}")))
    }
}

/// Any `volume.*` section: the grid and the data, whatever the quantity
/// (the Python `qvf.kinds.Volume` view). The per-kind models in
/// [`crate::kinds`] (`VolumeDensity`, `VolumeOrbital`, …) also type the
/// kind-specific section fields.
#[derive(Debug, Clone, PartialEq)]
pub struct Volume {
    /// The `volume.*` kind.
    pub kind: String,
    /// The section id.
    pub id: String,
    /// The section fields other than `id`, `kind` and `members`.
    pub fields: Map<String, Value>,
    pub grid: VolumeGridPayload,
    pub data: Tensor,
}

impl Volume {
    /// Decode a section of any `volume.*` kind.
    pub fn from_section(reader: &Reader, section: Section<'_>) -> Result<Self, QvfError> {
        if !section.kind.starts_with("volume.") {
            return Err(QvfError::new(
                Code::E_KIND_UNKNOWN,
                format!("section {:?} is a {:?} section, not a volume", section.id, section.kind),
            ));
        }
        let mut fields = section.raw.clone();
        for key in ["id", "kind", "members"] {
            fields.remove(key);
        }
        Ok(Volume {
            kind: section.kind.to_owned(),
            id: section.id.to_owned(),
            fields,
            grid: json_member(reader, section, "grid")?,
            data: tensor_member(reader, section, "data")?,
        })
    }
}

// --- helpers of the generated models --------------------------------------------

#[doc(hidden)]
pub fn check_kind(section: Section<'_>, kind: &str) -> Result<(), QvfError> {
    if section.kind == kind {
        Ok(())
    } else {
        Err(QvfError::new(
            Code::E_KIND_UNKNOWN,
            format!("section {:?} is a {:?} section, not {kind:?}", section.id, section.kind),
        ))
    }
}

#[doc(hidden)]
pub fn section_fields<T: DeserializeOwned>(section: Section<'_>) -> Result<T, QvfError> {
    let mut fields = section.raw.clone();
    fields.remove("kind");
    fields.remove("members");
    serde_json::from_value(Value::Object(fields))
        .map_err(|e| QvfError::new(Code::E_SCHEMA, format!("section {:?}: {e}", section.id)))
}

#[doc(hidden)]
pub fn json_member<T: DeserializeOwned>(reader: &Reader, section: Section<'_>, role: &str) -> Result<T, QvfError> {
    let value = reader.read_json(section.id, role)?;
    serde_json::from_value(value).map_err(|e| QvfError::new(Code::E_PAYLOAD_SCHEMA, format!("member {}/{role}: {e}", section.id)))
}

#[doc(hidden)]
pub fn tensor_member(reader: &Reader, section: Section<'_>, role: &str) -> Result<Tensor, QvfError> {
    reader.read_tensor(section.id, role)
}

#[doc(hidden)]
pub fn opaque_member(reader: &Reader, section: Section<'_>, role: &str) -> Result<OpaqueData, QvfError> {
    let media_type = reader.member_spec(section.id, role)?.get("media_type").and_then(Value::as_str).map(str::to_owned);
    Ok(OpaqueData { bytes: reader.read_bytes(section.id, role)?, media_type })
}

#[doc(hidden)]
pub fn optional<T>(
    section: Section<'_>,
    role: &str,
    read: impl FnOnce(&str) -> Result<T, QvfError>,
) -> Result<Option<T>, QvfError> {
    if section.members().contains_key(role) {
        read(role).map(Some)
    } else {
        Ok(None)
    }
}

#[doc(hidden)]
pub fn pattern_members<T>(
    section: Section<'_>,
    known: &[&str],
    mut read: impl FnMut(&str) -> Result<T, QvfError>,
) -> Result<BTreeMap<String, T>, QvfError> {
    let mut out = BTreeMap::new();
    for role in section.members().keys().filter(|r| !known.contains(&r.as_str())) {
        out.insert(role.clone(), read(role)?);
    }
    Ok(out)
}

fn to_json<T: Serialize>(value: &T) -> Result<Value, WriteError> {
    serde_json::to_value(value).map_err(|e| WriteError::Invalid(format!("cannot encode a typed model: {e}")))
}

#[doc(hidden)]
pub fn add_section<W: Write, T: Serialize>(writer: &mut Writer<W>, kind: &str, fields: &T) -> Result<String, WriteError> {
    let Value::Object(mut map) = to_json(fields)? else {
        return Err(WriteError::Invalid("section fields must encode to a JSON object".into()));
    };
    let Some(Value::String(id)) = map.remove("id") else {
        return Err(WriteError::Invalid("section fields need a string id".into()));
    };
    writer.add_section(&id, kind, map)?;
    Ok(id)
}

#[doc(hidden)]
pub fn put_json<W: Write, T: Serialize>(writer: &mut Writer<W>, id: &str, role: &str, value: &T) -> Result<(), WriteError> {
    writer.add_json(id, role, &to_json(value)?).map(drop)
}

#[doc(hidden)]
pub fn put_tensor<W: Write>(writer: &mut Writer<W>, id: &str, role: &str, t: &Tensor) -> Result<(), WriteError> {
    let shape = &t.shape;
    match &t.data {
        TensorData::Int8(d) => writer.add_tensor(id, role, d, shape),
        TensorData::Int16(d) => writer.add_tensor(id, role, d, shape),
        TensorData::Int32(d) => writer.add_tensor(id, role, d, shape),
        TensorData::Int64(d) => writer.add_tensor(id, role, d, shape),
        TensorData::UInt8(d) => writer.add_tensor(id, role, d, shape),
        TensorData::UInt16(d) => writer.add_tensor(id, role, d, shape),
        TensorData::UInt32(d) => writer.add_tensor(id, role, d, shape),
        TensorData::UInt64(d) => writer.add_tensor(id, role, d, shape),
        TensorData::Float32(d) => writer.add_tensor(id, role, d, shape),
        TensorData::Float64(d) => writer.add_tensor(id, role, d, shape),
    }
    .map(drop)
}

#[doc(hidden)]
pub fn put_opaque<W: Write>(writer: &mut Writer<W>, id: &str, role: &str, data: &OpaqueData) -> Result<(), WriteError> {
    writer.add_opaque(id, role, &data.bytes, data.media_type.as_deref()).map(drop)
}

/// serde helper for an optional nullable field: absent is `None`, `null` is
/// `Some(None)`.
#[doc(hidden)]
pub mod double_option {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<T: Serialize, S: Serializer>(value: &Option<Option<T>>, s: S) -> Result<S::Ok, S::Error> {
        match value {
            Some(inner) => inner.serialize(s),
            None => s.serialize_none(),
        }
    }

    pub fn deserialize<'de, T: Deserialize<'de>, D: Deserializer<'de>>(d: D) -> Result<Option<Option<T>>, D::Error> {
        Option::<T>::deserialize(d).map(Some)
    }
}
