//! In-place updates of an existing archive (spec § 2.3) and repacking.
//!
//! An [`Editor`] opens an archive, keeps the bytes of every entry where they
//! are, appends the new members after the last of them and writes a new
//! `manifest.json` (STORED, last) and a new central directory that lists the
//! old entries unchanged and then the new ones. The new manifest is complete:
//! it keeps `source`, the root blocks and the sections of the archive, and
//! records what the update added, replaced and removed. The previous manifest
//! and every entry it no longer references stay behind as orphans
//! (`W_ORPHAN_MEMBER`) until [`Editor::repack`] writes the archive anew from
//! its current manifest. This follows the Python `qvf.Editor` and the C++
//! `qvf::Editor`:
//!
//! * a new member never reuses an entry name: one whose path the archive
//!   already holds is written as `stem.N.suffix`, the first free N from 2;
//! * `metadata.modified_utc` is set to the time of the update unless the
//!   caller set one, and `created_utc` and `uuid` are kept;
//! * `job.spec` sections are kept verbatim (spec § 5.9);
//! * an archive of the withdrawn `qvf_version` 2 is written back as 1 with
//!   the v1 `schema_uri`.
//!
//! In place (the default), nothing of the archive changes before the first
//! new entry is written, and an editor dropped or aborted before
//! [`Editor::commit`] restores the old central directory, so the archive is
//! as it was. A process killed in between leaves an archive that cannot be
//! opened (spec § 2.3). `atomic` edits a copy and renames it over the
//! original on commit.
//!
//! ```no_run
//! use qvf::{Editor, Spectrum};
//! use serde_json::json;
//!
//! let mut ed = Editor::open("water.qvf")?;
//! ed.writer().add_spectrum(&Spectrum::new("spectra.ir", vec![1595.0, 3657.0], vec![70.0, 4.0]))?;
//! ed.replace_json("mol", "structure", &json!({"atoms": [], "pbc": [false, false, false]}))?;
//! ed.remove_section("old_volume")?;
//! ed.commit()?;
//! qvf::Editor::repack("water.qvf", "water-clean.qvf", None)?;
//! # Ok::<(), qvf::WriteError>(())
//! ```

use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use serde_json::{Map, Value};
use sha2::{Digest, Sha256};

use crate::error::{Code, QvfError};
use crate::limits::Limits;
use crate::reader::{check_entry_caps, declared_size_problem, Reader};
use crate::writer::{invalid, Compression, Element, SectionInfo, UpdateBase, WriteError, Writer, WriterOptions};
use crate::zip::{METHOD_DEFLATE, METHOD_STORE};
use crate::zip_writer::{Method, ZipWriter};

/// Options of an [`Editor`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EditorOptions {
    /// Edit a sibling copy and rename it over the archive on commit.
    pub atomic: bool,
    /// Compression of the new members.
    pub compression: Compression,
    /// ZIP64 fields for the new entries and end records.
    pub force_zip64: bool,
    /// The caps the archive is opened under (spec § 6.1).
    pub limits: Limits,
}

impl Default for EditorOptions {
    fn default() -> Self {
        EditorOptions { atomic: false, compression: Compression::Deflate, force_zip64: false, limits: Limits::default() }
    }
}

/// An update of an archive, see the module documentation.
pub struct Editor {
    writer: Option<Writer<BufWriter<File>>>,
    target: PathBuf,
    /// The file being written: `target`, or the sibling copy when atomic.
    working: PathBuf,
    atomic: bool,
    /// Where the old central directory began and the bytes from there to
    /// the end, to restore an aborted in-place update.
    tail: (u64, Vec<u8>),
}

fn read_error(e: QvfError) -> WriteError {
    WriteError::Read(e)
}

/// A new, exclusively created sibling of `target` for staged output.
fn sibling(target: &Path, prefix: &str) -> Result<(PathBuf, File), WriteError> {
    let dir = match target.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.subsec_nanos());
    for attempt in 0..64u32 {
        let path = dir.join(format!("{prefix}{}-{nanos:08x}-{attempt}.tmp", std::process::id()));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(f) => return Ok((path, f)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e.into()),
        }
    }
    invalid(format!("cannot reserve a temporary file beside {}", target.display()))
}

/// The path a replacement derives its fresh name from: the old path without
/// the `.N` a previous update added (`a/b.2.json` gives `a/b.json`, which
/// the writer then turns into the first free `a/b.N.json`).
fn base_path(path: &str) -> String {
    let (head, tail) = match path.rfind('/') {
        Some(i) => path.split_at(i + 1),
        None => ("", path),
    };
    let parts: Vec<&str> = tail.splitn(3, '.').collect();
    match parts.as_slice() {
        [stem, n, rest] if !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) => format!("{head}{stem}.{rest}"),
        [stem, n] if !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) => format!("{head}{stem}"),
        _ => path.to_owned(),
    }
}

/// What the builders need to know of the sections already in the archive.
fn section_info(reader: &Reader, section: &Map<String, Value>) -> SectionInfo {
    let kind = section.get("kind").and_then(Value::as_str).unwrap_or_default().to_owned();
    let id = section.get("id").and_then(Value::as_str).unwrap_or_default();
    let mut info = SectionInfo { kind, ..SectionInfo::default() };
    let members = section.get("members").and_then(Value::as_object);
    let rows = |role: &str| members?.get(role)?.get("shape")?.as_array()?.first()?.as_u64().map(|n| n as usize);
    match info.kind.as_str() {
        "structure" => {
            info.atoms = reader.read_json(id, "structure").ok().and_then(|p| p.get("atoms")?.as_array().map(Vec::len));
        }
        "trajectory" => {
            info.atoms = members.and_then(|m| m.get("coords")?.get("shape")?.as_array()?.get(1)?.as_u64()).map(|n| n as usize);
        }
        "wavefunction.gto" => {
            for (role, channel) in [("mo_coefficients", ""), ("mo_coefficients_alpha", "alpha"), ("mo_coefficients_beta", "beta")]
            {
                if let Some(n) = rows(role) {
                    info.mo_rows.insert(channel.to_owned(), n);
                }
            }
        }
        "excited_states" => {
            info.states = reader
                .read_json(id, "states")
                .ok()
                .and_then(|p| p.get("states")?.as_array()?.iter().map(|s| s.get("index")?.as_u64().map(|i| i as u32)).collect());
        }
        _ => {}
    }
    info
}

impl Editor {
    /// Open an archive for an in-place update with the default options.
    pub fn open(path: impl AsRef<Path>) -> Result<Editor, WriteError> {
        Self::open_with(path, EditorOptions::default())
    }

    /// Open an archive for an update.
    ///
    /// Refuses an archive the [`Reader`] refuses (a duplicated entry name
    /// other than `manifest.json`, no or an invalid manifest, the caps of
    /// spec § 6.1, a `byte_order` other than `"little"`), a `qvf_version`
    /// other than 1 or the withdrawn 2, a `source` without `program` and
    /// `version` strings, a `sections` array of anything but objects with a
    /// string `id` and object `members`, and an archive whose entries do not
    /// all lie before its central directory.
    pub fn open_with(path: impl AsRef<Path>, options: EditorOptions) -> Result<Editor, WriteError> {
        let target = path.as_ref().to_path_buf();
        let bytes = std::fs::read(&target)
            .map_err(|e| read_error(QvfError::new(Code::E_ZIP_UNREADABLE, format!("cannot read {}: {e}", target.display()))))?;
        let reader = Reader::from_bytes(bytes, options.limits).map_err(read_error)?;
        let manifest = reader.manifest().as_object().cloned().unwrap_or_default();
        let version = manifest.get("qvf_version");
        let v2 = match version {
            Some(Value::Number(n)) if n.as_f64() == Some(1.0) => false,
            Some(Value::Number(n)) if n.as_f64() == Some(2.0) => true,
            _ => return invalid("only qvf_version 1 archives (or the withdrawn 2, read as 1) can be updated"),
        };
        let source = match manifest.get("source") {
            Some(Value::Object(s))
                if s.get("program").is_some_and(Value::is_string) && s.get("version").is_some_and(Value::is_string) =>
            {
                s.clone()
            }
            _ => return invalid("manifest.json has no usable source block"),
        };
        let sections: Vec<Map<String, Value>> = match manifest.get("sections").and_then(Value::as_array) {
            Some(list)
                if list
                    .iter()
                    .all(|s| s.get("id").is_some_and(Value::is_string) && s.get("members").is_some_and(Value::is_object)) =>
            {
                list.iter().filter_map(|s| s.as_object().cloned()).collect()
            }
            _ => return invalid("manifest.json has no well-formed sections array"),
        };
        let zip = reader.zip();
        let start = zip.central_directory_offset();
        if let Some(e) = zip.entries().iter().find(|e| zip.data_end(e).map_or(true, |end| end > start)) {
            return invalid(format!("ZIP entry {:?} does not lie before the central directory; refusing to append", e.name));
        }
        let mut root = Map::new();
        for (k, v) in &manifest {
            if k == "qvf_version" || k == "source" || k == "sections" || (k == "schema_uri" && v2) {
                continue;
            }
            root.insert(k.clone(), v.clone());
        }
        let mut taken: HashSet<String> = zip.entries().iter().map(|e| e.name.clone()).collect();
        let mut info = HashMap::new();
        for s in &sections {
            for m in s.get("members").and_then(Value::as_object).into_iter().flat_map(|m| m.values()) {
                if let Some(p) = m.get("path").and_then(Value::as_str) {
                    taken.insert(p.to_owned());
                }
            }
            let id = s.get("id").and_then(Value::as_str).unwrap_or_default().to_owned();
            info.entry(id).or_insert_with(|| section_info(&reader, s));
        }
        let mut kept = Vec::new();
        for i in 0..zip.entries().len() {
            kept.extend_from_slice(zip.central_record(i));
        }
        let kept_count = zip.entries().len() as u64;
        let tail = (start, zip.bytes()[start as usize..].to_vec());
        let loaded_modified = manifest.get("metadata").and_then(|m| m.get("modified_utc")).cloned();
        let base = UpdateBase { source, root, sections, info, taken, loaded_modified, start, kept, kept_count };
        drop(reader);

        let working = if options.atomic {
            let (tmp, file) = sibling(&target, ".qvf-update-")?;
            drop(file);
            if let Err(e) = std::fs::copy(&target, &tmp) {
                let _ = std::fs::remove_file(&tmp);
                return Err(e.into());
            }
            tmp
        } else {
            target.clone()
        };
        let opened = OpenOptions::new().read(true).write(true).open(&working).and_then(|mut f| {
            f.seek(SeekFrom::Start(start))?;
            Ok(f)
        });
        let file = match opened {
            Ok(f) => f,
            Err(e) => {
                if options.atomic {
                    let _ = std::fs::remove_file(&working);
                }
                return Err(e.into());
            }
        };
        let writer_options =
            WriterOptions { compression: options.compression, force_zip64: options.force_zip64, created_utc: false };
        let writer = Writer::for_update(BufWriter::new(file), base, writer_options);
        Ok(Editor { writer: Some(writer), target, working, atomic: options.atomic, tail })
    }

    /// The writer of the update, with every member method and builder.
    /// Sections already in the archive can be referenced by the builders.
    pub fn writer(&mut self) -> &mut Writer<BufWriter<File>> {
        self.writer.as_mut().expect("an editor holds its writer until commit or abort")
    }

    /// Remove a section and return it. Its members stay in the archive as
    /// orphans until a repack. References to it from other sections dangle
    /// until a section with that id is added again. `job.spec` sections are
    /// refused (spec § 5.9).
    pub fn remove_section(&mut self, id: &str) -> Result<Map<String, Value>, WriteError> {
        self.writer().remove_section(id)
    }

    /// Remove one member of a section and return its spec. Its entry becomes
    /// an orphan.
    pub fn remove_member(&mut self, id: &str, role: &str) -> Result<Map<String, Value>, WriteError> {
        self.writer().remove_member(id, role)
    }

    fn old_path(&mut self, id: &str, role: &str) -> Result<String, WriteError> {
        let old = self.remove_member(id, role)?;
        match old.get("path").and_then(Value::as_str) {
            Some(p) => Ok(base_path(p)),
            None => invalid(format!("section {id:?} member {role:?} has no path")),
        }
    }

    /// Replace a JSON member. The new content is written under a fresh path
    /// derived from the old one (spec § 2.3).
    pub fn replace_json(&mut self, id: &str, role: &str, payload: &Value) -> Result<Value, WriteError> {
        let path = self.old_path(id, role)?;
        self.writer().add_json_at(id, role, payload, &path)
    }

    /// Replace a tensor member (little-endian, C order, `shape`).
    pub fn replace_tensor<T: Element>(&mut self, id: &str, role: &str, data: &[T], shape: &[usize]) -> Result<Value, WriteError> {
        let expected = crate::writer::tensor_bytes(T::DTYPE, shape);
        if expected != Some((data.len() * T::DTYPE.itemsize()) as u64) {
            return invalid(format!("E_BINARY_SIZE: {} elements do not fill shape {shape:?} of member {id}/{role}", data.len()));
        }
        let path = self.old_path(id, role)?;
        let mut sink = self.writer().begin_tensor_at(id, role, T::DTYPE, shape, &path)?;
        sink.write(data)?;
        sink.finish()
    }

    /// Replace an opaque member.
    pub fn replace_opaque(&mut self, id: &str, role: &str, data: &[u8], media_type: Option<&str>) -> Result<Value, WriteError> {
        let path = self.old_path(id, role)?;
        self.writer().add_opaque_at(id, role, data, media_type, &path)
    }

    /// Set the lifecycle status (spec § 3.2, § 5.8): `provenance.run_status`
    /// of the archive, or with `run_id` the `run_status` of that
    /// `run.record` section. `pending` is a root status only.
    pub fn set_run_status(&mut self, status: &str, run_id: Option<&str>) -> Result<(), WriteError> {
        if !["pending", "running", "converged", "failed"].contains(&status) || (run_id.is_some() && status == "pending") {
            return invalid(format!("invalid run status {status:?}"));
        }
        let w = self.writer();
        match run_id {
            None => {
                let provenance = w.root_mut().entry("provenance").or_insert_with(|| Value::Object(Map::new()));
                match provenance.as_object_mut() {
                    Some(p) => {
                        p.insert("run_status".into(), Value::from(status));
                    }
                    None => return invalid("provenance is not an object"),
                }
            }
            Some(id) => {
                let section = w.sections_mut().iter_mut().find(|s| s.get("id").and_then(Value::as_str) == Some(id));
                let Some(section) = section else { return invalid(format!("no section {id:?}")) };
                if section.get("kind").and_then(Value::as_str) != Some("run.record") {
                    return invalid("run_id must name a run.record section");
                }
                section.insert("run_status".into(), Value::from(status));
            }
        }
        Ok(())
    }

    /// Write the new manifest and central directory and, when atomic,
    /// replace the archive. Returns the archive's path.
    pub fn commit(mut self) -> Result<PathBuf, WriteError> {
        let writer = self.writer.take().expect("an editor holds its writer until commit or abort");
        let result = (|| -> Result<(), WriteError> {
            let sink = writer.finish()?;
            let mut file = sink.into_inner().map_err(|e| e.into_error())?;
            // An update shorter than the old directory leaves stale bytes after it.
            let end = file.stream_position()?;
            file.set_len(end)?;
            file.sync_all()?;
            Ok(())
        })();
        if let Err(e) = result {
            self.rollback();
            return Err(e);
        }
        if self.atomic {
            if let Err(e) = std::fs::rename(&self.working, &self.target) {
                self.rollback();
                return Err(e.into());
            }
        }
        self.atomic = false;
        self.tail.1.clear();
        self.tail.0 = u64::MAX;
        Ok(self.target.clone())
    }

    /// Give up the update: the archive is left as it was.
    pub fn abort(mut self) {
        self.rollback();
    }

    fn rollback(&mut self) {
        if let Some(w) = self.writer.take() {
            drop(w.abandon());
        }
        if self.atomic {
            let _ = std::fs::remove_file(&self.working);
            self.atomic = false;
        } else if self.tail.0 != u64::MAX {
            let (start, tail) = std::mem::take(&mut self.tail);
            self.tail.0 = u64::MAX;
            let restored = OpenOptions::new().write(true).open(&self.target).and_then(|mut f| {
                f.seek(SeekFrom::Start(start))?;
                f.write_all(&tail)?;
                f.set_len(start + tail.len() as u64)?;
                f.sync_all()
            });
            if let Err(e) = restored {
                eprintln!("qvf: could not restore {} after an aborted update: {e}", self.target.display());
            }
        }
    }

    /// Write `source` anew from its current manifest into `destination`:
    /// every referenced member, verified (size, CRC-32 and sha256) while it
    /// is streamed in chunks, in archive order, then the manifest (STORED),
    /// and nothing else, so every orphan and superseded manifest is dropped
    /// (spec § 2.3). `None` keeps each entry's compression method, `Some`
    /// sets one. The output is staged beside `destination` and renamed over
    /// it only when complete, so `destination` may be `source`.
    pub fn repack(
        source: impl AsRef<Path>,
        destination: impl AsRef<Path>,
        compression: Option<Compression>,
    ) -> Result<(), WriteError> {
        Self::repack_with(source, destination, compression, &Limits::default())
    }

    /// [`Editor::repack`] under explicit caps.
    pub fn repack_with(
        source: impl AsRef<Path>,
        destination: impl AsRef<Path>,
        compression: Option<Compression>,
        limits: &Limits,
    ) -> Result<(), WriteError> {
        let bytes = std::fs::read(source.as_ref()).map_err(|e| {
            read_error(QvfError::new(Code::E_ZIP_UNREADABLE, format!("cannot read {}: {e}", source.as_ref().display())))
        })?;
        let reader = Reader::from_bytes(bytes, *limits).map_err(read_error)?;
        let (tmp, file) = sibling(destination.as_ref(), ".qvf-repack-")?;
        let result = repack_into(&reader, file, compression);
        match result.and_then(|()| std::fs::rename(&tmp, destination.as_ref()).map_err(WriteError::from)) {
            Ok(()) => Ok(()),
            Err(e) => {
                let _ = std::fs::remove_file(&tmp);
                Err(e)
            }
        }
    }
}

impl Drop for Editor {
    fn drop(&mut self) {
        if self.writer.is_some() {
            self.rollback();
        }
    }
}

/// The member specs the manifest references, by path. Two different
/// declarations of one path are refused.
/// Member specs by path, in manifest order.
type Members = Vec<(String, Map<String, Value>)>;

pub(crate) fn referenced_members(manifest: &Value) -> Result<Members, WriteError> {
    let mut order = Vec::new();
    let mut members: HashMap<String, Map<String, Value>> = HashMap::new();
    for s in manifest.get("sections").and_then(Value::as_array).into_iter().flatten() {
        for m in s.get("members").and_then(Value::as_object).into_iter().flat_map(|m| m.values()) {
            let Some(spec) = m.as_object() else { continue };
            let Some(path) = spec.get("path").and_then(Value::as_str) else { continue };
            if path == "manifest.json" {
                return invalid("a member names manifest.json");
            }
            match members.get(path) {
                Some(seen) if seen != spec => {
                    return invalid(format!("inconsistent declarations for the shared member {path:?}"))
                }
                Some(_) => {}
                None => {
                    members.insert(path.to_owned(), spec.clone());
                    order.push(path.to_owned());
                }
            }
        }
    }
    Ok(order
        .into_iter()
        .map(|p| {
            let spec = members.remove(&p).unwrap_or_default();
            (p, spec)
        })
        .collect())
}

/// Stream one referenced member out of `reader`, verified against its spec,
/// to `sink`. The caller discards the output when this fails.
pub(crate) fn copy_member(
    reader: &Reader,
    path: &str,
    spec: &Map<String, Value>,
    mut sink: impl FnMut(&[u8]) -> std::io::Result<()>,
) -> Result<(), WriteError> {
    let zip = reader.zip();
    let limits = reader.limits();
    let entry = zip.find(path).ok_or_else(|| {
        read_error(QvfError::new(Code::E_MEMBER_MISSING, format!("member {path:?} declared but not in archive")))
    })?;
    if zip.total_size() as f64 > limits.max_total_bytes {
        return Err(read_error(QvfError::new(
            Code::E_TOTAL_TOO_LARGE,
            format!("archive expands to {} bytes, over the {}-byte total cap", zip.total_size(), limits.max_total_bytes),
        )));
    }
    check_entry_caps(entry, limits).map_err(read_error)?;
    if let Some((code, msg)) = declared_size_problem(spec, entry.size, path) {
        return Err(read_error(QvfError::new(code, msg)));
    }
    let mut hash = Sha256::new();
    let mut io_error = None;
    let streamed = zip.stream(entry, Some(limits.max_member_bytes), |chunk| {
        hash.update(chunk);
        sink(chunk).map_err(|e| {
            io_error = Some(e);
            QvfError::new(Code::E_ZIP_UNREADABLE, "write failed")
        })
    });
    if let Some(e) = io_error {
        return Err(e.into());
    }
    streamed.map_err(read_error)?;
    if let Some(declared) = spec.get("sha256") {
        let actual: String = hash.finalize().iter().map(|b| format!("{b:02x}")).collect();
        if declared.as_str() != Some(actual.as_str()) {
            return Err(read_error(QvfError::new(Code::E_SHA256_MISMATCH, format!("sha256 mismatch for {path}"))));
        }
    }
    Ok(())
}

fn repack_into(reader: &Reader, file: File, compression: Option<Compression>) -> Result<(), WriteError> {
    let zip = reader.zip();
    let mut wanted: HashMap<String, Map<String, Value>> = referenced_members(reader.manifest())?.into_iter().collect();
    let mut out = ZipWriter::new(BufWriter::new(file), false);
    for entry in zip.entries() {
        let Some(spec) = wanted.remove(&entry.name) else { continue };
        let method = match (compression, entry.method) {
            (Some(Compression::Stored), _) | (None, METHOD_STORE) => Method::Store,
            (Some(Compression::Deflate), _) | (None, METHOD_DEFLATE) => Method::Deflate,
            (None, m) => return invalid(format!("ZIP entry {:?}: unsupported compression method {m}", entry.name)),
        };
        out.start_entry(&entry.name, method, Some(entry.size))?;
        copy_member(reader, &entry.name, &spec, |chunk| out.write_data(chunk))?;
        out.finish_entry()?;
    }
    if !wanted.is_empty() {
        let mut missing: Vec<&String> = wanted.keys().collect();
        missing.sort();
        return Err(read_error(QvfError::new(Code::E_MEMBER_MISSING, format!("members not in the archive: {missing:?}"))));
    }
    let manifest = zip.find("manifest.json").expect("the reader found it");
    let raw = zip.read(manifest, Some(reader.limits().max_member_bytes)).map_err(read_error)?;
    out.start_entry("manifest.json", Method::Store, Some(raw.len() as u64))?;
    out.write_data(&raw)?;
    out.finish_entry()?;
    let sink = out.finish()?;
    let file = sink.into_inner().map_err(|e| e.into_error())?;
    file.sync_all()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn replacement_paths_drop_the_update_counter() {
        assert_eq!(super::base_path("structure/structure.2.json"), "structure/structure.json");
        assert_eq!(super::base_path("structure/structure.json"), "structure/structure.json");
        assert_eq!(super::base_path("volumes/rho.12.dat"), "volumes/rho.dat");
        assert_eq!(super::base_path("a.b/c.v2.json"), "a.b/c.v2.json");
        assert_eq!(super::base_path("x.3"), "x");
    }
}
