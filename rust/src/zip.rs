//! Container layer: a minimal, defensive ZIP reader (spec § 2.1).
//!
//! Reads the central directory of an archive held in memory, including ZIP64
//! (archives over 4 GiB, more than 65535 entries), and extracts STORE and
//! DEFLATE entries with CRC-32 verification. The archive is untrusted input:
//! every offset is bounds-checked, duplicate entry names are recorded rather
//! than silently resolved, and an entry is never inflated past the
//! uncompressed size its central-directory record declares. It is a port of
//! `js/src/zip.ts`.

use std::collections::HashMap;
use std::io::Read;

use crate::error::{Code, QvfError};
use crate::limits::Limits;

const SIG_EOCD: u32 = 0x0605_4b50;
const SIG_EOCD64: u32 = 0x0606_4b50;
const SIG_LOCATOR64: u32 = 0x0706_4b50;
const SIG_CENTRAL: u32 = 0x0201_4b50;
const SIG_LOCAL: u32 = 0x0403_4b50;

const EOCD_SIZE: usize = 22;
const CENTRAL_HEADER_SIZE: usize = 46;
const LOCAL_HEADER_SIZE: usize = 30;
const MAX_COMMENT: usize = 0xffff;
const SATURATED16: u64 = 0xffff;
const SATURATED32: u64 = 0xffff_ffff;

/// Compression method STORE.
pub const METHOD_STORE: u16 = 0;
/// Compression method DEFLATE.
pub const METHOD_DEFLATE: u16 = 8;

/// One central-directory record.
#[derive(Debug, Clone)]
pub struct ZipEntry {
    /// The entry name as stored (UTF-8 when it decodes as such, else Latin-1).
    pub name: String,
    /// Compression method: 0 = STORE, 8 = DEFLATE.
    pub method: u16,
    /// General-purpose bit flags.
    pub flags: u16,
    /// Declared CRC-32 of the uncompressed data.
    pub crc32: u32,
    /// Declared compressed size in bytes.
    pub compressed_size: u64,
    /// Declared uncompressed size in bytes.
    pub size: u64,
    /// Offset of the local header from the start of the archive.
    pub local_offset: u64,
    /// The raw name bytes, compared against the local header's copy on read.
    pub name_bytes: Vec<u8>,
}

fn unreadable(message: impl Into<String>) -> QvfError {
    QvfError::new(Code::E_ZIP_UNREADABLE, message)
}

/// CRC-32 (IEEE 802.3, the ZIP checksum) of `data`.
pub fn crc32(data: &[u8]) -> u32 {
    let mut crc = flate2::Crc::new();
    crc.update(data);
    crc.sum()
}

/// Why a member path or ZIP entry name is illegal (spec § 2.1), or `None`.
///
/// Illegal: a backslash, a NUL, a leading `/`, a drive-letter prefix, or a
/// `..` path segment.
pub fn path_problem(path: &str) -> Option<&'static str> {
    let b = path.as_bytes();
    if path.contains('\\') {
        Some("contains a backslash")
    } else if path.contains('\0') {
        Some("contains a NUL character")
    } else if path.starts_with('/') {
        Some("begins with '/'")
    } else if b.len() >= 2 && b[0].is_ascii_alphabetic() && b[1] == b':' {
        Some("begins with a drive letter")
    } else if path.split('/').any(|seg| seg == "..") {
        Some("contains a '..' segment")
    } else {
        None
    }
}

fn decode_name(raw: &[u8], flags: u16) -> String {
    if flags & 0x800 != 0 {
        return String::from_utf8_lossy(raw).into_owned();
    }
    match std::str::from_utf8(raw) {
        Ok(s) => s.to_owned(),
        Err(_) => raw.iter().map(|&b| char::from(b)).collect(),
    }
}

struct View<'a>(&'a [u8]);

impl View<'_> {
    fn u16(&self, at: usize) -> Result<u16, QvfError> {
        let s = self.0.get(at..at + 2).ok_or_else(|| unreadable("read past the end of the archive"))?;
        Ok(u16::from_le_bytes([s[0], s[1]]))
    }
    fn u32(&self, at: usize) -> Result<u32, QvfError> {
        let s = self.0.get(at..at + 4).ok_or_else(|| unreadable("read past the end of the archive"))?;
        Ok(u32::from_le_bytes([s[0], s[1], s[2], s[3]]))
    }
    fn u64(&self, at: usize) -> Result<u64, QvfError> {
        let s = self.0.get(at..at + 8).ok_or_else(|| unreadable("read past the end of the archive"))?;
        let mut b = [0u8; 8];
        b.copy_from_slice(s);
        Ok(u64::from_le_bytes(b))
    }
}

fn to_usize(v: u64, what: &str) -> Result<usize, QvfError> {
    usize::try_from(v).map_err(|_| unreadable(format!("{what} does not fit in memory")))
}

fn find_eocd(view: &View, length: usize) -> Result<usize, QvfError> {
    if length < EOCD_SIZE {
        return Err(unreadable("too small to be a ZIP archive"));
    }
    let lowest = length.saturating_sub(EOCD_SIZE + MAX_COMMENT);
    let mut p = length - EOCD_SIZE;
    loop {
        if view.u32(p)? == SIG_EOCD {
            let comment = usize::from(view.u16(p + 20)?);
            if p + EOCD_SIZE + comment <= length {
                return Ok(p);
            }
        }
        if p == lowest {
            break;
        }
        p -= 1;
    }
    Err(unreadable("no end-of-central-directory record (not a ZIP archive)"))
}

struct Directory {
    count: u64,
    size: u64,
    offset: u64,
}

fn read_directory_location(view: &View, length: usize, eocd: usize) -> Result<Directory, QvfError> {
    let mut count = u64::from(view.u16(eocd + 10)?);
    let mut size = u64::from(view.u32(eocd + 12)?);
    let mut offset = u64::from(view.u32(eocd + 16)?);
    let disk = view.u16(eocd + 4)?;
    let directory_disk = view.u16(eocd + 6)?;
    let saturated = count == SATURATED16 || size == SATURATED32 || offset == SATURATED32;
    let locator = eocd.checked_sub(20);
    let has_locator = match locator {
        Some(l) => view.u32(l)? == SIG_LOCATOR64,
        None => false,
    };
    if saturated && has_locator {
        let locator = locator.unwrap_or(0);
        let record = to_usize(view.u64(locator + 8)?, "the ZIP64 end record offset")?;
        if record.checked_add(56).map_or(true, |end| end > length) || view.u32(record)? != SIG_EOCD64 {
            return Err(unreadable("invalid ZIP64 end-of-central-directory record"));
        }
        if view.u32(locator + 4)? != 0 || view.u32(locator + 16)? != 1 {
            return Err(unreadable("multi-volume ZIP archives are not supported"));
        }
        if view.u32(record + 16)? != 0 || view.u32(record + 20)? != 0 {
            return Err(unreadable("multi-volume ZIP archives are not supported"));
        }
        if count == SATURATED16 {
            count = view.u64(record + 32)?;
        }
        if size == SATURATED32 {
            size = view.u64(record + 40)?;
        }
        if offset == SATURATED32 {
            offset = view.u64(record + 48)?;
        }
    } else if disk != 0 || directory_disk != 0 {
        return Err(unreadable("multi-volume ZIP archives are not supported"));
    }
    Ok(Directory { count, size, offset })
}

/// Resolve the 0x0001 (ZIP64) extra field of a central-directory record.
fn zip64_fields(view: &View, extra_start: usize, extra_length: usize, wanted: [bool; 3]) -> Result<[Option<u64>; 3], QvfError> {
    let end = extra_start + extra_length;
    let mut q = extra_start;
    while q + 4 <= end {
        let id = view.u16(q)?;
        let length = usize::from(view.u16(q + 2)?);
        let body = q + 4;
        if body + length > end {
            break;
        }
        if id == 0x0001 {
            let mut out = [None; 3];
            let mut o = body;
            for (k, want) in wanted.iter().enumerate() {
                if !*want {
                    continue;
                }
                if o + 8 > body + length {
                    return Err(unreadable("truncated ZIP64 extra field"));
                }
                out[k] = Some(view.u64(o)?);
                o += 8;
            }
            return Ok(out);
        }
        q = body + length;
    }
    Err(unreadable("a size or offset is saturated but there is no ZIP64 extra field"))
}

fn read_central_directory(bytes: &[u8], view: &View, dir: &Directory) -> Result<Vec<ZipEntry>, QvfError> {
    let length = bytes.len() as u64;
    if dir.offset > length || dir.size > length - dir.offset {
        return Err(unreadable("central directory out of range"));
    }
    if dir.count > dir.size / CENTRAL_HEADER_SIZE as u64 {
        return Err(unreadable("central directory too small for its entry count"));
    }
    let length = bytes.len();
    let mut entries = Vec::with_capacity(to_usize(dir.count, "the entry count")?);
    let mut p = to_usize(dir.offset, "the central directory offset")?;
    for _ in 0..dir.count {
        if p + CENTRAL_HEADER_SIZE > length || view.u32(p)? != SIG_CENTRAL {
            return Err(unreadable("bad central directory entry"));
        }
        let flags = view.u16(p + 8)?;
        let method = view.u16(p + 10)?;
        let crc = view.u32(p + 16)?;
        let mut compressed_size = u64::from(view.u32(p + 20)?);
        let mut size = u64::from(view.u32(p + 24)?);
        let name_length = usize::from(view.u16(p + 28)?);
        let extra_length = usize::from(view.u16(p + 30)?);
        let comment_length = usize::from(view.u16(p + 32)?);
        let mut local_offset = u64::from(view.u32(p + 42)?);
        let name_start = p + CENTRAL_HEADER_SIZE;
        let next = name_start + name_length + extra_length + comment_length;
        if next > length {
            return Err(unreadable("truncated central directory entry"));
        }
        let wanted = [size == SATURATED32, compressed_size == SATURATED32, local_offset == SATURATED32];
        if wanted.iter().any(|w| *w) {
            let z = zip64_fields(view, name_start + name_length, extra_length, wanted)?;
            if let Some(v) = z[0] {
                size = v;
            }
            if let Some(v) = z[1] {
                compressed_size = v;
            }
            if let Some(v) = z[2] {
                local_offset = v;
            }
        }
        let name_bytes = bytes[name_start..name_start + name_length].to_vec();
        entries.push(ZipEntry {
            name: decode_name(&name_bytes, flags),
            method,
            flags,
            crc32: crc,
            compressed_size,
            size,
            local_offset,
            name_bytes,
        });
        p = next;
    }
    Ok(entries)
}

/// Inflate `data` into exactly `size` bytes. The output is bounded: a stream
/// that expands past `size` is abandoned after one byte more, and memory grows
/// with what the stream really produces, not with the declared size.
fn inflate_exact(data: &[u8], size: u64, name: &str) -> Result<Vec<u8>, QvfError> {
    let mut out = Vec::new();
    let mut decoder = flate2::read::DeflateDecoder::new(data).take(size.saturating_add(1));
    decoder.read_to_end(&mut out).map_err(|e| unreadable(format!("ZIP entry {name:?} does not inflate ({e})")))?;
    let written = out.len() as u64;
    if written > size {
        return Err(unreadable(format!("ZIP entry {name:?} inflates to more than its declared {size} bytes")));
    }
    if written != size {
        return Err(unreadable(format!("ZIP entry {name:?} inflates to {written} bytes, not its declared {size}")));
    }
    Ok(out)
}

/// Whether a compressed entry expands past the ratio cap (spec § 6.1).
pub fn ratio_exceeded(entry: &ZipEntry, limits: &Limits) -> bool {
    limits.max_ratio > 0.0
        && entry.method != METHOD_STORE
        && (entry.size as f64) >= limits.ratio_exempt_bytes
        && (entry.compressed_size == 0 || (entry.size as f64) > limits.max_ratio * entry.compressed_size as f64)
}

/// A ZIP archive held in memory, opened through its central directory.
#[derive(Debug)]
pub struct ZipArchive {
    bytes: Vec<u8>,
    entries: Vec<ZipEntry>,
    /// name -> index of the last record of that name (like Python's zipfile)
    index: HashMap<String, usize>,
    duplicates: Vec<String>,
    total_size: u64,
}

impl ZipArchive {
    /// Open an archive. Fails with `E_ZIP_UNREADABLE` when the bytes are not a
    /// ZIP archive or its central directory is damaged. Entries are not
    /// decompressed here.
    pub fn open(bytes: Vec<u8>) -> Result<Self, QvfError> {
        let view = View(&bytes);
        let eocd = find_eocd(&view, bytes.len())?;
        let dir = read_directory_location(&view, bytes.len(), eocd)?;
        let entries = read_central_directory(&bytes, &view, &dir)?;
        let mut index = HashMap::new();
        let mut duplicates: Vec<String> = Vec::new();
        let mut total: u64 = 0;
        for (i, e) in entries.iter().enumerate() {
            if index.contains_key(&e.name) && !duplicates.contains(&e.name) {
                duplicates.push(e.name.clone());
            }
            index.insert(e.name.clone(), i);
            total = total.saturating_add(e.size);
        }
        Ok(ZipArchive { bytes, entries, index, duplicates, total_size: total })
    }

    /// Every central-directory record, in archive order, duplicates included.
    pub fn entries(&self) -> &[ZipEntry] {
        &self.entries
    }

    /// Names that occur more than once, in order of their second occurrence.
    pub fn duplicates(&self) -> &[String] {
        &self.duplicates
    }

    /// Sum of the declared uncompressed sizes of all entries.
    pub fn total_size(&self) -> u64 {
        self.total_size
    }

    /// The record a lookup by name resolves to (the last one when duplicated).
    pub fn find(&self, name: &str) -> Option<&ZipEntry> {
        self.index.get(name).map(|&i| &self.entries[i])
    }

    /// Unique entry names in order of first occurrence.
    pub fn names(&self) -> Vec<&str> {
        let mut seen = std::collections::HashSet::new();
        self.entries.iter().map(|e| e.name.as_str()).filter(|n| seen.insert(*n)).collect()
    }

    /// The uncompressed bytes of an entry.
    ///
    /// Fails with `E_MEMBER_TOO_LARGE` when the declared size exceeds
    /// `max_bytes` (nothing is decompressed) and with `E_ZIP_UNREADABLE` when
    /// the data cannot be read (truncated, encrypted, unsupported method, does
    /// not inflate to exactly the declared size, CRC-32 mismatch).
    pub fn read(&self, entry: &ZipEntry, max_bytes: Option<f64>) -> Result<Vec<u8>, QvfError> {
        let label = format!("{:?}", entry.name);
        if let Some(max) = max_bytes {
            if entry.size as f64 > max {
                return Err(QvfError::new(
                    Code::E_MEMBER_TOO_LARGE,
                    format!("ZIP entry {label} declares {} bytes, over the {max}-byte cap", entry.size),
                ));
            }
        }
        if entry.flags & 0x1 != 0 || entry.flags & 0x40 != 0 {
            return Err(unreadable(format!("ZIP entry {label} is encrypted")));
        }
        let bytes = &self.bytes;
        let view = View(bytes);
        let offset = usize::try_from(entry.local_offset).unwrap_or(usize::MAX);
        if offset.checked_add(LOCAL_HEADER_SIZE).map_or(true, |end| end > bytes.len()) {
            return Err(unreadable(format!("ZIP entry {label}: truncated local header")));
        }
        if view.u32(offset)? != SIG_LOCAL {
            return Err(unreadable(format!("ZIP entry {label}: bad local header signature")));
        }
        let name_length = usize::from(view.u16(offset + 26)?);
        let extra_length = usize::from(view.u16(offset + 28)?);
        let name_start = offset + LOCAL_HEADER_SIZE;
        let data_start = name_start + name_length + extra_length;
        if data_start > bytes.len() || entry.compressed_size > (bytes.len() - data_start) as u64 {
            return Err(unreadable(format!("ZIP entry {label}: truncated data")));
        }
        if bytes[name_start..name_start + name_length] != entry.name_bytes[..] {
            return Err(unreadable(format!("ZIP entry {label}: local header name differs from the directory")));
        }
        let data = &bytes[data_start..data_start + entry.compressed_size as usize];
        let out = match entry.method {
            METHOD_STORE => {
                if entry.compressed_size != entry.size {
                    return Err(unreadable(format!("ZIP entry {label}: stored size differs from uncompressed size")));
                }
                data.to_vec()
            }
            METHOD_DEFLATE => inflate_exact(data, entry.size, &entry.name)?,
            m => return Err(unreadable(format!("ZIP entry {label}: unsupported compression method {m}"))),
        };
        if crc32(&out) != entry.crc32 {
            return Err(unreadable(format!("ZIP entry {label}: CRC-32 mismatch")));
        }
        Ok(out)
    }
}
