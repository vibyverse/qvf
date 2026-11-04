//! Container layer for writing: a streaming ZIP writer (spec § 2.1).
//!
//! Entries are written front to back to any [`Write`] sink, which need not
//! be seekable. Every local header sets the data-descriptor flag, so the
//! CRC-32 and sizes follow the data, and the central directory at the end
//! carries the authoritative values. An entry whose size can reach 4 GiB, an
//! offset past 4 GiB or more than 65535 entries switch to ZIP64 fields and
//! end records automatically. `force_zip64` writes them for every entry,
//! which tests the ZIP64 path of readers without gigabytes of data.
//!
//! Entry names are UTF-8 (general-purpose bit 11). Methods are STORE and raw
//! DEFLATE (`flate2`, pure Rust backend). Nothing of an entry is retained
//! after it is written except its central-directory record.

use std::io::{self, Write};

use flate2::{Compress, Compression as Level, Crc, FlushCompress, Status};

use crate::zip::{METHOD_DEFLATE, METHOD_STORE};

const SIG_LOCAL: u32 = 0x0403_4b50;
const SIG_DESCRIPTOR: u32 = 0x0807_4b50;
const SIG_CENTRAL: u32 = 0x0201_4b50;
const SIG_EOCD: u32 = 0x0605_4b50;
const SIG_EOCD64: u32 = 0x0606_4b50;
const SIG_LOCATOR64: u32 = 0x0706_4b50;

const FLAG_DESCRIPTOR: u16 = 0x0008;
const FLAG_UTF8: u16 = 0x0800;
const VERSION_DEFAULT: u16 = 20;
const VERSION_ZIP64: u16 = 45;
/// Made by: UNIX (3), APPNOTE 4.5.
const VERSION_MADE_BY: u16 = (3 << 8) | VERSION_ZIP64;
/// A regular file, mode 0644.
const EXTERNAL_ATTR: u32 = 0o100_644 << 16;

const MAX32: u64 = 0xffff_ffff;
const MAX16: u64 = 0xffff;
const OUT_CHUNK: usize = 64 * 1024;

/// The ZIP compression method of an entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// No compression (method 0).
    Store,
    /// Raw DEFLATE (method 8).
    Deflate,
}

impl Method {
    fn code(self) -> u16 {
        match self {
            Method::Store => METHOD_STORE,
            Method::Deflate => METHOD_DEFLATE,
        }
    }
}

/// A finished entry, as its central-directory record describes it.
#[derive(Debug, Clone)]
struct Record {
    name: Vec<u8>,
    method: Method,
    crc32: u32,
    compressed_size: u64,
    size: u64,
    offset: u64,
    zip64: bool,
}

struct Active {
    record: Record,
    crc: Crc,
    deflate: Option<Compress>,
}

/// A sink that counts the bytes written through it.
struct Counting<W> {
    inner: W,
    pos: u64,
}

impl<W: Write> Counting<W> {
    fn put(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.inner.write_all(bytes)?;
        self.pos += bytes.len() as u64;
        Ok(())
    }
}

/// A streaming ZIP writer. Start an entry, write its data in any number of
/// pieces, finish it, and call [`ZipWriter::finish`] after the last one.
pub struct ZipWriter<W: Write> {
    out: Counting<W>,
    records: Vec<Record>,
    active: Option<Active>,
    force_zip64: bool,
    dos_time: u16,
    dos_date: u16,
}

/// The MS-DOS time and date fields of a UTC time stamp.
fn dos_stamp(t: crate::writer::UtcTime) -> (u16, u16) {
    let year = t.year.clamp(1980, 2107) - 1980;
    let time = (t.hour << 11) | (t.minute << 5) | (t.second / 2);
    let date = ((year as u16) << 9) | (t.month << 5) | t.day;
    (time, date)
}

impl<W: Write> ZipWriter<W> {
    /// A writer on `sink`. Entries get the current UTC time as their
    /// modification time.
    pub fn new(sink: W, force_zip64: bool) -> Self {
        let (dos_time, dos_date) = dos_stamp(crate::writer::UtcTime::now());
        ZipWriter { out: Counting { inner: sink, pos: 0 }, records: Vec::new(), active: None, force_zip64, dos_time, dos_date }
    }

    /// Bytes written so far.
    pub fn position(&self) -> u64 {
        self.out.pos
    }

    /// Number of finished entries.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether no entry has been finished yet.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Begin an entry. `size_hint` is its uncompressed size when known, which
    /// decides whether its local header carries ZIP64 fields. An entry
    /// without a hint that grows past 4 GiB fails at [`ZipWriter::finish_entry`].
    pub fn start_entry(&mut self, name: &str, method: Method, size_hint: Option<u64>) -> io::Result<()> {
        if self.active.is_some() {
            return Err(io::Error::other("a ZIP entry is still open"));
        }
        if name.len() > MAX16 as usize {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "ZIP entry name over 65535 bytes"));
        }
        // DEFLATE can expand incompressible data slightly; leave a margin.
        let bound = size_hint.map(|s| s.saturating_add(s / 1000).saturating_add(1024));
        let zip64 = self.force_zip64 || matches!(bound, Some(b) if b >= MAX32);
        let offset = self.out.pos;
        let mut header = Vec::with_capacity(30 + name.len() + 20);
        put_u32(&mut header, SIG_LOCAL);
        put_u16(&mut header, if zip64 { VERSION_ZIP64 } else { VERSION_DEFAULT });
        put_u16(&mut header, FLAG_DESCRIPTOR | FLAG_UTF8);
        put_u16(&mut header, method.code());
        put_u16(&mut header, self.dos_time);
        put_u16(&mut header, self.dos_date);
        put_u32(&mut header, 0); // CRC-32 in the data descriptor
        if zip64 {
            put_u32(&mut header, MAX32 as u32);
            put_u32(&mut header, MAX32 as u32);
        } else {
            put_u32(&mut header, 0);
            put_u32(&mut header, 0);
        }
        put_u16(&mut header, name.len() as u16);
        put_u16(&mut header, if zip64 { 20 } else { 0 });
        header.extend_from_slice(name.as_bytes());
        if zip64 {
            put_u16(&mut header, 0x0001);
            put_u16(&mut header, 16);
            put_u64(&mut header, 0); // sizes in the data descriptor
            put_u64(&mut header, 0);
        }
        self.out.put(&header)?;
        let deflate = match method {
            Method::Deflate => Some(Compress::new(Level::default(), false)),
            Method::Store => None,
        };
        self.active = Some(Active {
            record: Record { name: name.as_bytes().to_vec(), method, crc32: 0, compressed_size: 0, size: 0, offset, zip64 },
            crc: Crc::new(),
            deflate,
        });
        Ok(())
    }

    /// Append data to the open entry.
    pub fn write_data(&mut self, data: &[u8]) -> io::Result<()> {
        let active = self.active.as_mut().ok_or_else(|| io::Error::other("no ZIP entry is open"))?;
        active.crc.update(data);
        active.record.size += data.len() as u64;
        match active.deflate.as_mut() {
            None => {
                self.out.put(data)?;
                active.record.compressed_size += data.len() as u64;
            }
            Some(c) => {
                let n = deflate_into(c, data, FlushCompress::None, &mut self.out)?;
                active.record.compressed_size += n;
            }
        }
        Ok(())
    }

    /// Finish the open entry: flush the compressor and write its data
    /// descriptor. Returns the uncompressed size.
    pub fn finish_entry(&mut self) -> io::Result<u64> {
        let mut active = self.active.take().ok_or_else(|| io::Error::other("no ZIP entry is open"))?;
        if let Some(c) = active.deflate.as_mut() {
            active.record.compressed_size += deflate_into(c, &[], FlushCompress::Finish, &mut self.out)?;
        }
        let rec = &mut active.record;
        rec.crc32 = active.crc.sum();
        if !rec.zip64 && (rec.size >= MAX32 || rec.compressed_size >= MAX32) {
            return Err(io::Error::other(format!(
                "ZIP entry {:?} reached 4 GiB without ZIP64 fields; give its size when starting it",
                String::from_utf8_lossy(&rec.name)
            )));
        }
        let mut desc = Vec::with_capacity(24);
        put_u32(&mut desc, SIG_DESCRIPTOR);
        put_u32(&mut desc, rec.crc32);
        if rec.zip64 {
            put_u64(&mut desc, rec.compressed_size);
            put_u64(&mut desc, rec.size);
        } else {
            put_u32(&mut desc, rec.compressed_size as u32);
            put_u32(&mut desc, rec.size as u32);
        }
        self.out.put(&desc)?;
        let size = rec.size;
        self.records.push(active.record);
        Ok(size)
    }

    /// Abandon the open entry, if any. Its bytes stay in the sink, but no
    /// central-directory record names them.
    pub fn abandon_entry(&mut self) {
        self.active = None;
    }

    /// Write the central directory and end records and return the sink.
    pub fn finish(mut self) -> io::Result<W> {
        if self.active.is_some() {
            return Err(io::Error::other("a ZIP entry is still open"));
        }
        let cd_offset = self.out.pos;
        for rec in &self.records {
            let big_size = self.force_zip64 || rec.size >= MAX32 || rec.compressed_size >= MAX32;
            let big_offset = self.force_zip64 || rec.offset >= MAX32;
            let mut extra = Vec::new();
            if big_size || big_offset {
                let mut body = Vec::new();
                if big_size {
                    put_u64(&mut body, rec.size);
                    put_u64(&mut body, rec.compressed_size);
                }
                if big_offset {
                    put_u64(&mut body, rec.offset);
                }
                put_u16(&mut extra, 0x0001);
                put_u16(&mut extra, body.len() as u16);
                extra.extend_from_slice(&body);
            }
            let needs64 = rec.zip64 || !extra.is_empty();
            let mut h = Vec::with_capacity(46 + rec.name.len() + extra.len());
            put_u32(&mut h, SIG_CENTRAL);
            put_u16(&mut h, VERSION_MADE_BY);
            put_u16(&mut h, if needs64 { VERSION_ZIP64 } else { VERSION_DEFAULT });
            put_u16(&mut h, FLAG_DESCRIPTOR | FLAG_UTF8);
            put_u16(&mut h, rec.method.code());
            put_u16(&mut h, self.dos_time);
            put_u16(&mut h, self.dos_date);
            put_u32(&mut h, rec.crc32);
            put_u32(&mut h, if big_size { MAX32 as u32 } else { rec.compressed_size as u32 });
            put_u32(&mut h, if big_size { MAX32 as u32 } else { rec.size as u32 });
            put_u16(&mut h, rec.name.len() as u16);
            put_u16(&mut h, extra.len() as u16);
            put_u16(&mut h, 0); // comment
            put_u16(&mut h, 0); // disk
            put_u16(&mut h, 0); // internal attributes
            put_u32(&mut h, EXTERNAL_ATTR);
            put_u32(&mut h, if big_offset { MAX32 as u32 } else { rec.offset as u32 });
            h.extend_from_slice(&rec.name);
            h.extend_from_slice(&extra);
            self.out.put(&h)?;
        }
        let cd_size = self.out.pos - cd_offset;
        let count = self.records.len() as u64;
        let end64 = self.force_zip64 || count >= MAX16 || cd_offset >= MAX32 || cd_size >= MAX32;
        let mut tail = Vec::with_capacity(98);
        if end64 {
            let eocd64_offset = self.out.pos;
            put_u32(&mut tail, SIG_EOCD64);
            put_u64(&mut tail, 44);
            put_u16(&mut tail, VERSION_MADE_BY);
            put_u16(&mut tail, VERSION_ZIP64);
            put_u32(&mut tail, 0);
            put_u32(&mut tail, 0);
            put_u64(&mut tail, count);
            put_u64(&mut tail, count);
            put_u64(&mut tail, cd_size);
            put_u64(&mut tail, cd_offset);
            put_u32(&mut tail, SIG_LOCATOR64);
            put_u32(&mut tail, 0);
            put_u64(&mut tail, eocd64_offset);
            put_u32(&mut tail, 1);
        }
        put_u32(&mut tail, SIG_EOCD);
        put_u16(&mut tail, 0);
        put_u16(&mut tail, 0);
        let count16 = if end64 { MAX16 as u16 } else { count as u16 };
        put_u16(&mut tail, count16);
        put_u16(&mut tail, count16);
        put_u32(&mut tail, if end64 { MAX32 as u32 } else { cd_size as u32 });
        put_u32(&mut tail, if end64 { MAX32 as u32 } else { cd_offset as u32 });
        put_u16(&mut tail, 0);
        self.out.put(&tail)?;
        self.out.inner.flush()?;
        Ok(self.out.inner)
    }
}

/// Run `input` through the compressor and write its output. Returns the
/// number of compressed bytes written.
fn deflate_into<W: Write>(c: &mut Compress, input: &[u8], flush: FlushCompress, out: &mut Counting<W>) -> io::Result<u64> {
    let mut buf: Vec<u8> = Vec::with_capacity(OUT_CHUNK);
    let mut consumed = 0usize;
    let mut written = 0u64;
    loop {
        buf.clear();
        let before = c.total_in();
        let status = c.compress_vec(&input[consumed..], &mut buf, flush).map_err(io::Error::other)?;
        consumed += (c.total_in() - before) as usize;
        out.put(&buf)?;
        written += buf.len() as u64;
        let done = match flush {
            FlushCompress::Finish => status == Status::StreamEnd,
            _ => consumed == input.len() && buf.len() < OUT_CHUNK,
        };
        if done {
            return Ok(written);
        }
    }
}

fn put_u16(v: &mut Vec<u8>, x: u16) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn put_u32(v: &mut Vec<u8>, x: u32) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn put_u64(v: &mut Vec<u8>, x: u64) {
    v.extend_from_slice(&x.to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zip::ZipArchive;

    fn roundtrip(force_zip64: bool) {
        let mut z = ZipWriter::new(Vec::new(), force_zip64);
        let big: Vec<u8> = (0..300_000u32).map(|i| (i % 251) as u8).collect();
        for (name, method, data) in [
            ("a.json", Method::Deflate, &b"{}"[..]),
            ("b/c.dat", Method::Deflate, &big[..]),
            ("manifest.json", Method::Store, b"[]"),
        ] {
            z.start_entry(name, method, Some(data.len() as u64)).unwrap();
            for chunk in data.chunks(7777) {
                z.write_data(chunk).unwrap();
            }
            assert_eq!(z.finish_entry().unwrap(), data.len() as u64);
        }
        let bytes = z.finish().unwrap();
        let archive = ZipArchive::open(bytes).unwrap();
        assert_eq!(archive.names(), vec!["a.json", "b/c.dat", "manifest.json"]);
        assert_eq!(archive.read(archive.find("b/c.dat").unwrap(), None).unwrap(), big);
        assert_eq!(archive.read(archive.find("manifest.json").unwrap(), None).unwrap(), b"[]");
        assert_eq!(archive.find("manifest.json").unwrap().method, METHOD_STORE);
    }

    #[test]
    fn writes_readable_archives() {
        roundtrip(false);
    }

    #[test]
    fn writes_readable_zip64_archives() {
        roundtrip(true);
    }
}
