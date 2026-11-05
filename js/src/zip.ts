/**
 * Container layer: a minimal, defensive ZIP reader (spec § 2.1).
 *
 * Reads the central directory of an archive held in memory, including ZIP64
 * (archives over 4 GiB, more than 65535 entries) and extracts STORE and
 * DEFLATE entries with CRC-32 verification. The archive is untrusted input:
 * every offset is bounds-checked, duplicate entry names are recorded rather
 * than silently resolved, and an entry is never inflated past the
 * uncompressed size its central-directory record declares.
 *
 * Nothing here depends on Node: the module runs in browsers unchanged.
 */

import { Inflate } from "fflate";

import { QvfError } from "./errors.js";
import type { Limits } from "./limits.js";

const SIG_EOCD = 0x06054b50;
const SIG_EOCD64 = 0x06064b50;
const SIG_LOCATOR64 = 0x07064b50;
const SIG_CENTRAL = 0x02014b50;
const SIG_LOCAL = 0x04034b50;

const EOCD_SIZE = 22;
const CENTRAL_HEADER_SIZE = 46;
const LOCAL_HEADER_SIZE = 30;
const MAX_COMMENT = 0xffff;
const SATURATED16 = 0xffff;
const SATURATED32 = 0xffffffff;

/** Largest compressed input handed to the inflater in one step. */
const INFLATE_STEP_MAX = 64 * 1024;
/** DEFLATE cannot expand by more than this factor (RFC 1951 worst case is 1032:1). */
const DEFLATE_MAX_EXPANSION = 1032;
/** Allowance on top of that bound, so that it can never refuse a real stream. */
const DEFLATE_SLACK = 1024;

export const METHOD_STORE = 0;
export const METHOD_DEFLATE = 8;

/** One central-directory record. */
export interface ZipEntry {
  /** The entry name as stored (UTF-8 when it decodes as such, else Latin-1). */
  readonly name: string;
  /** Compression method: 0 = STORE, 8 = DEFLATE. */
  readonly method: number;
  /** General-purpose bit flags. */
  readonly flags: number;
  /** Declared CRC-32 of the uncompressed data. */
  readonly crc32: number;
  /** Declared compressed size in bytes. */
  readonly compressedSize: number;
  /** Declared uncompressed size in bytes. */
  readonly size: number;
  /** Offset of the local header from the start of the archive. */
  readonly localOffset: number;
  /** The raw name bytes, compared against the local header's copy on read. */
  readonly nameBytes: Uint8Array;
}

const unreadable = (message: string): QvfError => new QvfError("E_ZIP_UNREADABLE", message);

// --- CRC-32 -----------------------------------------------------------------

const CRC_TABLE: Uint32Array = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

/** CRC-32 (IEEE 802.3, the ZIP checksum) of `data`, as an unsigned integer. */
export function crc32(data: Uint8Array): number {
  let c = 0xffffffff;
  for (let i = 0; i < data.length; i++) c = CRC_TABLE[(c ^ data[i]!) & 0xff]! ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

// --- paths, names ------------------------------------------------------------

const DRIVE_LETTER = /^[A-Za-z]:/;

/**
 * Why a member path or ZIP entry name is illegal (spec § 2.1), or `null`.
 *
 * Illegal: a backslash, a NUL, a leading `/`, a drive-letter prefix, or a `..`
 * path segment.
 */
export function pathProblem(path: string): string | null {
  if (path.includes("\\")) return "contains a backslash";
  if (path.includes("\0")) return "contains a NUL character";
  if (path.startsWith("/")) return "begins with '/'";
  if (DRIVE_LETTER.test(path)) return "begins with a drive letter";
  if (path.split("/").includes("..")) return "contains a '..' segment";
  return null;
}

const strictUtf8 = new TextDecoder("utf-8", { fatal: true, ignoreBOM: true });
const lenientUtf8 = new TextDecoder("utf-8", { ignoreBOM: true });

function decodeName(raw: Uint8Array, flags: number): string {
  if (flags & 0x800) return lenientUtf8.decode(raw);
  try {
    return strictUtf8.decode(raw);
  } catch {
    let s = "";
    for (let i = 0; i < raw.length; i++) s += String.fromCharCode(raw[i]!);
    return s;
  }
}

function bytesEqual(a: Uint8Array, b: Uint8Array): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

// --- central directory ------------------------------------------------------

function u64(view: DataView, offset: number): number {
  const low = view.getUint32(offset, true);
  const high = view.getUint32(offset + 4, true);
  if (high > 0x1fffff) throw unreadable("a ZIP64 field exceeds 2^53 bytes");
  return high * 0x100000000 + low;
}

function findEndOfCentralDirectory(view: DataView, length: number): number {
  if (length < EOCD_SIZE) throw unreadable("too small to be a ZIP archive");
  const lowest = Math.max(0, length - EOCD_SIZE - MAX_COMMENT);
  for (let p = length - EOCD_SIZE; p >= lowest; p--) {
    if (view.getUint32(p, true) !== SIG_EOCD) continue;
    const commentLength = view.getUint16(p + 20, true);
    if (p + EOCD_SIZE + commentLength <= length) return p;
  }
  throw unreadable("no end-of-central-directory record (not a ZIP archive)");
}

interface Directory {
  count: number;
  size: number;
  offset: number;
}

function readDirectoryLocation(view: DataView, length: number, eocd: number): Directory {
  let count = view.getUint16(eocd + 10, true);
  let size = view.getUint32(eocd + 12, true);
  let offset = view.getUint32(eocd + 16, true);
  const disk = view.getUint16(eocd + 4, true);
  const directoryDisk = view.getUint16(eocd + 6, true);
  const saturated = count === SATURATED16 || size === SATURATED32 || offset === SATURATED32;
  const locator = eocd - 20;

  if (saturated && locator >= 0 && view.getUint32(locator, true) === SIG_LOCATOR64) {
    const record = u64(view, locator + 8);
    if (record + 56 > length || view.getUint32(record, true) !== SIG_EOCD64) {
      throw unreadable("invalid ZIP64 end-of-central-directory record");
    }
    if (view.getUint32(locator + 4, true) !== 0 || view.getUint32(locator + 16, true) !== 1) {
      throw unreadable("multi-volume ZIP archives are not supported");
    }
    if (view.getUint32(record + 16, true) !== 0 || view.getUint32(record + 20, true) !== 0) {
      throw unreadable("multi-volume ZIP archives are not supported");
    }
    if (count === SATURATED16) count = u64(view, record + 32);
    if (size === SATURATED32) size = u64(view, record + 40);
    if (offset === SATURATED32) offset = u64(view, record + 48);
  } else if (disk !== 0 || directoryDisk !== 0) {
    throw unreadable("multi-volume ZIP archives are not supported");
  }
  return { count, size, offset };
}

/** Resolve the 0x0001 (ZIP64) extra field of a central-directory record. */
function zip64Fields(
  view: DataView,
  extraStart: number,
  extraLength: number,
  wanted: { size: boolean; compressed: boolean; offset: boolean },
): { size?: number; compressed?: number; offset?: number } {
  const end = extraStart + extraLength;
  let q = extraStart;
  while (q + 4 <= end) {
    const id = view.getUint16(q, true);
    const length = view.getUint16(q + 2, true);
    const body = q + 4;
    if (body + length > end) break;
    if (id === 0x0001) {
      const out: { size?: number; compressed?: number; offset?: number } = {};
      let o = body;
      const take = (): number => {
        if (o + 8 > body + length) throw unreadable("truncated ZIP64 extra field");
        const v = u64(view, o);
        o += 8;
        return v;
      };
      if (wanted.size) out.size = take();
      if (wanted.compressed) out.compressed = take();
      if (wanted.offset) out.offset = take();
      return out;
    }
    q = body + length;
  }
  throw unreadable("a size or offset is saturated but there is no ZIP64 extra field");
}

function readCentralDirectory(bytes: Uint8Array, view: DataView, dir: Directory): ZipEntry[] {
  const length = bytes.length;
  if (dir.offset > length || dir.size > length - dir.offset) {
    throw unreadable("central directory out of range");
  }
  if (dir.count > Math.floor(dir.size / CENTRAL_HEADER_SIZE)) {
    throw unreadable("central directory too small for its entry count");
  }
  const entries: ZipEntry[] = [];
  let p = dir.offset;
  for (let i = 0; i < dir.count; i++) {
    if (p + CENTRAL_HEADER_SIZE > length || view.getUint32(p, true) !== SIG_CENTRAL) {
      throw unreadable("bad central directory entry");
    }
    const flags = view.getUint16(p + 8, true);
    const method = view.getUint16(p + 10, true);
    const crc = view.getUint32(p + 16, true);
    let compressedSize = view.getUint32(p + 20, true);
    let size = view.getUint32(p + 24, true);
    const nameLength = view.getUint16(p + 28, true);
    const extraLength = view.getUint16(p + 30, true);
    const commentLength = view.getUint16(p + 32, true);
    let localOffset = view.getUint32(p + 42, true);
    const nameStart = p + CENTRAL_HEADER_SIZE;
    const next = nameStart + nameLength + extraLength + commentLength;
    if (next > length) throw unreadable("truncated central directory entry");

    const wanted = {
      size: size === SATURATED32,
      compressed: compressedSize === SATURATED32,
      offset: localOffset === SATURATED32,
    };
    if (wanted.size || wanted.compressed || wanted.offset) {
      const z = zip64Fields(view, nameStart + nameLength, extraLength, wanted);
      if (z.size !== undefined) size = z.size;
      if (z.compressed !== undefined) compressedSize = z.compressed;
      if (z.offset !== undefined) localOffset = z.offset;
    }

    const nameBytes = bytes.subarray(nameStart, nameStart + nameLength);
    entries.push({
      name: decodeName(nameBytes, flags),
      method,
      flags,
      crc32: crc,
      compressedSize,
      size,
      localOffset,
      nameBytes,
    });
    p = next;
  }
  return entries;
}

// --- inflate ----------------------------------------------------------------

/**
 * Inflate `data` into a buffer of exactly `size` bytes.
 *
 * The compressed input is fed to the inflater in steps, and a step is sized so
 * that its output cannot exceed the declared size by more than about 1 MiB
 * (DEFLATE expands at most 1032:1). A stream that expands beyond `size` is
 * therefore abandoned almost at once instead of being inflated in full: a lying
 * header cannot make the reader produce, or allocate, much more than declared.
 */
function inflateExact(data: Uint8Array, size: number, name: string): Uint8Array {
  // The output buffer is allocated up front, so first refuse a declared size
  // that the stored bytes cannot reach: a tiny entry declaring gigabytes must
  // not make the reader allocate them (spec § 2.1).
  if (size > data.length * DEFLATE_MAX_EXPANSION + DEFLATE_SLACK) {
    throw unreadable(
      `ZIP entry ${JSON.stringify(name)} declares ${size} bytes, more than its ${data.length} compressed bytes can inflate to`,
    );
  }
  let out: Uint8Array;
  try {
    out = new Uint8Array(size);
  } catch {
    throw unreadable(`ZIP entry ${JSON.stringify(name)}: cannot allocate ${size} bytes`);
  }
  let written = 0;
  let overflow = false;
  const inflater = new Inflate((chunk) => {
    if (overflow) return;
    if (written + chunk.length > size) {
      overflow = true;
      return;
    }
    out.set(chunk, written);
    written += chunk.length;
  });
  const step = Math.min(INFLATE_STEP_MAX, Math.ceil(size / DEFLATE_MAX_EXPANSION) + 1024);
  try {
    let pos = 0;
    do {
      const end = Math.min(pos + step, data.length);
      inflater.push(data.subarray(pos, end), end === data.length);
      pos = end;
    } while (pos < data.length && !overflow);
  } catch (err) {
    const reason = err instanceof Error ? err.message : String(err);
    throw unreadable(`ZIP entry ${JSON.stringify(name)} does not inflate (${reason})`);
  }
  if (overflow) {
    throw unreadable(
      `ZIP entry ${JSON.stringify(name)} inflates to more than its declared ${size} bytes`,
    );
  }
  if (written !== size) {
    throw unreadable(
      `ZIP entry ${JSON.stringify(name)} inflates to ${written} bytes, not its declared ${size}`,
    );
  }
  return out;
}

// --- the archive ------------------------------------------------------------

/** Whether a compressed entry expands past the ratio cap (spec § 6.1). */
export function ratioExceeded(entry: ZipEntry, limits: Limits): boolean {
  return (
    limits.max_ratio > 0 &&
    entry.method !== METHOD_STORE &&
    entry.size >= limits.ratio_exempt_bytes &&
    (entry.compressedSize === 0 || entry.size > limits.max_ratio * entry.compressedSize)
  );
}

/** A ZIP archive held in memory, opened through its central directory. */
export class ZipArchive {
  /** Every central-directory record, in archive order, duplicates included. */
  readonly entries: readonly ZipEntry[];
  /** Names that occur more than once, in order of their second occurrence. */
  readonly duplicates: readonly string[];
  /** Sum of the declared uncompressed sizes of all entries. */
  readonly totalSize: number;

  private readonly bytes: Uint8Array;
  private readonly index: Map<string, ZipEntry>;

  private constructor(bytes: Uint8Array, entries: ZipEntry[]) {
    this.bytes = bytes;
    this.entries = entries;
    this.index = new Map();
    const duplicates: string[] = [];
    let total = 0;
    for (const entry of entries) {
      if (this.index.has(entry.name) && !duplicates.includes(entry.name)) {
        duplicates.push(entry.name);
      }
      this.index.set(entry.name, entry); // the last record wins, like Python's zipfile
      total += entry.size;
    }
    this.duplicates = duplicates;
    this.totalSize = total;
  }

  /**
   * Open an archive. Throws `QvfError` with code `E_ZIP_UNREADABLE` when the
   * bytes are not a ZIP archive or its central directory is damaged; entries
   * are not decompressed here.
   */
  static open(bytes: Uint8Array): ZipArchive {
    try {
      // Node's binary-data subclass of Uint8Array overrides slice() to share memory;
      // always work on a plain Uint8Array view of the same bytes.
      const data =
        Object.getPrototypeOf(bytes) === Uint8Array.prototype
          ? bytes
          : new Uint8Array(bytes.buffer, bytes.byteOffset, bytes.byteLength);
      const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
      const eocd = findEndOfCentralDirectory(view, data.length);
      const directory = readDirectoryLocation(view, data.length, eocd);
      return new ZipArchive(data, readCentralDirectory(data, view, directory));
    } catch (err) {
      if (err instanceof QvfError) throw err;
      const reason = err instanceof Error ? err.message : String(err);
      throw unreadable(`damaged ZIP structure (${reason})`);
    }
  }

  /** Whether an entry of that name exists. */
  has(name: string): boolean {
    return this.index.has(name);
  }

  /** The record a lookup by name resolves to (the last one when duplicated). */
  find(name: string): ZipEntry | undefined {
    return this.index.get(name);
  }

  /** Unique entry names in order of first occurrence. */
  names(): string[] {
    return [...new Set(this.entries.map((e) => e.name))];
  }

  /**
   * The uncompressed bytes of an entry, as a fresh buffer.
   *
   * Throws `QvfError`: `E_MEMBER_MISSING` for an unknown name,
   * `E_MEMBER_TOO_LARGE` when the declared size exceeds `maxBytes` (nothing is
   * decompressed), and `E_ZIP_UNREADABLE` when the data cannot be read
   * (truncated, encrypted, unsupported method, does not inflate to exactly the
   * declared size, CRC-32 mismatch).
   */
  read(target: ZipEntry | string, options: { maxBytes?: number } = {}): Uint8Array {
    const entry = typeof target === "string" ? this.index.get(target) : target;
    if (!entry) {
      throw new QvfError("E_MEMBER_MISSING", `no ZIP entry named ${JSON.stringify(target)}`);
    }
    const label = JSON.stringify(entry.name);
    if (options.maxBytes !== undefined && entry.size > options.maxBytes) {
      throw new QvfError(
        "E_MEMBER_TOO_LARGE",
        `ZIP entry ${label} declares ${entry.size} bytes, over the ${options.maxBytes}-byte cap`,
      );
    }
    if (entry.flags & 0x1 || entry.flags & 0x40) {
      throw unreadable(`ZIP entry ${label} is encrypted`);
    }

    const bytes = this.bytes;
    const offset = entry.localOffset;
    if (offset + LOCAL_HEADER_SIZE > bytes.length) {
      throw unreadable(`ZIP entry ${label}: truncated local header`);
    }
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    if (view.getUint32(offset, true) !== SIG_LOCAL) {
      throw unreadable(`ZIP entry ${label}: bad local header signature`);
    }
    const nameLength = view.getUint16(offset + 26, true);
    const extraLength = view.getUint16(offset + 28, true);
    const nameStart = offset + LOCAL_HEADER_SIZE;
    const dataStart = nameStart + nameLength + extraLength;
    if (dataStart > bytes.length || entry.compressedSize > bytes.length - dataStart) {
      throw unreadable(`ZIP entry ${label}: truncated data`);
    }
    if (!bytesEqual(bytes.subarray(nameStart, nameStart + nameLength), entry.nameBytes)) {
      throw unreadable(`ZIP entry ${label}: local header name differs from the directory`);
    }
    const data = bytes.subarray(dataStart, dataStart + entry.compressedSize);

    let out: Uint8Array;
    if (entry.method === METHOD_STORE) {
      if (entry.compressedSize !== entry.size) {
        throw unreadable(`ZIP entry ${label}: stored size differs from uncompressed size`);
      }
      out = data.slice();
    } else if (entry.method === METHOD_DEFLATE) {
      out = inflateExact(data, entry.size, entry.name);
    } else {
      throw unreadable(`ZIP entry ${label}: unsupported compression method ${entry.method}`);
    }
    if (crc32(out) !== entry.crc32) {
      throw unreadable(`ZIP entry ${label}: CRC-32 mismatch`);
    }
    return out;
  }
}
