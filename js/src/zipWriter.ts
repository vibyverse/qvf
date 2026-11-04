/**
 * Container layer, write side: a ZIP writer for whole entries (spec § 2.1).
 *
 * Each entry is CRC-checked and, for DEFLATE, compressed with `fflate` before
 * its local header is written, so no data descriptors are needed. ZIP64
 * records are written when a size, an offset or the entry count needs them,
 * or for every entry with `forceZip64`. Output goes to a callback as it is
 * produced, or is collected in memory.
 *
 * Nothing here depends on Node: the module runs in browsers unchanged.
 */

import { deflateSync } from "fflate";

import { QvfError } from "./errors.js";
import { crc32, METHOD_DEFLATE, METHOD_STORE } from "./zip.js";

/** Receives the archive bytes in order. A returned promise is awaited before the next chunk. */
export type ZipSink = (chunk: Uint8Array) => void | Promise<void>;

export interface ZipWriterOptions {
  /** Where the bytes go; without one they are kept and returned by `bytes()` / `blob()`. */
  sink?: ZipSink;
  /** Write ZIP64 fields for every entry and the ZIP64 end records, whatever the sizes. */
  forceZip64?: boolean;
  /** Modification time of every entry (its UTC fields). Default 1980-01-01T00:00:00Z, the DOS epoch. */
  modified?: Date;
}

interface CentralRecord {
  nameBytes: Uint8Array;
  flags: number;
  method: number;
  crc: number;
  compressedSize: number;
  size: number;
  offset: number;
}

const U32_MAX = 0xffffffff;
const U16_MAX = 0xffff;
/** Unix host (3), APPNOTE version 2.0 (20) or 4.5 (45, ZIP64). */
const MADE_BY_UNIX = 3 << 8;
/** Regular file, rw-r--r--, in the high half of the external attributes. */
const EXTERNAL_ATTRIBUTES = (0o100644 << 16) >>> 0;

const encoder = new TextEncoder();

function dosDateTime(date: Date): [number, number] {
  const year = date.getUTCFullYear();
  if (Number.isNaN(year) || year < 1980) return [0, (1 << 5) | 1];
  if (year > 2107) throw new RangeError(`ZIP time stamps end in 2107, got ${date.toISOString()}`);
  const time = (date.getUTCHours() << 11) | (date.getUTCMinutes() << 5) | (date.getUTCSeconds() >> 1);
  const day = ((year - 1980) << 9) | ((date.getUTCMonth() + 1) << 5) | date.getUTCDate();
  return [time, day];
}

class Bytes {
  readonly buf: Uint8Array;
  private readonly view: DataView;
  private o = 0;
  constructor(n: number) {
    this.buf = new Uint8Array(n);
    this.view = new DataView(this.buf.buffer);
  }
  u16(v: number): this {
    this.view.setUint16(this.o, v, true);
    this.o += 2;
    return this;
  }
  u32(v: number): this {
    this.view.setUint32(this.o, v >>> 0, true);
    this.o += 4;
    return this;
  }
  u64(v: number): this {
    this.u32(v % 0x100000000);
    return this.u32(Math.floor(v / 0x100000000));
  }
  raw(b: Uint8Array): this {
    this.buf.set(b, this.o);
    this.o += b.length;
    return this;
  }
}

/** A ZIP archive written entry by entry. */
export class ZipWriter {
  private readonly sink: ZipSink | undefined;
  private readonly forceZip64: boolean;
  private readonly time: number;
  private readonly date: number;
  private readonly parts: Uint8Array[] = [];
  private readonly records: CentralRecord[] = [];
  private readonly names = new Set<string>();
  private offset = 0;
  private finished = false;

  constructor(options: ZipWriterOptions = {}) {
    this.sink = options.sink;
    this.forceZip64 = options.forceZip64 ?? false;
    [this.time, this.date] = dosDateTime(options.modified ?? new Date(Date.UTC(1980, 0, 1)));
  }

  /** Bytes written so far. */
  get length(): number {
    return this.offset;
  }

  /** Whether an entry of that name was written. */
  has(name: string): boolean {
    return this.names.has(name);
  }

  private async emit(chunk: Uint8Array): Promise<void> {
    this.offset += chunk.length;
    if (this.sink) await this.sink(chunk);
    else this.parts.push(chunk);
  }

  /**
   * Write one entry. `method` 8 deflates `data` (at `level`, default 6); 0
   * stores it. Names must be unique (`E_ZIP_DUPLICATE_ENTRY`). `data` is
   * not retained after the call returns unless it is stored and no sink is
   * set, in which case a copy is kept.
   */
  async add(name: string, data: Uint8Array, method: number = METHOD_DEFLATE, level = 6): Promise<void> {
    if (this.finished) throw new Error("the ZIP archive is already finished");
    if (this.names.has(name)) throw new QvfError("E_ZIP_DUPLICATE_ENTRY", `duplicate ZIP entry name ${JSON.stringify(name)}`);
    if (method !== METHOD_STORE && method !== METHOD_DEFLATE) throw new RangeError(`unsupported ZIP method ${method}`);
    this.names.add(name);
    const nameBytes = encoder.encode(name);
    if (nameBytes.length > U16_MAX) throw new RangeError(`ZIP entry name of ${nameBytes.length} bytes is too long`);
    const ascii = nameBytes.every((b) => b < 0x80);
    const flags = ascii ? 0 : 0x800;
    const crc = crc32(data);
    const stored =
      method === METHOD_DEFLATE
        ? deflateSync(data, { level: level as 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 })
        : this.sink
          ? data
          : data.slice();
    const offset = this.offset;
    const zip64 = this.forceZip64 || data.length >= U32_MAX || stored.length >= U32_MAX;
    const header = new Bytes(30 + nameBytes.length + (zip64 ? 20 : 0))
      .u32(0x04034b50)
      .u16(zip64 ? 45 : 20)
      .u16(flags)
      .u16(method)
      .u16(this.time)
      .u16(this.date)
      .u32(crc)
      .u32(zip64 ? U32_MAX : stored.length)
      .u32(zip64 ? U32_MAX : data.length)
      .u16(nameBytes.length)
      .u16(zip64 ? 20 : 0)
      .raw(nameBytes);
    if (zip64) header.u16(0x0001).u16(16).u64(data.length).u64(stored.length);
    await this.emit(header.buf);
    await this.emit(stored);
    this.records.push({ nameBytes, flags, method, crc, compressedSize: stored.length, size: data.length, offset });
  }

  /** Write the central directory and the end records. No entry can be added afterwards. */
  async finish(): Promise<void> {
    if (this.finished) return;
    this.finished = true;
    const start = this.offset;
    for (const r of this.records) {
      const sizes64 = this.forceZip64 || r.size >= U32_MAX || r.compressedSize >= U32_MAX;
      const offset64 = this.forceZip64 || r.offset >= U32_MAX;
      const extraLength = (sizes64 || offset64 ? 4 : 0) + (sizes64 ? 16 : 0) + (offset64 ? 8 : 0);
      const version = sizes64 || offset64 ? 45 : 20;
      const rec = new Bytes(46 + r.nameBytes.length + extraLength)
        .u32(0x02014b50)
        .u16(MADE_BY_UNIX | version)
        .u16(version)
        .u16(r.flags)
        .u16(r.method)
        .u16(this.time)
        .u16(this.date)
        .u32(r.crc)
        .u32(sizes64 ? U32_MAX : r.compressedSize)
        .u32(sizes64 ? U32_MAX : r.size)
        .u16(r.nameBytes.length)
        .u16(extraLength)
        .u16(0) // comment
        .u16(0) // disk
        .u16(0) // internal attributes
        .u32(EXTERNAL_ATTRIBUTES)
        .u32(offset64 ? U32_MAX : r.offset)
        .raw(r.nameBytes);
      if (extraLength > 0) {
        rec.u16(0x0001).u16(extraLength - 4);
        if (sizes64) rec.u64(r.size).u64(r.compressedSize);
        if (offset64) rec.u64(r.offset);
      }
      await this.emit(rec.buf);
    }
    const size = this.offset - start;
    const count = this.records.length;
    const zip64 = this.forceZip64 || count >= U16_MAX || size >= U32_MAX || start >= U32_MAX;
    if (zip64) {
      const eocd64 = this.offset;
      await this.emit(
        new Bytes(56)
          .u32(0x06064b50)
          .u64(44)
          .u16(MADE_BY_UNIX | 45)
          .u16(45)
          .u32(0)
          .u32(0)
          .u64(count)
          .u64(count)
          .u64(size)
          .u64(start).buf,
      );
      await this.emit(new Bytes(20).u32(0x07064b50).u32(0).u64(eocd64).u32(1).buf);
    }
    const sat16 = (v: number): number => (zip64 ? U16_MAX : v);
    const sat32 = (v: number): number => (zip64 ? U32_MAX : v);
    await this.emit(
      new Bytes(22)
        .u32(0x06054b50)
        .u16(0)
        .u16(0)
        .u16(sat16(count))
        .u16(sat16(count))
        .u32(sat32(size))
        .u32(sat32(start))
        .u16(0).buf,
    );
  }

  /** The archive in one buffer (in-memory mode, after `finish`). */
  bytes(): Uint8Array {
    const parts = this.collected();
    const out = new Uint8Array(this.offset);
    let o = 0;
    for (const p of parts) {
      out.set(p, o);
      o += p.length;
    }
    return out;
  }

  /** The archive as a `Blob` (in-memory mode, after `finish`); it can exceed the largest single buffer. */
  blob(type = "application/zip"): Blob {
    return new Blob(this.collected() as BlobPart[], { type });
  }

  private collected(): Uint8Array[] {
    if (!this.finished) throw new Error("finish the archive first");
    if (this.sink) throw new Error("the archive went to the sink; it is not kept in memory");
    return this.parts;
  }
}
