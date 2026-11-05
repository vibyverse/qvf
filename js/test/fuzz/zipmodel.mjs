// An editable model of a ZIP archive for structure-aware fuzzing: the port of
// python/tests/fuzz/zipmodel.py. parse() turns archive bytes into entries
// (header fields plus the stored bytes); build() writes them back. Each header
// field can be overridden for the local or the central copy alone, an entry
// can carry a ZIP64 extra field with true or lying values, and a data
// descriptor.

import { deflateRawSync, inflateRawSync } from "node:zlib";

import { ZipArchive, crc32 } from "../../dist/index.js";

export const SAT16 = 0xffff;
export const SAT32 = 0xffffffff;
const SIG_LOCAL = 0x04034b50;
const SIG_CENTRAL = 0x02014b50;
const SIG_EOCD = 0x06054b50;
const SIG_EOCD64 = 0x06064b50;
const SIG_LOCATOR64 = 0x07064b50;
const SIG_DESCRIPTOR = 0x08074b50;

const big = (v) => (typeof v === "bigint" ? v : BigInt(Math.trunc(v)));

/** A growable little-endian byte writer. */
class Out {
  constructor() {
    this.parts = [];
    this.length = 0;
  }
  bytes(b) {
    this.parts.push(b);
    this.length += b.length;
  }
  u16(v) {
    const b = Buffer.alloc(2);
    b.writeUInt16LE(Number(BigInt.asUintN(16, big(v))));
    this.bytes(b);
  }
  u32(v) {
    const b = Buffer.alloc(4);
    b.writeUInt32LE(Number(BigInt.asUintN(32, big(v))));
    this.bytes(b);
  }
  u64(v) {
    const b = Buffer.alloc(8);
    b.writeBigUInt64LE(BigInt.asUintN(64, big(v)));
    this.bytes(b);
  }
  finish() {
    return Buffer.concat(this.parts, this.length);
  }
}

export function newEntry(fields) {
  return {
    name: Buffer.alloc(0), data: Buffer.alloc(0), method: 0, flags: 0, crc: 0, usize: 0,
    local: {}, central: {}, localExtra: Buffer.alloc(0), centralExtra: Buffer.alloc(0),
    comment: Buffer.alloc(0), zip64: [], zip64Values: {}, zip64Truncate: 0,
    descriptor: null, descriptorValues: {}, zeroLocalSizes: false, ...fields,
  };
}

export function cloneEntry(e) {
  return { ...e, local: { ...e.local }, central: { ...e.central }, zip64: [...e.zip64],
    zip64Values: { ...e.zip64Values }, descriptorValues: { ...e.descriptorValues } };
}

/** A fresh, valid entry holding `payload`. */
export function store(name, payload, deflate = true) {
  const data = deflate ? deflateRawSync(payload) : Buffer.from(payload);
  return newEntry({ name: Buffer.from(name), data, method: deflate ? 8 : 0,
    crc: crc32(payload), usize: payload.length });
}

/** The uncompressed payload of a STORE or DEFLATE entry, or null. */
export function payloadOf(e) {
  try {
    if (e.method === 0) return Buffer.from(e.data);
    if (e.method === 8) return inflateRawSync(e.data);
  } catch {
    return null;
  }
  return null;
}

export function parse(raw) {
  const zip = ZipArchive.open(new Uint8Array(raw));
  const view = Buffer.from(raw);
  const entries = [];
  for (const z of zip.entries) {
    const off = z.localOffset;
    const nlen = view.readUInt16LE(off + 26);
    const xlen = view.readUInt16LE(off + 28);
    const start = off + 30 + nlen + xlen;
    entries.push(newEntry({ name: Buffer.from(z.nameBytes), data: view.subarray(start, start + z.compressedSize),
      method: z.method, flags: z.flags & ~0x8, crc: z.crc32, usize: z.size }));
  }
  return { entries, comment: Buffer.alloc(0), prefix: Buffer.alloc(0), suffix: Buffer.alloc(0),
    zip64End: false, end: {}, end64: {} };
}

function zip64Extra(e, offset, local) {
  if (e.zip64.length === 0) return Buffer.alloc(0);
  const v = e.zip64Values;
  const body = new Out();
  const keys = local ? ["usize", "csize"] : ["usize", "csize", "offset"].filter((k) => e.zip64.includes(k));
  const truth = { usize: e.usize, csize: e.data.length, offset };
  for (const k of keys) body.u64(v[k] ?? truth[k]);
  let b = body.finish();
  if (e.zip64Truncate) b = b.subarray(0, Math.max(0, b.length - e.zip64Truncate));
  const out = new Out();
  out.u16(1);
  out.u16(!local && v.length !== undefined ? v.length : b.length);
  out.bytes(b);
  return out.finish();
}

function local(e) {
  const c = e.local;
  let usize = c.usize ?? e.usize;
  let csize = c.csize ?? e.data.length;
  let crc = c.crc ?? e.crc;
  if (e.zeroLocalSizes) usize = csize = crc = 0;
  let extra = e.localExtra;
  if (e.zip64.length) {
    extra = Buffer.concat([zip64Extra(e, 0, true), extra]);
    usize = csize = SAT32;
  }
  const name = c.name ?? e.name;
  const o = new Out();
  o.u32(SIG_LOCAL); o.u16(c.version ?? 20); o.u16(c.flags ?? e.flags); o.u16(c.method ?? e.method);
  o.u16(0); o.u16(0x21); o.u32(crc); o.u32(csize); o.u32(usize); o.u16(name.length);
  o.u16(c.extraLength ?? extra.length); o.bytes(name); o.bytes(extra);
  return o.finish();
}

function central(e, at) {
  const c = e.central;
  let offset = c.offset ?? at;
  let usize = c.usize ?? e.usize;
  let csize = c.csize ?? e.data.length;
  const extra = Buffer.concat([zip64Extra(e, offset, false), e.centralExtra]);
  if (e.zip64.includes("usize")) usize = SAT32;
  if (e.zip64.includes("csize")) csize = SAT32;
  if (e.zip64.includes("offset")) offset = SAT32;
  const name = c.name ?? e.name;
  const o = new Out();
  o.u32(SIG_CENTRAL); o.u16(e.zip64.length ? 45 : 20); o.u16(c.version ?? (e.zip64.length ? 45 : 20));
  o.u16(c.flags ?? e.flags); o.u16(c.method ?? e.method); o.u16(0); o.u16(0x21);
  o.u32(c.crc ?? e.crc); o.u32(csize); o.u32(usize); o.u16(name.length);
  o.u16(c.extraLength ?? extra.length); o.u16(e.comment.length); o.u16(0); o.u16(0); o.u32(0);
  o.u32(offset); o.bytes(name); o.bytes(extra); o.bytes(e.comment);
  return o.finish();
}

function descriptor(e) {
  if (!e.descriptor) return Buffer.alloc(0);
  const v = e.descriptorValues;
  const o = new Out();
  if (e.descriptor !== "nosig32") o.u32(SIG_DESCRIPTOR);
  o.u32(v.crc ?? e.crc);
  if (e.descriptor === "sig64") {
    o.u64(v.csize ?? e.data.length); o.u64(v.usize ?? e.usize);
  } else {
    o.u32(v.csize ?? e.data.length); o.u32(v.usize ?? e.usize);
  }
  return o.finish();
}

export function build(a) {
  const out = new Out();
  out.bytes(a.prefix);
  const cd = new Out();
  for (const e of a.entries) {
    const at = out.length;
    out.bytes(local(e));
    out.bytes(e.data);
    out.bytes(descriptor(e));
    cd.bytes(central(e, at));
  }
  const cdOffset = out.length;
  const cdBytes = cd.finish();
  out.bytes(cdBytes);
  const n = a.entries.length;
  let count = a.end.count ?? n;
  let cdSize = a.end.cdSize ?? cdBytes.length;
  let cdOff = a.end.cdOffset ?? cdOffset;
  if (a.zip64End) {
    const record = out.length;
    const x = a.end64;
    out.u32(SIG_EOCD64); out.u64(x.recordSize ?? 44); out.u16(45); out.u16(45);
    out.u32(x.disk ?? 0); out.u32(x.cdDisk ?? 0); out.u64(x.count ?? n); out.u64(x.count ?? n);
    out.u64(x.cdSize ?? cdBytes.length); out.u64(x.cdOffset ?? cdOffset);
    out.u32(SIG_LOCATOR64); out.u32(x.locatorDisk ?? 0); out.u64(x.locatorOffset ?? record);
    out.u32(x.disks ?? 1);
    count = SAT16; cdSize = SAT32; cdOff = SAT32;
  }
  out.u32(SIG_EOCD); out.u16(a.end.disk ?? 0); out.u16(a.end.cdDisk ?? 0); out.u16(count);
  out.u16(count); out.u32(cdSize); out.u32(cdOff); out.u16(a.end.commentLength ?? a.comment.length);
  out.bytes(a.comment);
  out.bytes(a.suffix);
  return out.finish();
}

export function find(a, name) {
  const key = Buffer.from(name);
  let hit = null;
  for (const e of a.entries) if (Buffer.compare(e.name, key) === 0) hit = e;
  return hit;
}
