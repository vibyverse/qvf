// Structure-aware mutations of QVF archives at the ZIP and the JSON level:
// the port of python/tests/fuzz/mutators.py, driven by a seeded PRNG so that
// every input is reproducible from (seed archive, plan, byte seed).

import { createHash } from "node:crypto";
import { readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { deflateRawSync } from "node:zlib";

import { crc32 as crc32Of } from "../../dist/index.js";

import { SAT32, build, cloneEntry, find, parse, payloadOf, store } from "./zipmodel.mjs";

const ROOT = fileURLToPath(new URL("../../../", import.meta.url));

/** mulberry32: a small, fast, seedable PRNG. */
export function rngFrom(seed) {
  let a = seed >>> 0;
  const next = () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  return {
    random: next,
    int: (n) => Math.floor(next() * n),
    randint: (lo, hi) => lo + Math.floor(next() * (hi - lo + 1)),
    choice: (xs) => xs[Math.floor(next() * xs.length)],
    bytes: (n) => Buffer.from(Array.from({ length: n }, () => Math.floor(next() * 256))),
  };
}

const INTERESTING_INTS = [0, 1, -1, 2, 3, 7, 255, 256, 0x7fff, 0xffff, 0x10000, 2 ** 31 - 1, 2 ** 31,
  2 ** 32 - 1, 2 ** 32, 2 ** 53, 2 ** 53 + 2, 2 ** 63, 2 ** 64, 1e30, -(2 ** 63)];
const SIZES = [0, 1, 2, 30, 0xffff, 0x10000, 2 ** 24, 2 ** 30, 2 ** 31 - 1, 2 ** 31, 0xfffffffe, SAT32,
  2 ** 33, 8 * 1024 ** 3, 8 * 1024 ** 3 + 1, 2 ** 53 - 1, 2 ** 53, 2n ** 62n, 2n ** 64n - 1n];
const EVIL_NAMES = ["../escape", "a/../../b", "/abs/path", "C:/drive", "c:relative", "a\\b",
  "back\\..\\slash", "nul\u0000byte", "dir/", "", ".", "..", "./manifest.json", "manifest.json",
  "MANIFEST.JSON", "manifest.json/", "a//b", "\u00e9t\u00e9", "x".repeat(300)];
const EVIL_STRINGS = ["", " ", "x", "\u0000", "\u{1F600}", "float128", "complex64", "<f8", "bool",
  "float64 ", "Float64", "json", "binary", "opaque", "../x", "/x", "0".repeat(64), "A".repeat(64),
  "g".repeat(64), "x".repeat(100000)];

export function seeds() {
  const out = [];
  for (const d of ["corpus", "negative", "legacy"]) {
    const dir = join(ROOT, "conformance", d);
    for (const f of readdirSync(dir).filter((x) => x.endsWith(".qvf")).sort()) {
      out.push([`${d}/${f}`, readFileSync(join(dir, f))]);
    }
  }
  return out;
}

const entry = (r, a) => (a.entries.length ? r.choice(a.entries) : null);
const which = (r) => r.choice(["local", "central", "both"]);
function set(e, where, key, value) {
  if (where !== "central") e.local[key] = value;
  if (where !== "local") e.central[key] = value;
}
function sizeValue(r, truth) {
  const x = r.random();
  if (x < 0.4) return Math.max(0, truth + r.choice([-1, 1, -2, 2, -1024, 1024, r.randint(-truth, 4096)]));
  if (x < 0.8) return Number(BigInt.asUintN(32, BigInt(r.choice(SIZES))));
  return r.int(2 ** 32);
}

// --- ZIP-level operators -----------------------------------------------------

const zipOps = {
  method(r, a) {
    const e = entry(r, a);
    if (e) set(e, which(r), "method", r.choice([0, 8, 1, 6, 9, 12, 14, 93, 95, 98, 99, r.int(0x10000)]));
  },
  flags(r, a) {
    const e = entry(r, a);
    if (e) set(e, which(r), "flags", e.flags ^ r.choice([0x1, 0x8, 0x40, 0x800, 0x2000, 0x20, 1 << r.int(16)]));
  },
  crc(r, a) {
    const e = entry(r, a);
    if (e) set(e, which(r), "crc", r.choice([0, SAT32, (e.crc ^ (1 << r.int(32))) >>> 0, r.int(2 ** 32)]));
  },
  csize(r, a) {
    const e = entry(r, a);
    if (e) set(e, which(r), "csize", sizeValue(r, e.data.length));
  },
  usize(r, a) {
    const e = entry(r, a);
    if (e) set(e, which(r), "usize", sizeValue(r, e.usize));
  },
  zip64(r, a) {
    const e = entry(r, a);
    if (!e) return;
    const keys = ["usize", "csize", "offset"].filter(() => r.random() < 0.6);
    e.zip64 = keys.length ? keys : ["usize"];
    if (r.random() < 0.5) e.zip64Values[r.choice(e.zip64)] = r.choice([...SIZES, 2n ** 63n]);
    if (r.random() < 0.15) e.zip64Truncate = r.randint(1, 8 * e.zip64.length);
    if (r.random() < 0.1) e.zip64Values.length = r.choice([0, 4, 7, 9, 0xffff]);
    if (r.random() < 0.3) a.zip64End = true;
  },
  zip64End(r, a) {
    a.zip64End = true;
    const x = r.random();
    if (x < 0.2) a.end64.count = r.choice(SIZES);
    else if (x < 0.4) a.end64.cdOffset = r.choice(SIZES);
    else if (x < 0.5) a.end64.cdSize = r.choice(SIZES);
    else if (x < 0.6) a.end64.locatorOffset = r.choice(SIZES);
    else if (x < 0.7) a.end64.disks = r.choice([0, 2, SAT32]);
    else if (x < 0.8) a.end64.disk = r.choice([1, SAT32]);
    else if (x < 0.9) a.end64.recordSize = r.choice([0, 12, 2n ** 63n]);
  },
  descriptor(r, a) {
    const e = entry(r, a);
    if (!e) return;
    e.flags |= 0x8;
    e.descriptor = r.choice(["sig32", "nosig32", "sig64"]);
    e.zeroLocalSizes = r.random() < 0.7;
    if (r.random() < 0.4) e.descriptorValues[r.choice(["crc", "csize", "usize"])] = r.choice(SIZES);
    if (r.random() < 0.1) e.flags &= ~0x8;
  },
  duplicate(r, a) {
    const e = entry(r, a);
    if (!e) return;
    const d = cloneEntry(e);
    if (r.random() < 0.3) d.data = d.data.subarray(0, d.data.length >> 1);
    a.entries.splice(r.int(a.entries.length + 1), 0, d);
  },
  duplicateManifest(r, a) {
    const m = find(a, "manifest.json");
    if (m) a.entries.splice(r.int(a.entries.length + 1), 0, cloneEntry(m));
  },
  rename(r, a) {
    const e = entry(r, a);
    if (!e) return;
    const x = r.random();
    let name;
    if (x < 0.6) name = Buffer.from(r.choice(EVIL_NAMES));
    else if (x < 0.8 && a.entries.length > 1) name = r.choice(a.entries).name;
    else name = r.bytes(r.randint(0, 12));
    const where = r.choice(["both", "both", "local", "central"]);
    if (where === "both") e.name = name;
    else set(e, where, "name", name);
    if (r.random() < 0.3) e.flags ^= 0x800;
  },
  offset(r, a) {
    const e = entry(r, a);
    if (e) e.central.offset = r.choice([...SIZES.filter((s) => typeof s === "number"), r.int(4096)]);
  },
  corruptData(r, a) {
    const e = entry(r, a);
    if (!e || !e.data.length) return;
    const d = Buffer.from(e.data);
    const x = r.random();
    if (x < 0.5) {
      for (let i = r.randint(1, 4); i > 0; i--) d[r.int(d.length)] ^= 1 << r.int(8);
      e.data = d;
    } else if (x < 0.75) e.data = d.subarray(0, r.int(d.length));
    else e.data = Buffer.concat([d, r.bytes(r.randint(1, 64))]);
  },
  recompress(r, a) {
    const e = entry(r, a);
    const payload = e && payloadOf(e);
    if (!payload) return;
    const fresh = store("x", payload, e.method === 0 || r.random() < 0.3);
    e.data = fresh.data;
    e.method = fresh.method;
  },
  drop(r, a) {
    if (a.entries.length) a.entries.splice(r.int(a.entries.length), 1);
  },
  reorder(r, a) {
    for (let i = a.entries.length - 1; i > 0; i--) {
      const j = r.int(i + 1);
      [a.entries[i], a.entries[j]] = [a.entries[j], a.entries[i]];
    }
  },
  extra(r, a) {
    const e = entry(r, a);
    if (!e) return;
    const parts = [];
    for (let i = r.randint(1, 3); i > 0; i--) {
      const body = r.bytes(r.randint(0, 24));
      const head = Buffer.alloc(4);
      head.writeUInt16LE(r.choice([0x0001, 0x5455, 0x7875, 0x000a, 0x9901, r.int(0x10000)]));
      head.writeUInt16LE(r.random() < 0.7 ? body.length : r.int(0x10000), 2);
      parts.push(head, body);
    }
    if (r.random() < 0.5) e.centralExtra = Buffer.concat(parts);
    else e.localExtra = Buffer.concat(parts);
    if (r.random() < 0.1) set(e, which(r), "extraLength", r.choice([0, 1, 3, 0xffff]));
  },
  end(r, a) {
    const x = r.random();
    const n = a.entries.length;
    if (x < 0.2) a.end.count = r.choice([0, Math.max(0, n - 1), n + 1, 0xffff]);
    else if (x < 0.35) a.end.cdSize = r.choice([0, 1, 45, 46, SAT32, r.int(1 << 20)]);
    else if (x < 0.5) a.end.cdOffset = r.choice([0, 1, SAT32, r.int(1 << 20)]);
    else if (x < 0.6) a.end.commentLength = r.choice([1, 100, 0xffff]);
    else if (x < 0.7) a.end.disk = r.choice([1, 0xffff]);
    else if (x < 0.8) a.comment = Buffer.concat([Buffer.from("PK\x05\x06", "latin1"), r.bytes(18)]);
    else if (x < 0.9) a.prefix = r.bytes(r.randint(1, 64));
    else a.suffix = r.bytes(r.randint(1, 64));
  },
  bomb(r, a) {
    const e = entry(r, a);
    if (!e) return;
    const n = r.choice([1 << 20, 1 << 24, 1 << 26]);
    const z = zeros(n);
    e.data = z.data;
    e.crc = z.crc;
    e.method = 8;
    const x = r.random();
    if (x < 0.4) e.usize = n;
    else if (x < 0.7) e.usize = r.choice([n >> 1, 1024, 0]);
    else e.usize = r.choice([n * 2, 2 ** 31, 0xfffffffe]);
  },
  declareHuge(r, a) {
    const e = entry(r, a);
    if (!e) return;
    const size = r.choice([2 ** 30, 2 ** 31, 0xfffffffe, 2 ** 33, 7 * 1024 ** 3, 8 * 1024 ** 3]);
    if (size > 0xfffffffe || r.random() < 0.3) {
      e.zip64 = ["usize"];
      e.zip64Values.usize = size;
      if (size <= SAT32) e.usize = size;
    } else e.usize = size;
    delete e.local.usize;
  },
};

const zeroCache = new Map();
function zeros(n) {
  if (!zeroCache.has(n)) {
    const z = Buffer.alloc(n);
    zeroCache.set(n, { data: deflateRawSync(z, { level: 9 }), crc: crc32Of(z) });
  }
  return zeroCache.get(n);
}

// --- manifest / JSON-level operators -----------------------------------------

function deep(depth) {
  let v = 0;
  for (let i = 0; i < depth; i++) v = i % 2 ? [v] : { k: v };
  return v;
}

function evil(r) {
  const x = r.random();
  if (x < 0.25) return r.choice(INTERESTING_INTS);
  if (x < 0.4) return r.choice([null, true, false, 0.5, -0, 1e308, -1e308, 5e-324, 1.5, 2 ** 53]);
  if (x < 0.55) return r.choice(EVIL_STRINGS);
  if (x < 0.7) return r.choice([[], {}, [null], [[]], { "": null }, new Array(1000).fill(0)]);
  if (x < 0.8) return deep(r.choice([10, 62, 63, 64, 65, 200]));
  if (x < 0.9) return Array.from({ length: r.randint(0, 6) }, () => r.choice(INTERESTING_INTS));
  return { path: r.choice(EVIL_STRINGS), format: r.choice(["json", "binary"]) };
}

function nodes(doc, path = [], out = []) {
  out.push([path, doc]);
  if (Array.isArray(doc)) doc.forEach((v, i) => nodes(v, [...path, i], out));
  else if (doc && typeof doc === "object") for (const [k, v] of Object.entries(doc)) nodes(v, [...path, k], out);
  return out;
}

function replace(doc, path, value) {
  if (!path.length) return value;
  let parent = doc;
  for (const p of path.slice(0, -1)) parent = parent[p];
  parent[path[path.length - 1]] = value;
  return doc;
}

function manifestOf(a) {
  const e = find(a, "manifest.json");
  if (!e) return [null, null];
  const raw = payloadOf(e);
  try {
    return [e, raw ? JSON.parse(raw.toString("utf8")) : null];
  } catch {
    return [e, null];
  }
}

function rewrite(e, raw) {
  const fresh = store("x", raw, e.method === 8);
  Object.assign(e, { data: fresh.data, method: fresh.method, crc: fresh.crc, usize: fresh.usize,
    local: {}, central: {} });
}

const writeJson = (e, doc) => rewrite(e, Buffer.from(JSON.stringify(doc)));

function members(doc) {
  const out = [];
  if (!doc || typeof doc !== "object" || !Array.isArray(doc.sections)) return out;
  for (const s of doc.sections) {
    if (s && typeof s === "object" && s.members && typeof s.members === "object" && !Array.isArray(s.members)) {
      for (const [role, m] of Object.entries(s.members)) {
        if (m && typeof m === "object" && !Array.isArray(m)) out.push([s, role, m]);
      }
    }
  }
  return out;
}

const jsonOps = {
  jsonNode(r, a) {
    let [e, doc] = manifestOf(a);
    if (doc === null || doc === undefined) return;
    const all = nodes(doc);
    const [path, node] = r.choice(all);
    const x = r.random();
    if (x < 0.6 || !path.length) doc = replace(doc, path, evil(r));
    else if (x < 0.8) {
      let parent = doc;
      for (const p of path.slice(0, -1)) parent = parent[p];
      if (Array.isArray(parent)) parent.splice(path[path.length - 1], 1);
      else delete parent[path[path.length - 1]];
    } else if (Array.isArray(node)) node.splice(r.int(node.length + 1), 0, evil(r));
    else if (node && typeof node === "object") {
      node[r.choice(["x", "critical", "sha256", "size", "dtype", "shape", "kind", "id", "structure_ref",
        "members", "qvf_version"])] = evil(r);
    }
    writeJson(e, doc);
  },
  memberField(r, a) {
    const [e, doc] = manifestOf(a);
    const ms = members(doc);
    if (!ms.length) return;
    const [, , m] = r.choice(ms);
    const field = r.choice(["dtype", "shape", "size", "sha256", "path", "format"]);
    if (field === "dtype") {
      m.dtype = r.choice(["int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64", "float32",
        "float64", ...EVIL_STRINGS.slice(0, 10), null, 8, ["float64"]]);
    } else if (field === "shape") {
      m.shape = r.choice([[], [0], [-1], [1.5], ["3"], [true], 3, "3", null, [2 ** 31, 2 ** 31],
        [2 ** 40, 2 ** 40, 2 ** 40], [2 ** 64], [1e30], [2 ** 40, 0], [2 ** 62, 2 ** 62, 0], new Array(100).fill(1), [0, 2 ** 62],
        Array.from({ length: r.randint(1, 5) }, () => r.randint(0, 100))]);
    } else if (field === "size") {
      m.size = r.choice([...INTERESTING_INTS, null, true, 1.0, "12", -0, 1e20]);
    } else if (field === "sha256") {
      const h = typeof m.sha256 === "string" ? m.sha256 : "0".repeat(64);
      m.sha256 = r.choice([h.toUpperCase(), h.slice(0, -1), h + "0", "", null, 0, "0".repeat(64),
        h.slice(0, -1) + (h.endsWith("0") ? "1" : "0")]);
    } else if (field === "path") {
      const others = ms.map(([, , x]) => x.path).filter((p) => typeof p === "string");
      m.path = r.choice([...EVIL_NAMES, ...others, null, 7, ["a"]]);
    } else m.format = r.choice(["json", "binary", "opaque", null, 1, ""]);
    writeJson(e, doc);
  },
  payload(r, a) {
    const [e, doc] = manifestOf(a);
    const ms = members(doc).filter(([, , m]) => m.format === "json" && typeof m.path === "string");
    if (!ms.length) return;
    const [, , m] = r.choice(ms);
    const target = find(a, m.path);
    const raw = target && payloadOf(target);
    if (!raw) return;
    let payload;
    try {
      payload = JSON.parse(raw.toString("utf8"));
    } catch {
      return;
    }
    const [path] = r.choice(nodes(payload));
    payload = replace(payload, path, evil(r));
    const text = Buffer.from(JSON.stringify(payload) ?? "null");
    rewrite(target, text);
    if (r.random() < 0.8) m.sha256 = createHash("sha256").update(text).digest("hex");
    if ("size" in m && r.random() < 0.8) m.size = text.length;
    writeJson(e, doc);
  },
  binaryPayload(r, a) {
    const [e, doc] = manifestOf(a);
    const ms = members(doc).filter(([, , m]) => m.format === "binary" && typeof m.path === "string");
    if (!ms.length) return;
    const [, , m] = r.choice(ms);
    const target = find(a, m.path);
    if (!target) return;
    let raw = payloadOf(target) ?? Buffer.alloc(0);
    const x = r.random();
    if (x < 0.3) raw = raw.subarray(0, r.int(raw.length + 1));
    else if (x < 0.5) raw = Buffer.concat([raw, r.bytes(r.randint(1, 16))]);
    else if (x < 0.75) {
      const vals = [NaN, Infinity, -Infinity, 5e-324, 1e308, -0];
      const out = Buffer.alloc(raw.length - (raw.length % 8));
      for (let i = 0; i < out.length; i += 8) out.writeDoubleLE(r.choice(vals), i);
      raw = out;
    } else raw = r.bytes(raw.length);
    rewrite(target, raw);
    m.sha256 = createHash("sha256").update(raw).digest("hex");
    if ("size" in m) m.size = raw.length;
    writeJson(e, doc);
  },
  jsonText(r, a) {
    const e = find(a, "manifest.json");
    let raw = e && payloadOf(e);
    if (!raw) return;
    let s = raw.toString("latin1");
    const x = r.random();
    if (x < 0.1) s = s.replace("1", "NaN");
    else if (x < 0.2) s = s.includes(": 1") ? s.replace(": 1", ": Infinity") : "-Infinity";
    else if (x < 0.3) s = s.replace("1", "1" + "0".repeat(r.choice([400, 5000])));
    else if (x < 0.4) {
      const text = raw.toString("utf8");
      const enc = r.choice(["utf16le", "utf16bom"]);
      raw = enc === "utf16le" ? Buffer.from(text, "utf16le") : Buffer.concat([Buffer.from([0xff, 0xfe]), Buffer.from(text, "utf16le")]);
      rewrite(e, raw);
      return;
    } else if (x < 0.45) s = "\xef\xbb\xbf" + s;
    else if (x < 0.55) s += r.choice([" x", "}", "\x00", "\n\n"]);
    else if (x < 0.65) {
      const i = r.int(s.length + 1);
      s = s.slice(0, i) + r.choice(["\xff", "\xc0\xaf", "\xed\xa0\x80", "\x80"]) + s.slice(i);
    } else if (x < 0.75) {
      const d = r.choice([63, 64, 65, 1000, 100000]);
      s = "[".repeat(d) + "]".repeat(d);
    } else if (x < 0.85) {
      const d = r.choice([63, 64, 65, 3000]);
      if (s.endsWith("}")) s = s.slice(0, -1) + ', "deep": ' + "[".repeat(d) + "]".repeat(d) + "}";
    } else if (x < 0.9) s = s.replace('"', '"\\ud800');
    else {
      const i = r.int(s.length + 1);
      s = s.slice(0, i) + s.slice(i + r.randint(1, 8));
    }
    rewrite(e, Buffer.from(s, "latin1"));
  },
};

export const OPS = { ...zipOps, ...jsonOps };
export const ZIP_OPS = Object.keys(zipOps);
export const JSON_OPS = Object.keys(jsonOps);

// --- byte level -----------------------------------------------------------------

export function byteMutate(r, raw) {
  const d = Buffer.from(raw);
  if (!d.length) return r.bytes(r.randint(1, 64));
  const x = r.random();
  if (x < 0.4) {
    for (let i = r.randint(1, 8); i > 0; i--) d[r.int(d.length)] = r.int(256);
    return d;
  }
  if (x < 0.6) return d.subarray(0, r.int(d.length));
  if (x < 0.8) {
    const sig = Buffer.from(r.choice(["PK\x03\x04", "PK\x01\x02", "PK\x05\x06", "PK\x06\x06"]), "latin1");
    const pos = [];
    for (let i = d.indexOf(sig); i !== -1; i = d.indexOf(sig, i + 1)) pos.push(i);
    if (pos.length) {
      const p = r.choice(pos) + r.randint(4, 45);
      if (p + 4 <= d.length) d.writeUInt32LE(Number(BigInt.asUintN(32, BigInt(r.choice(SIZES)))), p);
    }
    return d;
  }
  const i = r.int(d.length);
  const j = Math.min(d.length, i + r.randint(1, 256));
  return Buffer.concat([d.subarray(0, i), d.subarray(i, j), d.subarray(i)]);
}

/** Apply the operators of `plan` ([name, seed] pairs) to a seed archive, then optional byte noise. */
export function mutate(raw, plan, byteSeed = null) {
  let out = raw;
  let a = null;
  try {
    a = parse(raw);
  } catch {
    a = null; // a non-ZIP seed only gets byte noise
  }
  if (a) {
    for (const [name, seed] of plan) OPS[name](rngFrom(seed), a);
    out = build(a);
  }
  if (byteSeed !== null) out = byteMutate(rngFrom(byteSeed), out);
  return out;
}
