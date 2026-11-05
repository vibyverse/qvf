// The invariants every fuzz input is checked against (worker side), the
// counterpart of python/tests/fuzz/harness.py:
//
// - validate() never throws, every finding carries an Appendix B code of the
//   matching severity, and `ok` agrees with the findings;
// - QvfReader.fromBytes and the member reads refuse an input only with a
//   QvfError whose code is in Appendix B;
// - the typed arrays allocated while checking one input stay within what the
//   input justifies (largest allocation and sum, see below): a fixed base plus a small multiple of the bytes its
//   entries can really inflate to (the declared size, capped by the DEFLATE
//   limit of 1032:1). Time is bounded by the runner's watchdog.
//
// The allocation counter wraps the global typed-array constructors, so it
// must be installed before @qvf/core (and fflate, which captures Uint8Array
// at load) is imported: call install() first, then load().

export const MEMORY_BASE = 48 * 1024 ** 2;
// No single typed array may exceed the base plus twice what the input can
// inflate to: that catches allocation driven by a declared size that no
// stored byte backs. The sum over all allocations (two passes, plus fflate's
// short-lived scratch buffers per inflate step) only has to stay within a loose
// multiple, a guard against runaway churn.
export const LARGEST_FACTOR = 2;
export const TOTAL_FACTOR = 256;

export const SMALL = { max_member_bytes: 2 ** 24, max_total_bytes: 2 ** 26, max_json_bytes: 2 ** 22 };
export const DEFAULT = {};

let allocated = 0;
let largest = 0;

export function install() {
  for (const name of ["Uint8Array", "Int8Array", "Uint16Array", "Int16Array", "Uint32Array", "Int32Array",
    "Float32Array", "Float64Array", "BigInt64Array", "BigUint64Array"]) {
    const Base = globalThis[name];
    globalThis[name] = new Proxy(Base, {
      construct(target, args, newTarget) {
        if (typeof args[0] === "number") {
          const bytes = args[0] * Base.BYTES_PER_ELEMENT;
          allocated += bytes;
          if (bytes > largest) largest = bytes;
        }
        return Reflect.construct(target, args, newTarget === globalThis[name] ? target : newTarget);
      },
    });
  }
}

let lib;
export async function load() {
  lib = await import("../../dist/index.js");
  return lib;
}

/** Upper bound on the bytes the entries of `raw` can legitimately inflate to. */
export function inflatable(raw) {
  let total = raw.length;
  const sig = Buffer.from([0x50, 0x4b, 0x01, 0x02]);
  for (let p = raw.indexOf(sig); p !== -1 && p + 46 <= raw.length; p = raw.indexOf(sig, p + 4)) {
    const csize = raw.readUInt32LE(p + 20);
    let usize = raw.readUInt32LE(p + 24);
    if (usize === 0xffffffff || csize === 0xffffffff) usize = Math.min(usize, 2 ** 34);
    total += Math.min(usize, 1032 * Math.min(csize, raw.length) + 65536);
  }
  return total;
}

function where(err) {
  const line = String(err?.stack ?? "").split("\n").find((l) => l.includes("/dist/") || l.includes("fflate"));
  return line ? line.trim().replace(/^at /, "").replace(/\(.*\/(dist|node_modules)\//, "($1/") : "?";
}

const describe = (err) => `${err?.name ?? typeof err} at ${where(err)}: ${String(err?.message ?? err).slice(0, 200)}`;

async function checkValidator(raw, limits, problems) {
  let report;
  try {
    report = await lib.validate(raw, { limits });
  } catch (err) {
    problems.push(`validator threw ${describe(err)}`);
    return;
  }
  let errors = 0;
  for (const f of report.findings) {
    const severity = lib.FINDING_CODES[f.code];
    if (severity === undefined) problems.push(`validator: undocumented code ${f.code}`);
    else if (severity !== f.severity) problems.push(`validator: ${f.code} reported as ${f.severity}`);
    if (f.severity === "error") errors++;
  }
  if (report.ok !== (errors === 0)) problems.push("validator: ok disagrees with the findings");
}

async function checkReader(raw, limits, problems) {
  const guard = async (label, fn) => {
    try {
      return await fn();
    } catch (err) {
      if (!(err instanceof lib.QvfError) || lib.FINDING_CODES[err.code] === undefined) {
        problems.push(`${label} threw ${describe(err)}`);
      }
      return undefined;
    }
  };
  const reader = await guard("QvfReader.fromBytes", () => lib.QvfReader.fromBytes(raw, limits));
  if (!reader) return;
  const seen = new Set();
  for (const s of reader.sections) {
    // Members are addressed by section id, which resolves to the first section
    // of that id (duplicate ids are a validator finding, E_SECTION_ID_DUPLICATE).
    if (seen.has(s.id)) continue;
    seen.add(s.id);
    for (const [role, m] of Object.entries(s.members)) {
      if (m && m.format === "json") await guard("readJson", () => reader.readJson(s.id, role));
      else if (m && m.format === "binary" && "dtype" in m) await guard("readBinary", () => reader.readBinary(s.id, role));
      else await guard("readBytes", () => reader.readBytes(s.id, role));
    }
  }
}

/** Check one input; returns { problems, allocated, largest }. */
export async function check(raw, { limits = SMALL, reader = true } = {}) {
  const problems = [];
  allocated = 0;
  largest = 0;
  await checkValidator(raw, limits, problems);
  if (reader) await checkReader(raw, limits, problems);
  const backed = inflatable(raw);
  if (largest > MEMORY_BASE + LARGEST_FACTOR * backed) {
    problems.push(`memory: one allocation of ${largest} bytes over the bound ${MEMORY_BASE + LARGEST_FACTOR * backed}`);
  }
  if (allocated > MEMORY_BASE + TOTAL_FACTOR * backed) {
    problems.push(`memory: allocated ${allocated} bytes in all, over the bound ${MEMORY_BASE + TOTAL_FACTOR * backed}`);
  }
  return { problems, allocated, largest };
}
