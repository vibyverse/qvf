// The writer: typed builders, round trips through the reader, the validator,
// deterministic output, ZIP64, streaming output, and refusals.

import assert from "node:assert/strict";
import { test } from "node:test";

import { QvfError, QvfReader, Writer, ZipArchive, ZipWriter, validate } from "../dist/index.js";

import { EYE, FIXED, H2, SHELLS, SOURCE, fullArchive } from "./writerFixtures.mjs";

/** A predicate for assert.throws / rejects: a QvfError with this code. @param {string} code */
const hasCode = (code) => (/** @type {unknown} */ err) => err instanceof QvfError && err.code === code;

/** @param {Uint8Array} bytes */
async function assertValid(bytes) {
  const report = await validate(bytes);
  assert.deepEqual(
    report.findings.map((f) => `${f.code} ${f.where} ${f.message}`),
    [],
  );
  assert.equal(report.ok, true);
}

test("every builder produces an archive the validator accepts", async () => {
  const bytes = await (await fullArchive()).toBytes();
  await assertValid(bytes);
});

test("round trip through the reader", async () => {
  const bytes = await (await fullArchive()).toBytes();
  const q = QvfReader.fromBytes(bytes);
  /** @param {string} id @param {string} role @returns {Promise<any>} */
  const json = (id, role) => q.readJson(id, role);
  assert.equal(q.manifest.qvf_version, 1);
  assert.equal(q.manifest.byte_order, "little");
  assert.deepEqual(q.manifest.source, SOURCE);
  assert.equal(q.manifest.metadata?.created_utc, FIXED.timestamp);
  assert.equal(q.manifest.metadata?.uuid, FIXED.uuid);
  assert.equal(q.manifest.metadata?.title, "H2 test archive");
  assert.equal(q.manifest.provenance?.method, "rhf");

  const s = await json("structure", "structure");
  assert.deepEqual(s.atoms.map((/** @type {any} */ a) => [a.symbol, a.atomic_number]), [["H", 1], ["H", 1]]);
  assert.deepEqual(s.pbc, [false, false, false]);
  assert.equal(q.section("structure")?.label, "H2");
  assert.equal(q.section("structure")?.charge, 0);
  const nacl = await json("nacl", "structure");
  assert.equal(nacl.atoms[1].symbol, "Cl");
  assert.equal(nacl.dimensionality, 3);

  const coords = await q.readBinary("trajectory", "coords");
  assert.equal(coords.dtype, "float32");
  assert.deepEqual(coords.shape, [3, 2, 3]);
  assert.ok(Math.abs(coords.data[5] - 0.74) < 1e-6);
  const meta = await json("trajectory", "metadata");
  assert.equal(meta.precision, "float32");
  assert.deepEqual([...(await q.readBinary("trajectory", "times")).data], [0, 0.5, 1]);

  const C = await q.readBinary("wf", "mo_coefficients");
  assert.deepEqual([C.dtype, C.shape, [...C.data]], ["float64", [2, 2], [0.55, 0.55, 1.2, -1.2]]);
  const mo = await json("wf", "mo_metadata");
  assert.deepEqual([mo.spin, mo.n_mo, mo.n_ao, mo.occupations], ["restricted", 2, 2, [2, 0]]);
  const mou = await json("wf_u", "mo_metadata");
  assert.deepEqual([mou.spin, mou.alpha.n_mo, mou.beta.n_mo, mou.beta.symmetry_labels], ["unrestricted", 2, 1, ["sg"]]);
  assert.deepEqual((await q.readBinary("wf_u", "mo_coefficients_beta")).shape, [1, 2]);
  const basis = await json("wf", "basis");
  assert.deepEqual([basis.structure_ref, basis.n_ao, basis.n_electrons], ["structure", 2, 2]);

  const rho = await q.readBinary("rho", "data");
  assert.deepEqual([rho.dtype, rho.shape, rho.data[7]], ["float32", [2, 2, 2], 0.25]);
  const orb = await q.readBinary(q.sectionsOf("volume.orbital")[0].id, "data");
  assert.deepEqual([orb.dtype, [...orb.data.slice(0, 3)]], ["float64", [1, -1, 0.5]]);
  assert.deepEqual((await json("rho_t", "axis")).values, [0, 1, 2]);
  assert.deepEqual((await json("rho_t", "grid")).shape, [2, 2, 2]);
  assert.equal(q.section("drho")?.operand_a, "rho");
  assert.equal(q.section("rho")?.group, "densities");

  const states = await json("excited_states", "states");
  assert.deepEqual(states.states[0].excitation_energy, { value: 12.1, units: "ev" });
  assert.equal(q.section("excited_states")?.wavefunction_ref, "wf");
  const ir = await json("spectra_ir", "spectrum");
  assert.deepEqual([ir.frequencies, ir.x_units], [[4400], "cm-1"]);
  const scheme = await json("reaction_scheme", "scheme");
  assert.deepEqual(scheme.steps[0].reactants, [{ structure_ref: "structure" }]);

  assert.equal(new TextDecoder().decode(await q.readBytes("notes", "log")), "line one\nline two\n");
  assert.equal(q.section("notes")?.members.log.media_type, "text/plain");
  assert.deepEqual((await json("notes", "extra")).values, [1, 2]);
  const ids = await q.readBinary("notes", "ids");
  assert.deepEqual([ids.dtype, [...ids.data]], ["int64", [1n, -2n, 2n ** 62n]]);
});

test("member specs carry sha256 and size; manifest.json is stored and last", async () => {
  const bytes = await (await fullArchive()).toBytes();
  const zip = ZipArchive.open(bytes);
  const last = zip.entries[zip.entries.length - 1];
  assert.equal(last.name, "manifest.json");
  assert.equal(last.method, 0);
  for (const e of zip.entries.slice(0, -1)) assert.equal(e.method, 8, e.name);
  const q = QvfReader.fromBytes(bytes);
  for (const s of q.sections) {
    for (const m of Object.values(s.members)) {
      assert.match(m.sha256, /^[0-9a-f]{64}$/);
      assert.equal(m.size, zip.find(m.path)?.size, m.path);
    }
  }
});

test("output is deterministic for a fixed timestamp and uuid", async () => {
  const a = await (await fullArchive()).toBytes();
  const b = await (await fullArchive()).toBytes();
  assert.deepEqual(a, b);
  const c = await (await fullArchive({ timestamp: "2026-01-02T03:04:07Z" })).toBytes();
  assert.notDeepEqual(a, c);
  // Without automatic metadata the entries carry the DOS epoch and no time stamp is written.
  const make = async () => {
    const w = new Writer(SOURCE, { autoMetadata: false });
    await w.addStructure(H2);
    return w.toBytes();
  };
  const d = await make();
  assert.deepEqual(d, await make());
  assert.equal(QvfReader.fromBytes(d).manifest.metadata, undefined);
});

test("automatic metadata is a v4 uuid and the current time", async () => {
  const w = new Writer(SOURCE);
  await w.addStructure(H2);
  const m = QvfReader.fromBytes(await w.toBytes()).manifest.metadata;
  assert.match(m?.uuid ?? "", /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/);
  assert.ok(Math.abs(Date.parse(m?.created_utc ?? "") - Date.now()) < 60_000);
});

test("forced ZIP64 archives read back and validate", async () => {
  const bytes = await (await fullArchive({ forceZip64: true })).toBytes();
  await assertValid(bytes);
  const view = new DataView(bytes.buffer, bytes.byteOffset);
  assert.equal(view.getUint32(bytes.length - 22 - 20, true), 0x07064b50, "ZIP64 locator before the EOCD");
  const q = QvfReader.fromBytes(bytes);
  assert.deepEqual([...(await q.readBinary("wf", "mo_coefficients")).data], [0.55, 0.55, 1.2, -1.2]);
});

test("more than 65535 entries switch to ZIP64 end records", async () => {
  const zw = new ZipWriter();
  const n = 65536;
  const one = Uint8Array.of(7);
  for (let i = 0; i < n; i++) await zw.add(`e/${i}`, one, 0);
  await zw.finish();
  const bytes = zw.bytes();
  const view = new DataView(bytes.buffer, bytes.byteOffset);
  assert.equal(view.getUint16(bytes.length - 22 + 10, true), 0xffff);
  const zip = ZipArchive.open(bytes);
  assert.equal(zip.entries.length, n);
  assert.deepEqual(zip.read(`e/${n - 1}`), one);
});

test("streaming to a sink gives the same bytes as the in-memory archive", async () => {
  /** @type {Uint8Array[]} */
  const chunks = [];
  let calls = 0;
  const w = await fullArchive({
    sink: async (chunk) => {
      calls++;
      await new Promise((r) => setTimeout(r, 0));
      chunks.push(chunk.slice());
    },
  });
  await w.finish();
  const streamed = new Uint8Array(chunks.reduce((n, c) => n + c.length, 0));
  let o = 0;
  for (const c of chunks) {
    streamed.set(c, o);
    o += c.length;
  }
  assert.ok(calls > 10);
  assert.equal(w.length, streamed.length);
  assert.deepEqual(streamed, await (await fullArchive()).toBytes());
  assert.throws(() => w.bytes(), /sink/);
});

test("a Blob holds the same bytes", async () => {
  const w = await fullArchive();
  const blob = await w.toBlob();
  assert.equal(blob.type, "application/zip");
  assert.deepEqual(new Uint8Array(await blob.arrayBuffer()), w.bytes());
});

test("store compression and unawaited calls keep call order", async () => {
  const w = new Writer(SOURCE, { ...FIXED, compression: "store" });
  w.addSection({ id: "v", kind: "x_test.v" });
  const jobs = [];
  for (let i = 0; i < 20; i++) jobs.push(w.addBinary("v", `m${i}`, new Float64Array([i])));
  await Promise.all(jobs);
  const bytes = await w.toBytes();
  const zip = ZipArchive.open(bytes);
  assert.deepEqual(
    zip.entries.map((e) => e.name),
    [...Array.from({ length: 20 }, (_, i) => `v/m${i}.bin`), "manifest.json"],
  );
  assert.ok(zip.entries.every((e) => e.method === 0));
  await assertValid(bytes);
});

test("binary members are little-endian and converted to the dtype", async () => {
  const w = new Writer(SOURCE, FIXED);
  w.addSection({ id: "b", kind: "x_test.b" });
  await w.addBinary("b", "u16", [1, 258], { dtype: "uint16" });
  await w.addBinary("b", "f32", new Float64Array([1.5, -2]), { dtype: "float32", shape: [2, 1] });
  await w.addBinary("b", "i32", new Int32Array([-1, 2, 3, 4, 5, 6]), { shape: [2, 3] });
  const bytes = await w.toBytes();
  const zip = ZipArchive.open(bytes);
  assert.deepEqual([...zip.read("b/u16.bin")], [1, 0, 2, 1]);
  const q = QvfReader.fromBytes(bytes);
  const f32 = await q.readBinary("b", "f32");
  assert.deepEqual([f32.dtype, f32.shape, [...f32.data]], ["float32", [2, 1], [1.5, -2]]);
  assert.deepEqual((await q.readBinary("b", "i32")).shape, [2, 3]);
});

test("payloads that fail the schema are refused before anything is written", async () => {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(H2);
  await assert.rejects(
    w.addJson("structure", "bonds", { pairs: [{ i: 0 }] }),
    (err) => err instanceof QvfError && err.code === "E_PAYLOAD_SCHEMA" && /BondsPayload/.test(err.message) && /\/pairs\/0/.test(err.message),
  );
  await assert.rejects(
    w.addSpectrum("spectra.ir", { payload: { frequencies: [1], intensities: ["x"] } }),
    (err) => err instanceof QvfError && err.code === "E_PAYLOAD_SCHEMA" && /\/intensities\/0/.test(err.message),
  );
  assert.ok(!w.sectionIds.includes("spectra_ir"), "a refused builder adds no section");
  // The writer is still usable, and the refused path is free again.
  await w.addJson("structure", "bonds", { pairs: [{ i: 0, j: 1, order: 1 }] });
  await assertValid(await w.toBytes());
});

test("a manifest that fails the schema is refused by finish, which can be retried", async () => {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(H2);
  w.setProvenance({ n_threads: 0 });
  await assert.rejects(w.finish(), (err) => err instanceof QvfError && err.code === "E_SCHEMA" && /\/provenance\/n_threads/.test(err.message));
  w.setProvenance({ n_threads: 4 });
  await assertValid(await w.toBytes());
});

test("builder arguments are checked", async () => {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(H2);
  await assert.rejects(w.addStructure(H2), hasCode("E_SECTION_ID_DUPLICATE"));
  await assert.rejects(w.addStructure(H2, { sectionId: "camera" }), hasCode("E_RESERVED_ID"));
  await assert.rejects(w.addStructure(H2, { sectionId: "p", pbc: [true, false, false] }), hasCode("E_STRUCTURE_PERIODICITY"));
  await assert.rejects(w.addStructure(H2, { sectionId: "q", charge: 0, multiplicity: 2 }), hasCode("E_CHARGE_MULTIPLICITY"));
  await assert.rejects(w.addStructure([{ symbol: "H", position: [0, 0] }], { sectionId: "r" }), /position/);
  await assert.rejects(w.addTrajectory({ symbols: ["H", "H"] }, [[[0, 0, 0], [0, 0]]]), /ragged/);
  await assert.rejects(
    w.addTrajectory({ symbols: ["H", "H"], energies: [-1] }, [H2.map((a) => a.position)], { energies: [-1.5] }),
    /agree within 1e-10/,
  );
  await assert.rejects(w.addWavefunctionGto(SHELLS, { moCoefficients: [[1, 0, 0]] }), /n_mo, 2/);
  await assert.rejects(w.addWavefunctionGto(SHELLS, { moCoefficients: [[1, 0]], occupations: [2], nElectrons: 1 }), /nElectrons/);
  await assert.rejects(w.addDensity(new Float32Array(8), [0, 0, 0], EYE), /shape/);
  await assert.rejects(w.addVolume("volume.nope", [[[0]]], [0, 0, 0], EYE), hasCode("E_KIND_UNKNOWN"));
  await assert.rejects(w.addDensity([[[0]]], [0, 0, 0], EYE, { structureRef: "nowhere" }), hasCode("E_REF_UNRESOLVED"));
  await assert.rejects(w.addSpectrum("spectra.ir", { frequencies: [1, 2], intensities: [1] }), /entries of frequencies/);
  await assert.rejects(w.addSpectrum("spectra.ir", { frequencies: [1], intensities: [1], xUnits: "furlong" }), /unit name/);
  await assert.rejects(w.addReactionScheme([{ id: "a", reactants: ["structure"], products: ["missing"] }]), hasCode("E_REF_UNRESOLVED"));
  await assert.rejects(w.addExcitedStates([{ index: 2, excitation_energy: 1 }, { index: 1, excitation_energy: 2 }]), /ascend/);
  w.addSection({ id: "x", kind: "x_test.x" });
  await assert.rejects(w.addBinary("x", "big", [300], { dtype: "uint8" }), /range/);
  await assert.rejects(w.addBinary("x", "nan", [1, 2, 3], { shape: [2] }), hasCode("E_BINARY_SIZE"));
  await assert.rejects(w.addJson("x", "nan", { v: NaN }), /NaN/);
  await assert.rejects(w.addJson("x", "p", {}, { path: "../evil.json" }), hasCode("E_MEMBER_PATH"));
  await assert.rejects(w.addJson("x", "p", {}, { path: "structure/structure.json" }), hasCode("E_ZIP_DUPLICATE_ENTRY"));
  assert.throws(() => w.addSection({ id: "k", kind: "nonsense" }), hasCode("E_KIND_UNKNOWN"));
  assert.throws(() => new Writer({ program: "", version: "1" }), hasCode("E_SOURCE"));
  // Every refusal above left the archive consistent.
  await assertValid(await w.toBytes());
  await assert.rejects(w.addStructure(H2, { sectionId: "late" }), /finished/);
});

test("a refused member leaves its path and role free", async () => {
  const w = new Writer(SOURCE, FIXED);
  w.addSection({ id: "x", kind: "x_test.x" });
  await assert.rejects(w.addJson("x", "doc", { v: Infinity }, { path: "x/doc.json" }), /infinity/);
  await assert.rejects(w.addBinary("x", "arr", [0.5], { dtype: "int32", path: "x/arr.bin" }), /integer/);
  await assert.rejects(w.addOpaque("x", "log", "text", { path: "x/log.txt", mediaType: "text" }), /media type/);
  await w.addJson("x", "doc", { v: 1 }, { path: "x/doc.json" });
  await w.addBinary("x", "arr", [1], { dtype: "int32", path: "x/arr.bin" });
  await w.addOpaque("x", "log", "text", { path: "x/log.txt", mediaType: "text/plain" });
  await assert.rejects(w.addJson("x", "doc", {}), /already has a member/);
  await assertValid(await w.toBytes());
});
