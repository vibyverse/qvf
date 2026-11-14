// The energy_decomposition kind (spec § 5.16): addEnergyDecomposition, the volume.difference fields that name
// an orbital pair (§ 5.2), and the validator rules of § 6.

import assert from "node:assert/strict";
import { test } from "node:test";

import { QvfReader, Writer, validate } from "../dist/index.js";

import { FIXED, SOURCE } from "./writerFixtures.mjs";

const DIMER = [
  { symbol: "O", position: [0, 0, 0] }, { symbol: "H", position: [0, 0, 0.96] },
  { symbol: "H", position: [0, 0.929, -0.24] }, { symbol: "O", position: [0, 0, 2.95] },
  { symbol: "H", position: [0, 0.759, 3.538] }, { symbol: "H", position: [0, -0.759, 3.538] },
];

/** A Ziegler-Rauk payload with one NOCV pair. @returns {any} */
const ZR = () => ({
  scheme: "ziegler-rauk",
  variant: "eda-nocv",
  fragments: [{ label: "donor", atoms: [0, 1, 2] }, { label: "acceptor", atoms: [3, 4, 5] }],
  terms: [
    { name: "interaction", energy: -0.008 },
    { name: "electrostatic", energy: -0.0125 },
    { name: "pauli", energy: 0.0132 },
    { name: "orbital", energy: -0.0073 },
    { name: "dispersion", energy: { value: -0.0014 * 627.5094740631, units: "kcal_per_mol" } },
  ],
  sum_rules: [{ total: 0, parts: [1, 2, 3, 4] }],
  orbital_pairs: {
    type: "nocv", eigenvalue_convention: "per_spin",
    pairs: [{ index: 1, eigenvalue: 0.045, energy: -0.0049, density_ref: "nocv_1" }],
  },
});

/** @param {any} [payload] */
async function archive(payload = ZR()) {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(DIMER);
  const eye = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]];
  await w.addVolume("volume.difference", { data: new Float32Array(27), shape: [3, 3, 3] }, [0, 0, 0], eye, {
    sectionId: "nocv_1", decompositionRef: "energy_decomposition", decompositionRole: "orbital_pair", pairIndex: 1,
  });
  await w.addEnergyDecomposition(payload);
  return w.toBytes();
}

/** An archive whose energy_decomposition payload skipped addEnergyDecomposition's checks (written as a raw section). */
/** @param {(p: any) => void} edit @param {number} [pairIndex] */
async function raw(edit, pairIndex = 1) {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(DIMER);
  const eye = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]];
  await w.addVolume("volume.difference", { data: new Float32Array(27), shape: [3, 3, 3] }, [0, 0, 0], eye, {
    sectionId: "nocv_1", decompositionRef: "energy_decomposition", decompositionRole: "orbital_pair", pairIndex,
  });
  const p = ZR();
  for (const t of p.terms) if (typeof t.energy === "number") t.energy = { value: t.energy, units: "hartree" };
  for (const q of p.orbital_pairs.pairs) q.energy = { value: q.energy, units: "hartree" };
  edit(p);
  w.addSection({ id: "energy_decomposition", kind: "energy_decomposition" });
  await w.addJson("energy_decomposition", "decomposition", p, { path: "energy_decomposition/decomposition.json" });
  return w.toBytes();
}

/** @param {Uint8Array} bytes */
const codes = async (bytes) => new Set((await validate(bytes)).findings.map((f) => f.code));

test("a Ziegler-Rauk decomposition with a NOCV density validates", async () => {
  const bytes = await archive();
  const report = await validate(bytes);
  assert.deepEqual(report.findings.map((f) => `${f.code} ${f.message}`), []);
  const q = QvfReader.fromBytes(bytes);
  const p = /** @type {any} */ (await q.readJson("energy_decomposition", "decomposition"));
  assert.deepEqual(p.terms[0].energy, { value: -0.008, units: "hartree" });
  assert.equal(q.section("nocv_1")?.pair_index, 1);
});

test("addEnergyDecomposition refuses what the validator would reject", async () => {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(DIMER);
  const broken = ZR();
  delete broken.orbital_pairs;
  broken.sum_rules = [{ total: 0, parts: [1, 2, 3] }];
  await assert.rejects(w.addEnergyDecomposition(broken), /add up/);
  const overlap = ZR();
  delete overlap.orbital_pairs;
  overlap.fragments[1].atoms.push(0);
  await assert.rejects(w.addEnergyDecomposition(overlap), /disjoint/);
  await assert.rejects(w.addEnergyDecomposition(ZR()), /density_ref/);
});

test("the validator reports the § 6 rules", async () => {
  /** @type {Array<[(p: any) => void, import("../dist/index.js").FindingCode]>} */
  const cases = [
    [(p) => p.sum_rules[0].parts.pop(), "E_SUM_RULE"],
    [(p) => p.fragments[1].atoms.push(9), "E_MEMBER_MALFORMED"],
    [(p) => p.terms.push({ name: "interaction", energy: { value: 0, units: "hartree" } }), "E_MEMBER_MALFORMED"],
    [(p) => { p.terms[1].energy.units = "bohr"; }, "E_UNIT_CONTRADICTED"],
    [(p) => { p.orbital_pairs.pairs[0].density_ref = "missing"; }, "E_REF_UNRESOLVED"],
  ];
  for (const [edit, code] of cases) {
    assert.ok((await codes(await raw(edit))).has(code), code);
  }
  assert.deepEqual([...(await codes(await raw(() => {})))], []);
  assert.ok((await codes(await raw(() => {}, 4))).has("E_REF_UNRESOLVED"));
});
