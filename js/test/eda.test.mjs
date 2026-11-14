// The eda kind (spec § 5.16): addEda, the volume.difference fields that name
// an orbital pair (§ 5.2), and the validator rules of § 6.

import assert from "node:assert/strict";
import { test } from "node:test";

import { QvfReader, Writer, validate } from "../dist/index.js";

import { FIXED, SOURCE } from "./writerFixtures.mjs";

const DIMER = [
  ["O", [0, 0, 0]], ["H", [0, 0, 0.96]], ["H", [0, 0.929, -0.24]],
  ["O", [0, 0, 2.95]], ["H", [0, 0.759, 3.538]], ["H", [0, -0.759, 3.538]],
].map(([symbol, position]) => ({ symbol, position }));

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

async function archive(payload = ZR()) {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(DIMER);
  const eye = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]];
  await w.addVolume("volume.difference", { data: new Float32Array(27), shape: [3, 3, 3] }, [0, 0, 0], eye, {
    sectionId: "nocv_1", edaRef: "eda", edaRole: "orbital_pair", pairIndex: 1,
  });
  await w.addEda(payload);
  return w.toBytes();
}

/** An archive whose eda payload skipped addEda's checks (written as a raw section). */
async function raw(edit, pairIndex = 1) {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(DIMER);
  const eye = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]];
  await w.addVolume("volume.difference", { data: new Float32Array(27), shape: [3, 3, 3] }, [0, 0, 0], eye, {
    sectionId: "nocv_1", edaRef: "eda", edaRole: "orbital_pair", pairIndex,
  });
  const p = ZR();
  for (const t of p.terms) if (typeof t.energy === "number") t.energy = { value: t.energy, units: "hartree" };
  for (const q of p.orbital_pairs.pairs) q.energy = { value: q.energy, units: "hartree" };
  edit(p);
  w.addSection({ id: "eda", kind: "eda" });
  await w.addJson("eda", "decomposition", p, { path: "eda/decomposition.json" });
  return w.toBytes();
}

const codes = async (bytes) => new Set((await validate(bytes)).findings.map((f) => f.code));

test("a Ziegler-Rauk decomposition with a NOCV density validates", async () => {
  const bytes = await archive();
  const report = await validate(bytes);
  assert.deepEqual(report.findings.map((f) => `${f.code} ${f.message}`), []);
  const q = QvfReader.fromBytes(bytes);
  const p = await q.readJson("eda", "decomposition");
  assert.deepEqual(p.terms[0].energy, { value: -0.008, units: "hartree" });
  assert.equal(q.section("nocv_1")?.pair_index, 1);
});

test("addEda refuses what the validator would reject", async () => {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(DIMER);
  const broken = ZR();
  delete broken.orbital_pairs;
  broken.sum_rules = [{ total: 0, parts: [1, 2, 3] }];
  await assert.rejects(w.addEda(broken), /add up/);
  const overlap = ZR();
  delete overlap.orbital_pairs;
  overlap.fragments[1].atoms.push(0);
  await assert.rejects(w.addEda(overlap), /disjoint/);
  await assert.rejects(w.addEda(ZR()), /density_ref/);
});

test("the validator reports the § 6 rules", async () => {
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
