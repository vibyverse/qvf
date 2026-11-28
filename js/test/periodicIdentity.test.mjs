// Periodic bond images and the symmetry site map (spec § 5.1, § 5.7), with
// the codes and JSON pointers of the Python reference validator
// (python/tests/test_periodic_identity.py).

import assert from "node:assert/strict";
import { test } from "node:test";

import { validate } from "../dist/index.js";
import { jsonBytes, jsonMember, pack } from "./helpers.mjs";

/** `[code, where]` of every error. */
function errors(/** @type {{ findings: { code: string, severity: string, where: string }[] }} */ report) {
  return report.findings.filter((f) => f.severity === "error").map((f) => [f.code, f.where]);
}

/** The set of error codes. */
const codes = (/** @type {Parameters<typeof errors>[0]} */ report) => [...new Set(errors(report).map(([c]) => c))].sort();

/**
 * An archive of `sections`, each `{ id, kind, payloads: { role: value }, ...extra }`.
 * @param {Array<{ id: string, kind: string, payloads: Record<string, unknown>, [k: string]: unknown }>} sections
 */
function archive(sections) {
  /** @type {{ name: string, data: Uint8Array }[]} */
  const entries = [];
  const manifestSections = sections.map(({ payloads, ...section }) => {
    /** @type {Record<string, unknown>} */
    const members = {};
    for (const [role, raw] of Object.entries(payloads)) {
      const value = role === "structure" ? withNumbers(raw) : raw;
      const name = `${section.id}/${role}.json`;
      const data = jsonBytes(value);
      entries.push({ name, data });
      members[role] = jsonMember(name, data);
    }
    return { ...section, members };
  });
  return pack({
    manifest: { qvf_version: 1, source: { program: "qvf-test", version: "0.0.1", calculation: "periodic identity" }, sections: manifestSections },
    entries,
  });
}

const Z = /** @type {Record<string, number>} */ ({ H: 1, C: 6, O: 8, Ti: 22, Po: 84 });

/** A structure payload whose atoms carry their `atomic_number`. */
const withNumbers = (/** @type {any} */ s) => ({ ...s, atoms: s.atoms.map((/** @type {any} */ a) => ({ atomic_number: Z[a.symbol], ...a })) });

/** @param {number[]} f @param {number[][]} lattice */
const cart = (f, lattice) => [0, 1, 2].map((r) => f[0] * lattice[0][r] + f[1] * lattice[1][r] + f[2] * lattice[2][r]);

// --- fixtures -------------------------------------------------------------------

const A_GRAPHENE = 2.46;
const GRAPHENE_LATTICE = [
  [A_GRAPHENE, 0, 0],
  [-A_GRAPHENE / 2, (A_GRAPHENE * Math.sqrt(3)) / 2, 0],
  [0, 0, 15],
];
/** Each carbon bonds to the other through three images. */
const GRAPHENE_BONDS = [
  { i: 0, j: 1, image: [0, 0, 0] },
  { i: 0, j: 1, image: [-1, 0, 0] },
  { i: 0, j: 1, image: [0, 1, 0] },
];
const graphene = () => ({
  atoms: [
    [1 / 3, 2 / 3, 0.5],
    [2 / 3, 1 / 3, 0.5],
  ].map((f) => ({ symbol: "C", position: cart(f, GRAPHENE_LATTICE) })),
  pbc: [true, true, false],
  lattice_vectors: GRAPHENE_LATTICE,
});

const A_CUBIC = 3.35;
const CUBIC = {
  atoms: [{ symbol: "Po", position: [0, 0, 0] }],
  pbc: [true, true, true],
  lattice_vectors: [
    [A_CUBIC, 0, 0],
    [0, A_CUBIC, 0],
    [0, 0, A_CUBIC],
  ],
};

/**
 * Bonds `pairs` on `structure`: the structure's own `bonds` member (`"own"`),
 * a standalone `bonds` section (`"bonds"`) or a `bond_orders` section.
 * @param {unknown} structure @param {Record<string, unknown>[]} pairs
 * @param {"own" | "bonds" | "bond_orders"} [where]
 */
function bondsArchive(structure, pairs, where = "own") {
  if (where === "own") return archive([{ id: "structure", kind: "structure", payloads: { structure, bonds: { pairs } } }]);
  const second =
    where === "bonds"
      ? { id: "bonds", kind: "bonds", structure_ref: "structure", payloads: { bonds: { pairs } } }
      : { id: "bond_orders", kind: "bond_orders", payloads: { bond_orders: { method: "mayer", pairs: pairs.map((p) => ({ order: 1.0, ...p })) } } };
  return archive([{ id: "structure", kind: "structure", payloads: { structure } }, second]);
}

// --- periodic bond images: archives that validate ----------------------------------

for (const where of /** @type {const} */ (["own", "bonds", "bond_orders"])) {
  test(`graphene bonds through three images validate (${where})`, async () => {
    const report = await validate(bondsArchive(graphene(), GRAPHENE_BONDS, where));
    assert.deepEqual(errors(report), []);
  });
}

test("a cubic cell bonds to its own images", async () => {
  const pairs = [
    { i: 0, j: 0, image: [1, 0, 0] },
    { i: 0, j: 0, image: [0, 1, 0] },
    { i: 0, j: 0, image: [0, 0, 1] },
  ];
  assert.deepEqual(errors(await validate(bondsArchive(CUBIC, pairs))), []);
});

test("molecule bonds without images keep their meaning, repeats included", async () => {
  const water = {
    atoms: [
      { symbol: "O", position: [0, 0, 0.117] },
      { symbol: "H", position: [0.757, 0, -0.47] },
      { symbol: "H", position: [-0.757, 0, -0.47] },
    ],
  };
  const pairs = [
    { i: 0, j: 1 },
    { i: 0, j: 2 },
    { i: 1, j: 0 },
  ];
  assert.deepEqual(errors(await validate(bondsArchive(water, pairs))), []);
});

test("the periodic axes need not lead", async () => {
  const slab = {
    atoms: [{ symbol: "H", position: [0, 0, 0] }],
    pbc: [false, true, true],
    lattice_vectors: [
      [20, 0, 0],
      [0, 3, 0],
      [0, 0, 3],
    ],
  };
  assert.deepEqual(errors(await validate(bondsArchive(slab, [{ i: 0, j: 0, image: [0, 0, 1] }]))), []);
});

// --- periodic bond images: archives that must fail -------------------------------

test("a self-bond with a zero image is refused", async () => {
  const report = await validate(bondsArchive(graphene(), [{ i: 1, j: 1, image: [0, 0, 0] }]));
  assert.deepEqual(errors(report), [["E_MEMBER_MALFORMED", "/sections/0/members/bonds"]]);
  assert.match(report.findings[0]?.message ?? "", /itself in the same cell/);
});

test("an image that is not three integers is the schema's finding only", async () => {
  for (const image of [[0.5, 0, 0], [1, 0]]) {
    assert.deepEqual(codes(await validate(bondsArchive(graphene(), [{ i: 0, j: 1, image }]))), ["E_PAYLOAD_SCHEMA"]);
  }
});

test("an image on a molecule is refused", async () => {
  const h2 = {
    atoms: [
      { symbol: "H", position: [0, 0, 0] },
      { symbol: "H", position: [0, 0, 0.74] },
    ],
  };
  const report = await validate(bondsArchive(h2, [{ i: 0, j: 1, image: [0, 0, 0] }], "bonds"));
  assert.deepEqual(errors(report), [["E_MEMBER_MALFORMED", "/sections/1/members/bonds"]]);
  assert.match(report.findings[0]?.message ?? "", /no periodic axis/);
});

test("an image along an aperiodic axis is refused", async () => {
  const report = await validate(bondsArchive(graphene(), [{ i: 0, j: 1, image: [0, 0, 1] }]));
  assert.deepEqual(errors(report), [["E_MEMBER_MALFORMED", "/sections/0/members/bonds"]]);
  assert.match(report.findings[0]?.message ?? "", /axis 2, which is not periodic/);
});

test("one bond written twice, or as its reversal, is refused", async () => {
  for (const second of [
    { i: 0, j: 1, image: [-1, 0, 0] },
    { i: 1, j: 0, image: [1, 0, 0] },
  ]) {
    const report = await validate(bondsArchive(graphene(), [GRAPHENE_BONDS[1], second]));
    assert.deepEqual(errors(report), [["E_MEMBER_MALFORMED", "/sections/0/members/bonds"]]);
    assert.match(report.findings[0]?.message ?? "", /pairs\[1\]\.image .* bond of pairs\[0\] again/);
  }
});

test("a self-bond and its reversal are one bond", async () => {
  const pairs = [
    { i: 0, j: 0, image: [1, 0, 0] },
    { i: 0, j: 0, image: [-1, 0, 0] },
  ];
  assert.deepEqual(errors(await validate(bondsArchive(CUBIC, pairs))), [["E_MEMBER_MALFORMED", "/sections/0/members/bonds"]]);
});

test("bond orders follow the same rules", async () => {
  const pairs = [GRAPHENE_BONDS[0], { ...GRAPHENE_BONDS[0], i: 1, j: 0 }];
  assert.deepEqual(errors(await validate(bondsArchive(graphene(), pairs, "bond_orders"))), [
    ["E_MEMBER_MALFORMED", "/sections/1/members/bond_orders"],
  ]);
});

test("bond images are checked at the references level only", async () => {
  const bytes = bondsArchive(graphene(), [{ i: 1, j: 1, image: [0, 0, 0] }]);
  assert.deepEqual(errors(await validate(bytes, { level: "semantics" })), []);
});

test("bonds resolve a structure that comes later; an unresolved structure_ref checks no image", async () => {
  const bonds = { pairs: [{ i: 1, j: 1, image: [0, 0, 0] }] };
  const later = archive([
    { id: "bonds", kind: "bonds", structure_ref: "structure", payloads: { bonds } },
    { id: "structure", kind: "structure", payloads: { structure: graphene() } },
  ]);
  assert.deepEqual(errors(await validate(later)), [["E_MEMBER_MALFORMED", "/sections/0/members/bonds"]]);
  const unresolved = archive([
    { id: "structure", kind: "structure", payloads: { structure: graphene() } },
    { id: "bonds", kind: "bonds", structure_ref: "missing", payloads: { bonds } },
  ]);
  assert.deepEqual(codes(await validate(unresolved)), ["E_REF_UNRESOLVED"]);
});

// --- the site map ------------------------------------------------------------------

/** An operation such as `x+1/2,-y,z` as `[R, t]` on fractional coordinates. */
function parseOp(/** @type {string} */ text) {
  /** @type {number[][]} */
  const rot = [];
  /** @type {number[]} */
  const shift = [];
  for (const part of text.split(",")) {
    const row = [0, 0, 0];
    let t = 0;
    let sign = 1;
    for (const m of part.matchAll(/([+-]?)(?:([xyz])|(\d+)\/(\d+))/g)) {
      sign = m[1] === "-" ? -1 : 1;
      if (m[2]) row["xyz".indexOf(m[2])] = sign;
      else t += (sign * Number(m[3])) / Number(m[4]);
    }
    rot.push(row);
    shift.push(t);
  }
  return /** @type {[number[][], number[]]} */ ([rot, shift]);
}

/**
 * Expand `sites` (label, symbol, fractional position) with the operations
 * into the reference cell, removing copies that coincide, and record each
 * atom's site, operation and translation.
 * @param {string[]} opTexts @param {[string, string, number[]][]} sites @param {number[][]} lattice
 */
function expand(opTexts, sites, lattice) {
  const ops = opTexts.map(parseOp);
  /** @type {{ symbol: string, position: number[] }[]} */
  const atoms = [];
  /** @type {{ site: number, operation: number, translation: number[] }[]} */
  const atomSites = [];
  sites.forEach(([, symbol, pos], s) => {
    /** @type {number[][]} */
    const seen = [];
    ops.forEach(([rot, shift], o) => {
      const raw = [0, 1, 2].map((r) => rot[r][0] * pos[0] + rot[r][1] * pos[1] + rot[r][2] * pos[2] + shift[r]);
      const n = raw.map((x) => -Math.floor(x + 1e-9) || 0);
      const frac = raw.map((x, r) => x + n[r]);
      if (seen.some((f) => Math.max(...[0, 1, 2].map((r) => Math.abs(frac[r] - f[r]))) < 1e-6)) return;
      seen.push(frac);
      atoms.push({ symbol, position: cart(frac, lattice) });
      atomSites.push({ site: s, operation: o, translation: n });
    });
  });
  return { atoms, operations: ops.map(([rotation, translation]) => ({ rotation, translation })), atomSites };
}

/** P 42/m n m (No. 136), the sixteen operations of the standard setting. */
const P42MNM = [
  "x,y,z", "-x,-y,z", "-y+1/2,x+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2",
  "-x+1/2,y+1/2,-z+1/2", "x+1/2,-y+1/2,-z+1/2", "y,x,-z", "-y,-x,-z",
  "-x,-y,-z", "x,y,-z", "y+1/2,-x+1/2,-z+1/2", "-y+1/2,x+1/2,-z+1/2",
  "x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,z+1/2", "-y,-x,z", "y,x,z",
];
const RUTILE_LATTICE = [
  [4.594, 0, 0],
  [0, 4.594, 0],
  [0, 0, 2.959],
];
/** Ti on 2a; O on 4f, deposited outside the reference cell so that its images carry translations. */
const RUTILE_SITES = /** @type {[string, string, number[]][]} */ ([
  ["Ti1", "Ti", [0, 0, 0]],
  ["O1", "O", [1.3049, 0.3049, 0]],
]);

/** P 1 21/n 1, a non-default setting of P 21/c (No. 14), in a skew cell. */
const P21N = ["x,y,z", "-x+1/2,y+1/2,-z+1/2", "-x,-y,-z", "x+1/2,-y+1/2,z+1/2"];
const BETA = (97 * Math.PI) / 180;
const P21N_LATTICE = [
  [5, 0, 0],
  [0, 6, 0],
  [7 * Math.cos(BETA), 0, 7 * Math.sin(BETA)],
];

/**
 * A structure and its `structure.symmetry` site map.
 * @param {string[]} opTexts @param {[string, string, number[]][]} sites @param {number[][]} lattice
 * @param {{ edit?: (d: any) => void, atomsEdit?: (a: any[]) => void, pbc?: boolean[] }} [options]
 */
function symmetryArchive(opTexts, sites, lattice, { edit, atomsEdit, pbc = [true, true, true] } = {}) {
  const { atoms, operations, atomSites } = expand(opTexts, sites, lattice);
  /** @type {Record<string, any>} */
  const data = { source: "deposited", operations, sites: sites.map(([label, , position]) => ({ label, position })), atom_sites: atomSites };
  edit?.(data);
  atomsEdit?.(atoms);
  return archive([
    { id: "structure", kind: "structure", payloads: { structure: { atoms, pbc, lattice_vectors: lattice } } },
    { id: "symmetry", kind: "structure.symmetry", payloads: { data } },
  ]);
}

const rutile = (/** @type {Parameters<typeof symmetryArchive>[3]} */ options = {}) => symmetryArchive(P42MNM, RUTILE_SITES, RUTILE_LATTICE, options);
const DATA = "/sections/1/members/data";

test("the rutile site map, a special position and nonzero translations, validates", async () => {
  const { atoms, atomSites } = expand(P42MNM, RUTILE_SITES, RUTILE_LATTICE);
  assert.equal(atoms.length, 6);
  assert.deepEqual(atomSites.map((e) => e.site), [0, 0, 1, 1, 1, 1]);
  assert.ok(atomSites.some((e) => e.translation.some((x) => x !== 0)));
  assert.deepEqual(errors(await validate(rutile())), []);
});

test("a non-default setting in a skew cell validates", async () => {
  const bytes = symmetryArchive(P21N, [["C1", "C", [0.11, 0.22, 0.33]]], P21N_LATTICE);
  assert.deepEqual(errors(await validate(bytes)), []);
});

test("a molecule's site map is Cartesian and has no translation", async () => {
  const water = {
    atoms: [
      { symbol: "O", position: [0, 0, 0.117] },
      { symbol: "H", position: [0.757, 0, -0.47] },
      { symbol: "H", position: [-0.757, 0, -0.47] },
    ],
  };
  const data = (/** @type {number[]} */ translation) => ({
    point_group: "C2v",
    source: "detected",
    symprec: 1e-3,
    operations: [{ rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]] }, { rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]] }],
    sites: [
      { label: "O1", position: [0, 0, 0.117] },
      { label: "H1", position: [0.757, 0, -0.47] },
    ],
    atom_sites: [{ site: 0, operation: 0 }, { site: 1, operation: 0 }, { site: 1, operation: 1, translation }],
  });
  const bytes = (/** @type {number[]} */ t) =>
    archive([
      { id: "structure", kind: "structure", payloads: { structure: water } },
      { id: "symmetry", kind: "structure.symmetry", payloads: { data: data(t) } },
    ]);
  assert.deepEqual(errors(await validate(bytes([0, 0, 0]))), []);
  assert.deepEqual(errors(await validate(bytes([1, 0, 0]))), [["E_MEMBER_MALFORMED", DATA]]);
});

test("symmetry without a site map is unchanged", async () => {
  const bytes = rutile({ edit: (d) => { delete d.atom_sites; delete d.sites; } });
  assert.deepEqual(errors(await validate(bytes)), []);
});

test("site_tolerance widens the position check", async () => {
  const atomsEdit = (/** @type {any[]} */ atoms) => { atoms[3].position[0] += 0.05; };
  assert.deepEqual(errors(await validate(rutile({ atomsEdit }))), [["E_MEMBER_MALFORMED", DATA]]);
  assert.deepEqual(errors(await validate(rutile({ atomsEdit, edit: (d) => { d.site_tolerance = 0.1; } }))), []);
});

/** @param {Record<string, unknown>} change */
const badEntry = (change) => (/** @type {any} */ d) => Object.assign(d.atom_sites[2], change);

test("site map rules", async () => {
  /** @type {[string, (d: any) => void, string][]} */
  const cases = [
    ["operation out of range", badEntry({ operation: 16 }), "E_MEMBER_MALFORMED"],
    ["site out of range", badEntry({ site: 2 }), "E_MEMBER_MALFORMED"],
    ["a wrong position", badEntry({ operation: 3 }), "E_MEMBER_MALFORMED"],
    ["fractional translation", badEntry({ translation: [0, 0, 0.5] }), "E_PAYLOAD_SCHEMA"],
    ["a wrong translation", badEntry({ translation: [1, 0, 0] }), "E_MEMBER_MALFORMED"],
    ["one entry short", (d) => d.atom_sites.pop(), "E_MEMBER_MALFORMED"],
    ["no sites", (d) => { delete d.sites; }, "E_MEMBER_MALFORMED"],
    ["no operations", (d) => { delete d.operations; }, "E_MEMBER_MALFORMED"],
    ["unknown source", (d) => { d.source = "guessed"; }, "E_PAYLOAD_SCHEMA"],
    ["a site without position", (d) => { delete d.sites[0].position; }, "E_PAYLOAD_SCHEMA"],
  ];
  for (const [name, edit, code] of cases) {
    assert.deepEqual(codes(await validate(rutile({ edit }))), [code], name);
  }
});

test("a translation along an aperiodic axis is refused", async () => {
  const bytes = rutile({
    pbc: [true, true, false],
    edit: (d) => { d.atom_sites[0].translation = [0, 0, 1]; },
    atomsEdit: (atoms) => { atoms[0].position[2] += RUTILE_LATTICE[2][2]; },
  });
  const report = await validate(bytes);
  assert.deepEqual(errors(report), [["E_MEMBER_MALFORMED", DATA]]);
  assert.match(report.findings.find((f) => f.severity === "error")?.message ?? "", /translation \[0, 0, 1\] is nonzero along axis 2/);
});

test("an edited structure disagrees with its site map", async () => {
  const report = await validate(rutile({ atomsEdit: (atoms) => { atoms[0].position = [1, 1, 1]; } }));
  assert.deepEqual(errors(report), [["E_MEMBER_MALFORMED", DATA]]);
  assert.match(report.findings.find((f) => f.severity === "error")?.message ?? "", /atom_sites\[0\] generates .* from atom 0, more than the site tolerance 0\.01/);
});

test("site maps are checked at the references level only", async () => {
  const bytes = rutile({ edit: badEntry({ operation: 16 }) });
  assert.deepEqual(errors(await validate(bytes, { level: "semantics" })), []);
});
