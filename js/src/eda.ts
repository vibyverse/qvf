/**
 * Rules of an `energy_decomposition` payload that JSON Schema cannot express (spec § 5.16),
 * shared by the validator and the writer. Mirrors `python/qvf/_eda.py`.
 *
 * License: Apache-2.0.
 */

type Obj = Record<string, unknown>;

/** The finding codes the payload rules report. */
export type EdaCode = "E_MEMBER_MALFORMED" | "E_UNIT_CONTRADICTED" | "E_SUM_RULE";

const isObject = (v: unknown): v is Obj => typeof v === "object" && v !== null && !Array.isArray(v);
const isInt = (v: unknown): v is number => typeof v === "number" && Number.isInteger(v);
const isNum = (v: unknown): v is number => typeof v === "number" && Number.isFinite(v);
const list = (v: unknown): unknown[] => (Array.isArray(v) ? v : []);
const has = (o: Obj, k: string): boolean => Object.hasOwn(o, k);
const show = (v: unknown): string => (typeof v === "string" ? `'${v}'` : JSON.stringify(v) ?? String(v));

/** Hartree per unit of every energy unit a decomposition may use (§ 4.3). */
export const ENERGY_UNITS: ReadonlyMap<string, number> = new Map([
  ["hartree", 1.0],
  ["Eh", 1.0],
  ["Ha", 1.0],
  ["Hartree", 1.0],
  ["ev", 1.0 / 27.211386245988],
  ["eV", 1.0 / 27.211386245988],
  ["kcal_per_mol", 1.0 / 627.5094740631],
  ["kj_per_mol", 1.0 / 2625.4996394799],
]);

/** Default sum-rule tolerance, hartree (§ 5.16). */
export const DEFAULT_SUM_TOLERANCE_EH = 1e-6;

const SECTION = "spec § 5.16";

/** An `{value, units}` energy in hartree, or null when it is not one. */
export function energyEh(energy: unknown): number | null {
  if (!isObject(energy) || !isNum(energy.value)) return null;
  const factor = typeof energy.units === "string" ? ENERGY_UNITS.get(energy.units) : undefined;
  return factor === undefined ? null : energy.value * factor;
}

function energies(payload: Obj): Array<[string, unknown]> {
  const found: Array<[string, unknown]> = [];
  list(payload.fragments).forEach((f, k) => {
    if (isObject(f) && has(f, "energy")) found.push([`fragments[${k}].energy`, f.energy]);
  });
  if (isObject(payload.energies) && has(payload.energies, "system")) {
    found.push(["energies.system", payload.energies.system]);
  }
  list(payload.terms).forEach((t, k) => {
    if (isObject(t) && has(t, "energy")) found.push([`terms[${k}].energy`, t.energy]);
  });
  list(payload.sum_rules).forEach((r, k) => {
    if (isObject(r) && has(r, "tolerance")) found.push([`sum_rules[${k}].tolerance`, r.tolerance]);
  });
  const pairs = payload.orbital_pairs;
  if (isObject(pairs)) {
    if (has(pairs, "total")) found.push(["orbital_pairs.total", pairs.total]);
    list(pairs.pairs).forEach((p, k) => {
      if (isObject(p) && has(p, "energy")) found.push([`orbital_pairs.pairs[${k}].energy`, p.energy]);
    });
  }
  return found;
}

/** `"<index>/<spin>"` of every orbital pair, or null without a pair table. */
export function pairKeys(payload: unknown): Set<string> | null {
  if (!isObject(payload) || !isObject(payload.orbital_pairs)) return null;
  const keys = new Set<string>();
  for (const p of list(payload.orbital_pairs.pairs)) {
    if (isObject(p) && isInt(p.index)) keys.add(`${p.index}/${typeof p.spin === "string" ? p.spin : "restricted"}`);
  }
  return keys;
}

/**
 * `[code, message]` for every rule of § 5.16 an `energy_decomposition` payload breaks.
 * `nAtoms` is the atom count of the section's structure when it resolves.
 */
export function edaPayloadErrors(payload: unknown, nAtoms: number | null = null): Array<[EdaCode, string]> {
  if (!isObject(payload)) return [];
  const errors: Array<[EdaCode, string]> = [];
  const bad = (msg: string): void => {
    errors.push(["E_MEMBER_MALFORMED", `${msg} (${SECTION})`]);
  };

  // --- fragments --------------------------------------------------------
  const fragments = list(payload.fragments);
  const nFrag = fragments.length;
  const atomsOf: Array<Set<number> | null> = fragments.map((f) =>
    isObject(f) && Array.isArray(f.atoms) ? new Set(f.atoms.filter(isInt)) : null,
  );
  const top: number[] = [];
  fragments.forEach((f, k) => {
    if (isObject(f) && !has(f, "parent")) top.push(k);
  });
  if (Array.isArray(payload.fragments) && top.length < 2) {
    bad(`fragments has ${top.length} interacting fragment(s) (no parent); a decomposition needs at least two`);
  }
  const owner = new Map<number, number>();
  for (const k of top) {
    for (const a of [...(atomsOf[k] ?? [])].sort((x, y) => x - y)) {
      const prev = owner.get(a);
      if (prev !== undefined) bad(`fragments[${k}] and fragments[${prev}] both contain atom ${a}; interacting fragments are disjoint`);
      else owner.set(a, k);
    }
  }
  fragments.forEach((f, k) => {
    if (!isObject(f)) return;
    if (nAtoms !== null) {
      for (const a of [...(atomsOf[k] ?? [])].sort((x, y) => x - y)) {
        if (a < 0 || a >= nAtoms) bad(`fragments[${k}].atoms names atom ${a}, but the structure has ${nAtoms} atoms`);
      }
    }
    if (has(f, "parent")) {
      const p = f.parent;
      if (!isInt(p) || p < 0 || p >= nFrag || p === k) {
        bad(`fragments[${k}].parent ${show(p)} is not the index of another fragment`);
      } else if (!top.includes(p)) {
        bad(`fragments[${k}].parent ${p} is itself a sub-fragment; a sub-fragment subdivides an interacting fragment`);
      } else {
        const mine = atomsOf[k];
        const theirs = atomsOf[p];
        if (mine && theirs) {
          const extra = [...mine].filter((a) => !theirs.has(a)).sort((x, y) => x - y);
          if (extra.length > 0) bad(`fragments[${k}] has atoms [${extra.join(", ")}] that its parent fragments[${p}] does not contain`);
        }
      }
    }
  });

  // --- terms ------------------------------------------------------------
  const terms = list(payload.terms);
  const nTerms = terms.length;
  const keys = new Map<string, number>();
  let interactions = 0;
  terms.forEach((t, k) => {
    if (!isObject(t)) return;
    const frs = t.fragments;
    let scope: string | null = null;
    if (Array.isArray(frs)) {
      for (const f of frs) {
        if (isInt(f) && (f < 0 || f >= nFrag)) bad(`terms[${k}].fragments names fragment ${f}, but there are ${nFrag} fragments`);
      }
      if (t.directed === true && frs.length !== 2) {
        bad(`terms[${k}] is directed but names ${frs.length} fragments; a directed term names [donor, acceptor]`);
      }
      const ints = frs.filter(isInt);
      scope = JSON.stringify(t.directed === true ? ints : [...ints].sort((x, y) => x - y));
    }
    if (has(t, "parent")) {
      const p = t.parent;
      if (!isInt(p) || p < 0 || p >= nTerms || p === k) bad(`terms[${k}].parent ${show(p)} is not the index of another term`);
    }
    const key = JSON.stringify([t.name ?? null, scope, t.spin ?? null, has(t, "parent") ? t.parent : null]);
    const prev = keys.get(key);
    if (prev !== undefined) bad(`terms[${k}] repeats terms[${prev}] (same name, fragments, spin and parent)`);
    else keys.set(key, k);
    if (t.name === "interaction" && !has(t, "fragments") && !has(t, "spin") && !has(t, "parent")) interactions += 1;
  });
  terms.forEach((t, k) => {
    const seen = new Set<number>([k]);
    let node: unknown = isObject(t) ? t.parent : undefined;
    while (isInt(node) && node >= 0 && node < nTerms) {
      if (seen.has(node)) {
        bad(`terms[${k}] is its own ancestor through parent`);
        break;
      }
      seen.add(node);
      const parent = terms[node];
      node = isObject(parent) ? parent.parent : undefined;
    }
  });
  if (Array.isArray(payload.terms) && interactions !== 1) {
    bad(`terms has ${interactions} interaction terms (name 'interaction' without fragments, spin or parent); exactly one is required`);
  }

  // --- units ------------------------------------------------------------
  for (const [where, energy] of energies(payload)) {
    if (isObject(energy) && typeof energy.units === "string" && !ENERGY_UNITS.has(energy.units)) {
      errors.push([
        "E_UNIT_CONTRADICTED",
        `${where} is in '${energy.units}', which is not an energy unit (hartree, ev, kcal_per_mol, kj_per_mol) (${SECTION})`,
      ]);
    }
  }

  // --- sum rules --------------------------------------------------------
  list(payload.sum_rules).forEach((r, k) => {
    if (!isObject(r)) return;
    const total = r.total;
    const parts = r.parts;
    const idx = [...(isInt(total) ? [total] : []), ...list(parts).filter(isInt)];
    const outOfRange = idx.filter((i) => i < 0 || i >= nTerms);
    if (outOfRange.length > 0) {
      bad(`sum_rules[${k}] names term(s) [${outOfRange.join(", ")}], but there are ${nTerms} terms`);
      return;
    }
    if (!isInt(total) || !Array.isArray(parts)) return;
    if (parts.includes(total)) {
      bad(`sum_rules[${k}] lists its total terms[${total}] among its parts`);
      return;
    }
    const values = [total, ...(parts as number[])].map((i) => {
      const term = terms[i];
      return isObject(term) ? energyEh(term.energy) : null;
    });
    const tol = has(r, "tolerance") ? energyEh(r.tolerance) : DEFAULT_SUM_TOLERANCE_EH;
    if (tol === null || values.some((v) => v === null)) return;
    const nums = values as number[];
    const sum = nums.slice(1).reduce((a, b) => a + b, 0);
    const residual = (nums[0] as number) - sum;
    if (Math.abs(residual) > tol) {
      errors.push([
        "E_SUM_RULE",
        `sum_rules[${k}]: terms[${total}] is ${nums[0]} hartree but its parts [${parts.join(", ")}] add up to ${sum} hartree (residual ${residual.toPrecision(3)}, tolerance ${tol}) (${SECTION})`,
      ]);
    }
  });

  // --- orbital pairs ----------------------------------------------------
  let nK: number | null = null;
  const periodic = payload.periodic;
  if (isObject(periodic) && isObject(periodic.kpoints) && Array.isArray(periodic.kpoints.kpoints)) {
    nK = periodic.kpoints.kpoints.length;
    const w = periodic.kpoints.weights;
    if (Array.isArray(w) && w.length !== nK) bad(`periodic.kpoints has ${nK} k-points but ${w.length} weights`);
  }
  const pairs = payload.orbital_pairs;
  if (isObject(pairs)) {
    if (pairs.type === "nocv" && !has(pairs, "eigenvalue_convention")) {
      bad("orbital_pairs of type 'nocv' must state eigenvalue_convention");
    }
    const seen = new Map<string, number>();
    list(pairs.pairs).forEach((p, k) => {
      if (!isObject(p)) return;
      const spin = typeof p.spin === "string" ? p.spin : "restricted";
      const key = `${JSON.stringify(p.index)}/${spin}`;
      const prev = seen.get(key);
      if (prev !== undefined) bad(`orbital_pairs.pairs[${k}] repeats index ${show(p.index)} of spin '${spin}' (pairs[${prev}])`);
      else seen.set(key, k);
      for (const field of ["donor", "acceptor"]) {
        const f = p[field];
        if (isInt(f) && (f < 0 || f >= nFrag)) bad(`orbital_pairs.pairs[${k}].${field} names fragment ${f}, but there are ${nFrag} fragments`);
      }
      for (const field of ["kpoint_eigenvalues", "kpoint_energies"]) {
        if (!has(p, field)) continue;
        const arr = p[field];
        if (nK === null) bad(`orbital_pairs.pairs[${k}].${field} needs periodic.kpoints`);
        else if (Array.isArray(arr) && arr.length !== nK) bad(`orbital_pairs.pairs[${k}].${field} has ${arr.length} entries for ${nK} k-points`);
      }
    });
  }
  return errors;
}
