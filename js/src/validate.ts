/**
 * The validator (spec § 6, Appendix B).
 *
 * `validate` treats the archive as untrusted input and never throws for what
 * is wrong with it: every problem is a finding with a stable code. It performs
 * the semantic checks of the Python and C++ reference validators, in the same
 * order and with the same codes. The JSON Schema of `manifest.json`
 * (`E_SCHEMA`) and of every JSON member with a payload definition
 * (`E_PAYLOAD_SCHEMA`) is checked by the checker generated from the schema
 * (`schemaCheck.ts`), with no JSON Schema engine; `schemaCheck: false`
 * disables it, and the validator then reports `W_SCHEMA_SKIPPED` once.
 */

import { declaredSizeProblem, expectedByteLength, isDType, isShape } from "./binary.js";
import { Findings, type Finding } from "./findings.js";
import { decodeUtf8, jsonDepth, isValidUtf8, parseJsonBytes } from "./json.js";
import { CANONICAL_KINDS } from "./kinds.js";
import { resolveLimits, type Limits } from "./limits.js";
import { DTYPE_ITEMSIZE } from "./model.js";
import { sha256Hex } from "./sha256.js";
import { BOHR_PER_ANGSTROM, canonicalUnit, isUnitName } from "./units.js";
import { pathProblem, ratioExceeded, ZipArchive, type ZipEntry } from "./zip.js";
import { checkDefinition, checkManifest, jsonEqual, payloadDefinition } from "./schemaCheck.js";

/**
 * How much of spec § 6 to check. Levels are cumulative; each runs everything
 * the previous ones do:
 *
 * - `container`: ZIP readability, unique and legal entry names (several
 *   `manifest.json` entries are an archive updated in place: the last is the
 *   manifest, each earlier one a `W_ORPHAN_MEMBER`, spec § 2.3), size and
 *   compression-ratio caps (central directory only; nothing is decompressed).
 * - `manifest`: `manifest.json` (caps, JSON, version, source, the JSON Schema),
 *   sections (ids, kinds, critical governance), member specs (shape, path,
 *   existence, a declared `size` against the central directory), orphan
 *   entries. The withdrawn `qvf_version` 2 is checked as 1, the JSON Schema
 *   check included, and reported as `W_VERSION_WITHDRAWN` (spec § 3.1).
 * - `payloads`: reads every member: sha256, JSON caps and syntax, the payload
 *   definition of a JSON member, dtype and byte length, UTF-8 of `run.record`
 *   input/log, the bytes of a `preview` member against its `media_type`.
 * - `semantics`: kind-specific payload rules (`structure` periodicity, the
 *   `bands` spin axis, the per-frame members, `coords` precision, binary
 *   `energies` and `lattice` / `pbc` of `trajectory` and `reaction.path`,
 *   `forces` and `hessian` shapes, the `meta` entry of every
 *   `atom_properties` custom role, the state indices and `n_states` of
 *   `excited_states`, the `wavefunction.gto` rules of spec § 5.3 and the
 *   Appendix A.1 diagonal of its `overlap`, `E_NORMALIZATION`, the data shapes
 *   of the vector, complex, series and Laplacian volumes and a series `axis`
 *   against its slabs, a grid `units.data` against the unit § 4.3 fixes for
 *   its volume kind, `E_UNIT_CONTRADICTED`, `wavefunction.bloch` member shapes and k-point
 *   weights, the property kinds' symmetry and sum rules and the
 *   agreement of a `properties.electric` dipole with the root
 *   `dipole_moment`, the shared array length of a 1-D spectrum, the step ids
 *   and repeated participants of a `reaction.scheme`, a `structure` section's
 *   `charge` and `multiplicity` against its atoms, `E_CHARGE_MULTIPLICITY`), unit names (`W_UNIT_UNKNOWN`) and the agreement of the
 *   root `energies` block with `provenance` and `thermochemistry`, and of the
 *   only `structure` section's `charge` and `multiplicity` with `provenance`.
 * - `references`: cross-references (`volume.difference` operands,
 *   `reaction.waypoints` `trajectory_ref`, any section's `structure_ref` and
 *   `run_ref`, `job.spec` payload `structure_ref`, `run.record` `job_ref` and
 *   `files` index, the participants of a `reaction.scheme`,
 *   `vibrations` `hessian_ref`, the structure of `forces` and `hessian`
 *   sections and their atom counts, the references and orbital rows of
 *   `excited_states`, a spectrum's `states_ref`, a transition density's
 *   `excited_states_ref`, the `wavefunction_ref`, `mo_index` and `spin` of
 *   a density, spin density or orbital volume, a `volume.series` axis `trajectory_ref`, a grid
 *   `structure_ref`, and the span of a periodic grid against the cell of its
 *   structure, the Born charge count of `properties.dielectric`, and the
 *   periodic structure of a `wavefunction.bloch` section).
 * - `all`: everything (the default; currently the same as `references`).
 */
export type ValidationLevel = "container" | "manifest" | "payloads" | "semantics" | "references" | "all";

const LEVELS: readonly ValidationLevel[] = ["container", "manifest", "payloads", "semantics", "references", "all"];

const AT_MANIFEST = LEVELS.indexOf("manifest");
const AT_PAYLOADS = LEVELS.indexOf("payloads");
const AT_SEMANTICS = LEVELS.indexOf("semantics");
const AT_REFERENCES = LEVELS.indexOf("references");

export interface ValidateOptions {
  /** How much to check; `"all"` by default. */
  level?: ValidationLevel;
  /** Overrides of the resource caps of spec § 6.1. */
  limits?: Partial<Limits>;
  /**
   * Run the JSON Schema check generated from the schema (`E_SCHEMA`,
   * `E_PAYLOAD_SCHEMA`); `true` by default. `false` skips it and reports
   * `W_SCHEMA_SKIPPED`, the behaviour of a validator without a schema engine.
   */
  schemaCheck?: boolean;
}

export interface ValidationReport {
  /** `true` when there is no finding of severity `error`. */
  ok: boolean;
  findings: Finding[];
  /** Number of entries in `manifest.sections` (0 when the manifest is unusable). */
  nSections: number;
}

/** Keys of `viewer_defaults` that are not section ids (spec § 3.1, § 3.4). */
const RESERVED_SECTION_IDS: ReadonlySet<string> = new Set(["auto_open", "bookmarks", "camera", "crossfade"]);

const CANONICAL_KIND_SET: ReadonlySet<string> = new Set(CANONICAL_KINDS);

/** Generic references any section may carry, and the kind each must name (spec § 4.1). */
const SECTION_REFS: readonly (readonly [string, string])[] = [
  ["structure_ref", "structure"],
  ["run_ref", "run.record"],
];

type Obj = Record<string, unknown>;

function isObject(value: unknown): value is Obj {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

const has = (obj: Obj, key: string): boolean => Object.hasOwn(obj, key);

/** Python-style truthiness, so loosely typed flags behave as in the reference validators. */
function truthy(value: unknown): boolean {
  if (value === null || value === undefined || value === false || value === 0 || value === "") return false;
  if (Array.isArray(value)) return value.length > 0;
  if (isObject(value)) return Object.keys(value).length > 0;
  return true;
}

/** A compact rendering of a JSON value for messages. */
function show(value: unknown): string {
  if (value === undefined) return "nothing";
  let text: string;
  try {
    text = JSON.stringify(value) ?? String(value);
  } catch {
    return "<value>";
  }
  return text.length > 80 ? `${text.slice(0, 77)}...` : text;
}

const quote = (s: string): string => `'${s}'`;

function escapePointer(part: string): string {
  return part.replaceAll("~", "~0").replaceAll("/", "~1");
}

const sectionPointer = (i: number): string => `/sections/${i}`;
const memberPointer = (i: number, role: string): string => `${sectionPointer(i)}/members/${escapePointer(role)}`;

function errorText(err: unknown): string {
  return err instanceof Error ? err.message : String(err);
}

/**
 * Semantic checks on a `structure` payload (spec § 5.1): `pbc` is the
 * normative carrier of periodicity and `dimensionality` is derived from it.
 */
/** A member spec's `shape` when it is a list of non-negative integers, else null. */
function shapeOf(member: unknown): number[] | null {
  if (!isObject(member)) return null;
  const shape = member.shape;
  if (!Array.isArray(shape) || !shape.every((d) => Number.isInteger(d) && (d as number) >= 0)) return null;
  return shape as number[];
}

const sameExtents = (a: number[], b: number[]): boolean => a.length === b.length && a.every((v, i) => v === b[i]);

/**
 * Spin-axis agreement of a `bands` section's binary members (spec § 5.6), from the
 * manifest alone: a rank-4 `projections` and an `occupations` repeat the
 * `[n_spin, n_k, n_bands]` extents of `eigenvalues`; a legacy rank-3
 * `projections` repeats `[n_k, n_bands]`. Other ranks are left to the schema.
 */
function bandsShapeErrors(sid: string, members: Obj): Array<[string, string]> {
  const eig = shapeOf(members.eigenvalues);
  if (eig === null || eig.length !== 3) return [];
  const errors: Array<[string, string]> = [];
  const proj = shapeOf(members.projections);
  if (proj !== null && proj.length === 4 && !sameExtents(proj.slice(0, 3), eig)) {
    errors.push(["projections", `section ${quote(sid)}: projections shape [${proj}] has [n_spin, n_k, n_bands] = [${proj.slice(0, 3)}], but eigenvalues is [${eig}] (spec § 5.6)`]);
  } else if (proj !== null && proj.length === 3 && !sameExtents(proj.slice(0, 2), eig.slice(1))) {
    errors.push(["projections", `section ${quote(sid)}: rank-3 projections shape [${proj}] has [n_k, n_bands] = [${proj.slice(0, 2)}], but eigenvalues has [${eig.slice(1)}] (spec § 5.6)`]);
  }
  const occ = shapeOf(members.occupations);
  if (occ !== null && occ.length === 3 && !sameExtents(occ, eig)) {
    errors.push(["occupations", `section ${quote(sid)}: occupations shape [${occ}] differs from eigenvalues [${eig}]; both are [n_spin, n_k, n_bands] (spec § 5.6)`]);
  }
  return errors;
}

/** The kinds whose sections are frames of one set of atoms (spec § 5.5). */
const TRAJECTORY_KINDS: ReadonlySet<string> = new Set(["trajectory", "reaction.path"]);

/** The `[n_frames]` members of each trajectory kind (spec § 5.5). */
const TRAJECTORY_FRAME_SCALARS: Readonly<Record<string, readonly string[]>> = {
  trajectory: ["times", "energies", "cell_volumes", "temperatures"],
  "reaction.path": ["times", "energies"],
};

/**
 * Largest allowed difference, in hartree, between the binary `energies` member
 * of a trajectory or reaction path and its `metadata.energies` (spec § 5.5).
 */
const TRAJECTORY_ENERGY_TOLERANCE_EH = 1e-10;

/** The values of a rank-1 float64 member whose byte length matches its shape, little-endian; else null. */
function float64Values(bytes: Uint8Array | undefined, member: unknown): number[] | null {
  const shape = shapeOf(member);
  if (bytes === undefined || !isObject(member) || member.dtype !== "float64" || shape === null) return null;
  const [n] = shape;
  if (shape.length !== 1 || n === undefined || bytes.length !== 8 * n) return null;
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  return Array.from({ length: n }, (_, k) => view.getFloat64(8 * k, true));
}

/**
 * Shape and consistency rules of a `trajectory` or `reaction.path` section
 * (spec § 5.5): `coords` float32 only with `metadata.precision` `"float32"`
 * (and `precision`, when present, naming the dtype); against a
 * `[n_frames, n_atoms, 3]` `coords` (another rank is left to the schema),
 * `velocities` and `gradient` `[n_frames, n_atoms, 3]`, `times` and
 * `energies` (and, for a trajectory, `cell_volumes` and `temperatures`)
 * `[n_frames]`, `lattice` `[3, 3]` or `[n_frames, 3, 3]` (each kind's lattice
 * read in its own convention, spec § 4.3, but checked here by shape only); the
 * binary `energies` equal to `metadata.energies` within 1e-10 hartree; a
 * trajectory `lattice` needs `metadata.pbc`, and a true `pbc` entry needs the
 * `lattice`. Rules that need the metadata are skipped when it could not be
 * read. Returns `[role, message]` pairs.
 */
function trajectoryErrors(
  sid: string,
  kind: string,
  members: Obj,
  metadata: unknown,
  energiesBytes: Uint8Array | undefined,
): Array<[string, string]> {
  const errors: Array<[string, string]> = [];
  const sec = `section ${quote(sid)}: `;
  const meta = isObject(metadata) ? metadata : null;
  const coordsSpec = members.coords;
  const dtype = isObject(coordsSpec) ? coordsSpec.dtype : undefined;
  if (meta !== null && (dtype === "float32" || dtype === "float64")) {
    const f32 = meta.precision === "float32";
    if (dtype === "float32" && !f32) {
      errors.push([
        "coords",
        `${sec}coords are float32 but metadata.precision is ${has(meta, "precision") ? show(meta.precision) : "absent"}; float32 coords require precision 'float32' (spec § 5.5)`,
      ]);
    } else if (dtype === "float64" && f32) {
      errors.push([
        "coords",
        `${sec}metadata.precision is 'float32' but coords are float64; precision must name the dtype of coords (spec § 5.5)`,
      ]);
    }
  }
  const coords = shapeOf(coordsSpec);
  if (coords !== null && coords.length === 3 && coords[2] === 3) {
    const [nFrames, nAtoms] = coords as [number, number, number];
    const vectors = [nFrames, nAtoms, 3];
    for (const role of ["velocities", "gradient"]) {
      const shape = shapeOf(members[role]);
      if (shape !== null && !sameExtents(shape, vectors)) {
        errors.push([role, `${sec}${role} shape [${shape}] is not [n_frames, n_atoms, 3] = [${vectors}], the extents of coords (spec § 5.5)`]);
      }
    }
    for (const role of TRAJECTORY_FRAME_SCALARS[kind] ?? []) {
      const shape = shapeOf(members[role]);
      if (shape !== null && !sameExtents(shape, [nFrames])) {
        errors.push([role, `${sec}${role} shape [${shape}] is not [n_frames] = [${nFrames}], the frames of coords (spec § 5.5)`]);
      }
    }
    const lattice = shapeOf(members.lattice);
    if (lattice !== null && !sameExtents(lattice, [3, 3]) && !sameExtents(lattice, [nFrames, 3, 3])) {
      errors.push(["lattice", `${sec}lattice shape [${lattice}] is neither [3, 3] nor [n_frames, 3, 3] = [${nFrames},3,3] (spec § 5.5)`]);
    }
  }
  if (meta !== null) {
    const binary = float64Values(energiesBytes, members.energies);
    const listed = meta.energies;
    if (binary !== null && Array.isArray(listed) && listed.every((x) => typeof x === "number")) {
      const values = listed as number[];
      if (values.length !== binary.length) {
        errors.push([
          "energies",
          `${sec}the energies member has ${binary.length} values but metadata.energies has ${values.length}; they must agree (spec § 5.5)`,
        ]);
      } else {
        const bad = binary.flatMap((b, k) => (Math.abs(b - (values[k] as number)) <= TRAJECTORY_ENERGY_TOLERANCE_EH ? [] : [k]));
        const first = bad[0];
        if (first !== undefined) {
          errors.push([
            "energies",
            `${sec}the energies member disagrees with metadata.energies at ${bad.length} frame(s), first frame ${first}: ${binary[first]} vs ${values[first]} hartree; they must agree within 1e-10 hartree (spec § 5.5)`,
          ]);
        }
      }
    }
  }
  if (kind === "trajectory" && meta !== null) {
    const pbc = meta.pbc;
    const hasLattice = has(members, "lattice");
    if (hasLattice && (pbc === undefined || pbc === null)) {
      errors.push([
        "metadata",
        `${sec}the section has a lattice member but its metadata has no pbc; pbc is REQUIRED when lattice is present (spec § 5.5)`,
      ]);
    } else if (Array.isArray(pbc) && pbc.some((x) => x === true) && !hasLattice) {
      errors.push(["metadata", `${sec}metadata.pbc ${show(pbc)} marks a periodic axis but the section has no lattice member (spec § 5.5)`]);
    }
  }
  return errors;
}

/** A one-dimensional spectrum kind: its `spectrum` member is a `$defs/SpectrumPayload` (spec § 5.4). */
const isSpectrum1d = (kind: string): boolean => payloadDefinition(kind, "spectrum") === "SpectrumPayload";

/**
 * The equal-length rule of a 1-D `spectrum` payload (spec § 5.4):
 * `frequencies`, `intensities` and, when present, `x` have one shared length.
 * JSON Schema cannot say so, so this runs with or without the schema check.
 */
function spectrumPayloadErrors(sid: string, payload: unknown): string[] {
  if (!isObject(payload) || !Array.isArray(payload.frequencies)) return [];
  const n = payload.frequencies.length;
  const errors: string[] = [];
  for (const key of ["intensities", "x"]) {
    const other = payload[key];
    if (Array.isArray(other) && other.length !== n) {
      errors.push(
        `section ${quote(sid)}: spectrum ${key} has ${other.length} entries but frequencies has ${n}; they must have the same length (spec § 5.4)`,
      );
    }
  }
  return errors;
}

/**
 * Every declared unit of a member payload (spec § 4.3, "Declared units") as
 * `[location inside the payload, value]`: each value of its `units` object,
 * for a one-dimensional spectrum `x_units`, `y_units` and
 * `broadening.units`, and for `excited_states` each value of `meta.units`.
 */
function declaredUnits(kind: string, role: string, payload: unknown): Array<[string, unknown]> {
  if (!isObject(payload)) return [];
  const out: Array<[string, unknown]> = [];
  if (isObject(payload.units)) {
    for (const [key, value] of Object.entries(payload.units)) out.push([`units/${key}`, value]);
  } else if (has(payload, "units") && kind === "volume.series" && role === "axis") {
    out.push(["units", payload.units]);
  }
  if (kind === "excited_states" && isObject(payload.meta) && isObject(payload.meta.units)) {
    for (const [key, value] of Object.entries(payload.meta.units)) out.push([`meta/units/${key}`, value]);
  }
  if (isSpectrum1d(kind)) {
    for (const key of ["x_units", "y_units"]) if (has(payload, key)) out.push([key, payload[key]]);
    const broadening = payload.broadening;
    if (isObject(broadening) && has(broadening, "units")) out.push(["broadening/units", broadening.units]);
  }
  return out;
}

/** Hartree per unit for the energy units the root check converts (spec § 3.6, CODATA 2018). */
const HARTREE_PER_UNIT: ReadonlyMap<string, number> = new Map([
  ["hartree", 1],
  ["eh", 1],
  ["ha", 1],
  ["ev", 1 / 27.211386245988],
  ["kcal_per_mol", 1 / 627.5094740631],
  ["kj_per_mol", 1 / 2625.4996394799],
]);

/** Largest allowed difference, in hartree, between two root statements of one energy. */
const ENERGY_TOLERANCE_EH = 1e-8;

const asciiLower = (s: string): string => s.replace(/[A-Z]/g, (c) => c.toLowerCase());

const isNumber = (v: unknown): v is number => typeof v === "number" && Number.isFinite(v);

/** A `{value, units}` energy in hartree, or null when malformed or in an unrecognised unit. */
function energyEh(energy: unknown): number | null {
  if (!isObject(energy)) return null;
  const { value, units } = energy;
  if (!isNumber(value) || typeof units !== "string") return null;
  const factor = HARTREE_PER_UNIT.get(asciiLower(units));
  return factor === undefined ? null : value * factor;
}

/**
 * Root blocks that state one energy twice must agree (spec § 3.6): every
 * `energies` total stated for the method `provenance.method` names against
 * `provenance.scf_energy`, and `energies.components.zpve` against
 * `thermochemistry.zpve_eh`. Returns `[where, message]` pairs.
 */
function rootConsistencyErrors(manifest: Obj): Array<[string, string]> {
  const energies = manifest.energies;
  if (!isObject(energies)) return [];
  const errors: Array<[string, string]> = [];
  const compare = (where: string, label: string, got: number | null, ref: number, refLabel: string): void => {
    if (got !== null && Math.abs(got - ref) > ENERGY_TOLERANCE_EH) {
      errors.push([
        where,
        `${label} is ${got.toFixed(10)} hartree but ${refLabel} is ${ref.toFixed(10)} hartree; they state one energy and must agree within 1e-8 hartree (spec § 3.6)`,
      ]);
    }
  };
  const prov = isObject(manifest.provenance) ? manifest.provenance : null;
  const method = prov?.method;
  const scf = prov ? energyEh(prov.scf_energy) : null;
  if (typeof method === "string" && scf !== null) {
    const key = asciiLower(method);
    const refLabel = `provenance.scf_energy (method ${quote(method)})`;
    if (typeof energies.method === "string" && asciiLower(energies.method) === key) {
      compare("/energies/total", "energies.total", energyEh(energies.total), scf, refLabel);
    }
    const byMethod = Array.isArray(energies.by_method) ? energies.by_method : [];
    byMethod.forEach((entry, k) => {
      if (isObject(entry) && typeof entry.method === "string" && asciiLower(entry.method) === key) {
        compare(`/energies/by_method/${k}/total`, `energies.by_method[${k}].total`, energyEh(entry.total), scf, refLabel);
      }
    });
  }
  const thermo = manifest.thermochemistry;
  const components = energies.components;
  if (isObject(thermo) && isNumber(thermo.zpve_eh) && isObject(components)) {
    compare("/energies/components/zpve", "energies.components.zpve", energyEh(components.zpve), thermo.zpve_eh, "thermochemistry.zpve_eh");
  }
  return errors;
}

/** What the forces / hessian checks know about the archive's structures. */
interface StructureIndex {
  kindsById: ReadonlyMap<string, string>;
  nAtomsOf: ReadonlyMap<string, number>;
  onlyStructure: string | null;
}

/**
 * Shape and reference rules of a `forces` or `hessian` section (spec § 5.10,
 * § 5.11): `gradient` `[n_atoms, 3]`, `stress` `[3, 3]`, `matrix` square with an
 * extent 3N divisible by 3, `meta.atom_masses_amu` with N entries. With
 * `structures` (the references level), `meta.structure_ref` must name a
 * `structure` section, and n_atoms / N must equal the atom count of the
 * section's structure (the section-level `structure_ref`, else
 * `meta.structure_ref`, else the only `structure` section) when that is known.
 * Returns `[code, role, message]` triples.
 */
function derivativeErrors(
  sid: string,
  kind: string,
  section: Obj,
  members: Obj,
  metaDoc: unknown,
  structures: StructureIndex | null,
): Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> {
  const errors: Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> = [];
  const sec = `section ${quote(sid)}: `;
  const meta: Obj = isObject(metaDoc) ? metaDoc : {};
  let structure: string | null = null;
  if (structures) {
    const metaRef = meta.structure_ref;
    if (has(meta, "structure_ref") && (typeof metaRef !== "string" || structures.kindsById.get(metaRef) !== "structure")) {
      errors.push(["E_REF_UNRESOLVED", "meta", `${sec}meta.structure_ref ${show(metaRef)} does not name a structure section (spec § 5.10)`]);
    }
    if (has(section, "structure_ref")) {
      structure = typeof section.structure_ref === "string" ? section.structure_ref : null;
    } else if (typeof metaRef === "string") {
      structure = metaRef;
    } else {
      structure = structures.onlyStructure;
    }
  }
  const nStruct = structure !== null && structures ? (structures.nAtomsOf.get(structure) ?? null) : null;

  if (kind === "forces") {
    const grad = shapeOf(members.gradient);
    if (grad !== null && (grad.length !== 2 || grad[1] !== 3)) {
      errors.push(["E_MEMBER_MALFORMED", "gradient", `${sec}gradient shape [${grad}] is not [n_atoms, 3] (spec § 5.10)`]);
    } else if (grad !== null && nStruct !== null && grad[0] !== nStruct) {
      errors.push([
        "E_MEMBER_MALFORMED",
        "gradient",
        `${sec}gradient has ${grad[0]} rows but structure ${quote(structure ?? "")} has ${nStruct} atoms (spec § 5.10)`,
      ]);
    }
    const stress = shapeOf(members.stress);
    if (stress !== null && !sameExtents(stress, [3, 3])) {
      errors.push(["E_MEMBER_MALFORMED", "stress", `${sec}stress shape [${stress}] is not [3, 3] (spec § 5.10)`]);
    }
    return errors;
  }

  const mat = shapeOf(members.matrix);
  if (mat === null) return errors;
  if (mat.length !== 2 || mat[0] !== mat[1]) {
    errors.push(["E_MEMBER_MALFORMED", "matrix", `${sec}matrix shape [${mat}] is not square (spec § 5.11)`]);
    return errors;
  }
  const extent = mat[0] as number;
  if (extent % 3 !== 0) {
    errors.push(["E_MEMBER_MALFORMED", "matrix", `${sec}matrix extent ${extent} is not divisible by 3, so the matrix is not [3N, 3N] (spec § 5.11)`]);
    return errors;
  }
  const n = extent / 3;
  if (nStruct !== null && n !== nStruct) {
    errors.push([
      "E_MEMBER_MALFORMED",
      "matrix",
      `${sec}matrix [${mat}] is for ${n} atoms but structure ${quote(structure ?? "")} has ${nStruct} atoms (spec § 5.11)`,
    ]);
  }
  const masses = meta.atom_masses_amu;
  if (Array.isArray(masses) && masses.length !== n) {
    errors.push(["E_MEMBER_MALFORMED", "meta", `${sec}meta.atom_masses_amu has ${masses.length} entries but the matrix is for ${n} atoms (spec § 5.11)`]);
  }
  return errors;
}

/**
 * Atom indices of a standalone `bonds` section against the structure its
 * `structure_ref` names (spec § 5.1). Only a section whose `structure_ref`
 * names a `structure` section with a known atom count is checked; a reference
 * that does not resolve is the generic reference check's finding, and a
 * section without one is not checked. Every integer `i` and `j` of the
 * payload's `pairs` must be an atom of that structure. Returns messages,
 * located at the `bonds` member.
 */
function bondsErrors(sid: string, section: Obj, payload: unknown, structures: StructureIndex): string[] {
  const ref = section.structure_ref;
  if (typeof ref !== "string" || structures.kindsById.get(ref) !== "structure") return [];
  const n = structures.nAtomsOf.get(ref);
  const pairs = isObject(payload) ? payload.pairs : undefined;
  if (n === undefined || !Array.isArray(pairs)) return [];
  const errors: string[] = [];
  pairs.forEach((pair, k) => {
    if (!isObject(pair)) return;
    for (const key of ["i", "j"]) {
      const atom = pair[key];
      if (isInt(atom) && !(atom >= 0 && atom < n)) {
        errors.push(
          `section ${quote(sid)}: pairs[${k}].${key} ${atom} is not one of the ${n} atoms of structure ${quote(ref)}, which structure_ref names (spec § 5.1)`,
        );
      }
    }
  });
  return errors;
}

/** Kinds a state's `density_ref` may name (spec § 5.12). */
const STATE_DENSITY_KINDS: ReadonlySet<string> = new Set(["volume.transition_density", "volume.difference"]);

const isInt = (v: unknown): v is number => typeof v === "number" && Number.isInteger(v);

/**
 * Orbital rows of a `wavefunction.gto` section's coefficient members, from
 * their manifest shapes: key `""` for `mo_coefficients`, `"alpha"` and
 * `"beta"` for the unrestricted pair.
 */
function moRows(members: unknown): Map<string, number> {
  const rows = new Map<string, number>();
  if (!isObject(members)) return rows;
  for (const [role, key] of [
    ["mo_coefficients", ""],
    ["mo_coefficients_alpha", "alpha"],
    ["mo_coefficients_beta", "beta"],
  ] as const) {
    const shape = shapeOf(members[role]);
    if (shape !== null && shape.length > 0) rows.set(key, shape[0] as number);
  }
  return rows;
}

/**
 * The rows a contribution with `spin` may index: those of `mo_coefficients`,
 * else of the channel `spin` selects, else (no `spin`) of both channels; null
 * when unknown.
 */
function contributionRows(rows: Map<string, number>, spin: unknown): number | null {
  const single = rows.get("");
  if (single !== undefined) return single;
  if (spin === "alpha" || spin === "beta") return rows.get(spin) ?? null;
  const known = ["alpha", "beta"].map((c) => rows.get(c)).filter((n): n is number => n !== undefined);
  return known.length > 0 ? Math.min(...known) : null;
}

/** The `index` values of an `excited_states` payload's states, or null without a `states` array. */
function stateIndices(payload: unknown): Set<number> | null {
  if (!isObject(payload) || !Array.isArray(payload.states)) return null;
  const out = new Set<number>();
  for (const s of payload.states) if (isObject(s) && isInt(s.index)) out.add(s.index);
  return out;
}

/** What the excited-state reference checks know about the archive. */
interface ExcitedStatesIndex {
  kindsById: ReadonlyMap<string, string>;
  moRowsOf: ReadonlyMap<string, Map<string, number>>;
}

/**
 * Rules of an `excited_states` section that JSON Schema cannot express (spec
 * § 5.12): the `index` values of `states` are unique and ascending, and
 * `meta.n_states` equals the number of states. With `refs` (the references
 * level), the section's `wavefunction_ref` and every state's `nto_ref` name
 * `wavefunction.gto` sections, every `density_ref` a
 * `volume.transition_density` or `volume.difference` section, and, when
 * `wavefunction_ref` names a `wavefunction.gto` section, every contribution's
 * `from_mo` and `to_mo` is one of its orbital rows. Returns `[code, field,
 * message]` triples; `field` is the section field a finding is about, or `""`
 * for the `states` member.
 */
function excitedStatesErrors(
  sid: string,
  section: Obj,
  payload: unknown,
  refs: ExcitedStatesIndex | null,
): Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> {
  const errors: Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> = [];
  const sec = `section ${quote(sid)}: `;
  const kindOf = (ref: unknown): string | undefined =>
    refs && typeof ref === "string" ? refs.kindsById.get(ref) : undefined;
  let wfRef: string | null = null;
  if (refs && has(section, "wavefunction_ref")) {
    const ref = section.wavefunction_ref;
    if (kindOf(ref) !== "wavefunction.gto") {
      errors.push(["E_REF_UNRESOLVED", "wavefunction_ref", `${sec}wavefunction_ref ${show(ref)} does not name a wavefunction.gto section (spec § 5.12)`]);
    } else {
      wfRef = ref as string;
    }
  }
  if (!isObject(payload) || !Array.isArray(payload.states)) return errors;
  const states = payload.states;
  const rows = wfRef !== null && refs ? (refs.moRowsOf.get(wfRef) ?? new Map<string, number>()) : new Map<string, number>();
  const seen = new Set<number>();
  let prev: number | null = null;
  states.forEach((state, k) => {
    if (!isObject(state)) return;
    const idx = state.index;
    if (isInt(idx)) {
      if (seen.has(idx)) {
        errors.push(["E_MEMBER_MALFORMED", "", `${sec}states[${k}].index ${idx} repeats an earlier index; state indices are unique (spec § 5.12)`]);
      } else if (prev !== null && idx < prev) {
        errors.push(["E_MEMBER_MALFORMED", "", `${sec}states[${k}].index ${idx} follows index ${prev}; states are listed in ascending order of index (spec § 5.12)`]);
      }
      seen.add(idx);
      prev = idx;
    }
    if (!refs) return;
    if (has(state, "nto_ref") && kindOf(state.nto_ref) !== "wavefunction.gto") {
      errors.push(["E_REF_UNRESOLVED", "", `${sec}states[${k}].nto_ref ${show(state.nto_ref)} does not name a wavefunction.gto section (spec § 5.12)`]);
    }
    if (has(state, "density_ref") && !STATE_DENSITY_KINDS.has(kindOf(state.density_ref) ?? "")) {
      errors.push([
        "E_REF_UNRESOLVED",
        "",
        `${sec}states[${k}].density_ref ${show(state.density_ref)} does not name a volume.transition_density or volume.difference section (spec § 5.12)`,
      ]);
    }
    const contributions = state.contributions;
    if (wfRef === null || !Array.isArray(contributions)) return;
    contributions.forEach((contrib, c) => {
      if (!isObject(contrib)) return;
      const n = contributionRows(rows, contrib.spin);
      if (n === null) return;
      for (const key of ["from_mo", "to_mo"]) {
        const mo = contrib[key];
        if (isInt(mo) && (mo < 0 || mo >= n)) {
          errors.push([
            "E_MEMBER_MALFORMED",
            "",
            `${sec}states[${k}].contributions[${c}].${key} ${mo} is not one of the ${n} orbital rows of wavefunction ${quote(wfRef as string)} (spec § 5.12)`,
          ]);
        }
      }
    });
  });
  const meta = payload.meta;
  if (isObject(meta) && isInt(meta.n_states) && meta.n_states !== states.length) {
    errors.push(["E_MEMBER_MALFORMED", "", `${sec}meta.n_states is ${meta.n_states} but states has ${states.length} entries (spec § 5.12)`]);
  }
  return errors;
}

/** Largest allowed difference between 1 and the sum of a `wavefunction.bloch` section's k-point weights. */
const BLOCH_WEIGHT_TOLERANCE = 1e-6;

/** A JSON number with an integral value. */
const isIntegral = (x: unknown): x is number => typeof x === "number" && Number.isInteger(x);

/** What the `wavefunction.bloch` structure rule knows about the archive's structures. */
interface PeriodicIndex {
  kindsById: ReadonlyMap<string, string>;
  /** structure id -> whether its `pbc` marks at least one periodic axis */
  periodicOf: ReadonlyMap<string, boolean>;
}

/**
 * Shape and structure rules of a `wavefunction.bloch` section that JSON Schema
 * cannot express (spec § 5.3.1): `mo_coefficients` `[n_k, n_spin, n_mo, n_ao, 2]`
 * with n_spin 1 or 2; `energies` and `occupations` with its `[n_k, n_spin, n_mo]`
 * (each other's when the coefficients are malformed); n_k k-points, `labels` and
 * `weights` in the `kpoints` payload, the weights summing to 1 within 1e-6;
 * `mo_metadata` `n_spin` (1 for `restricted`, 2 for `unrestricted`), `n_mo`,
 * `n_ao` and `n_mo_per_k`, and `basis.n_ao`, agreeing with those axes. With
 * `structures` (the references level), the section's structure (the
 * section-level `structure_ref`, else `basis.structure_ref`) must be a
 * `structure` section with at least one periodic axis. Mirrors `_bloch_errors`
 * of the Python reference validator. Returns `[code, location, message]`
 * triples, the location relative to the section (`members/<role>` or
 * `structure_ref`).
 */
function blochErrors(
  sid: string,
  section: Obj,
  members: Obj,
  docs: Map<string, unknown>,
  structures: PeriodicIndex | null,
): Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> {
  const errors: Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> = [];
  const err = (code: "E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", where: string, msg: string): void => {
    errors.push([code, where, `section ${quote(sid)}: ${msg} (spec § 5.3.1)`]);
  };

  let coeff = shapeOf(members.mo_coefficients);
  if (coeff !== null && (coeff.length !== 5 || coeff[4] !== 2 || (coeff[1] !== 1 && coeff[1] !== 2))) {
    err("E_MEMBER_MALFORMED", "members/mo_coefficients", `mo_coefficients shape [${coeff}] is not [n_k, n_spin, n_mo, n_ao, 2] with n_spin 1 or 2`);
    coeff = null;
  }
  let ref: number[] | null = coeff !== null ? coeff.slice(0, 3) : null;
  let refRole = "mo_coefficients";
  for (const role of ["energies", "occupations"]) {
    const shape = shapeOf(members[role]);
    if (shape === null) continue;
    if (shape.length !== 3) {
      err("E_MEMBER_MALFORMED", `members/${role}`, `${role} shape [${shape}] is not [n_k, n_spin, n_mo]`);
    } else if (ref === null) {
      ref = shape;
      refRole = role;
    } else if (!sameExtents(shape, ref)) {
      err("E_MEMBER_MALFORMED", `members/${role}`, `${role} shape [${shape}] differs from the [n_k, n_spin, n_mo] = [${ref}] of ${refRole}`);
    }
  }

  const kp = docs.get("kpoints");
  if (isObject(kp)) {
    const nK = Array.isArray(kp.kpoints) ? kp.kpoints.length : null;
    if (nK !== null && ref !== null && nK !== ref[0]) {
      err("E_MEMBER_MALFORMED", "members/kpoints", `kpoints lists ${nK} k-points but axis 0 of ${refRole} is ${ref[0]}`);
    }
    for (const key of ["weights", "labels"]) {
      const arr = kp[key];
      if (Array.isArray(arr) && nK !== null && arr.length !== nK) {
        err("E_MEMBER_MALFORMED", "members/kpoints", `kpoints.${key} has ${arr.length} entries for ${nK} k-points`);
      }
    }
    const weights = kp.weights;
    if (Array.isArray(weights) && weights.every((w) => typeof w === "number")) {
      const total = (weights as number[]).reduce((a, b) => a + b, 0);
      if (Math.abs(total - 1) > BLOCH_WEIGHT_TOLERANCE) {
        err("E_MEMBER_MALFORMED", "members/kpoints", `kpoints.weights sum to ${total}, not to 1 within 1e-6`);
      }
    }
  }

  const meta = docs.get("mo_metadata");
  if (isObject(meta)) {
    const want = meta.spin === "restricted" ? 1 : meta.spin === "unrestricted" ? 2 : null;
    if (want !== null && isIntegral(meta.n_spin) && meta.n_spin !== want) {
      err("E_MEMBER_MALFORMED", "members/mo_metadata", `mo_metadata.spin ${show(meta.spin)} has n_spin ${want}, not ${meta.n_spin}`);
    }
    const axes: Array<[string, number, number[] | null, string]> = [
      ["n_spin", 1, ref, refRole],
      ["n_mo", 2, ref, refRole],
      ["n_ao", 3, coeff, "mo_coefficients"],
    ];
    for (const [key, axis, dims, name] of axes) {
      const value = meta[key];
      if (isIntegral(value) && dims !== null && value !== dims[axis]) {
        err("E_MEMBER_MALFORMED", "members/mo_metadata", `mo_metadata.${key} is ${value} but axis ${axis} of ${name} is ${dims[axis]}`);
      }
    }
    const perK = meta.n_mo_per_k;
    if (Array.isArray(perK)) {
      if (ref !== null && perK.length !== ref[0]) {
        err("E_MEMBER_MALFORMED", "members/mo_metadata", `mo_metadata.n_mo_per_k has ${perK.length} entries but axis 0 of ${refRole} is ${ref[0]}`);
      }
      const nMo = ref !== null ? ref[2] : meta.n_mo;
      if (isIntegral(nMo) && perK.some((v) => typeof v === "number" && v > nMo)) {
        err("E_MEMBER_MALFORMED", "members/mo_metadata", `mo_metadata.n_mo_per_k has an entry above n_mo = ${nMo}`);
      }
    }
  }

  const basis = docs.get("basis");
  if (isObject(basis) && coeff !== null && isIntegral(basis.n_ao) && basis.n_ao !== coeff[3]) {
    err("E_MEMBER_MALFORMED", "members/basis", `basis.n_ao is ${basis.n_ao} but axis 3 of mo_coefficients is ${coeff[3]}`);
  }

  if (!structures) return errors;
  const namesStructure = (r: unknown): r is string => typeof r === "string" && structures.kindsById.get(r) === "structure";
  let structure: string;
  let where: string;
  if (has(section, "structure_ref")) {
    if (!namesStructure(section.structure_ref)) return errors; // reported by the generic rule (spec § 4.1)
    structure = section.structure_ref;
    where = "structure_ref";
  } else if (isObject(basis)) {
    where = "members/basis";
    if (!namesStructure(basis.structure_ref)) {
      err("E_REF_UNRESOLVED", where, `basis.structure_ref ${show(basis.structure_ref)} does not name a structure section`);
      return errors;
    }
    structure = basis.structure_ref;
  } else {
    return errors; // no basis payload: its own finding says why
  }
  if (structures.periodicOf.get(structure) === false) {
    err("E_MEMBER_MALFORMED", where, `structure ${quote(structure)} has no periodic axis; Bloch orbitals need a periodic structure`);
  }
  return errors;
}

/** The property kinds (spec § 5.13). */
const PROPERTY_KINDS: ReadonlySet<string> = new Set(["properties.electric", "properties.elastic", "properties.dielectric"]);

/**
 * Largest allowed disagreement between two statements of one property (spec
 * § 5.13): matrix symmetry (relative to the largest entry), `static` against
 * `electronic + ionic`, and the electric dipole against the root dipole (debye).
 */
const PROPERTY_TOLERANCE = 1e-6;

/**
 * `value` as a `rows` x `cols` matrix of finite numbers, row-major, or null
 * when it is anything else (the payload schema reports that).
 */
function matrixOf(value: unknown, rows: number, cols: number): number[] | null {
  if (!Array.isArray(value) || value.length !== rows) return null;
  const out: number[] = [];
  for (const row of value) {
    if (!Array.isArray(row) || row.length !== cols) return null;
    for (const x of row) {
      if (!isNumber(x)) return null;
      out.push(x);
    }
  }
  return out;
}

/** The first `[i, j]`, i < j, of an n x n matrix that breaks symmetry, or null. */
function firstAsymmetry(m: number[], n: number): [number, number] | null {
  const scale = m.reduce((s, x) => Math.max(s, Math.abs(x)), 0);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (Math.abs((m[i * n + j] as number) - (m[j * n + i] as number)) > PROPERTY_TOLERANCE * scale) return [i, j];
    }
  }
  return null;
}

/**
 * Rules of the property kinds that JSON Schema cannot express (spec § 5.13):
 * symmetric `properties.elastic` matrices, `static` = `electronic` + `ionic`
 * and one Born tensor per atom of the section's structure (the section-level
 * `structure_ref`, else the only `structure` section) in `properties.dielectric`,
 * and a `properties.electric` dipole that agrees with the root `dipole_moment`
 * within 1e-6 debye. With `structures` null (below the references level) the
 * Born charges are not counted. Returns `[code, message]` pairs.
 */
function propertyErrors(
  sid: string,
  kind: string,
  section: Obj,
  payload: unknown,
  manifest: Obj,
  structures: StructureIndex | null,
): Array<["E_MEMBER_MALFORMED" | "E_ROOT_INCONSISTENT", string]> {
  if (!isObject(payload)) return [];
  const sec = `section ${quote(sid)}: `;
  const errors: Array<["E_MEMBER_MALFORMED" | "E_ROOT_INCONSISTENT", string]> = [];
  if (kind === "properties.elastic") {
    for (const name of ["stiffness", "compliance"]) {
      const m = matrixOf(payload[name], 6, 6);
      const bad = m === null ? null : firstAsymmetry(m, 6);
      if (m !== null && bad !== null) {
        const [i, j] = bad;
        errors.push([
          "E_MEMBER_MALFORMED",
          `${sec}${name} is not symmetric: entry (${i}, ${j}) is ${m[i * 6 + j]} but entry (${j}, ${i}) is ${m[j * 6 + i]}; they must agree within 1e-6 of the largest entry (spec § 5.13)`,
        ]);
      }
    }
  } else if (kind === "properties.dielectric") {
    const e = matrixOf(payload.electronic, 3, 3);
    const ion = matrixOf(payload.ionic, 3, 3);
    const s = matrixOf(payload.static, 3, 3);
    if (e !== null && ion !== null && s !== null) {
      for (let k = 0; k < 9; k++) {
        const sum = (e[k] as number) + (ion[k] as number);
        if (Math.abs((s[k] as number) - sum) > PROPERTY_TOLERANCE) {
          errors.push([
            "E_MEMBER_MALFORMED",
            `${sec}static[${Math.floor(k / 3)}][${k % 3}] is ${s[k]} but electronic + ionic is ${sum}; static must equal their sum within 1e-6 (spec § 5.13)`,
          ]);
          break;
        }
      }
    }
    const born = payload.born_charges;
    if (Array.isArray(born) && structures) {
      const ref = has(section, "structure_ref") ? section.structure_ref : structures.onlyStructure;
      const nStruct = typeof ref === "string" ? structures.nAtomsOf.get(ref) : undefined;
      if (nStruct !== undefined && born.length !== nStruct) {
        errors.push([
          "E_MEMBER_MALFORMED",
          `${sec}born_charges has ${born.length} tensors but structure ${quote(String(ref))} has ${nStruct} atoms; there must be one per atom (spec § 5.13)`,
        ]);
      }
    }
  } else {
    const dipole = payload.dipole;
    const root = manifest.dipole_moment;
    if (isObject(dipole) && isObject(root)) {
      const vec = dipole.vector;
      const rootVec = root.vector_debye;
      if (
        Array.isArray(vec) &&
        Array.isArray(rootVec) &&
        vec.length === 3 &&
        rootVec.length === 3 &&
        vec.every(isNumber) &&
        rootVec.every(isNumber) &&
        vec.some((x, k) => Math.abs((x as number) - (rootVec[k] as number)) > PROPERTY_TOLERANCE)
      ) {
        errors.push([
          "E_ROOT_INCONSISTENT",
          `${sec}dipole.vector ${show(vec)} differs from the root dipole_moment.vector_debye ${show(rootVec)}; the two state one dipole and must agree within 1e-6 debye (spec § 5.13)`,
        ]);
      }
      if (isNumber(dipole.total) && isNumber(root.total_debye) && Math.abs(dipole.total - root.total_debye) > PROPERTY_TOLERANCE) {
        errors.push([
          "E_ROOT_INCONSISTENT",
          `${sec}dipole.total ${show(dipole.total)} differs from the root dipole_moment.total_debye ${show(root.total_debye)}; the two state one dipole and must agree within 1e-6 debye (spec § 5.13)`,
        ]);
      }
    }
  }
  return errors;
}

/**
 * A one-dimensional spectrum's `states_ref` names an `excited_states` section,
 * and every `peaks[].assignment.state_index` is the `index` of one of its
 * states (spec § 5.4). Returns messages.
 */
function spectrumStatesErrors(
  sid: string,
  payload: Obj,
  kindsById: ReadonlyMap<string, string>,
  indicesOf: ReadonlyMap<string, Set<number> | null>,
): string[] {
  if (!has(payload, "states_ref")) return [];
  const ref = payload.states_ref;
  const sec = `section ${quote(sid)}: `;
  if (typeof ref !== "string" || kindsById.get(ref) !== "excited_states") {
    return [`${sec}spectrum states_ref ${show(ref)} does not name an excited_states section (spec § 5.4)`];
  }
  const indices = indicesOf.get(ref) ?? null;
  const peaks = payload.peaks;
  if (indices === null || !Array.isArray(peaks)) return [];
  const errors: string[] = [];
  peaks.forEach((peak, k) => {
    const assignment = isObject(peak) ? peak.assignment : undefined;
    const state = isObject(assignment) ? assignment.state_index : undefined;
    if (isInt(state) && !indices.has(state)) {
      errors.push(`${sec}peaks[${k}].assignment.state_index ${state} names no state of excited_states section ${quote(ref)} (spec § 5.4)`);
    }
  });
  return errors;
}

/**
 * A `volume.transition_density` section's `excited_states_ref` names an
 * `excited_states` section, and its `state_index` the `index` of one of its
 * states (spec § 5.2). Returns `[field, message]` pairs.
 */
function transitionDensityErrors(
  sid: string,
  section: Obj,
  kindsById: ReadonlyMap<string, string>,
  indicesOf: ReadonlyMap<string, Set<number> | null>,
): Array<[string, string]> {
  if (!has(section, "excited_states_ref")) return [];
  const ref = section.excited_states_ref;
  const sec = `section ${quote(sid)}: `;
  if (typeof ref !== "string" || kindsById.get(ref) !== "excited_states") {
    return [["excited_states_ref", `${sec}excited_states_ref ${show(ref)} does not name an excited_states section (spec § 5.2)`]];
  }
  const state = section.state_index;
  const indices = indicesOf.get(ref) ?? null;
  if (isInt(state) && indices !== null && !indices.has(state)) {
    return [["state_index", `${sec}state_index ${state} names no state of excited_states section ${quote(ref)} (spec § 5.2)`]];
  }
  return [];
}

/** Volume kinds that may name the wavefunction they were computed from (spec § 5.2). */
const WAVEFUNCTION_VOLUME_KINDS: ReadonlySet<string> = new Set(["volume.density", "volume.spin", "volume.orbital"]);

/**
 * A `volume.density`, `volume.spin` or `volume.orbital` section's
 * `wavefunction_ref` names a `wavefunction.gto` section. For a
 * `volume.orbital`, `spin` is absent when that section carries
 * `mo_coefficients` and present, whenever `mo_index` is, when it carries the
 * alpha and beta coefficients, and `mo_index` is a row of the member they
 * select, from the manifest shapes (spec § 5.2, Volumes of a wavefunction).
 * Returns `[field, message]` pairs, all `E_REF_UNRESOLVED`.
 */
function wavefunctionVolumeErrors(
  sid: string,
  kind: string,
  section: Obj,
  kindsById: ReadonlyMap<string, string>,
  moRowsOf: ReadonlyMap<string, Map<string, number>>,
): Array<[string, string]> {
  if (!has(section, "wavefunction_ref")) return [];
  const ref = section.wavefunction_ref;
  const sec = `section ${quote(sid)}: `;
  if (typeof ref !== "string" || kindsById.get(ref) !== "wavefunction.gto") {
    return [["wavefunction_ref", `${sec}wavefunction_ref ${show(ref)} does not name a wavefunction.gto section (spec § 5.2)`]];
  }
  if (kind !== "volume.orbital") return [];
  const errors: Array<[string, string]> = [];
  const rows = moRowsOf.get(ref) ?? new Map<string, number>();
  const mo = section.mo_index;
  const spin = section.spin;
  let n: number | undefined;
  const single = rows.get("");
  if (single !== undefined) {
    if (has(section, "spin")) {
      errors.push(["spin", `${sec}spin ${show(spin)} selects a spin channel, but wavefunction ${quote(ref)} carries one set of orbitals in mo_coefficients (spec § 5.2)`]);
    }
    n = single;
  } else if (rows.has("alpha") || rows.has("beta")) {
    if (!has(section, "spin")) {
      if (has(section, "mo_index")) {
        errors.push(["mo_index", `${sec}mo_index ${show(mo)} needs spin: wavefunction ${quote(ref)} carries alpha and beta coefficients (spec § 5.2)`]);
      }
      return errors;
    }
    if (spin === "alpha" || spin === "beta") n = rows.get(spin);
  }
  if (isInt(mo) && n !== undefined && !(mo >= 0 && mo < n)) {
    const member = single !== undefined ? "mo_coefficients" : `mo_coefficients_${String(spin)}`;
    errors.push(["mo_index", `${sec}mo_index ${mo} is not one of the ${n} rows of ${member} of wavefunction ${quote(ref)} (spec § 5.2)`]);
  }
  return errors;
}

/**
 * Custom roles of an `atom_properties` section need a `meta` entry (spec § 5.7):
 * every `custom.<name>` member needs an object entry of that name in the
 * section's `meta` member, and a section without a `meta` member has none. A
 * `meta` member that is declared but could not be read (its own finding says
 * why) is not judged. Returns `[role, message]` pairs.
 */
function atomPropertiesMetaErrors(sid: string, members: Obj, docs: Map<string, unknown>): Array<[string, string]> {
  const custom = Object.keys(members).filter((role) => role.startsWith("custom."));
  if (custom.length === 0) return [];
  const sec = `section ${quote(sid)}: custom role `;
  if (!has(members, "meta")) {
    return custom.map((r): [string, string] => [r, `${sec}${quote(r)} has no meta member to describe it (spec § 5.7)`]);
  }
  if (!docs.has("meta")) return [];
  const meta = docs.get("meta");
  if (!isObject(meta)) {
    return custom.map((r): [string, string] => [r, `${sec}${quote(r)}: the meta member is not a JSON object (spec § 5.7)`]);
  }
  return custom
    .filter((r) => !(has(meta, r) && isObject(meta[r])))
    .map((r): [string, string] => [r, `${sec}${quote(r)} has no entry in the meta member (spec § 5.7)`]);
}

// --- wavefunction.gto (spec § 5.3, Appendix A.1) -------------------------------

/** Tolerances of the electron count and of the Appendix A.1 normalization check. */
const N_ELECTRONS_TOLERANCE = 1e-6;
const NORMALIZATION_RTOL = 1e-6;

/** `%.10g`, as the reference validators print numbers in these messages. */
const g10 = (x: number): string => String(Number(x.toPrecision(10)));

/** One shell of a `basis` payload; `exponents` is null when its primitive arrays are unusable. */
interface ShellInfo {
  l: number;
  pure: boolean;
  exponents: number[] | null;
  coefficients: number[] | null;
}

/**
 * The shells of a `basis` payload, or null when the AO layout cannot be read
 * (no shell list, or a shell whose `l` or `pure` is malformed). Mirrors
 * `_basis_shells()` of the Python reference validator.
 */
function basisShells(basis: unknown): ShellInfo[] | null {
  if (!isObject(basis) || !Array.isArray(basis.shells)) return null;
  const defaultPure = has(basis, "pure") ? basis.pure : true;
  if (typeof defaultPure !== "boolean") return null;
  const out: ShellInfo[] = [];
  for (const sh of basis.shells) {
    if (!isObject(sh)) return null;
    const l = sh.l;
    const pure = has(sh, "pure") ? sh.pure : defaultPure;
    if (typeof l !== "number" || !Number.isInteger(l) || l < 0 || l > 7 || typeof pure !== "boolean") return null;
    const e = sh.exponents;
    const c = sh.coefficients;
    const usable =
      Array.isArray(e) &&
      Array.isArray(c) &&
      e.length === c.length &&
      e.length > 0 &&
      e.every((a) => isNumber(a) && a > 0) &&
      c.every((x) => isNumber(x));
    out.push({
      l,
      pure,
      exponents: usable ? (e as number[]) : null,
      coefficients: usable ? (c as number[]) : null,
    });
  }
  return out;
}

const shellSize = (l: number, pure: boolean): number => (pure ? 2 * l + 1 : ((l + 1) * (l + 2)) / 2);

/** n_ao of a `basis` payload: counted from its shells, else its `n_ao`, else null. */
function basisNAo(basis: unknown): number | null {
  const shells = basisShells(basis);
  if (shells) return shells.reduce((n, s) => n + shellSize(s.l, s.pure), 0);
  const n = isObject(basis) ? basis.n_ao : undefined;
  return typeof n === "number" && Number.isInteger(n) && n >= 0 ? n : null;
}

/**
 * Self-overlap of the axial component (l, 0, 0) of a contracted shell on the
 * N_i-normalized primitives of Appendix A.1:
 * sum_ij c_i c_j (2 sqrt(a_i a_j) / (a_i + a_j))^(l + 3/2).
 */
function axialSelfOverlap(l: number, exponents: number[], coefficients: number[]): number {
  let total = 0;
  exponents.forEach((a, i) => {
    exponents.forEach((b, j) => {
      total += (coefficients[i] as number) * (coefficients[j] as number) * ((2 * Math.sqrt(a * b)) / (a + b)) ** (l + 1.5);
    });
  });
  return total;
}

/** (2n - 1)!!, with (-1)!! = 1. */
function oddDoubleFactorial(n: number): number {
  let out = 1;
  for (let k = 1; k < 2 * n; k += 2) out *= k;
  return out;
}

/** `[name, factor]` of each component of a shell in Appendix A.2 order: its diagonal overlap is factor * S_axial. */
function shellComponents(l: number, pure: boolean): Array<[string, number]> {
  const out: Array<[string, number]> = [];
  if (pure) {
    for (let m = -l; m <= l; m++) out.push([m > 0 ? `m = +${m}` : `m = ${m}`, 1]);
    return out;
  }
  for (let lx = l; lx >= 0; lx--) {
    for (let ly = l - lx; ly >= 0; ly--) {
      const lz = l - lx - ly;
      const name = "x".repeat(lx) + "y".repeat(ly) + "z".repeat(lz);
      out.push([
        `Cartesian ${name || "s"} [${lx}, ${ly}, ${lz}]`,
        (oddDoubleFactorial(lx) * oddDoubleFactorial(ly) * oddDoubleFactorial(lz)) / oddDoubleFactorial(l),
      ]);
    }
  }
  return out;
}

/**
 * The Appendix A.1 check of a `wavefunction.gto` `overlap` member: each diagonal
 * entry agrees with the self-overlap the basis gives for that AO, to a relative
 * tolerance of 1e-6. Runs only on a well-formed float64 [n_ao, n_ao] overlap of
 * a basis whose AO layout can be read; returns at most one message.
 */
function normalizationErrors(sid: string, basis: unknown, overlap: unknown, bytes: Uint8Array | undefined): string[] {
  const shells = basisShells(basis);
  const shape = shapeOf(overlap);
  if (!shells || !shape || !bytes || !isObject(overlap)) return [];
  const nAo = shells.reduce((n, s) => n + shellSize(s.l, s.pure), 0);
  if (!sameExtents(shape, [nAo, nAo]) || overlap.dtype !== "float64" || bytes.length !== 8 * nAo * nAo) return [];
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let nBad = 0;
  let first = "";
  let ao = 0;
  shells.forEach((s, k) => {
    const components = shellComponents(s.l, s.pure);
    if (s.exponents && s.coefficients) {
      const sAxial = axialSelfOverlap(s.l, s.exponents, s.coefficients);
      components.forEach(([name, factor], j) => {
        const i = ao + j;
        const want = sAxial * factor;
        const got = view.getFloat64(8 * (i * nAo + i), true);
        if (!(Math.abs(got - want) <= NORMALIZATION_RTOL * Math.abs(want))) {
          if (nBad++ === 0) {
            first = `AO ${i} (shell ${k}, l = ${s.l}, ${name}) is ${g10(got)} but the basis gives ${g10(want)}`;
          }
        }
      });
    }
    ao += components.length;
  });
  if (!nBad) return [];
  return [
    `section ${quote(sid)}: ${nBad} of ${nAo} overlap diagonal entries disagree with the self-overlaps Appendix A.1 gives for the basis; first: ${first}. The basis.json coefficients are not on N_i-normalized primitives, or the AO order disagrees with Appendix A.2 (spec § 5.3)`,
  ];
}

/** Whether `mo_metadata` declares complex coefficients: `complex`, else `coefficient_encoding`, else real; null when malformed. */
function coefficientsComplex(meta: unknown): boolean | null {
  if (!isObject(meta)) return null;
  if (has(meta, "complex")) return typeof meta.complex === "boolean" ? meta.complex : null;
  if (has(meta, "coefficient_encoding")) {
    if (meta.coefficient_encoding === "complex_split_last_axis") return true;
    return meta.coefficient_encoding === "real" ? false : null;
  }
  return false;
}

type NtoSection = { i: number; meta: unknown; members: Obj };

/** Unrestricted NTO rules (spec § 5.3); normalization is conditional, not imposed. */
function unrestrictedNtoErrors(meta: unknown, members: Obj): string[] {
  if (!isObject(meta) || !has(meta, "nto_role") || meta.orbital_kind !== "nto" || meta.spin !== "unrestricted") return [];
  const errors: string[] = [];
  const rows = moRows(members);
  for (const ch of ["alpha", "beta"]) {
    const block = meta[ch];
    if (!isObject(block)) continue;
    const n = rows.get(ch);
    if (n === undefined) {
      errors.push(`unrestricted NTOs need mo_coefficients_${ch} (spec § 5.3)`);
    } else if (["energies", "occupations"].some((key) => Array.isArray(block[key]) && block[key].length !== n) ||
               (isInt(block.n_mo) && block.n_mo !== n)) {
      errors.push(`NTO ${ch} metadata must have one entry per coefficient row (spec § 5.3)`);
    }
    const occ = block.occupations;
    if (Array.isArray(occ) && occ.every(isNumber) &&
        (occ.some((x) => x < 0) || occ.some((x, i) => i > 0 && occ[i - 1] < x))) {
      errors.push(`NTO ${ch} weights must be nonnegative and non-increasing (spec § 5.3)`);
    }
  }
  return errors;
}

function ntoPairErrors(meta: unknown, members: Obj, ntos: NtoSection[],
                       kinds: Map<string, string>, indices: Map<string, Set<number> | null>): Array<["E_REF_UNRESOLVED" | "E_MEMBER_MALFORMED", string]> {
  if (!isObject(meta) || !has(meta, "nto_role") || meta.orbital_kind !== "nto" || meta.spin !== "unrestricted") return [];
  const ref = meta.state_ref, index = meta.state_index, role = meta.nto_role;
  if (typeof ref !== "string" || !isInt(index) || (role !== "hole" && role !== "particle")) return [];
  const errors: Array<["E_REF_UNRESOLVED" | "E_MEMBER_MALFORMED", string]> = [];
  const states = indices.get(ref);
  if (kinds.get(ref) !== "excited_states" || (states != null && !states.has(index))) {
    errors.push(["E_REF_UNRESOLVED", "NTO state_ref/state_index must name an excited state (spec § 5.3)"]);
  }
  const opposite = role === "hole" ? "particle" : "hole";
  const partners = ntos.filter(({ meta: m }) => isObject(m) && m.orbital_kind === "nto" &&
    m.state_ref === ref && m.state_index === index && m.nto_role === opposite);
  if (partners.length !== 1) {
    errors.push(["E_REF_UNRESOLVED", "NTOs must have one opposite-role section for the same state (spec § 5.3)"]);
    return errors;
  }
  const other = partners[0]!.meta as Obj;
  if (other.spin !== "unrestricted") {
    errors.push(["E_MEMBER_MALFORMED", "NTO partners must have the same spin layout (spec § 5.3)"]);
    return errors;
  }
  const rows = moRows(members), otherRows = moRows(partners[0]!.members);
  for (const ch of ["alpha", "beta"]) {
    if (rows.has(ch) && otherRows.has(ch) && rows.get(ch) !== otherRows.get(ch)) {
      errors.push(["E_MEMBER_MALFORMED", `NTO partners must have the same ${ch} row count (spec § 5.3)`]);
    }
    const a = meta[ch], b = other[ch];
    if (isObject(a) && isObject(b)) {
      const av = a.occupations, bv = b.occupations;
      if (Array.isArray(av) && Array.isArray(bv) && av.every(isNumber) && bv.every(isNumber) &&
          (av.length !== bv.length || av.some((x, i) => x !== bv[i]))) {
        errors.push(["E_MEMBER_MALFORMED", `NTO partners must have identical ${ch} weights (spec § 5.3)`]);
      }
    }
  }
  return errors;
}

type WfCode = "E_MEMBER_MALFORMED" | "E_ROOT_INCONSISTENT";

/**
 * Cross-member rules of a `wavefunction.gto` section (spec § 5.3): coefficient
 * rank against `complex`, `symmetry_labels` lengths, the [n_ao, n_ao] shape and
 * pairing of the AO-basis matrices, `basis.n_electrons` against the occupations
 * (`E_MEMBER_MALFORMED`) and against `provenance.n_electrons`
 * (`E_ROOT_INCONSISTENT`). Returns `[code, role, message]` triples.
 */
function wavefunctionErrors(sid: string, members: Obj, docs: Map<string, unknown>, provenance: unknown): Array<[WfCode, string, string]> {
  const errors: Array<[WfCode, string, string]> = [];
  const sec = `section ${quote(sid)}: `;
  const basis = docs.get("basis");
  const meta = docs.get("mo_metadata");

  const cplx = coefficientsComplex(meta);
  if (cplx !== null) {
    for (const role of ["mo_coefficients", "mo_coefficients_alpha", "mo_coefficients_beta"]) {
      const shape = has(members, role) ? shapeOf(members[role]) : null;
      if (!shape) continue;
      if (cplx && (shape.length !== 3 || shape[2] !== 2)) {
        errors.push(["E_MEMBER_MALFORMED", role, `${sec}${role} shape [${shape}] has no trailing [re, im] axis of 2, but mo_metadata declares complex coefficients (spec § 5.3)`]);
      } else if (!cplx && shape.length !== 2) {
        errors.push(["E_MEMBER_MALFORMED", role, `${sec}${role} shape [${shape}] is not [n_mo, n_ao], and mo_metadata does not declare complex coefficients (spec § 5.3)`]);
      }
    }
  }

  if (isObject(meta)) {
    const sets: Array<[string, Obj]> = [["", meta]];
    for (const ch of ["alpha", "beta"]) {
      const block = meta[ch];
      if (isObject(block)) sets.push([`${ch}.`, block]);
    }
    for (const [prefix, block] of sets) {
      const labels = block.symmetry_labels;
      const energies = block.energies;
      if (Array.isArray(labels) && Array.isArray(energies) && labels.length !== energies.length) {
        errors.push([
          "E_MEMBER_MALFORMED",
          "mo_metadata",
          `${sec}mo_metadata ${prefix}symmetry_labels has ${labels.length} entries but ${prefix}energies has ${energies.length}; there is one label per MO (spec § 5.3)`,
        ]);
      }
    }
  }

  const nAo = basisNAo(basis);
  for (const role of ["overlap", "density_matrix", "density_matrix_alpha", "density_matrix_beta"]) {
    const shape = has(members, role) ? shapeOf(members[role]) : null;
    if (shape && nAo !== null && !sameExtents(shape, [nAo, nAo])) {
      errors.push(["E_MEMBER_MALFORMED", role, `${sec}${role} shape [${shape}] is not [n_ao, n_ao] = [${nAo}, ${nAo}] for the basis (spec § 5.3)`]);
    }
  }
  const unrestricted = has(members, "mo_coefficients_alpha") && has(members, "mo_coefficients_beta");
  if (has(members, "density_matrix") && !has(members, "mo_coefficients")) {
    errors.push([
      "E_MEMBER_MALFORMED",
      "density_matrix",
      `${sec}density_matrix accompanies mo_coefficients only; an unrestricted wavefunction carries density_matrix_alpha and density_matrix_beta (spec § 5.3)`,
    ]);
  }
  for (const [role, other] of [
    ["density_matrix_alpha", "density_matrix_beta"],
    ["density_matrix_beta", "density_matrix_alpha"],
  ] as const) {
    if (!has(members, role)) continue;
    if (!unrestricted) {
      errors.push([
        "E_MEMBER_MALFORMED",
        role,
        `${sec}${role} accompanies mo_coefficients_alpha and mo_coefficients_beta only; a single-set wavefunction carries density_matrix (spec § 5.3)`,
      ]);
    } else if (!has(members, other)) {
      errors.push(["E_MEMBER_MALFORMED", role, `${sec}${role} without ${other}; the two come together (spec § 5.3)`]);
    }
  }

  const nEl = isObject(basis) ? basis.n_electrons : undefined;
  if (!isNumber(nEl)) return errors;
  if (isObject(meta) && meta.occupation_semantics !== "transition_weight") {
    const lists: unknown[] = has(meta, "occupations") ? [meta.occupations] : [];
    for (const ch of ["alpha", "beta"]) {
      const block = meta[ch];
      if (isObject(block) && has(block, "occupations")) lists.push(block.occupations);
    }
    if (lists.length && lists.every((occ) => Array.isArray(occ) && occ.every((x) => isNumber(x)))) {
      const total = (lists as number[][]).flat().reduce((a, b) => a + b, 0);
      if (!(Math.abs(total - nEl) <= N_ELECTRONS_TOLERANCE)) {
        errors.push([
          "E_MEMBER_MALFORMED",
          "basis",
          `${sec}basis.n_electrons is ${nEl} but the mo_metadata occupations sum to ${g10(total)}; they must agree within 1e-6 (spec § 5.3)`,
        ]);
      }
    }
  }
  const provN = isObject(provenance) ? provenance.n_electrons : undefined;
  if (isNumber(provN) && provN !== nEl) {
    errors.push([
      "E_ROOT_INCONSISTENT",
      "basis",
      `${sec}basis.n_electrons is ${nEl} but provenance.n_electrons is ${provN}; they state one electron count and must be equal (spec § 3.2, § 5.3)`,
    ]);
  }
  return errors;
}

/** The volume kinds whose `data` shape is checked (spec § 5.2): rank, fixed trailing extent, form. */
const VOLUME_DATA_SHAPES: ReadonlyMap<string, readonly [number, number | null, string]> = new Map([
  ["volume.laplacian", [3, null, "[nx, ny, nz]"]],
  ["volume.vector", [4, 3, "[nx, ny, nz, 3]"]],
  ["volume.complex", [4, 2, "[nx, ny, nz, 2]"]],
  ["volume.series", [4, null, "[n, nx, ny, nz]"]],
] as const);

/**
 * Data shape and series axis of a volume section (spec § 5.2): `volume.laplacian`
 * data is `[nx, ny, nz]`, `volume.vector` `[nx, ny, nz, 3]`, `volume.complex`
 * `[nx, ny, nz, 2]` and `volume.series` `[n, nx, ny, nz]`; a series `axis` has
 * one `values` entry, and one `labels` entry when it has labels, per slab. With
 * `kindsById` (the references level), its `trajectory_ref` must name a
 * `trajectory` or `reaction.path` section. Returns `[code, role, message]`.
 */
function volumeErrors(
  sid: string,
  kind: string,
  members: Obj,
  docs: Map<string, unknown>,
  kindsById: ReadonlyMap<string, string> | null,
): Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> {
  const rule = VOLUME_DATA_SHAPES.get(kind);
  if (!rule) return [];
  const [rank, trailing, form] = rule;
  const sec = `section ${quote(sid)}: `;
  const errors: Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string, string]> = [];
  let shape = shapeOf(members.data);
  if (shape !== null && (shape.length !== rank || (trailing !== null && shape[shape.length - 1] !== trailing))) {
    errors.push(["E_MEMBER_MALFORMED", "data", `${sec}${kind} data shape [${shape}] is not ${form} (spec § 5.2)`]);
    shape = null;
  }
  const axis = docs.get("axis");
  if (kind !== "volume.series" || !isObject(axis)) return errors;
  for (const key of ["values", "labels"]) {
    const seq = axis[key];
    if (shape !== null && Array.isArray(seq) && seq.length !== shape[0]) {
      errors.push([
        "E_MEMBER_MALFORMED",
        "axis",
        `${sec}axis.${key} has ${seq.length} entries but data [${shape}] has ${shape[0]} slabs; it needs one per slab (spec § 5.2)`,
      ]);
    }
  }
  if (kindsById && has(axis, "trajectory_ref")) {
    const ref = axis.trajectory_ref;
    const target = typeof ref === "string" ? kindsById.get(ref) : undefined;
    if (target !== "trajectory" && target !== "reaction.path") {
      errors.push([
        "E_REF_UNRESOLVED",
        "axis",
        `${sec}axis.trajectory_ref ${show(ref)} does not name a trajectory or reaction.path section (spec § 5.2)`,
      ]);
    }
  }
  return errors;
}

/**
 * True for a member whose payload is a volume grid (`$defs/VolumeGridPayload`):
 * the `grid` of a canonical `volume.*` or `basis.ao` section, and the
 * `volume_grid` of a `reaction.path` section.
 */
function isGridPayload(kind: string, role: string): boolean {
  if (role === "volume_grid") return kind === "reaction.path";
  return role === "grid" && CANONICAL_KIND_SET.has(kind) && (kind.startsWith("volume.") || kind === "basis.ao");
}

const isBool3 = (v: unknown): v is boolean[] => Array.isArray(v) && v.length === 3 && v.every((x) => typeof x === "boolean");

function vec3(v: unknown): number[] | null {
  return Array.isArray(v) && v.length === 3 && v.every(isNumber) ? (v as number[]) : null;
}

const showVec = (v: number[]): string => `[${v.map((x) => Number(x.toPrecision(6))).join(", ")}]`;

/**
 * The unit name spec § 4.3 fixes for the grid `units.data` of a volume kind.
 * One grid payload serves every `volume.*` kind, so the schema cannot fix it
 * per kind.
 */
const VOLUME_DATA_UNITS: ReadonlyMap<string, string> = new Map([
  ["volume.density", "e_per_bohr3"],
  ["volume.transition_density", "e_per_bohr3"],
  ["volume.potential", "hartree_per_e"],
  ["volume.rdg", "dimensionless"],
  ["volume.laplacian", "e_per_bohr5"],
]);

/**
 * A grid `units.data` that is a unit name other than the one § 4.3 fixes for
 * the volume kind (`E_UNIT_CONTRADICTED`), as a message, or null. A string
 * that is not a unit name is `W_UNIT_UNKNOWN` alone.
 */
function gridDataUnitError(sid: string, kind: string, grid: Obj): string | null {
  const fixed = VOLUME_DATA_UNITS.get(kind);
  const units = grid["units"];
  if (fixed === undefined || !isObject(units)) return null;
  const unit = units["data"];
  if (typeof unit !== "string" || unit === fixed || !isUnitName(unit)) return null;
  return `section ${quote(sid)}: grid units.data ${quote(unit)} contradicts the unit ${quote(fixed)} that spec § 4.3 fixes for ${kind} data`;
}

/** Largest relative difference between a periodic grid span and its lattice vector (spec § 5.2). */
const GRID_SPAN_TOLERANCE = 1e-6;

/**
 * Reference and periodic-span rules of a grid payload (spec § 5.2). A grid
 * `structure_ref` names a `structure` section. The grid's structure is the
 * section-level `structure_ref`, else the grid's, else the archive's only
 * `structure`; when it resolves to a structure whose payload was read
 * (`cellsOf`), every axis a with `periodic[a]` true has `pbc[a]` true and a
 * span `voxel_vectors[a] × (shape[a] − 1)` (`includes_endpoint` true, the
 * default) or `voxel_vectors[a] × shape[a]` (false) equal to lattice row a in
 * bohr within 1e-6 relative. Malformed fields are left to the schema. Returns
 * `[code, message]` pairs.
 */
function gridErrors(
  sid: string,
  section: Obj,
  grid: Obj,
  kindsById: ReadonlyMap<string, string>,
  cellsOf: ReadonlyMap<string, Obj>,
  onlyStructure: string | null,
): Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string]> {
  const sec = `section ${quote(sid)}: `;
  const errors: Array<["E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED", string]> = [];
  const ref = grid.structure_ref;
  if (has(grid, "structure_ref") && (typeof ref !== "string" || kindsById.get(ref) !== "structure")) {
    errors.push(["E_REF_UNRESOLVED", `${sec}grid structure_ref ${show(ref)} does not name a structure section (spec § 5.2)`]);
  }
  const periodic = grid.periodic;
  if (!isBool3(periodic) || !periodic.some(Boolean)) return errors;
  let structure: string | null;
  if (has(section, "structure_ref")) structure = typeof section.structure_ref === "string" ? section.structure_ref : null;
  else if (typeof ref === "string") structure = ref;
  else structure = onlyStructure;
  const cell = structure !== null ? cellsOf.get(structure) : undefined;
  if (structure === null || cell === undefined) return errors;
  const pbc = cell.pbc ?? [false, false, false];
  const endpoint = has(grid, "includes_endpoint") ? grid.includes_endpoint : true;
  const shape = grid.shape;
  const vox = grid.voxel_vectors;
  if (
    !isBool3(pbc) ||
    typeof endpoint !== "boolean" ||
    !Array.isArray(shape) ||
    shape.length !== 3 ||
    !shape.every((n) => Number.isInteger(n)) ||
    !Array.isArray(vox) ||
    vox.length !== 3
  ) {
    return errors;
  }
  const lattice = cell.lattice_vectors;
  for (let a = 0; a < 3; a++) {
    if (!periodic[a]) continue;
    if (!pbc[a]) {
      errors.push([
        "E_MEMBER_MALFORMED",
        `${sec}grid axis ${a} is periodic but structure ${quote(structure)} has pbc[${a}] false; a grid axis can be periodic only along a periodic axis of its structure (spec § 5.2)`,
      ]);
      continue;
    }
    const step = vec3(vox[a]);
    const row = Array.isArray(lattice) && lattice.length === 3 ? vec3(lattice[a]) : null;
    if (step === null || row === null) continue;
    const n = (shape[a] as number) - (endpoint ? 1 : 0);
    const span = step.map((x) => x * n);
    const target = row.map((x) => x * BOHR_PER_ANGSTROM);
    const diff = Math.hypot(...span.map((s, c) => s - (target[c] as number)));
    if (diff > GRID_SPAN_TOLERANCE * Math.hypot(...target)) {
      errors.push([
        "E_MEMBER_MALFORMED",
        `${sec}periodic grid axis ${a} spans ${showVec(span)} bohr (voxel_vectors[${a}] × ${n}, includes_endpoint ${endpoint}) but lattice vector ${a} of structure ${quote(structure)} is ${showVec(target)} bohr; a periodic grid spans its cell within 1e-6 relative (spec § 5.2)`,
      ]);
    }
  }
  return errors;
}

function structurePayloadErrors(sid: string, payload: unknown): string[] {
  const sec = `section ${quote(sid)}: `;
  if (!isObject(payload)) return [`${sec}structure payload must be a JSON object`];

  const errors: string[] = [];
  const pbc = payload.pbc ?? null;
  const lattice = payload.lattice_vectors ?? null;
  let dim: unknown = payload.dimensionality ?? null;

  if (dim !== null && typeof dim !== "number") {
    errors.push(`${sec}dimensionality must be an integer, got ${show(dim)} (spec § 5.1)`);
    dim = null;
  }

  if (pbc === null) {
    if (lattice !== null) {
      errors.push(
        `${sec}structure payload has lattice_vectors but no 'pbc'; pbc is REQUIRED for periodic structures (spec § 5.1)`,
      );
    }
    if (dim !== null && dim !== 0) {
      errors.push(
        `${sec}dimensionality ${show(dim)} without 'pbc'; absent pbc means a molecule, so dimensionality must be 0 (spec § 5.1)`,
      );
    }
    return errors;
  }

  if (!Array.isArray(pbc) || pbc.length !== 3 || !pbc.every((x) => typeof x === "boolean")) {
    return [...errors, `${sec}pbc must be three booleans, got ${show(pbc)}`];
  }

  const nPeriodic = pbc.filter(Boolean).length;
  if (nPeriodic > 0 && !Array.isArray(lattice)) {
    errors.push(
      `${sec}pbc marks ${nPeriodic} periodic axis/axes but lattice_vectors is ${lattice === null ? "null or absent" : "not an array"}`,
    );
  }
  if (dim !== null && dim !== nPeriodic) {
    errors.push(`${sec}dimensionality ${show(dim)} != sum(pbc) ${nPeriodic} (spec § 5.1 invariant)`);
  }
  return errors;
}

/**
 * The sum of the `atomic_number` values of a `structure` payload's atoms, or
 * null when the payload has no `atoms` array or an atomic number is missing or
 * not a non-negative integer (spec § 5.1, Chemical identity).
 */
function sumAtomicNumbers(payload: unknown): number | null {
  if (!isObject(payload) || !Array.isArray(payload.atoms)) return null;
  let sum = 0;
  for (const atom of payload.atoms) {
    const z = isObject(atom) ? atom.atomic_number : undefined;
    if (!isInt(z) || z < 0) return null;
    sum += z;
  }
  return sum;
}

type IdentityCode = "E_CHARGE_MULTIPLICITY" | "E_ROOT_INCONSISTENT";

/**
 * The chemical identity of a `structure` section (spec § 5.1, Chemical
 * identity). With an integer `charge` and the atoms' summed atomic numbers
 * `sumZ`, the electron count N = ΣZ − `charge` must not be negative, and an
 * integer `multiplicity` (≥ 1) must have `multiplicity` − 1 at most N and of
 * the parity of N (`E_CHARGE_MULTIPLICITY`). In an archive with exactly one
 * `structure` section (`onlyStructure`), `charge` and `multiplicity` must
 * equal an integer `provenance.charge` / `provenance.multiplicity`
 * (`E_ROOT_INCONSISTENT`). Returns `[code, field, message]` triples.
 */
function identityErrors(
  sid: string,
  section: Obj,
  sumZ: number | null,
  provenance: unknown,
  onlyStructure: boolean,
): Array<[IdentityCode, "charge" | "multiplicity", string]> {
  const errors: Array<[IdentityCode, "charge" | "multiplicity", string]> = [];
  const charge = isInt(section.charge) ? section.charge : null;
  const mult = isInt(section.multiplicity) && section.multiplicity >= 1 ? section.multiplicity : null;
  if (charge !== null && sumZ !== null) {
    const electrons = sumZ - charge;
    if (electrons < 0) {
      errors.push([
        "E_CHARGE_MULTIPLICITY",
        "charge",
        `structure ${quote(sid)}: charge ${charge} leaves ${electrons} electrons for atoms whose atomic numbers sum to ${sumZ} (spec § 5.1)`,
      ]);
    } else if (mult !== null && (mult - 1 > electrons || (electrons - (mult - 1)) % 2 !== 0)) {
      errors.push([
        "E_CHARGE_MULTIPLICITY",
        "multiplicity",
        `structure ${quote(sid)}: multiplicity ${mult} is impossible with ${electrons} electrons (charge ${charge}): the unpaired electrons must not exceed the electron count and must have its parity (spec § 5.1)`,
      ]);
    }
  }
  if (onlyStructure && isObject(provenance)) {
    for (const [field, value] of [
      ["charge", charge],
      ["multiplicity", mult],
    ] as const) {
      const stated = provenance[field];
      if (value !== null && isInt(stated) && stated !== value) {
        errors.push([
          "E_ROOT_INCONSISTENT",
          field,
          `structure ${quote(sid)}: ${field} ${value} differs from provenance.${field} ${stated}, which describes the archive's only structure (spec § 5.1)`,
        ]);
      }
    }
  }
  return errors;
}

/** The participant lists of a `reaction.scheme` step (spec § 5.14). */
const SCHEME_ROLES = ["reactants", "products", "reagents"] as const;

type SchemeCode = "E_MEMBER_MALFORMED" | "E_REF_UNRESOLVED";

/**
 * The rules of a `reaction.scheme` payload that JSON Schema cannot state
 * (spec § 5.14): the step `id` values are unique and no `structure_ref` occurs
 * twice within one role of one step (`E_MEMBER_MALFORMED`); with `kindsById`
 * (the references level), every participant's `structure_ref` names a
 * `structure` section (`E_REF_UNRESOLVED`). Returns `[code, message]` pairs.
 */
function reactionSchemeErrors(
  sid: string,
  payload: unknown,
  kindsById: ReadonlyMap<string, string> | null,
): Array<[SchemeCode, string]> {
  const errors: Array<[SchemeCode, string]> = [];
  const steps = isObject(payload) ? payload.steps : undefined;
  if (!Array.isArray(steps)) return errors;
  const sec = `section ${quote(sid)}: `;
  const seen = new Set<string>();
  for (let k = 0; k < steps.length; k++) {
    const step: unknown = steps[k];
    if (!isObject(step)) continue;
    const stepId = step.id;
    if (typeof stepId === "string") {
      if (seen.has(stepId)) {
        errors.push(["E_MEMBER_MALFORMED", `${sec}step id ${quote(stepId)} is used by more than one step (spec § 5.14)`]);
      }
      seen.add(stepId);
    }
    for (const role of SCHEME_ROLES) {
      const items = step[role];
      if (!Array.isArray(items)) continue;
      const counts = new Map<string, number>();
      for (const item of items) {
        const ref = isObject(item) ? item.structure_ref : undefined;
        if (typeof ref !== "string") continue;
        counts.set(ref, (counts.get(ref) ?? 0) + 1);
        if (kindsById !== null && kindsById.get(ref) !== "structure") {
          errors.push([
            "E_REF_UNRESOLVED",
            `${sec}step ${k} ${role} structure_ref ${quote(ref)} does not name a structure section (spec § 5.14)`,
          ]);
        }
      }
      const repeated = [...counts].filter(([, n]) => n > 1).map(([ref]) => ref);
      for (const ref of repeated.sort()) {
        errors.push(["E_MEMBER_MALFORMED", `${sec}step ${k} names ${quote(ref)} more than once among its ${role} (spec § 5.14)`]);
      }
    }
    if (isObject(step.energies)) {
      errors.push(...reactionEnergiesErrors(`${sec}step ${k}`, step, step.energies, kindsById));
    }
  }
  return errors;
}

/** Reaction energy quantities and the term field of each absolute value (spec § 5.14). */
const REACTION_QUANTITIES = [
  ["energy", "energy_eh"],
  ["enthalpy", "enthalpy_eh"],
  ["gibbs_free_energy", "gibbs_free_energy_eh"],
] as const;
const HARTREE_TO_KCAL_MOL = 627.509474; // CODATA 2018
const REACTION_DELTA_TOL_EH = 1e-8;
const REACTION_KCAL_REL_TOL = 1e-6;
const REACTION_COEFFICIENT_TOL = 1e-9;

/** Sum without cancellation loss beyond one rounding, as Python's math.fsum. */
function fsum(values: readonly number[]): number {
  const partials: number[] = [];
  for (let x of values) {
    let i = 0;
    for (let y of partials) {
      if (Math.abs(x) < Math.abs(y)) [x, y] = [y, x];
      const hi = x + y;
      const lo = y - (hi - x);
      if (lo !== 0) partials[i++] = lo;
      x = hi;
    }
    partials.length = i;
    partials.push(x);
  }
  return partials.reduce((a, b) => a + b, 0);
}

/**
 * A step's `energies` (spec § 5.14, *Reaction energies*): no energies on a
 * resonance arrow; one term per reactant and product with the signed
 * coefficient; each delta equals the sum of its terms and its kcal/mol value
 * the conversion; one method (`E_MEMBER_MALFORMED`); with `kindsById`, every
 * `run_ref` names a `run.record` (`E_REF_UNRESOLVED`). Mirrors
 * _reaction_energies_errors() in python/qvf/validate.py.
 */
function reactionEnergiesErrors(
  where: string,
  step: Obj,
  energies: Obj,
  kindsById: ReadonlyMap<string, string> | null,
): Array<[SchemeCode, string]> {
  const errors: Array<[SchemeCode, string]> = [];
  const bad = (msg: string): void => {
    errors.push(["E_MEMBER_MALFORMED", `${where} energies: ${msg} (spec § 5.14)`]);
  };
  if (step.arrow === "resonance") bad("a resonance arrow has no reaction energy");
  if (!Array.isArray(energies.terms)) return errors;
  const expected = new Map<string, number>();
  for (const [role, sign] of [["reactants", -1], ["products", 1]] as const) {
    const items = step[role];
    if (!Array.isArray(items)) continue;
    for (const item of items) {
      if (!isObject(item) || typeof item.structure_ref !== "string") continue;
      const c = item.coefficient === undefined ? 1 : item.coefficient;
      if (isNumber(c)) expected.set(item.structure_ref, (expected.get(item.structure_ref) ?? 0) + sign * c);
    }
  }
  for (const [ref, c] of [...expected]) if (c === 0) expected.delete(ref);
  const terms = (energies.terms as unknown[]).filter(isObject);
  const seen = new Set<string>();
  terms.forEach((term, t) => {
    const ref = term.structure_ref;
    const c = term.coefficient;
    if (typeof ref !== "string") return;
    if (seen.has(ref)) {
      bad(`structure ${quote(ref)} has more than one term`);
      return;
    }
    seen.add(ref);
    const want = expected.get(ref);
    if (want === undefined) {
      bad(`term ${t} names ${quote(ref)}, which is not a reactant or product`);
    } else if (isNumber(c) && Math.abs(c - want) > REACTION_COEFFICIENT_TOL * Math.max(1, Math.abs(c))) {
      bad(`term ${t} coefficient ${c} for ${quote(ref)}, but the step gives ${want}`);
    }
    const runRef = term.run_ref;
    if (kindsById !== null && typeof runRef === "string" && kindsById.get(runRef) !== "run.record") {
      errors.push([
        "E_REF_UNRESOLVED",
        `${where} energies: term ${t} run_ref ${quote(runRef)} does not name a run.record section (spec § 5.14)`,
      ]);
    }
  });
  for (const ref of [...expected.keys()].filter((r) => !seen.has(r)).sort()) {
    bad(`no term for participant ${quote(ref)}`);
  }
  const methods = terms.map((term) => term.method).filter(isObject);
  if (methods.some((m) => !jsonEqual(m, methods[0]))) bad("the terms have different methods");
  const deltas = energies.deltas;
  if (!isObject(deltas)) return errors;
  for (const [quantity, field] of REACTION_QUANTITIES) {
    const delta = deltas[quantity];
    if (!isObject(delta) || !isNumber(delta.value_eh)) continue;
    const value = delta.value_eh;
    const missing = terms.findIndex((term) => !isNumber(term[field]));
    if (missing >= 0) {
      bad(`${quantity} is given but term ${missing} has no ${field}`);
    } else if (terms.every((term) => isNumber(term.coefficient))) {
      const total = fsum(terms.map((term) => (term.coefficient as number) * (term[field] as number)));
      if (Math.abs(total - value) > REACTION_DELTA_TOL_EH) {
        bad(`${quantity} value_eh ${value} is not the sum of its terms, ${total}`);
      }
    }
    const kcal = delta.kcal_mol;
    if (isNumber(kcal)) {
      const converted = value * HARTREE_TO_KCAL_MOL;
      if (Math.abs(kcal - converted) > REACTION_KCAL_REL_TOL * Math.max(1, Math.abs(converted))) {
        bad(`${quantity} kcal_mol ${kcal} is not value_eh x 627.509474 = ${converted}`);
      }
    }
  }
  return errors;
}

/** The first bytes of a PNG and of a JPEG file (spec § 5.15). */
const PNG_SIGNATURE = Uint8Array.of(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a);
const JPEG_SIGNATURE = Uint8Array.of(0xff, 0xd8, 0xff);

const startsWith = (bytes: Uint8Array, prefix: Uint8Array): boolean =>
  bytes.length >= prefix.length && prefix.every((b, k) => bytes[k] === b);

/**
 * The bytes of a `preview` member against its `media_type` (spec § 5.15): an
 * `image/png` member begins with the PNG signature and an `image/jpeg` member
 * with FF D8 FF (`E_MEMBER_MALFORMED`); an `image/svg+xml` member decodes as
 * UTF-8 (`E_UTF8`) and contains `<svg` (`E_MEMBER_MALFORMED`). The
 * `media_type` values themselves are left to the JSON Schema. Returns
 * `[code, message]` pairs.
 */
function previewMemberErrors(
  sid: string,
  role: string,
  member: Obj,
  bytes: Uint8Array,
): Array<["E_MEMBER_MALFORMED" | "E_UTF8", string]> {
  const mediaType = member.media_type;
  const where = `section ${quote(sid)} member ${quote(role)}`;
  if (mediaType === "image/png" && !startsWith(bytes, PNG_SIGNATURE)) {
    return [["E_MEMBER_MALFORMED", `${where}: image/png bytes do not begin with the PNG signature (spec § 5.15)`]];
  }
  if (mediaType === "image/jpeg" && !startsWith(bytes, JPEG_SIGNATURE)) {
    return [["E_MEMBER_MALFORMED", `${where}: image/jpeg bytes do not begin with FF D8 FF (spec § 5.15)`]];
  }
  if (mediaType === "image/svg+xml") {
    let text: string;
    try {
      text = decodeUtf8(bytes);
    } catch (err) {
      return [["E_UTF8", `${where}: image/svg+xml bytes do not decode as UTF-8 (${errorText(err)}; spec § 5.15)`]];
    }
    if (!text.includes("<svg")) {
      return [["E_MEMBER_MALFORMED", `${where}: image/svg+xml bytes contain no <svg element (spec § 5.15)`]];
    }
  }
  return [];
}

/**
 * Validate a `.qvf` archive held in memory (spec § 6).
 *
 * Returns `{ ok, findings, nSections }`. It does not throw for anything that
 * is wrong with the archive: an unreadable ZIP is the finding
 * `E_ZIP_UNREADABLE`. A `TypeError` reports an invalid `options` argument.
 */
export async function validate(
  bytes: Uint8Array | ArrayBuffer,
  options: ValidateOptions = {},
): Promise<ValidationReport> {
  const limits = resolveLimits(options.limits);
  const level = options.level ?? "all";
  if (!LEVELS.includes(level)) {
    throw new TypeError(`unknown validation level ${JSON.stringify(level)}`);
  }
  const rank = LEVELS.indexOf(level);
  if (rank >= AT_PAYLOADS && !globalThis.crypto?.subtle) {
    // An environment problem, not a property of the archive: do not report it as a finding.
    throw new Error("WebCrypto (globalThis.crypto.subtle) is required to verify sha256 checksums");
  }
  const out = new Findings();
  let nSections = 0;
  try {
    const data = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
    let zip: ZipArchive | undefined;
    try {
      zip = ZipArchive.open(data);
    } catch (err) {
      out.error("E_ZIP_UNREADABLE", "archive", `not a readable zip: ${errorText(err)}`);
    }
    if (zip) nSections = await run(zip, limits, rank, options.schemaCheck ?? true, out);
  } catch (err) {
    out.error("E_ZIP_UNREADABLE", "archive", `the archive could not be processed: ${errorText(err)}`);
  }
  return {
    ok: !out.items.some((f) => f.severity === "error"),
    findings: out.items,
    nSections,
  };
}

async function run(zip: ZipArchive, lim: Limits, rank: number, schemaCheck: boolean, out: Findings): Promise<number> {
  const entries = zip.entries;

  // --- container: central-directory checks, before decompressing anything ---
  const counts = new Map<string, number>();
  for (const e of entries) counts.set(e.name, (counts.get(e.name) ?? 0) + 1);
  for (const [name, n] of counts) {
    if (n > 1 && name === "manifest.json") {
      // An archive updated in place: the last manifest.json is the manifest,
      // every earlier one is superseded (spec § 2.3).
      for (let copy = 1; copy < n; copy++) {
        out.warning(
          "W_ORPHAN_MEMBER",
          `entry:${name}`,
          `ZIP entry ${quote(name)} (${copy} of ${n}) is superseded by the last ${quote(name)} entry (archive updated in place; consumers ignore it)`,
        );
      }
    } else if (n > 1) {
      out.error(
        "E_ZIP_DUPLICATE_ENTRY",
        `entry:${name}`,
        `ZIP entry name ${quote(name)} occurs ${n} times; readers disagree on which entry wins`,
      );
    }
  }
  for (const name of counts.keys()) {
    const why = pathProblem(name);
    if (why) out.error("E_MEMBER_PATH", `entry:${name}`, `ZIP entry ${quote(name)}: illegal path (${why})`);
  }

  const unsafe = new Set<string>(); // entries that must not be decompressed
  let total = 0;
  for (const e of entries) {
    total += e.size;
    if (e.size > lim.max_member_bytes) {
      out.error(
        "E_MEMBER_TOO_LARGE",
        `entry:${e.name}`,
        `ZIP entry ${quote(e.name)} exceeds the per-member size cap (${e.size} > ${lim.max_member_bytes} bytes uncompressed)`,
      );
      unsafe.add(e.name);
    }
    if (ratioExceeded(e, lim)) {
      const ratio = e.compressedSize ? String(Math.round(e.size / e.compressedSize)) : "inf";
      out.error(
        "E_COMPRESSION_RATIO",
        `entry:${e.name}`,
        `ZIP entry ${quote(e.name)} expands ${e.size} bytes from ${e.compressedSize} (ratio ${ratio} > ${lim.max_ratio})`,
      );
      unsafe.add(e.name);
    }
  }
  const payloadsBlocked = total > lim.max_total_bytes;
  if (payloadsBlocked) {
    out.error(
      "E_TOTAL_TOO_LARGE",
      "archive",
      `archive expands to ${total} bytes uncompressed, over the ${lim.max_total_bytes}-byte total cap; member payloads were not read`,
    );
  }
  if (rank < AT_MANIFEST) return 0;

  // --- manifest ---
  const manifestEntry = zip.find("manifest.json");
  if (!manifestEntry) {
    out.error("E_MANIFEST_MISSING", "manifest.json", "missing manifest.json");
    return 0;
  }
  if (manifestEntry.size > lim.max_json_bytes) {
    out.error(
      "E_JSON_TOO_LARGE",
      "manifest.json",
      `manifest.json is ${manifestEntry.size} bytes, over the ${lim.max_json_bytes}-byte JSON cap`,
    );
    return 0;
  }
  if (unsafe.has("manifest.json")) return 0;
  let raw: Uint8Array;
  try {
    raw = zip.read(manifestEntry, { maxBytes: lim.max_member_bytes });
  } catch (err) {
    out.error("E_ZIP_UNREADABLE", "entry:manifest.json", `ZIP entry 'manifest.json' cannot be read: ${errorText(err)}`);
    return 0;
  }
  if (jsonDepth(raw, lim.max_json_depth) > lim.max_json_depth) {
    out.error("E_JSON_TOO_DEEP", "manifest.json", `manifest.json nests deeper than ${lim.max_json_depth} levels`);
    return 0;
  }
  let manifest: unknown;
  try {
    manifest = parseJsonBytes(raw);
  } catch (err) {
    out.error("E_MANIFEST_JSON", "manifest.json", `manifest.json is not valid JSON: ${errorText(err)}`);
    return 0;
  }
  if (!isObject(manifest)) {
    out.error("E_MANIFEST_JSON", "manifest.json", "manifest.json must be a JSON object");
    return 0;
  }

  // Root sanity. The withdrawn version 2 is read as 1 for every check, the
  // JSON Schema check included (spec § 3.1).
  const version = manifest.qvf_version;
  let checkedManifest: Obj = manifest;
  if (version === 2) {
    out.warning(
      "W_VERSION_WITHDRAWN",
      "/qvf_version",
      "qvf_version 2 was withdrawn: the archive is checked as version 1, which producers must write (spec § 3.1)",
    );
    checkedManifest = { ...manifest, qvf_version: 1 };
  } else if (typeof version !== "number" || version !== 1) {
    out.error("E_VERSION", "/qvf_version", `qvf_version must be 1, got ${show(version)}`);
  }
  let source: Obj = {};
  if (manifest.source !== undefined) {
    if (isObject(manifest.source)) source = manifest.source;
    else out.error("E_SOURCE", "/source", "source must be a JSON object");
  }
  for (const key of ["program", "version", "calculation"]) {
    if (!has(source, key)) out.error("E_SOURCE", "/source", `source.${key} is missing`);
  }
  if (!truthy(source.program) || !truthy(source.version)) {
    out.error("E_SOURCE", "/source", "source.program and source.version must be non-empty");
  }
  if (rank >= AT_SEMANTICS) {
    for (const [where, msg] of rootConsistencyErrors(manifest)) out.error("E_ROOT_INCONSISTENT", where, msg);
  }

  // The JSON Schema, through the checker generated from it.
  if (schemaCheck) {
    for (const e of checkManifest(checkedManifest)) {
      out.error("E_SCHEMA", e.pointer || "manifest.json", `${e.location || "(root)"}: ${e.message}`);
    }
  } else {
    out.warning("W_SCHEMA_SKIPPED", "manifest.json", "the JSON Schema check was disabled; structural schema check skipped");
  }

  // Declared units of {value, units} objects in provenance (spec § 4.3).
  if (rank >= AT_SEMANTICS && isObject(manifest.provenance)) {
    for (const [key, value] of Object.entries(manifest.provenance)) {
      if (!isObject(value) || !has(value, "value") || !has(value, "units")) continue;
      if (canonicalUnit(value.units) === null) {
        out.warning(
          "W_UNIT_UNKNOWN",
          `/provenance/${escapePointer(key)}/units`,
          `provenance.${key}.units ${show(value.units)} is not a unit name of spec § 4.3 or a legacy spelling`,
        );
      }
    }
  }

  let sections: unknown[] = [];
  if (manifest.sections !== undefined) {
    if (Array.isArray(manifest.sections)) sections = manifest.sections;
    else out.error("E_SECTION_MALFORMED", "/sections", "sections must be an array");
  }
  const sectionIds = new Set<string>();
  const kindsById = new Map<string, string>(); // first occurrence wins
  const structureIds: string[] = [];
  // Section id -> kind (undefined when not a string), for reference resolution
  // (spec § 4.1, § 5.5). With duplicate ids the first section wins.
  const sectionKinds = new Map<string, string | undefined>();
  for (const s of sections) {
    if (isObject(s) && typeof s.id === "string") {
      sectionIds.add(s.id);
      if (typeof s.kind === "string" && !kindsById.has(s.id)) kindsById.set(s.id, s.kind);
      if (s.kind === "structure") structureIds.push(s.id);
      if (!sectionKinds.has(s.id)) sectionKinds.set(s.id, typeof s.kind === "string" ? s.kind : undefined);
    }
  }
  const nAtomsOf = new Map<string, number>(); // structure id -> atom count
  const sumZOf = new Map<string, number>(); // structure id -> sum of atomic numbers
  // structure sections, whose chemical identity is checked once every section has been read
  const identitySections: Array<{ i: number; sid: string; s: Obj }> = [];
  // Excited states (spec § 5.12): orbital rows of every wavefunction.gto, state
  // indices of every excited_states payload, and the sections checked once
  // every payload has been read.
  const moRowsOf = new Map<string, Map<string, number>>();
  for (const s of sections) {
    if (isObject(s) && s.kind === "wavefunction.gto" && typeof s.id === "string" && !moRowsOf.has(s.id)) {
      moRowsOf.set(s.id, moRows(s.members));
    }
  }
  const stateIndicesOf = new Map<string, Set<number> | null>();
  const ntoSections: NtoSection[] = [];
  const excitedSections: Array<{ i: number; sid: string; s: Obj; payload: unknown }> = [];
  const statesRefSpectra: Array<{ i: number; sid: string; payload: Obj }> = [];
  const transitionDensities: Array<{ i: number; sid: string; s: Obj }> = [];
  const wavefunctionVolumes: Array<{ i: number; sid: string; kind: string; s: Obj }> = [];
  const cellsOf = new Map<string, Obj>(); // structure id -> structure payload
  const periodicOf = new Map<string, boolean>(); // structure id -> any(pbc)
  // wavefunction.bloch sections, checked once every structure has been read
  const blochSections: Array<{ i: number; sid: string; s: Obj; members: Obj; docs: Map<string, unknown> }> = [];
  // forces / hessian sections, checked once every structure has been read
  const derivativeSections: Array<{ i: number; sid: string; kind: string; s: Obj; members: Obj; meta: unknown }> = [];
  // standalone bonds sections with a structure_ref, checked the same way (spec § 5.1)
  const bondsSections: Array<{ i: number; sid: string; s: Obj; payload: unknown }> = [];
  // grid payloads, checked once every structure has been read (spec § 5.2)
  const gridMembers: Array<{ where: string; sid: string; s: Obj; grid: Obj }> = [];
  // properties.* sections, checked the same way (the Born charges need the structures)
  const propertySections: Array<{ i: number; sid: string; kind: string; s: Obj; payload: unknown }> = [];
  const declaredExtensions: Obj = isObject(manifest.extensions) ? manifest.extensions : {};
  const declaredPaths = new Set<string>();
  const seenIds = new Set<string>();

  for (let i = 0; i < sections.length; i++) {
    const s = sections[i];
    const here = sectionPointer(i);
    if (!isObject(s)) {
      out.error("E_SECTION_MALFORMED", here, `section ${i} must be a JSON object`);
      continue;
    }

    let sid: string;
    if (typeof s.id !== "string") {
      out.error("E_SECTION_MALFORMED", here, `section ${i}: id must be a string, got ${show(s.id)}`);
      sid = `<section ${i}>`;
    } else {
      sid = s.id;
      if (seenIds.has(sid)) out.error("E_SECTION_ID_DUPLICATE", `${here}/id`, `duplicate section id ${quote(sid)}`);
      seenIds.add(sid);
      if (RESERVED_SECTION_IDS.has(sid)) {
        out.error(
          "E_RESERVED_ID",
          `${here}/id`,
          `section id ${quote(sid)} is reserved (a viewer_defaults key, spec § 3.1, § 3.4)`,
        );
      }
    }

    let kind = "";
    if (typeof s.kind !== "string") {
      out.error("E_SECTION_MALFORMED", here, `section ${quote(sid)}: kind must be a string, got ${show(s.kind)}`);
    } else {
      kind = s.kind;
      if (!CANONICAL_KIND_SET.has(kind) && !kind.startsWith("x_")) {
        out.error("E_KIND_UNKNOWN", `${here}/kind`, `section ${quote(sid)}: unknown non-vendor kind ${quote(kind)}`);
      }
    }

    // critical governance (spec § 7.4).
    if (truthy(s.critical) && kind.startsWith("x_")) {
      const dot = kind.indexOf(".");
      const vendor = dot < 0 ? kind : kind.slice(0, dot);
      if (!has(declaredExtensions, vendor)) {
        out.error(
          "E_CRITICAL_UNDECLARED",
          `${here}/critical`,
          `section ${quote(sid)} is critical but its namespace ${quote(vendor)} is not declared in root extensions`,
        );
      }
    }

    let members: Obj = {};
    if (isObject(s.members)) {
      members = s.members;
    } else {
      out.error("E_SECTION_MALFORMED", here, `section ${quote(sid)}: members must be a JSON object`);
    }

    const docs = new Map<string, unknown>(); // parsed JSON members of this section, by role
    const binaries = new Map<string, Uint8Array>(); // trajectory-family energies and wavefunction.gto overlap bytes
    for (const [role, m] of Object.entries(members)) {
      const mwhere = memberPointer(i, role);
      if (!isObject(m) || typeof m.path !== "string") {
        out.error("E_MEMBER_MALFORMED", mwhere, `section ${quote(sid)} member ${quote(role)}: malformed spec`);
        continue;
      }
      const mpath = m.path;
      declaredPaths.add(mpath);
      const why = pathProblem(mpath);
      if (why) {
        out.error("E_MEMBER_PATH", mwhere, `member ${quote(mpath)}: illegal path (${why})`);
        continue;
      }
      const entry = zip.find(mpath);
      if (!entry) {
        out.error("E_MEMBER_MISSING", mwhere, `member ${quote(mpath)} declared but not in archive`);
        continue;
      }
      // Declared size against the central directory, before reading.
      const sizeProblem = declaredSizeProblem(m, entry.size, mpath);
      if (sizeProblem) out.error(sizeProblem[0], mwhere, sizeProblem[1]);
      if (rank < AT_PAYLOADS) continue;
      if (payloadsBlocked || unsafe.has(mpath)) continue; // the container finding already explains why
      await checkPayload({ zip, entry, m, mpath, mwhere, sid, kind, role, lim, rank, schemaCheck, docs, binaries, out });
    }

    // 1-D spectra: one shared array length (spec § 5.4).
    if (rank >= AT_SEMANTICS && isSpectrum1d(kind)) {
      for (const msg of spectrumPayloadErrors(sid, docs.get("spectrum"))) {
        out.error("E_PAYLOAD_SCHEMA", memberPointer(i, "spectrum"), msg);
      }
    }

    if (rank >= AT_SEMANTICS && kind === "bands" && isObject(members)) {
      for (const [role, msg] of bandsShapeErrors(sid, members)) {
        out.error("E_MEMBER_MALFORMED", `${here}/members/${escapePointer(role)}`, msg);
      }
    }
    if (rank >= AT_SEMANTICS && TRAJECTORY_KINDS.has(kind) && isObject(members)) {
      for (const [role, msg] of trajectoryErrors(sid, kind, members, docs.get("metadata"), binaries.get("energies"))) {
        out.error("E_MEMBER_MALFORMED", memberPointer(i, role), msg);
      }
    }
    if (kind === "structure") {
      const doc = docs.get("structure");
      if (isObject(doc) && Array.isArray(doc.atoms) && !nAtomsOf.has(sid)) nAtomsOf.set(sid, doc.atoms.length);
      if (!sumZOf.has(sid)) {
        const sumZ = sumAtomicNumbers(doc);
        if (sumZ !== null) sumZOf.set(sid, sumZ);
      }
      if (isObject(doc) && !cellsOf.has(sid)) cellsOf.set(sid, doc);
      if (isObject(doc) && !periodicOf.has(sid)) {
        periodicOf.set(sid, Array.isArray(doc.pbc) && doc.pbc.some((x) => x === true));
      }
      if (rank >= AT_SEMANTICS) identitySections.push({ i, sid, s });
    }
    if (rank >= AT_SEMANTICS && (kind === "forces" || kind === "hessian")) {
      derivativeSections.push({ i, sid, kind, s, members, meta: docs.get("meta") });
    }
    if (rank >= AT_REFERENCES && kind === "bonds" && has(s, "structure_ref")) {
      bondsSections.push({ i, sid, s, payload: docs.get("bonds") });
    }
    if (rank >= AT_SEMANTICS) {
      for (const [code, role, msg] of volumeErrors(sid, kind, members, docs, rank >= AT_REFERENCES ? kindsById : null)) {
        out.error(code, memberPointer(i, role), msg);
      }
      for (const [role, doc] of docs) {
        if (!isGridPayload(kind, role) || !isObject(doc)) continue;
        const msg = gridDataUnitError(sid, kind, doc);
        if (msg !== null) out.error("E_UNIT_CONTRADICTED", memberPointer(i, role), msg);
      }
    }
    if (rank >= AT_REFERENCES) {
      for (const [role, doc] of docs) {
        if (isGridPayload(kind, role) && isObject(doc)) gridMembers.push({ where: memberPointer(i, role), sid, s, grid: doc });
      }
    }
    if (rank >= AT_SEMANTICS && PROPERTY_KINDS.has(kind) && docs.has("properties")) {
      propertySections.push({ i, sid, kind, s, payload: docs.get("properties") });
    }
    if (rank >= AT_SEMANTICS && kind === "wavefunction.bloch") {
      blochSections.push({ i, sid, s, members, docs });
    }

    if (kind === "excited_states") {
      if (!stateIndicesOf.has(sid)) stateIndicesOf.set(sid, stateIndices(docs.get("states")));
      if (rank >= AT_SEMANTICS) excitedSections.push({ i, sid, s, payload: docs.get("states") });
    }
    if (rank >= AT_REFERENCES && isSpectrum1d(kind)) {
      const spectrum = docs.get("spectrum");
      if (isObject(spectrum) && has(spectrum, "states_ref")) statesRefSpectra.push({ i, sid, payload: spectrum });
    }
    if (rank >= AT_REFERENCES && kind === "volume.transition_density") transitionDensities.push({ i, sid, s });
    if (rank >= AT_REFERENCES && WAVEFUNCTION_VOLUME_KINDS.has(kind) && has(s, "wavefunction_ref")) {
      wavefunctionVolumes.push({ i, sid, kind, s });
    }
    if (rank >= AT_SEMANTICS && kind === "atom_properties" && isObject(members)) {
      for (const [role, msg] of atomPropertiesMetaErrors(sid, members, docs)) {
        out.error("E_MEMBER_MALFORMED", `${here}/members/${escapePointer(role)}`, msg);
      }
    }
    if (rank >= AT_SEMANTICS && kind === "wavefunction.gto" && isObject(members)) {
      ntoSections.push({ i, meta: docs.get("mo_metadata"), members });
      for (const msg of unrestrictedNtoErrors(docs.get("mo_metadata"), members)) {
        out.error("E_MEMBER_MALFORMED", memberPointer(i, "mo_metadata"), msg);
      }
      for (const [code, role, msg] of wavefunctionErrors(sid, members, docs, manifest.provenance)) {
        out.error(code, memberPointer(i, role), msg);
      }
      for (const msg of normalizationErrors(sid, docs.get("basis"), members.overlap, binaries.get("overlap"))) {
        out.error("E_NORMALIZATION", memberPointer(i, "overlap"), msg);
      }
    }
    // reaction.scheme (spec § 5.14): step ids and repeated participants, and,
    // from the references level, the structures the participants name.
    if (rank >= AT_SEMANTICS && kind === "reaction.scheme") {
      for (const [code, msg] of reactionSchemeErrors(sid, docs.get("scheme"), rank >= AT_REFERENCES ? kindsById : null)) {
        out.error(code, memberPointer(i, "scheme"), msg);
      }
    }

    if (rank < AT_REFERENCES) continue;

    // Cross-references.
    if (kind === "volume.difference") {
      for (const op of ["operand_a", "operand_b"]) {
        const ref = s[op];
        if (ref !== undefined && ref !== null && (typeof ref !== "string" || !sectionIds.has(ref))) {
          out.error("E_REF_UNRESOLVED", `${here}/${op}`, `section ${quote(sid)}: ${op} ${show(ref)} does not resolve`);
        }
      }
    }
    if (kind === "reaction.waypoints") {
      const ref = s.trajectory_ref;
      if (typeof ref !== "string" || !sectionIds.has(ref)) {
        out.error(
          "E_REF_UNRESOLVED",
          `${here}/trajectory_ref`,
          `section ${quote(sid)}: trajectory_ref ${show(ref)} does not resolve`,
        );
      } else if (sectionKinds.get(ref) !== "trajectory") {
        out.error(
          "E_REF_UNRESOLVED",
          `${here}/trajectory_ref`,
          `section ${quote(sid)}: trajectory_ref ${show(ref)} names a section of kind ${show(sectionKinds.get(ref))}, not 'trajectory' (spec § 5.5)`,
        );
      }
    }
    // Generic section references, on any kind (spec § 4.1).
    for (const [field, target] of SECTION_REFS) {
      const ref = s[field];
      if (ref === undefined || ref === null) continue;
      if (typeof ref !== "string" || !sectionIds.has(ref)) {
        out.error("E_REF_UNRESOLVED", `${here}/${field}`, `section ${quote(sid)}: ${field} ${show(ref)} does not resolve`);
      } else if (sectionKinds.get(ref) !== target) {
        out.error(
          "E_REF_UNRESOLVED",
          `${here}/${field}`,
          `section ${quote(sid)}: ${field} ${show(ref)} names a section of kind ${show(sectionKinds.get(ref))}, not ${quote(target)} (spec § 4.1)`,
        );
      }
    }
    if (kind === "vibrations" && has(s, "hessian_ref")) {
      const ref = s.hessian_ref;
      if (typeof ref !== "string" || kindsById.get(ref) !== "hessian") {
        out.error(
          "E_REF_UNRESOLVED",
          `${here}/hessian_ref`,
          `section ${quote(sid)}: hessian_ref ${show(ref)} does not name a hessian section (spec § 5.5)`,
        );
      }
    }
    if (kind === "job.spec") {
      // The payload's structure_ref names a structure section (spec § 5.9).
      const spec = docs.get("spec");
      if (isObject(spec) && has(spec, "structure_ref")) {
        const ref = spec.structure_ref;
        if (typeof ref !== "string" || kindsById.get(ref) !== "structure") {
          out.error(
            "E_REF_UNRESOLVED",
            memberPointer(i, "spec"),
            `section ${quote(sid)}: spec.structure_ref ${show(ref)} does not name a structure section (spec § 5.9)`,
          );
        }
      }
    }
    if (kind === "run.record") {
      // job_ref names a job.spec section (spec § 5.8).
      if (has(s, "job_ref")) {
        const ref = s.job_ref;
        if (typeof ref !== "string" || kindsById.get(ref) !== "job.spec") {
          out.error(
            "E_REF_UNRESOLVED",
            `${here}/job_ref`,
            `section ${quote(sid)}: job_ref ${show(ref)} does not name a job.spec section (spec § 5.8)`,
          );
        }
      }
      const files = docs.get("files");
      if (isObject(files)) {
        for (const frole of Object.keys(files)) {
          if (!has(members, frole)) {
            out.error(
              "E_FILES_INDEX",
              memberPointer(i, "files"),
              `section ${quote(sid)}: files index names role ${quote(frole)} with no matching member`,
            );
          }
        }
      }
    }
  }

  // Chemical identity of the structures (spec § 5.1).
  for (const d of identitySections) {
    const sumZ = sumZOf.get(d.sid) ?? null;
    for (const [code, field, msg] of identityErrors(d.sid, d.s, sumZ, manifest.provenance, structureIds.length === 1)) {
      out.error(code, `${sectionPointer(d.i)}/${field}`, msg);
    }
  }

  // forces / hessian: shapes against the structure, which may come later.
  const structures: StructureIndex | null =
    rank >= AT_REFERENCES
      ? { kindsById, nAtomsOf, onlyStructure: structureIds.length === 1 ? (structureIds[0] as string) : null }
      : null;
  for (const d of derivativeSections) {
    for (const [code, role, msg] of derivativeErrors(d.sid, d.kind, d.s, d.members, d.meta, structures)) {
      out.error(code, memberPointer(d.i, role), msg);
    }
  }
  // Standalone bonds: atom indices against the structure they name, which may
  // come later (spec § 5.1).
  if (structures) {
    for (const d of bondsSections) {
      for (const msg of bondsErrors(d.sid, d.s, d.payload, structures)) {
        out.error("E_MEMBER_MALFORMED", memberPointer(d.i, "bonds"), msg);
      }
    }
  }
  // Excited states, the spectra that name them, and transition densities
  // (spec § 5.2, § 5.4, § 5.12): the sections they name may come later.
  if (rank >= AT_REFERENCES) {
    for (const d of ntoSections) {
      for (const [code, msg] of ntoPairErrors(d.meta, d.members, ntoSections, kindsById, stateIndicesOf)) {
        out.error(code, memberPointer(d.i, "mo_metadata"), msg);
      }
    }
  }
  const excitedRefs: ExcitedStatesIndex | null = rank >= AT_REFERENCES ? { kindsById, moRowsOf } : null;
  for (const d of excitedSections) {
    for (const [code, field, msg] of excitedStatesErrors(d.sid, d.s, d.payload, excitedRefs)) {
      out.error(code, field ? `${sectionPointer(d.i)}/${field}` : memberPointer(d.i, "states"), msg);
    }
  }
  for (const d of statesRefSpectra) {
    for (const msg of spectrumStatesErrors(d.sid, d.payload, kindsById, stateIndicesOf)) {
      out.error("E_REF_UNRESOLVED", memberPointer(d.i, "spectrum"), msg);
    }
  }
  for (const d of transitionDensities) {
    for (const [field, msg] of transitionDensityErrors(d.sid, d.s, kindsById, stateIndicesOf)) {
      out.error("E_REF_UNRESOLVED", `${sectionPointer(d.i)}/${field}`, msg);
    }
  }
  for (const d of wavefunctionVolumes) {
    for (const [field, msg] of wavefunctionVolumeErrors(d.sid, d.kind, d.s, kindsById, moRowsOf)) {
      out.error("E_REF_UNRESOLVED", `${sectionPointer(d.i)}/${field}`, msg);
    }
  }

  // Grids: references and periodic spans against the structure (spec § 5.2).
  if (structures) {
    for (const g of gridMembers) {
      for (const [code, msg] of gridErrors(g.sid, g.s, g.grid, kindsById, cellsOf, structures.onlyStructure)) {
        out.error(code, g.where, msg);
      }
    }
  }
  // properties.*: rules JSON Schema cannot express (spec § 5.13).
  for (const p of propertySections) {
    for (const [code, msg] of propertyErrors(p.sid, p.kind, p.s, p.payload, manifest, structures)) {
      out.error(code, memberPointer(p.i, "properties"), msg);
    }
  }

  // wavefunction.bloch: shapes, and a periodic structure, which may come later.
  const periodic: PeriodicIndex | null = rank >= AT_REFERENCES ? { kindsById, periodicOf } : null;
  for (const b of blochSections) {
    for (const [code, where, msg] of blochErrors(b.sid, b.s, b.members, b.docs, periodic)) {
      out.error(code, `${sectionPointer(b.i)}/${where}`, msg);
    }
  }

  // --- orphans: entries no manifest member references (spec § 2.1) ---
  for (const name of counts.keys()) {
    if (name !== "manifest.json" && !name.endsWith("/") && !declaredPaths.has(name)) {
      out.warning(
        "W_ORPHAN_MEMBER",
        `entry:${name}`,
        `ZIP entry ${quote(name)} is not referenced by any manifest member (orphan; consumers ignore it)`,
      );
    }
  }

  return sections.length;
}

interface PayloadContext {
  zip: ZipArchive;
  entry: ZipEntry;
  m: Obj;
  mpath: string;
  mwhere: string;
  sid: string;
  kind: string;
  role: string;
  lim: Limits;
  rank: number;
  schemaCheck: boolean;
  docs: Map<string, unknown>;
  binaries: Map<string, Uint8Array>;
  out: Findings;
}

/** Read one member and run the payload-level checks on its bytes. */
async function checkPayload(c: PayloadContext): Promise<void> {
  const { zip, entry, m, mpath, mwhere, sid, kind, role, lim, rank, schemaCheck, docs, binaries, out } = c;
  const fmt = typeof m.format === "string" ? m.format : "";
  if (fmt === "json" && entry.size > lim.max_json_bytes) {
    out.error(
      "E_JSON_TOO_LARGE",
      mwhere,
      `member ${quote(mpath)}: ${entry.size} bytes, over the ${lim.max_json_bytes}-byte JSON cap`,
    );
    return;
  }
  let bytes: Uint8Array;
  try {
    bytes = zip.read(entry, { maxBytes: lim.max_member_bytes });
  } catch (err) {
    out.error("E_ZIP_UNREADABLE", `entry:${mpath}`, `ZIP entry ${quote(mpath)} cannot be read: ${errorText(err)}`);
    return;
  }

  if (has(m, "sha256")) {
    const declared = m.sha256;
    if (typeof declared !== "string" || declared !== (await sha256Hex(bytes))) {
      out.error("E_SHA256_MISMATCH", mwhere, `member ${quote(mpath)}: sha256 mismatch`);
    }
  }

  if (fmt === "json") {
    if (jsonDepth(bytes, lim.max_json_depth) > lim.max_json_depth) {
      out.error(
        "E_JSON_TOO_DEEP",
        mwhere,
        `member ${quote(mpath)}: JSON nests deeper than ${lim.max_json_depth} levels`,
      );
      return;
    }
    let doc: unknown;
    let parsed = false;
    try {
      doc = parseJsonBytes(bytes);
      parsed = true;
    } catch (err) {
      out.error("E_JSON_INVALID", mwhere, `member ${quote(mpath)}: not valid JSON (${errorText(err)})`);
    }
    if (parsed) {
      if (rank >= AT_SEMANTICS && kind === "structure" && role === "structure") {
        for (const msg of structurePayloadErrors(sid, doc)) out.error("E_STRUCTURE_PERIODICITY", mwhere, msg);
      }
      // The member's payload definition (spec § 5, § 6).
      const def = payloadDefinition(kind, role);
      if (schemaCheck && def !== undefined) {
        for (const e of checkDefinition(def, doc)) {
          out.error("E_PAYLOAD_SCHEMA", mwhere, `schema: ${mpath}: ${e.location}: ${e.message}`);
        }
      }
      if (rank >= AT_SEMANTICS && CANONICAL_KIND_SET.has(kind)) {
        for (const [loc, unit] of declaredUnits(kind, role, doc)) {
          if (!isUnitName(unit)) {
            out.warning(
              "W_UNIT_UNKNOWN",
              mwhere,
              `member ${quote(mpath)}: declared unit ${show(unit)} at ${loc} is not a unit name of spec § 4.3`,
            );
          }
        }
      }
      docs.set(role, doc);
    }
  } else if (fmt === "binary" && has(m, "dtype")) {
    const dtype = m.dtype;
    const itemsize = isDType(dtype) ? DTYPE_ITEMSIZE[dtype] : 0;
    const shape = m.shape === undefined ? [] : m.shape;
    if (!itemsize) {
      out.error("E_DTYPE", mwhere, `member ${quote(mpath)}: bad dtype ${show(dtype)}`);
    } else if (!isShape(shape)) {
      out.error(
        "E_MEMBER_MALFORMED",
        mwhere,
        `member ${quote(mpath)}: shape must be a list of non-negative integers, got ${show(shape)}`,
      );
    } else {
      const want = expectedByteLength(itemsize, shape);
      if (BigInt(bytes.length) !== want) {
        out.error(
          "E_BINARY_SIZE",
          mwhere,
          `member ${quote(mpath)}: byte length ${bytes.length} != itemsize*prod(shape) ${want}`,
        );
      } else if (kind === "wavefunction.gto" && role === "overlap") {
        binaries.set(role, bytes);
      }
    }
  }

  if (TRAJECTORY_KINDS.has(kind) && role === "energies" && fmt === "binary") binaries.set(role, bytes);

  if (kind === "run.record" && (role === "input" || role === "log") && !isValidUtf8(bytes)) {
    out.error("E_UTF8", mwhere, `member ${quote(mpath)}: run.record ${role} must decode as UTF-8`);
  }
  if (kind === "preview" && fmt === "binary") {
    for (const [code, msg] of previewMemberErrors(sid, role, m, bytes)) out.error(code, mwhere, msg);
  }
}
