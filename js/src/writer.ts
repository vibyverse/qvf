/**
 * `Writer`: build a QVF archive in memory or stream it to a sink (spec § 2,
 * § 3, § 4), in Node and in browsers.
 *
 * Members are written to the ZIP in the order they are added, each compressed
 * and checksummed once (SHA-256 through WebCrypto), and `manifest.json` is
 * written last and STOREd (spec § 2.1). Every JSON payload with a schema
 * definition is checked against it before it is written, and the manifest is
 * checked against the root schema before it is written, so an archive the
 * writer finishes passes the schema stage of the validator.
 *
 * The method names follow the Python `qvf.Writer` (`add_structure` is
 * `addStructure`, keyword `lattice_vectors` is option `latticeVectors`); the
 * payload JSON keeps the snake_case field names of the specification.
 */

import { QvfError } from "./errors.js";
import { CANONICAL_KINDS } from "./kinds.js";
import { type DType, type MemberSpec, type Section, type Source } from "./model.js";
import { checkDefinition, checkManifest, payloadDefinition, type SchemaError } from "./schemaCheck.js";
import { sha256Hex } from "./sha256.js";
import { dtypeOf, encodeTensor, flatten, type NestedNumbers, type TensorInput } from "./tensor.js";
import { isUnitName } from "./units.js";
import { pathProblem, METHOD_DEFLATE, METHOD_STORE } from "./zip.js";
import { ZipWriter, type ZipSink } from "./zipWriter.js";

export const SCHEMA_URI = "https://vibe-qc.org/spec/qvf/1/manifest.schema.json";

const ID_RE = /^[A-Za-z0-9_.-]+$/;
const VENDOR_RE = /^x_[A-Za-z0-9_]+(\.[A-Za-z0-9_.-]+)?$/;
const RESERVED_IDS = new Set(["auto_open", "bookmarks", "camera", "crossfade"]);
const UTC_RE = /^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])T([01][0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9]Z$/;
const UUID_RE = /^[0-9a-f]{8}-[0-9a-f]{4}-[1-8][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/;
const MEDIA_TYPE_RE = /^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}\/[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$/;
const COMPRESSED_MEDIA = new Set(["image/png", "image/jpeg"]);
const ROOT_KEYS = [
  "provenance",
  "thermochemistry",
  "dipole_moment",
  "constraints",
  "energies",
  "extensions",
  "viewer_defaults",
] as const;

/** A root block other than `source`, `metadata` and `sections` (spec § 3.1). */
export type RootBlock = (typeof ROOT_KEYS)[number];

const ELEMENTS = (
  "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr " +
  "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb " +
  "Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf " +
  "Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split(" ");

type Obj = Record<string, unknown>;
const isObject = (v: unknown): v is Obj => typeof v === "object" && v !== null && !Array.isArray(v);

// --- options and inputs --------------------------------------------------------

/** A binary array with an explicit shape, such as `QvfReader.readBinary` returns. */
export interface Tensor {
  data: TensorInput;
  shape: readonly number[];
  dtype?: DType;
}

/** Array data for a builder: nested arrays, a flat typed array (shape derived where it is unambiguous), or a {@link Tensor}. */
export type ArrayInput = TensorInput | Tensor;

export interface WriterOptions {
  /**
   * Add `created_utc` and a random version 4 `uuid` to the root `metadata`
   * unless `setMetadata` supplies them (spec § 3.5). Default `true`.
   */
  autoMetadata?: boolean;
  /**
   * The time of the archive: the automatic `created_utc` and the modification
   * time of every ZIP entry. Default: now with `autoMetadata`, else the DOS
   * epoch 1980-01-01T00:00:00Z. Fix it (and `uuid`) for byte-identical output.
   */
  timestamp?: Date | string;
  /** The automatic `uuid`, instead of a random one. */
  uuid?: string;
  /** `"deflate"` (the default) or `"store"` for the members; `manifest.json` is always stored. */
  compression?: "deflate" | "store";
  /** DEFLATE level 0 to 9, default 6. */
  level?: number;
  /** ZIP64 records for every entry, whatever the sizes (they are written anyway when needed). */
  forceZip64?: boolean;
  /**
   * Stream the archive: every chunk goes to `sink` as soon as it is written
   * and nothing is kept, so memory is bounded by the largest member. Without
   * a sink, `bytes()` and `blob()` return the archive after `finish()`.
   */
  sink?: ZipSink;
  /** Check payloads and the manifest against the JSON Schema before writing them. Default `true`. */
  checkSchema?: boolean;
}

export interface MemberOptions {
  /** The ZIP entry name; default `<section id>/<role>.<ext>`. */
  path?: string;
  /** IANA media type of the member's bytes (spec § 4.2). */
  mediaType?: string;
}

export interface BinaryOptions extends MemberOptions {
  /** The stored dtype; default the typed array's own, `float64` for plain arrays. Values are converted. */
  dtype?: DType;
  /** The shape; default the nesting of a plain array, `[length]` for a typed array. */
  shape?: readonly number[];
}

/** One atom for {@link Writer.addStructure}; positions in Å. */
export interface AtomInput {
  symbol?: string;
  atomic_number?: number;
  position: readonly number[];
  [biomoleculeField: string]: unknown;
}

const ATOM_EXTRA_FIELDS = [
  "atom_name",
  "residue_name",
  "residue_seq",
  "chain_id",
  "b_factor",
  "insertion_code",
  "alt_loc",
  "occupancy",
  "hetero",
];

export interface StructureOptions {
  /** Which axes are periodic; `latticeVectors` is required when one is. */
  pbc?: readonly boolean[];
  /** 3x3 rows in Å, row i paired with `pbc[i]`. */
  latticeVectors?: readonly (readonly number[])[];
  /** Written when given; must equal the number of periodic axes. */
  dimensionality?: number;
  /** `{i, j, order}` pairs, written as the `bonds` member. */
  bonds?: readonly Obj[];
  /** `chains`, `residues`, `secondary_structure`, `b_factors` peer fields. */
  biomolecule?: Obj;
  charge?: number;
  multiplicity?: number;
  name?: string;
  smiles?: string;
  inchi?: string;
  inchikey?: string;
  label?: string;
  /** Default `"structure"`. */
  sectionId?: string;
}

export interface TrajectoryOptions {
  /** Rows in Å: `[3, 3]` or `[n_frames, 3, 3]`. Needs `metadata.pbc`. */
  lattice?: ArrayInput;
  velocities?: ArrayInput;
  gradient?: ArrayInput;
  times?: ArrayInput;
  energies?: ArrayInput;
  cellVolumes?: ArrayInput;
  temperatures?: ArrayInput;
  /** `coords` dtype; also written as `metadata.precision`. */
  precision?: "float64" | "float32";
  /** The metadata `units` declaration (spec § 4.3). */
  units?: Record<string, string>;
  label?: string;
  /** Default `"trajectory"`. */
  sectionId?: string;
}

/** One shell for {@link Writer.addWavefunctionGto}. */
export interface ShellInput {
  center: number;
  l: number;
  exponents: readonly number[];
  coefficients: readonly number[];
  pure?: boolean;
}

export interface WavefunctionGtoOptions {
  /** `[n_mo, n_ao]` (rows are MOs) for a single set of orbitals. */
  moCoefficients?: ArrayInput;
  /** `[n_mo_alpha, n_ao]` and `[n_mo_beta, n_ao]` of an unrestricted wavefunction. */
  moCoefficientsAlpha?: ArrayInput;
  moCoefficientsBeta?: ArrayInput;
  energies?: readonly number[];
  occupations?: readonly number[];
  energiesAlpha?: readonly number[];
  occupationsAlpha?: readonly number[];
  energiesBeta?: readonly number[];
  occupationsBeta?: readonly number[];
  /** Default `"structure"`. */
  structureRef?: string;
  /** Default `"canonical"`. */
  orbitalKind?: string;
  /** Default for shells without `pure`; default `true`. */
  pure?: boolean;
  kPoint?: readonly number[];
  /** Divide the contraction coefficients by the axial primitive norm (spec Appendix A.1). */
  coeffsAreLibintNormalized?: boolean;
  /** Default `"restricted"` with `moCoefficients`, `"unrestricted"` with alpha and beta. */
  spin?: "restricted" | "unrestricted" | "restricted_open" | "generalized";
  ecp?: readonly Obj[];
  nElectrons?: number;
  /** One label per MO; `{alpha, beta}` for an unrestricted wavefunction. */
  symmetryLabels?: readonly string[] | { alpha?: readonly string[]; beta?: readonly string[] };
  /** `[n_ao, n_ao]`, the AO overlap of the AOs as written; producers should pass it. */
  overlap?: ArrayInput;
  densityMatrix?: ArrayInput;
  densityMatrixAlpha?: ArrayInput;
  densityMatrixBeta?: ArrayInput;
  /** Coefficients are `[..., 2]` = `[re, im]` pairs. */
  complexCoefficients?: boolean;
  stateRef?: string;
  stateIndex?: number;
  ntoRole?: "hole" | "particle";
  localizationReferenceBasis?: string;
  label?: string;
  /** Default `"wf"`. */
  sectionId?: string;
}

export interface VolumeOptions {
  /** Stored dtype of `data`, default `float32`. */
  dtype?: "float32" | "float64";
  label?: string;
  component?: string;
  /** Default `<kind with _>_<section count>`. */
  sectionId?: string;
  /** The two operands of a `volume.difference`. */
  operandA?: string;
  operandB?: string;
  description?: string;
  /** The grid `units` declaration (spec § 4.3). */
  units?: Record<string, string>;
  periodic?: readonly boolean[];
  includesEndpoint?: boolean;
  structureRef?: string;
  /** The `axis` of a `volume.series`: `{kind, values, units?, labels?, trajectory_ref?}`. */
  axis?: Obj;
  excitedStatesRef?: string;
  stateIndex?: number;
}

export interface ExcitedStatesOptions {
  meta?: Obj;
  wavefunctionRef?: string;
  label?: string;
  /** Default `"excited_states"`. */
  sectionId?: string;
}

export interface SpectrumOptions {
  /** x values of a 1-D spectrum, in the kind's unit unless `xUnits` declares another. */
  frequencies?: readonly number[];
  intensities?: readonly number[];
  /** The whole payload instead (the only form for `spectra.nmr` and `spectra.epr`). */
  payload?: Obj;
  xUnits?: string;
  yUnits?: string;
  peaks?: readonly Obj[];
  broadening?: Obj;
  label?: string;
  /** Default the kind with `.` replaced by `_`. */
  sectionId?: string;
}

/** A participant of a `reaction.scheme` step: a structure section id or `{structure_ref, coefficient?}`. */
export type Participant = string | { structure_ref: string; coefficient?: number };

/** One step of {@link Writer.addReactionScheme}. */
export interface SchemeStep {
  id: string;
  reactants: readonly Participant[];
  products: readonly Participant[];
  reagents?: readonly Participant[];
  label?: string;
  arrow?: "forward" | "equilibrium" | "retrosynthetic" | "resonance";
  conditions?: readonly string[];
  yield_percent?: number;
  yield_text?: string;
  energies?: Obj;
  [vendorKey: `x_${string}`]: unknown;
}

export interface ReactionSchemeOptions {
  label?: string;
  /** Default `"reaction_scheme"`. */
  sectionId?: string;
}

export interface AnnotateOptions {
  structureRef?: string;
  runRef?: string;
  group?: string;
  label?: string;
}

// --- helpers --------------------------------------------------------------

const encoder = new TextEncoder();

/**
 * UTF-8 JSON bytes of a member (2-space indentation, as the Python writer
 * writes). Typed arrays become arrays; NaN, infinities, bigints, functions and
 * symbols are refused.
 */
export function jsonBytes(value: unknown, what = "payload"): Uint8Array {
  const text = JSON.stringify(
    value,
    function (this: unknown, key: string, v: unknown): unknown {
      if (typeof v === "number" && !Number.isFinite(v)) {
        throw new RangeError(`${what}: ${key || "value"} is ${v}; JSON has no NaN or infinity`);
      }
      if (typeof v === "bigint" || typeof v === "function" || typeof v === "symbol") {
        throw new TypeError(`${what}: ${key || "value"} is a ${typeof v}, not JSON`);
      }
      if (ArrayBuffer.isView(v) && !(v instanceof DataView)) {
        const own = dtypeOf(v);
        if (own === "int64" || own === "uint64") {
          return Array.from(v as BigInt64Array, (b) => {
            const n = Number(b);
            if (!Number.isSafeInteger(n)) throw new RangeError(`${what}: ${key} holds ${b}, beyond 2^53`);
            return n;
          });
        }
        return Array.from(v as unknown as ArrayLike<number>);
      }
      return v;
    },
    2,
  );
  if (text === undefined) throw new TypeError(`${what} is not JSON`);
  return encoder.encode(text);
}

function schemaMessage(errors: readonly SchemaError[], max = 5): string {
  const shown = errors.slice(0, max).map((e) => `${e.pointer || "/"}: ${e.message}`);
  if (errors.length > max) shown.push(`... and ${errors.length - max} more`);
  return shown.join("; ");
}

function utcString(value: Date | string, what: string): string {
  if (value instanceof Date) {
    if (Number.isNaN(value.getTime())) throw new RangeError(`${what} is an invalid date`);
    return value.toISOString().replace(/\.\d{3}Z$/, "Z");
  }
  if (typeof value !== "string" || !UTC_RE.test(value)) {
    throw new RangeError(`${what} must be YYYY-MM-DDThh:mm:ssZ (spec § 3.5), got ${JSON.stringify(value)}`);
  }
  return value;
}

function randomUuid(): string {
  const c = globalThis.crypto;
  if (c && typeof c.randomUUID === "function") return c.randomUUID();
  if (!c || typeof c.getRandomValues !== "function") throw new Error("no WebCrypto random source for the uuid");
  const b = c.getRandomValues(new Uint8Array(16));
  b[6] = (b[6]! & 0x0f) | 0x40;
  b[8] = (b[8]! & 0x3f) | 0x80;
  const h = Array.from(b, (x) => x.toString(16).padStart(2, "0")).join("");
  return `${h.slice(0, 8)}-${h.slice(8, 12)}-${h.slice(12, 16)}-${h.slice(16, 20)}-${h.slice(20)}`;
}

function finite(value: unknown, what: string): number {
  if (typeof value !== "number" || !Number.isFinite(value)) {
    throw new RangeError(`${what} must be a finite number, got ${JSON.stringify(value) ?? String(value)}`);
  }
  return value;
}

function integer(value: unknown, what: string, minimum = -Infinity): number {
  if (typeof value !== "number" || !Number.isInteger(value)) {
    throw new RangeError(`${what} must be an integer, got ${JSON.stringify(value) ?? String(value)}`);
  }
  if (value < minimum) throw new RangeError(`${what} must be >= ${minimum}, got ${value}`);
  return value;
}

function numbers(values: readonly unknown[] | undefined, what: string, n?: number): number[] {
  if (values === undefined) return new Array<number>(n ?? 0).fill(0);
  if (!Array.isArray(values) && !ArrayBuffer.isView(values)) throw new TypeError(`${what} must be an array of numbers`);
  const out = Array.from(values as ArrayLike<unknown>, (v, i) => finite(v, `${what}[${i}]`));
  if (n !== undefined && out.length !== n) throw new RangeError(`${what} has ${out.length} entries, ${n} expected`);
  return out;
}

function vector(value: unknown, n: number, what: string): number[] {
  if (!Array.isArray(value) && !ArrayBuffer.isView(value)) throw new TypeError(`${what} must be ${n} numbers`);
  return numbers(Array.from(value as ArrayLike<unknown>), what, n);
}

function matrix3(value: unknown, what: string): number[][] {
  if (!Array.isArray(value) || value.length !== 3) throw new RangeError(`${what} must be a 3x3 array`);
  return value.map((row, i) => vector(row, 3, `${what}[${i}]`));
}

function unitName(value: unknown, what: string): string {
  if (!isUnitName(value)) throw new RangeError(`${what}: ${JSON.stringify(value)} is not a QVF unit name (spec § 4.3)`);
  return value;
}

function units(value: Record<string, string>, what: string): Record<string, string> {
  if (!isObject(value)) throw new TypeError(`${what}: units must map names to unit names`);
  const out: Record<string, string> = {};
  for (const [k, v] of Object.entries(value)) out[k] = unitName(v, `${what} units.${k}`);
  return out;
}

function isTensor(value: unknown): value is Tensor {
  return isObject(value) && "data" in value && Array.isArray(value.shape);
}

/** Split an {@link ArrayInput} into data and shape, giving a flat typed array the shape `[..., trailing]` where that is unambiguous. */
function resolveArray(input: ArrayInput, what: string, trailing: readonly number[] = []): Tensor {
  if (isTensor(input)) return input;
  const own = dtypeOf(input);
  if (own !== undefined) {
    const n = (input as ArrayLike<unknown>).length;
    const inner = trailing.reduce((a, b) => a * b, 1);
    if (trailing.length === 0) return { data: input, shape: [n] };
    if (inner === 0 || n % inner !== 0) {
      throw new RangeError(`${what}: ${n} elements do not fill rows of shape ${JSON.stringify(trailing)}`);
    }
    return { data: input, shape: [n / inner, ...trailing] };
  }
  return { data: input, shape: flatten(input as NestedNumbers, what).shape };
}

function sameShape(shape: readonly number[], want: readonly number[], what: string, form: string): void {
  if (shape.length !== want.length || shape.some((d, i) => d !== want[i])) {
    throw new RangeError(`${what} must be ${form} = ${JSON.stringify(want)}, got shape ${JSON.stringify(shape)}`);
  }
}

/** Float64 values of an array input, flat (for checks that read the numbers). */
function float64Values(t: Tensor, what: string): Float64Array {
  const { bytes } = encodeTensor(t.data, { dtype: "float64", shape: t.shape, what });
  return new Float64Array(bytes.buffer, bytes.byteOffset, bytes.length / 8);
}

/** Norm of the axial Cartesian Gaussian primitive of total angular momentum l (spec Appendix A.1). */
export function primitiveNorm(alpha: number, l: number): number {
  let df = 1;
  for (let k = 1; k < 2 * l; k += 2) df *= k;
  return ((2 * alpha) / Math.PI) ** 0.75 * (4 * alpha) ** (l / 2) / Math.sqrt(df);
}

// --- the writer ------------------------------------------------------------

/**
 * Builds one QVF archive. Every `add*` method returns a promise; members are
 * written in call order even when calls are not awaited one by one, and a
 * failed member write makes `finish()` refuse.
 */
export class Writer {
  /** Whether `created_utc` and `uuid` are added automatically (see {@link WriterOptions.autoMetadata}). */
  readonly autoMetadata: boolean;

  private readonly source: Source;
  private readonly zip: ZipWriter;
  private readonly method: number;
  private readonly level: number;
  private readonly schema: boolean;
  private readonly autoValues: { created_utc: string; uuid: string } | null;
  private readonly sections: Section[] = [];
  private readonly ids = new Set<string>();
  private readonly paths = new Set<string>(["manifest.json"]);
  private readonly root: Obj = {};
  /** Members queued but not yet written, as `<section id>\0<role>`. */
  private readonly pendingRoles = new Set<string>();
  private metadata: Obj = {};
  private chain: Promise<void> = Promise.resolve();
  private failure: unknown = null;
  private state: "open" | "finishing" | "finished" = "open";

  constructor(source: { program: string; version: string; calculation?: string }, options: WriterOptions = {}) {
    if (!isObject(source) || typeof source.program !== "string" || !source.program) {
      throw new QvfError("E_SOURCE", "source.program must be a non-empty string");
    }
    if (typeof source.version !== "string" || !source.version) {
      throw new QvfError("E_SOURCE", "source.version must be a non-empty string");
    }
    if (source.calculation !== undefined && typeof source.calculation !== "string") {
      throw new QvfError("E_SOURCE", "source.calculation must be a string");
    }
    this.source = { program: source.program, version: source.version, calculation: source.calculation ?? "" };
    this.autoMetadata = options.autoMetadata ?? true;
    let stamp: string | null = options.timestamp !== undefined ? utcString(options.timestamp, "timestamp") : null;
    if (this.autoMetadata) {
      stamp ??= utcString(new Date(), "timestamp");
      const uuid = options.uuid ?? randomUuid();
      if (!UUID_RE.test(uuid)) throw new RangeError(`uuid ${JSON.stringify(uuid)} is not a lowercase RFC 4122 UUID`);
      this.autoValues = { created_utc: stamp, uuid };
    } else {
      this.autoValues = null;
    }
    if (options.compression !== undefined && options.compression !== "deflate" && options.compression !== "store") {
      throw new RangeError(`compression must be "deflate" or "store"`);
    }
    this.method = options.compression === "store" ? METHOD_STORE : METHOD_DEFLATE;
    this.level = integer(options.level ?? 6, "level", 0);
    if (this.level > 9) throw new RangeError("level must be 0 to 9");
    this.schema = options.checkSchema ?? true;
    const zipOptions: ConstructorParameters<typeof ZipWriter>[0] = {
      forceZip64: options.forceZip64 ?? false,
      modified: stamp === null ? new Date(Date.UTC(1980, 0, 1)) : new Date(stamp),
    };
    if (options.sink) zipOptions.sink = options.sink;
    this.zip = new ZipWriter(zipOptions);
  }

  // -- root blocks -------------------------------------------------------------

  /**
   * Merge fields into the root `metadata` block (spec § 3.5); `null` removes a
   * key. `created_utc` / `modified_utc` take a `Date` or a UTC string, `uuid`
   * is lowercased. Fields set here win over the automatic ones.
   */
  setMetadata(fields: Obj): void {
    this.ensureOpen();
    const next = { ...this.metadata };
    for (const [key, value] of Object.entries(fields)) {
      if (value === null || value === undefined) {
        delete next[key];
      } else if (key === "created_utc" || key === "modified_utc") {
        next[key] = utcString(value as Date | string, `metadata ${key}`);
      } else if (key === "uuid") {
        const text = typeof value === "string" ? value.toLowerCase() : value;
        if (typeof text !== "string" || !UUID_RE.test(text)) {
          throw new RangeError(`metadata uuid must be an RFC 4122 UUID, got ${JSON.stringify(value)}`);
        }
        next[key] = text;
      } else {
        next[key] = value;
      }
    }
    this.checkPayload("Metadata", next, "metadata");
    this.metadata = next;
  }

  /** Merge fields into the root `provenance` block (spec § 3.2). */
  setProvenance(fields: Obj): void {
    this.setRoot("provenance", { ...((this.root.provenance as Obj | undefined) ?? {}), ...fields });
  }

  /**
   * Set a root block (spec § 3.1): `provenance`, `thermochemistry`,
   * `dipole_moment`, `constraints`, `energies`, `extensions` or
   * `viewer_defaults`; `null` removes it. The whole manifest is checked
   * against the schema in `finish()`.
   */
  setRoot(key: RootBlock, value: unknown): void {
    this.ensureOpen();
    if (!(ROOT_KEYS as readonly string[]).includes(key)) throw new RangeError(`${JSON.stringify(key)} is not a settable root block`);
    if (value === null || value === undefined) delete this.root[key];
    else this.root[key] = JSON.parse(new TextDecoder().decode(jsonBytes(value, key))) as unknown;
  }

  /** Root `energies` block (spec § 3.6). */
  setEnergies(block: Obj): void {
    this.setRoot("energies", block);
  }

  /** Root `viewer_defaults` block (spec § 3.4). */
  setViewerDefaults(block: Obj): void {
    this.setRoot("viewer_defaults", block);
  }

  // -- generic sections and members -------------------------------------------

  /** The ids of the sections added so far, in order. */
  get sectionIds(): string[] {
    return this.sections.map((s) => s.id);
  }

  /**
   * Add a section: `{id, kind, label?, ...peer fields}` with no members yet;
   * attach them with `addJson`, `addBinary` and `addOpaque`. Returns the id.
   */
  addSection(section: { id: string; kind: string; [field: string]: unknown }): string {
    this.ensureOpen();
    if (!isObject(section)) throw new TypeError("section must be an object");
    const { id, kind, members, ...rest } = section;
    if (members !== undefined && !(isObject(members) && Object.keys(members).length === 0)) {
      throw new RangeError("add members through addJson, addBinary or addOpaque, not in addSection");
    }
    this.checkKind(kind);
    this.reserveId(id);
    const value = { id, kind, ...JSON.parse(new TextDecoder().decode(jsonBytes(rest, `section ${id}`))), members: {} } as Section;
    this.sections.push(value);
    return id;
  }

  /** Set the generic fields of an added section (spec § 4.1); references must name added sections of the right kind. */
  annotateSection(sectionId: string, options: AnnotateOptions): string {
    this.ensureOpen();
    const section = this.sectionOrThrow(sectionId);
    if (options.structureRef !== undefined) {
      this.checkRef("structure_ref", options.structureRef, ["structure"]);
      section.structure_ref = options.structureRef;
    }
    if (options.runRef !== undefined) {
      this.checkRef("run_ref", options.runRef, ["run.record"]);
      section.run_ref = options.runRef;
    }
    if (options.group !== undefined) section.group = String(options.group);
    if (options.label !== undefined) section.label = String(options.label);
    return sectionId;
  }

  /** Write a JSON member; a payload with a schema definition is checked first (`E_PAYLOAD_SCHEMA`). */
  async addJson(sectionId: string, role: string, payload: unknown, options: MemberOptions = {}): Promise<MemberSpec> {
    const what = `section ${JSON.stringify(sectionId)} member ${JSON.stringify(role)}`;
    const section = this.sectionOrThrow(sectionId);
    const extra = this.mediaType(options.mediaType);
    const bytes = jsonBytes(payload, what);
    const def = payloadDefinition(section.kind, role);
    if (def !== undefined) this.checkPayload(def, JSON.parse(new TextDecoder().decode(bytes)), what, "E_PAYLOAD_SCHEMA");
    const path = this.prepareMember(section, role, options.path ?? `${sectionId}/${role}.json`);
    return this.write(section, role, path, bytes, this.method, (sha256) => ({
      path,
      format: "json",
      sha256,
      size: bytes.length,
      ...extra,
    }));
  }

  /**
   * Write a binary array member: little-endian, C-ordered, of `dtype` and
   * `shape` (spec § 4.2). Values are converted to the dtype; an integer dtype
   * refuses values it cannot hold.
   */
  async addBinary(sectionId: string, role: string, data: TensorInput, options: BinaryOptions = {}): Promise<MemberSpec> {
    const section = this.sectionOrThrow(sectionId);
    const extra = this.mediaType(options.mediaType);
    const opts: Parameters<typeof encodeTensor>[1] = { what: `section ${JSON.stringify(sectionId)} member ${JSON.stringify(role)}` };
    if (options.dtype !== undefined) opts.dtype = options.dtype;
    if (options.shape !== undefined) opts.shape = options.shape;
    const { bytes, dtype, shape } = encodeTensor(data, opts);
    const path = this.prepareMember(section, role, options.path ?? `${sectionId}/${role}.bin`);
    return this.write(section, role, path, bytes, this.method, (sha256) => ({
      path,
      format: "binary",
      dtype,
      shape,
      sha256,
      size: bytes.length,
      ...extra,
    }));
  }

  /**
   * Write an opaque binary member (no dtype or shape, spec § 4.2): text or
   * bytes such as a log, a BibTeX file or a PNG image. PNG and JPEG are
   * stored uncompressed.
   */
  async addOpaque(sectionId: string, role: string, data: Uint8Array | string, options: MemberOptions = {}): Promise<MemberSpec> {
    const section = this.sectionOrThrow(sectionId);
    const extra = this.mediaType(options.mediaType);
    const bytes = typeof data === "string" ? encoder.encode(data) : data;
    if (!(bytes instanceof Uint8Array)) throw new TypeError("opaque member data must be a Uint8Array or a string");
    const path = this.prepareMember(section, role, options.path ?? `${sectionId}/${role}.bin`);
    const method = extra.media_type !== undefined && COMPRESSED_MEDIA.has(extra.media_type) ? METHOD_STORE : this.method;
    return this.write(section, role, path, bytes, method, (sha256) => ({
      path,
      format: "binary",
      sha256,
      size: bytes.length,
      ...extra,
    }));
  }

  // -- typed builders (the names of the Python writer) -----------------------------

  /** `add_structure`: a `structure` section (spec § 5.1); positions in Å. */
  async addStructure(atoms: readonly AtomInput[], options: StructureOptions = {}): Promise<string> {
    const sid = options.sectionId ?? "structure";
    this.checkNewId(sid);
    if (!Array.isArray(atoms)) throw new TypeError("atoms must be an array");
    const outAtoms = atoms.map((raw, k) => {
      if (!isObject(raw)) throw new TypeError(`atoms[${k}] must be an object`);
      const a = raw as AtomInput;
      let z: number | undefined = a.atomic_number;
      let sym: string | undefined = a.symbol;
      if (z === undefined && sym !== undefined) z = ELEMENTS.indexOf(sym) >= 0 ? ELEMENTS.indexOf(sym) : 0;
      if (sym === undefined && z !== undefined) sym = ELEMENTS[integer(z, `atoms[${k}].atomic_number`, 0)] ?? "X";
      const atom: Obj = {
        symbol: sym,
        position: vector(a.position, 3, `atoms[${k}].position`),
        atomic_number: z === undefined ? 0 : integer(z, `atoms[${k}].atomic_number`, 0),
      };
      for (const f of ATOM_EXTRA_FIELDS) if (a[f] !== undefined && a[f] !== null) atom[f] = a[f];
      return atom;
    });
    let flags: boolean[];
    if (options.pbc !== undefined) {
      if (!Array.isArray(options.pbc) || options.pbc.length !== 3 || !options.pbc.every((x) => typeof x === "boolean")) {
        throw new RangeError("pbc must be three booleans");
      }
      flags = [...options.pbc];
    } else if (options.dimensionality !== undefined) {
      const d = integer(options.dimensionality, "dimensionality", 0);
      flags = [0, 1, 2].map((i) => i < d);
    } else {
      flags = [false, false, false];
    }
    const nPeriodic = flags.filter(Boolean).length;
    if (options.dimensionality !== undefined && options.dimensionality !== nPeriodic) {
      throw new QvfError(
        "E_STRUCTURE_PERIODICITY",
        `dimensionality ${options.dimensionality} != sum(pbc) ${nPeriodic} (spec § 5.1)`,
      );
    }
    const lattice = options.latticeVectors !== undefined ? matrix3(options.latticeVectors, "latticeVectors") : null;
    if (nPeriodic > 0 && lattice === null) {
      throw new QvfError("E_STRUCTURE_PERIODICITY", `pbc ${JSON.stringify(flags)} marks a periodic axis, so latticeVectors is required (spec § 5.1)`);
    }
    const payload: Obj = { atoms: outAtoms, pbc: flags, lattice_vectors: lattice };
    if (options.dimensionality !== undefined) payload.dimensionality = options.dimensionality;

    const identity: Obj = {};
    if (options.charge !== undefined) identity.charge = integer(options.charge, "charge");
    if (options.multiplicity !== undefined) identity.multiplicity = integer(options.multiplicity, "multiplicity", 1);
    for (const key of ["name", "smiles", "inchi", "inchikey"] as const) {
      const v = options[key];
      if (v !== undefined) {
        if (typeof v !== "string" || !v) throw new RangeError(`${key} must be a non-empty string`);
        identity[key] = v;
      }
    }
    if (identity.charge !== undefined) {
      const sumZ = outAtoms.reduce((n, a) => n + (a.atomic_number as number), 0);
      const electrons = sumZ - (identity.charge as number);
      const m = identity.multiplicity as number | undefined;
      if (electrons < 0 || (m !== undefined && (m - 1 > electrons || (electrons - (m - 1)) % 2 !== 0))) {
        throw new QvfError(
          "E_CHARGE_MULTIPLICITY",
          `structure ${JSON.stringify(sid)}: charge ${String(identity.charge)} and multiplicity ${String(m)} are impossible for ${sumZ} protons (spec § 5.1)`,
        );
      }
    }
    const peers: Obj = {};
    if (options.label !== undefined) peers.label = options.label;
    if (options.biomolecule !== undefined) {
      for (const key of ["chains", "residues", "secondary_structure", "b_factors"]) {
        if (key in options.biomolecule) peers[key] = options.biomolecule[key];
      }
    }
    Object.assign(peers, identity);
    this.precheck("structure", "structure", payload, sid);
    const bonds = options.bonds !== undefined ? { pairs: options.bonds.map((b) => ({ ...b })) } : undefined;
    if (bonds) this.precheck("structure", "bonds", bonds, sid);
    this.addSection({ id: sid, kind: "structure", ...peers });
    const jobs = [this.addJson(sid, "structure", payload, { path: `${sid}/structure.json` })];
    if (bonds) jobs.push(this.addJson(sid, "bonds", bonds, { path: `${sid}/bonds.json` }));
    await Promise.all(jobs);
    return sid;
  }

  /** `add_trajectory`: a `trajectory` section (spec § 5.5); `coords` is `[n_frames, n_atoms, 3]` in Å. */
  async addTrajectory(metadata: Obj, coords: ArrayInput, options: TrajectoryOptions = {}): Promise<string> {
    const sid = options.sectionId ?? "trajectory";
    this.checkNewId(sid);
    if (!isObject(metadata)) throw new TypeError("trajectory metadata must be an object");
    const meta: Obj = { ...metadata };
    const atomCount = this.atomCount(meta);
    const c = resolveArray(coords, "trajectory coords", atomCount ? [atomCount, 3] : []);
    if (c.shape.length !== 3 || c.shape[2] !== 3) {
      throw new RangeError(`trajectory coords must be [n_frames, n_atoms, 3], got shape ${JSON.stringify(c.shape)}`);
    }
    const [nFrames, nAtoms] = c.shape as [number, number, number];
    const declared = meta.precision;
    for (const [v, what] of [
      [options.precision, "precision"],
      [declared, "metadata precision"],
    ] as const) {
      if (v !== undefined && v !== "float64" && v !== "float32") throw new RangeError(`trajectory ${what} must be "float64" or "float32"`);
    }
    if (options.precision !== undefined && declared !== undefined && declared !== options.precision) {
      throw new RangeError(`trajectory precision ${options.precision} contradicts metadata precision ${String(declared)}`);
    }
    if (options.precision !== undefined) meta.precision = options.precision;
    const coordsDtype = (meta.precision as DType | undefined) ?? "float64";
    const pbc = meta.pbc;
    if (pbc !== undefined && !(Array.isArray(pbc) && pbc.length === 3 && pbc.every((x) => typeof x === "boolean"))) {
      throw new RangeError("trajectory metadata pbc must be three booleans");
    }
    let lattice: Tensor | undefined;
    if (options.lattice !== undefined) {
      lattice = resolveArray(options.lattice, "trajectory lattice", [3, 3]);
      const s = lattice.shape;
      const ok = (s.length === 2 && s[0] === 3 && s[1] === 3) || (s.length === 3 && s[0] === nFrames && s[1] === 3 && s[2] === 3);
      if (!ok) throw new RangeError(`trajectory lattice must be [3, 3] or [${nFrames}, 3, 3], got ${JSON.stringify(s)}`);
      if (pbc === undefined) throw new QvfError("E_STRUCTURE_PERIODICITY", "trajectory lattice needs metadata.pbc (spec § 5.5)");
    } else if (Array.isArray(pbc) && pbc.some(Boolean)) {
      throw new QvfError("E_STRUCTURE_PERIODICITY", `trajectory metadata pbc ${JSON.stringify(pbc)} marks a periodic axis, so lattice is required (spec § 5.5)`);
    }
    const frames: [string, Tensor][] = [];
    const frameInputs: [string, ArrayInput | undefined, boolean][] = [
      ["velocities", options.velocities, true],
      ["gradient", options.gradient, true],
      ["times", options.times, false],
      ["energies", options.energies, false],
      ["cell_volumes", options.cellVolumes, false],
      ["temperatures", options.temperatures, false],
    ];
    for (const [role, input, isVector] of frameInputs) {
      if (input === undefined) continue;
      const t = resolveArray(input, `trajectory ${role}`, isVector ? [nAtoms, 3] : []);
      sameShape(t.shape, isVector ? [nFrames, nAtoms, 3] : [nFrames], `trajectory ${role}`, isVector ? "[n_frames, n_atoms, 3]" : "[n_frames]");
      frames.push([role, t]);
    }
    const binaryEnergies = frames.find(([r]) => r === "energies");
    if (binaryEnergies && meta.energies !== undefined) {
      const listed = numbers(meta.energies as number[], "trajectory metadata energies");
      const values = float64Values(binaryEnergies[1], "trajectory energies");
      if (listed.length !== values.length || listed.some((v, i) => Math.abs(v - values[i]!) > 1e-10)) {
        throw new RangeError("trajectory energies and metadata energies must agree within 1e-10 hartree (spec § 5.5)");
      }
    }
    if (options.units !== undefined) {
      if (meta.units !== undefined) throw new RangeError("trajectory units given both as an option and in metadata");
      meta.units = units(options.units, "trajectory");
    }
    this.precheck("trajectory", "metadata", meta, sid);
    this.addSection({ id: sid, kind: "trajectory", ...(options.label !== undefined ? { label: options.label } : {}) });
    const jobs: Promise<MemberSpec>[] = [
      this.addJson(sid, "metadata", meta, { path: `${sid}/metadata.json` }),
      this.addBinary(sid, "coords", c.data, { dtype: coordsDtype, shape: c.shape, path: `${sid}/coords.bin` }),
    ];
    if (lattice) jobs.push(this.addBinary(sid, "lattice", lattice.data, { dtype: "float64", shape: lattice.shape, path: `${sid}/lattice.bin` }));
    for (const [role, t] of frames) {
      jobs.push(this.addBinary(sid, role, t.data, { dtype: "float64", shape: t.shape, path: `${sid}/${role}.bin` }));
    }
    await Promise.all(jobs);
    return sid;
  }

  /**
   * `add_wavefunction_gto`: a `wavefunction.gto` section (spec § 5.3,
   * Appendix A), restricted (`moCoefficients`) or unrestricted
   * (`moCoefficientsAlpha` + `moCoefficientsBeta`). Coefficients are float64
   * `[n_mo, n_ao]`, or `[n_mo, n_ao, 2]` with `complexCoefficients`.
   */
  async addWavefunctionGto(shells: readonly ShellInput[], options: WavefunctionGtoOptions = {}): Promise<string> {
    const sid = options.sectionId ?? "wf";
    this.checkNewId(sid);
    const pure = options.pure ?? true;
    if (!Array.isArray(shells)) throw new TypeError("shells must be an array");
    const outShells = shells.map((raw, k) => {
      if (!isObject(raw)) throw new TypeError(`shells[${k}] must be an object`);
      const sh = raw as unknown as ShellInput;
      const l = integer(sh.l, `shells[${k}].l`, 0);
      const exps = numbers(sh.exponents, `shells[${k}].exponents`);
      let coeffs = numbers(sh.coefficients, `shells[${k}].coefficients`, exps.length);
      if (options.coeffsAreLibintNormalized) coeffs = coeffs.map((cf, i) => cf / primitiveNorm(exps[i]!, l));
      return {
        center: integer(sh.center, `shells[${k}].center`, 0),
        l,
        exponents: exps,
        coefficients: coeffs,
        pure: typeof sh.pure === "boolean" ? sh.pure : pure,
      };
    });
    const nAo = outShells.reduce((n, sh) => n + (sh.pure ? 2 * sh.l + 1 : ((sh.l + 1) * (sh.l + 2)) / 2), 0);

    const haveAlphaBeta = options.moCoefficientsAlpha !== undefined && options.moCoefficientsBeta !== undefined;
    if (options.moCoefficients === undefined && !haveAlphaBeta) {
      throw new RangeError("pass moCoefficients, or both moCoefficientsAlpha and moCoefficientsBeta");
    }
    const spin = options.spin ?? (options.moCoefficients !== undefined ? "restricted" : "unrestricted");
    if (!["restricted", "unrestricted", "restricted_open", "generalized"].includes(spin)) {
      throw new RangeError(`spin ${JSON.stringify(spin)} is not restricted, unrestricted, restricted_open or generalized`);
    }
    const unrestricted = spin === "unrestricted";
    if (unrestricted && !haveAlphaBeta) throw new RangeError("an unrestricted wavefunction needs moCoefficientsAlpha and moCoefficientsBeta");
    if (!unrestricted && options.moCoefficients === undefined) throw new RangeError(`a ${spin} wavefunction needs moCoefficients`);
    const complex = options.complexCoefficients ?? false;
    if (spin === "generalized" && !complex) throw new RangeError("a generalized wavefunction needs complexCoefficients (spec § 5.3)");
    const aoExtent = spin === "generalized" ? 2 * nAo : nAo;
    const coefficients = (input: ArrayInput, what: string): Tensor => {
      const t = resolveArray(input, what, complex ? [aoExtent, 2] : [aoExtent]);
      const want = complex ? 3 : 2;
      if (t.shape.length !== want || t.shape[1] !== aoExtent || (complex && t.shape[2] !== 2)) {
        throw new RangeError(
          `${what} must be [n_mo, ${aoExtent}${complex ? ", 2" : ""}] for ${nAo} AOs, got shape ${JSON.stringify(t.shape)}`,
        );
      }
      return t;
    };
    const labels = (v: unknown, what: string, n: number): string[] => {
      if (!Array.isArray(v) || v.length !== n || !v.every((x) => typeof x === "string")) {
        throw new RangeError(`${what} must be ${n} strings, one per MO`);
      }
      return [...(v as string[])];
    };

    const matrices: [string, Tensor][] = [];
    let moMeta: Obj;
    let occupationLists: number[][];
    const coefficientMembers: [string, Tensor][] = [];
    if (!unrestricted) {
      const C = coefficients(options.moCoefficients!, "moCoefficients");
      const nMo = C.shape[0]!;
      moMeta = {
        n_mo: nMo,
        n_ao: nAo,
        spin,
        orbital_kind: options.orbitalKind ?? "canonical",
        energies: numbers(options.energies, "energies", nMo),
        occupations: numbers(options.occupations, "occupations", nMo),
      };
      if (options.symmetryLabels !== undefined) moMeta.symmetry_labels = labels(options.symmetryLabels, "symmetryLabels", nMo);
      coefficientMembers.push(["mo_coefficients", C]);
      occupationLists = [moMeta.occupations as number[]];
    } else {
      const Ca = coefficients(options.moCoefficientsAlpha!, "moCoefficientsAlpha");
      const Cb = coefficients(options.moCoefficientsBeta!, "moCoefficientsBeta");
      const channel = (n: number, e: readonly number[] | undefined, o: readonly number[] | undefined, name: string): Obj => ({
        n_mo: n,
        energies: numbers(e, `energies${name}`, n),
        occupations: numbers(o, `occupations${name}`, n),
      });
      const alpha = channel(Ca.shape[0]!, options.energiesAlpha, options.occupationsAlpha, "Alpha");
      const beta = channel(Cb.shape[0]!, options.energiesBeta, options.occupationsBeta, "Beta");
      const sl = options.symmetryLabels;
      if (sl !== undefined) {
        if (Array.isArray(sl) || !isObject(sl)) throw new RangeError("symmetryLabels of an unrestricted wavefunction must be {alpha, beta}");
        const map = sl as { alpha?: readonly string[]; beta?: readonly string[] };
        if (map.alpha !== undefined) alpha.symmetry_labels = labels(map.alpha, "symmetryLabels.alpha", alpha.n_mo as number);
        if (map.beta !== undefined) beta.symmetry_labels = labels(map.beta, "symmetryLabels.beta", beta.n_mo as number);
      }
      moMeta = { n_ao: nAo, spin: "unrestricted", orbital_kind: options.orbitalKind ?? "canonical", alpha, beta };
      coefficientMembers.push(["mo_coefficients_alpha", Ca], ["mo_coefficients_beta", Cb]);
      occupationLists = [alpha.occupations as number[], beta.occupations as number[]];
    }
    if (options.kPoint !== undefined) moMeta.k_point = vector(options.kPoint, 3, "kPoint");
    if (complex) moMeta.complex = true;
    if (options.localizationReferenceBasis !== undefined) moMeta.localization_reference_basis = options.localizationReferenceBasis;
    if (options.stateRef !== undefined) moMeta.state_ref = options.stateRef;
    if (options.stateIndex !== undefined) {
      if (options.stateRef === undefined) throw new RangeError("stateIndex needs stateRef (spec § 5.3)");
      moMeta.state_index = integer(options.stateIndex, "stateIndex", 1);
    }
    if (options.ntoRole !== undefined) {
      if (options.orbitalKind !== "nto" || !["hole", "particle"].includes(options.ntoRole)) {
        throw new RangeError('ntoRole is "hole" or "particle", for orbitalKind "nto" (spec § 5.3)');
      }
      moMeta.nto_role = options.ntoRole;
    }
    if (options.orbitalKind === "nto") moMeta.occupation_semantics = "transition_weight";

    const basis: Obj = { structure_ref: options.structureRef ?? "structure", pure, n_ao: nAo, shells: outShells };
    if (options.ecp !== undefined) basis.ecp = options.ecp.map((e) => ({ ...e }));
    if (options.nElectrons !== undefined) {
      const n = integer(options.nElectrons, "nElectrons", 0);
      const total = occupationLists.flat().reduce((a, b) => a + b, 0);
      if (options.orbitalKind !== "nto" && Math.abs(total - n) > 1e-6) {
        throw new RangeError(`nElectrons is ${n} but the occupations sum to ${total}; they must agree (spec § 5.3)`);
      }
      basis.n_electrons = n;
    }
    const aoMatrix = (input: ArrayInput, what: string): Tensor => {
      const t = resolveArray(input, what, [nAo]);
      sameShape(t.shape, [nAo, nAo], what, "[n_ao, n_ao]");
      return t;
    };
    if (options.overlap !== undefined) matrices.push(["overlap", aoMatrix(options.overlap, "overlap")]);
    if (options.densityMatrix !== undefined) {
      if (unrestricted) throw new RangeError("densityMatrix accompanies moCoefficients; pass densityMatrixAlpha and densityMatrixBeta");
      matrices.push(["density_matrix", aoMatrix(options.densityMatrix, "densityMatrix")]);
    }
    if (options.densityMatrixAlpha !== undefined || options.densityMatrixBeta !== undefined) {
      if (!unrestricted) throw new RangeError("densityMatrixAlpha and densityMatrixBeta accompany the alpha and beta coefficients");
      if (options.densityMatrixAlpha === undefined || options.densityMatrixBeta === undefined) {
        throw new RangeError("pass densityMatrixAlpha and densityMatrixBeta together");
      }
      matrices.push(["density_matrix_alpha", aoMatrix(options.densityMatrixAlpha, "densityMatrixAlpha")]);
      matrices.push(["density_matrix_beta", aoMatrix(options.densityMatrixBeta, "densityMatrixBeta")]);
    }
    this.precheck("wavefunction.gto", "basis", basis, sid);
    this.precheck("wavefunction.gto", "mo_metadata", moMeta, sid);
    this.addSection({ id: sid, kind: "wavefunction.gto", ...(options.label !== undefined ? { label: options.label } : {}) });
    const jobs: Promise<MemberSpec>[] = [
      this.addJson(sid, "basis", basis, { path: `${sid}/basis.json` }),
      this.addJson(sid, "mo_metadata", moMeta, { path: `${sid}/mo_metadata.json` }),
    ];
    for (const [role, t] of [...coefficientMembers, ...matrices]) {
      jobs.push(this.addBinary(sid, role, t.data, { dtype: "float64", shape: t.shape, path: `${sid}/${role}.dat` }));
    }
    await Promise.all(jobs);
    return sid;
  }

  /**
   * `add_volume`: any `volume.*` section (spec § 5.2). `origin` and
   * `voxelVectors` are in bohr. `data` is `[nx, ny, nz]` for the scalar
   * kinds, `[nx, ny, nz, 3]` for `volume.vector`, `[nx, ny, nz, 2]` for
   * `volume.complex` and `[n, nx, ny, nz]` for `volume.series` (with `axis`).
   * A flat typed array needs a {@link Tensor} with its shape.
   */
  async addVolume(
    kind: string,
    data: ArrayInput,
    origin: readonly number[],
    voxelVectors: readonly (readonly number[])[],
    options: VolumeOptions = {},
  ): Promise<string> {
    const four: Record<string, [number | null, string]> = {
      "volume.vector": [3, "[nx, ny, nz, 3]"],
      "volume.complex": [2, "[nx, ny, nz, 2]"],
      "volume.series": [null, "[n, nx, ny, nz]"],
    };
    const scalar = [
      "volume.density",
      "volume.orbital",
      "volume.spin",
      "volume.elf",
      "volume.generic",
      "volume.potential",
      "volume.rdg",
      "volume.difference",
      "volume.transition_density",
      "volume.laplacian",
    ];
    if (!scalar.includes(kind) && !(kind in four)) throw new QvfError("E_KIND_UNKNOWN", `${JSON.stringify(kind)} is not a volume kind`);
    if (options.sectionId !== undefined) this.checkNewId(options.sectionId);
    if (dtypeOf(data) !== undefined && !isTensor(data)) {
      throw new RangeError(`${kind}: a flat typed array needs its shape; pass {data, shape}`);
    }
    const t = resolveArray(data, `${kind} data`);
    let spatial: number[];
    const form = four[kind];
    if (form === undefined) {
      if (t.shape.length !== 3) throw new RangeError(`${kind}: data must be 3-D, got shape ${JSON.stringify(t.shape)}`);
      spatial = [...t.shape];
    } else {
      if (t.shape.length !== 4 || (form[0] !== null && t.shape[3] !== form[0])) {
        throw new RangeError(`${kind}: data must be ${form[1]}, got shape ${JSON.stringify(t.shape)}`);
      }
      spatial = kind === "volume.series" ? t.shape.slice(1) : t.shape.slice(0, 3);
    }
    if ((options.axis !== undefined) !== (kind === "volume.series")) {
      throw new RangeError(`volume.series needs axis, and only volume.series takes one (${kind})`);
    }
    let axis: Obj | undefined;
    if (options.axis !== undefined) {
      const n = t.shape[0]!;
      axis = { ...options.axis };
      const values = numbers(axis.values as number[], "volume.series axis values");
      if (values.length !== n) throw new RangeError(`volume.series axis has ${values.length} values but data has ${n} slabs`);
      if (Array.isArray(axis.labels) && axis.labels.length !== n) {
        throw new RangeError(`volume.series axis has ${axis.labels.length} labels but data has ${n} slabs`);
      }
      if (axis.units !== undefined) unitName(axis.units, "volume.series axis units");
      if (axis.trajectory_ref !== undefined) this.checkRef("volume.series axis trajectory_ref", axis.trajectory_ref, ["trajectory", "reaction.path"]);
    }
    const grid: Obj = {
      origin: vector(origin, 3, `${kind} origin`),
      voxel_vectors: matrix3(voxelVectors, `${kind} voxelVectors`),
      shape: spatial,
    };
    if (options.periodic !== undefined) grid.periodic = options.periodic;
    if (options.includesEndpoint !== undefined) grid.includes_endpoint = options.includesEndpoint;
    if (options.structureRef !== undefined) {
      this.checkRef(`${kind} structureRef`, options.structureRef, ["structure"]);
      grid.structure_ref = options.structureRef;
    }
    if (options.units !== undefined) grid.units = units(options.units, kind);
    const peers: Obj = {};
    if (options.label !== undefined) peers.label = options.label;
    if (options.component !== undefined) peers.component = options.component;
    if (options.operandA !== undefined || options.operandB !== undefined) {
      if (kind !== "volume.difference") throw new RangeError("operandA and operandB apply to volume.difference only");
      if (!options.operandA || !options.operandB) throw new RangeError("volume.difference needs both operandA and operandB");
      peers.operand_a = options.operandA;
      peers.operand_b = options.operandB;
      if (options.description !== undefined) peers.description = options.description;
    }
    if (options.excitedStatesRef !== undefined || options.stateIndex !== undefined) {
      if (kind !== "volume.transition_density") throw new RangeError("excitedStatesRef and stateIndex apply to volume.transition_density only");
      if (options.excitedStatesRef === undefined) throw new RangeError("stateIndex needs excitedStatesRef (spec § 5.2)");
      peers.excited_states_ref = options.excitedStatesRef;
      if (options.stateIndex !== undefined) peers.state_index = integer(options.stateIndex, "stateIndex", 1);
    }
    const sid = options.sectionId ?? `${kind.replace(/\./g, "_")}_${this.sections.length}`;
    this.checkNewId(sid);
    this.precheck(kind, "grid", grid, sid);
    if (axis) this.precheck(kind, "axis", axis, sid);
    this.addSection({ id: sid, kind, ...peers });
    const jobs = [
      this.addJson(sid, "grid", grid, { path: `volumes/${sid}_grid.json` }),
      this.addBinary(sid, "data", t.data, { dtype: options.dtype ?? "float32", shape: t.shape, path: `volumes/${sid}.dat` }),
    ];
    if (axis) jobs.push(this.addJson(sid, "axis", axis, { path: `volumes/${sid}_axis.json` }));
    await Promise.all(jobs);
    return sid;
  }

  /** `add_density`: `addVolume("volume.density", ...)`. */
  addDensity(data: ArrayInput, origin: readonly number[], voxelVectors: readonly (readonly number[])[], options: VolumeOptions = {}): Promise<string> {
    return this.addVolume("volume.density", data, origin, voxelVectors, options);
  }

  /** `add_orbital`: `addVolume("volume.orbital", ...)`. */
  addOrbital(data: ArrayInput, origin: readonly number[], voxelVectors: readonly (readonly number[])[], options: VolumeOptions = {}): Promise<string> {
    return this.addVolume("volume.orbital", data, origin, voxelVectors, options);
  }

  /**
   * `add_excited_states`: an `excited_states` section (spec § 5.12). States
   * ascend by `index`; a numeric `excitation_energy` is written in `ev` and a
   * numeric `total_energy` in `hartree`. `nto_ref`, `density_ref` and
   * `wavefunctionRef` must name added sections of the right kind.
   */
  async addExcitedStates(states: readonly Obj[], options: ExcitedStatesOptions = {}): Promise<string> {
    const sid = options.sectionId ?? "excited_states";
    this.checkNewId(sid);
    if (!Array.isArray(states)) throw new TypeError("excited_states: states must be an array");
    let rows: Map<string, number> | null = null;
    if (options.wavefunctionRef !== undefined) {
      const wf = this.checkRef("excited_states wavefunctionRef", options.wavefunctionRef, ["wavefunction.gto"]);
      rows = new Map();
      for (const [role, key] of [
        ["mo_coefficients", ""],
        ["mo_coefficients_alpha", "alpha"],
        ["mo_coefficients_beta", "beta"],
      ] as const) {
        const m = wf.members[role] as { shape?: number[] } | undefined;
        if (m?.shape) rows.set(key, m.shape[0]!);
      }
    }
    const energy = (v: unknown, unit: string, what: string): unknown => (typeof v === "number" ? { value: finite(v, what), units: unit } : v);
    const out = states.map((state, k) => {
      const what = `excited_states: states[${k}]`;
      if (!isObject(state)) throw new TypeError(`${what} must be an object`);
      if (state.index === undefined || state.excitation_energy === undefined) throw new RangeError(`${what} needs index and excitation_energy (spec § 5.12)`);
      const s: Obj = { ...state };
      s.excitation_energy = energy(state.excitation_energy, "ev", `${what}.excitation_energy`);
      if (state.total_energy !== undefined) s.total_energy = energy(state.total_energy, "hartree", `${what}.total_energy`);
      if (state.nto_ref !== undefined) this.checkRef(`${what}.nto_ref`, state.nto_ref, ["wavefunction.gto"]);
      if (state.density_ref !== undefined) this.checkRef(`${what}.density_ref`, state.density_ref, ["volume.transition_density", "volume.difference"]);
      if (rows && Array.isArray(state.contributions)) {
        state.contributions.forEach((c: unknown, j: number) => {
          if (!isObject(c)) return;
          const n = rows.has("") ? rows.get("") : typeof c.spin === "string" ? rows.get(c.spin) : undefined;
          for (const key of ["from_mo", "to_mo"]) {
            const mo = c[key];
            if (n !== undefined && typeof mo === "number" && mo >= n) {
              throw new QvfError("E_REF_UNRESOLVED", `${what}.contributions[${j}].${key} ${mo} is not one of the ${n} orbital rows of ${options.wavefunctionRef}`);
            }
          }
        });
      }
      return s;
    });
    for (let k = 1; k < out.length; k++) {
      if (!((out[k]!.index as number) > (out[k - 1]!.index as number))) {
        throw new RangeError(`excited_states: states[${k}].index does not ascend; indices are unique and ascending (spec § 5.12)`);
      }
    }
    const payload: Obj = { states: out };
    if (options.meta !== undefined) {
      const n = options.meta.n_states;
      if (n !== undefined && n !== out.length) throw new RangeError(`excited_states: meta.n_states ${String(n)} but ${out.length} states`);
      payload.meta = { ...options.meta };
    }
    this.precheck("excited_states", "states", payload, sid);
    const peers: Obj = {};
    if (options.label !== undefined) peers.label = options.label;
    if (options.wavefunctionRef !== undefined) peers.wavefunction_ref = options.wavefunctionRef;
    this.addSection({ id: sid, kind: "excited_states", ...peers });
    await this.addJson(sid, "states", payload, { path: `${sid}/states.json` });
    return sid;
  }

  /**
   * `add_spectrum`: a `spectra.*` section (spec § 5.4). For the 1-D kinds
   * pass `frequencies` and `intensities` (plus `xUnits`, `yUnits`, `peaks`,
   * `broadening`) or the whole `payload`; `spectra.nmr` and `spectra.epr`
   * take `payload` only.
   */
  async addSpectrum(kind: string, options: SpectrumOptions = {}): Promise<string> {
    const kinds = ["spectra.ir", "spectra.raman", "spectra.uvvis", "spectra.ecd", "spectra.vcd", "spectra.nmr", "spectra.epr", "spectra.generic"];
    if (!kinds.includes(kind)) throw new QvfError("E_KIND_UNKNOWN", `${JSON.stringify(kind)} is not a spectra kind`);
    const sid = options.sectionId ?? kind.replace(/\./g, "_");
    this.checkNewId(sid);
    const oneD = kind !== "spectra.nmr" && kind !== "spectra.epr";
    let body: Obj;
    if (options.payload !== undefined) {
      body = { ...options.payload };
    } else {
      if (options.frequencies === undefined || options.intensities === undefined) {
        throw new RangeError(`${kind}: pass frequencies and intensities, or payload`);
      }
      body = { frequencies: numbers(options.frequencies, "frequencies"), intensities: numbers(options.intensities, "intensities") };
    }
    const extras: [string, unknown][] = [
      ["x_units", options.xUnits],
      ["y_units", options.yUnits],
      ["peaks", options.peaks],
      ["broadening", options.broadening],
    ];
    for (const [key, value] of extras) {
      if (value === undefined) continue;
      if (!oneD) throw new RangeError(`${kind}: ${key} applies to the 1-D spectra kinds only`);
      if (key in body) throw new RangeError(`${kind}: ${key} given both as an option and in payload`);
      body[key] = value;
    }
    if (oneD) {
      for (const key of ["x_units", "y_units"]) if (body[key] !== undefined) unitName(body[key], `${kind} ${key}`);
      if (isObject(body.broadening) && body.broadening.units !== undefined) unitName(body.broadening.units, `${kind} broadening units`);
      const f = body.frequencies;
      if (!Array.isArray(f) || !Array.isArray(body.intensities)) throw new RangeError(`${kind}: the payload needs frequencies and intensities arrays (spec § 5.4)`);
      for (const key of ["intensities", "x"]) {
        const other = body[key];
        if (other !== undefined && (!Array.isArray(other) || other.length !== f.length)) {
          throw new RangeError(`${kind}: ${key} must have the ${f.length} entries of frequencies (spec § 5.4)`);
        }
      }
      if (body.states_ref !== undefined) this.checkRef(`${kind} states_ref`, body.states_ref, ["excited_states"]);
    }
    this.precheck(kind, "spectrum", body, sid);
    this.addSection({ id: sid, kind, ...(options.label !== undefined ? { label: options.label } : {}) });
    await this.addJson(sid, "spectrum", body, { path: `spectra/${sid}.json` });
    return sid;
  }

  /**
   * `add_reaction_scheme`: a `reaction.scheme` section (spec § 5.14). Every
   * participant names an added `structure` section, once per role of a step.
   */
  async addReactionScheme(steps: readonly SchemeStep[], options: ReactionSchemeOptions = {}): Promise<string> {
    const sid = options.sectionId ?? "reaction_scheme";
    this.checkNewId(sid);
    if (!Array.isArray(steps) || steps.length === 0) throw new RangeError("reaction.scheme needs at least one step");
    const seen = new Set<string>();
    const out = steps.map((step, k) => {
      const what = `reaction.scheme step ${k}`;
      if (!isObject(step)) throw new TypeError(`${what} must be an object`);
      if (typeof step.id !== "string" || !ID_RE.test(step.id)) throw new RangeError(`${what}: id must match [A-Za-z0-9_.-]+`);
      if (seen.has(step.id)) throw new RangeError(`${what}: id ${JSON.stringify(step.id)} is used by another step`);
      seen.add(step.id);
      const entry: Obj = { ...step };
      for (const role of ["reactants", "products", "reagents"] as const) {
        const items = step[role];
        if (items === undefined) {
          if (role !== "reagents") throw new RangeError(`${what} ${role} is required`);
          continue;
        }
        if (!Array.isArray(items) || (role !== "reagents" && items.length === 0)) {
          throw new RangeError(`${what} ${role} must be a non-empty list of participants`);
        }
        const parts = items.map((p, j) => {
          const e = typeof p === "string" ? { structure_ref: p } : isObject(p) ? { ...p } : null;
          if (e === null) throw new TypeError(`${what} ${role}[${j}] must be a section id or {structure_ref, coefficient?}`);
          this.checkRef(`${what} ${role}[${j}] structure_ref`, e.structure_ref, ["structure"]);
          return e;
        });
        const refs = parts.map((p) => p.structure_ref);
        const repeated = refs.filter((r, i) => refs.indexOf(r) !== i);
        if (repeated.length > 0) throw new RangeError(`${what} ${role} names ${JSON.stringify(repeated[0])} more than once; a coefficient says how many`);
        entry[role] = parts;
      }
      return entry;
    });
    const payload = { steps: out };
    this.precheck("reaction.scheme", "scheme", payload, sid);
    this.addSection({ id: sid, kind: "reaction.scheme", ...(options.label !== undefined ? { label: options.label } : {}) });
    await this.addJson(sid, "scheme", payload, { path: `${sid}/scheme.json` });
    return sid;
  }

  // -- output ----------------------------------------------------------------

  /** The manifest as it would be written now (nothing is written). */
  buildManifest(): Obj {
    const manifest: Obj = {
      qvf_version: 1,
      schema_uri: SCHEMA_URI,
      byte_order: "little",
      source: { ...this.source },
    };
    const metadata: Obj = { ...this.metadata };
    if (this.autoValues) {
      for (const [k, v] of Object.entries(this.autoValues)) if (!(k in metadata)) metadata[k] = v;
    }
    if (Object.keys(metadata).length > 0) manifest.metadata = metadata;
    for (const key of ROOT_KEYS) if (key in this.root) manifest[key] = this.root[key];
    manifest.sections = JSON.parse(JSON.stringify(this.sections)) as unknown;
    return manifest;
  }

  /**
   * Wait for every member, check the manifest against the schema (`E_SCHEMA`),
   * then write `manifest.json` and the ZIP central directory. A schema
   * refusal leaves the writer open, so the root blocks can be fixed and
   * `finish()` called again.
   */
  async finish(): Promise<void> {
    if (this.state === "finished") return;
    this.ensureOpen();
    await this.chain;
    if (this.failure !== null) throw new Error(`a member write failed; the archive cannot be finished (${errorText(this.failure)})`);
    const manifest = this.buildManifest();
    const bytes = jsonBytes(manifest, "manifest");
    if (this.schema) {
      const errors = checkManifest(manifest);
      if (errors.length > 0) throw new QvfError("E_SCHEMA", `manifest.json: ${schemaMessage(errors)}`);
    }
    this.state = "finishing";
    try {
      await this.zip.add("manifest.json", bytes, METHOD_STORE);
      await this.zip.finish();
    } catch (err) {
      this.failure = err;
      throw err;
    }
    this.state = "finished";
  }

  /** `finish()`, then the archive bytes (in-memory mode). */
  async toBytes(): Promise<Uint8Array> {
    await this.finish();
    return this.zip.bytes();
  }

  /** `finish()`, then the archive as a `Blob` (in-memory mode). */
  async toBlob(type?: string): Promise<Blob> {
    await this.finish();
    return this.zip.blob(type);
  }

  /** The finished archive (in-memory mode). */
  bytes(): Uint8Array {
    return this.zip.bytes();
  }

  /** The finished archive as a `Blob` (in-memory mode); it can exceed the largest single buffer. */
  blob(type?: string): Blob {
    return this.zip.blob(type);
  }

  /** Bytes written to the archive so far. */
  get length(): number {
    return this.zip.length;
  }

  // -- internals ---------------------------------------------------------------

  private ensureOpen(): void {
    if (this.state !== "open") throw new Error("the writer is finished");
    if (this.failure !== null) throw new Error(`a member write failed; the writer cannot continue (${errorText(this.failure)})`);
  }

  private checkKind(kind: unknown): void {
    if (typeof kind !== "string" || (!CANONICAL_KINDS.includes(kind) && !VENDOR_RE.test(kind))) {
      throw new QvfError("E_KIND_UNKNOWN", `unknown section kind ${JSON.stringify(kind)}: neither canonical nor x_<vendor>.*`);
    }
  }

  private checkNewId(id: unknown): void {
    if (typeof id !== "string" || !ID_RE.test(id)) throw new QvfError("E_SECTION_MALFORMED", `section id ${JSON.stringify(id)} must match [A-Za-z0-9_.-]+`);
    if (RESERVED_IDS.has(id)) throw new QvfError("E_RESERVED_ID", `section id ${JSON.stringify(id)} is reserved (a viewer_defaults key, spec § 3.4)`);
    if (this.ids.has(id)) throw new QvfError("E_SECTION_ID_DUPLICATE", `duplicate section id ${JSON.stringify(id)}`);
  }

  private reserveId(id: unknown): void {
    this.checkNewId(id);
    this.ids.add(id as string);
  }

  private sectionOrThrow(id: string): Section {
    this.ensureOpen();
    const s = this.sections.find((x) => x.id === id);
    if (!s) throw new RangeError(`no section ${JSON.stringify(id)} was added`);
    return s;
  }

  private checkRef(what: string, ref: unknown, kinds: readonly string[]): Section {
    const target = typeof ref === "string" ? this.sections.find((s) => s.id === ref) : undefined;
    if (!target || !kinds.includes(target.kind)) {
      const found = target ? `a ${target.kind} section` : "no such section";
      throw new QvfError("E_REF_UNRESOLVED", `${what} ${JSON.stringify(ref)} must name an added ${kinds.join(" or ")} section, found ${found}`);
    }
    return target;
  }

  private atomCount(meta: Obj): number | undefined {
    for (const key of ["symbols", "atomic_numbers", "atoms"]) {
      const v = meta[key];
      if (Array.isArray(v)) return v.length;
    }
    return undefined;
  }

  /** Check a builder payload before anything of the section is written. */
  private precheck(kind: string, role: string, payload: unknown, sid: string): void {
    const def = payloadDefinition(kind, role);
    if (def === undefined) return;
    const parsed: unknown = JSON.parse(new TextDecoder().decode(jsonBytes(payload, `section ${sid} member ${role}`)));
    this.checkPayload(def, parsed, `section ${JSON.stringify(sid)} member ${JSON.stringify(role)}`, "E_PAYLOAD_SCHEMA");
  }

  private checkPayload(def: string, doc: unknown, what: string, code: "E_SCHEMA" | "E_PAYLOAD_SCHEMA" = "E_SCHEMA"): void {
    if (!this.schema) return;
    const errors = checkDefinition(def, doc);
    if (errors.length > 0) throw new QvfError(code, `${what} (${def}): ${schemaMessage(errors)}`);
  }

  private mediaType(value: string | undefined): { media_type?: string } {
    if (value === undefined) return {};
    if (typeof value !== "string" || !MEDIA_TYPE_RE.test(value)) {
      throw new RangeError(`media type ${JSON.stringify(value)} is not type/subtype (spec § 4.2)`);
    }
    return { media_type: value };
  }

  /** Check the role and the path, then reserve the path. Nothing that can fail follows the reservation. */
  private prepareMember(section: Section, role: string, path: string): string {
    this.ensureOpen();
    if (typeof role !== "string" || !role) throw new RangeError("a member role must be a non-empty string");
    if (Object.hasOwn(section.members, role) || this.pendingRoles.has(`${section.id}\u0000${role}`)) {
      throw new RangeError(`section ${JSON.stringify(section.id)} already has a member ${JSON.stringify(role)}`);
    }
    if (typeof path !== "string" || !path || path.endsWith("/")) throw new QvfError("E_MEMBER_PATH", `illegal member path ${JSON.stringify(path)}`);
    const problem = pathProblem(path);
    if (problem) throw new QvfError("E_MEMBER_PATH", `member path ${JSON.stringify(path)} ${problem} (spec § 2.1)`);
    if (this.paths.has(path)) throw new QvfError("E_ZIP_DUPLICATE_ENTRY", `member path ${JSON.stringify(path)} is already used`);
    this.paths.add(path);
    return path;
  }

  /** Queue one entry: hash, compress and append in call order, then record its spec. */
  private write(
    section: Section,
    role: string,
    path: string,
    bytes: Uint8Array,
    method: number,
    spec: (sha256: string) => MemberSpec,
  ): Promise<MemberSpec> {
    const key = `${section.id}\u0000${role}`;
    this.pendingRoles.add(key);
    const digest = sha256Hex(bytes);
    digest.catch(() => undefined); // awaited in the chain; keep an early rejection from being reported as unhandled
    const job = this.chain.then(async () => {
      if (this.failure !== null) throw new Error("an earlier member write failed");
      try {
        const sha256 = await digest;
        await this.zip.add(path, bytes, method, this.level);
        const member = spec(sha256);
        section.members[role] = member;
        return member;
      } catch (err) {
        this.failure ??= err;
        throw err;
      } finally {
        this.pendingRoles.delete(key);
      }
    });
    this.chain = job.then(
      () => undefined,
      () => undefined,
    );
    return job;
  }
}

function errorText(err: unknown): string {
  return err instanceof Error ? err.message : String(err);
}
