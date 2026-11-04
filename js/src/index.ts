/**
 * `@qvf/core`: a reader, validator and writer for QVF (Quantum Visualization Format)
 * archives, in pure TypeScript for Node and browsers.
 */

export { QvfReader } from "./reader.js";
export { validate } from "./validate.js";
export type { ValidateOptions, ValidationLevel, ValidationReport } from "./validate.js";
export { FINDING_CODES } from "./findings.js";
export type { Finding, FindingCode, Severity } from "./findings.js";
export { QvfError } from "./errors.js";
export { DEFAULT_LIMITS, resolveLimits } from "./limits.js";
export type { Limits } from "./limits.js";
export { CANONICAL_KINDS, READABLE_VERSIONS } from "./kinds.js";
export { BOHR_PER_ANGSTROM, LEGACY_UNIT_NAMES, UNIT_NAMES, canonicalUnit, isUnitName } from "./units.js";
export type { UnitName } from "./units.js";
export { DTYPE_ITEMSIZE } from "./model.js";
export type {
  BinaryArray,
  BinaryMemberSpec,
  Bookmark,
  Camera,
  Constraints,
  CrossfadeHint,
  DType,
  DTypeArrays,
  DipoleMoment,
  Energies,
  EnergyByMethod,
  EnergyComponents,
  Extension,
  JsonMemberSpec,
  Manifest,
  MemberSpec,
  Metadata,
  MetadataAuthor,
  Provenance,
  ProvenanceSoftware,
  Section,
  Source,
  Thermochemistry,
  TypedArray,
  ValueWithUnits,
  Vec3,
  ViewerDefaults,
} from "./model.js";
export { ZipArchive, crc32, pathProblem } from "./zip.js";
export type { ZipEntry } from "./zip.js";
export { sha256Hex } from "./sha256.js";
export { Writer, SCHEMA_URI, jsonBytes, primitiveNorm } from "./writer.js";
export type {
  AnnotateOptions,
  ArrayInput,
  AtomInput,
  BinaryOptions,
  ExcitedStatesOptions,
  MemberOptions,
  Participant,
  ReactionSchemeOptions,
  RootBlock,
  SchemeStep,
  ShellInput,
  SpectrumOptions,
  StructureOptions,
  Tensor,
  TrajectoryOptions,
  VolumeOptions,
  WavefunctionGtoOptions,
  WriterOptions,
} from "./writer.js";
export { ZipWriter } from "./zipWriter.js";
export type { ZipSink, ZipWriterOptions } from "./zipWriter.js";
export { encodeTensor } from "./tensor.js";
export type { EncodedTensor, NestedNumbers, TensorInput } from "./tensor.js";
