/**
 * Validation findings and their stable codes (spec Appendix B).
 *
 * A code is never renamed or reused; tools match on the code, not on the
 * message. `E_...` codes are errors and `W_...` codes are warnings.
 */

/**
 * Every code of Appendix B and its severity. The validator reports `E_SCHEMA`
 * and `E_PAYLOAD_SCHEMA` through the JSON Schema checker generated from the
 * schema (`schemaCheck.ts`), or `W_SCHEMA_SKIPPED` when that check is
 * disabled, exactly like the C++ validator. (The reader throws `E_SCHEMA` for
 * the one schema rule it applies itself: a `byte_order` other than
 * `"little"`, spec § 2.2.)
 */
export const FINDING_CODES = {
  // Container (ZIP) level.
  E_ZIP_UNREADABLE: "error",
  E_ZIP_DUPLICATE_ENTRY: "error",
  E_MEMBER_PATH: "error",
  E_MEMBER_TOO_LARGE: "error",
  E_TOTAL_TOO_LARGE: "error",
  E_COMPRESSION_RATIO: "error",
  W_ORPHAN_MEMBER: "warning",
  // Manifest.
  E_MANIFEST_MISSING: "error",
  E_MANIFEST_JSON: "error",
  E_JSON_TOO_LARGE: "error",
  E_JSON_TOO_DEEP: "error",
  E_SCHEMA: "error",
  W_SCHEMA_SKIPPED: "warning",
  E_VERSION: "error",
  W_VERSION_WITHDRAWN: "warning",
  E_SOURCE: "error",
  E_SECTION_MALFORMED: "error",
  E_SECTION_ID_DUPLICATE: "error",
  E_RESERVED_ID: "error",
  E_KIND_UNKNOWN: "error",
  E_CRITICAL_UNDECLARED: "error",
  // Members and payloads.
  E_MEMBER_MALFORMED: "error",
  E_MEMBER_MISSING: "error",
  E_SHA256_MISMATCH: "error",
  E_JSON_INVALID: "error",
  E_DTYPE: "error",
  E_BINARY_SIZE: "error",
  E_UTF8: "error",
  E_PAYLOAD_SCHEMA: "error",
  E_STRUCTURE_PERIODICITY: "error",
  E_CHARGE_MULTIPLICITY: "error",
  E_REF_UNRESOLVED: "error",
  E_FILES_INDEX: "error",
  W_UNIT_UNKNOWN: "warning",
  // Root blocks.
  E_ROOT_INCONSISTENT: "error",
  // Wavefunction normalization (spec § 5.3, Appendix A.1).
  E_NORMALIZATION: "error",
  // A grid units.data against the unit of its volume kind (spec § 4.3).
  E_UNIT_CONTRADICTED: "error",
} as const;

export type Severity = "error" | "warning";

/** A finding code of Appendix B. */
export type FindingCode = keyof typeof FINDING_CODES;

type CodesOfSeverity<S extends Severity> = {
  [C in FindingCode]: (typeof FINDING_CODES)[C] extends S ? C : never;
}[FindingCode];

/** One validation finding. */
export interface Finding {
  severity: Severity;
  /** A stable code of spec Appendix B. */
  code: FindingCode;
  /**
   * Where it is: `archive` (the ZIP as a whole), `entry:<name>` (one ZIP
   * entry), `manifest.json`, or a JSON Pointer (RFC 6901) into `manifest.json`
   * such as `/sections/2/members/data`.
   */
  where: string;
  /** Free text for people. */
  message: string;
}

/** Collects findings; the types make a code's severity impossible to get wrong. */
export class Findings {
  readonly items: Finding[] = [];

  error(code: CodesOfSeverity<"error">, where: string, message: string): void {
    this.items.push({ severity: "error", code, where, message });
  }

  warning(code: CodesOfSeverity<"warning">, where: string, message: string): void {
    this.items.push({ severity: "warning", code, where, message });
  }
}
