/**
 * The QVF unit vocabulary (spec § 4.3, "Unit names").
 *
 * The schema's `$defs/UnitName` and `$defs/LegacyUnitName` are the source of
 * truth; this table repeats them in schema order, and a test compares the two.
 */

/** The unit names of spec § 4.3, in the order of `$defs/UnitName`. */
export const UNIT_NAMES = [
  "angstrom",
  "bohr",
  "hartree",
  "ev",
  "kcal_per_mol",
  "kj_per_mol",
  "cm-1",
  "nm",
  "debye",
  "mhz",
  "hz",
  "ppm",
  "km_per_mol",
  "a3",
  "bohr3",
  "gpa",
  "fs",
  "ps",
  "kelvin",
  "atm",
  "e",
  "amu",
  "radian",
  "degree",
  "dimensionless",
  "angstrom-1",
  "gpa-1",
  "angstrom_per_fs",
  "hartree_per_bohr",
  "hartree_per_bohr2",
  "hartree_per_bohr3",
  "ev_per_angstrom",
  "hartree_per_e",
  "hartree_per_e_bohr",
  "e_per_bohr3",
  "e_per_bohr4",
  "e_per_bohr5",
  "e_bohr",
  "e_bohr2",
  "debye_angstrom",
  "debye_angstrom2",
  "e2_bohr2_per_hartree",
  "e3_bohr3_per_hartree2",
  "c_per_m2",
  "states_per_ev",
] as const;

/** A unit name of spec § 4.3. */
export type UnitName = (typeof UNIT_NAMES)[number];

/**
 * Bohr per ångström: the reciprocal of CODATA 2018's bohr radius
 * (1 bohr = 0.529177210903 Å), rounded to ten decimals. It is the factor `f`
 * of spec § 4.3, "Lattice conventions", and the value of the Python
 * `BOHR_PER_ANGSTROM` and the C++ `qvf::kBohrPerAngstrom`.
 */
export const BOHR_PER_ANGSTROM = 1.8897261246;

/**
 * The deprecated spellings a `{value, units}` object still accepts
 * (`$defs/LegacyUnitName`, spec § 4.3 "Legacy spellings"), each with the unit
 * name it stands for.
 */
export const LEGACY_UNIT_NAMES: Readonly<Record<string, UnitName>> = {
  Eh: "hartree",
  eV: "ev",
  Hartree: "hartree",
};

const UNIT_NAME_SET: ReadonlySet<string> = new Set(UNIT_NAMES);

/** True when `value` is a unit name (case-sensitive; a legacy spelling is not one). */
export function isUnitName(value: unknown): value is UnitName {
  return typeof value === "string" && UNIT_NAME_SET.has(value);
}

/** The unit name for a unit name or a legacy spelling, else `null`. */
export function canonicalUnit(value: unknown): UnitName | null {
  if (isUnitName(value)) return value;
  if (typeof value === "string" && Object.hasOwn(LEGACY_UNIT_NAMES, value)) {
    return LEGACY_UNIT_NAMES[value] ?? null;
  }
  return null;
}
