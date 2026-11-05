//! The validator (spec § 6, Appendix B).
//!
//! [`validate_bytes`] treats the archive as untrusted input and never fails
//! for what is wrong with it: every problem is a finding with a stable code.
//! It is a port of `js/src/validate.ts`, which performs the semantic checks of
//! the Python and C++ reference validators in the same order and with the
//! same codes. The levels are staged as in Python and C++ (see [`Level`]). The JSON Schema of `manifest.json` (`E_SCHEMA`) and of every
//! JSON member with a payload definition (`E_PAYLOAD_SCHEMA`) is checked by
//! the checker generated from the schema ([`crate::schema_check`]);
//! `schema_check: false` disables it, and the validator then reports
//! `W_SCHEMA_SKIPPED` once.
//!
//! The finding codes and their locations follow the reference validators.
//! Messages are worded like the JavaScript ones but are not guaranteed to be
//! byte-identical, and findings of one rule may come in a different order
//! (JSON object keys iterate sorted here, in insertion order there).

use std::collections::{HashMap, HashSet};
use std::path::Path;

use serde_json::Value;

use crate::error::{Code, Severity};
use crate::generated::{CANONICAL_KINDS, LEGACY_UNITS, UNIT_NAMES};
use crate::json::{
    as_obj, decode_utf8, escape_pointer, g10, int_of, join, json_depth, json_equal, num, num_of, parse_json_bytes, quote, show,
    str_of, truthy, Obj,
};
use crate::limits::Limits;
use crate::reader::{declared_size_problem, sha256_hex, DType};
use crate::schema_check::{check_definition, check_manifest, payload_definition};
use crate::zip::{path_problem, ratio_exceeded, ZipArchive, ZipEntry};

/// Bohr per ångström (spec § 4.3, CODATA 2018 rounded to ten decimals).
pub const BOHR_PER_ANGSTROM: f64 = 1.8897261246;

/// How much of spec § 6 to check. Levels are cumulative and staged as in
/// the Python and C++ reference validators: `Container` reads only the ZIP
/// central directory (and reports a missing `manifest.json`), `Manifest`
/// adds `manifest.json`, `Payloads` the members, `Semantics` every semantic
/// check (cross-section ones included) without `E_REF_UNRESOLVED`, and
/// `References` and `All` everything.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Level {
    Container,
    Manifest,
    Payloads,
    Semantics,
    References,
    All,
}

/// Validation options.
#[derive(Debug, Clone)]
pub struct ValidateOptions {
    /// How much to check; [`Level::All`] by default.
    pub level: Level,
    /// The resource caps of spec § 6.1.
    pub limits: Limits,
    /// Run the JSON Schema check generated from the schema; `true` by default.
    pub schema_check: bool,
}

impl Default for ValidateOptions {
    fn default() -> Self {
        ValidateOptions { level: Level::All, limits: Limits::default(), schema_check: true }
    }
}

/// One validation finding.
#[derive(Debug, Clone)]
pub struct Finding {
    pub severity: Severity,
    /// A stable code of spec Appendix B.
    pub code: Code,
    /// `archive`, `entry:<name>`, `manifest.json`, or a JSON pointer into
    /// `manifest.json` such as `/sections/2/members/data`.
    pub location: String,
    /// Free text for people.
    pub message: String,
}

/// The result of a validation.
#[derive(Debug, Clone)]
pub struct Report {
    /// `true` when there is no finding of severity `error`.
    pub ok: bool,
    pub findings: Vec<Finding>,
    /// Number of entries in `manifest.sections` (0 when the manifest is unusable).
    pub n_sections: usize,
}

impl Report {
    /// Whether a finding of that code was reported.
    pub fn has(&self, code: Code) -> bool {
        self.findings.iter().any(|f| f.code == code)
    }
}

#[derive(Default)]
struct Findings(Vec<Finding>);

impl Findings {
    fn push(&mut self, code: Code, location: impl Into<String>, message: impl Into<String>) {
        self.0.push(Finding { severity: code.severity(), code, location: location.into(), message: message.into() });
    }
}

/// Validate an archive file. An unreadable file is the finding `E_ZIP_UNREADABLE`.
pub fn validate_file(path: impl AsRef<Path>, options: &ValidateOptions) -> Report {
    match std::fs::read(path.as_ref()) {
        Ok(bytes) => validate_bytes(bytes, options),
        Err(e) => {
            let mut out = Findings::default();
            out.push(Code::E_ZIP_UNREADABLE, "archive", format!("cannot read {}: {e}", path.as_ref().display()));
            finish(out, 0)
        }
    }
}

/// Validate a `.qvf` archive held in memory (spec § 6).
pub fn validate_bytes(bytes: Vec<u8>, options: &ValidateOptions) -> Report {
    let mut out = Findings::default();
    // The reference validators run every check from the semantics level on
    // and drop unresolved references below the references level.
    let level = if options.level >= Level::Semantics { Level::All } else { options.level };
    let n = match ZipArchive::open(bytes) {
        Ok(zip) => run(&zip, &options.limits, level, options.schema_check, &mut out),
        Err(e) => {
            out.push(Code::E_ZIP_UNREADABLE, "archive", format!("not a readable zip: {}", e.message));
            0
        }
    };
    if options.level < Level::References {
        out.0.retain(|f| f.code != Code::E_REF_UNRESOLVED);
    }
    finish(out, n)
}

fn finish(out: Findings, n_sections: usize) -> Report {
    let ok = !out.0.iter().any(|f| f.severity == Severity::Error);
    Report { ok, findings: out.0, n_sections }
}

// --- small helpers -----------------------------------------------------------------

/// Keys of `viewer_defaults` that are not section ids (spec § 3.1, § 3.4).
const RESERVED_SECTION_IDS: [&str; 4] = ["auto_open", "bookmarks", "camera", "crossfade"];

/// Generic references any section may carry, and the kind each must name (spec § 4.1).
const SECTION_REFS: [(&str, &str); 2] = [("structure_ref", "structure"), ("run_ref", "run.record")];

fn is_canonical(kind: &str) -> bool {
    CANONICAL_KINDS.contains(&kind)
}

fn is_unit_name(v: Option<&Value>) -> bool {
    str_of(v).is_some_and(|s| UNIT_NAMES.contains(&s))
}

fn canonical_unit(v: Option<&Value>) -> Option<&'static str> {
    let s = str_of(v)?;
    UNIT_NAMES.iter().find(|u| **u == s).copied().or_else(|| LEGACY_UNITS.iter().find(|(l, _)| *l == s).map(|(_, u)| *u))
}

fn section_pointer(i: usize) -> String {
    format!("/sections/{i}")
}

fn member_pointer(i: usize, role: &str) -> String {
    format!("{}/members/{}", section_pointer(i), escape_pointer(role))
}

fn sec(sid: &str) -> String {
    format!("section {}: ", quote(sid))
}

/// A member spec's `shape` when it is a list of non-negative integers.
fn shape_of(member: Option<&Value>) -> Option<Vec<f64>> {
    let shape = as_obj(member)?.get("shape")?.as_array()?;
    shape.iter().map(|d| d.as_f64().filter(|x| x.fract() == 0.0 && *x >= 0.0)).collect()
}

fn same(a: &[f64], b: &[f64]) -> bool {
    a == b
}

fn is_bool3(v: Option<&Value>) -> Option<[bool; 3]> {
    let a = v?.as_array()?;
    if a.len() != 3 {
        return None;
    }
    Some([a[0].as_bool()?, a[1].as_bool()?, a[2].as_bool()?])
}

fn vec3(v: Option<&Value>) -> Option<[f64; 3]> {
    let a = v?.as_array()?;
    if a.len() != 3 {
        return None;
    }
    Some([a[0].as_f64()?, a[1].as_f64()?, a[2].as_f64()?])
}

fn numbers(v: Option<&Value>) -> Option<Vec<f64>> {
    v?.as_array()?.iter().map(Value::as_f64).collect()
}

/// `x` with `p` significant digits, as JavaScript `Number(x.toPrecision(p))` prints it.
fn gprec(x: f64, p: i32) -> String {
    if x == 0.0 || !x.is_finite() {
        return num(x);
    }
    let digits = p - 1 - x.abs().log10().floor() as i32;
    if (0..=20).contains(&digits) {
        let s = format!("{x:.*}", digits as usize);
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        format!("{x:e}")
    }
}

fn show_vec(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|x| gprec(*x, 6)).collect::<Vec<_>>().join(", "))
}

type Docs = HashMap<String, Value>;

/// `|d| <= tol`, false for NaN (a NaN in binary data fails the comparison, as
/// in the reference validators).
fn within(d: f64, tol: f64) -> bool {
    d.abs() <= tol
}

// --- bands ----------------------------------------------------------------------------

fn bands_shape_errors(sid: &str, members: &Obj) -> Vec<(String, String)> {
    let Some(eig) = shape_of(members.get("eigenvalues")) else { return vec![] };
    if eig.len() != 3 {
        return vec![];
    }
    let mut errors = Vec::new();
    let s = sec(sid);
    if let Some(proj) = shape_of(members.get("projections")) {
        if proj.len() == 4 && !same(&proj[..3], &eig) {
            errors.push((
                "projections".into(),
                format!(
                    "{s}projections shape [{}] has [n_spin, n_k, n_bands] = [{}], but eigenvalues is [{}] (spec § 5.6)",
                    join(&proj),
                    join(&proj[..3]),
                    join(&eig)
                ),
            ));
        } else if proj.len() == 3 && !same(&proj[..2], &eig[1..]) {
            errors.push((
                "projections".into(),
                format!(
                    "{s}rank-3 projections shape [{}] has [n_k, n_bands] = [{}], but eigenvalues has [{}] (spec § 5.6)",
                    join(&proj),
                    join(&proj[..2]),
                    join(&eig[1..])
                ),
            ));
        }
    }
    if let Some(occ) = shape_of(members.get("occupations")) {
        if occ.len() == 3 && !same(&occ, &eig) {
            errors.push((
                "occupations".into(),
                format!(
                    "{s}occupations shape [{}] differs from eigenvalues [{}]; both are [n_spin, n_k, n_bands] (spec § 5.6)",
                    join(&occ),
                    join(&eig)
                ),
            ));
        }
    }
    errors
}

// --- trajectory / reaction.path -------------------------------------------------------

fn is_trajectory_kind(kind: &str) -> bool {
    kind == "trajectory" || kind == "reaction.path"
}

fn trajectory_frame_scalars(kind: &str) -> &'static [&'static str] {
    match kind {
        "trajectory" => &["times", "energies", "cell_volumes", "temperatures"],
        "reaction.path" => &["times", "energies"],
        _ => &[],
    }
}

const TRAJECTORY_ENERGY_TOLERANCE_EH: f64 = 1e-10;

fn float64_values(bytes: Option<&Vec<u8>>, member: Option<&Value>) -> Option<Vec<f64>> {
    let shape = shape_of(member)?;
    let bytes = bytes?;
    if str_of(as_obj(member)?.get("dtype")) != Some("float64") || shape.len() != 1 || bytes.len() as f64 != 8.0 * shape[0] {
        return None;
    }
    Some(
        bytes
            .chunks_exact(8)
            .map(|c| {
                let mut b = [0u8; 8];
                b.copy_from_slice(c);
                f64::from_le_bytes(b)
            })
            .collect(),
    )
}

fn trajectory_errors(
    sid: &str,
    kind: &str,
    members: &Obj,
    metadata: Option<&Value>,
    energies_bytes: Option<&Vec<u8>>,
) -> Vec<(String, String)> {
    let mut errors: Vec<(String, String)> = Vec::new();
    let s = sec(sid);
    let meta = as_obj(metadata);
    let dtype = str_of(as_obj(members.get("coords")).and_then(|c| c.get("dtype")));
    if let (Some(meta), Some(dtype)) = (meta, dtype) {
        let f32_ = str_of(meta.get("precision")) == Some("float32");
        if dtype == "float32" && !f32_ {
            let p = if meta.contains_key("precision") { show(meta.get("precision")) } else { "absent".into() };
            errors.push((
                "coords".into(),
                format!("{s}coords are float32 but metadata.precision is {p}; float32 coords require precision 'float32' (spec § 5.5)"),
            ));
        } else if dtype == "float64" && f32_ {
            errors.push((
                "coords".into(),
                format!("{s}metadata.precision is 'float32' but coords are float64; precision must name the dtype of coords (spec § 5.5)"),
            ));
        }
    }
    if let Some(coords) = shape_of(members.get("coords")) {
        if coords.len() == 3 && coords[2] == 3.0 {
            let (n_frames, n_atoms) = (coords[0], coords[1]);
            let vectors = [n_frames, n_atoms, 3.0];
            for role in ["velocities", "gradient"] {
                if let Some(shape) = shape_of(members.get(role)) {
                    if !same(&shape, &vectors) {
                        errors.push((
                            role.into(),
                            format!(
                                "{s}{role} shape [{}] is not [n_frames, n_atoms, 3] = [{}], the extents of coords (spec § 5.5)",
                                join(&shape),
                                join(&vectors)
                            ),
                        ));
                    }
                }
            }
            for role in trajectory_frame_scalars(kind) {
                if let Some(shape) = shape_of(members.get(*role)) {
                    if !same(&shape, &[n_frames]) {
                        errors.push((
                            (*role).into(),
                            format!(
                                "{s}{role} shape [{}] is not [n_frames] = [{}], the frames of coords (spec § 5.5)",
                                join(&shape),
                                num(n_frames)
                            ),
                        ));
                    }
                }
            }
            if let Some(lattice) = shape_of(members.get("lattice")) {
                if !same(&lattice, &[3.0, 3.0]) && !same(&lattice, &[n_frames, 3.0, 3.0]) {
                    errors.push((
                        "lattice".into(),
                        format!(
                            "{s}lattice shape [{}] is neither [3, 3] nor [n_frames, 3, 3] = [{},3,3] (spec § 5.5)",
                            join(&lattice),
                            num(n_frames)
                        ),
                    ));
                }
            }
        }
    }
    if let Some(meta) = meta {
        if let (Some(binary), Some(values)) =
            (float64_values(energies_bytes, members.get("energies")), numbers(meta.get("energies")))
        {
            if values.len() != binary.len() {
                errors.push((
                    "energies".into(),
                    format!(
                        "{s}the energies member has {} values but metadata.energies has {}; they must agree (spec § 5.5)",
                        binary.len(),
                        values.len()
                    ),
                ));
            } else {
                let bad: Vec<usize> =
                    (0..binary.len()).filter(|&k| !within(binary[k] - values[k], TRAJECTORY_ENERGY_TOLERANCE_EH)).collect();
                if let Some(&first) = bad.first() {
                    errors.push((
                        "energies".into(),
                        format!(
                            "{s}the energies member disagrees with metadata.energies at {} frame(s), first frame {first}: {} vs {} hartree; they must agree within 1e-10 hartree (spec § 5.5)",
                            bad.len(),
                            num(binary[first]),
                            num(values[first])
                        ),
                    ));
                }
            }
        }
        if kind == "trajectory" {
            let pbc = meta.get("pbc");
            let has_lattice = members.contains_key("lattice");
            if has_lattice && matches!(pbc, None | Some(Value::Null)) {
                errors.push((
                    "metadata".into(),
                    format!("{s}the section has a lattice member but its metadata has no pbc; pbc is REQUIRED when lattice is present (spec § 5.5)"),
                ));
            } else if pbc.and_then(Value::as_array).is_some_and(|a| a.iter().any(|x| x == &Value::Bool(true))) && !has_lattice {
                errors.push((
                    "metadata".into(),
                    format!(
                        "{s}metadata.pbc {} marks a periodic axis but the section has no lattice member (spec § 5.5)",
                        show(pbc)
                    ),
                ));
            }
        }
    }
    errors
}

// --- spectra and units -----------------------------------------------------------------

fn is_spectrum_1d(kind: &str) -> bool {
    payload_definition(kind, "spectrum") == Some("SpectrumPayload")
}

fn spectrum_payload_errors(sid: &str, payload: Option<&Value>) -> Vec<String> {
    let Some(p) = as_obj(payload) else { return vec![] };
    let Some(freq) = p.get("frequencies").and_then(Value::as_array) else { return vec![] };
    let n = freq.len();
    let mut errors = Vec::new();
    for key in ["intensities", "x"] {
        if let Some(other) = p.get(key).and_then(Value::as_array) {
            if other.len() != n {
                errors.push(format!(
                    "section {}: spectrum {key} has {} entries but frequencies has {n}; they must have the same length (spec § 5.4)",
                    quote(sid),
                    other.len()
                ));
            }
        }
    }
    errors
}

fn declared_units<'a>(kind: &str, role: &str, payload: &'a Value) -> Vec<(String, &'a Value)> {
    let Some(p) = payload.as_object() else { return vec![] };
    let mut out = Vec::new();
    if let Some(units) = p.get("units").and_then(Value::as_object) {
        for (k, v) in units {
            out.push((format!("units/{k}"), v));
        }
    } else if let Some(units) = p.get("units") {
        if kind == "volume.series" && role == "axis" {
            out.push(("units".into(), units));
        }
    }
    if kind == "excited_states" {
        if let Some(units) = as_obj(p.get("meta")).and_then(|m| m.get("units")).and_then(Value::as_object) {
            for (k, v) in units {
                out.push((format!("meta/units/{k}"), v));
            }
        }
    }
    if is_spectrum_1d(kind) {
        for key in ["x_units", "y_units"] {
            if let Some(v) = p.get(key) {
                out.push((key.into(), v));
            }
        }
        if let Some(v) = as_obj(p.get("broadening")).and_then(|b| b.get("units")) {
            out.push(("broadening/units".into(), v));
        }
    }
    out
}

// --- root consistency ----------------------------------------------------------------------

const ENERGY_TOLERANCE_EH: f64 = 1e-8;

fn hartree_per_unit(units: &str) -> Option<f64> {
    Some(match units.to_ascii_lowercase().as_str() {
        "hartree" | "eh" | "ha" => 1.0,
        "ev" => 1.0 / 27.211386245988,
        "kcal_per_mol" => 1.0 / 627.5094740631,
        "kj_per_mol" => 1.0 / 2625.4996394799,
        _ => return None,
    })
}

fn energy_eh(energy: Option<&Value>) -> Option<f64> {
    let e = as_obj(energy)?;
    let value = num_of(e.get("value"))?;
    let units = str_of(e.get("units"))?;
    Some(value * hartree_per_unit(units)?)
}

fn root_consistency_errors(manifest: &Obj) -> Vec<(String, String)> {
    let Some(energies) = as_obj(manifest.get("energies")) else { return vec![] };
    let mut errors = Vec::new();
    let mut compare = |where_: String, label: String, got: Option<f64>, reference: f64, ref_label: &str| {
        if let Some(got) = got {
            if (got - reference).abs() > ENERGY_TOLERANCE_EH {
                errors.push((
                    where_,
                    format!("{label} is {got:.10} hartree but {ref_label} is {reference:.10} hartree; they state one energy and must agree within 1e-8 hartree (spec § 3.6)"),
                ));
            }
        }
    };
    let prov = as_obj(manifest.get("provenance"));
    let method = prov.and_then(|p| str_of(p.get("method")));
    let scf = prov.and_then(|p| energy_eh(p.get("scf_energy")));
    if let (Some(method), Some(scf)) = (method, scf) {
        let key = method.to_ascii_lowercase();
        let ref_label = format!("provenance.scf_energy (method {})", quote(method));
        if str_of(energies.get("method")).is_some_and(|m| m.to_ascii_lowercase() == key) {
            compare("/energies/total".into(), "energies.total".into(), energy_eh(energies.get("total")), scf, &ref_label);
        }
        if let Some(by_method) = energies.get("by_method").and_then(Value::as_array) {
            for (k, entry) in by_method.iter().enumerate() {
                if let Some(e) = entry.as_object() {
                    if str_of(e.get("method")).is_some_and(|m| m.to_ascii_lowercase() == key) {
                        compare(
                            format!("/energies/by_method/{k}/total"),
                            format!("energies.by_method[{k}].total"),
                            energy_eh(e.get("total")),
                            scf,
                            &ref_label,
                        );
                    }
                }
            }
        }
    }
    let thermo = as_obj(manifest.get("thermochemistry"));
    let zpve = thermo.and_then(|t| num_of(t.get("zpve_eh")));
    if let (Some(zpve), Some(components)) = (zpve, as_obj(energies.get("components"))) {
        compare(
            "/energies/components/zpve".into(),
            "energies.components.zpve".into(),
            energy_eh(components.get("zpve")),
            zpve,
            "thermochemistry.zpve_eh",
        );
    }
    errors
}

// --- forces / hessian / bonds -------------------------------------------------------------

struct StructureIndex<'a> {
    kinds_by_id: &'a HashMap<String, String>,
    n_atoms_of: &'a HashMap<String, f64>,
    only_structure: Option<&'a str>,
}

impl StructureIndex<'_> {
    fn is_structure(&self, r: &str) -> bool {
        self.kinds_by_id.get(r).map(String::as_str) == Some("structure")
    }
}

fn derivative_errors(
    sid: &str,
    kind: &str,
    section: &Obj,
    members: &Obj,
    meta_doc: Option<&Value>,
    structures: Option<&StructureIndex>,
) -> Vec<(Code, String, String)> {
    let mut errors = Vec::new();
    let s = sec(sid);
    let empty = Obj::new();
    let meta = as_obj(meta_doc).unwrap_or(&empty);
    let mut structure: Option<String> = None;
    if let Some(st) = structures {
        let meta_ref = meta.get("structure_ref");
        if meta.contains_key("structure_ref") && !str_of(meta_ref).is_some_and(|r| st.is_structure(r)) {
            errors.push((
                Code::E_REF_UNRESOLVED,
                "meta".into(),
                format!("{s}meta.structure_ref {} does not name a structure section (spec § 5.10)", show(meta_ref)),
            ));
        }
        structure = if section.contains_key("structure_ref") {
            str_of(section.get("structure_ref")).map(String::from)
        } else if let Some(r) = str_of(meta_ref) {
            Some(r.to_string())
        } else {
            st.only_structure.map(String::from)
        };
    }
    let n_struct = match (&structure, structures) {
        (Some(id), Some(st)) => st.n_atoms_of.get(id).copied(),
        _ => None,
    };
    let sname = structure.clone().unwrap_or_default();

    if kind == "forces" {
        if let Some(grad) = shape_of(members.get("gradient")) {
            if grad.len() != 2 || grad[1] != 3.0 {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    "gradient".into(),
                    format!("{s}gradient shape [{}] is not [n_atoms, 3] (spec § 5.10)", join(&grad)),
                ));
            } else if let Some(n) = n_struct.filter(|n| *n != grad[0]) {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    "gradient".into(),
                    format!(
                        "{s}gradient has {} rows but structure {} has {} atoms (spec § 5.10)",
                        num(grad[0]),
                        quote(&sname),
                        num(n)
                    ),
                ));
            }
        }
        if let Some(stress) = shape_of(members.get("stress")) {
            if !same(&stress, &[3.0, 3.0]) {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    "stress".into(),
                    format!("{s}stress shape [{}] is not [3, 3] (spec § 5.10)", join(&stress)),
                ));
            }
        }
        return errors;
    }

    let Some(mat) = shape_of(members.get("matrix")) else { return errors };
    if mat.len() != 2 || mat[0] != mat[1] {
        errors.push((
            Code::E_MEMBER_MALFORMED,
            "matrix".into(),
            format!("{s}matrix shape [{}] is not square (spec § 5.11)", join(&mat)),
        ));
        return errors;
    }
    let extent = mat[0];
    if extent % 3.0 != 0.0 {
        errors.push((
            Code::E_MEMBER_MALFORMED,
            "matrix".into(),
            format!("{s}matrix extent {} is not divisible by 3, so the matrix is not [3N, 3N] (spec § 5.11)", num(extent)),
        ));
        return errors;
    }
    let n = extent / 3.0;
    if let Some(ns) = n_struct.filter(|ns| *ns != n) {
        errors.push((
            Code::E_MEMBER_MALFORMED,
            "matrix".into(),
            format!(
                "{s}matrix [{}] is for {} atoms but structure {} has {} atoms (spec § 5.11)",
                join(&mat),
                num(n),
                quote(&sname),
                num(ns)
            ),
        ));
    }
    if let Some(masses) = meta.get("atom_masses_amu").and_then(Value::as_array) {
        if masses.len() as f64 != n {
            errors.push((
                Code::E_MEMBER_MALFORMED,
                "meta".into(),
                format!(
                    "{s}meta.atom_masses_amu has {} entries but the matrix is for {} atoms (spec § 5.11)",
                    masses.len(),
                    num(n)
                ),
            ));
        }
    }
    errors
}

fn bonds_errors(sid: &str, section: &Obj, payload: Option<&Value>, st: &StructureIndex) -> Vec<String> {
    let Some(r) = str_of(section.get("structure_ref")) else { return vec![] };
    if !st.is_structure(r) {
        return vec![];
    }
    let (Some(&n), Some(pairs)) = (st.n_atoms_of.get(r), as_obj(payload).and_then(|p| p.get("pairs")).and_then(Value::as_array))
    else {
        return vec![];
    };
    let mut errors = Vec::new();
    for (k, pair) in pairs.iter().enumerate() {
        let Some(pair) = pair.as_object() else { continue };
        for key in ["i", "j"] {
            if let Some(atom) = int_of(pair.get(key)) {
                if !(atom >= 0.0 && atom < n) {
                    errors.push(format!(
                        "section {}: pairs[{k}].{key} {} is not one of the {} atoms of structure {}, which structure_ref names (spec § 5.1)",
                        quote(sid),
                        num(atom),
                        num(n),
                        quote(r)
                    ));
                }
            }
        }
    }
    errors
}

// --- excited states -----------------------------------------------------------------------

type MoRows = HashMap<&'static str, f64>;

fn mo_rows(members: Option<&Value>) -> MoRows {
    let mut rows = MoRows::new();
    let Some(m) = as_obj(members) else { return rows };
    for (role, key) in [("mo_coefficients", ""), ("mo_coefficients_alpha", "alpha"), ("mo_coefficients_beta", "beta")] {
        if let Some(shape) = shape_of(m.get(role)) {
            if !shape.is_empty() {
                rows.insert(key, shape[0]);
            }
        }
    }
    rows
}

fn contribution_rows(rows: &MoRows, spin: Option<&Value>) -> Option<f64> {
    if let Some(n) = rows.get("") {
        return Some(*n);
    }
    match str_of(spin) {
        Some("alpha") => return rows.get("alpha").copied(),
        Some("beta") => return rows.get("beta").copied(),
        _ => {}
    }
    ["alpha", "beta"].iter().filter_map(|c| rows.get(c).copied()).reduce(f64::min)
}

type StateSet = HashSet<i64>;

fn state_indices(payload: Option<&Value>) -> Option<StateSet> {
    let states = as_obj(payload)?.get("states")?.as_array()?;
    Some(states.iter().filter_map(|s| int_of(as_obj(Some(s))?.get("index")).map(|x| x as i64)).collect())
}

struct ExcitedRefs<'a> {
    kinds_by_id: &'a HashMap<String, String>,
    mo_rows_of: &'a HashMap<String, MoRows>,
}

fn excited_states_errors(
    sid: &str,
    section: &Obj,
    payload: Option<&Value>,
    refs: Option<&ExcitedRefs>,
) -> Vec<(Code, String, String)> {
    let mut errors = Vec::new();
    let s = sec(sid);
    let kind_of = |r: Option<&Value>| -> Option<&str> {
        let refs = refs?;
        refs.kinds_by_id.get(str_of(r)?).map(String::as_str)
    };
    let mut wf_ref: Option<&str> = None;
    if refs.is_some() && section.contains_key("wavefunction_ref") {
        let r = section.get("wavefunction_ref");
        if kind_of(r) != Some("wavefunction.gto") {
            errors.push((
                Code::E_REF_UNRESOLVED,
                "wavefunction_ref".into(),
                format!("{s}wavefunction_ref {} does not name a wavefunction.gto section (spec § 5.12)", show(r)),
            ));
        } else {
            wf_ref = str_of(r);
        }
    }
    let Some(states) = as_obj(payload).and_then(|p| p.get("states")).and_then(Value::as_array) else { return errors };
    let empty = MoRows::new();
    let rows = match (wf_ref, refs) {
        (Some(w), Some(r)) => r.mo_rows_of.get(w).unwrap_or(&empty),
        _ => &empty,
    };
    let mut seen: HashSet<i64> = HashSet::new();
    let mut prev: Option<f64> = None;
    for (k, state) in states.iter().enumerate() {
        let Some(state) = state.as_object() else { continue };
        if let Some(idx) = int_of(state.get("index")) {
            if seen.contains(&(idx as i64)) {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    String::new(),
                    format!("{s}states[{k}].index {} repeats an earlier index; state indices are unique (spec § 5.12)", num(idx)),
                ));
            } else if let Some(p) = prev.filter(|p| idx < *p) {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    String::new(),
                    format!(
                        "{s}states[{k}].index {} follows index {}; states are listed in ascending order of index (spec § 5.12)",
                        num(idx),
                        num(p)
                    ),
                ));
            }
            seen.insert(idx as i64);
            prev = Some(idx);
        }
        if refs.is_none() {
            continue;
        }
        if state.contains_key("nto_ref") && kind_of(state.get("nto_ref")) != Some("wavefunction.gto") {
            errors.push((
                Code::E_REF_UNRESOLVED,
                String::new(),
                format!(
                    "{s}states[{k}].nto_ref {} does not name a wavefunction.gto section (spec § 5.12)",
                    show(state.get("nto_ref"))
                ),
            ));
        }
        if state.contains_key("density_ref")
            && !matches!(kind_of(state.get("density_ref")), Some("volume.transition_density") | Some("volume.difference"))
        {
            errors.push((
                Code::E_REF_UNRESOLVED,
                String::new(),
                format!(
                    "{s}states[{k}].density_ref {} does not name a volume.transition_density or volume.difference section (spec § 5.12)",
                    show(state.get("density_ref"))
                ),
            ));
        }
        let (Some(w), Some(contributions)) = (wf_ref, state.get("contributions").and_then(Value::as_array)) else { continue };
        for (c, contrib) in contributions.iter().enumerate() {
            let Some(contrib) = contrib.as_object() else { continue };
            let Some(n) = contribution_rows(rows, contrib.get("spin")) else { continue };
            for key in ["from_mo", "to_mo"] {
                if let Some(mo) = int_of(contrib.get(key)) {
                    if mo < 0.0 || mo >= n {
                        errors.push((
                            Code::E_MEMBER_MALFORMED,
                            String::new(),
                            format!(
                                "{s}states[{k}].contributions[{c}].{key} {} is not one of the {} orbital rows of wavefunction {} (spec § 5.12)",
                                num(mo),
                                num(n),
                                quote(w)
                            ),
                        ));
                    }
                }
            }
        }
    }
    if let Some(meta) = as_obj(payload).and_then(|p| as_obj(p.get("meta"))) {
        if let Some(n) = int_of(meta.get("n_states")) {
            if n != states.len() as f64 {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    String::new(),
                    format!("{s}meta.n_states is {} but states has {} entries (spec § 5.12)", num(n), states.len()),
                ));
            }
        }
    }
    errors
}

// --- wavefunction.bloch --------------------------------------------------------------------

const BLOCH_WEIGHT_TOLERANCE: f64 = 1e-6;

struct PeriodicIndex<'a> {
    kinds_by_id: &'a HashMap<String, String>,
    periodic_of: &'a HashMap<String, bool>,
}

fn bloch_errors(
    sid: &str,
    section: &Obj,
    members: &Obj,
    docs: &Docs,
    structures: Option<&PeriodicIndex>,
) -> Vec<(Code, String, String)> {
    let mut errors: Vec<(Code, String, String)> = Vec::new();
    let mut err = |code: Code, where_: &str, msg: String| {
        errors.push((code, where_.to_string(), format!("section {}: {msg} (spec § 5.3.1)", quote(sid))));
    };
    let mut coeff = shape_of(members.get("mo_coefficients"));
    if let Some(c) = &coeff {
        if c.len() != 5 || c[4] != 2.0 || (c[1] != 1.0 && c[1] != 2.0) {
            err(
                Code::E_MEMBER_MALFORMED,
                "members/mo_coefficients",
                format!("mo_coefficients shape [{}] is not [n_k, n_spin, n_mo, n_ao, 2] with n_spin 1 or 2", join(c)),
            );
            coeff = None;
        }
    }
    let mut reference: Option<Vec<f64>> = coeff.as_ref().map(|c| c[..3].to_vec());
    let mut ref_role = "mo_coefficients";
    for role in ["energies", "occupations"] {
        let Some(shape) = shape_of(members.get(role)) else { continue };
        if shape.len() != 3 {
            err(
                Code::E_MEMBER_MALFORMED,
                &format!("members/{role}"),
                format!("{role} shape [{}] is not [n_k, n_spin, n_mo]", join(&shape)),
            );
        } else if let Some(r) = &reference {
            if !same(&shape, r) {
                err(
                    Code::E_MEMBER_MALFORMED,
                    &format!("members/{role}"),
                    format!("{role} shape [{}] differs from the [n_k, n_spin, n_mo] = [{}] of {ref_role}", join(&shape), join(r)),
                );
            }
        } else {
            reference = Some(shape);
            ref_role = role;
        }
    }

    if let Some(kp) = as_obj(docs.get("kpoints")) {
        let n_k = kp.get("kpoints").and_then(Value::as_array).map(Vec::len);
        if let (Some(nk), Some(r)) = (n_k, &reference) {
            if nk as f64 != r[0] {
                err(
                    Code::E_MEMBER_MALFORMED,
                    "members/kpoints",
                    format!("kpoints lists {nk} k-points but axis 0 of {ref_role} is {}", num(r[0])),
                );
            }
        }
        for key in ["weights", "labels"] {
            if let (Some(arr), Some(nk)) = (kp.get(key).and_then(Value::as_array), n_k) {
                if arr.len() != nk {
                    err(
                        Code::E_MEMBER_MALFORMED,
                        "members/kpoints",
                        format!("kpoints.{key} has {} entries for {nk} k-points", arr.len()),
                    );
                }
            }
        }
        if let Some(w) = numbers(kp.get("weights")) {
            let total: f64 = w.iter().sum();
            if (total - 1.0).abs() > BLOCH_WEIGHT_TOLERANCE {
                err(
                    Code::E_MEMBER_MALFORMED,
                    "members/kpoints",
                    format!("kpoints.weights sum to {}, not to 1 within 1e-6", num(total)),
                );
            }
        }
    }

    if let Some(meta) = as_obj(docs.get("mo_metadata")) {
        let want = match str_of(meta.get("spin")) {
            Some("restricted") => Some(1.0),
            Some("unrestricted") => Some(2.0),
            _ => None,
        };
        if let (Some(want), Some(ns)) = (want, int_of(meta.get("n_spin"))) {
            if ns != want {
                err(
                    Code::E_MEMBER_MALFORMED,
                    "members/mo_metadata",
                    format!("mo_metadata.spin {} has n_spin {}, not {}", show(meta.get("spin")), num(want), num(ns)),
                );
            }
        }
        let axes: [(&str, usize, Option<&Vec<f64>>, &str); 3] = [
            ("n_spin", 1, reference.as_ref(), ref_role),
            ("n_mo", 2, reference.as_ref(), ref_role),
            ("n_ao", 3, coeff.as_ref(), "mo_coefficients"),
        ];
        for (key, axis, dims, name) in axes {
            if let (Some(value), Some(dims)) = (int_of(meta.get(key)), dims) {
                if value != dims[axis] {
                    err(
                        Code::E_MEMBER_MALFORMED,
                        "members/mo_metadata",
                        format!("mo_metadata.{key} is {} but axis {axis} of {name} is {}", num(value), num(dims[axis])),
                    );
                }
            }
        }
        if let Some(per_k) = meta.get("n_mo_per_k").and_then(Value::as_array) {
            if let Some(r) = &reference {
                if per_k.len() as f64 != r[0] {
                    err(
                        Code::E_MEMBER_MALFORMED,
                        "members/mo_metadata",
                        format!("mo_metadata.n_mo_per_k has {} entries but axis 0 of {ref_role} is {}", per_k.len(), num(r[0])),
                    );
                }
            }
            let n_mo = match &reference {
                Some(r) => Some(r[2]),
                None => int_of(meta.get("n_mo")),
            };
            if let Some(n_mo) = n_mo.filter(|x| x.fract() == 0.0) {
                if per_k.iter().any(|v| v.as_f64().is_some_and(|x| x > n_mo)) {
                    err(
                        Code::E_MEMBER_MALFORMED,
                        "members/mo_metadata",
                        format!("mo_metadata.n_mo_per_k has an entry above n_mo = {}", num(n_mo)),
                    );
                }
            }
        }
    }

    let basis = as_obj(docs.get("basis"));
    if let (Some(b), Some(c)) = (basis, &coeff) {
        if let Some(n) = int_of(b.get("n_ao")) {
            if n != c[3] {
                err(
                    Code::E_MEMBER_MALFORMED,
                    "members/basis",
                    format!("basis.n_ao is {} but axis 3 of mo_coefficients is {}", num(n), num(c[3])),
                );
            }
        }
    }

    let Some(st) = structures else { return errors };
    let names_structure = |r: Option<&Value>| -> Option<String> {
        let r = str_of(r)?;
        (st.kinds_by_id.get(r).map(String::as_str) == Some("structure")).then(|| r.to_string())
    };
    let (structure, where_) = if section.contains_key("structure_ref") {
        match names_structure(section.get("structure_ref")) {
            Some(s) => (s, "structure_ref"),
            None => return errors, // reported by the generic rule (spec § 4.1)
        }
    } else if let Some(b) = basis {
        match names_structure(b.get("structure_ref")) {
            Some(s) => (s, "members/basis"),
            None => {
                err(
                    Code::E_REF_UNRESOLVED,
                    "members/basis",
                    format!("basis.structure_ref {} does not name a structure section", show(b.get("structure_ref"))),
                );
                return errors;
            }
        }
    } else {
        return errors;
    };
    if st.periodic_of.get(&structure) == Some(&false) {
        err(
            Code::E_MEMBER_MALFORMED,
            where_,
            format!("structure {} has no periodic axis; Bloch orbitals need a periodic structure", quote(&structure)),
        );
    }
    errors
}

// --- properties -------------------------------------------------------------------------------

fn is_property_kind(kind: &str) -> bool {
    matches!(kind, "properties.electric" | "properties.elastic" | "properties.dielectric")
}

const PROPERTY_TOLERANCE: f64 = 1e-6;

fn matrix_of(value: Option<&Value>, rows: usize, cols: usize) -> Option<Vec<f64>> {
    let v = value?.as_array()?;
    if v.len() != rows {
        return None;
    }
    let mut out = Vec::with_capacity(rows * cols);
    for row in v {
        let row = row.as_array()?;
        if row.len() != cols {
            return None;
        }
        for x in row {
            out.push(x.as_f64()?);
        }
    }
    Some(out)
}

fn first_asymmetry(m: &[f64], n: usize) -> Option<(usize, usize)> {
    let scale = m.iter().fold(0.0f64, |s, x| s.max(x.abs()));
    for i in 0..n {
        for j in i + 1..n {
            if (m[i * n + j] - m[j * n + i]).abs() > PROPERTY_TOLERANCE * scale {
                return Some((i, j));
            }
        }
    }
    None
}

fn property_errors(
    sid: &str,
    kind: &str,
    section: &Obj,
    payload: Option<&Value>,
    manifest: &Obj,
    structures: Option<&StructureIndex>,
) -> Vec<(Code, String)> {
    let Some(p) = as_obj(payload) else { return vec![] };
    let s = sec(sid);
    let mut errors = Vec::new();
    if kind == "properties.elastic" {
        for name in ["stiffness", "compliance"] {
            if let Some(m) = matrix_of(p.get(name), 6, 6) {
                if let Some((i, j)) = first_asymmetry(&m, 6) {
                    errors.push((
                        Code::E_MEMBER_MALFORMED,
                        format!(
                            "{s}{name} is not symmetric: entry ({i}, {j}) is {} but entry ({j}, {i}) is {}; they must agree within 1e-6 of the largest entry (spec § 5.13)",
                            num(m[i * 6 + j]),
                            num(m[j * 6 + i])
                        ),
                    ));
                }
            }
        }
    } else if kind == "properties.dielectric" {
        if let (Some(e), Some(ion), Some(st)) =
            (matrix_of(p.get("electronic"), 3, 3), matrix_of(p.get("ionic"), 3, 3), matrix_of(p.get("static"), 3, 3))
        {
            for k in 0..9 {
                let sum = e[k] + ion[k];
                if (st[k] - sum).abs() > PROPERTY_TOLERANCE {
                    errors.push((
                        Code::E_MEMBER_MALFORMED,
                        format!(
                            "{s}static[{}][{}] is {} but electronic + ionic is {}; static must equal their sum within 1e-6 (spec § 5.13)",
                            k / 3,
                            k % 3,
                            num(st[k]),
                            num(sum)
                        ),
                    ));
                    break;
                }
            }
        }
        if let (Some(born), Some(st)) = (p.get("born_charges").and_then(Value::as_array), structures) {
            let r: Option<&str> =
                if section.contains_key("structure_ref") { str_of(section.get("structure_ref")) } else { st.only_structure };
            let r_shown = if section.contains_key("structure_ref") {
                match section.get("structure_ref") {
                    Some(Value::String(x)) => x.clone(),
                    Some(v) => v.to_string(),
                    None => String::new(),
                }
            } else {
                st.only_structure.unwrap_or("null").to_string()
            };
            if let Some(&n) = r.and_then(|r| st.n_atoms_of.get(r)) {
                if born.len() as f64 != n {
                    errors.push((
                        Code::E_MEMBER_MALFORMED,
                        format!(
                            "{s}born_charges has {} tensors but structure {} has {} atoms; there must be one per atom (spec § 5.13)",
                            born.len(),
                            quote(&r_shown),
                            num(n)
                        ),
                    ));
                }
            }
        }
    } else if let (Some(dipole), Some(root)) = (as_obj(p.get("dipole")), as_obj(manifest.get("dipole_moment"))) {
        if let (Some(vec), Some(root_vec)) = (numbers(dipole.get("vector")), numbers(root.get("vector_debye"))) {
            if vec.len() == 3 && root_vec.len() == 3 && vec.iter().zip(&root_vec).any(|(a, b)| (a - b).abs() > PROPERTY_TOLERANCE)
            {
                errors.push((
                    Code::E_ROOT_INCONSISTENT,
                    format!(
                        "{s}dipole.vector {} differs from the root dipole_moment.vector_debye {}; the two state one dipole and must agree within 1e-6 debye (spec § 5.13)",
                        show(dipole.get("vector")),
                        show(root.get("vector_debye"))
                    ),
                ));
            }
        }
        if let (Some(t), Some(rt)) = (num_of(dipole.get("total")), num_of(root.get("total_debye"))) {
            if (t - rt).abs() > PROPERTY_TOLERANCE {
                errors.push((
                    Code::E_ROOT_INCONSISTENT,
                    format!(
                        "{s}dipole.total {} differs from the root dipole_moment.total_debye {}; the two state one dipole and must agree within 1e-6 debye (spec § 5.13)",
                        show(dipole.get("total")),
                        show(root.get("total_debye"))
                    ),
                ));
            }
        }
    }
    errors
}

fn spectrum_states_errors(
    sid: &str,
    payload: &Obj,
    kinds_by_id: &HashMap<String, String>,
    indices_of: &HashMap<String, Option<StateSet>>,
) -> Vec<String> {
    if !payload.contains_key("states_ref") {
        return vec![];
    }
    let r = payload.get("states_ref");
    let s = sec(sid);
    let Some(rs) = str_of(r).filter(|x| kinds_by_id.get(*x).map(String::as_str) == Some("excited_states")) else {
        return vec![format!("{s}spectrum states_ref {} does not name an excited_states section (spec § 5.4)", show(r))];
    };
    let (Some(Some(indices)), Some(peaks)) = (indices_of.get(rs), payload.get("peaks").and_then(Value::as_array)) else {
        return vec![];
    };
    let mut errors = Vec::new();
    for (k, peak) in peaks.iter().enumerate() {
        let state = as_obj(Some(peak)).and_then(|p| as_obj(p.get("assignment"))).and_then(|a| int_of(a.get("state_index")));
        if let Some(st) = state {
            if !indices.contains(&(st as i64)) {
                errors.push(format!(
                    "{s}peaks[{k}].assignment.state_index {} names no state of excited_states section {} (spec § 5.4)",
                    num(st),
                    quote(rs)
                ));
            }
        }
    }
    errors
}

fn transition_density_errors(
    sid: &str,
    section: &Obj,
    kinds_by_id: &HashMap<String, String>,
    indices_of: &HashMap<String, Option<StateSet>>,
) -> Vec<(String, String)> {
    if !section.contains_key("excited_states_ref") {
        return vec![];
    }
    let r = section.get("excited_states_ref");
    let s = sec(sid);
    let Some(rs) = str_of(r).filter(|x| kinds_by_id.get(*x).map(String::as_str) == Some("excited_states")) else {
        return vec![(
            "excited_states_ref".into(),
            format!("{s}excited_states_ref {} does not name an excited_states section (spec § 5.2)", show(r)),
        )];
    };
    if let (Some(state), Some(Some(indices))) = (int_of(section.get("state_index")), indices_of.get(rs)) {
        if !indices.contains(&(state as i64)) {
            return vec![(
                "state_index".into(),
                format!("{s}state_index {} names no state of excited_states section {} (spec § 5.2)", num(state), quote(rs)),
            )];
        }
    }
    vec![]
}

/// Volume kinds that may name the wavefunction they were computed from (spec § 5.2).
fn is_wavefunction_volume_kind(kind: &str) -> bool {
    matches!(kind, "volume.density" | "volume.spin" | "volume.orbital")
}

/// A `volume.density`, `volume.spin` or `volume.orbital` section's
/// `wavefunction_ref` names a `wavefunction.gto` section. For a
/// `volume.orbital`, `spin` is absent when that section carries
/// `mo_coefficients` and present, whenever `mo_index` is, when it carries the
/// alpha and beta coefficients, and `mo_index` is a row of the member they
/// select, from the manifest shapes (spec § 5.2, Volumes of a wavefunction).
/// Returns `(field, message)` pairs, all `E_REF_UNRESOLVED`.
fn wavefunction_volume_errors(
    sid: &str,
    kind: &str,
    section: &Obj,
    kinds_by_id: &HashMap<String, String>,
    mo_rows_of: &HashMap<String, MoRows>,
) -> Vec<(String, String)> {
    if !section.contains_key("wavefunction_ref") {
        return vec![];
    }
    let r = section.get("wavefunction_ref");
    let s = sec(sid);
    let Some(wf) = str_of(r).filter(|x| kinds_by_id.get(*x).map(String::as_str) == Some("wavefunction.gto")) else {
        return vec![(
            "wavefunction_ref".into(),
            format!("{s}wavefunction_ref {} does not name a wavefunction.gto section (spec § 5.2)", show(r)),
        )];
    };
    if kind != "volume.orbital" {
        return vec![];
    }
    let mut errors = Vec::new();
    let empty = MoRows::new();
    let rows = mo_rows_of.get(wf).unwrap_or(&empty);
    let mo = section.get("mo_index");
    let spin = section.get("spin");
    let mut n: Option<f64> = None;
    let mut member = "mo_coefficients".to_string();
    if let Some(single) = rows.get("") {
        if section.contains_key("spin") {
            errors.push((
                "spin".into(),
                format!(
                    "{s}spin {} selects a spin channel, but wavefunction {} carries one set of orbitals in mo_coefficients (spec § 5.2)",
                    show(spin),
                    quote(wf)
                ),
            ));
        }
        n = Some(*single);
    } else if rows.contains_key("alpha") || rows.contains_key("beta") {
        if !section.contains_key("spin") {
            if section.contains_key("mo_index") {
                errors.push((
                    "mo_index".into(),
                    format!(
                        "{s}mo_index {} needs spin: wavefunction {} carries alpha and beta coefficients (spec § 5.2)",
                        show(mo),
                        quote(wf)
                    ),
                ));
            }
            return errors;
        }
        if let Some(channel @ ("alpha" | "beta")) = str_of(spin) {
            n = rows.get(channel).copied();
            member = format!("mo_coefficients_{channel}");
        }
    }
    if let (Some(m), Some(n)) = (int_of(mo), n) {
        if !(m >= 0.0 && m < n) {
            errors.push((
                "mo_index".into(),
                format!(
                    "{s}mo_index {} is not one of the {} rows of {member} of wavefunction {} (spec § 5.2)",
                    num(m),
                    num(n),
                    quote(wf)
                ),
            ));
        }
    }
    errors
}

fn atom_properties_meta_errors(sid: &str, members: &Obj, docs: &Docs) -> Vec<(String, String)> {
    let custom: Vec<&String> = members.keys().filter(|r| r.starts_with("custom.")).collect();
    if custom.is_empty() {
        return vec![];
    }
    let s = format!("section {}: custom role ", quote(sid));
    if !members.contains_key("meta") {
        return custom
            .iter()
            .map(|r| ((*r).clone(), format!("{s}{} has no meta member to describe it (spec § 5.7)", quote(r))))
            .collect();
    }
    let Some(meta) = docs.get("meta") else { return vec![] };
    let Some(meta) = meta.as_object() else {
        return custom
            .iter()
            .map(|r| ((*r).clone(), format!("{s}{}: the meta member is not a JSON object (spec § 5.7)", quote(r))))
            .collect();
    };
    custom
        .iter()
        .filter(|r| !meta.get(r.as_str()).is_some_and(Value::is_object))
        .map(|r| ((*r).clone(), format!("{s}{} has no entry in the meta member (spec § 5.7)", quote(r))))
        .collect()
}

// --- wavefunction.gto (spec § 5.3, Appendix A.1) ---------------------------------------------

const N_ELECTRONS_TOLERANCE: f64 = 1e-6;
const NORMALIZATION_RTOL: f64 = 1e-6;

struct ShellInfo {
    l: u32,
    pure: bool,
    prims: Option<(Vec<f64>, Vec<f64>)>,
}

fn basis_shells(basis: Option<&Value>) -> Option<Vec<ShellInfo>> {
    let b = as_obj(basis)?;
    let shells = b.get("shells")?.as_array()?;
    let default_pure = match b.get("pure") {
        Some(v) => v.as_bool()?,
        None => true,
    };
    let mut out = Vec::new();
    for sh in shells {
        let sh = sh.as_object()?;
        let l = int_of(sh.get("l")).filter(|l| (0.0..=7.0).contains(l))? as u32;
        let pure = match sh.get("pure") {
            Some(v) => v.as_bool()?,
            None => default_pure,
        };
        let e = numbers(sh.get("exponents"));
        let c = numbers(sh.get("coefficients"));
        let prims = match (e, c) {
            (Some(e), Some(c)) if e.len() == c.len() && !e.is_empty() && e.iter().all(|a| *a > 0.0) => Some((e, c)),
            _ => None,
        };
        out.push(ShellInfo { l, pure, prims });
    }
    Some(out)
}

fn shell_size(l: u32, pure: bool) -> usize {
    let l = l as usize;
    if pure {
        2 * l + 1
    } else {
        (l + 1) * (l + 2) / 2
    }
}

fn basis_n_ao(basis: Option<&Value>) -> Option<f64> {
    if let Some(shells) = basis_shells(basis) {
        return Some(shells.iter().map(|s| shell_size(s.l, s.pure)).sum::<usize>() as f64);
    }
    int_of(as_obj(basis)?.get("n_ao")).filter(|n| *n >= 0.0)
}

fn axial_self_overlap(l: u32, e: &[f64], c: &[f64]) -> f64 {
    let mut total = 0.0;
    for (i, a) in e.iter().enumerate() {
        for (j, b) in e.iter().enumerate() {
            total += c[i] * c[j] * ((2.0 * (a * b).sqrt()) / (a + b)).powf(f64::from(l) + 1.5);
        }
    }
    total
}

fn odd_double_factorial(n: u32) -> f64 {
    let mut out = 1.0;
    let mut k = 1;
    while k < 2 * n {
        out *= f64::from(k);
        k += 2;
    }
    out
}

fn shell_components(l: u32, pure: bool) -> Vec<(String, f64)> {
    let mut out = Vec::new();
    if pure {
        let l = l as i32;
        for m in -l..=l {
            out.push((if m > 0 { format!("m = +{m}") } else { format!("m = {m}") }, 1.0));
        }
        return out;
    }
    for lx in (0..=l).rev() {
        for ly in (0..=l - lx).rev() {
            let lz = l - lx - ly;
            let name = "x".repeat(lx as usize) + &"y".repeat(ly as usize) + &"z".repeat(lz as usize);
            let name = if name.is_empty() { "s".to_string() } else { name };
            out.push((
                format!("Cartesian {name} [{lx}, {ly}, {lz}]"),
                odd_double_factorial(lx) * odd_double_factorial(ly) * odd_double_factorial(lz) / odd_double_factorial(l),
            ));
        }
    }
    out
}

fn normalization_errors(sid: &str, basis: Option<&Value>, overlap: Option<&Value>, bytes: Option<&Vec<u8>>) -> Vec<String> {
    let (Some(shells), Some(shape), Some(bytes), Some(ov)) = (basis_shells(basis), shape_of(overlap), bytes, as_obj(overlap))
    else {
        return vec![];
    };
    let n_ao: usize = shells.iter().map(|s| shell_size(s.l, s.pure)).sum();
    if !same(&shape, &[n_ao as f64, n_ao as f64]) || str_of(ov.get("dtype")) != Some("float64") || bytes.len() != 8 * n_ao * n_ao
    {
        return vec![];
    }
    let get = |k: usize| {
        let mut b = [0u8; 8];
        b.copy_from_slice(&bytes[8 * k..8 * k + 8]);
        f64::from_le_bytes(b)
    };
    let mut n_bad = 0;
    let mut first = String::new();
    let mut ao = 0;
    for (k, s) in shells.iter().enumerate() {
        let components = shell_components(s.l, s.pure);
        if let Some((e, c)) = &s.prims {
            let s_axial = axial_self_overlap(s.l, e, c);
            for (j, (name, factor)) in components.iter().enumerate() {
                let i = ao + j;
                let want = s_axial * factor;
                let got = get(i * n_ao + i);
                if !within(got - want, NORMALIZATION_RTOL * want.abs()) {
                    if n_bad == 0 {
                        first =
                            format!("AO {i} (shell {k}, l = {}, {name}) is {} but the basis gives {}", s.l, g10(got), g10(want));
                    }
                    n_bad += 1;
                }
            }
        }
        ao += components.len();
    }
    if n_bad == 0 {
        return vec![];
    }
    vec![format!(
        "section {}: {n_bad} of {n_ao} overlap diagonal entries disagree with the self-overlaps Appendix A.1 gives for the basis; first: {first}. The basis.json coefficients are not on N_i-normalized primitives, or the AO order disagrees with Appendix A.2 (spec § 5.3)",
        quote(sid)
    )]
}

fn coefficients_complex(meta: Option<&Value>) -> Option<bool> {
    let m = as_obj(meta)?;
    if let Some(c) = m.get("complex") {
        return c.as_bool();
    }
    if let Some(enc) = m.get("coefficient_encoding") {
        return match enc.as_str() {
            Some("complex_split_last_axis") => Some(true),
            Some("real") => Some(false),
            _ => None,
        };
    }
    Some(false)
}

fn is_unrestricted_nto(meta: &Obj) -> bool {
    meta.contains_key("nto_role")
        && str_of(meta.get("orbital_kind")) == Some("nto")
        && str_of(meta.get("spin")) == Some("unrestricted")
}

fn unrestricted_nto_errors(meta: Option<&Value>, members: &Obj) -> Vec<String> {
    let Some(meta) = as_obj(meta).filter(|m| is_unrestricted_nto(m)) else { return vec![] };
    let mut errors = Vec::new();
    let rows = mo_rows(Some(&Value::Object(members.clone())));
    for ch in ["alpha", "beta"] {
        let Some(block) = as_obj(meta.get(ch)) else { continue };
        match rows.get(ch) {
            None => errors.push(format!("unrestricted NTOs need mo_coefficients_{ch} (spec § 5.3)")),
            Some(&n) => {
                let lens_bad = ["energies", "occupations"]
                    .iter()
                    .any(|k| block.get(*k).and_then(Value::as_array).is_some_and(|a| a.len() as f64 != n));
                if lens_bad || int_of(block.get("n_mo")).is_some_and(|m| m != n) {
                    errors.push(format!("NTO {ch} metadata must have one entry per coefficient row (spec § 5.3)"));
                }
            }
        }
        if let Some(occ) = numbers(block.get("occupations")) {
            if occ.iter().any(|x| *x < 0.0) || occ.windows(2).any(|w| w[0] < w[1]) {
                errors.push(format!("NTO {ch} weights must be nonnegative and non-increasing (spec § 5.3)"));
            }
        }
    }
    errors
}

struct NtoSection {
    i: usize,
}

fn nto_pair_errors(
    meta: Option<&Value>,
    members: &Obj,
    ntos: &[(Option<&Value>, &Obj)],
    kinds: &HashMap<String, String>,
    indices: &HashMap<String, Option<StateSet>>,
) -> Vec<(Code, String)> {
    let Some(meta) = as_obj(meta).filter(|m| is_unrestricted_nto(m)) else { return vec![] };
    let (Some(r), Some(index), Some(role)) =
        (str_of(meta.get("state_ref")), int_of(meta.get("state_index")), str_of(meta.get("nto_role")))
    else {
        return vec![];
    };
    if role != "hole" && role != "particle" {
        return vec![];
    }
    let mut errors = Vec::new();
    let states = indices.get(r);
    let state_missing = matches!(states, Some(Some(set)) if !set.contains(&(index as i64)));
    if kinds.get(r).map(String::as_str) != Some("excited_states") || state_missing {
        errors.push((Code::E_REF_UNRESOLVED, "NTO state_ref/state_index must name an excited state (spec § 5.3)".into()));
    }
    let opposite = if role == "hole" { "particle" } else { "hole" };
    let partners: Vec<&(Option<&Value>, &Obj)> = ntos
        .iter()
        .filter(|(m, _)| {
            as_obj(*m).is_some_and(|m| {
                str_of(m.get("orbital_kind")) == Some("nto")
                    && str_of(m.get("state_ref")) == Some(r)
                    && num_of(m.get("state_index")) == Some(index)
                    && str_of(m.get("nto_role")) == Some(opposite)
            })
        })
        .collect();
    if partners.len() != 1 {
        errors.push((Code::E_REF_UNRESOLVED, "NTOs must have one opposite-role section for the same state (spec § 5.3)".into()));
        return errors;
    }
    let (other_meta, other_members) = partners[0];
    let other = as_obj(*other_meta).expect("filtered above");
    if str_of(other.get("spin")) != Some("unrestricted") {
        errors.push((Code::E_MEMBER_MALFORMED, "NTO partners must have the same spin layout (spec § 5.3)".into()));
        return errors;
    }
    let rows = mo_rows(Some(&Value::Object(members.clone())));
    let other_rows = mo_rows(Some(&Value::Object((*other_members).clone())));
    for ch in ["alpha", "beta"] {
        if let (Some(a), Some(b)) = (rows.get(ch), other_rows.get(ch)) {
            if a != b {
                errors.push((Code::E_MEMBER_MALFORMED, format!("NTO partners must have the same {ch} row count (spec § 5.3)")));
            }
        }
        if let (Some(a), Some(b)) = (as_obj(meta.get(ch)), as_obj(other.get(ch))) {
            if let (Some(av), Some(bv)) = (numbers(a.get("occupations")), numbers(b.get("occupations"))) {
                if av != bv {
                    errors
                        .push((Code::E_MEMBER_MALFORMED, format!("NTO partners must have identical {ch} weights (spec § 5.3)")));
                }
            }
        }
    }
    errors
}

fn wavefunction_errors(sid: &str, members: &Obj, docs: &Docs, provenance: Option<&Value>) -> Vec<(Code, String, String)> {
    let mut errors = Vec::new();
    let s = sec(sid);
    let basis = docs.get("basis");
    let meta = docs.get("mo_metadata");

    if let Some(cplx) = coefficients_complex(meta) {
        for role in ["mo_coefficients", "mo_coefficients_alpha", "mo_coefficients_beta"] {
            let Some(shape) = shape_of(members.get(role)) else { continue };
            if cplx && (shape.len() != 3 || shape[2] != 2.0) {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    role.to_string(),
                    format!("{s}{role} shape [{}] has no trailing [re, im] axis of 2, but mo_metadata declares complex coefficients (spec § 5.3)", join(&shape)),
                ));
            } else if !cplx && shape.len() != 2 {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    role.to_string(),
                    format!("{s}{role} shape [{}] is not [n_mo, n_ao], and mo_metadata does not declare complex coefficients (spec § 5.3)", join(&shape)),
                ));
            }
        }
    }

    if let Some(m) = as_obj(meta) {
        let mut sets: Vec<(String, &Obj)> = vec![(String::new(), m)];
        for ch in ["alpha", "beta"] {
            if let Some(block) = as_obj(m.get(ch)) {
                sets.push((format!("{ch}."), block));
            }
        }
        for (prefix, block) in sets {
            if let (Some(labels), Some(energies)) =
                (block.get("symmetry_labels").and_then(Value::as_array), block.get("energies").and_then(Value::as_array))
            {
                if labels.len() != energies.len() {
                    errors.push((
                        Code::E_MEMBER_MALFORMED,
                        "mo_metadata".into(),
                        format!(
                            "{s}mo_metadata {prefix}symmetry_labels has {} entries but {prefix}energies has {}; there is one label per MO (spec § 5.3)",
                            labels.len(),
                            energies.len()
                        ),
                    ));
                }
            }
        }
    }

    let n_ao = basis_n_ao(basis);
    for role in ["overlap", "density_matrix", "density_matrix_alpha", "density_matrix_beta"] {
        if let (Some(shape), Some(n)) = (shape_of(members.get(role)), n_ao) {
            if !same(&shape, &[n, n]) {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    role.to_string(),
                    format!(
                        "{s}{role} shape [{}] is not [n_ao, n_ao] = [{}, {}] for the basis (spec § 5.3)",
                        join(&shape),
                        num(n),
                        num(n)
                    ),
                ));
            }
        }
    }
    let unrestricted = members.contains_key("mo_coefficients_alpha") && members.contains_key("mo_coefficients_beta");
    if members.contains_key("density_matrix") && !members.contains_key("mo_coefficients") {
        errors.push((
            Code::E_MEMBER_MALFORMED,
            "density_matrix".into(),
            format!("{s}density_matrix accompanies mo_coefficients only; an unrestricted wavefunction carries density_matrix_alpha and density_matrix_beta (spec § 5.3)"),
        ));
    }
    for (role, other) in [("density_matrix_alpha", "density_matrix_beta"), ("density_matrix_beta", "density_matrix_alpha")] {
        if !members.contains_key(role) {
            continue;
        }
        if !unrestricted {
            errors.push((
                Code::E_MEMBER_MALFORMED,
                role.into(),
                format!("{s}{role} accompanies mo_coefficients_alpha and mo_coefficients_beta only; a single-set wavefunction carries density_matrix (spec § 5.3)"),
            ));
        } else if !members.contains_key(other) {
            errors.push((
                Code::E_MEMBER_MALFORMED,
                role.into(),
                format!("{s}{role} without {other}; the two come together (spec § 5.3)"),
            ));
        }
    }

    let Some(n_el) = as_obj(basis).and_then(|b| num_of(b.get("n_electrons"))) else { return errors };
    if let Some(m) = as_obj(meta) {
        if str_of(m.get("occupation_semantics")) != Some("transition_weight") {
            let mut lists: Vec<Option<&Value>> = Vec::new();
            if m.contains_key("occupations") {
                lists.push(m.get("occupations"));
            }
            for ch in ["alpha", "beta"] {
                if let Some(block) = as_obj(m.get(ch)) {
                    if block.contains_key("occupations") {
                        lists.push(block.get("occupations"));
                    }
                }
            }
            let parsed: Option<Vec<Vec<f64>>> = lists.iter().map(|l| numbers(*l)).collect();
            if let Some(parsed) = parsed.filter(|p| !p.is_empty()) {
                let total: f64 = parsed.iter().flatten().sum();
                if !within(total - n_el, N_ELECTRONS_TOLERANCE) {
                    errors.push((
                        Code::E_MEMBER_MALFORMED,
                        "basis".into(),
                        format!(
                            "{s}basis.n_electrons is {} but the mo_metadata occupations sum to {}; they must agree within 1e-6 (spec § 5.3)",
                            num(n_el),
                            g10(total)
                        ),
                    ));
                }
            }
        }
    }
    if let Some(prov_n) = as_obj(provenance).and_then(|p| num_of(p.get("n_electrons"))) {
        if prov_n != n_el {
            errors.push((
                Code::E_ROOT_INCONSISTENT,
                "basis".into(),
                format!(
                    "{s}basis.n_electrons is {} but provenance.n_electrons is {}; they state one electron count and must be equal (spec § 3.2, § 5.3)",
                    num(n_el),
                    num(prov_n)
                ),
            ));
        }
    }
    errors
}

// --- volumes and grids ----------------------------------------------------------------------

fn volume_data_shape(kind: &str) -> Option<(usize, Option<f64>, &'static str)> {
    Some(match kind {
        "volume.laplacian" => (3, None, "[nx, ny, nz]"),
        "volume.vector" => (4, Some(3.0), "[nx, ny, nz, 3]"),
        "volume.complex" => (4, Some(2.0), "[nx, ny, nz, 2]"),
        "volume.series" => (4, None, "[n, nx, ny, nz]"),
        _ => return None,
    })
}

fn volume_errors(
    sid: &str,
    kind: &str,
    members: &Obj,
    docs: &Docs,
    kinds_by_id: Option<&HashMap<String, String>>,
) -> Vec<(Code, String, String)> {
    let Some((rank, trailing, form)) = volume_data_shape(kind) else { return vec![] };
    let s = sec(sid);
    let mut errors = Vec::new();
    let mut shape = shape_of(members.get("data"));
    if let Some(sh) = &shape {
        if sh.len() != rank || trailing.is_some_and(|t| sh.last() != Some(&t)) {
            errors.push((
                Code::E_MEMBER_MALFORMED,
                "data".into(),
                format!("{s}{kind} data shape [{}] is not {form} (spec § 5.2)", join(sh)),
            ));
            shape = None;
        }
    }
    let Some(axis) = as_obj(docs.get("axis")).filter(|_| kind == "volume.series") else { return errors };
    for key in ["values", "labels"] {
        if let (Some(sh), Some(seq)) = (&shape, axis.get(key).and_then(Value::as_array)) {
            if seq.len() as f64 != sh[0] {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    "axis".into(),
                    format!(
                        "{s}axis.{key} has {} entries but data [{}] has {} slabs; it needs one per slab (spec § 5.2)",
                        seq.len(),
                        join(sh),
                        num(sh[0])
                    ),
                ));
            }
        }
    }
    if let Some(kinds) = kinds_by_id {
        if axis.contains_key("trajectory_ref") {
            let r = axis.get("trajectory_ref");
            let target = str_of(r).and_then(|x| kinds.get(x)).map(String::as_str);
            if target != Some("trajectory") && target != Some("reaction.path") {
                errors.push((
                    Code::E_REF_UNRESOLVED,
                    "axis".into(),
                    format!(
                        "{s}axis.trajectory_ref {} does not name a trajectory or reaction.path section (spec § 5.2)",
                        show(r)
                    ),
                ));
            }
        }
    }
    errors
}

/// The unit name spec § 4.3 fixes for the grid `units.data` of a volume kind.
/// One grid payload serves every `volume.*` kind, so the schema cannot fix it
/// per kind.
fn volume_data_unit(kind: &str) -> Option<&'static str> {
    Some(match kind {
        "volume.density" | "volume.transition_density" => "e_per_bohr3",
        "volume.potential" => "hartree_per_e",
        "volume.rdg" => "dimensionless",
        "volume.laplacian" => "e_per_bohr5",
        _ => return None,
    })
}

/// A grid `units.data` that is a unit name other than the one § 4.3 fixes for
/// the volume kind (`E_UNIT_CONTRADICTED`), as a message. A string that is not
/// a unit name is `W_UNIT_UNKNOWN` alone.
fn grid_data_unit_error(sid: &str, kind: &str, grid: &Obj) -> Option<String> {
    let fixed = volume_data_unit(kind)?;
    let unit = as_obj(grid.get("units"))?.get("data");
    let name = str_of(unit).filter(|u| *u != fixed && is_unit_name(unit))?;
    Some(format!(
        "{}grid units.data {} contradicts the unit {} that spec § 4.3 fixes for {kind} data",
        sec(sid),
        quote(name),
        quote(fixed)
    ))
}

fn is_grid_payload(kind: &str, role: &str) -> bool {
    if role == "volume_grid" {
        return kind == "reaction.path";
    }
    role == "grid" && is_canonical(kind) && (kind.starts_with("volume.") || kind == "basis.ao")
}

const GRID_SPAN_TOLERANCE: f64 = 1e-6;

fn grid_errors(
    sid: &str,
    section: &Obj,
    grid: &Obj,
    kinds_by_id: &HashMap<String, String>,
    cells_of: &HashMap<String, &Obj>,
    only_structure: Option<&str>,
) -> Vec<(Code, String)> {
    let s = sec(sid);
    let mut errors = Vec::new();
    let r = grid.get("structure_ref");
    if grid.contains_key("structure_ref")
        && !str_of(r).is_some_and(|x| kinds_by_id.get(x).map(String::as_str) == Some("structure"))
    {
        errors.push((
            Code::E_REF_UNRESOLVED,
            format!("{s}grid structure_ref {} does not name a structure section (spec § 5.2)", show(r)),
        ));
    }
    let Some(periodic) = is_bool3(grid.get("periodic")) else { return errors };
    if !periodic.iter().any(|p| *p) {
        return errors;
    }
    let structure: Option<&str> = if section.contains_key("structure_ref") {
        str_of(section.get("structure_ref"))
    } else if let Some(x) = str_of(r) {
        Some(x)
    } else {
        only_structure
    };
    let Some(structure) = structure else { return errors };
    let Some(cell) = cells_of.get(structure) else { return errors };
    let pbc = match cell.get("pbc") {
        None | Some(Value::Null) => Some([false; 3]),
        v => is_bool3(v),
    };
    let endpoint = match grid.get("includes_endpoint") {
        None => Some(true),
        Some(v) => v.as_bool(),
    };
    let shape: Option<Vec<f64>> = grid
        .get("shape")
        .and_then(Value::as_array)
        .filter(|a| a.len() == 3)
        .and_then(|a| a.iter().map(|n| n.as_f64().filter(|x| x.fract() == 0.0)).collect());
    let vox = grid.get("voxel_vectors").and_then(Value::as_array).filter(|a| a.len() == 3);
    let (Some(pbc), Some(endpoint), Some(shape), Some(vox)) = (pbc, endpoint, shape, vox) else { return errors };
    let lattice = cell.get("lattice_vectors").and_then(Value::as_array).filter(|l| l.len() == 3);
    for a in 0..3 {
        if !periodic[a] {
            continue;
        }
        if !pbc[a] {
            errors.push((
                Code::E_MEMBER_MALFORMED,
                format!(
                    "{s}grid axis {a} is periodic but structure {} has pbc[{a}] false; a grid axis can be periodic only along a periodic axis of its structure (spec § 5.2)",
                    quote(structure)
                ),
            ));
            continue;
        }
        let (Some(step), Some(row)) = (vec3(vox.get(a)), lattice.and_then(|l| vec3(l.get(a)))) else { continue };
        let n = shape[a] - if endpoint { 1.0 } else { 0.0 };
        let span: Vec<f64> = step.iter().map(|x| x * n).collect();
        let target: Vec<f64> = row.iter().map(|x| x * BOHR_PER_ANGSTROM).collect();
        let diff = span.iter().zip(&target).map(|(p, q)| (p - q) * (p - q)).sum::<f64>().sqrt();
        let norm = target.iter().map(|q| q * q).sum::<f64>().sqrt();
        if diff > GRID_SPAN_TOLERANCE * norm {
            errors.push((
                Code::E_MEMBER_MALFORMED,
                format!(
                    "{s}periodic grid axis {a} spans {} bohr (voxel_vectors[{a}] × {}, includes_endpoint {endpoint}) but lattice vector {a} of structure {} is {} bohr; a periodic grid spans its cell within 1e-6 relative (spec § 5.2)",
                    show_vec(&span),
                    num(n),
                    quote(structure),
                    show_vec(&target)
                ),
            ));
        }
    }
    errors
}

// --- structure -----------------------------------------------------------------------------------

fn structure_payload_errors(sid: &str, payload: &Value) -> Vec<String> {
    let s = sec(sid);
    let Some(p) = payload.as_object() else { return vec![format!("{s}structure payload must be a JSON object")] };
    let mut errors = Vec::new();
    let nonnull = |k: &str| p.get(k).filter(|v| !v.is_null());
    let pbc = nonnull("pbc");
    let lattice = nonnull("lattice_vectors");
    let mut dim = nonnull("dimensionality");
    if let Some(d) = dim {
        if !d.is_number() {
            errors.push(format!("{s}dimensionality must be an integer, got {} (spec § 5.1)", show(Some(d))));
            dim = None;
        }
    }
    let dim = dim.and_then(Value::as_f64);
    let Some(pbc) = pbc else {
        if lattice.is_some() {
            errors.push(format!(
                "{s}structure payload has lattice_vectors but no 'pbc'; pbc is REQUIRED for periodic structures (spec § 5.1)"
            ));
        }
        if let Some(d) = dim.filter(|d| *d != 0.0) {
            errors.push(format!(
                "{s}dimensionality {} without 'pbc'; absent pbc means a molecule, so dimensionality must be 0 (spec § 5.1)",
                num(d)
            ));
        }
        return errors;
    };
    let Some(flags) = is_bool3(Some(pbc)) else {
        errors.push(format!("{s}pbc must be three booleans, got {}", show(Some(pbc))));
        return errors;
    };
    let n_periodic = flags.iter().filter(|x| **x).count();
    if n_periodic > 0 && !lattice.is_some_and(Value::is_array) {
        errors.push(format!(
            "{s}pbc marks {n_periodic} periodic axis/axes but lattice_vectors is {}",
            if lattice.is_none() { "null or absent" } else { "not an array" }
        ));
    }
    if let Some(d) = dim.filter(|d| *d != n_periodic as f64) {
        errors.push(format!("{s}dimensionality {} != sum(pbc) {n_periodic} (spec § 5.1 invariant)", num(d)));
    }
    errors
}

fn sum_atomic_numbers(payload: Option<&Value>) -> Option<f64> {
    let atoms = as_obj(payload)?.get("atoms")?.as_array()?;
    let mut sum = 0.0;
    for atom in atoms {
        let z = int_of(as_obj(Some(atom)).and_then(|a| a.get("atomic_number")))?;
        if z < 0.0 {
            return None;
        }
        sum += z;
    }
    Some(sum)
}

fn identity_errors(
    sid: &str,
    section: &Obj,
    sum_z: Option<f64>,
    provenance: Option<&Value>,
    only_structure: bool,
) -> Vec<(Code, &'static str, String)> {
    let mut errors = Vec::new();
    let charge = int_of(section.get("charge"));
    let mult = int_of(section.get("multiplicity")).filter(|m| *m >= 1.0);
    if let (Some(charge), Some(sum_z)) = (charge, sum_z) {
        let electrons = sum_z - charge;
        if electrons < 0.0 {
            errors.push((
                Code::E_CHARGE_MULTIPLICITY,
                "charge",
                format!(
                    "structure {}: charge {} leaves {} electrons for atoms whose atomic numbers sum to {} (spec § 5.1)",
                    quote(sid),
                    num(charge),
                    num(electrons),
                    num(sum_z)
                ),
            ));
        } else if let Some(m) = mult.filter(|m| m - 1.0 > electrons || (electrons - (m - 1.0)) % 2.0 != 0.0) {
            errors.push((
                Code::E_CHARGE_MULTIPLICITY,
                "multiplicity",
                format!(
                    "structure {}: multiplicity {} is impossible with {} electrons (charge {}): the unpaired electrons must not exceed the electron count and must have its parity (spec § 5.1)",
                    quote(sid),
                    num(m),
                    num(electrons),
                    num(charge)
                ),
            ));
        }
    }
    if only_structure {
        if let Some(prov) = as_obj(provenance) {
            for (field, value) in [("charge", charge), ("multiplicity", mult)] {
                if let (Some(v), Some(stated)) = (value, int_of(prov.get(field))) {
                    if stated != v {
                        errors.push((
                            Code::E_ROOT_INCONSISTENT,
                            field,
                            format!(
                                "structure {}: {field} {} differs from provenance.{field} {}, which describes the archive's only structure (spec § 5.1)",
                                quote(sid),
                                num(v),
                                num(stated)
                            ),
                        ));
                    }
                }
            }
        }
    }
    errors
}

// --- reaction.scheme ------------------------------------------------------------------------------

fn reaction_scheme_errors(
    sid: &str,
    payload: Option<&Value>,
    kinds_by_id: Option<&HashMap<String, String>>,
) -> Vec<(Code, String)> {
    let mut errors = Vec::new();
    let Some(steps) = as_obj(payload).and_then(|p| p.get("steps")).and_then(Value::as_array) else { return errors };
    let s = sec(sid);
    let mut seen: HashSet<&str> = HashSet::new();
    for (k, step) in steps.iter().enumerate() {
        let Some(step) = step.as_object() else { continue };
        if let Some(step_id) = str_of(step.get("id")) {
            if seen.contains(step_id) {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    format!("{s}step id {} is used by more than one step (spec § 5.14)", quote(step_id)),
                ));
            }
            seen.insert(step_id);
        }
        for role in ["reactants", "products", "reagents"] {
            let Some(items) = step.get(role).and_then(Value::as_array) else { continue };
            let mut counts: HashMap<&str, usize> = HashMap::new();
            for item in items {
                let Some(r) = as_obj(Some(item)).and_then(|i| str_of(i.get("structure_ref"))) else { continue };
                *counts.entry(r).or_insert(0) += 1;
                if let Some(kinds) = kinds_by_id {
                    if kinds.get(r).map(String::as_str) != Some("structure") {
                        errors.push((
                            Code::E_REF_UNRESOLVED,
                            format!(
                                "{s}step {k} {role} structure_ref {} does not name a structure section (spec § 5.14)",
                                quote(r)
                            ),
                        ));
                    }
                }
            }
            let mut repeated: Vec<&str> = counts.iter().filter(|(_, n)| **n > 1).map(|(r, _)| *r).collect();
            repeated.sort_unstable();
            for r in repeated {
                errors.push((
                    Code::E_MEMBER_MALFORMED,
                    format!("{s}step {k} names {} more than once among its {role} (spec § 5.14)", quote(r)),
                ));
            }
        }
        if let Some(energies) = as_obj(step.get("energies")) {
            errors.extend(reaction_energies_errors(&format!("{s}step {k}"), step, energies, kinds_by_id));
        }
    }
    errors
}

const HARTREE_TO_KCAL_MOL: f64 = 627.509474;
const REACTION_DELTA_TOL_EH: f64 = 1e-8;
const REACTION_KCAL_REL_TOL: f64 = 1e-6;
const REACTION_COEFFICIENT_TOL: f64 = 1e-9;

/// Sum without cancellation loss beyond one rounding, as Python's `math.fsum`.
fn fsum(values: &[f64]) -> f64 {
    let mut partials: Vec<f64> = Vec::new();
    for &v in values {
        let mut x = v;
        let mut i = 0;
        for j in 0..partials.len() {
            let mut y = partials[j];
            if x.abs() < y.abs() {
                std::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            let lo = y - (hi - x);
            if lo != 0.0 {
                partials[i] = lo;
                i += 1;
            }
            x = hi;
        }
        partials.truncate(i);
        partials.push(x);
    }
    partials.iter().sum()
}

fn reaction_energies_errors(
    where_: &str,
    step: &Obj,
    energies: &Obj,
    kinds_by_id: Option<&HashMap<String, String>>,
) -> Vec<(Code, String)> {
    let mut errors: Vec<(Code, String)> = Vec::new();
    let bad = |errors: &mut Vec<(Code, String)>, msg: String| {
        errors.push((Code::E_MEMBER_MALFORMED, format!("{where_} energies: {msg} (spec § 5.14)")));
    };
    if str_of(step.get("arrow")) == Some("resonance") {
        bad(&mut errors, "a resonance arrow has no reaction energy".into());
    }
    let Some(all_terms) = energies.get("terms").and_then(Value::as_array) else { return errors };
    let mut expected: HashMap<String, f64> = HashMap::new();
    for (role, sign) in [("reactants", -1.0), ("products", 1.0)] {
        let Some(items) = step.get(role).and_then(Value::as_array) else { continue };
        for item in items {
            let Some(item) = item.as_object() else { continue };
            let Some(r) = str_of(item.get("structure_ref")) else { continue };
            let c = match item.get("coefficient") {
                None => Some(1.0),
                Some(v) => v.as_f64(),
            };
            if let Some(c) = c {
                *expected.entry(r.to_string()).or_insert(0.0) += sign * c;
            }
        }
    }
    expected.retain(|_, c| *c != 0.0);
    let terms: Vec<&Obj> = all_terms.iter().filter_map(Value::as_object).collect();
    let mut seen: HashSet<&str> = HashSet::new();
    for (t, term) in terms.iter().enumerate() {
        let Some(r) = str_of(term.get("structure_ref")) else { continue };
        if seen.contains(r) {
            bad(&mut errors, format!("structure {} has more than one term", quote(r)));
            continue;
        }
        seen.insert(r);
        match expected.get(r) {
            None => bad(&mut errors, format!("term {t} names {}, which is not a reactant or product", quote(r))),
            Some(&want) => {
                if let Some(c) = num_of(term.get("coefficient")) {
                    if (c - want).abs() > REACTION_COEFFICIENT_TOL * 1f64.max(c.abs()) {
                        bad(
                            &mut errors,
                            format!("term {t} coefficient {} for {}, but the step gives {}", num(c), quote(r), num(want)),
                        );
                    }
                }
            }
        }
        if let (Some(kinds), Some(run_ref)) = (kinds_by_id, str_of(term.get("run_ref"))) {
            if kinds.get(run_ref).map(String::as_str) != Some("run.record") {
                errors.push((
                    Code::E_REF_UNRESOLVED,
                    format!(
                        "{where_} energies: term {t} run_ref {} does not name a run.record section (spec § 5.14)",
                        quote(run_ref)
                    ),
                ));
            }
        }
    }
    let mut missing: Vec<&String> = expected.keys().filter(|r| !seen.contains(r.as_str())).collect();
    missing.sort();
    for r in missing {
        bad(&mut errors, format!("no term for participant {}", quote(r)));
    }
    let methods: Vec<&Value> = terms.iter().filter_map(|t| t.get("method")).filter(|m| m.is_object()).collect();
    if methods.iter().any(|m| !json_equal(m, methods[0])) {
        bad(&mut errors, "the terms have different methods".into());
    }
    let Some(deltas) = as_obj(energies.get("deltas")) else { return errors };
    for (quantity, field) in [("energy", "energy_eh"), ("enthalpy", "enthalpy_eh"), ("gibbs_free_energy", "gibbs_free_energy_eh")]
    {
        let Some(delta) = as_obj(deltas.get(quantity)) else { continue };
        let Some(value) = num_of(delta.get("value_eh")) else { continue };
        if let Some(missing) = terms.iter().position(|t| num_of(t.get(field)).is_none()) {
            bad(&mut errors, format!("{quantity} is given but term {missing} has no {field}"));
        } else if terms.iter().all(|t| num_of(t.get("coefficient")).is_some()) {
            let products: Vec<f64> =
                terms.iter().map(|t| num_of(t.get("coefficient")).unwrap_or(0.0) * num_of(t.get(field)).unwrap_or(0.0)).collect();
            let total = fsum(&products);
            if (total - value).abs() > REACTION_DELTA_TOL_EH {
                bad(&mut errors, format!("{quantity} value_eh {} is not the sum of its terms, {}", num(value), num(total)));
            }
        }
        if let Some(kcal) = num_of(delta.get("kcal_mol")) {
            let converted = value * HARTREE_TO_KCAL_MOL;
            if (kcal - converted).abs() > REACTION_KCAL_REL_TOL * 1f64.max(converted.abs()) {
                bad(&mut errors, format!("{quantity} kcal_mol {} is not value_eh x 627.509474 = {}", num(kcal), num(converted)));
            }
        }
    }
    errors
}

// --- preview -----------------------------------------------------------------------------------

fn preview_member_errors(sid: &str, role: &str, member: &Obj, bytes: &[u8]) -> Vec<(Code, String)> {
    let media = str_of(member.get("media_type"));
    let w = format!("section {} member {}", quote(sid), quote(role));
    match media {
        Some("image/png") if !bytes.starts_with(&[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]) => {
            vec![(Code::E_MEMBER_MALFORMED, format!("{w}: image/png bytes do not begin with the PNG signature (spec § 5.15)"))]
        }
        Some("image/jpeg") if !bytes.starts_with(&[0xff, 0xd8, 0xff]) => {
            vec![(Code::E_MEMBER_MALFORMED, format!("{w}: image/jpeg bytes do not begin with FF D8 FF (spec § 5.15)"))]
        }
        Some("image/svg+xml") => match decode_utf8(bytes) {
            Err(e) => vec![(Code::E_UTF8, format!("{w}: image/svg+xml bytes do not decode as UTF-8 ({e}; spec § 5.15)"))],
            Ok(text) if !text.contains("<svg") => {
                vec![(Code::E_MEMBER_MALFORMED, format!("{w}: image/svg+xml bytes contain no <svg element (spec § 5.15)"))]
            }
            Ok(_) => vec![],
        },
        _ => vec![],
    }
}

// --- the run -----------------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Deferred<'a> {
    i: usize,
    sid: &'a str,
    s: &'a Obj,
}

fn run(zip: &ZipArchive, lim: &Limits, level: Level, schema_check: bool, out: &mut Findings) -> usize {
    let entries = zip.entries();

    // --- container: central-directory checks, before decompressing anything ---
    let names = zip.names();
    let count_of = |name: &str| entries.iter().filter(|e| e.name == name).count();
    for name in &names {
        let n = count_of(name);
        if n > 1 && *name == "manifest.json" {
            // An archive updated in place: the last manifest.json is the manifest (spec § 2.3).
            for copy in 1..n {
                out.push(
                    Code::W_ORPHAN_MEMBER,
                    format!("entry:{name}"),
                    format!(
                        "ZIP entry {} ({copy} of {n}) is superseded by the last {} entry (archive updated in place; consumers ignore it)",
                        quote(name),
                        quote(name)
                    ),
                );
            }
        } else if n > 1 {
            out.push(
                Code::E_ZIP_DUPLICATE_ENTRY,
                format!("entry:{name}"),
                format!("ZIP entry name {} occurs {n} times; readers disagree on which entry wins", quote(name)),
            );
        }
    }
    for name in &names {
        if let Some(why) = path_problem(name) {
            out.push(Code::E_MEMBER_PATH, format!("entry:{name}"), format!("ZIP entry {}: illegal path ({why})", quote(name)));
        }
    }

    let mut unsafe_entries: HashSet<&str> = HashSet::new();
    let mut total: f64 = 0.0;
    for e in entries {
        total += e.size as f64;
        if e.size as f64 > lim.max_member_bytes {
            out.push(
                Code::E_MEMBER_TOO_LARGE,
                format!("entry:{}", e.name),
                format!(
                    "ZIP entry {} exceeds the per-member size cap ({} > {} bytes uncompressed)",
                    quote(&e.name),
                    e.size,
                    num(lim.max_member_bytes)
                ),
            );
            unsafe_entries.insert(&e.name);
        }
        if ratio_exceeded(e, lim) {
            let ratio = if e.compressed_size > 0 {
                format!("{}", (e.size as f64 / e.compressed_size as f64).round())
            } else {
                "inf".into()
            };
            out.push(
                Code::E_COMPRESSION_RATIO,
                format!("entry:{}", e.name),
                format!(
                    "ZIP entry {} expands {} bytes from {} (ratio {ratio} > {})",
                    quote(&e.name),
                    e.size,
                    e.compressed_size,
                    num(lim.max_ratio)
                ),
            );
            unsafe_entries.insert(&e.name);
        }
    }
    let payloads_blocked = total > lim.max_total_bytes;
    if payloads_blocked {
        out.push(
            Code::E_TOTAL_TOO_LARGE,
            "archive",
            format!(
                "archive expands to {} bytes uncompressed, over the {}-byte total cap; member payloads were not read",
                num(total),
                num(lim.max_total_bytes)
            ),
        );
    }
    // --- manifest (its presence is a container check) ---
    let Some(manifest_entry) = zip.find("manifest.json") else {
        out.push(Code::E_MANIFEST_MISSING, "manifest.json", "missing manifest.json");
        return 0;
    };
    if level < Level::Manifest {
        return 0;
    }
    if manifest_entry.size as f64 > lim.max_json_bytes {
        out.push(
            Code::E_JSON_TOO_LARGE,
            "manifest.json",
            format!("manifest.json is {} bytes, over the {}-byte JSON cap", manifest_entry.size, num(lim.max_json_bytes)),
        );
        return 0;
    }
    if unsafe_entries.contains("manifest.json") {
        return 0;
    }
    let raw = match zip.read(manifest_entry, Some(lim.max_member_bytes)) {
        Ok(r) => r,
        Err(e) => {
            out.push(
                Code::E_ZIP_UNREADABLE,
                "entry:manifest.json",
                format!("ZIP entry 'manifest.json' cannot be read: {}", e.message),
            );
            return 0;
        }
    };
    if json_depth(&raw, lim.max_json_depth) > lim.max_json_depth {
        out.push(
            Code::E_JSON_TOO_DEEP,
            "manifest.json",
            format!("manifest.json nests deeper than {} levels", num(lim.max_json_depth)),
        );
        return 0;
    }
    let manifest_value = match parse_json_bytes(&raw) {
        Ok(v) => v,
        Err(e) => {
            out.push(Code::E_MANIFEST_JSON, "manifest.json", format!("manifest.json is not valid JSON: {e}"));
            return 0;
        }
    };
    let Some(manifest) = manifest_value.as_object() else {
        out.push(Code::E_MANIFEST_JSON, "manifest.json", "manifest.json must be a JSON object");
        return 0;
    };

    // Root sanity. The withdrawn version 2 is read as 1 for every check (spec § 3.1).
    let version = manifest.get("qvf_version");
    let mut checked_manifest: Option<Value> = None;
    if num_of(version) == Some(2.0) {
        out.push(
            Code::W_VERSION_WITHDRAWN,
            "/qvf_version",
            "qvf_version 2 was withdrawn: the archive is checked as version 1, which producers must write (spec § 3.1)",
        );
        let mut m = manifest.clone();
        m.insert("qvf_version".into(), Value::from(1));
        checked_manifest = Some(Value::Object(m));
    } else if num_of(version) != Some(1.0) {
        out.push(Code::E_VERSION, "/qvf_version", format!("qvf_version must be 1, got {}", show(version)));
    }
    let empty = Obj::new();
    let mut source: &Obj = &empty;
    if let Some(src) = manifest.get("source") {
        match src.as_object() {
            Some(o) => source = o,
            None => out.push(Code::E_SOURCE, "/source", "source must be a JSON object"),
        }
    }
    for key in ["program", "version", "calculation"] {
        if !source.contains_key(key) {
            out.push(Code::E_SOURCE, "/source", format!("source.{key} is missing"));
        }
    }
    if !truthy(source.get("program")) || !truthy(source.get("version")) {
        out.push(Code::E_SOURCE, "/source", "source.program and source.version must be non-empty");
    }
    if level >= Level::Semantics {
        for (w, msg) in root_consistency_errors(manifest) {
            out.push(Code::E_ROOT_INCONSISTENT, w, msg);
        }
    }

    // The JSON Schema, through the checker generated from it.
    if schema_check {
        for e in check_manifest(checked_manifest.as_ref().unwrap_or(&manifest_value)) {
            let w = if e.pointer.is_empty() { "manifest.json".to_string() } else { e.pointer.clone() };
            let loc = if e.location.is_empty() { "(root)".to_string() } else { e.location.clone() };
            out.push(Code::E_SCHEMA, w, format!("{loc}: {}", e.message));
        }
    } else {
        out.push(Code::W_SCHEMA_SKIPPED, "manifest.json", "the JSON Schema check was disabled; structural schema check skipped");
    }

    // Declared units of {value, units} objects in provenance (spec § 4.3).
    if level >= Level::Semantics {
        if let Some(prov) = as_obj(manifest.get("provenance")) {
            for (key, value) in prov {
                let Some(v) = value.as_object() else { continue };
                if !v.contains_key("value") || !v.contains_key("units") {
                    continue;
                }
                if canonical_unit(v.get("units")).is_none() {
                    out.push(
                        Code::W_UNIT_UNKNOWN,
                        format!("/provenance/{}/units", escape_pointer(key)),
                        format!(
                            "provenance.{key}.units {} is not a unit name of spec § 4.3 or a legacy spelling",
                            show(v.get("units"))
                        ),
                    );
                }
            }
        }
    }

    let no_sections: Vec<Value> = Vec::new();
    let mut sections: &Vec<Value> = &no_sections;
    if let Some(v) = manifest.get("sections") {
        match v.as_array() {
            Some(a) => sections = a,
            None => out.push(Code::E_SECTION_MALFORMED, "/sections", "sections must be an array"),
        }
    }
    let mut section_ids: HashSet<&str> = HashSet::new();
    let mut kinds_by_id: HashMap<String, String> = HashMap::new();
    let mut structure_ids: Vec<&str> = Vec::new();
    let mut section_kinds: HashMap<&str, Option<&str>> = HashMap::new();
    for s in sections {
        let Some(o) = s.as_object() else { continue };
        let Some(id) = str_of(o.get("id")) else { continue };
        section_ids.insert(id);
        let kind = str_of(o.get("kind"));
        if let Some(k) = kind {
            kinds_by_id.entry(id.to_string()).or_insert_with(|| k.to_string());
        }
        if kind == Some("structure") {
            structure_ids.push(id);
        }
        section_kinds.entry(id).or_insert(kind);
    }
    let mut mo_rows_of: HashMap<String, MoRows> = HashMap::new();
    for s in sections {
        let Some(o) = s.as_object() else { continue };
        if str_of(o.get("kind")) == Some("wavefunction.gto") {
            if let Some(id) = str_of(o.get("id")) {
                mo_rows_of.entry(id.to_string()).or_insert_with(|| mo_rows(o.get("members")));
            }
        }
    }
    let mut n_atoms_of: HashMap<String, f64> = HashMap::new();
    let mut sum_z_of: HashMap<String, f64> = HashMap::new();
    let mut cells_of_owned: HashMap<String, usize> = HashMap::new(); // structure id -> section index
    let mut periodic_of: HashMap<String, bool> = HashMap::new();
    let mut state_indices_of: HashMap<String, Option<StateSet>> = HashMap::new();
    let mut identity_sections: Vec<Deferred> = Vec::new();
    let mut nto_sections: Vec<NtoSection> = Vec::new();
    let mut excited_sections: Vec<Deferred> = Vec::new();
    let mut states_ref_spectra: Vec<Deferred> = Vec::new();
    let mut transition_densities: Vec<Deferred> = Vec::new();
    let mut wavefunction_volumes: Vec<(Deferred, &str)> = Vec::new();
    let mut bloch_sections: Vec<Deferred> = Vec::new();
    let mut derivative_sections: Vec<(Deferred, &str)> = Vec::new();
    let mut bonds_sections: Vec<Deferred> = Vec::new();
    let mut grid_members: Vec<(Deferred, String)> = Vec::new(); // (section, role)
    let mut property_sections: Vec<(Deferred, &str)> = Vec::new();
    let declared_extensions = as_obj(manifest.get("extensions")).unwrap_or(&empty);
    let mut declared_paths: HashSet<&str> = HashSet::new();
    let mut seen_ids: HashSet<&str> = HashSet::new();
    let mut all_docs: Vec<Docs> = Vec::with_capacity(sections.len());
    let placeholder_ids: Vec<String> = (0..sections.len()).map(|i| format!("<section {i}>")).collect();

    for (i, sv) in sections.iter().enumerate() {
        let here = section_pointer(i);
        let Some(s) = sv.as_object() else {
            out.push(Code::E_SECTION_MALFORMED, here, format!("section {i} must be a JSON object"));
            all_docs.push(Docs::new());
            continue;
        };
        let sid: &str = match str_of(s.get("id")) {
            None => {
                out.push(
                    Code::E_SECTION_MALFORMED,
                    here.clone(),
                    format!("section {i}: id must be a string, got {}", show(s.get("id"))),
                );
                &placeholder_ids[i]
            }
            Some(id) => {
                if seen_ids.contains(id) {
                    out.push(Code::E_SECTION_ID_DUPLICATE, format!("{here}/id"), format!("duplicate section id {}", quote(id)));
                }
                seen_ids.insert(id);
                if RESERVED_SECTION_IDS.contains(&id) {
                    out.push(
                        Code::E_RESERVED_ID,
                        format!("{here}/id"),
                        format!("section id {} is reserved (a viewer_defaults key, spec § 3.1, § 3.4)", quote(id)),
                    );
                }
                id
            }
        };

        let kind: &str = match str_of(s.get("kind")) {
            None => {
                out.push(
                    Code::E_SECTION_MALFORMED,
                    here.clone(),
                    format!("section {}: kind must be a string, got {}", quote(sid), show(s.get("kind"))),
                );
                ""
            }
            Some(k) => {
                if !is_canonical(k) && !k.starts_with("x_") {
                    out.push(
                        Code::E_KIND_UNKNOWN,
                        format!("{here}/kind"),
                        format!("section {}: unknown non-vendor kind {}", quote(sid), quote(k)),
                    );
                }
                k
            }
        };

        // critical governance (spec § 7.4).
        if truthy(s.get("critical")) && kind.starts_with("x_") {
            let vendor = kind.split('.').next().unwrap_or(kind);
            if !declared_extensions.contains_key(vendor) {
                out.push(
                    Code::E_CRITICAL_UNDECLARED,
                    format!("{here}/critical"),
                    format!(
                        "section {} is critical but its namespace {} is not declared in root extensions",
                        quote(sid),
                        quote(vendor)
                    ),
                );
            }
        }

        let members: &Obj = match s.get("members").and_then(Value::as_object) {
            Some(m) => m,
            None => {
                out.push(
                    Code::E_SECTION_MALFORMED,
                    here.clone(),
                    format!("section {}: members must be a JSON object", quote(sid)),
                );
                &empty
            }
        };

        let mut docs = Docs::new();
        let mut binaries: HashMap<String, Vec<u8>> = HashMap::new();
        for (role, m) in members {
            let mwhere = member_pointer(i, role);
            let Some((m, mpath)) = m.as_object().and_then(|o| str_of(o.get("path")).map(|p| (o, p))) else {
                out.push(
                    Code::E_MEMBER_MALFORMED,
                    mwhere,
                    format!("section {} member {}: malformed spec", quote(sid), quote(role)),
                );
                continue;
            };
            declared_paths.insert(mpath);
            if let Some(why) = path_problem(mpath) {
                out.push(Code::E_MEMBER_PATH, mwhere, format!("member {}: illegal path ({why})", quote(mpath)));
                continue;
            }
            let Some(entry) = zip.find(mpath) else {
                out.push(Code::E_MEMBER_MISSING, mwhere, format!("member {} declared but not in archive", quote(mpath)));
                continue;
            };
            if let Some((code, msg)) = declared_size_problem(m, entry.size, mpath) {
                out.push(code, mwhere.clone(), msg);
            }
            if level < Level::Payloads || payloads_blocked || unsafe_entries.contains(mpath) {
                continue;
            }
            let ctx = PayloadContext { zip, entry, m, mpath, mwhere: &mwhere, sid, kind, role, lim, level, schema_check };
            check_payload(&ctx, &mut docs, &mut binaries, out);
        }

        if level >= Level::Semantics && is_spectrum_1d(kind) {
            for msg in spectrum_payload_errors(sid, docs.get("spectrum")) {
                out.push(Code::E_PAYLOAD_SCHEMA, member_pointer(i, "spectrum"), msg);
            }
        }
        if level >= Level::Semantics && kind == "bands" {
            for (role, msg) in bands_shape_errors(sid, members) {
                out.push(Code::E_MEMBER_MALFORMED, member_pointer(i, &role), msg);
            }
        }
        if level >= Level::Semantics && is_trajectory_kind(kind) {
            for (role, msg) in trajectory_errors(sid, kind, members, docs.get("metadata"), binaries.get("energies")) {
                out.push(Code::E_MEMBER_MALFORMED, member_pointer(i, &role), msg);
            }
        }
        let d = Deferred { i, sid, s };
        if kind == "structure" {
            let doc = docs.get("structure");
            if let Some(atoms) = as_obj(doc).and_then(|o| o.get("atoms")).and_then(Value::as_array) {
                n_atoms_of.entry(sid.to_string()).or_insert(atoms.len() as f64);
            }
            if !sum_z_of.contains_key(sid) {
                if let Some(z) = sum_atomic_numbers(doc) {
                    sum_z_of.insert(sid.to_string(), z);
                }
            }
            if let Some(o) = as_obj(doc) {
                cells_of_owned.entry(sid.to_string()).or_insert(i);
                periodic_of.entry(sid.to_string()).or_insert_with(|| {
                    o.get("pbc").and_then(Value::as_array).is_some_and(|a| a.iter().any(|x| x == &Value::Bool(true)))
                });
            }
            if level >= Level::Semantics {
                identity_sections.push(d);
            }
        }
        if level >= Level::Semantics && (kind == "forces" || kind == "hessian") {
            derivative_sections.push((d, kind));
        }
        if level >= Level::References && kind == "bonds" && s.contains_key("structure_ref") {
            bonds_sections.push(d);
        }
        if level >= Level::Semantics {
            let kinds = if level >= Level::References { Some(&kinds_by_id) } else { None };
            for (code, role, msg) in volume_errors(sid, kind, members, &docs, kinds) {
                out.push(code, member_pointer(i, &role), msg);
            }
            let mut roles: Vec<&String> = docs.keys().filter(|r| is_grid_payload(kind, r)).collect();
            roles.sort();
            for r in roles {
                if let Some(msg) = as_obj(docs.get(r)).and_then(|g| grid_data_unit_error(sid, kind, g)) {
                    out.push(Code::E_UNIT_CONTRADICTED, member_pointer(i, r), msg);
                }
            }
        }
        if level >= Level::References {
            let mut roles: Vec<&String> = docs.keys().filter(|r| is_grid_payload(kind, r) && docs[*r].is_object()).collect();
            roles.sort();
            for r in roles {
                grid_members.push((d, r.clone()));
            }
        }
        if level >= Level::Semantics && is_property_kind(kind) && docs.contains_key("properties") {
            property_sections.push((d, kind));
        }
        if level >= Level::Semantics && kind == "wavefunction.bloch" {
            bloch_sections.push(d);
        }
        if kind == "excited_states" {
            state_indices_of.entry(sid.to_string()).or_insert_with(|| state_indices(docs.get("states")));
            if level >= Level::Semantics {
                excited_sections.push(d);
            }
        }
        if level >= Level::References
            && is_spectrum_1d(kind)
            && as_obj(docs.get("spectrum")).is_some_and(|p| p.contains_key("states_ref"))
        {
            states_ref_spectra.push(d);
        }
        if level >= Level::References && kind == "volume.transition_density" {
            transition_densities.push(d);
        }
        if level >= Level::References && is_wavefunction_volume_kind(kind) && s.contains_key("wavefunction_ref") {
            wavefunction_volumes.push((d, kind));
        }
        if level >= Level::Semantics && kind == "atom_properties" {
            for (role, msg) in atom_properties_meta_errors(sid, members, &docs) {
                out.push(Code::E_MEMBER_MALFORMED, member_pointer(i, &role), msg);
            }
        }
        if level >= Level::Semantics && kind == "wavefunction.gto" {
            nto_sections.push(NtoSection { i });
            for msg in unrestricted_nto_errors(docs.get("mo_metadata"), members) {
                out.push(Code::E_MEMBER_MALFORMED, member_pointer(i, "mo_metadata"), msg);
            }
            for (code, role, msg) in wavefunction_errors(sid, members, &docs, manifest.get("provenance")) {
                out.push(code, member_pointer(i, &role), msg);
            }
            for msg in normalization_errors(sid, docs.get("basis"), members.get("overlap"), binaries.get("overlap")) {
                out.push(Code::E_NORMALIZATION, member_pointer(i, "overlap"), msg);
            }
        }
        if level >= Level::Semantics && kind == "reaction.scheme" {
            let kinds = if level >= Level::References { Some(&kinds_by_id) } else { None };
            for (code, msg) in reaction_scheme_errors(sid, docs.get("scheme"), kinds) {
                out.push(code, member_pointer(i, "scheme"), msg);
            }
        }

        if level >= Level::References {
            cross_references(i, sid, kind, s, members, &docs, &section_ids, &section_kinds, &kinds_by_id, out);
        }
        all_docs.push(docs);
    }

    let only_structure = if structure_ids.len() == 1 { Some(structure_ids[0]) } else { None };
    let provenance = manifest.get("provenance");

    // Chemical identity of the structures (spec § 5.1).
    for d in &identity_sections {
        for (code, field, msg) in identity_errors(d.sid, d.s, sum_z_of.get(d.sid).copied(), provenance, structure_ids.len() == 1)
        {
            out.push(code, format!("{}/{field}", section_pointer(d.i)), msg);
        }
    }

    let structures = StructureIndex { kinds_by_id: &kinds_by_id, n_atoms_of: &n_atoms_of, only_structure };
    let structures_opt = if level >= Level::References { Some(&structures) } else { None };
    for (d, kind) in &derivative_sections {
        let members = d.s.get("members").and_then(Value::as_object).unwrap_or(&empty);
        for (code, role, msg) in derivative_errors(d.sid, kind, d.s, members, all_docs[d.i].get("meta"), structures_opt) {
            out.push(code, member_pointer(d.i, &role), msg);
        }
    }
    if let Some(st) = structures_opt {
        for d in &bonds_sections {
            for msg in bonds_errors(d.sid, d.s, all_docs[d.i].get("bonds"), st) {
                out.push(Code::E_MEMBER_MALFORMED, member_pointer(d.i, "bonds"), msg);
            }
        }
    }
    if level >= Level::References {
        let ntos: Vec<(Option<&Value>, &Obj)> = nto_sections
            .iter()
            .map(|n| {
                let m = sections[n.i].get("members").and_then(Value::as_object).unwrap_or(&empty);
                (all_docs[n.i].get("mo_metadata"), m)
            })
            .collect();
        for (k, n) in nto_sections.iter().enumerate() {
            let (meta, members) = ntos[k];
            for (code, msg) in nto_pair_errors(meta, members, &ntos, &kinds_by_id, &state_indices_of) {
                out.push(code, member_pointer(n.i, "mo_metadata"), msg);
            }
        }
    }
    let excited_refs = ExcitedRefs { kinds_by_id: &kinds_by_id, mo_rows_of: &mo_rows_of };
    let excited_opt = if level >= Level::References { Some(&excited_refs) } else { None };
    for d in &excited_sections {
        for (code, field, msg) in excited_states_errors(d.sid, d.s, all_docs[d.i].get("states"), excited_opt) {
            let w = if field.is_empty() { member_pointer(d.i, "states") } else { format!("{}/{field}", section_pointer(d.i)) };
            out.push(code, w, msg);
        }
    }
    for d in &states_ref_spectra {
        if let Some(p) = as_obj(all_docs[d.i].get("spectrum")) {
            for msg in spectrum_states_errors(d.sid, p, &kinds_by_id, &state_indices_of) {
                out.push(Code::E_REF_UNRESOLVED, member_pointer(d.i, "spectrum"), msg);
            }
        }
    }
    for d in &transition_densities {
        for (field, msg) in transition_density_errors(d.sid, d.s, &kinds_by_id, &state_indices_of) {
            out.push(Code::E_REF_UNRESOLVED, format!("{}/{field}", section_pointer(d.i)), msg);
        }
    }
    for (d, kind) in &wavefunction_volumes {
        for (field, msg) in wavefunction_volume_errors(d.sid, kind, d.s, &kinds_by_id, &mo_rows_of) {
            out.push(Code::E_REF_UNRESOLVED, format!("{}/{field}", section_pointer(d.i)), msg);
        }
    }

    // Grids: references and periodic spans against the structure (spec § 5.2).
    if let Some(st) = structures_opt {
        let cells_of: HashMap<String, &Obj> =
            cells_of_owned.iter().filter_map(|(id, &i)| as_obj(all_docs[i].get("structure")).map(|o| (id.clone(), o))).collect();
        for (d, role) in &grid_members {
            if let Some(grid) = as_obj(all_docs[d.i].get(role)) {
                for (code, msg) in grid_errors(d.sid, d.s, grid, &kinds_by_id, &cells_of, st.only_structure) {
                    out.push(code, member_pointer(d.i, role), msg);
                }
            }
        }
    }
    for (d, kind) in &property_sections {
        for (code, msg) in property_errors(d.sid, kind, d.s, all_docs[d.i].get("properties"), manifest, structures_opt) {
            out.push(code, member_pointer(d.i, "properties"), msg);
        }
    }
    let periodic = PeriodicIndex { kinds_by_id: &kinds_by_id, periodic_of: &periodic_of };
    let periodic_opt = if level >= Level::References { Some(&periodic) } else { None };
    for d in &bloch_sections {
        let members = d.s.get("members").and_then(Value::as_object).unwrap_or(&empty);
        for (code, w, msg) in bloch_errors(d.sid, d.s, members, &all_docs[d.i], periodic_opt) {
            out.push(code, format!("{}/{w}", section_pointer(d.i)), msg);
        }
    }

    // --- orphans: entries no manifest member references (spec § 2.1) ---
    for name in &names {
        if *name != "manifest.json" && !name.ends_with('/') && !declared_paths.contains(name) {
            out.push(
                Code::W_ORPHAN_MEMBER,
                format!("entry:{name}"),
                format!("ZIP entry {} is not referenced by any manifest member (orphan; consumers ignore it)", quote(name)),
            );
        }
    }
    sections.len()
}

#[allow(clippy::too_many_arguments)]
fn cross_references(
    i: usize,
    sid: &str,
    kind: &str,
    s: &Obj,
    members: &Obj,
    docs: &Docs,
    section_ids: &HashSet<&str>,
    section_kinds: &HashMap<&str, Option<&str>>,
    kinds_by_id: &HashMap<String, String>,
    out: &mut Findings,
) {
    let here = section_pointer(i);
    fn resolves_in<'v>(r: Option<&'v Value>, ids: &HashSet<&str>) -> Option<&'v str> {
        str_of(r).filter(|x| ids.contains(x))
    }
    let resolves = |r| resolves_in(r, section_ids);
    let kind_shown = |r: &str| match section_kinds.get(r).copied().flatten() {
        Some(k) => format!("\"{k}\""),
        None => "nothing".to_string(),
    };
    if kind == "volume.difference" {
        for op in ["operand_a", "operand_b"] {
            let r = s.get(op);
            if !matches!(r, None | Some(Value::Null)) && resolves(r).is_none() {
                out.push(
                    Code::E_REF_UNRESOLVED,
                    format!("{here}/{op}"),
                    format!("section {}: {op} {} does not resolve", quote(sid), show(r)),
                );
            }
        }
    }
    if kind == "reaction.waypoints" {
        let r = s.get("trajectory_ref");
        match resolves(r) {
            None => out.push(
                Code::E_REF_UNRESOLVED,
                format!("{here}/trajectory_ref"),
                format!("section {}: trajectory_ref {} does not resolve", quote(sid), show(r)),
            ),
            Some(rs) if section_kinds.get(rs).copied().flatten() != Some("trajectory") => out.push(
                Code::E_REF_UNRESOLVED,
                format!("{here}/trajectory_ref"),
                format!(
                    "section {}: trajectory_ref {} names a section of kind {}, not 'trajectory' (spec § 5.5)",
                    quote(sid),
                    show(r),
                    kind_shown(rs)
                ),
            ),
            Some(_) => {}
        }
    }
    // Generic section references, on any kind (spec § 4.1).
    for (field, target) in SECTION_REFS {
        let r = s.get(field);
        if matches!(r, None | Some(Value::Null)) {
            continue;
        }
        match resolves(r) {
            None => out.push(
                Code::E_REF_UNRESOLVED,
                format!("{here}/{field}"),
                format!("section {}: {field} {} does not resolve", quote(sid), show(r)),
            ),
            Some(rs) if section_kinds.get(rs).copied().flatten() != Some(target) => out.push(
                Code::E_REF_UNRESOLVED,
                format!("{here}/{field}"),
                format!(
                    "section {}: {field} {} names a section of kind {}, not {} (spec § 4.1)",
                    quote(sid),
                    show(r),
                    kind_shown(rs),
                    quote(target)
                ),
            ),
            Some(_) => {}
        }
    }
    let names_kind = |r: Option<&Value>, want: &str| str_of(r).and_then(|x| kinds_by_id.get(x)).map(String::as_str) == Some(want);
    if kind == "vibrations" && s.contains_key("hessian_ref") && !names_kind(s.get("hessian_ref"), "hessian") {
        out.push(
            Code::E_REF_UNRESOLVED,
            format!("{here}/hessian_ref"),
            format!(
                "section {}: hessian_ref {} does not name a hessian section (spec § 5.5)",
                quote(sid),
                show(s.get("hessian_ref"))
            ),
        );
    }
    if kind == "job.spec" {
        if let Some(spec) = as_obj(docs.get("spec")) {
            if spec.contains_key("structure_ref") && !names_kind(spec.get("structure_ref"), "structure") {
                out.push(
                    Code::E_REF_UNRESOLVED,
                    member_pointer(i, "spec"),
                    format!(
                        "section {}: spec.structure_ref {} does not name a structure section (spec § 5.9)",
                        quote(sid),
                        show(spec.get("structure_ref"))
                    ),
                );
            }
        }
    }
    if kind == "run.record" {
        if s.contains_key("job_ref") && !names_kind(s.get("job_ref"), "job.spec") {
            out.push(
                Code::E_REF_UNRESOLVED,
                format!("{here}/job_ref"),
                format!(
                    "section {}: job_ref {} does not name a job.spec section (spec § 5.8)",
                    quote(sid),
                    show(s.get("job_ref"))
                ),
            );
        }
        if let Some(files) = as_obj(docs.get("files")) {
            for frole in files.keys() {
                if !members.contains_key(frole) {
                    out.push(
                        Code::E_FILES_INDEX,
                        member_pointer(i, "files"),
                        format!("section {}: files index names role {} with no matching member", quote(sid), quote(frole)),
                    );
                }
            }
        }
    }
}

struct PayloadContext<'a> {
    zip: &'a ZipArchive,
    entry: &'a ZipEntry,
    m: &'a Obj,
    mpath: &'a str,
    mwhere: &'a str,
    sid: &'a str,
    kind: &'a str,
    role: &'a str,
    lim: &'a Limits,
    level: Level,
    schema_check: bool,
}

/// Read one member and run the payload-level checks on its bytes.
fn check_payload(c: &PayloadContext, docs: &mut Docs, binaries: &mut HashMap<String, Vec<u8>>, out: &mut Findings) {
    let (m, mpath, mwhere, kind, role, lim) = (c.m, c.mpath, c.mwhere, c.kind, c.role, c.lim);
    let fmt = str_of(m.get("format")).unwrap_or("");
    if fmt == "json" && c.entry.size as f64 > lim.max_json_bytes {
        out.push(
            Code::E_JSON_TOO_LARGE,
            mwhere,
            format!("member {}: {} bytes, over the {}-byte JSON cap", quote(mpath), c.entry.size, num(lim.max_json_bytes)),
        );
        return;
    }
    let bytes = match c.zip.read(c.entry, Some(lim.max_member_bytes)) {
        Ok(b) => b,
        Err(e) => {
            out.push(
                Code::E_ZIP_UNREADABLE,
                format!("entry:{mpath}"),
                format!("ZIP entry {} cannot be read: {}", quote(mpath), e.message),
            );
            return;
        }
    };
    if let Some(declared) = m.get("sha256") {
        if declared.as_str() != Some(sha256_hex(&bytes).as_str()) {
            out.push(Code::E_SHA256_MISMATCH, mwhere, format!("member {}: sha256 mismatch", quote(mpath)));
        }
    }

    if fmt == "json" {
        if json_depth(&bytes, lim.max_json_depth) > lim.max_json_depth {
            out.push(
                Code::E_JSON_TOO_DEEP,
                mwhere,
                format!("member {}: JSON nests deeper than {} levels", quote(mpath), num(lim.max_json_depth)),
            );
            return;
        }
        match parse_json_bytes(&bytes) {
            Err(e) => out.push(Code::E_JSON_INVALID, mwhere, format!("member {}: not valid JSON ({e})", quote(mpath))),
            Ok(doc) => {
                if c.level >= Level::Semantics && kind == "structure" && role == "structure" {
                    for msg in structure_payload_errors(c.sid, &doc) {
                        out.push(Code::E_STRUCTURE_PERIODICITY, mwhere, msg);
                    }
                }
                if c.schema_check {
                    if let Some(def) = payload_definition(kind, role) {
                        for e in check_definition(def, &doc) {
                            out.push(Code::E_PAYLOAD_SCHEMA, mwhere, format!("schema: {mpath}: {}: {}", e.location, e.message));
                        }
                    }
                }
                if c.level >= Level::Semantics && is_canonical(kind) {
                    for (loc, unit) in declared_units(kind, role, &doc) {
                        if !is_unit_name(Some(unit)) {
                            out.push(
                                Code::W_UNIT_UNKNOWN,
                                mwhere,
                                format!(
                                    "member {}: declared unit {} at {loc} is not a unit name of spec § 4.3",
                                    quote(mpath),
                                    show(Some(unit))
                                ),
                            );
                        }
                    }
                }
                docs.insert(role.to_string(), doc);
            }
        }
    } else if fmt == "binary" && m.contains_key("dtype") {
        let dtype = m.get("dtype");
        let itemsize = str_of(dtype).and_then(DType::parse).map(DType::itemsize);
        let shape = m.get("shape");
        let parsed_shape: Option<Vec<f64>> = match shape {
            None => Some(Vec::new()),
            Some(v) => {
                v.as_array().and_then(|a| a.iter().map(|d| d.as_f64().filter(|x| x.fract() == 0.0 && *x >= 0.0)).collect())
            }
        };
        match (itemsize, parsed_shape) {
            (None, _) => out.push(Code::E_DTYPE, mwhere, format!("member {}: bad dtype {}", quote(mpath), show(dtype))),
            (Some(_), None) => out.push(
                Code::E_MEMBER_MALFORMED,
                mwhere,
                format!("member {}: shape must be a list of non-negative integers, got {}", quote(mpath), show(shape)),
            ),
            (Some(size), Some(dims)) => {
                // Exact for any shape: compare in 128-bit integers, saturating.
                let want = dims.iter().fold(size as u128, |n, d| n.saturating_mul(*d as u128));
                if bytes.len() as u128 != want {
                    out.push(
                        Code::E_BINARY_SIZE,
                        mwhere,
                        format!("member {}: byte length {} != itemsize*prod(shape) {want}", quote(mpath), bytes.len()),
                    );
                } else if kind == "wavefunction.gto" && role == "overlap" {
                    binaries.insert(role.to_string(), bytes.clone());
                }
            }
        }
    }

    if is_trajectory_kind(kind) && role == "energies" && fmt == "binary" {
        binaries.insert(role.to_string(), bytes.clone());
    }
    if kind == "run.record" && (role == "input" || role == "log") && decode_utf8(&bytes).is_err() {
        out.push(Code::E_UTF8, mwhere, format!("member {}: run.record {role} must decode as UTF-8", quote(mpath)));
    }
    if kind == "preview" && fmt == "binary" {
        for (code, msg) in preview_member_errors(c.sid, role, m, &bytes) {
            out.push(code, mwhere, msg);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fsum_is_exact_on_cancellation() {
        assert_eq!(fsum(&[1e16, 1.0, -1e16]), 1.0);
    }

    #[test]
    fn cartesian_d_factors() {
        let comps = shell_components(2, false);
        let factors: Vec<f64> = comps.iter().map(|c| c.1).collect();
        let third = 1.0 / 3.0;
        assert_eq!(factors, vec![1.0, third, third, 1.0, third, 1.0]);
    }

    #[test]
    fn not_a_zip() {
        let r = validate_bytes(b"not a zip at all, definitely not".to_vec(), &ValidateOptions::default());
        assert!(!r.ok);
        assert!(r.has(Code::E_ZIP_UNREADABLE));
    }

    #[test]
    fn unit_vocabulary() {
        assert!(is_unit_name(Some(&Value::from("hartree"))));
        assert!(!is_unit_name(Some(&Value::from("Eh"))));
        assert_eq!(canonical_unit(Some(&Value::from("Eh"))), Some("hartree"));
    }
}
