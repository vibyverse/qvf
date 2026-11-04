//! Section builders for trajectories, band structures, densities of states
//! and spectra, the continuation of [`crate::sections`] (and re-exported
//! there).
//!
//! * [`Trajectory`]: `trajectory` (spec § 5.5): frames in Å, optionally
//!   periodic with a fixed or per-frame cell (rows in Å, the structure
//!   convention of § 4.3), velocities, gradients, times, energies, cell
//!   volumes and temperatures, coordinates in float64 or float32.
//! * [`Bands`]: `bands` (§ 5.6) with its [`Kpath`], fat-band projections and
//!   occupations.
//! * [`DosTotal`], [`DosProjected`] and [`DosBonding`]: `dos.total`,
//!   `dos.projected`, `dos.coop` and `dos.cohp` (§ 5.6), each with or
//!   without the spin axis.
//! * [`Spectrum`]: the six one-dimensional `spectra.*` kinds on
//!   `$defs/SpectrumPayload` (§ 5.4). `spectra.nmr` and `spectra.epr` take
//!   their payload as JSON through [`Writer::add_spectrum_payload`].
//!
//! As in [`crate::sections`], everything that can be checked is checked
//! before anything is written: extents, the spin axis, unit names,
//! references to sections added earlier, and each JSON payload against its
//! schema definition.

use std::io::Write;

use serde_json::{json, Map, Value};

use crate::generated;
use crate::reader::DType;
use crate::schema_check::check_definition;
use crate::sections::ELEMENTS;
use crate::writer::{invalid, WriteError, Writer};

fn finite(what: &str, values: &[f64]) -> Result<(), WriteError> {
    match values.iter().position(|x| !x.is_finite()) {
        Some(i) => invalid(format!("{what}[{i}] is not finite")),
        None => Ok(()),
    }
}

fn unit_name(what: &str, unit: &str) -> Result<(), WriteError> {
    if generated::UNIT_NAMES.contains(&unit) {
        Ok(())
    } else {
        invalid(format!("{what}: {unit:?} is not a unit name of spec § 4.3"))
    }
}

/// Refuse a payload that does not fit its schema definition.
fn schema(what: &str, def: &str, payload: &Value) -> Result<(), WriteError> {
    match check_definition(def, payload).first() {
        Some(e) => {
            invalid(format!("{what}: the payload does not fit $defs/{def} at {:?}: {} (E_PAYLOAD_SCHEMA)", e.location, e.message))
        }
        None => Ok(()),
    }
}

/// Insert the optional keys that are set.
fn put(m: &mut Map<String, Value>, pairs: Vec<(&str, Option<Value>)>) {
    for (k, v) in pairs {
        if let Some(v) = v {
            m.insert(k.to_owned(), v);
        }
    }
}

/// Merge the caller's extra keys, refusing those the builder writes.
fn merge_extra(what: &str, m: &mut Map<String, Value>, extra: &Map<String, Value>) -> Result<(), WriteError> {
    for (k, v) in extra {
        if m.contains_key(k) {
            return invalid(format!("{what}: extra key {k:?} is written by the builder"));
        }
        m.insert(k.clone(), v.clone());
    }
    Ok(())
}

// -- trajectory --------------------------------------------------------------

/// The cell of a periodic trajectory, rows in Å (the structure convention,
/// spec § 4.3, *Lattice conventions*).
#[derive(Debug, Clone, PartialEq)]
pub enum TrajectoryCell {
    /// One cell for every frame, written `[3, 3]`.
    Fixed([[f64; 3]; 3]),
    /// One cell per frame (variable cell, NPT), written `[n_frames, 3, 3]`.
    PerFrame(Vec<[[f64; 3]; 3]>),
}

/// The precision of `coords` (spec § 5.5). Every other member is float64.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Precision {
    #[default]
    Float64,
    /// Halves a long run; sets `metadata.precision` to `"float32"`.
    Float32,
}

/// A `trajectory` section (spec § 5.5).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Trajectory {
    pub id: String,
    /// Element symbol of each atom, in the atom order of `coords`. The
    /// metadata names the atoms by `symbols` and `atomic_numbers`.
    pub symbols: Vec<String>,
    /// `[n_frames, n_atoms, 3]` in Å, C order.
    pub coords: Vec<f64>,
    pub precision: Precision,
    /// `optimization`, `irc`, `md`, `neb`, `scan` or `custom`.
    pub kind: Option<String>,
    /// Required with `cell`, and a periodic axis requires `cell`.
    pub pbc: Option<[bool; 3]>,
    pub cell: Option<TrajectoryCell>,
    /// `[n_frames, n_atoms, 3]` in Å/fs.
    pub velocities: Option<Vec<f64>>,
    /// The energy gradient dE/dR `[n_frames, n_atoms, 3]` in hartree/bohr
    /// (the forces are its negative).
    pub gradient: Option<Vec<f64>>,
    /// `[n_frames]` in fs.
    pub times: Option<Vec<f64>>,
    /// `[n_frames]` in hartree, the binary `energies` member.
    pub energies: Option<Vec<f64>>,
    /// Also write the energies as `metadata.energies` (identical values).
    pub energies_in_metadata: bool,
    /// `[n_frames]` in Å³.
    pub cell_volumes: Option<Vec<f64>>,
    /// `[n_frames]` in kelvin.
    pub temperatures: Option<Vec<f64>>,
    /// `nve`, `nvt`, `npt`, `nph` or `custom`.
    pub ensemble: Option<String>,
    /// Time between frames in fs (> 0).
    pub timestep_fs: Option<f64>,
    /// The metadata `units` declaration (spec § 4.3), restating the fixed units.
    pub units: Option<Map<String, Value>>,
    /// Further metadata keys (typed-open payload).
    pub extra: Map<String, Value>,
    pub label: Option<String>,
}

impl Trajectory {
    /// A trajectory of the atoms `symbols` with frames `coords`
    /// (`[n_frames, n_atoms, 3]` in Å).
    pub fn new(id: &str, symbols: Vec<String>, coords: Vec<f64>) -> Self {
        Trajectory { id: id.to_owned(), symbols, coords, ..Trajectory::default() }
    }

    /// A molecular-dynamics run (`kind` `md`) with a time step.
    pub fn md(id: &str, symbols: Vec<String>, coords: Vec<f64>, timestep_fs: f64) -> Self {
        let mut t = Self::new(id, symbols, coords);
        t.kind = Some("md".into());
        t.timestep_fs = Some(timestep_fs);
        t
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.symbols.len()
    }

    /// Number of frames (0 when `coords` does not fill whole frames).
    pub fn n_frames(&self) -> usize {
        let per = 3 * self.n_atoms();
        if per == 0 || self.coords.len() % per != 0 {
            0
        } else {
            self.coords.len() / per
        }
    }

    /// The `metadata` member payload.
    pub fn metadata_payload(&self) -> Result<Value, WriteError> {
        let what = format!("trajectory {:?}", self.id);
        let mut m = Map::new();
        let z: Vec<u64> = self.symbols.iter().map(|s| ELEMENTS.iter().position(|e| e == s).unwrap_or(0) as u64).collect();
        m.insert("symbols".into(), json!(self.symbols));
        m.insert("atomic_numbers".into(), json!(z));
        m.insert("n_frames".into(), Value::from(self.n_frames()));
        put(
            &mut m,
            vec![
                ("kind", self.kind.clone().map(Value::from)),
                ("pbc", self.pbc.map(|p| json!(p))),
                ("ensemble", self.ensemble.clone().map(Value::from)),
                ("timestep_fs", self.timestep_fs.map(Value::from)),
                ("units", self.units.clone().map(Value::Object)),
            ],
        );
        if self.precision == Precision::Float32 {
            m.insert("precision".into(), Value::from("float32"));
        }
        if self.energies_in_metadata {
            match &self.energies {
                Some(e) => {
                    m.insert("energies".into(), json!(e));
                }
                None => return invalid(format!("{what}: energies_in_metadata needs energies")),
            }
        }
        merge_extra(&what, &mut m, &self.extra)?;
        Ok(Value::Object(m))
    }

    fn check(&self) -> Result<(usize, usize), WriteError> {
        let what = format!("trajectory {:?}", self.id);
        let n_atoms = self.n_atoms();
        if n_atoms == 0 {
            return invalid(format!("{what}: needs at least one atom"));
        }
        if let Some(s) = self.symbols.iter().find(|s| s.is_empty() || s.chars().count() > 3) {
            return invalid(format!("{what}: {s:?} is not an element symbol"));
        }
        if self.coords.len() % (3 * n_atoms) != 0 {
            return invalid(format!("{what}: {} coords are not [n_frames, {n_atoms}, 3]", self.coords.len()));
        }
        let n_frames = self.n_frames();
        finite(&format!("{what} coords"), &self.coords)?;
        if self.precision == Precision::Float32 && self.coords.iter().any(|&x| !(x as f32).is_finite()) {
            return invalid(format!("{what}: coords overflow float32"));
        }
        for (role, v) in [("velocities", &self.velocities), ("gradient", &self.gradient)] {
            if let Some(v) = v {
                if v.len() != self.coords.len() {
                    return invalid(format!(
                        "{what}: {role} must be [n_frames, n_atoms, 3] = [{n_frames}, {n_atoms}, 3] (spec § 5.5)"
                    ));
                }
                finite(&format!("{what} {role}"), v)?;
            }
        }
        for (role, v) in [
            ("times", &self.times),
            ("energies", &self.energies),
            ("cell_volumes", &self.cell_volumes),
            ("temperatures", &self.temperatures),
        ] {
            if let Some(v) = v {
                if v.len() != n_frames {
                    return invalid(format!(
                        "{what}: {role} must be [n_frames] = [{n_frames}], got {} values (spec § 5.5)",
                        v.len()
                    ));
                }
                finite(&format!("{what} {role}"), v)?;
            }
        }
        match &self.cell {
            Some(cell) => {
                if self.pbc.is_none() {
                    return invalid(format!("{what}: a cell needs pbc (spec § 5.5)"));
                }
                let cells: Vec<&[[f64; 3]; 3]> = match cell {
                    TrajectoryCell::Fixed(c) => vec![c],
                    TrajectoryCell::PerFrame(list) => {
                        if list.len() != n_frames {
                            return invalid(format!("{what}: a per-frame cell needs {n_frames} cells, got {}", list.len()));
                        }
                        list.iter().collect()
                    }
                };
                for (k, c) in cells.iter().enumerate() {
                    finite(&format!("{what} cell {k}"), &c.concat())?;
                }
            }
            None => {
                if self.pbc.is_some_and(|p| p.iter().any(|&x| x)) {
                    return invalid(format!("{what}: a periodic axis in pbc needs a cell (spec § 5.5)"));
                }
            }
        }
        if let Some(dt) = self.timestep_fs {
            if !(dt.is_finite() && dt > 0.0) {
                return invalid(format!("{what}: timestep_fs must be positive"));
            }
        }
        if let Some(units) = &self.units {
            for (k, v) in units {
                unit_name(&format!("{what}: units.{k}"), v.as_str().unwrap_or_default())?;
            }
        }
        Ok((n_frames, n_atoms))
    }
}

// -- bands -------------------------------------------------------------------

/// One projection channel (`$defs/DosChannel`) of fat bands or a projected
/// DOS.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DosChannel {
    pub label: String,
    /// Atom index into the structure, from 0. `None`: summed over atoms.
    pub atom: Option<usize>,
    /// Shell angular momentum. `None`: summed over shells.
    pub l: Option<u32>,
    /// Component −l..l (Appendix A.2). `None`: summed over the shell.
    pub m: Option<i32>,
    /// 0 (α) or 1 (β) for a channel that is itself spin-resolved.
    pub spin: Option<u8>,
    /// Further keys, such as `symbol`.
    pub extra: Map<String, Value>,
}

impl DosChannel {
    pub fn new(label: &str) -> Self {
        DosChannel { label: label.to_owned(), ..DosChannel::default() }
    }

    /// A channel of one shell of one atom.
    pub fn shell(label: &str, atom: usize, l: u32) -> Self {
        DosChannel { label: label.to_owned(), atom: Some(atom), l: Some(l), ..DosChannel::default() }
    }

    fn payload(&self, what: &str) -> Result<Value, WriteError> {
        if let (Some(l), Some(m)) = (self.l, self.m) {
            if m.unsigned_abs() > l {
                return invalid(format!("{what}: channel {:?} has m = {m} outside -l..l for l = {l}", self.label));
            }
        }
        if self.spin.is_some_and(|s| s > 1) {
            return invalid(format!("{what}: channel {:?} spin must be 0 or 1", self.label));
        }
        let mut o = Map::new();
        o.insert("label".into(), Value::from(self.label.as_str()));
        put(
            &mut o,
            vec![
                ("atom", self.atom.map(Value::from)),
                ("l", self.l.map(Value::from)),
                ("m", self.m.map(Value::from)),
                ("spin", self.spin.map(Value::from)),
            ],
        );
        merge_extra(what, &mut o, &self.extra)?;
        Ok(Value::Object(o))
    }
}

/// A straight segment of a k-path (`$defs/KpathSegment`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct KpathSegment {
    pub start: Option<usize>,
    pub end: Option<usize>,
    pub n_points: Option<usize>,
    pub label_start: Option<String>,
    pub label_end: Option<String>,
    pub k_start: Option<[f64; 3]>,
    pub k_end: Option<[f64; 3]>,
}

impl KpathSegment {
    /// A segment of `n_points` points between two labelled points.
    pub fn between(label_start: &str, k_start: [f64; 3], label_end: &str, k_end: [f64; 3], n_points: usize) -> Self {
        KpathSegment {
            n_points: Some(n_points),
            label_start: Some(label_start.to_owned()),
            label_end: Some(label_end.to_owned()),
            k_start: Some(k_start),
            k_end: Some(k_end),
            ..KpathSegment::default()
        }
    }

    fn payload(&self) -> Value {
        let mut o = Map::new();
        put(
            &mut o,
            vec![
                ("start", self.start.map(Value::from)),
                ("end", self.end.map(Value::from)),
                ("n_points", self.n_points.map(Value::from)),
                ("label_start", self.label_start.clone().map(Value::from)),
                ("label_end", self.label_end.clone().map(Value::from)),
                ("k_start", self.k_start.map(|k| json!(k))),
                ("k_end", self.k_end.map(|k| json!(k))),
            ],
        );
        Value::Object(o)
    }
}

/// The `kpath` payload of a `bands` section (`$defs/KpathPayload`, strict).
/// The builder writes `n_spin`, `n_kpoints` and `n_bands` from the
/// eigenvalues.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Kpath {
    /// Fractional reciprocal-lattice coordinates, one per k-point (or empty).
    pub kpoints: Vec<[f64; 3]>,
    /// One label per k-point (`""` unlabelled), or empty.
    pub labels: Vec<String>,
    pub segments: Vec<KpathSegment>,
    /// Fermi level in eV on the scale of the eigenvalues, written as both
    /// `fermi` and `fermi_energy_ev`.
    pub fermi_ev: Option<f64>,
    /// Two-component (spin-orbit) bands.
    pub spinor: Option<bool>,
    /// One per last-axis slice of `projections`.
    pub channels: Vec<DosChannel>,
    pub title: Option<String>,
    pub units: Option<Map<String, Value>>,
    /// Vendor keys (`x_*`).
    pub extra: Map<String, Value>,
}

/// The fat-band weights of a `bands` section.
#[derive(Debug, Clone, PartialEq)]
pub struct BandProjections {
    /// `[n_spin, n_k, n_bands, n_channels]`, or `[n_k, n_bands, n_channels]`
    /// when `spin_axis` is false.
    pub data: Vec<f64>,
    pub n_channels: usize,
    pub spin_axis: bool,
}

/// A `bands` section (spec § 5.6).
#[derive(Debug, Clone, PartialEq)]
pub struct Bands {
    pub id: String,
    pub kpath: Kpath,
    /// `[n_spin, n_k, n_bands]` in eV, absolute scale.
    pub eigenvalues: Vec<f64>,
    /// `[n_spin, n_k, n_bands]`; `n_spin` is 1 or 2.
    pub shape: [usize; 3],
    pub projections: Option<BandProjections>,
    /// `[n_spin, n_k, n_bands]` electrons per state.
    pub occupations: Option<Vec<f64>>,
    pub label: Option<String>,
}

impl Bands {
    pub fn new(id: &str, kpath: Kpath, eigenvalues: Vec<f64>, shape: [usize; 3]) -> Self {
        Bands { id: id.to_owned(), kpath, eigenvalues, shape, projections: None, occupations: None, label: None }
    }

    /// The `kpath` member payload.
    pub fn kpath_payload(&self) -> Result<Value, WriteError> {
        let what = format!("bands {:?}", self.id);
        let k = &self.kpath;
        let [n_spin, n_k, n_bands] = self.shape;
        let mut m = Map::new();
        if !k.kpoints.is_empty() {
            m.insert("kpoints".into(), json!(k.kpoints));
        }
        if !k.labels.is_empty() {
            m.insert("labels".into(), json!(k.labels));
        }
        if !k.segments.is_empty() {
            m.insert("segments".into(), Value::Array(k.segments.iter().map(KpathSegment::payload).collect()));
        }
        if let Some(f) = k.fermi_ev {
            m.insert("fermi".into(), Value::from(f));
            m.insert("fermi_energy_ev".into(), Value::from(f));
        }
        put(
            &mut m,
            vec![
                ("spinor", k.spinor.map(Value::from)),
                ("title", k.title.clone().map(Value::from)),
                ("units", k.units.clone().map(Value::Object)),
            ],
        );
        m.insert("n_spin".into(), Value::from(n_spin));
        m.insert("n_kpoints".into(), Value::from(n_k));
        m.insert("n_bands".into(), Value::from(n_bands));
        if !k.channels.is_empty() {
            let channels = k.channels.iter().map(|c| c.payload(&what)).collect::<Result<Vec<_>, _>>()?;
            m.insert("channels".into(), Value::Array(channels));
        }
        merge_extra(&what, &mut m, &k.extra)?;
        Ok(Value::Object(m))
    }

    fn check(&self) -> Result<(), WriteError> {
        let what = format!("bands {:?}", self.id);
        let [n_spin, n_k, n_bands] = self.shape;
        if !(n_spin == 1 || n_spin == 2) || n_k == 0 || n_bands == 0 {
            return invalid(format!("{what}: shape {:?} needs n_spin 1 or 2 and n_k, n_bands at least 1", self.shape));
        }
        if self.kpath.spinor == Some(true) && n_spin != 1 {
            return invalid(format!("{what}: spinor bands have n_spin 1 (spec § 5.6)"));
        }
        let n = n_spin * n_k * n_bands;
        if self.eigenvalues.len() != n {
            return invalid(format!("{what}: {} eigenvalues do not fill [{n_spin}, {n_k}, {n_bands}]", self.eigenvalues.len()));
        }
        finite(&format!("{what} eigenvalues"), &self.eigenvalues)?;
        let k = &self.kpath;
        if !k.kpoints.is_empty() && k.kpoints.len() != n_k {
            return invalid(format!("{what}: {} kpoints for {n_k} k-points", k.kpoints.len()));
        }
        if !k.labels.is_empty() && k.labels.len() != n_k {
            return invalid(format!("{what}: {} labels for {n_k} k-points", k.labels.len()));
        }
        finite(&format!("{what} kpoints"), &k.kpoints.concat())?;
        if let Some(f) = k.fermi_ev {
            finite(&format!("{what} fermi"), &[f])?;
        }
        if let Some(seg) = k.segments.iter().find(|s| s.n_points == Some(0)) {
            return invalid(format!("{what}: segment {seg:?} has no points"));
        }
        if let Some(units) = &k.units {
            for (key, v) in units {
                unit_name(&format!("{what}: kpath units.{key}"), v.as_str().unwrap_or_default())?;
            }
        }
        match &self.projections {
            Some(p) => {
                if p.n_channels == 0 {
                    return invalid(format!("{what}: projections need at least one channel"));
                }
                let expected = if p.spin_axis { n } else { n_k * n_bands } * p.n_channels;
                if p.data.len() != expected {
                    return invalid(format!("{what}: {} projection values for {expected}", p.data.len()));
                }
                if !k.channels.is_empty() && k.channels.len() != p.n_channels {
                    return invalid(format!(
                        "{what}: {} kpath channels for {} projection channels",
                        k.channels.len(),
                        p.n_channels
                    ));
                }
                finite(&format!("{what} projections"), &p.data)?;
            }
            None if !k.channels.is_empty() => return invalid(format!("{what}: kpath channels need projections")),
            None => {}
        }
        if let Some(o) = &self.occupations {
            if o.len() != n {
                return invalid(format!("{what}: occupations must be [{n_spin}, {n_k}, {n_bands}], got {} values", o.len()));
            }
            finite(&format!("{what} occupations"), o)?;
        }
        Ok(())
    }
}

// -- densities of states -----------------------------------------------------

fn spin_form(what: &str, n_spin: Option<usize>) -> Result<usize, WriteError> {
    match n_spin {
        None => Ok(1),
        Some(n @ (1 | 2)) => Ok(n),
        Some(n) => invalid(format!("{what}: n_spin must be 1 or 2, got {n}")),
    }
}

/// Broadening and Fermi level fields shared by the DOS kinds.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DosInfo {
    /// Broadening width in eV (≥ 0).
    pub smearing: Option<f64>,
    pub smearing_type: Option<String>,
    /// Absolute Fermi level the energies are referenced to.
    pub fermi_energy_ev: Option<f64>,
    /// Electrons (≥ 0); a section field of `dos.total` and `dos.projected`.
    pub n_electrons: Option<f64>,
}

impl DosInfo {
    fn fields(&self, what: &str) -> Result<Vec<(&'static str, Option<Value>)>, WriteError> {
        for (name, v) in [("smearing", self.smearing), ("n_electrons", self.n_electrons)] {
            if v.is_some_and(|x| !(x.is_finite() && x >= 0.0)) {
                return invalid(format!("{what}: {name} must be finite and non-negative"));
            }
        }
        if let Some(f) = self.fermi_energy_ev {
            finite(&format!("{what} fermi_energy_ev"), &[f])?;
        }
        Ok(vec![
            ("smearing", self.smearing.map(Value::from)),
            ("smearing_type", self.smearing_type.clone().map(Value::from)),
            ("fermi_energy_ev", self.fermi_energy_ev.map(Value::from)),
            ("n_electrons", self.n_electrons.map(Value::from)),
        ])
    }
}

/// A `dos.total` section (spec § 5.6).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DosTotal {
    pub id: String,
    /// `[n_points]` in eV, referenced to E_F at 0.
    pub energies: Vec<f64>,
    /// `[n_spin, n_points]` with `n_spin` set, `[n_points]` without, states/eV.
    pub dos: Vec<f64>,
    /// `Some` writes the spin axis; `None` the lower-rank form (n_spin 1).
    pub n_spin: Option<usize>,
    pub info: DosInfo,
    pub label: Option<String>,
}

impl DosTotal {
    pub fn new(id: &str, energies: Vec<f64>, dos: Vec<f64>) -> Self {
        DosTotal { id: id.to_owned(), energies, dos, ..DosTotal::default() }
    }
}

/// A `dos.projected` section (spec § 5.6).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DosProjected {
    pub id: String,
    pub energies: Vec<f64>,
    /// `[n_spin, n_channels, n_points]` with `n_spin` set, else
    /// `[n_channels, n_points]`; row `c` belongs to `channels[c]`.
    pub projections: Vec<f64>,
    pub channels: Vec<DosChannel>,
    pub n_spin: Option<usize>,
    pub info: DosInfo,
    pub label: Option<String>,
}

impl DosProjected {
    pub fn new(id: &str, energies: Vec<f64>, projections: Vec<f64>, channels: Vec<DosChannel>) -> Self {
        DosProjected { id: id.to_owned(), energies, projections, channels, ..DosProjected::default() }
    }
}

/// An atom pair of a COOP or COHP analysis.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DosPair {
    pub i: usize,
    pub j: usize,
    pub label: Option<String>,
    pub distance_ang: Option<f64>,
}

/// The bonding-analysis kind of a [`DosBonding`] section.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondingAnalysis {
    /// `dos.coop`, crystal orbital overlap populations.
    #[default]
    Coop,
    /// `dos.cohp`, crystal orbital Hamilton populations.
    Cohp,
}

/// A `dos.coop` or `dos.cohp` section (spec § 5.6).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DosBonding {
    pub id: String,
    pub analysis: BondingAnalysis,
    pub energies: Vec<f64>,
    /// `[n_spin, n_pairs, n_points]` with `n_spin` set, else `[n_pairs, n_points]`.
    pub projections: Vec<f64>,
    /// `[n_spin, n_pairs]` with `n_spin` set, else `[n_pairs]`.
    pub integrated: Vec<f64>,
    pub pairs: Vec<DosPair>,
    pub n_spin: Option<usize>,
    /// `smearing`, `smearing_type` and `fermi_energy_ev` go to the `meta`
    /// payload (`n_electrons` is not a field of it and is refused).
    pub info: DosInfo,
    /// Further `meta` keys (typed-open payload).
    pub extra: Map<String, Value>,
    pub label: Option<String>,
}

impl DosBonding {
    pub fn new(id: &str, analysis: BondingAnalysis, energies: Vec<f64>, pairs: Vec<DosPair>) -> Self {
        DosBonding { id: id.to_owned(), analysis, energies, pairs, ..DosBonding::default() }
    }

    /// `dos.coop` or `dos.cohp`.
    pub fn kind(&self) -> &'static str {
        match self.analysis {
            BondingAnalysis::Coop => "dos.coop",
            BondingAnalysis::Cohp => "dos.cohp",
        }
    }
}

fn energy_grid(what: &str, energies: &[f64]) -> Result<usize, WriteError> {
    if energies.is_empty() {
        return invalid(format!("{what}: the energy grid is empty"));
    }
    finite(&format!("{what} energies"), energies)?;
    Ok(energies.len())
}

// -- spectra -----------------------------------------------------------------

/// The one-dimensional spectrum kinds (`$defs/SpectrumPayload`, spec § 5.4).
pub const SPECTRUM_1D_KINDS: &[&str] =
    &["spectra.ir", "spectra.raman", "spectra.uvvis", "spectra.ecd", "spectra.vcd", "spectra.generic"];

/// An annotated peak of a spectrum.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Peak {
    pub x: f64,
    pub y: Option<f64>,
    /// Full width at half maximum (≥ 0), in the x unit.
    pub width: Option<f64>,
    pub label: Option<String>,
    /// 0-based index into the modes of the archive's `vibrations` section.
    pub mode_index: Option<u32>,
    /// The electronic state the transition reaches (an `index` of the
    /// `states_ref` section when the spectrum has one).
    pub state_index: Option<u32>,
}

/// The line shape of a spectrum.
#[derive(Debug, Clone, PartialEq)]
pub struct Broadening {
    /// `none`, `gaussian` or `lorentzian`.
    pub kind: String,
    /// Full width at half maximum (> 0).
    pub fwhm: Option<f64>,
    /// Unit of `fwhm` (spec § 4.3), else the x unit.
    pub units: Option<String>,
}

/// A one-dimensional `spectra.*` section (spec § 5.4).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Spectrum {
    /// The section id; empty takes the kind with `_` for `.` (`spectra_ir`).
    pub id: String,
    /// One of [`SPECTRUM_1D_KINDS`].
    pub kind: String,
    /// The x values, in the kind's § 4.3 unit unless `x_units` is set.
    pub frequencies: Vec<f64>,
    pub intensities: Vec<f64>,
    /// Also write `x`, the quantity-neutral alias of `frequencies`.
    pub write_x_alias: bool,
    pub x_units: Option<String>,
    pub y_units: Option<String>,
    pub peaks: Vec<Peak>,
    pub broadening: Option<Broadening>,
    /// Temperature in K (> 0).
    pub temperature_k: Option<f64>,
    /// `spectra.raman`: excitation wavelength in nm (> 0).
    pub raman_laser_nm: Option<f64>,
    /// `spectra.raman`: one depolarization ratio per line.
    pub raman_depolarization_ratios: Option<Vec<f64>>,
    /// `spectra.uvvis`: `energy_ev` or `wavelength_nm`.
    pub uvvis_x_axis: Option<String>,
    /// An `excited_states` section added earlier, whose state indices the
    /// peak assignments name.
    pub states_ref: Option<String>,
    /// Further payload keys (typed-open payload).
    pub extra: Map<String, Value>,
    pub label: Option<String>,
}

impl Spectrum {
    /// A spectrum of `kind` with its x and y values.
    pub fn new(kind: &str, frequencies: Vec<f64>, intensities: Vec<f64>) -> Self {
        Spectrum { kind: kind.to_owned(), frequencies, intensities, ..Spectrum::default() }
    }

    /// The section id this spectrum is written under.
    pub fn section_id(&self) -> String {
        if self.id.is_empty() {
            self.kind.replace('.', "_")
        } else {
            self.id.clone()
        }
    }

    /// The `spectrum` member payload.
    pub fn payload(&self) -> Result<Value, WriteError> {
        let what = format!("{} {:?}", self.kind, self.section_id());
        let mut m = Map::new();
        m.insert("frequencies".into(), json!(self.frequencies));
        m.insert("intensities".into(), json!(self.intensities));
        if self.write_x_alias {
            m.insert("x".into(), json!(self.frequencies));
        }
        put(
            &mut m,
            vec![
                ("x_units", self.x_units.clone().map(Value::from)),
                ("y_units", self.y_units.clone().map(Value::from)),
                ("temperature_k", self.temperature_k.map(Value::from)),
                ("states_ref", self.states_ref.clone().map(Value::from)),
            ],
        );
        if !self.peaks.is_empty() {
            let peaks: Vec<Value> = self
                .peaks
                .iter()
                .map(|p| {
                    let mut o = Map::new();
                    o.insert("x".into(), Value::from(p.x));
                    put(
                        &mut o,
                        vec![
                            ("y", p.y.map(Value::from)),
                            ("width", p.width.map(Value::from)),
                            ("label", p.label.clone().map(Value::from)),
                        ],
                    );
                    if p.mode_index.is_some() || p.state_index.is_some() {
                        let mut a = Map::new();
                        put(
                            &mut a,
                            vec![("mode_index", p.mode_index.map(Value::from)), ("state_index", p.state_index.map(Value::from))],
                        );
                        o.insert("assignment".into(), Value::Object(a));
                    }
                    Value::Object(o)
                })
                .collect();
            m.insert("peaks".into(), Value::Array(peaks));
        }
        if let Some(b) = &self.broadening {
            let mut o = Map::new();
            o.insert("type".into(), Value::from(b.kind.as_str()));
            put(&mut o, vec![("fwhm", b.fwhm.map(Value::from)), ("units", b.units.clone().map(Value::from))]);
            m.insert("broadening".into(), Value::Object(o));
        }
        if self.raman_laser_nm.is_some() || self.raman_depolarization_ratios.is_some() {
            let mut o = Map::new();
            put(
                &mut o,
                vec![
                    ("laser_nm", self.raman_laser_nm.map(Value::from)),
                    ("depolarization_ratios", self.raman_depolarization_ratios.as_ref().map(|r| json!(r))),
                ],
            );
            m.insert("raman".into(), Value::Object(o));
        }
        if let Some(axis) = &self.uvvis_x_axis {
            m.insert("uvvis".into(), json!({ "x_axis": axis }));
        }
        merge_extra(&what, &mut m, &self.extra)?;
        Ok(Value::Object(m))
    }

    fn check(&self) -> Result<(), WriteError> {
        let what = format!("{} {:?}", self.kind, self.section_id());
        if !SPECTRUM_1D_KINDS.contains(&self.kind.as_str()) {
            return invalid(format!("{:?} is not a one-dimensional spectrum kind (spec § 5.4)", self.kind));
        }
        if self.intensities.len() != self.frequencies.len() {
            return invalid(format!(
                "{what}: intensities has {} entries but frequencies has {}; they must have the same length (spec § 5.4)",
                self.intensities.len(),
                self.frequencies.len()
            ));
        }
        finite(&format!("{what} frequencies"), &self.frequencies)?;
        finite(&format!("{what} intensities"), &self.intensities)?;
        for (key, u) in [("x_units", &self.x_units), ("y_units", &self.y_units)] {
            if let Some(u) = u {
                unit_name(&format!("{what}: {key}"), u)?;
            }
        }
        for (k, p) in self.peaks.iter().enumerate() {
            finite(&format!("{what} peaks[{k}]"), &[p.x, p.y.unwrap_or(0.0), p.width.unwrap_or(0.0)])?;
            if p.width.is_some_and(|w| w < 0.0) {
                return invalid(format!("{what}: peaks[{k}].width must be non-negative"));
            }
        }
        if let Some(b) = &self.broadening {
            if !["none", "gaussian", "lorentzian"].contains(&b.kind.as_str()) {
                return invalid(format!("{what}: broadening type {:?} is not none, gaussian or lorentzian", b.kind));
            }
            if b.fwhm.is_some_and(|w| !(w.is_finite() && w > 0.0)) {
                return invalid(format!("{what}: broadening fwhm must be positive"));
            }
            if let Some(u) = &b.units {
                unit_name(&format!("{what}: broadening units"), u)?;
            }
        }
        if self.temperature_k.is_some_and(|t| !(t.is_finite() && t > 0.0)) {
            return invalid(format!("{what}: temperature_k must be positive"));
        }
        let raman = self.raman_laser_nm.is_some() || self.raman_depolarization_ratios.is_some();
        if raman && self.kind != "spectra.raman" {
            return invalid(format!("{what}: the raman block belongs to spectra.raman (spec § 5.4)"));
        }
        if self.raman_laser_nm.is_some_and(|l| !(l.is_finite() && l > 0.0)) {
            return invalid(format!("{what}: raman laser_nm must be positive"));
        }
        if let Some(r) = &self.raman_depolarization_ratios {
            if r.len() != self.frequencies.len() {
                return invalid(format!("{what}: one depolarization ratio per line ({})", self.frequencies.len()));
            }
            finite(&format!("{what} depolarization_ratios"), r)?;
        }
        if let Some(axis) = &self.uvvis_x_axis {
            if self.kind != "spectra.uvvis" {
                return invalid(format!("{what}: the uvvis block belongs to spectra.uvvis (spec § 5.4)"));
            }
            if axis != "energy_ev" && axis != "wavelength_nm" {
                return invalid(format!("{what}: uvvis x_axis {axis:?} is not energy_ev or wavelength_nm"));
            }
        }
        Ok(())
    }
}

// -- writing -----------------------------------------------------------------

impl<W: Write> Writer<W> {
    /// Add a `trajectory` section: `metadata`, `coords` and the per-frame
    /// members that are set (`lattice`, `velocities`, `gradient`, `times`,
    /// `energies`, `cell_volumes`, `temperatures`).
    pub fn add_trajectory(&mut self, t: &Trajectory) -> Result<String, WriteError> {
        let what = format!("trajectory {:?}", t.id);
        let (n_frames, n_atoms) = t.check()?;
        let meta = t.metadata_payload()?;
        schema(&what, "TrajectoryMetadataPayload", &meta)?;
        let fields = match &t.label {
            Some(l) => Map::from_iter([("label".to_owned(), Value::from(l.as_str()))]),
            None => Map::new(),
        };
        self.add_section(&t.id, "trajectory", fields)?;
        self.add_json(&t.id, "metadata", &meta)?;
        let frame_shape = [n_frames, n_atoms, 3];
        match t.precision {
            Precision::Float64 => {
                self.add_tensor(&t.id, "coords", &t.coords, &frame_shape)?;
            }
            Precision::Float32 => {
                let mut sink = self.begin_tensor(&t.id, "coords", DType::Float32, &frame_shape)?;
                for chunk in t.coords.chunks(8192) {
                    let c: Vec<f32> = chunk.iter().map(|&x| x as f32).collect();
                    sink.write(&c)?;
                }
                sink.finish()?;
            }
        }
        match &t.cell {
            Some(TrajectoryCell::Fixed(c)) => {
                self.add_tensor(&t.id, "lattice", &c.concat(), &[3, 3])?;
            }
            Some(TrajectoryCell::PerFrame(list)) => {
                let flat: Vec<f64> = list.iter().flat_map(|c| c.concat()).collect();
                self.add_tensor(&t.id, "lattice", &flat, &[n_frames, 3, 3])?;
            }
            None => {}
        }
        for (role, v) in [("velocities", &t.velocities), ("gradient", &t.gradient)] {
            if let Some(v) = v {
                self.add_tensor(&t.id, role, v, &frame_shape)?;
            }
        }
        for (role, v) in
            [("times", &t.times), ("energies", &t.energies), ("cell_volumes", &t.cell_volumes), ("temperatures", &t.temperatures)]
        {
            if let Some(v) = v {
                self.add_tensor(&t.id, role, v, &[n_frames])?;
            }
        }
        if let Some(info) = self.info.get_mut(&t.id) {
            info.atoms = Some(n_atoms);
        }
        Ok(t.id.clone())
    }

    /// Add a `bands` section: `kpath`, `eigenvalues` and the optional
    /// `projections` and `occupations`.
    pub fn add_bands(&mut self, b: &Bands) -> Result<String, WriteError> {
        let what = format!("bands {:?}", b.id);
        b.check()?;
        let kpath = b.kpath_payload()?;
        schema(&what, "KpathPayload", &kpath)?;
        let fields = match &b.label {
            Some(l) => Map::from_iter([("label".to_owned(), Value::from(l.as_str()))]),
            None => Map::new(),
        };
        let [n_spin, n_k, n_bands] = b.shape;
        self.add_section(&b.id, "bands", fields)?;
        self.add_json(&b.id, "kpath", &kpath)?;
        self.add_tensor(&b.id, "eigenvalues", &b.eigenvalues, &b.shape)?;
        if let Some(p) = &b.projections {
            let shape: Vec<usize> =
                if p.spin_axis { vec![n_spin, n_k, n_bands, p.n_channels] } else { vec![n_k, n_bands, p.n_channels] };
            self.add_tensor(&b.id, "projections", &p.data, &shape)?;
        }
        if let Some(o) = &b.occupations {
            self.add_tensor(&b.id, "occupations", o, &b.shape)?;
        }
        Ok(b.id.clone())
    }

    fn dos_fields(what: &str, n_spin: usize, info: &DosInfo, label: &Option<String>) -> Result<Map<String, Value>, WriteError> {
        let mut fields = Map::new();
        fields.insert("n_spin".into(), Value::from(n_spin));
        put(&mut fields, info.fields(what)?);
        put(&mut fields, vec![("label", label.clone().map(Value::from))]);
        Ok(fields)
    }

    /// Add a `dos.total` section: `energies` and `dos`, with `n_spin` and the
    /// broadening fields on the section.
    pub fn add_dos_total(&mut self, d: &DosTotal) -> Result<String, WriteError> {
        let what = format!("dos.total {:?}", d.id);
        let n_points = energy_grid(&what, &d.energies)?;
        let n_spin = spin_form(&what, d.n_spin)?;
        if d.dos.len() != n_spin * n_points {
            return invalid(format!("{what}: {} dos values for n_spin {n_spin} and {n_points} points", d.dos.len()));
        }
        finite(&format!("{what} dos"), &d.dos)?;
        let fields = Self::dos_fields(&what, n_spin, &d.info, &d.label)?;
        self.add_section(&d.id, "dos.total", fields)?;
        self.add_tensor(&d.id, "energies", &d.energies, &[n_points])?;
        let shape: Vec<usize> = if d.n_spin.is_some() { vec![n_spin, n_points] } else { vec![n_points] };
        self.add_tensor(&d.id, "dos", &d.dos, &shape)?;
        Ok(d.id.clone())
    }

    /// Add a `dos.projected` section: `energies` and `projections`, with
    /// `n_spin`, `channels` and the broadening fields on the section.
    pub fn add_dos_projected(&mut self, d: &DosProjected) -> Result<String, WriteError> {
        let what = format!("dos.projected {:?}", d.id);
        let n_points = energy_grid(&what, &d.energies)?;
        let n_spin = spin_form(&what, d.n_spin)?;
        let n_channels = d.channels.len();
        if n_channels == 0 {
            return invalid(format!("{what}: needs at least one channel"));
        }
        if d.projections.len() != n_spin * n_channels * n_points {
            return invalid(format!(
                "{what}: {} projection values for n_spin {n_spin}, {n_channels} channels and {n_points} points",
                d.projections.len()
            ));
        }
        finite(&format!("{what} projections"), &d.projections)?;
        let mut fields = Self::dos_fields(&what, n_spin, &d.info, &d.label)?;
        let channels = d.channels.iter().map(|c| c.payload(&what)).collect::<Result<Vec<_>, _>>()?;
        fields.insert("channels".into(), Value::Array(channels));
        self.add_section(&d.id, "dos.projected", fields)?;
        self.add_tensor(&d.id, "energies", &d.energies, &[n_points])?;
        let shape: Vec<usize> = if d.n_spin.is_some() { vec![n_spin, n_channels, n_points] } else { vec![n_channels, n_points] };
        self.add_tensor(&d.id, "projections", &d.projections, &shape)?;
        Ok(d.id.clone())
    }

    /// Add a `dos.coop` or `dos.cohp` section: `energies`, `projections`,
    /// `integrated` and the `meta` payload (pairs, `n_spin`, Fermi level).
    pub fn add_dos_bonding(&mut self, d: &DosBonding) -> Result<String, WriteError> {
        let kind = d.kind();
        let what = format!("{kind} {:?}", d.id);
        let n_points = energy_grid(&what, &d.energies)?;
        let n_spin = spin_form(&what, d.n_spin)?;
        let n_pairs = d.pairs.len();
        if n_pairs == 0 {
            return invalid(format!("{what}: needs at least one atom pair"));
        }
        if d.projections.len() != n_spin * n_pairs * n_points {
            return invalid(format!(
                "{what}: {} projection values for n_spin {n_spin}, {n_pairs} pairs and {n_points} points",
                d.projections.len()
            ));
        }
        if d.integrated.len() != n_spin * n_pairs {
            return invalid(format!("{what}: {} integrated values for n_spin {n_spin} and {n_pairs} pairs", d.integrated.len()));
        }
        finite(&format!("{what} projections"), &d.projections)?;
        finite(&format!("{what} integrated"), &d.integrated)?;
        if d.info.n_electrons.is_some() {
            return invalid(format!("{what}: n_electrons is not a field of the meta payload"));
        }
        let mut meta = Map::new();
        meta.insert("method".into(), Value::from(&kind[4..]));
        meta.insert("n_spin".into(), Value::from(n_spin));
        put(&mut meta, d.info.fields(&what)?);
        let mut pairs = Vec::with_capacity(n_pairs);
        for p in &d.pairs {
            if p.distance_ang.is_some_and(|x| !(x.is_finite() && x >= 0.0)) {
                return invalid(format!("{what}: pair ({}, {}) distance must be non-negative", p.i, p.j));
            }
            if p.label.is_none() && p.distance_ang.is_none() {
                pairs.push(json!([p.i, p.j]));
            } else {
                let mut o = Map::from_iter([("i".to_owned(), Value::from(p.i)), ("j".to_owned(), Value::from(p.j))]);
                put(&mut o, vec![("label", p.label.clone().map(Value::from)), ("distance_ang", p.distance_ang.map(Value::from))]);
                pairs.push(Value::Object(o));
            }
        }
        meta.insert("pairs".into(), Value::Array(pairs));
        merge_extra(&what, &mut meta, &d.extra)?;
        let meta = Value::Object(meta);
        schema(&what, "DosMetaPayload", &meta)?;
        let fields = match &d.label {
            Some(l) => Map::from_iter([("label".to_owned(), Value::from(l.as_str()))]),
            None => Map::new(),
        };
        self.add_section(&d.id, kind, fields)?;
        self.add_tensor(&d.id, "energies", &d.energies, &[n_points])?;
        let (proj, integ): (Vec<usize>, Vec<usize>) = if d.n_spin.is_some() {
            (vec![n_spin, n_pairs, n_points], vec![n_spin, n_pairs])
        } else {
            (vec![n_pairs, n_points], vec![n_pairs])
        };
        self.add_tensor(&d.id, "projections", &d.projections, &proj)?;
        self.add_tensor(&d.id, "integrated", &d.integrated, &integ)?;
        self.add_json(&d.id, "meta", &meta)?;
        Ok(d.id.clone())
    }

    /// Add a one-dimensional `spectra.*` section with its `spectrum` member
    /// (path `spectra/<id>.json`, as the Python writer names it).
    pub fn add_spectrum(&mut self, s: &Spectrum) -> Result<String, WriteError> {
        s.check()?;
        let id = s.section_id();
        let what = format!("{} {id:?}", s.kind);
        if let Some(r) = &s.states_ref {
            let states = self.require_kind(&format!("{what}: states_ref"), r, &["excited_states"])?.states.clone();
            for (k, p) in s.peaks.iter().enumerate() {
                if let (Some(i), Some(states)) = (p.state_index, &states) {
                    if !states.contains(&i) {
                        return invalid(format!("{what}: peaks[{k}] state_index {i} is not a state of {r:?} (spec § 5.4)"));
                    }
                }
            }
        }
        let payload = s.payload()?;
        schema(&what, "SpectrumPayload", &payload)?;
        let fields = match &s.label {
            Some(l) => Map::from_iter([("label".to_owned(), Value::from(l.as_str()))]),
            None => Map::new(),
        };
        self.add_section(&id, &s.kind, fields)?;
        self.add_json_at(&id, "spectrum", &payload, &format!("spectra/{id}.json"))?;
        Ok(id)
    }

    /// Add a `spectra.*` section of any kind (`spectra.nmr` and
    /// `spectra.epr` included) from its whole payload, checked against the
    /// kind's payload definition.
    pub fn add_spectrum_payload(
        &mut self,
        id: &str,
        kind: &str,
        payload: &Value,
        label: Option<&str>,
    ) -> Result<String, WriteError> {
        let what = format!("{kind} {id:?}");
        let Some(def) = crate::schema_check::payload_definition(kind, "spectrum").filter(|_| kind.starts_with("spectra.")) else {
            return invalid(format!("{kind:?} is not a spectra kind (spec § 5.4)"));
        };
        schema(&what, def, payload)?;
        if SPECTRUM_1D_KINDS.contains(&kind) {
            let len = |k: &str| payload.get(k).and_then(Value::as_array).map(Vec::len);
            for key in ["intensities", "x"] {
                if len(key).is_some_and(|n| Some(n) != len("frequencies")) {
                    return invalid(format!("{what}: {key} and frequencies must have the same length (spec § 5.4)"));
                }
            }
            if let Some(r) = payload.get("states_ref").and_then(Value::as_str) {
                self.require_kind(&format!("{what}: states_ref"), r, &["excited_states"])?;
            }
        }
        let fields = match label {
            Some(l) => Map::from_iter([("label".to_owned(), Value::from(l))]),
            None => Map::new(),
        };
        self.add_section(id, kind, fields)?;
        self.add_json_at(id, "spectrum", payload, &format!("spectra/{id}.json"))?;
        Ok(id.to_owned())
    }
}
