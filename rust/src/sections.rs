//! Section builders: typed records for the common kinds, written through a
//! [`Writer`] with their payload layout (spec § 5).
//!
//! * [`Structure`]: `structure` (§ 5.1), positions in Å.
//! * [`WavefunctionGto`]: `wavefunction.gto` (§ 5.3, Appendix A), restricted,
//!   restricted open-shell or unrestricted real coefficients, with the
//!   optional `overlap` and density matrices.
//! * [`Volume`]: every `volume.*` kind (§ 5.2), grids in bohr, streamed in
//!   slabs if wanted.
//! * [`ExcitedStates`]: `excited_states` (§ 5.12).
//!
//! The builders check what can be checked before anything is written: shapes
//! and lengths, references to sections added earlier, atom indices, orbital
//! rows. Checks that need the whole archive (references to sections added
//! later, the overlap diagonal of Appendix A.1) are left to the validator.

use std::io::Write;

use serde_json::{json, Map, Value};

use crate::reader::DType;
use crate::writer::{invalid, Element, TensorSink, WriteError, Writer};

pub use crate::builders::{
    BandProjections, Bands, BondingAnalysis, Broadening, DosBonding, DosChannel, DosInfo, DosPair, DosProjected, DosTotal, Kpath,
    KpathSegment, Peak, Precision, Spectrum, Trajectory, TrajectoryCell, SPECTRUM_1D_KINDS,
};

/// Element symbols by atomic number (index 0 is a dummy atom).
pub const ELEMENTS: [&str; 119] = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
    "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

fn finite(what: &str, values: &[f64]) -> Result<(), WriteError> {
    match values.iter().position(|x| !x.is_finite()) {
        Some(i) => invalid(format!("{what}[{i}] is not finite")),
        None => Ok(()),
    }
}

fn fields_from(pairs: Vec<(&str, Option<Value>)>) -> Map<String, Value> {
    pairs.into_iter().filter_map(|(k, v)| v.map(|v| (k.to_owned(), v))).collect()
}

// -- structure ---------------------------------------------------------------

/// One atom of a structure.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub symbol: String,
    pub atomic_number: u32,
    /// Cartesian position in Å.
    pub position: [f64; 3],
}

impl Atom {
    /// An atom of the element `symbol`. An unknown symbol gets atomic number 0.
    pub fn new(symbol: &str, position: [f64; 3]) -> Self {
        let z = ELEMENTS.iter().position(|s| *s == symbol).unwrap_or(0);
        Atom { symbol: symbol.to_owned(), atomic_number: z as u32, position }
    }

    /// An atom of atomic number `z`.
    pub fn from_z(z: u32, position: [f64; 3]) -> Self {
        let symbol = ELEMENTS.get(z as usize).copied().unwrap_or("X");
        Atom { symbol: symbol.to_owned(), atomic_number: z, position }
    }
}

/// A bond of the `bonds` member of a structure (atom indices from 0).
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub order: Option<f64>,
}

/// A `structure` section (spec § 5.1).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Structure {
    pub id: String,
    pub atoms: Vec<Atom>,
    /// Periodic axes. A periodic axis needs `lattice_vectors`.
    pub pbc: [bool; 3],
    /// Lattice vectors as rows, in Å. Row `i` pairs with `pbc[i]`.
    pub lattice_vectors: Option<[[f64; 3]; 3]>,
    pub bonds: Vec<Bond>,
    pub charge: Option<i64>,
    pub multiplicity: Option<u32>,
    pub name: Option<String>,
    pub label: Option<String>,
}

impl Structure {
    /// A molecule (no periodic axis).
    pub fn molecule(id: &str, atoms: Vec<Atom>) -> Self {
        Structure { id: id.to_owned(), atoms, ..Structure::default() }
    }

    /// A periodic structure.
    pub fn periodic(id: &str, atoms: Vec<Atom>, lattice_vectors: [[f64; 3]; 3], pbc: [bool; 3]) -> Self {
        Structure { id: id.to_owned(), atoms, pbc, lattice_vectors: Some(lattice_vectors), ..Structure::default() }
    }

    fn check(&self) -> Result<(), WriteError> {
        if self.pbc.iter().any(|&p| p) && self.lattice_vectors.is_none() {
            return invalid(format!("structure {:?}: a periodic axis needs lattice_vectors (spec § 5.1)", self.id));
        }
        for (k, a) in self.atoms.iter().enumerate() {
            finite(&format!("structure {:?}: atoms[{k}].position", self.id), &a.position)?;
        }
        if let Some(lv) = &self.lattice_vectors {
            finite(&format!("structure {:?}: lattice_vectors", self.id), &lv.concat())?;
        }
        if let Some(b) = self.bonds.iter().find(|b| b.i >= self.atoms.len() || b.j >= self.atoms.len() || b.i == b.j) {
            return invalid(format!("structure {:?}: bond ({}, {}) does not join two of its atoms", self.id, b.i, b.j));
        }
        if self.multiplicity == Some(0) {
            return invalid(format!("structure {:?}: multiplicity must be at least 1", self.id));
        }
        if let (Some(charge), Some(mult)) = (self.charge, self.multiplicity) {
            let electrons = self.atoms.iter().map(|a| i64::from(a.atomic_number)).sum::<i64>() - charge;
            let unpaired = i64::from(mult) - 1;
            if electrons < 0 || unpaired > electrons || (electrons - unpaired) % 2 != 0 {
                return invalid(format!(
                    "structure {:?}: charge {charge} and multiplicity {mult} are impossible for {electrons} electrons (spec § 5.1)",
                    self.id
                ));
            }
        }
        Ok(())
    }

    /// The `structure` member payload.
    pub fn payload(&self) -> Value {
        let atoms: Vec<Value> = self
            .atoms
            .iter()
            .map(|a| json!({"symbol": a.symbol, "position": a.position, "atomic_number": a.atomic_number}))
            .collect();
        json!({"atoms": atoms, "pbc": self.pbc, "lattice_vectors": self.lattice_vectors})
    }
}

// -- wavefunction.gto --------------------------------------------------------

/// A contracted shell of a GTO basis (Appendix A): contraction coefficients
/// apply to `N_i`-normalized primitives.
#[derive(Debug, Clone, PartialEq)]
pub struct Shell {
    /// Atom index into the structure, from 0.
    pub center: usize,
    pub l: u32,
    pub exponents: Vec<f64>,
    pub coefficients: Vec<f64>,
    /// Spherical (`true`) or Cartesian; `None` takes the basis default.
    pub pure: Option<bool>,
}

impl Shell {
    pub fn new(center: usize, l: u32, exponents: Vec<f64>, coefficients: Vec<f64>) -> Self {
        Shell { center, l, exponents, coefficients, pure: None }
    }
}

/// The `basis` payload of a wavefunction.
#[derive(Debug, Clone, PartialEq)]
pub struct Basis {
    /// The `structure` section the shell centres index.
    pub structure_ref: String,
    /// Spherical (`true`) or Cartesian shells by default.
    pub pure: bool,
    pub shells: Vec<Shell>,
    /// Electrons the wavefunction describes; the occupations must sum to it.
    pub n_electrons: Option<u64>,
}

impl Basis {
    pub fn new(structure_ref: &str, pure: bool, shells: Vec<Shell>) -> Self {
        Basis { structure_ref: structure_ref.to_owned(), pure, shells, n_electrons: None }
    }

    /// The number of AOs: 2l + 1 per spherical and (l + 1)(l + 2)/2 per
    /// Cartesian shell.
    pub fn n_ao(&self) -> usize {
        self.shells
            .iter()
            .map(|s| {
                let l = s.l as usize;
                if s.pure.unwrap_or(self.pure) {
                    2 * l + 1
                } else {
                    (l + 1) * (l + 2) / 2
                }
            })
            .sum()
    }

    /// The `basis` member payload.
    pub fn payload(&self) -> Value {
        let shells: Vec<Value> = self
            .shells
            .iter()
            .map(|s| {
                json!({"center": s.center, "l": s.l, "exponents": s.exponents,
                       "coefficients": s.coefficients, "pure": s.pure.unwrap_or(self.pure)})
            })
            .collect();
        let mut out = json!({"structure_ref": self.structure_ref, "pure": self.pure, "n_ao": self.n_ao(), "shells": shells});
        if let Some(n) = self.n_electrons {
            out["n_electrons"] = Value::from(n);
        }
        out
    }
}

/// One set of orbitals: `coefficients` is `[n_mo, n_ao]` in C order (rows
/// are MOs). Empty `energies` or `occupations` are written as zeros.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Orbitals {
    pub coefficients: Vec<f64>,
    /// Orbital energies in hartree.
    pub energies: Vec<f64>,
    pub occupations: Vec<f64>,
    pub symmetry_labels: Option<Vec<String>>,
}

impl Orbitals {
    pub fn new(coefficients: Vec<f64>, energies: Vec<f64>, occupations: Vec<f64>) -> Self {
        Orbitals { coefficients, energies, occupations, symmetry_labels: None }
    }

    fn rows(&self, n_ao: usize, what: &str) -> Result<usize, WriteError> {
        if n_ao == 0 || self.coefficients.len() % n_ao != 0 {
            return invalid(format!("{what}: {} coefficients are not rows of n_ao = {n_ao}", self.coefficients.len()));
        }
        finite(what, &self.coefficients)?;
        let n = self.coefficients.len() / n_ao;
        for (name, v) in [("energies", &self.energies), ("occupations", &self.occupations)] {
            if !v.is_empty() && v.len() != n {
                return invalid(format!("{what}: {} {name} for {n} orbitals", v.len()));
            }
            finite(&format!("{what} {name}"), v)?;
        }
        if self.symmetry_labels.as_ref().is_some_and(|s| s.len() != n) {
            return invalid(format!("{what}: symmetry_labels needs one label per orbital ({n})"));
        }
        Ok(n)
    }

    fn block(&self, n: usize) -> Map<String, Value> {
        let or_zeros = |v: &Vec<f64>| if v.is_empty() { vec![0.0; n] } else { v.clone() };
        let mut m = Map::new();
        m.insert("n_mo".into(), Value::from(n));
        m.insert("energies".into(), Value::from(or_zeros(&self.energies)));
        m.insert("occupations".into(), Value::from(or_zeros(&self.occupations)));
        if let Some(s) = &self.symmetry_labels {
            m.insert("symmetry_labels".into(), Value::from(s.clone()));
        }
        m
    }
}

/// The orbitals of a wavefunction and their spin treatment.
#[derive(Debug, Clone, PartialEq)]
pub enum SpinOrbitals {
    /// One set, doubly occupied up to the frontier (`spin: "restricted"`).
    Restricted(Orbitals),
    /// One set of a restricted open-shell wavefunction.
    RestrictedOpen(Orbitals),
    /// Separate α and β sets (`spin: "unrestricted"`).
    Unrestricted { alpha: Orbitals, beta: Orbitals },
}

/// The hole or particle orbitals of a natural-transition-orbital pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NtoRole {
    Hole,
    Particle,
}

/// Links a section of natural transition orbitals to its state (spec § 5.3).
#[derive(Debug, Clone, PartialEq)]
pub struct NtoLink {
    /// The `excited_states` section (it may be added later).
    pub state_ref: String,
    /// The state's `index` there.
    pub state_index: u32,
    pub role: NtoRole,
}

/// A `wavefunction.gto` section (spec § 5.3) with real coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct WavefunctionGto {
    pub id: String,
    pub basis: Basis,
    pub orbitals: SpinOrbitals,
    /// `canonical`, `natural`, `localized`, `nto`, … (spec § 5.3).
    pub orbital_kind: String,
    /// The AO overlap `[n_ao, n_ao]` of the AOs as written; producers
    /// SHOULD write it so that validators check the normalization.
    pub overlap: Option<Vec<f64>>,
    /// The total density matrix of a single set.
    pub density_matrix: Option<Vec<f64>>,
    /// The α and β density matrices of an unrestricted wavefunction.
    pub density_matrices_ab: Option<(Vec<f64>, Vec<f64>)>,
    /// Set for natural transition orbitals (`orbital_kind` `nto`).
    pub nto: Option<NtoLink>,
    pub label: Option<String>,
}

impl WavefunctionGto {
    pub fn new(id: &str, basis: Basis, orbitals: SpinOrbitals) -> Self {
        WavefunctionGto {
            id: id.to_owned(),
            basis,
            orbitals,
            orbital_kind: "canonical".into(),
            overlap: None,
            density_matrix: None,
            density_matrices_ab: None,
            nto: None,
            label: None,
        }
    }
}

// -- volumes -----------------------------------------------------------------

/// The scalar volume kinds, whose data is `[nx, ny, nz]` (spec § 5.2).
pub const SCALAR_VOLUME_KINDS: &[&str] = &[
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

/// The `axis` of a `volume.series` (one value per slab).
#[derive(Debug, Clone, PartialEq)]
pub struct SeriesAxis {
    /// `time`, `frame`, `state` or `custom`.
    pub kind: String,
    pub values: Vec<f64>,
    pub units: Option<String>,
    pub labels: Option<Vec<String>>,
    pub trajectory_ref: Option<String>,
}

/// A `volume.*` section (spec § 5.2). The grid is in bohr.
#[derive(Debug, Clone, PartialEq)]
pub struct Volume {
    pub id: String,
    pub kind: String,
    pub origin: [f64; 3],
    /// Step vectors of the three grid axes, rows, in bohr.
    pub voxel_vectors: [[f64; 3]; 3],
    /// Points along the three axes.
    pub shape: [usize; 3],
    pub periodic: Option<[bool; 3]>,
    pub includes_endpoint: Option<bool>,
    pub structure_ref: Option<String>,
    /// The grid's `units` declaration (spec § 4.3), keyed by field or role.
    pub units: Option<Map<String, Value>>,
    /// Required for `volume.series`, refused for the others.
    pub axis: Option<SeriesAxis>,
    pub label: Option<String>,
    /// Free text of the represented component (spin channel, orbital name).
    pub component: Option<String>,
    /// The two operand sections of a `volume.difference`.
    pub operands: Option<(String, String)>,
    /// `excited_states` section and state index of a `volume.transition_density`.
    pub excited_state: Option<(String, Option<u32>)>,
}

impl Volume {
    /// A volume of `kind` on a grid.
    pub fn new(id: &str, kind: &str, origin: [f64; 3], voxel_vectors: [[f64; 3]; 3], shape: [usize; 3]) -> Self {
        Volume {
            id: id.to_owned(),
            kind: kind.to_owned(),
            origin,
            voxel_vectors,
            shape,
            periodic: None,
            includes_endpoint: None,
            structure_ref: None,
            units: None,
            axis: None,
            label: None,
            component: None,
            operands: None,
            excited_state: None,
        }
    }

    /// The shape of the `data` tensor of this kind.
    pub fn data_shape(&self) -> Result<Vec<usize>, WriteError> {
        let [nx, ny, nz] = self.shape;
        let shape = match self.kind.as_str() {
            "volume.vector" => vec![nx, ny, nz, 3],
            "volume.complex" => vec![nx, ny, nz, 2],
            "volume.series" => match &self.axis {
                Some(a) => vec![a.values.len(), nx, ny, nz],
                None => return invalid(format!("volume.series {:?} needs an axis (spec § 5.2)", self.id)),
            },
            k if SCALAR_VOLUME_KINDS.contains(&k) => vec![nx, ny, nz],
            k => return invalid(format!("{k:?} is not a volume kind")),
        };
        Ok(shape)
    }

    /// The `grid` member payload.
    pub fn grid_payload(&self) -> Value {
        let mut g = json!({"origin": self.origin, "voxel_vectors": self.voxel_vectors, "shape": self.shape});
        if let Some(p) = self.periodic {
            g["periodic"] = json!(p);
        }
        if let Some(e) = self.includes_endpoint {
            g["includes_endpoint"] = Value::from(e);
        }
        if let Some(s) = &self.structure_ref {
            g["structure_ref"] = Value::from(s.as_str());
        }
        if let Some(u) = &self.units {
            g["units"] = Value::Object(u.clone());
        }
        g
    }

    fn axis_payload(&self) -> Option<Value> {
        let a = self.axis.as_ref()?;
        let values: Vec<Value> = if a.kind == "frame" {
            a.values.iter().map(|&v| Value::from(v as u64)).collect()
        } else {
            a.values.iter().map(|&v| Value::from(v)).collect()
        };
        let mut out = json!({"kind": a.kind, "values": values});
        if let Some(u) = &a.units {
            out["units"] = Value::from(u.as_str());
        }
        if let Some(l) = &a.labels {
            out["labels"] = json!(l);
        }
        if let Some(t) = &a.trajectory_ref {
            out["trajectory_ref"] = Value::from(t.as_str());
        }
        Some(out)
    }
}

// -- excited states ----------------------------------------------------------

/// The spin channel of a contribution of an unrestricted reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Channel {
    Alpha,
    Beta,
}

/// One orbital excitation of a state: `from_mo` to `to_mo` (orbital rows of
/// the reference wavefunction, from 0) with its weight.
#[derive(Debug, Clone, PartialEq)]
pub struct Contribution {
    pub from_mo: usize,
    pub to_mo: usize,
    pub spin: Option<Channel>,
    pub weight: f64,
}

/// One excited state (`$defs/ExcitedState`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ExcitedState {
    /// From 1, ascending and unique within the section.
    pub index: u32,
    /// Excitation energy in eV.
    pub excitation_energy_ev: f64,
    /// Total energy of the state in hartree.
    pub total_energy_hartree: Option<f64>,
    pub multiplicity: Option<u32>,
    pub symmetry: Option<String>,
    pub oscillator_strength: Option<f64>,
    /// Rotatory strength in 10^-40 esu² cm².
    pub rotatory_strength_cgs: Option<f64>,
    /// Transition electric dipole, atomic units.
    pub transition_dipole: Option<[f64; 3]>,
    /// Transition magnetic dipole, atomic units.
    pub transition_magnetic_dipole: Option<[f64; 3]>,
    pub contributions: Vec<Contribution>,
    /// A `wavefunction.gto` section of NTOs, added earlier.
    pub nto_ref: Option<String>,
    /// A `volume.transition_density` or `volume.difference` section, added earlier.
    pub density_ref: Option<String>,
}

impl ExcitedState {
    pub fn new(index: u32, excitation_energy_ev: f64) -> Self {
        ExcitedState { index, excitation_energy_ev, ..ExcitedState::default() }
    }
}

/// An `excited_states` section (spec § 5.12).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ExcitedStates {
    pub id: String,
    pub states: Vec<ExcitedState>,
    /// The `$defs/ExcitedStatesMeta` object (`method`, `n_states`, `tda`, …).
    pub meta: Option<Map<String, Value>>,
    /// The `wavefunction.gto` of the reference state, added earlier.
    pub wavefunction_ref: Option<String>,
    pub label: Option<String>,
}

// -- writing -----------------------------------------------------------------

impl<W: Write> Writer<W> {
    /// Add a `structure` section with its `structure` and, when it has bonds,
    /// `bonds` members. Returns the section id.
    pub fn add_structure(&mut self, s: &Structure) -> Result<String, WriteError> {
        s.check()?;
        let fields = fields_from(vec![
            ("label", s.label.clone().map(Value::from)),
            ("charge", s.charge.map(Value::from)),
            ("multiplicity", s.multiplicity.map(Value::from)),
            ("name", s.name.clone().map(Value::from)),
        ]);
        self.add_section(&s.id, "structure", fields)?;
        self.add_json(&s.id, "structure", &s.payload())?;
        if !s.bonds.is_empty() {
            let pairs: Vec<Value> = s
                .bonds
                .iter()
                .map(|b| match b.order {
                    Some(o) => json!({"i": b.i, "j": b.j, "order": o}),
                    None => json!({"i": b.i, "j": b.j}),
                })
                .collect();
            self.add_json(&s.id, "bonds", &json!({ "pairs": pairs }))?;
        }
        if let Some(info) = self.info.get_mut(&s.id) {
            info.atoms = Some(s.atoms.len());
        }
        Ok(s.id.clone())
    }

    /// Add a `wavefunction.gto` section: `basis`, `mo_metadata`, the
    /// coefficient tensors (float64) and the optional AO matrices.
    pub fn add_wavefunction_gto(&mut self, wf: &WavefunctionGto) -> Result<String, WriteError> {
        let what = format!("wavefunction.gto {:?}", wf.id);
        let atoms = self.require_kind(&format!("{what}: basis structure_ref"), &wf.basis.structure_ref, &["structure"])?.atoms;
        let n_ao = wf.basis.n_ao();
        for (k, s) in wf.basis.shells.iter().enumerate() {
            if s.exponents.is_empty() || s.exponents.len() != s.coefficients.len() {
                return invalid(format!("{what}: shell {k} needs as many coefficients as exponents (at least one)"));
            }
            if s.exponents.iter().any(|&a| !(a.is_finite() && a > 0.0)) {
                return invalid(format!("{what}: shell {k} exponents must be positive"));
            }
            finite(&format!("{what}: shell {k} coefficients"), &s.coefficients)?;
            if atoms.is_some_and(|n| s.center >= n) {
                return invalid(format!("{what}: shell {k} centre {} is not an atom of {:?}", s.center, wf.basis.structure_ref));
            }
        }
        let spin = match &wf.orbitals {
            SpinOrbitals::Restricted(_) => "restricted",
            SpinOrbitals::RestrictedOpen(_) => "restricted_open",
            SpinOrbitals::Unrestricted { .. } => "unrestricted",
        };
        let mut meta = Map::new();
        let mut coefficients: Vec<(&str, &Orbitals, usize)> = Vec::new();
        let mut rows = std::collections::HashMap::new();
        match &wf.orbitals {
            SpinOrbitals::Restricted(o) | SpinOrbitals::RestrictedOpen(o) => {
                let n = o.rows(n_ao, &format!("{what} mo_coefficients"))?;
                meta.extend(o.block(n));
                meta.insert("n_ao".into(), Value::from(n_ao));
                coefficients.push(("mo_coefficients", o, n));
                rows.insert(String::new(), n);
            }
            SpinOrbitals::Unrestricted { alpha, beta } => {
                meta.insert("n_ao".into(), Value::from(n_ao));
                for (name, o) in [("alpha", alpha), ("beta", beta)] {
                    let n = o.rows(n_ao, &format!("{what} mo_coefficients_{name}"))?;
                    meta.insert(name.into(), Value::Object(o.block(n)));
                    rows.insert(name.to_owned(), n);
                }
                coefficients.push(("mo_coefficients_alpha", alpha, rows["alpha"]));
                coefficients.push(("mo_coefficients_beta", beta, rows["beta"]));
            }
        }
        meta.insert("spin".into(), Value::from(spin));
        meta.insert("orbital_kind".into(), Value::from(wf.orbital_kind.as_str()));
        if let Some(nto) = &wf.nto {
            if wf.orbital_kind != "nto" {
                return invalid(format!("{what}: an NTO link needs orbital_kind \"nto\""));
            }
            if nto.state_index == 0 {
                return invalid(format!("{what}: state_index counts from 1"));
            }
            meta.insert("state_ref".into(), Value::from(nto.state_ref.as_str()));
            meta.insert("state_index".into(), Value::from(nto.state_index));
            let role = match nto.role {
                NtoRole::Hole => "hole",
                NtoRole::Particle => "particle",
            };
            meta.insert("nto_role".into(), Value::from(role));
        }
        if wf.orbital_kind == "nto" {
            meta.insert("occupation_semantics".into(), Value::from("transition_weight"));
        }
        if let (Some(n), true) = (wf.basis.n_electrons, wf.orbital_kind != "nto") {
            let total: f64 = coefficients.iter().flat_map(|(_, o, _)| o.occupations.iter()).sum();
            if (total - n as f64).abs() > 1e-6 {
                return invalid(format!("{what}: n_electrons is {n} but the occupations sum to {total} (spec § 5.3)"));
            }
        }
        let mut matrices: Vec<(&str, &Vec<f64>)> = Vec::new();
        if let Some(s) = &wf.overlap {
            matrices.push(("overlap", s));
        }
        match (&wf.orbitals, &wf.density_matrix, &wf.density_matrices_ab) {
            (SpinOrbitals::Unrestricted { .. }, Some(_), _) => {
                return invalid(format!("{what}: an unrestricted wavefunction takes density_matrices_ab, not density_matrix"))
            }
            (SpinOrbitals::Restricted(_) | SpinOrbitals::RestrictedOpen(_), _, Some(_)) => {
                return invalid(format!("{what}: density_matrices_ab needs unrestricted orbitals"))
            }
            (_, Some(p), _) => matrices.push(("density_matrix", p)),
            (_, _, Some((a, b))) => {
                matrices.push(("density_matrix_alpha", a));
                matrices.push(("density_matrix_beta", b));
            }
            _ => {}
        }
        for (role, m) in &matrices {
            if m.len() != n_ao * n_ao {
                return invalid(format!("{what}: {role} must be [{n_ao}, {n_ao}], got {} values", m.len()));
            }
            finite(&format!("{what} {role}"), m)?;
        }
        let fields = fields_from(vec![("label", wf.label.clone().map(Value::from))]);
        self.add_section(&wf.id, "wavefunction.gto", fields)?;
        self.add_json(&wf.id, "basis", &wf.basis.payload())?;
        self.add_json(&wf.id, "mo_metadata", &Value::Object(meta))?;
        for (role, o, n) in coefficients {
            self.add_tensor(&wf.id, role, &o.coefficients, &[n, n_ao])?;
        }
        for (role, m) in matrices {
            self.add_tensor(&wf.id, role, m, &[n_ao, n_ao])?;
        }
        if let Some(info) = self.info.get_mut(&wf.id) {
            info.mo_rows = rows;
        }
        Ok(wf.id.clone())
    }

    fn start_volume(&mut self, v: &Volume, dtype: DType) -> Result<Vec<usize>, WriteError> {
        let what = format!("{} {:?}", v.kind, v.id);
        if !matches!(dtype, DType::Float32 | DType::Float64) {
            return invalid(format!("{what}: volume data is float32 or float64"));
        }
        let shape = v.data_shape()?;
        finite(&format!("{what} origin"), &v.origin)?;
        finite(&format!("{what} voxel_vectors"), &v.voxel_vectors.concat())?;
        if v.axis.is_some() != (v.kind == "volume.series") {
            return invalid(format!("{what}: an axis belongs to volume.series only"));
        }
        if let Some(a) = &v.axis {
            if !["time", "frame", "state", "custom"].contains(&a.kind.as_str()) {
                return invalid(format!("{what}: axis kind {:?} is not time, frame, state or custom", a.kind));
            }
            finite(&format!("{what} axis values"), &a.values)?;
            if a.kind == "frame" && a.values.iter().any(|&x| x < 0.0 || x.fract() != 0.0) {
                return invalid(format!("{what}: frame axis values are 0-based frame indices"));
            }
            if a.labels.as_ref().is_some_and(|l| l.len() != a.values.len()) {
                return invalid(format!("{what}: axis needs one label per slab"));
            }
            if let Some(t) = &a.trajectory_ref {
                self.require_kind(&format!("{what}: axis trajectory_ref"), t, &["trajectory", "reaction.path"])?;
            }
        }
        if let Some(s) = &v.structure_ref {
            self.require_kind(&format!("{what}: structure_ref"), s, &["structure"])?;
        }
        if v.operands.is_some() && v.kind != "volume.difference" {
            return invalid(format!("{what}: operands belong to volume.difference"));
        }
        if v.excited_state.is_some() && v.kind != "volume.transition_density" {
            return invalid(format!("{what}: an excited state belongs to volume.transition_density"));
        }
        let mut fields =
            fields_from(vec![("label", v.label.clone().map(Value::from)), ("component", v.component.clone().map(Value::from))]);
        if let Some((a, b)) = &v.operands {
            fields.insert("operand_a".into(), Value::from(a.as_str()));
            fields.insert("operand_b".into(), Value::from(b.as_str()));
        }
        if let Some((es, index)) = &v.excited_state {
            fields.insert("excited_states_ref".into(), Value::from(es.as_str()));
            if let Some(i) = index {
                if *i == 0 {
                    return invalid(format!("{what}: state_index counts from 1"));
                }
                fields.insert("state_index".into(), Value::from(*i));
            }
        }
        self.add_section(&v.id, &v.kind, fields)?;
        self.add_json_at(&v.id, "grid", &v.grid_payload(), &format!("volumes/{}_grid.json", v.id))?;
        if let Some(axis) = v.axis_payload() {
            self.add_json_at(&v.id, "axis", &axis, &format!("volumes/{}_axis.json", v.id))?;
        }
        Ok(shape)
    }

    /// Add a `volume.*` section with its `grid` (and `axis`) and its `data`
    /// (`f32` or `f64`, C order, shape [`Volume::data_shape`]).
    pub fn add_volume<T: Element>(&mut self, v: &Volume, data: &[T]) -> Result<String, WriteError> {
        let shape = v.data_shape()?;
        if data.len() as u64 != shape.iter().map(|&n| n as u64).product::<u64>() {
            return invalid(format!("{} {:?}: {} values do not fill the data shape {shape:?}", v.kind, v.id, data.len()));
        }
        let mut sink = self.begin_volume(v, T::DTYPE)?;
        sink.write(data)?;
        sink.finish()?;
        Ok(v.id.clone())
    }

    /// Add a `volume.*` section and return the sink of its `data`, to be
    /// written in slabs (C order) and finished.
    pub fn begin_volume(&mut self, v: &Volume, dtype: DType) -> Result<TensorSink<'_, W>, WriteError> {
        let shape = self.start_volume(v, dtype)?;
        let path = format!("volumes/{}.dat", v.id);
        self.begin_tensor_at(&v.id, "data", dtype, &shape, &path)
    }

    /// Add an `excited_states` section with its `states` member.
    pub fn add_excited_states(&mut self, es: &ExcitedStates) -> Result<String, WriteError> {
        let what = format!("excited_states {:?}", es.id);
        let rows = match &es.wavefunction_ref {
            Some(r) => Some(self.require_kind(&format!("{what}: wavefunction_ref"), r, &["wavefunction.gto"])?.mo_rows.clone()),
            None => None,
        };
        let mut states = Vec::with_capacity(es.states.len());
        let mut last = 0u32;
        for (k, s) in es.states.iter().enumerate() {
            let w = format!("{what}: states[{k}]");
            if s.index == 0 || s.index <= last {
                return invalid(format!("{w}: index {} must be at least 1 and above the previous one (spec § 5.12)", s.index));
            }
            last = s.index;
            let mut o = Map::new();
            o.insert("index".into(), Value::from(s.index));
            finite(&w, &[s.excitation_energy_ev])?;
            o.insert("excitation_energy".into(), json!({"value": s.excitation_energy_ev, "units": "ev"}));
            if let Some(e) = s.total_energy_hartree {
                finite(&w, &[e])?;
                o.insert("total_energy".into(), json!({"value": e, "units": "hartree"}));
            }
            if let Some(m) = s.multiplicity {
                if m == 0 {
                    return invalid(format!("{w}: multiplicity must be at least 1"));
                }
                o.insert("multiplicity".into(), Value::from(m));
            }
            if let Some(sym) = &s.symmetry {
                o.insert("symmetry".into(), Value::from(sym.as_str()));
            }
            if let Some(f) = s.oscillator_strength {
                if !(f.is_finite() && f >= 0.0) {
                    return invalid(format!("{w}: oscillator_strength must be finite and non-negative"));
                }
                o.insert("oscillator_strength".into(), Value::from(f));
            }
            if let Some(r) = s.rotatory_strength_cgs {
                finite(&w, &[r])?;
                o.insert("rotatory_strength_cgs".into(), Value::from(r));
            }
            for (key, v) in
                [("transition_dipole", s.transition_dipole), ("transition_magnetic_dipole", s.transition_magnetic_dipole)]
            {
                if let Some(v) = v {
                    finite(&format!("{w}.{key}"), &v)?;
                    o.insert(key.into(), json!(v));
                }
            }
            if !s.contributions.is_empty() {
                let mut list = Vec::with_capacity(s.contributions.len());
                for (j, c) in s.contributions.iter().enumerate() {
                    finite(&format!("{w}.contributions[{j}].weight"), &[c.weight])?;
                    let channel = match c.spin {
                        Some(Channel::Alpha) => "alpha",
                        Some(Channel::Beta) => "beta",
                        None => "",
                    };
                    if let Some(rows) = &rows {
                        let n = rows.get(channel).or_else(|| if channel.is_empty() { None } else { rows.get("") });
                        match n {
                            Some(&n) if c.from_mo < n && c.to_mo < n => {}
                            Some(&n) => {
                                return invalid(format!(
                                    "{w}.contributions[{j}]: orbitals {} -> {} outside the {n} rows",
                                    c.from_mo, c.to_mo
                                ))
                            }
                            None => {
                                return invalid(format!(
                                    "{w}.contributions[{j}]: the reference has no {} orbitals",
                                    if channel.is_empty() { "single-set" } else { channel }
                                ))
                            }
                        }
                    }
                    let mut co = json!({"from_mo": c.from_mo, "to_mo": c.to_mo, "weight": c.weight});
                    if !channel.is_empty() {
                        co["spin"] = Value::from(channel);
                    }
                    list.push(co);
                }
                o.insert("contributions".into(), Value::Array(list));
            }
            if let Some(r) = &s.nto_ref {
                self.require_kind(&format!("{w}.nto_ref"), r, &["wavefunction.gto"])?;
                o.insert("nto_ref".into(), Value::from(r.as_str()));
            }
            if let Some(r) = &s.density_ref {
                self.require_kind(&format!("{w}.density_ref"), r, &["volume.transition_density", "volume.difference"])?;
                o.insert("density_ref".into(), Value::from(r.as_str()));
            }
            states.push(Value::Object(o));
        }
        let mut payload = json!({ "states": states });
        if let Some(meta) = &es.meta {
            if let Some(n) = meta.get("n_states") {
                if n.as_u64() != Some(es.states.len() as u64) {
                    return invalid(format!("{what}: meta.n_states {n} but {} states (spec § 5.12)", es.states.len()));
                }
            }
            payload["meta"] = Value::Object(meta.clone());
        }
        let fields = fields_from(vec![
            ("label", es.label.clone().map(Value::from)),
            ("wavefunction_ref", es.wavefunction_ref.clone().map(Value::from)),
        ]);
        self.add_section(&es.id, "excited_states", fields)?;
        self.add_json(&es.id, "states", &payload)?;
        if let Some(info) = self.info.get_mut(&es.id) {
            info.states = Some(es.states.iter().map(|s| s.index).collect());
        }
        Ok(es.id.clone())
    }
}
