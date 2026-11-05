# QVF v1 review — gaps, defects, and proposed changes

> Status: planning document (informative). Nothing here changes the contract
> by itself. Each item is either a **defect** to fix in place (D-nn), a
> **proposal** for an additive v1 change (P-nn), or a candidate for a future
> **major version** (V2-nn). Work packages in `work_packages.md` reference
> these ids; the roadmap in `ROADMAP.md` orders them.

The review was made against the repository state of 2026-10-06: spec
`spec/qvf-format-spec.md`, schema `spec/qvf_manifest.schema.json`, the Python
and C++ reference implementations, the C ABI, the conformance corpus, and the
docs. Baseline at review time: Python tests pass, 12/12 corpus archives pass,
the C++ library builds and its tests pass on macOS/clang.

## 1. Coverage matrix — what QVF v1 can carry today

Legend: ✅ canonical and usable · ⚠️ possible but under-specified or defective ·
❌ no canonical carrier (vendor namespace only) · — out of scope by design.

Refreshed 2026-10-10, after Phases 0–2 and WP-40 / WP-41 landed. The defect and
proposal sections below are kept as the record of what was found and decided.

Current reconciliation (2026-10-19): the scheduled format changes are
implemented, but this historical list is not a blanket completion claim.
D09 hosting remains open (the `$id` retention decision is settled); broader
packaging, Fortran/platform qualification and fuzzing remain in the
[work packages](work_packages.md). P13 and P16 retain demand-driven kinds.

| Data category | Status | Notes / id |
|---|---|---|
| **Geometry** | | |
| Molecule, periodic 0–3-D, slabs/wires | ✅ | `structure`, `pbc` contract is strong |
| Biomolecule metadata (chains, residues, SS, B-factors) | ✅ | |
| Several structures in one archive (ensembles, workflow steps) | ✅ | generic `structure_ref` / `run_ref` / `group` (P04, WP-20) |
| Forces / gradient, stress | ✅ | `forces` (P06, WP-21) |
| Velocities | ✅ | `trajectory` / `reaction.path` members (P12, WP-26) |
| Symmetry summary | ✅ | payload schema'd (P01, WP-11) |
| Constraints | ✅ | root block |
| **Energies** | | |
| Total SCF energy, Fermi energy | ✅ | `provenance` |
| Energy components, several levels of theory | ✅ | root `energies` (P05, WP-21) |
| Thermochemistry | ✅ | root block |
| **Wavefunction** | | |
| GTO molecular, Γ-point periodic, R/U | ✅ | Appendix A is the strongest part of the spec |
| ROHF / GHF / two-component / complex coefficients | ✅ | spin layouts and `complex` (P09, WP-23) |
| ECP / effective core electrons | ✅ | `basis.ecp`, `n_electrons` (P09, WP-23) |
| AO overlap, density matrix, natural occupations | ✅ | `overlap` and density members; normalization checkable (P09, WP-23) |
| NTOs, transition densities | ✅ | NTO pairs in `wavefunction.gto`, `volume.transition_density` (P08, WP-22) |
| k-resolved Bloch orbitals | ✅ | `wavefunction.bloch`, provisional (P15, WP-30) |
| Slater (STO) and numeric (NAO) bases | ❌ | P16, demand-driven |
| Plane waves | — | by design: visualise via `volume.*` |
| **Scalar fields** | | |
| Density, orbital, spin, ELF, potential, RDG, difference, generic | ✅ | 8 kinds |
| Vector fields (current density, field, gradient) | ✅ | `volume.vector` (P11, WP-25) |
| Complex fields | ✅ | `volume.complex` (P11, WP-25) |
| Field series (time / states / path) | ✅ | `volume.series` (P11, WP-25) |
| Periodic grid semantics (endpoint included or not) | ✅ | `grid.periodic`, `includes_endpoint`, span checked (P11, WP-25) |
| Huge grids (chunked, streamed, better compression) | ⚠️ | ZIP64 and streaming write in both writers (WP-40); chunked members stay V2-02 |
| **Dynamics** | | |
| Optimisation / IRC / MD frames | ✅ | coords only |
| Periodic MD (per-frame cell) | ✅ | `trajectory.lattice`, rows in Å (P12, D06, WP-26) |
| Velocities, forces, time axis, per-frame energies as binary | ✅ | P12, WP-26 |
| NEB / reaction paths with waypoints | ✅ | |
| 2-D relaxed scans | ✅ | units declared (D05, WP-13) |
| **Vibrations** | | |
| Normal modes + frequencies | ✅ | metadata schema'd (P01, WP-10) |
| Hessian | ✅ | `hessian` (P07, WP-21) |
| IR / Raman / VCD spectra | ✅ | typed peaks and broadening (P13, WP-11) |
| **Excited states** | | |
| UV-Vis / ECD spectra | ✅ | |
| State table (energies, f, R, transition dipoles, contributions) | ✅ | `excited_states` (P08, WP-22) |
| **Properties** | | |
| Dipole | ✅ | root block |
| Quadrupole and higher multipoles, polarisabilities | ✅ | `properties.electric` (P17, WP-29) |
| Atomic charges / spin populations | ✅ | one role list, vector roles, `custom.<name>` (D02, P10, WP-24) |
| Bond orders, QTAIM, NCI | ✅ | QTAIM positions in Å, declared (D05, WP-11) |
| NMR, EPR | ✅ | typed payloads (P13, WP-11) |
| Mössbauer, XAS/XES | ✅ | `spectra.xas`, `spectra.xes`, `spectra.mossbauer` (P13, #19, Unreleased) |
| Other spectroscopies | ❌ | `spectra.generic` or vendor; P13 adds typed kinds on demand |
| **Periodic analysis** | | |
| Bands, DOS, PDOS, COOP/COHP | ✅ | spin axis on projections, `occupations` (D04, WP-12) |
| Fermi surface | ✅ | mesh payload schema'd, `shape` optional (P14, WP-12) |
| Phonon bands / DOS | ✅ | optional members schema'd (D03, WP-04) |
| Equation of state | ✅ | |
| Elastic, dielectric, Born charges, stress | ✅ | `properties.elastic`, `properties.dielectric`, `forces.stress` (P06, P14, WP-21, WP-29) |
| Non-collinear magnetic moments | ✅ | `atom_properties` vector roles (P10, WP-24) |
| **Provenance and container** | | |
| Source, provenance, citations, run record, job spec | ✅ | |
| FAIR metadata (title, authors, DOI, licence, uuid, created) | ✅ | root `metadata` (P03, WP-20) |
| Privacy of shared archives | ✅ | `provenance.hostname` deprecated, sanitisation guidance (D11, WP-07) |
| Random access, per-member sha256 | ✅ | |
| ZIP64 (> 4 GiB) | ✅ | both readers and writers (D13, WP-40) |
| Streaming write, in-place update | ✅ | C++ streaming writer, `open_for_update`, `repack`; atomic editor with WP-42 (D19, WP-40) |
| Byte-order tag | ✅ | `byte_order` (P20, WP-28) |

## 2. Defects and inconsistencies (fix in v1, no version bump)

### D01 — Schema descriptions reference internal documents and stale section numbers
Several `description` strings in the schema cite `design_qvf_format.md § 4.x`,
`vibeqc.PeriodicSystem.lattice`, and "vibe-view loads it from a symlink". That
document is not in this repository and the implementation detail is not part
of the contract. Other descriptions cite spec sections by numbers that no
longer match (`§ 4.7` for root metadata, which is § 3.3; `§ 5.4` / `§ 5.5` for
extensions and `critical`, which are § 7.4).
**Fix.** Scrub all descriptions to self-contained text; cite the spec by
section *title* and a stable anchor; add a test that every `§ n.n` in the
schema resolves to a heading in the spec.

### D02 — `atom_properties` role list disagrees across spec, schema, and writers
Spec § 5.7 lists three roles, the § 5 table four (`hirshfeld_charge`), the
schema five (`iao_charge`). The Python writer and the C ABI expose three.
**Fix.** One list, in the schema, mirrored verbatim in spec § 5.7 and both
writers; extend per P10.

### D03 — Writers emit schema-invalid optional phonon members
`QvfWriter.add_phonon_bands(eigenvectors=…)` writes an `eigenvectors` member
and `add_phonon_dos(projected=…)` a `projected` member. Both section schemas
set `additionalProperties: false` without those roles, so the archive fails
schema validation (verified 2026-10-06). The C++ writer has the same
parameters.
**Fix.** Add the two optional members to the schema with dtype/shape
constraints (`eigenvectors`: float64 `[n_q, n_modes, n_atoms, 3]` real, or
`[…, 2]` complex per P09's convention; `projected`: float64
`[n_channels, n_points]`), document them in spec § 5.6, add a corpus archive.

### D04 — Band and DOS projections carry no spin axis
`bands.eigenvalues` is `[n_spin, n_k, n_bands]` but `projections` is
`[n_k, n_bands, n_channels]`; `dos.projected.projections` is
`[n_channels, n_points]` with spin only as a peer `n_spin` field.
**Fix (additive).** Accept `[n_spin, n_k, n_bands, n_channels]` and
`[n_spin, n_channels, n_points]`; the existing lower-rank shapes remain legal
and mean `n_spin = 1`. The validator cross-checks `n_spin` against
`eigenvalues`.

### D05 — Undeclared or mixed units in analysis payloads
`topology.qtaim` positions have no declared unit. `scan.surface` axes are
"bohr for bonds, radians for angles" with geometries in Å, and nothing in the
payload says which.
**Fix.** Declare the QTAIM unit (decision needed; recommended: Å, matching the
structure the points are overlaid on) and add an optional `units` object to
both payloads under P02.

### D06 — Two lattice conventions in one format
`structure.lattice_vectors` are row vectors in Å; `reaction.path.lattice` is a
column matrix in bohr. `trajectory` has no lattice at all.
**Fix.** Units are frozen within v1, so: document the discrepancy loudly in
§ 4.3, give `trajectory` a `lattice` member that follows the *structure*
convention (rows, Å) so new carriers converge, provide a `Lattice` helper in
the library that converts, and schedule unification for v2 (V2-05).

### D07 — Only two of about twenty JSON member payloads are schema'd
`StructurePayload` and `JobSpecPayload` are specified. Everything else —
volume `grid`, wavefunction `basis` and `mo_metadata`, every `metadata`,
`waypoints`, every `spectrum`, `kpath`, `qpath`, `mesh`, `meta`, `fit`,
`data` (symmetry), `iterations`, `bonds`, `bond_orders`, `critical_points`,
run-record `files` — is prose only. Two producers can diverge without any
validator noticing; this is the largest obstacle to "a standard external codes
can adopt".
**Fix.** P01.

### D08 — Section peer fields emitted but unschema'd
`dos.total` / `dos.projected` carry `n_spin`, `smearing`, `smearing_type`,
`fermi_energy_ev`, `n_electrons`, `channels` on the section object. The
`Section` definition is open, so they pass, but they are not a contract.
**Fix.** Name them in `SectionDOSTotal` / `SectionDOSProjected` (or move
`channels` into a `meta` member for symmetry with `dos.coop`); decide whether
`Section` should become closed per branch (it should, once every peer field is
named).

### D09 — Schema `$id` domain does not match the project domain and does not resolve
`$id` is `https://vibe-qc.org/…`; the project, registry, and packaging use
`vibe-qc.com`. Neither URL serves the schema.
**Decision update (2026-10-15).** Retain the existing `$id`; do not change
identifiers to implement the original domain proposal. Hosting remains open
under WP-06. The preceding URL observation is historical, not a fresh
availability probe.

### D10 — Referenced artefacts that do not exist in the repository
`GOVERNANCE.md` and `registry.json` state that `qvf_manifest_v2.schema.json`
"is retained as a frozen artifact" and that decisions are "recorded in the
CHANGELOG". Neither file exists here.
**Fix.** Add `CHANGELOG.md` (Keep-a-Changelog style, seeded from the decision
log) and `spec/legacy/qvf_manifest_v2.schema.json` (a strict superset of v1,
generated once), or amend the two documents.

### D11 — `provenance.hostname` is a recognised field
A format meant for sharing should not encourage machine names in the file;
this repository's own policy forbids them in examples.
**Fix.** Deprecate `hostname` in spec § 3.2 (consumers ignore; producers
SHOULD NOT emit), writers never emit it by default, and the toolkit gains a
`scrub` operation that also treats `run.record` text (absolute paths, user
names) with a reviewable report. See P18.

### D12 — `ValueWithUnits.units` is free text
`"Eh"`, `"hartree"`, `"Hartree"` are all accepted and all different to a
consumer.
**Fix.** P02 defines a small controlled vocabulary; the validator warns on
unknown strings.

### D13 — C++ ZIP writer and reader lack ZIP64
The spec permits members up to 8 GiB; Python's `zipfile` silently switches to
ZIP64 for members ≥ 4 GiB (and for > 65 535 entries). The C++ reader cannot
open those archives and the C++ writer cannot produce them. This is an
interoperability break between the two reference implementations.
**Fix.** ZIP64 in both (extra field 0x0001, ZIP64 end-of-central-directory
record and locator); corpus gets a synthetic ZIP64 archive (small payload,
ZIP64 structures forced).

### D14 — Portability and CI
`primitive_norm` uses `M_PI`, which is not standard C++ (fails on MSVC without
`_USE_MATH_DEFINES`). CI runs only on Debian/gcc and one Python; the public
mirror on GitHub has no workflow at all.
**Fix.** A `constexpr` pi; GitHub Actions matrix (Linux gcc/clang, macOS,
Windows MSVC; Python 3.9–3.14); sanitizers on one Linux job.

### D15 — Validator hardening gaps for untrusted input
Not checked today: duplicate ZIP entry names (readers disagree on which wins),
backslashes and drive letters in member paths, JSON nesting depth and size,
total uncompressed expansion across members (zip-bomb by count), and orphan
members (ZIP entries no manifest member references).
**Fix.** Add all five to spec § 6 and both validators; orphans are a warning
(they are the by-product of in-place updates, P19), the rest are errors.

### D16 — The C ABI covers a third of the format
`qvf_add_structure` has no `pbc`/lattice; the wavefunction entry is restricted
only; no volumes, trajectories, vibrations, bands, DOS, periodic kinds, no
binary vendor members, and no reader or validator at all. Most production
quantum-chemistry codes are Fortran and can only use this surface.
**Fix.** Library design § 6 (complete ABI) and WP-44 / WP-45.

### D17 — `viewer_defaults` key space collision
Per-section hints are keyed by section id in the same object that holds
`auto_open` and `bookmarks`; a section with either id is ambiguous.
**Fix.** Reserve the two names (validator error on a section with those ids);
v2 moves hints under `sections: {…}` (V2-08).

### D18 — Spectrum x-axis is called `frequencies` even when it is an energy
`spectra.uvvis` uses `frequencies` for eV values; the units table also mentions
`energies_ev`, which no writer emits.
**Fix.** Keep `frequencies` for compatibility and add `x_units` / `y_units`
and an optional explicit `x` alias in the spectrum payload schema (P13); v2
renames (V2-06).

### D19 — The writers are in-memory only; there is no update path
Both writers accumulate every member in memory and write once. The job-spec
lifecycle ("the runner updates the same container") and live checkpoints have
no toolkit support, and a 1-gigavoxel grid must exist twice in memory.
**Fix.** Library design § 4 (streaming sinks, ZIP data descriptors, editor).

### D20 — Python distribution shape
Two flat modules (`qvf_writer`, `qvf_reader`), the distribution is named
`qvf-writer` although it also reads and validates, the reader has no typed
accessors, and the only CLI is `qvf-validate`.
**Fix.** Library design § 7; compatibility shims keep the old import names.

## 3. Proposals — additive v1 changes

Each proposal is additive (new optional members, kinds, or fields), so by the
governance rules none bumps `qvf_version`. Each lands as spec + schema +
registry + both writers + validator + corpus archive + docs, with a decision
record in `GOVERNANCE.md`.

### P01 — Payload schemas for every JSON member
Add a `$defs/<Kind><Role>Payload` definition for every JSON member role, and
make the validator load and check each one (the mechanism `StructurePayload`
already uses). Three strictness levels, chosen per payload:

| Level | Applies to | Rule |
|---|---|---|
| strict | `grid`, `basis`, `mo_metadata`, `kpath`, `qpath`, `mesh`, `fit`, `iterations`, `bonds`, `bond_orders`, `waypoints`, `files` | `additionalProperties: false` except an `x_*` escape key |
| typed-open | `metadata` of trajectory/reaction/scan/vibrations, DOS `meta`, phonon `meta`, 1-D `spectrum` | known fields typed, root open |
| permissive | `spectra.nmr`, `spectra.epr`, `structure.symmetry.data`, `critical_points` | typed known fields, open, tightened later by P13 |

Where the reference producer already emits a field the schema must accept it
unchanged; the corpus is regenerated and must still pass. Deliverable also
includes a generated "payload catalogue" page in the docs.

### P02 — Units vocabulary and self-describing payloads
Spec § 4.3 gains a controlled list of unit strings:
`angstrom`, `bohr`, `hartree`, `ev`, `kcal_per_mol`, `kj_per_mol`,
`cm-1`, `nm`, `debye`, `mhz`, `ppm`, `km_per_mol`, `a3` (Å³), `gpa`, `fs`,
`ps`, `kelvin`, `atm`, `e` (elementary charge), `dimensionless`.
`ValueWithUnits.units` becomes an enum of these; payloads that carry
geometry-like data (`grid`, `critical_points`, `scan.surface.metadata`,
`trajectory.metadata`) gain an optional `units` object naming the unit of each
numeric field, defaulting to the § 4.3 convention. Producers MAY declare,
consumers MUST honour a declared unit. This is how v1 becomes
self-describing without changing any default.

### P03 — Root `metadata` block (FAIR)
Optional root object: `title`, `description`, `created_utc` (ISO 8601),
`modified_utc`, `uuid` (RFC 4122), `authors[]` (`{name, orcid?, affiliation?}`),
`license` (SPDX id), `doi`, `keywords[]`, `derived_from[]` (uuids or DOIs),
`tags`. Purely descriptive; writers set `created_utc` and `uuid` by default.

### P04 — Generic section references and grouping
Every section object MAY carry `structure_ref` (id of the `structure` it
describes), `run_ref` (id of the `run.record` that produced it), and `group`
(free string for UI grouping, e.g. `"step 2: frequencies"`). The validator
resolves refs. Existing specific refs (`operand_a`, `trajectory_ref`, and
`basis.json.structure_ref`) stay; `structure_ref` on the section object takes
precedence when both are present. This makes multi-structure and multi-step
archives unambiguous.

### P05 — Root `energies` block
`{ "total": {value, units}, "components": { "nuclear_repulsion", "electronic",
"one_electron", "two_electron", "exchange_correlation", "correlation",
"dispersion", "solvation", "zpve", … }, "by_method": [ {method, total,
components?} ] }`. `provenance.scf_energy` remains and must equal
`energies.by_method[method == provenance.method].total` when both exist
(validator check).

### P06 — `forces` kind
Members: `gradient` float64 `[n_atoms, 3]` in Eh/bohr (sign convention:
gradient, not force; `forces = -gradient`); optional `stress` float64 `[3, 3]`
in Eh/bohr³ (periodic); JSON `meta` with `structure_ref`, `method`, optional
`units`. Also used per frame via P12.

### P07 — `hessian` kind
Members: `matrix` float64 `[3N, 3N]` in Eh/bohr²; JSON `meta` with
`structure_ref`, `mass_weighted` (bool), `projected` (translations/rotations
removed), `atom_masses_amu[]`. `vibrations` MAY carry `hessian_ref`.

### P08 — `excited_states` kind
JSON `states` member: `[{index, energy_ev, excitation_energy_ev,
multiplicity, symmetry?, oscillator_strength?, rotatory_strength?,
transition_dipole_au?, transition_magnetic_dipole_au?, contributions:
[{from_mo, to_mo, spin?, weight}], nto_ref?, density_ref?}]` plus `meta`
(`method`, `reference_state`, `n_states`, `tda`). `spectra.uvvis` /
`spectra.ecd` payloads MAY reference states by index (P13). Covers TD-DFT,
CIS, ADC, EOM-CC outputs.

### P09 — Wavefunction extensions
`basis.json`: optional `ecp` per centre (`{n_core_electrons, label,
projectors?: [{l, exponents, coefficients, powers}]}`); `n_electrons`.
`mo_metadata.json`: `orbital_kind` enum grows to `canonical | natural |
localized | nto | nbo | iao | ibo | boys | pipek_mezey | custom`;
`symmetry_labels[]`; `state_ref` (for NTOs, P08); `spin` grows to
`restricted | unrestricted | restricted_open | generalized`;
`complex: true` documents the `[n_mo, n_ao, 2]` layout the schema already
allows. Optional binary members: `overlap` float64 `[n_ao, n_ao]` (lets a
consumer verify the Appendix A normalization numerically),
`density_matrix` / `density_matrix_alpha` / `density_matrix_beta`
`[n_ao, n_ao]`. Spec § 5.3 gains the two-component description.

### P10 — `atom_properties` expansion
Charge roles: add `cm5_charge`, `chelpg_charge`, `resp_charge`, `esp_charge`,
`nbo_charge`, `mbis_charge`, `bader_charge`, `ddec_charge`. Other scalar roles:
`valence`, `bonded_valence`, `fukui_plus`, `fukui_minus`, `fukui_zero`,
`dual_descriptor`, `atomic_energy_eh`, `atomic_volume_bohr3`,
`nmr_shielding_iso_ppm`. Vector roles (`AtomVectorBinary`, `[n_atoms, 3]`):
`magnetic_moment`, `force_eh_bohr` (alias of P06 for per-atom display).
Escape: `custom.<name>` roles allowed by `patternProperties` with a required
`meta` JSON member describing `{units, description}` per custom role.

### P11 — Volume additions and periodic grid semantics
New kinds: `volume.vector` (`[nx, ny, nz, 3]`), `volume.complex`
(`[nx, ny, nz, 2]`), `volume.series` (`[n, nx, ny, nz]` + `axis` JSON:
`{kind: time|frame|state|custom, values[], units}`),
`volume.transition_density`, `volume.laplacian`. `grid` JSON gains optional
`units` (`bohr` default), `periodic: [bool, bool, bool]`,
`includes_endpoint: bool` (default `true`, cube-style; CHGCAR-style producers
set `false`), `structure_ref`. The validator rejects `periodic` grids whose
voxel vectors do not span the referenced cell within tolerance.

### P12 — Trajectory extensions
Optional members on `trajectory` (and `reaction.path`): `lattice` float64
`[3, 3]` or `[n_frames, 3, 3]` **rows, Å** (structure convention; see D06),
`velocities` `[n_frames, n_atoms, 3]` Å/fs, `forces` `[n_frames, n_atoms, 3]`
Eh/bohr, `times` `[n_frames]` fs, `energies` binary `[n_frames]` Eh (preferred
over the JSON array for long runs), `cell_volumes`, `temperatures`.
`coords` MAY be float32 when `metadata.precision: "float32"` is set (MD).
`metadata.kind` enum: `optimization | irc | md | neb | scan | custom`;
`metadata.pbc`.

### P13 — Spectrum payload schema and typed NMR/EPR
Common 1-D payload: `frequencies[]`, `intensities[]`, `x_units`, `y_units`,
optional `x[]` (alias), `peaks[]` (`{x, y, width?, label?, assignment?:
{mode_index | state_index}}`), `broadening` (`{type: none|gaussian|lorentzian,
fwhm, units}`), `temperature_k`, kind-specific blocks (`raman: {laser_nm,
depolarization_ratios[]}`, `uvvis: {x_axis: energy_ev|wavelength_nm}`).
`spectra.nmr`: `chemical_shifts[]` typed, `shielding_tensors[]` `[3, 3]`,
`j_couplings[]` `{i, j, value_hz, isotope_i, isotope_j}`, `reference`,
`solvent`. `spectra.epr`: `g_tensor`, `hyperfine[]`, `zero_field_splitting`
typed with units fixed (MHz). New kinds on demand: `spectra.xas`,
`spectra.mossbauer`.

### P14 — Periodic analysis completeness
`bands`: optional `occupations` `[n_spin, n_k, n_bands]`, `kpath` schema
(`kpoints[]` fractional, `labels[]`, `segments[]`, `fermi`, `spinor: bool`),
spin-axis fix (D04). `dos.*`: `channels[]` typed (`{atom?, l?, m?, spin?,
label}`). `fermi_surface.mesh`: `shift[3]`, `symmetry_expanded`, `kpoints` for
non-uniform meshes. New demand-driven kinds: `properties.elastic`
(6×6 GPa), `properties.dielectric` (`[3, 3]` tensors, Born charges
`[n_atoms, 3, 3]`).

### P15 — `wavefunction.bloch` (provisional canonical)
Promote the shape of `x_vibeqc.bloch_wavefunction` once a second producer is in
sight (CP2K, CRYSTAL, or a PySCF-PBC adapter): `kpoints` JSON (fractional +
weights), `mo_coefficients` float64 `[n_k, n_spin, n_mo, n_ao, 2]`,
`energies` `[n_k, n_spin, n_mo]`, `occupations`, sharing `basis.json` with
`wavefunction.gto`. Marked "provisional" in the registry until the second
producer ships.

### P16 — `wavefunction.nao` and `wavefunction.sto` (demand-driven)
Numeric atomic orbitals (FHI-aims, SIESTA, DMol3): `basis` JSON (per function:
centre, l, radial table ref) + `radial_grids` / `radial_tables` binary.
Slater (ADF): shells with `{n, l, zeta}`. Both reuse `mo_metadata` and the
MO coefficient layout. Not scheduled until a producer asks.

### P17 — `properties.electric` kind
JSON: `quadrupole` (`[3, 3]`, Buckingham, traceless flag), `octupole`,
`polarizability` (`[3, 3]` au, frequency), `hyperpolarizability`, `origin`
(same convention as `dipole_moment.origin`). The root `dipole_moment` stays.

### P18 — Provenance and job-spec tightening
`provenance`: deprecate `hostname` (D11); add `software[]`
(`{name, version}`), `git_commit`, `n_threads`, `n_mpi_ranks`,
`started_utc`, `finished_utc`. `job.spec`: `structure_ref`, `solvation`
(`{model, solvent}`), `dispersion`, `relativistic`, `scf` (`{convergence,
max_iterations}`), `geometry` (`{convergence, max_steps}`), `resources`
(hints: `{n_cores, memory_mb, walltime_s}`).

### P19 — Container rules
Optional member field `size` (uncompressed byte count; validator checks it).
Duplicate ZIP entry names are an error. Orphan ZIP entries are a warning and
`repack` removes them. In-place update semantics are specified: an updating
producer MAY append members and MUST rewrite `manifest.json` last; the
previous manifest is superseded, not kept. `run.record.sequence` ordering is
reaffirmed for multi-step updates.

### P20 — Byte-order declaration
Optional root `byte_order` with `const: "little"` now (documents what every
archive already is; a validator rejects any other value). v2 moves the tag
into dtype names (V2-01).

## 4. v2 candidate ledger (breaking; collect, do not schedule)

| Id | Candidate | Why it cannot be v1 |
|---|---|---|
| V2-01 | Byte order in dtype names (`<f8`) | changes the dtype enum every consumer validates |
| V2-02 | Chunked members and alternative compression (zstd, byte-shuffle) | a v1 consumer cannot decode them; `critical` would be abused |
| V2-03 | dtype expansion: `complex64`, `complex128`, `float16`, `bool` | enum change |
| V2-04 | One length unit everywhere (Å *or* bohr) | units are frozen within a major |
| V2-05 | One lattice convention (rows, one unit) | same |
| V2-06 | Spectrum x-axis renamed `x` with mandatory `x_units` | required-member change |
| V2-07 | `structure_ref` required on every structure-bound section | required-member change |
| V2-08 | `viewer_defaults.sections` sub-object | key-space change |
| V2-09 | Manifest self-hash or detached signature | new required root field |

## 5. Decisions needed from the maintainer

Recorded decisions: retain the current `$id` until permanent hosting is
settled (D09); keep the Python distribution name `qvf-writer`; prepare an
IANA registration draft, with submission reserved to the maintainer.
These decisions do not authorize a domain migration or imply IANA registration.

Still open:

1. Permanent schema hosting (D09), while preserving the current `$id`.
2. QTAIM position unit (D05).
3. Whether `Section` becomes closed per branch once peer fields are named (D08).

The `wavefunction.bloch` contract (P15) is already provisional in § 5.3.1.
Second-producer acceptance remains implementation/validation work, rather
than a decision about whether to specify its current shape.
