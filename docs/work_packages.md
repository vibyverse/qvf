# Work packages

Implementation briefs for the roadmap in `ROADMAP.md`. Each brief is meant to
be handed to one implementation chat as-is. Ids are stable; D-nn / P-nn / V2-nn
refer to `format_review.md`; "design § n" refers to `library_design.md`.

## Status (2026-10-22)

Done and on `main`: WP-01 to WP-05, WP-07, WP-10 to WP-15, WP-20 to WP-30,
WP-40 to WP-44, WP-46, WP-47 and WP-60. Work continues with WP-45 / WP-48 / WP-49.
Toolkit 0.5.0 includes the earlier released slices. Toolkit 0.6.0 adds
XSF/PDB conversion, the local Conan recipe and drawing-preserving pending-job
preparation. See [the release boundary](../ROADMAP.md#release-boundary-and-current-main).
WP-42 supplies every canonical typed kind, Writer and Editor; WP-43 supplies
the Reader and `get<K>` protocol. Writer shares the compatibility backend
rather than duplicating ZIP/update machinery.
Deviations from the briefs that later packages should know: WP-41 emits
`python/qvf/_kinds.py` and `python/qvf/model.py` in the WP-46 package;
`qvf_kinds` and `qvf_model` remain compatibility imports. WP-40's
`qvf::Limits` is an alias of
`qvf::container::Limits`; both writers now write `manifest.json` last; the
C++ validator and the JavaScript validator run the generated schema checker,
so tests that isolate one semantic rule validate with the check disabled.
WP-06 retains the existing schema `$id` pending hosting; WP-61 needs a Rust toolchain.

## How to run a work package

Paste this into a fresh chat, replacing the id:

> Implement **WP-nn** from `docs/work_packages.md` in this repository. Read
> `AGENTS.md`, the brief, the Definition of Done below it, and the files the
> brief names before changing anything. Keep the change scoped to the brief.

**Suggested model.** *Opus* for anything that touches normative semantics,
numerical conventions, schema architecture, ZIP internals, or public API
shape. *Sonnet* for mechanical, well-specified work: docs, CI, codegen
application, bindings that follow a fixed API, corpus additions.

**Sizes.** S = under a day, M = one to three days, L = a week or more of
focused work. Split an L package into the sub-steps its brief lists if the
chat's context runs long.

## Definition of Done (applies to every package)

1. The change is complete across every artefact it touches: spec, schema,
   `registry.json`, Python writer and reader/validator, C++ writer and reader,
   C ABI, conformance corpus (archive + `expected.json`), docs, `CHANGELOG.md`.
   A format change without all of these is not done.
2. A contract change gets a dated decision record in `GOVERNANCE.md` (same
   style as the existing ones) and, for new kinds, a `promotions` entry in
   `registry.json`. Additive changes never bump `qvf_version`.
3. Tests: `python -m pytest python/tests`, `cmake --build && ctest`,
   `python conformance/run_conformance.py`, `python scripts/amalgamate.py
   --check` all pass. New behaviour has a test; a fixed defect has a
   regression test.
4. `AGENTS.md` public-safety rules hold; the staged-content privacy hook
   passes; no private paths, hosts, or addresses anywhere, including fixtures.
5. The commit message names the WP id and the D/P ids it resolves.

---

## Phase 0 — hygiene and defects

### WP-01 · CHANGELOG, legacy v2 schema, dangling references — S · Sonnet · D10
**Goal.** Make every artefact the governance documents mention exist.
**Scope.** Add `CHANGELOG.md` (Keep-a-Changelog; seed `Unreleased` and a
`0.1.0` section from the `GOVERNANCE.md` decision log). Generate
`spec/legacy/qvf_manifest_v2.schema.json` from the v1 schema with
`qvf_version` `const: 2` and the v2 `$id`, with a header comment that it is
frozen. Fix the wording in `GOVERNANCE.md` / `registry.json` that points at
these files. Add a test that every file path mentioned in `GOVERNANCE.md`,
`registry.json`, and `README.md` exists.
**Acceptance.** Test passes; `registry.json` `version_decisions[0]` names the
legacy path.

### WP-02 · Schema description scrub and section-number sync — S · Sonnet · D01
**Scope.** Rewrite every schema `description` that cites
`design_qvf_format.md`, `vibeqc.*`, or implementation details (symlinks) into
self-contained text. Replace `§ n.n` citations with the spec's current numbers
*and* the heading title. Add `python/tests/test_schema_descriptions.py`: no
forbidden tokens; every `§ n.n` in the schema matches a heading in
`spec/qvf-format-spec.md`.
**Acceptance.** Test passes; corpus unchanged (descriptions are not contract).

### WP-03 · Portable pi, GitHub Actions matrix, sanitizers — M · Sonnet · D14
**Scope.** Replace `M_PI` with a `constexpr double` in `cpp/src/qvf.cpp`
(and anywhere else; grep). Add `.github/workflows/ci.yml` mirroring
`.gitlab-ci.yml` jobs with a matrix: ubuntu (gcc, clang), macos (clang),
windows (MSVC) for C++; Python 3.9, 3.11, 3.14 for the package; one Linux job
with `-fsanitize=address,undefined`; the amalgamation drift check; the
conformance run. Keep `.gitlab-ci.yml` in sync (same steps).
**Acceptance.** All jobs green on a PR; MSVC compiles with `/W4 /WX`.

### WP-04 · Writer/schema mismatches — M · Opus · D02, D03, D08
**Scope.** (a) Schema: add `phonon_bands.eigenvectors` (float64
`[n_q, n_modes, n_atoms, 3]` or `[n_q, n_modes, n_atoms, 3, 2]` complex) and
`phonon_dos.projected` (float64 `[n_channels, n_points]`); document in spec
§ 5.6. (b) One `atom_properties` role list in schema, spec § 5 table, § 5.7,
Python writer, C++ writer, C ABI (five roles today; P10 extends later).
(c) Name the DOS peer fields (`n_spin`, `smearing`, `smearing_type`,
`fermi_energy_ev`, `n_electrons`, `channels`) in `SectionDOSTotal` /
`SectionDOSProjected`. (d) Regression test: a full-writer archive using every
optional argument validates with `jsonschema`.
**Acceptance.** The reproduction in `format_review.md` D03 validates; corpus
gains `phonons.qvf` and `atom_properties_all.qvf`.

### WP-05 · Validator hardening and negative corpus — M · Opus · D15, D17
**Scope.** Spec § 6 and both validators: duplicate ZIP entry names (error),
backslash or drive-letter member paths (error), JSON size/depth limits
(error), total expansion and compression-ratio caps (error), orphan ZIP
entries (warning), reserved `viewer_defaults` ids `auto_open`/`bookmarks`
used as section ids (error). Introduce stable finding codes
(`E_…`/`W_…`, design § 5) in Python; C++ `validate()` returns the same codes.
Add `conformance/negative/` with hand-built must-fail archives and a runner
mode that asserts the expected code.
**Acceptance.** Every negative archive fails with its expected code in both
validators.

### WP-06 · Schema `$id` domain and hosting — S · Sonnet · D09
**Status.** The 2026-10-15 decision retains the existing `$id`; a domain
migration is not pending authorization. Hosting remains open.
**Remaining scope.** Document and provide hosting for the retained identifier
and a check that fetched content matches the committed schema (allowed to
skip offline). Do not rename identifiers or package URLs merely to implement
the superseded domain-change proposal.

### WP-07 · Privacy: deprecate `hostname`, scrub semantics in the spec — S · Opus · D11
**Scope.** Spec § 3.2: `hostname` deprecated (producers SHOULD NOT emit,
consumers MUST ignore). Writers drop it unless `allow_hostname=True`. Spec
§ 5.8 gains a "sharing and sanitisation" paragraph. The `scrub` tool itself is
WP-47.

## Phase 1 — contract completeness

### WP-10 · Payload schemas, batch A (geometry and wavefunction) — L · Opus · P01, D07
**Scope.** `$defs` for volume `grid`, `basis` (shells, `structure_ref`,
`pure`, `n_ao`, `ecp` placeholder per P09), `mo_metadata` (restricted and
unrestricted shapes, `k_point`, `orbital_kind`, `spin`), trajectory
`metadata`, `reaction.path` `metadata` + `volume_grid`, `waypoints`,
`scan.surface` `metadata`, `vibrations` `metadata`. Validator loads each by
kind/role (generalise `_load_payload_validator`). Spec § 5 subsections gain a
"Payload" table per member. Corpus: every existing archive still passes;
add `trajectory.qvf`, `vibrations.qvf`, `scan_surface.qvf`,
`reaction_path.qvf`.
**Sub-steps.** 1 grid+basis+mo_metadata · 2 trajectory family · 3 vibrations.

### WP-11 · Payload schemas, batch B (spectra and analysis) — L · Opus · P01, P13, D05, D18
**Scope.** `SpectrumPayload` (common 1-D shape with `x_units`, `y_units`,
`peaks`, `broadening`, kind blocks), typed `spectra.nmr` and `spectra.epr`
(permissive root), `bonds`, `bond_orders`, `scf_history.iterations`,
`topology.qtaim.critical_points` (with the decided unit and `units` object),
`structure.symmetry.data` (permissive: `spacegroup_number`, `hall`,
`international`, `operations[]`, `wyckoffs[]`). Spec § 5.4 / § 5.7 rewritten
around the payload tables. Corpus: `spectra_all.qvf`, `nmr.qvf`,
`qtaim.qvf`, `symmetry.qvf`.

### WP-12 · Payload schemas, batch C (periodic) + spin axis — L · Opus · P01, P14, D04
**Scope.** `kpath` (`kpoints`, `labels`, `segments`, `fermi`, `spinor`,
`channels`), DOS `meta`/`channels` typed, `fermi_surface.mesh` (+`shift`,
`symmetry_expanded`), `qpath`, phonon `meta`, `equation_of_state.fit`
(`model` enum, parameters), run-record `files`. Spin axis: accept
`[n_spin, n_k, n_bands, n_channels]` and `[n_spin, n_channels, n_points]`;
lower ranks mean `n_spin = 1`; validator cross-checks against `eigenvalues`.
Optional `bands.occupations`. Writers gain the new shapes. Corpus:
`bands_spin.qvf`, `dos.qvf`, `fermi_surface.qvf`, `eos.qvf`.

### WP-13 · Units vocabulary — M · Opus · P02, D12
**Scope.** Spec § 4.3 vocabulary table; `ValueWithUnits.units` enum; optional
`units` objects on the payloads named in P02; `qvf.units` Python module and
`qvf::Unit` C++ enum (generated, WP-41 makes it so; hand-written first);
validator warns on unknown unit strings and errors on a declared unit that
contradicts a fixed-unit field. Writers accept `units=` and convert nothing
(declare only).

### WP-14 · Conformance corpus v2 and C++ runner — M · Sonnet · (after WP-10..12)
**Scope.** One archive per canonical kind (coverage test asserts it), edge
cases: unrestricted wavefunction, Cartesian `pure: false` d shell with the
Appendix A.1 overlap values in `expected.json`, 1-D periodic wire, variable
cell reaction path, `qvf_version: 2` legacy archive (read-only), a ZIP64
synthetic archive (small members with ZIP64 structures forced), a
data-descriptor archive. Add `cpp/examples/run_conformance.cpp` implementing
the check format and run it in CI over the corpus.

### WP-15 · Lattice conventions documented + helper — S · Sonnet · D06
**Scope.** Spec § 4.3 and § 5.5: an explicit side-by-side of the two lattice
conventions with a worked numeric example; `Lattice` helper in both
implementations (design § 4) with tests converting between them; note that
P12's `trajectory.lattice` follows the structure convention.

## Phase 2 — format expansion (all additive)

### WP-20 · Root `metadata` + generic section refs — M · Opus · P03, P04
**Scope.** Schema root `metadata` block; `structure_ref`, `run_ref`, `group`
on `Section`; validator resolves refs and checks `structure_ref` kind; writers
set `created_utc` + `uuid` by default (opt-out) and accept the refs on every
`add_*`; spec § 3.4 (new) and § 4.1; corpus `multi_structure.qvf`.

### WP-21 · `energies`, `forces`, `hessian` — M · Opus · P05, P06, P07
**Scope.** Root `energies` block with the consistency check against
`provenance.scf_energy`; kinds `forces` (gradient + optional stress + meta)
and `hessian` (matrix + meta); registry entries; writers, readers, C ABI
entries (`qvf_writer_add_forces`, `_add_hessian`); spec § 5.10–5.11; corpus
`forces_hessian.qvf`.

### WP-22 · `excited_states` + transition densities — M · Opus · P08, parts of P09/P11
**Scope.** Kind `excited_states` (`states` + `meta` payload schemas);
`mo_metadata.orbital_kind: nto` + `state_ref`; `volume.transition_density`
kind; `SpectrumPayload.peaks[].assignment.state_index` (from WP-11). Writers:
`add_excited_states(states, meta)`; corpus `tddft.qvf` with a consistent
UV-Vis spectrum derived from the state table.

### WP-23 · Wavefunction extensions — M · Opus · P09
**Scope.** `basis.ecp`, `n_electrons`, `symmetry_labels`, `spin` enum growth,
`complex` flag + two-component description in spec § 5.3, optional `overlap`
and `density_matrix*` members. Validator: when `overlap` is present, compute
the Appendix A.1 self-overlap of each shell's axial component from the
basis and compare (tolerance) — this turns the normalization trap into a
checkable error. Corpus `wavefunction_ecp.qvf`, `wavefunction_uhf.qvf`.

### WP-24 · `atom_properties` expansion — S · Sonnet · P10
**Scope.** New scalar roles, `AtomVectorBinary` roles, `custom.<name>` with
required `meta`; writers accept `**roles`; C ABI `add_atom_properties(role,
…)`; spec § 5.7; corpus `atom_properties_all.qvf` extended.

### WP-25 · Volume additions and periodic grid semantics — M · Opus · P11
**Scope.** Kinds `volume.vector`, `volume.complex`, `volume.series` (+`axis`),
`volume.laplacian`; `grid.units`, `grid.periodic`, `grid.includes_endpoint`,
`grid.structure_ref`; validator cell-span check for periodic grids; writers
(`add_volume` accepts 4-D for vector/complex/series); spec § 5.2 rewritten with
the endpoint rule and a cube-vs-CHGCAR example; corpus `volume_periodic.qvf`,
`volume_vector.qvf`.

### WP-26 · Trajectory extensions — M · Opus · P12
**Scope.** Optional `lattice` (rows, Å), `velocities`, `forces`, `times`,
binary `energies`, `cell_volumes`, `temperatures`; float32 coords with
`metadata.precision`; `metadata.kind` enum and `pbc`; validator shape
consistency across members; writers; C ABI `add_trajectory` full signature;
spec § 5.5; corpus `md_periodic.qvf`.

### WP-27 · Provenance and job-spec tightening — S · Sonnet · P18
**Scope.** `provenance.software[]`, `git_commit`, `n_threads`, `n_mpi_ranks`,
`started_utc`, `finished_utc`; `job.spec` `structure_ref`, `solvation`,
`dispersion`, `relativistic`, `scf`, `geometry`, `resources`; writers; spec
§ 3.2 / § 5.9; corpus `job_spec.qvf` regenerated.

### WP-28 · Container rules — M · Opus · P19, P20
**Scope.** Optional member `size` (writers emit, validators check), duplicate
entries forbidden (from WP-05), orphan rule, in-place update semantics in a
new spec § 2.3, root `byte_order` const; spec § 6 updated; corpus
`updated_in_place.qvf` (an archive with one orphan member, validates with a
warning).

### WP-29 · Periodic extras — M · Sonnet · P14, P17
**Scope.** `fermi_surface.mesh.shift` / `symmetry_expanded` / explicit
`kpoints`; `properties.electric` kind; `properties.elastic` and
`properties.dielectric` kinds (JSON payloads with fixed units); magnetic
moment vector role (from WP-24); writers; spec; corpus.

### WP-30 · `wavefunction.bloch` provisional — M · Opus · P15 (needs maintainer decision 5)
**Scope.** Promote the `x_vibeqc.bloch_wavefunction` layout to a provisional
canonical kind sharing `basis.json`; registry `promotions` entry marked
provisional; writers; spec § 5.3 companion section; corpus `bloch.qvf` with a
two-k-point toy system.

## Phase 3 — libqvf

### WP-40 · Container core: ZIP64, streaming, append — L · Opus · D13, D19 · design § 3
**Scope.** C++ `container::ZipWriter` (sequential, data descriptors on
non-seekable sinks, local-header patch on seekable, ZIP64 when any size or
offset ≥ 0xFFFFFFFF or on request), `ZipReader` (ZIP64 EOCD + locator, extra
field 0x0001, duplicate-name rejection, `Limits`), streaming inflate/deflate,
incremental SHA-256 and CRC. Append mode: open existing, truncate at the
central directory, add entries, rewrite. Python `_container.py` with the same
semantics on `zipfile` (ZIP64 forced on request; append via `mode="a"` plus
manifest rewrite). Tests: round-trip with Python ↔ C++ on synthetic ZIP64 and
descriptor archives; fuzz target for the reader.
**Sub-steps.** 1 ZIP64 read · 2 ZIP64 write · 3 streaming sinks · 4 append.

### WP-41 · Typed model + schema-driven codegen — L · Opus · design § 4–5
**Scope.** `scripts/gen_kinds.py` reads the schema and emits
`cpp/include/qvf/kinds_generated.hpp` (KindInfo/RoleInfo tables, unit enum)
and `python/qvf/_kinds.py` (with a `qvf_kinds` compatibility import); the
table-driven payload checker for the JSON
Schema subset the QVF schema uses; `Manifest`/`Section`/`MemberSpec` model
types in both languages; drift test that regenerated output equals committed
output; parity test: generated checker ≡ `jsonschema` on corpus + negative
corpus.

### WP-42 · C++ `Writer` (streaming) and `Editor` — L · Opus · design § 7
**Scope.** `qvf::Writer` over `container::ZipWriter` with typed `add()`
overloads for every kind in `kinds/`, `begin_binary()` member sinks,
`TensorView` with on-the-fly byte-order/dtype conversion, manifest written
last; `qvf::Editor` (open, append, remove_section, set_run_status, commit,
repack); `QvfWriter` kept as a compatibility layer; amalgamation updated.
Tests: streaming and in-memory produce identical member bytes; editor
round-trip through the Python validator.

### WP-43 · C++ `Reader` typed accessors + validator — M · Opus · design § 5, § 7
**Scope.** `qvf::Reader` with `sections_of`, `json`, `tensor`,
`tensor_into` (streaming), `get<K>()`; `validate::run` with levels and the
shared finding codes; `read_qvf` example becomes `qvf-inspect`; conformance
C++ runner (WP-14) switched to the new API.

### WP-44 · Complete C ABI — L · Sonnet (after 42/43 fix the API) · D16 · design § 8
**Status.** Included in toolkit 0.5.0. ABI 1 keeps all
36 existing C symbols and adds all-kind writers, periodic structures, lazy
readers, staged validation, optional atomic editing and shared-library exports.
The [C API guide](c_api.md) documents descriptor inputs, ownership and errors.
**Scope.** Every entry listed in design § 8 for writer, reader, validator,
editor; `QVF_API` export macro and `QVF_BUILD_SHARED`; `qvf_abi_version()`,
`qvf_toolkit_version()`; existing 16 symbols unchanged; C example extended
(`write_all_kinds_c.c`, `read_qvf_c.c`); Python-validated in CI.
**Sub-steps.** 1 writer · 2 reader · 3 validator/editor.

### WP-45 · Fortran module — M · Sonnet · design § 8
**Status.** Source/API slice included in toolkit 0.5.0.
Generated raw interfaces cover all 161 C exports; checked writers cover
structure, scalar volume, trajectory and restricted GTO. Focused gfortran
checks verify descriptor sizes, refusal/ownership and Python-validated array
layouts. Installed CMake consumption is checked. The 0.5.0 Linux GNU 12
candidate passed all 12 CTests, including Fortran layout/ownership checks
with runtime bounds checking, alongside prior focused GNU 16 checks.
Focused development checks on Ubuntu 24.04 x86_64 pass with gfortran 11, 12,
13 and 14 and LLVM flang 18, 19 and 20 (warnings as errors, Release and
Debug, CTest, Python validation, installed consumer). They fixed two
wrapper defects: truncated or garbage `message` text on gfortran 11/12, and
an IEEE invalid trap on NaN GTO exponents under halting mode with gfortran.
The matrix is in [the Fortran guide](fortran.md#checked-compilers). macOS,
Windows and other vendor compilers remain unchecked. WP-45 stays in progress.
**Scope.** `fortran/qvf.f90` (`module qvf`, `bind(C)` interfaces + wrappers
handling column-major → row-major), CMake `QVF_FORTRAN` option, example
writing structure + grid + trajectory + wavefunction, CI job with gfortran
validating the output with the Python reader; docs section in
`library_guide.md`.

### WP-46 · Python package restructure — L · Opus · D20 · design § 8
**Status.** Included in toolkit 0.5.0; distribution remains
`qvf-writer`, with optional `numpy` and `validate` extras.
**Scope.** `python/qvf/` package per the design layout; dataclass model;
streaming `Writer` (members written as added), lazy `Reader` with typed
accessors and a numpy-free `memoryview` path, `Editor`; `validate` with
finding codes; schema as package data; shims `qvf_writer.py` /
`qvf_reader.py` with `DeprecationWarning`; all existing tests pass unchanged
through the shims; new tests for the package API. Distribution name per
maintainer decision 4.
**Sub-steps.** 1 package + model + shims · 2 reader · 3 writer/editor.

### WP-47 · CLI `qvf` — ✅ implemented · design § 8
**Scope.** Subcommands `validate` (levels, JSON output of findings), `info`
(tree of sections, members, sizes, kinds), `cat` (JSON member), `extract`,
`diff` (two archives, numeric tolerance, exit code), `repack`, `scrub`
(hostname, `run.record` path/user patterns, dry-run report, never silent),
`conformance` (run the corpus), `schema` (print). `qvf-validate` alias kept.
**Delivered.** All commands share the packaged library; output paths must be
new, scrub reports changes/skipped content and supports dry runs, and early
validation stages avoid later reads/checks. See [CLI guide](cli.md).

### WP-48 · Converters — L · Sonnet (after 46) · design § 8
**Status.** First bounded slice landed on main: scalar cube and geometry-only
XYZ/extXYZ import/export, including trajectories, explicit periodic endpoint
handling, validation and reported conversion losses. See [converter guide](converters.md).
Cube/XYZ are included in toolkit 0.5.0. The 0.6.0 static XSF slice adds
primitive geometries and one scalar grid with explicit field periodicity,
unit/layout conversion and preserved endpoint samples. The 0.6.0 PDB
geometry slice handles explicit element/coordinate columns, refuses
ensembles and alternate locations, and reports metadata and precision losses.
The Unreleased Molden slice supports real molecular s/p/sp/d orbitals with
explicit spin, AO normalization/order conversion and loss reports. Independent
PySCF import/export checks compare orbital values. Broader Molden dialects,
richer extXYZ and the other periodic formats remain open. PDB atom-site identity
(names, residues, chains, insertion codes, alternate locations, occupancy,
B-factors, `HETATM`) and model ensembles are on main (Unreleased, #27), and
so are `HELIX`/`SHEET` secondary-structure ranges (#28), `CONECT` bonds
(#29), `CRYST1`/`SCALE` cells with the space group (#30), `SSBOND`/`LINK`
bonds (#31) and `REMARK 350` biological assemblies with `MTRIX` accounting
(#32), and opt-in crystal symmetry expansion to the unit cell (#33); PDB NCS
copy generation remains open. Core crystal CIF import with unit-cell expansion
and P 1 CIF export (#34), VASP POSCAR/CONTCAR (#35) and CHGCAR (#36)
import/export, mmCIF import/export (#37) and opt-in PDB NCS copy generation
(#38) extXYZ energy/forces/stress (#39), XSF forces and several grids (#40) and
Molden f/g shells (#41) are on main (Unreleased).
**Scope.** `qvf.convert`: import/export Gaussian cube (↔ `volume.*`, with
the endpoint rule), XYZ/extXYZ (↔ `structure`/`trajectory`), Molden
(↔ `structure` + `wavefunction.gto`, with Appendix A normalization and
ordering handled), PDB (→ `structure` with biomolecule metadata), XSF
(↔ periodic structure + grid), POSCAR/CHGCAR (↔ periodic structure + grid,
`includes_endpoint: false`), CIF (→ structure). Each converter is
round-trip-tested on small fixtures and validated.
**Sub-steps.** 1 cube/xyz · 2 molden · 3 periodic formats · 4 pdb/cif.

### WP-49 · Packaging and release — M · Sonnet · ROADMAP "frictionless packaging"
**Status.** In progress. ABI-versioned shared builds, CMake library-only build
options and relocatable `qvf.pc` metadata are included in toolkit 0.5.0.
The 0.6.0 local Conan 2 recipe and installed C/C++ consumers are in
`packaging/conan`; the release/adoption checklist is in `CONTRIBUTING.md`.
Registry publishing, vcpkg and broader platform acceptance remain open.
**Scope.** PyPI release workflow (trusted publishing from GitHub Actions,
tag-driven), `qvf.pc`, `QVF_BUILD_SHARED` with SONAME = ABI version, vcpkg
port and conan recipe in `packaging/`, release checklist in
`CONTRIBUTING.md`, `scripts/make_archive.sh` updated for the new layout.

### WP-50 · Fuzzing and big-endian CI — M · Sonnet · design § 10
**Scope.** libFuzzer targets for `ZipReader`, inflate, `Json::parse`,
`Manifest::from_json`, `validate::run`; a GitHub Actions job running each
target for a bounded time; OSS-Fuzz-ready layout; one QEMU s390x job that
writes the corpus with both writers on a big-endian host and validates it.

## Phase 4 — ecosystem

### WP-60 · JavaScript/TypeScript reader — L · Sonnet · design § 8
**Scope.** `js/` package `@qvf/core`: ZIP (central directory, ZIP64) + inflate
via `fflate`, SHA-256 via WebCrypto, manifest model, typed accessors, validator
with the shared finding codes, conformance runner over the corpus in CI
(Node), browser smoke test, `README` with a 20-line viewer example.

### WP-61 · Rust crate — L · Sonnet · design § 8
**Scope.** `rust/` crate `qvf`: reader, writer, validator (finding codes),
typed kinds for structure/volume/wavefunction/trajectory/bands/spectra,
conformance runner as an integration test, `qvf-cli` feature. Certified by
the corpus; no C ABI dependency.

### WP-62 · PySCF adapter (second producer) — M · Opus · ROADMAP Tier 2
**Scope.** `python/qvf/adapters/pyscf.py` (optional dependency): molecule →
`structure`, SCF/DFT → `wavefunction.gto` with `coeffs_are_libint_normalized`
handled, `overlap` member emitted so the WP-23 check runs, populations →
`atom_properties`, TD-DFT → `excited_states` + `spectra.uvvis`, densities on a
grid → `volume.density`, PBC Γ-point → periodic structure + Γ wavefunction.
Tests gated on `pyscf` being importable; a notebook-style example in docs.
This is the producer that earns the first real vendor→canonical promotions.
**Status.** Tranche 1 (#23) landed the molecular SCF/DFT adapter,
`python/qvf/adapters/pyscf.py` with the `qvf-writer[pyscf]` extra: RHF, UHF,
ROHF, RKS, UKS and ROKS results become `structure`, `wavefunction.gto`
(spherical or Cartesian `basis` on `N_i`-normalized primitives, Appendix A.2
order, `overlap` and density members), Mulliken and Löwdin `atom_properties`,
`job.spec` and the `provenance`/`energies` blocks, with no format change.
`python/tests/test_adapter_pyscf.py` checks it against PySCF and an
Appendix A evaluator, and `scripts/pyscf_cross_producer.py` against vibe-qc
archives (energies within 3e-12 hartree, `max |ΔP|` about 1e-8). See the
[PySCF adapter guide](adapters/pyscf.md). Tranche 2 (#25) adds closed-shell
RHF/RKS TDA and TDHF/TDDFT `excited_states` + `spectra.uvvis`, with paired
singlet-TDA NTOs and source-amplitude preservation. Frozen MO indices map
to the full reference wavefunction. Tranche 3 (#26) adds singlet full-response
NTOs from PySCF's `sqrt(2)(X+Y)`, retaining non-unit squared singular-value
weights and reconstructing the length-gauge transition moments. Triplet NTOs and
unrestricted/open-shell response remain deferred. Grid output adds
`volume.density`, `volume.spin` and `volume.orbital` (`add_density_grid`,
`add_orbital_grid`). PBC Γ-point → periodic structure + Γ wavefunction
remains. Scalar ECP projectors
(`basis.ecp[].projectors`) are written from `mol._ecpbas`. Spin-orbit terms are omitted.

### WP-63 · vibe-qc links libqvf — L · Opus · maintainer note (low priority)
**Scope.** Replace vibe-qc's internal QVF writer/reader with the standalone
library (C++ via `find_package(qvf)` or FetchContent at a pinned tag; Python
via the `qvf` package). Keep vibe-qc's schema symlink pointing at the
installed schema. Acceptance is vibe-qc's own test suite plus the corpus
passing through vibe-view. Scheduled after WP-42/43/46 stabilise the API;
until then vibe-qc keeps its copy and the cross-repo drift test.

### WP-64 · Docs site and kind catalogue — M · Sonnet
**Scope.** mkdocs site from `docs/` + generated API reference (Doxygen,
mkdocstrings) + a kind catalogue page generated from the schema (one entry per
kind: members, dtypes, shapes, units, payload fields, example JSON) + updated
guides for the new entry points; built in CI, published from the public
mirror.
**Status.** Tranche 1 (#19) landed the kind catalogue: `docs/kinds.md`,
rendered by `scripts/gen_kind_docs.py` from the schema, `registry.json` and
the must-pass corpus, and covered by `gen_kinds.py --check`. Tranche 2 (#20)
landed `mkdocs.yml`: a strict local build over `docs/` plus the top-level
documents, with a pinned toolchain (`scripts/requirements-docs.txt`) and link
hooks (`scripts/mkdocs_hooks.py`). Tranche 3 (#21) landed the Python API
reference (`docs/api/python.md`, mkdocstrings). Tranche 4 (#22) landed the
C++ and C reference: `cpp/Doxyfile`, run by the mkdocs hook into
`api/cpp-reference/`, with `scripts/doxygen_filter.py` presenting the
headers' plain comments as documentation. Publishing remains; it is a
decision for the release owner.

### WP-65 · Julia binding — M · Sonnet (demand-driven)
**Scope.** `julia/QVF.jl` over the C ABI via `ccall`, reader + writer
subset, conformance runner as a test.

---

## Dependency summary

```
WP-01..07  independent (Phase 0)
WP-10/11/12 → WP-13 → WP-14 ; WP-15 independent
WP-20..30  depend on WP-10..12 (payload-schema mechanism) ; WP-22 after WP-11 ; WP-29 after WP-24
WP-40 → WP-42 → WP-44 ; WP-41 → WP-42, WP-43, WP-46
WP-46 → WP-47, WP-48 ; WP-44 → WP-45, WP-65
WP-49, WP-50 after WP-42..46
WP-60, WP-61 need only the corpus (after WP-14) ; WP-62 after WP-22/23/46 ; WP-63 after WP-42/43/46
```
