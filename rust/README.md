# qvf (Rust)

An independent Rust implementation of the QVF archive format (WP-61): a
streaming writer with section builders, an editor for in-place updates and
repacking, a reader with typed section models, a validator and a `qvf`
command-line tool, certified by the shared
conformance corpus in
[`../conformance`](../conformance). It does not use the C++ library or its C
ABI. Status: in development (`0.1.0-dev`, `publish = false`).

## Build and test

```sh
cd rust
cargo build
cargo test                      # unit, conformance, robustness, writer, builder, editor, typed-kind, regex and CLI tests
cargo test --test conformance -- --nocapture   # prints the pass/fail/skip counts
cargo clippy --all-targets -- -D warnings
cargo run --example validate -- ../conformance/corpus/structure_periodic.qvf
cargo run --bin qvf -- validate ../conformance/corpus/structure_periodic.qvf
```

Edition 2021, declared MSRV 1.74 (`rust-version` in `Cargo.toml`). The crate
is developed with current stable Rust, and `cargo +1.74 test --locked` passes. The release-candidate CI (`rust` job)
runs the tests on stable and on 1.74, and fmt and clippy on stable. Build output
stays in `rust/target/` (ignored).

## Dependencies

| Crate | Why |
|---|---|
| `flate2` (pure Rust `miniz_oxide` backend) | DEFLATE and CRC-32 |
| `serde_json` | `manifest.json` and the JSON members |
| `serde` (with `derive`) | the typed models of `src/kinds.rs` |
| `sha2` | member checksums (spec § 4.2) |

The ZIP container layer is our own, a reader ported from `js/src/zip.ts`
(about 400 lines) and a streaming writer (`src/zip_writer.rs`, about 300
lines), instead of the `zip` crate. The validator must see every central-directory
record, duplicate names included, and check sizes and compression ratios
before anything is decompressed (spec § 2.1, § 2.3, § 6.1). It reads ZIP64
(extra fields, the ZIP64 end record and locator) and data-descriptor
archives, STORE and DEFLATE, and never inflates past the declared size.
The writer needs no seekable sink: it writes data descriptors after each
entry and switches to ZIP64 fields and end records when a size, an offset or
the entry count needs them. The writer adds no dependency.

The JSON Schema check uses no JSON Schema crate. `scripts/gen_kinds.py`
compiles `spec/qvf_manifest.schema.json` into an instruction table, the same
one the C++ and TypeScript validators interpret, and writes it with the kind,
unit and finding-code tables to `src/generated.rs`. `src/schema_check.rs`
interprets it (a port of `js/src/schemaCheck.ts`, with the same linear-time
Pike VM for patterns). The same script writes the typed models to `src/kinds.rs`
(rendered by `scripts/gen_typed_kinds.py`, the mapping of the C++
`cpp/include/qvf/kinds/all.hpp`). Do not edit `src/generated.rs` or
`src/kinds.rs` by hand: regenerate with `python3 scripts/gen_kinds.py` and
check freshness with `python3 scripts/gen_kinds.py --check`.

## API

```rust
use qvf::{validate_file, Reader, ValidateOptions};

let report = validate_file("water.qvf", &ValidateOptions::default());
assert!(report.ok);

let reader = Reader::open("water.qvf")?;
let structure = reader.read_json("mol", "structure")?;       // serde_json::Value
let coeffs = reader.read_tensor("wf", "mo_coefficients")?;   // dtype, shape, typed data
let first = coeffs.at(&[0, 0]);                              // Option<f64>
let log = reader.read_bytes("run", "log")?;                  // raw bytes, sha256 verified
```

* `Reader`: opens an archive (the last `manifest.json` of an archive updated
  in place, spec § 2.3), refuses duplicated entry names, applies the
  resource caps of spec § 6.1, and reads members as raw bytes, JSON or
  little-endian tensors of the manifest `dtype` and `shape`, verifying a
  declared `sha256` and `size`.
* `Reader::get::<K>(id)`: typed access, see *Typed kinds* below.
* `Writer`: see *Writing* below.
* `Editor`: see *Editing* below. `qvf::operations` holds `diff` and
  `scrub`, `qvf::conformance` the conformance runner, and `qvf::regex` the
  pattern engine of `scrub`.
* `validate_bytes` / `validate_file`: every finding of spec § 6 with its
  stable Appendix B code (`qvf::Code`), its severity and its location
  (`archive`, `entry:<name>`, `manifest.json` or a JSON pointer). Levels
  (`Level::Container` up to `Level::All`) and limits work as in the other
  validators. `schema_check: false` skips the JSON Schema check and reports
  `W_SCHEMA_SKIPPED`.

## Typed kinds

```rust
use qvf::kinds::{Structure, WavefunctionGTO, VolumeDensity, ExcitedStates, SpectraIR, Trajectory, AnyKind};

let reader = qvf::Reader::open("tddft.qvf")?;
let mol: Structure = reader.get("structure")?;            // fields, structure, bonds
let x = mol.structure.atoms[0].position;                   // [f64; 3], Å
let wf: WavefunctionGTO = reader.get("wf")?;               // basis, mo_metadata, tensors
let c = wf.mo_coefficients.as_ref().or(wf.mo_coefficients_alpha.as_ref());
let rho: VolumeDensity = reader.get("rho")?;               // grid (VolumeGridPayload), data (Tensor)
let any = reader.get_volume("rho")?;                       // any volume.* kind
let section = reader.get_any("excited_states")?;           // Option<AnyKind>, None for a vendor kind
```

`src/kinds.rs` is generated from the schema `$defs` by
`scripts/gen_kinds.py`: one serde model per payload definition and per
nested object, and one struct per canonical kind (all 56, named as in the
C++ `qvf::kinds`, such as `Structure`, `WavefunctionGTO`, `VolumeDensity`,
`ExcitedStates`, `SpectraIR`, `Trajectory`) that implements
`qvf::TypedKind`. A kind struct holds `fields` (the section fields of
`$defs/Section` and of the kind, typed), each JSON member as its payload
model, each binary member as a `Tensor`, each opaque member as `OpaqueData`
and pattern roles as a map. Required members are plain fields, optional ones
`Option`. In every model the keys the schema leaves open (`x_*` and the
like) are kept in a flattened `extra` map, an optional nullable field is an
`Option<Option<T>>` (absent and `null` differ), and fixed-length arrays are
Rust arrays.

`Reader::get::<K>` fails with `E_REF_UNRESOLVED` for an unknown id,
`E_KIND_UNKNOWN` for a section of another kind, `E_SCHEMA` when the section
fields do not fit the model and `E_PAYLOAD_SCHEMA` when a JSON member does
not. The models check shapes and types, not value ranges, patterns or
cross-references: that stays with the validator. `TypedKind::to_writer` (and
`AnyKind::to_writer`) writes a section back through a `Writer` at the
writer's default member paths.

`tests/typed.rs` decodes every section of the positive and legacy corpora
(81 archives, 237 typed sections) with the model of its kind, writes each
archive back through the models, and requires the rewrite to validate
without errors, to hold the same section fields, JSON payloads, tensors and
opaque bytes (numbers compared by value, as a typed `f64` writes `2` back as
`2.0`), and to decode to equal models.

## Command line

The `qvf` binary target (no extra dependencies, no feature gate) ports
`validate`, `info`, `cat`, `extract`, `diff`, `repack`, `scrub`,
`conformance` and `schema` of the Python `qvf` CLI:

```sh
qvf validate ARCHIVE [--level container|manifest|payloads|semantics|references|all] [--json]
qvf info ARCHIVE [--json]
qvf cat ARCHIVE [SECTION MEMBER]          # the manifest, or one JSON member
qvf extract ARCHIVE SECTION MEMBER -o OUTPUT
qvf diff LEFT RIGHT [--atol X] [--rtol X] [--ignore-metadata] [--max-differences N] [--json]
qvf repack ARCHIVE -o OUTPUT [--compression preserve|store|deflate]
qvf scrub ARCHIVE (-o OUTPUT | --dry-run) [--path-pattern REGEX] [--user-pattern REGEX] [--replacement TEXT]
qvf conformance CORPUS [--negative DIR | --no-negative] [--legacy DIR | --no-legacy]
qvf schema
```

Exit status: 0 success (or equal), 1 validation failed (or different, or a
conformance failure), 2 usage or I/O error, as in Python. `validate --json`
prints the Python report (`ok`, `errors`, `warnings`, `n_sections`,
`findings` with `severity`, `code`, `where`, `message`, and
`validation_context` with `reader: "qvf-rust"`), the human form the same
verdict, context and finding lines. `info`, `cat`, `extract`, `diff`,
`repack` and `scrub` open an archive with the checks the Python `Reader`
makes (legal entry paths, caps, `qvf_version`). `extract`, `repack` and
`scrub` verify each member's size, CRC-32 and sha256 while they copy it and
never overwrite an existing output (the output is staged beside it and
linked into place when complete).

* `diff` compares the logical manifests (member `path`, `size` and `sha256`
  left out, `metadata` too with `--ignore-metadata`) and the decoded
  payloads: JSON members by value, float tensors element by element
  (`math.isclose` with `--atol` and `--rtol`), other members byte by byte in
  1 MiB chunks. Integers compare exactly and NaN differs. Members are walked
  in manifest file order, as in Python, so the listed (and truncated)
  differences are the same.
* `repack` is `Editor::repack` (below).
* `scrub` removes `provenance.hostname` and the authors' `email` and
  `affiliation`, redacts the matches of the patterns in each `run.record`
  command and in its UTF-8 `input`, `log`, `attachment.*` and JSON `files`
  members, and writes a repacked archive with a new `metadata.uuid` (the old
  one in `derived_from`) and `created_utc`. Patterns are the commonly used
  subset of Python's `re` syntax (`qvf::regex`: literals, escapes, classes,
  `\d \w \s`, anchors and `\b`, groups, lookaheads, alternation, greedy and
  lazy quantifiers, leading `(?i)`, `(?m)`, `(?s)`), matched with Python's
  semantics by a linear-time Pike VM. Backreferences, lookbehinds, scoped
  flags and possessive quantifiers are refused with exit status 2.
* `conformance` runs `qvf::conformance` and prints the Python runner's
  `PASS`/`FAIL` lines and summaries.
* `schema` prints `rust/qvf_manifest.schema.json`, a byte-identical copy of
  `spec/qvf_manifest.schema.json` embedded at build time (a test checks the
  copy and its sha256 against `generated::SCHEMA_SHA256`; refresh it with
  `cp ../spec/qvf_manifest.schema.json .` when the schema changes).

`tests/cli.rs` runs the binary and the Python CLI (`qvf.cli.main`, in one
interpreter, with `PYTHONPATH=../python` or `$PYTHON`) on all 201 archives
of the positive, negative and legacy corpora. `validate --json` at each of
the six levels must give the same exit code, verdict, set of finding codes,
set of (code, severity, location) findings and validation context, the human
`validate` the same verdict line and finding lines, and `info --json` and
`cat` the same exit code and JSON (1809 commands). `diff` (each archive
against itself and against the next one, in JSON and human form), `repack`
(preserving and as STORE), `scrub` (dry run and written, with and without
patterns) and `schema` must give the same exit code and output, and the
written archives the same entries, manifest (but for the renewed `uuid` and
`created_utc`) and findings (1810 commands). A perturbed archive checks the
`diff` tolerances and truncation, and seven `conformance` invocations the
runner. Everything is skipped with a `SKIP` line when Python or the package
is unavailable. Where the Python CLI raises an uncaught exception (`diff`
and `scrub` on manifests whose `sections` hold a non-object, 9 commands)
there is nothing to compare, and where its repack drops every member of
such a manifest (2 commands) the Rust repack keeps them.
`tests/regex.rs` compares 1072 pattern and text cases with `re.subn`.

## Writing

```rust
use qvf::{Atom, Basis, Orbitals, Shell, Source, SpinOrbitals, Structure, Volume, WavefunctionGto, Writer};

let mut w = Writer::create("h2.qvf", Source::new("my-code", "1.0").with_calculation("h2-rhf"))?;
w.set_root("provenance", serde_json::json!({"method": "rhf", "basis": "sto-3g"}))?;
w.add_structure(&Structure::molecule("mol", vec![Atom::new("H", [0.0, 0.0, 0.0]), Atom::new("H", [0.0, 0.0, 0.74])]))?;
let shell = |c| Shell::new(c, 0, vec![3.42525091, 0.62391373, 0.1688554], vec![0.15432897, 0.53532814, 0.44463454]);
let orbitals = Orbitals::new(mo_coefficients, energies, vec![2.0, 0.0]);   // [n_mo, n_ao], C order
let mut wf = WavefunctionGto::new("wf", Basis::new("mol", true, vec![shell(0), shell(1)]), SpinOrbitals::Restricted(orbitals));
wf.overlap = Some(overlap);                                                 // lets validators check Appendix A.1
w.add_wavefunction_gto(&wf)?;
let grid = Volume::new("rho", "volume.density", origin, voxel_vectors, [nx, ny, nz]);
let mut sink = w.begin_volume(&grid, qvf::reader::DType::Float32)?;        // stream slabs of any size
for slab in slabs { sink.write(&slab)?; }
sink.finish()?;
w.finish()?;                                                                // manifest.json, STORED, last
```

* `Writer::new(sink, source)` writes to any `std::io::Write`, `Writer::create`
  to a file. `WriterOptions` choose DEFLATE (default) or STORE for members,
  forced ZIP64 and the `metadata.created_utc` stamp. `set_root` sets the
  other root blocks (`metadata`, `provenance`, `energies`, …).
* Members: `add_json`, `add_tensor` (any manifest dtype, written
  little-endian with `dtype` and `shape`), `begin_tensor` (a `TensorSink` to
  stream slabs, checked against the shape), `add_opaque` (with a media type,
  images STORED). Every member records `sha256` and `size`. `add_section`
  adds a section of any canonical or vendor kind for them.
* Builders (`qvf::sections`): `Structure` (molecules and periodic cells,
  bonds, charge and multiplicity), `WavefunctionGto` (restricted, restricted
  open-shell and unrestricted real coefficients, `overlap`, density
  matrices, NTO links), `Volume` (every `volume.*` kind, including vector,
  complex and series data, periodic grids, difference operands and
  transition-density states), `ExcitedStates` (states, contributions,
  `nto_ref`, `density_ref`, `meta`), and in `src/builders.rs` (re-exported
  by `qvf::sections`):
  * `Trajectory`: `trajectory` with `symbols` and `coords`
    `[n_frames, n_atoms, 3]` (float64, or float32 with
    `metadata.precision`), a fixed or per-frame cell (`TrajectoryCell`,
    rows in Å) with `pbc`, `velocities`, `gradient`, `times`, `energies`
    (optionally also as `metadata.energies`), `cell_volumes`,
    `temperatures`, `kind`, `ensemble`, `timestep_fs` and a `units`
    declaration.
  * `Bands` with its `Kpath` (k-points, labels, segments, Fermi level,
    spinor flag, channels), rank-4 or rank-3 `projections` and
    `occupations`.
  * `DosTotal`, `DosProjected` and `DosBonding` (`dos.coop`, `dos.cohp`),
    each with or without the spin axis, with `DosInfo` (smearing, Fermi
    level, electrons) and `DosChannel` and `DosPair` labels.
  * `Spectrum`: the six one-dimensional `spectra.*` kinds (frequencies,
    intensities, the `x` alias, `x_units`, `y_units`, `peaks` with mode and
    state assignments, `broadening`, `temperature_k`, the `raman` and
    `uvvis` blocks, `states_ref`), and `Writer::add_spectrum_payload` for any
    `spectra.*` kind, `spectra.nmr` and `spectra.epr` included.

  The new builders check every JSON payload against its schema definition
  (`schema_check::check_definition`) before anything is written, and refuse
  unknown unit names, a `states_ref` to a section not added yet and peak
  state indices that section does not have.
* Refused before anything is written: illegal or duplicate ids and paths,
  unknown kinds, impossible charge and multiplicity, shapes and lengths that
  do not match, references to sections not added yet, shell centres outside
  the structure, contributions outside the reference orbitals. A tensor sink
  that gets too many or too few bytes, or an I/O error, poisons the writer so
  that `finish` cannot write a manifest naming incomplete data.

`tests/builders.rs` writes an archive with every new builder (22 sections,
periodic NPT and slab trajectories among them) that has no findings from
either validator and reads back through the reader and the typed models,
and checks that misuse is refused before anything is written.
`tests/writer.rs` writes archives with every other builder, as DEFLATE, STORE and
forced ZIP64. Each must have no findings from the Rust validator and from the
Python reference validator (`python3 -m qvf validate --json`, run with
`PYTHONPATH=../python`, or `$PYTHON`), and must read back unchanged through
the Rust reader. The Python check is skipped with a `SKIP` line when Python
or the `qvf` package is unavailable, and accepts only `W_SCHEMA_SKIPPED`
when `jsonschema` is missing.

## Editing

```rust
use qvf::{Editor, EditorOptions, Spectrum};

let mut ed = Editor::open("water.qvf")?;                      // or open_with(path, EditorOptions { atomic: true, .. })
ed.writer().add_spectrum(&Spectrum::new("spectra.ir", freqs, intensities))?;   // any builder or member method
ed.replace_json("mol", "structure", &payload)?;              // written as mol/structure.2.json
ed.replace_tensor("rho", "data", &values, &shape)?;
ed.remove_section("old")?;
ed.set_run_status("converged", None)?;
ed.commit()?;
Editor::repack("water.qvf", "water-clean.qvf", None)?;       // None keeps each entry's method
```

`Editor` follows spec § 2.3 and the Python and C++ editors: it keeps every
entry's bytes where they are, appends the new members where the old central
directory began, and writes a complete new `manifest.json` (STORED, last)
and a central directory that lists the old records unchanged, then the new
ones. `source`, the root blocks and the sections are kept; builders can
reference the archive's sections (structure atom counts, orbital rows and
excited-state indices are read at open). A new member never reuses an entry
name: a path the archive holds becomes `stem.N.suffix`, the first free N
from 2 (`replace_*` derive it from the old path). `metadata.modified_utc` is
stamped unless the caller set it, `created_utc` and `uuid` are kept,
`job.spec` sections are kept verbatim, and an archive of the withdrawn
`qvf_version` 2 is written back as 1 with the v1 `schema_uri`. The editor
refuses what the `Reader` refuses and an archive whose entries do not all
lie before its central directory. In place, nothing changes before the
first new entry is written, and `abort` or dropping the editor restores the
old central directory. `atomic` edits a sibling copy and renames it over
the archive on commit. `Editor::repack` streams every referenced member
(verified while copied), then the manifest, into a staged sibling file and
renames it over the destination, which may be the source.

`tests/editor.rs` edits `tddft.qvf` in place and atomically (append with
references to the archive's sections, replace a JSON and a tensor member,
remove a section) and requires: the Rust and Python validators to report
the same findings, all of them `W_ORPHAN_MEMBER` (3, then 6); untouched
members to keep their spec and bytes; the orphans to count towards the
total-size cap (the edited archive exceeds a cap that its repack meets);
the repack to hold exactly the referenced entries and the manifest, with no
findings from either validator, and to equal the Python `Editor.repack` of
the same archive. Aborted and dropped edits leave the archive byte for byte
as it was, `job.spec` sections cannot be removed, and a version 2 archive
is written back as version 1.

## Supported checks

The validator is a port of the JavaScript validator (`js/src/validate.ts`)
and runs all of its checks in the same order, with the same codes:

* container: ZIP readability, duplicate entries, illegal paths, member and
  total size caps, compression ratio, orphan entries.
* manifest: size and depth caps, JSON, `qvf_version` (2 read as 1 with
  `W_VERSION_WITHDRAWN`), `source`, sections, ids, reserved ids, kinds,
  critical vendor sections, the JSON Schema (`E_SCHEMA`).
* members: paths, existence, declared `size`, sha256, JSON caps and syntax,
  payload definitions (`E_PAYLOAD_SCHEMA`), dtype and byte length, UTF-8 of
  `run.record` logs, `preview` signatures.
* semantics: structure periodicity, chemical identity, bands, trajectories
  and reaction paths, forces and hessians, atom properties, excited states
  and NTOs, `wavefunction.gto` including the Appendix A.1 overlap diagonal
  (`E_NORMALIZATION`), `wavefunction.bloch`, volumes and series, periodic
  grid spans, property kinds, 1-D spectra, reaction schemes and reaction
  energies, declared units, root consistency.
* references: every cross-reference the other validators resolve.

Current results (`cargo test --test conformance -- --nocapture`):

| Corpus | Passed | Failed | Skipped |
|---|---|---|---|
| `conformance/corpus` | 80 | 0 | 0 |
| `conformance/negative` | 120 | 0 | 0 |
| `conformance/legacy` | 1 | 0 | 0 |

On every archive, at every level, the set of (code, severity, location)
findings and the verdict equal those of the Python reference validator (with
`jsonschema` installed). `tests/cli.rs` checks this through the two CLIs.
The levels are staged as in the Python and C++ validators: `container`
reports a missing `manifest.json`, and `semantics` runs the cross-section
checks too and leaves out only `E_REF_UNRESOLVED`. The JavaScript validator
this one was ported from stages them differently (see `Level`).

## Known gaps

* Complex and `generalized` wavefunctions, `wavefunction.bloch`,
  `reaction.path`, vibrations, phonons and the other kinds without a
  builder are written through the generic member API (or
  `TypedKind::to_writer`).
* No `metadata.uuid` is generated (the crate has no random source). Set it
  with `set_root("metadata", …)`.
* ZIP64 is tested with forced ZIP64 records. No archive over 4 GiB has been
  written in the tests.
* `convert`, `prepare-job`, `export-xtb` and `run-xtb` of the Python CLI
  are not ported (format converters and program adapters). JSON output
  sorts object keys (the Python CLI keeps file order). `repack` and `scrub`
  stamp new entries with the current time instead of keeping each entry's
  time and attributes.
* `scrub` patterns are the `qvf::regex` subset of Python `re`; its `\w`,
  `\d` and `\s` use the Rust Unicode properties, which can differ from
  Python's for rare characters. The new `uuid` comes from the OS-seeded
  standard hasher keys, unique but not cryptographically strong.
* Typed models keep the keys of a member spec that the writer sets
  (`path`, `format`, `dtype`, `shape`, `sha256`, `size`, `media_type`) only
  through the `Tensor` and `OpaqueData` they decode to: `to_writer` writes
  members at default paths and does not carry other member-spec keys.
* `AtomPropertiesMetaPayload` entries other than the known metadata entry
  fall back to `serde_json::Value`, as in the C++ records.
* Finding messages follow the JavaScript wording but are not guaranteed to
  be byte-identical. Findings of one rule can come in another order, since
  JSON objects iterate with sorted keys here.
* JSON is parsed with `serde_json`: a duplicate key keeps the last value, as
  in the other validators, but a lone surrogate escape or a number outside
  the `f64` range is a JSON syntax error here, where JavaScript accepts it.
* The reader reads the whole archive into memory.
