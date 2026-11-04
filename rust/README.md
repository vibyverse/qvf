# qvf (Rust)

An independent Rust implementation of the QVF archive format (WP-61): a
streaming writer with section builders, a reader and a validator, certified
by the shared conformance corpus in
[`../conformance`](../conformance). It does not use the C++ library or its C
ABI. Status: in development (`0.1.0-dev`, `publish = false`).

## Build and test

```sh
cd rust
cargo build
cargo test                      # unit tests, conformance runner, robustness and writer tests
cargo test --test conformance -- --nocapture   # prints the pass/fail/skip counts
cargo clippy --all-targets -- -D warnings
cargo run --example validate -- ../conformance/corpus/structure_periodic.qvf
```

Edition 2021, declared MSRV 1.74 (`rust-version` in `Cargo.toml`). The crate
is developed with current stable Rust. The release-candidate CI (`rust` job)
runs the tests on stable and on 1.74, and fmt and clippy on stable. Build output
stays in `rust/target/` (ignored).

## Dependencies

| Crate | Why |
|---|---|
| `flate2` (pure Rust `miniz_oxide` backend) | DEFLATE and CRC-32 |
| `serde_json` | `manifest.json` and the JSON members |
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
Pike VM for patterns). Do not edit `src/generated.rs` by hand:
regenerate with `python3 scripts/gen_kinds.py` and check freshness with
`python3 scripts/gen_kinds.py --check`.

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
* `Writer`: see *Writing* below.
* `validate_bytes` / `validate_file`: every finding of spec § 6 with its
  stable Appendix B code (`qvf::Code`), its severity and its location
  (`archive`, `entry:<name>`, `manifest.json` or a JSON pointer). Levels
  (`Level::Container` up to `Level::All`) and limits work as in the other
  validators. `schema_check: false` skips the JSON Schema check and reports
  `W_SCHEMA_SKIPPED`.

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
  transition-density states) and `ExcitedStates` (states, contributions,
  `nto_ref`, `density_ref`, `meta`).
* Refused before anything is written: illegal or duplicate ids and paths,
  unknown kinds, impossible charge and multiplicity, shapes and lengths that
  do not match, references to sections not added yet, shell centres outside
  the structure, contributions outside the reference orbitals. A tensor sink
  that gets too many or too few bytes, or an I/O error, poisons the writer so
  that `finish` cannot write a manifest naming incomplete data.

`tests/writer.rs` writes archives with every builder, as DEFLATE, STORE and
forced ZIP64. Each must have no findings from the Rust validator and from the
Python reference validator (`python3 -m qvf validate --json`, run with
`PYTHONPATH=../python`, or `$PYTHON`), and must read back unchanged through
the Rust reader. The Python check is skipped with a `SKIP` line when Python
or the `qvf` package is unavailable, and accepts only `W_SCHEMA_SKIPPED`
when `jsonschema` is missing.

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
| `conformance/corpus` | 78 | 0 | 0 |
| `conformance/negative` | 118 | 0 | 0 |
| `conformance/legacy` | 1 | 0 | 0 |

On all 197 archives the set of (code, severity, location) findings and the
verdict equal those of the Python reference validator (`validate_qvf`, with
`jsonschema` installed). That comparison was run once for this change and is
not an automated test.

## Known gaps

* No editor (in-place update, spec § 2.3) and no typed reader records yet.
  The writer has builders for structure, wavefunction.gto, volumes and
  excited states. Trajectory, bands and spectra builders (the WP-61 brief)
  remain. Complex and `generalized` wavefunctions and `wavefunction.bloch`
  are written only through the generic member API.
* No `metadata.uuid` is generated (the crate has no random source). Set it
  with `set_root("metadata", …)`.
* ZIP64 is tested with forced ZIP64 records. No archive over 4 GiB has been
  written in the tests.
* No `qvf-cli` binary or feature. `examples/validate.rs` prints JSON reports.
* Finding messages follow the JavaScript wording but are not guaranteed to
  be byte-identical. Findings of one rule can come in another order, since
  JSON objects iterate with sorted keys here.
* JSON is parsed with `serde_json`: a duplicate key keeps the last value, as
  in the other validators, but a lone surrogate escape or a number outside
  the `f64` range is a JSON syntax error here, where JavaScript accepts it.
* The reader reads the whole archive into memory.
