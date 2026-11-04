# qvf (Rust)

An independent Rust implementation of the QVF archive format (WP-61): a
reader and a validator, certified by the shared conformance corpus in
[`../conformance`](../conformance). It does not use the C++ library or its C
ABI. Status: in development (`0.1.0-dev`, `publish = false`).

## Build and test

```sh
cd rust
cargo build
cargo test                      # unit tests, conformance runner, robustness test
cargo test --test conformance -- --nocapture   # prints the pass/fail/skip counts
cargo clippy --all-targets -- -D warnings
cargo run --example validate -- ../conformance/corpus/structure_periodic.qvf
```

Edition 2021, declared MSRV 1.74 (`rust-version` in `Cargo.toml`). The crate
is developed and tested with Rust 1.97. The MSRV has not been checked with an
older toolchain yet, and the crate is not part of the GitHub CI. Build output
stays in `rust/target/` (ignored).

## Dependencies

| Crate | Why |
|---|---|
| `flate2` (pure Rust `miniz_oxide` backend) | DEFLATE and CRC-32 |
| `serde_json` | `manifest.json` and the JSON members |
| `sha2` | member checksums (spec § 4.2) |

The ZIP container layer is our own port of `js/src/zip.ts` (about 400 lines)
instead of the `zip` crate. The validator must see every central-directory
record, duplicate names included, and check sizes and compression ratios
before anything is decompressed (spec § 2.1, § 2.3, § 6.1). It reads ZIP64
(extra fields, the ZIP64 end record and locator) and data-descriptor
archives, STORE and DEFLATE, and never inflates past the declared size.

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
* `validate_bytes` / `validate_file`: every finding of spec § 6 with its
  stable Appendix B code (`qvf::Code`), its severity and its location
  (`archive`, `entry:<name>`, `manifest.json` or a JSON pointer). Levels
  (`Level::Container` up to `Level::All`) and limits work as in the other
  validators. `schema_check: false` skips the JSON Schema check and reports
  `W_SCHEMA_SKIPPED`.

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

* No writer, no editor and no typed kind records yet (the WP-61 brief asks
  for structure, volume, wavefunction, trajectory, bands and spectra types).
* No `qvf-cli` binary or feature. `examples/validate.rs` prints JSON reports.
* Finding messages follow the JavaScript wording but are not guaranteed to
  be byte-identical. Findings of one rule can come in another order, since
  JSON objects iterate with sorted keys here.
* JSON is parsed with `serde_json`: a duplicate key keeps the last value, as
  in the other validators, but a lone surrogate escape or a number outside
  the `f64` range is a JSON syntax error here, where JavaScript accepts it.
* The reader reads the whole archive into memory.
* Not wired into GitHub CI, and the MSRV is declared but not verified.
