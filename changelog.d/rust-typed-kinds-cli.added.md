Add typed section models and a `qvf` command-line tool to the `rust/` crate
(WP-61, third tranche). `scripts/gen_kinds.py` now also writes
`rust/src/kinds.rs` (rendered by `scripts/gen_typed_kinds.py` with the
mapping of the C++ typed records): serde models of the schema `$defs`
payloads and one struct per canonical kind, read with
`Reader::get::<K>(id)` and written back with `TypedKind::to_writer`, plus
`Reader::get_volume` for any `volume.*` kind. Every section of the positive
and legacy corpora round-trips through the models (`rust/tests/typed.rs`).
The `qvf` binary has `validate` (human and `--json`, every `--level`),
`info`, `cat` and `extract`, with the report layout and exit codes of the
Python CLI. `rust/tests/cli.rs` runs both CLIs on all 201 corpus archives
and requires the same exit codes, finding codes and findings. New
dependency: `serde` with `derive`. An editor and the other Python CLI
commands remain. The crate is not in CI yet.
