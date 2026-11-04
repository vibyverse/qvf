Complete the WP-61 scope of the `rust/` crate (fourth tranche). The writer
gains builders for `trajectory` (including periodic MD with a fixed or
per-frame cell, velocities, gradients, times, energies, cell volumes and
temperatures, and float32 coordinates), `bands`, `dos.total`,
`dos.projected`, `dos.coop`, `dos.cohp` and the six one-dimensional
`spectra.*` kinds, plus `add_spectrum_payload` for any spectra kind. Each
JSON payload is checked against its schema definition before it is written.
`qvf::Editor` updates an archive in place or atomically per spec § 2.3
(append, replace members under fresh paths, remove sections or members, set
the run status, roll back on abort) and `Editor::repack` streams the
referenced members into a new archive without orphans. The `qvf` binary
adds `diff`, `repack`, `scrub` (with `--path-pattern` and `--user-pattern`
through `qvf::regex`, a Python `re` subset), `conformance` and `schema`.
`rust/tests/cli.rs` matches the Python CLI on all 201 corpus archives for
these commands (1810 commands), and the editor tests check orphan, size and
repack rules against the Python validator and editor. `cargo +1.74 test`
passes. Remaining Python CLI commands not ported: `convert`, `prepare-job`,
`export-xtb` and `run-xtb`.
