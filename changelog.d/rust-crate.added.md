Add `rust/`, an independent Rust reader and validator (WP-61, in
development). It reads manifests, raw members and little-endian tensors, and
validates with the stable finding codes of Appendix B, including the JSON
Schema through a table that `scripts/gen_kinds.py` now also writes to
`rust/src/generated.rs`. Its conformance test passes the positive (78),
negative (118) and legacy (1) corpora with no skips. The crate is not in CI yet.
