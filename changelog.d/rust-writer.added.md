Add a writer to the `rust/` crate (WP-61, second tranche). `qvf::Writer`
streams members into a ZIP that needs no seekable sink (data descriptors,
ZIP64 fields and end records when needed or forced) and writes
`manifest.json` STORED and last. Members are JSON, little-endian tensors of
every manifest dtype with `dtype`, `shape`, `sha256` and `size` (also
streamed in slabs), and opaque bytes. Section builders cover `structure`,
`wavefunction.gto` (restricted, restricted open-shell and unrestricted, with
`overlap`, density matrices and NTO links), every `volume.*` kind and
`excited_states`. `rust/tests/writer.rs` checks that the written archives
have no findings from the Rust validator and the Python reference validator
and read back unchanged. No new dependency. Trajectory, bands and spectra
builders and an editor remain.
