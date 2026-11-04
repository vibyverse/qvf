Add a writer to the JavaScript package `@qvf/core` (in development). `Writer`
runs in Node and browsers with no Node-only API: it builds the manifest, writes
JSON members and little-endian binary tensors with dtype and shape, computes
`sha256` through WebCrypto and `size` for every member, deflates through
`fflate`, writes `manifest.json` last and stored, and switches to ZIP64 when a
size, an offset or the entry count needs it. Output is a `Uint8Array`, a
`Blob`, or a stream of chunks to a sink. Typed builders follow the Python
writer: `addStructure`, `addTrajectory`, `addWavefunctionGto` (restricted and
unrestricted), `addVolume` with `addDensity` and `addOrbital` (every
`volume.*` kind), `addExcitedStates`, `addSpectrum` (every `spectra.*` kind)
and `addReactionScheme`. Payloads are checked against the JSON Schema before
they are written, the manifest before `finish()` writes it. A fixed
`timestamp` and `uuid` give byte-identical output. Tests round-trip every
builder through the JS reader and validator, and the Python validator accepts
the archives (a Node test that is skipped without Python). Slab-wise input of
a single large member is not available yet, because WebCrypto has no
incremental SHA-256.
