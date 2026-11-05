Add four jobs to the release-candidate CI workflow, with the existing
`release-candidate/**` push and manual dispatch triggers only: `rust` (the
WP-61 crate on stable and on its MSRV 1.74, with fmt and clippy on stable),
`parity` (the schema-checker parity test with Python, Node and the C++ tool
in one job, failing on any skip), `fuzz` (both C++ libFuzzer targets for
180 s each under ASan + UBSan, crashers uploaded as an artifact) and `s390x`
(cross-build, CTest under QEMU and the cross-host round trip, WP-50).
Development pushes do not run them. Each job's shell steps ran locally (the
fuzz targets for 30 s each) and passed. The Rust tests also pass on 1.74, so
the declared MSRV is now checked. The Python and JavaScript fuzz harnesses
were not on main and are not included.
