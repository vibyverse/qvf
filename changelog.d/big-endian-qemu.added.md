Added a local big-endian check of the C++ library (WP-50, no CI job): a
CMake toolchain file for s390x under QEMU user-mode emulation, the
`qvf_endian_probe` CTest test that writes all ten dtypes through every write
path and compares them with little-endian bytes, and
`scripts/check_endian_roundtrip.py`, which validates archives written on one
host with the Python validator and compares their tensors bit for bit with
those of the other host. Under emulated s390x all 13 CTest tests pass,
including the conformance runner (78 positive, 118 negative, 1 legacy
archive). No endianness defect was found, so the C++ sources are unchanged.
