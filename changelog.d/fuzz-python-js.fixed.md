Fuzzed the Python and JavaScript readers and validators with structure-aware
harnesses (`python/tests/fuzz/`: Hypothesis and atheris, `js/test/fuzz/`: a
seeded Node runner, both opt-in, see CONTRIBUTING.md) and fixed what they
found. Python no longer inflates a DEFLATE entry past its declared size
(`zipfile.read` used to inflate a lying stream in full before truncating
it), refuses an entry that inflates to fewer or more bytes than declared,
reads only STORE and DEFLATE members (it used to accept bzip2 and LZMA),
refuses an end of central directory record that names another disk or counts
a different number of entries than the central directory holds, and
parses JSON as RFC 8259 UTF-8 JSON (it used to accept `NaN`, `Infinity`,
UTF-16 and UTF-32). Damaged, encrypted or unsupported entries are now a
`ValueError` with `E_ZIP_UNREADABLE` from `Reader`, `QvfReader` and the
streaming API instead of `zlib.error`, `RuntimeError`, `NotImplementedError`
or `EOFError`. The error is also a `zipfile.BadZipFile`. `Tensor.tolist()`
refuses a zero-extent shape that would need more than 2^20 nested empty
lists. JavaScript refuses an entry whose declared size exceeds what DEFLATE
can expand its stored bytes to before allocating the output buffer (a tiny
archive could make it allocate 8 GiB). Seven negative archives pin the
container and JSON rules: `entry_inflates_past_size`, `entry_size_overstated`,
`entry_size_impossible`, `compression_method_unsupported`, `zip_multi_volume`
(`E_ZIP_UNREADABLE`), `manifest_json_nan` and `manifest_json_utf16`
(`E_MANIFEST_JSON`). Python, JavaScript, C++ and Rust agree on all seven.
