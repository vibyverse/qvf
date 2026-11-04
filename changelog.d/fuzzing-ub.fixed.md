Fix two undefined-behaviour defects in the C++ library found by a
coverage-guided fuzzing campaign (WP-50). `ZipReader::local_header` passed a
null pointer to `memcmp` for an entry with an empty name. `Json::as_int`
converted a double outside the int64 range unchecked (a manifest shape
integer beyond int64 reached it through `QvfReader::read_binary_member`); it
now throws. Both have regression tests.
