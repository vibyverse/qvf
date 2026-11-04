Add a second C++ libFuzzer target, `manifest_fuzzer` (`-DQVF_FUZZ=ON`). It
builds a consistent ZIP from a manifest.json and member payloads, with
optionally corrected member sha256 values, so mutations reach the JSON
parser, the typed manifest model, the reader and every validation stage.
`scripts/fuzz_manifest_seeds.py` turns .qvf archives into its seeds.
CONTRIBUTING.md documents both fuzz targets.
