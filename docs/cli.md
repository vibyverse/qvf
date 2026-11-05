# QVF command line

The `qvf` command uses the reference Python library to inspect, compare and
maintain archives. Install from the toolkit checkout with
`pip install './python[validate]'`. The `validate` extra enables full JSON
Schema checks. NumPy is unnecessary for these commands. `python -m qvf`
works identically, and the original `qvf-validate [--json] FILE` stays available.

```sh
qvf info result.qvf
qvf validate result.qvf --json
qvf cat result.qvf
qvf cat result.qvf structure structure
qvf extract result.qvf rho data --output density.bin
qvf diff before.qvf after.qvf --atol 1e-8 --rtol 1e-6
qvf repack result.qvf --output compact.qvf --compression deflate
qvf schema > manifest.schema.json
qvf conformance conformance/corpus
```

`info` displays section IDs and kinds, member roles and archive paths, formats,
tensor shapes and uncompressed sizes. `--json` also includes compressed sizes.
It reads the manifest and central directory only. `cat` prints the manifest
when no section/member pair is given, or decodes one JSON member when both
are present. `extract` writes one declared member's exact bytes to an explicit
new file, including opaque logs, images and little-endian tensors. It checks
sizes, CRC and SHA-256 before publishing that file.

Exit status is **0** for success/equal, **1** for failed validation, differing
archives or failed conformance, and **2** for usage or operational errors.

Validation reports include `validation_context`: the archive format and stated
schema URI, executing reader/toolkit version, selected bundled schema ID and
SHA-256, and whether schema checking ran. The archive URI is informational;
it does not select or download a schema. Format version 1 is not toolkit
release 0.1: compatible toolkit releases can carry different schema additions.
Compare the producer's pinned validator and schema digest when investigating
cross-version rejections. An unsupported archive format is an `E_VERSION`
error; other schema violations remain errors, not presumed version mismatches.
Null `job.spec` basis/functional values mean the method takes neither input
and are allowed by the current schema. This does not allow null in arbitrary
typed fields. Format validation establishes neither numerical accuracy nor
electronic stability.
A partially verified extraction or archive rewrite is discarded. Output files
must be new: `extract`, `repack`, `scrub` and `prepare-job` refuse existing paths, including
the input and symlinks. They stage a sibling temporary file and publish it
without replacing a concurrently created destination. No command edits an
input archive in place.

## Export xTB results (toolkit 0.7.0)

The `export-xtb` command packages molecular GFN0/1/2-xTB output as a new QVF:

```sh
qvf export-xtb xtbout.json --geometry xtbopt.xyz --charge 0 \
  --trajectory xtbopt.log --wbo wbo -o result.qvf
```

JSON output must be requested before running xTB. Optional `--molden` plus
`--spin`, `--modes` and repeatable `--cube FILE KIND` arguments add orbitals,
normal modes and scalar grids. The command validates before publishing and
refuses existing destinations. It neither runs xTB nor infers convergence.
See the [xTB exporter guide](adapters/xtb.md) for installation, accepted
formats and data that cannot be recovered from JSON alone.

`qvf run-xtb xtb-job.json --binary /opt/xtb/bin/xtb` executes a prepared
serial bundle, checks xTB termination/convergence evidence and exports
`result.qvf`. See the same guide for the request contract and failure behavior.

## Prepare one pending job (toolkit 0.6.0)

An existing one-species drawing or imported geometry can become a pending
calculation request while retaining its drawing, bonds and preview:

```sh
qvf prepare-job drawing.qvf --output job.qvf --job-type molecular \
  --method rhf --task single_point --basis sto-3g --no-functional
qvf validate job.qvf --json
```

This command is included in toolkit 0.6.0.
It requires `qvf-writer[validate]` for full source and output validation;
NumPy is unnecessary. `--job-type`, `--method` and `--task` are required.
Tasks are `single_point`, `optimize` or `hessian`. Method, basis and functional
names remain engine-dependent; preparing a request neither checks scientific
suitability nor runs or submits a calculation. Arbitrary engine options are
outside this command's scope.

Omitted `--basis` and `--functional` remain absent, leaving the engine's choice.
`--no-basis` and `--no-functional` write explicit JSON null, meaning the method
takes none; each conflicts with its named value option. `--charge` and
`--multiplicity` default to 0 and 1 as defined by `job.spec`. They are never
inferred from the geometry: a charged species or radical requires matching
explicit values. Periodic geometry requires `--job-type periodic`, with an
optional three-integer `--kpoints N1 N2 N3` mesh.

The source must contain exactly one non-empty structure, optionally bonds,
preview, citations and vendor sections. Existing jobs, runs, other canonical
sections, result metadata and calculation constraints are refused. Root
provenance may contain charge, multiplicity, dimensionality and vendor keys;
other calculation/lifecycle provenance is refused. Conflicting identity or
periodicity, invalid source bytes and incomplete schema validation fail before
publishing output. No section is selected or silently dropped.

Every original section and non-manifest member is preserved, including opaque
vendor documents and orphan entries; their declared hashes stay unchanged.
Superseded manifests are consolidated. The new manifest adds one `job.spec`,
sets the pending lifecycle, and records a new UUID derived from the original
UUID when available. JSON output reports the job/structure IDs and validation
findings, including inherited warnings.
The withdrawn legacy version-2 spelling is read with its existing warning;
the prepared archive always writes format version 1.

The equivalent API is `qvf.jobs.prepare_job(source, destination,
job_type="molecular", method="rhf", task="single_point", basis="sto-3g",
functional=None)`. See the [family workflow](family_workflow.md) for submission.

## Choose and record the validator

For acceptance against a producer's QVF pin, install that exact toolkit release
with the `validate` extra and run it in the intended Python environment.
`qvf_version: 1` identifies the archive format, not a toolkit release; the
schema `$id` also stays the same across compatible additions. Neither selects
the validator version. A later-toolkit compatibility check is a separate result.

Record the source tag and commit, or the installed distribution's artifact
hash, alongside the toolkit version, packaged schema hash and validation report.
A development checkout can still report the previous toolkit version, so the
version string alone is insufficient. The existing commands expose the version,
schema and findings:

```sh
python -m qvf --version
python -m qvf schema > validator.schema.json
python -c 'import hashlib, pathlib; print(hashlib.sha256(pathlib.Path("validator.schema.json").read_bytes()).hexdigest())'
python -m qvf validate result.qvf --json
```

Use a new output path for `validator.schema.json`. `schema` prints the schema
bundled with that interpreter's QVF package. Keep the original archive and its
hash unchanged while investigating a failure. Confirm the report has no
`W_SCHEMA_SKIPPED` before claiming full schema validation.

For example, since toolkit **0.3.0**, a `job.spec` `basis` or `functional` of
`null` means the method takes none; omission asks the engine to choose, and an
empty string is a deprecated spelling of `null`. An older validator can reject
these valid fields. For the [frozen family interfaces](family_workflow.md), use
**v0.6.0** and retain the original scientific results; do not invent a basis or
functional, remove fields or rerun a calculation to satisfy an older schema.

## Validation levels

`qvf validate FILE --level LEVEL` selects a cumulative stage. The library
accepts the same names in `qvf.validate.run(FILE, level="manifest")` and
`reader.validate(level="manifest")`.

| Level | Checks added |
|---|---|
| `container` | ZIP directory, duplicate and illegal paths, resource caps, presence of manifest; no entries decompressed |
| `manifest` | Manifest JSON/depth/schema, section IDs/kinds, declared member existence and size; only the manifest is read |
| `payloads` | Member CRC/SHA-256, binary sizes/dtypes, JSON parsing and payload schemas, run text UTF-8, preview signatures |
| `semantics` | Scientific consistency, shapes across members, units and normalization |
| `references` | Cross-section reference findings |
| `all` | All stages; the default, identical to the full reference validator |

Findings use the existing shared codes, with `severity`, `code`, `where` and
`message` in `--json` output. `W_SCHEMA_SKIPPED` means the optional schema engine
was unavailable; it is not evidence that structural schema checks passed.
The early stages avoid later reads and checks. A successful early stage does
not certify stages that have not run. Container-only reports have zero parsed
sections. Stage boundaries follow the C++ facade; JavaScript retains its full
validation entry point.

## Compare and repack

`diff` compares logical manifest metadata and common member contents. ZIP
compression, entry ordering, superseded manifests, orphans, and member transport
paths/sizes/hashes do not themselves cause a difference. Section order, kinds,
roles, dtypes and shapes are significant. JSON object key order is immaterial.
Members being compared are checksum-verified; use `validate` to check an entire
archive independently of a comparison.

Floating JSON numbers and tensor values use
`abs(a-b) <= max(atol, rtol * max(abs(a), abs(b)))`, with both tolerances zero
by default. Integers, booleans and opaque bytes compare exactly. NaNs differ;
equal infinities compare equal. `--ignore-metadata` excludes only the root
`metadata` object (including UUID, timestamps and authors). Source information
and other root fields still count. `--json` emits the verdict and locations;
`--max-differences N` limits reported locations, without stopping comparison
or checksum verification. Binary comparisons use bounded chunks.

`repack` copies only members referenced by the current manifest and writes one
manifest, removing orphan entries and old manifests. It preserves scientific
content, source metadata and UUID. `--compression preserve` is the default;
`store` and `deflate` select a new ZIP method. This is useful after append-mode
library edits, and `diff` can confirm that the logical content is unchanged.

## Scrub for review

Scrubbing always reports changes and requires either `--dry-run` or an explicit
new `--output`. The original archive remains available to undo or review the
operation. First inspect the proposed changes:

```sh
qvf scrub result.qvf --dry-run \
  --path-pattern '/scratch/example' --user-pattern 'user=example'
qvf scrub result.qvf --output shared.qvf \
  --path-pattern '/scratch/example' --user-pattern 'user=example'
qvf validate shared.qvf
```

The command removes deprecated root `provenance.hostname` and author `email`
and `affiliation` fields. Repeated `--path-pattern REGEX` and
`--user-pattern REGEX` options select additional replacements in `run.record`
commands, UTF-8 inputs/logs/attachments, and string values of the `files` index.
They use Python regular expressions. There is no built-in list of people,
paths or hosts. The replacement is the literal `[redacted]` by default;
`--replacement TEXT` changes it. Pattern matches and removals are reported by
location and count, without echoing matched text.

Text processing is bounded by the reader's JSON/text cap per member; changed
payloads are staged on disk, keeping memory bounded to one text member.
Shared payloads that would change are rejected rather than silently changing
another section's data. The new archive has updated member sizes/hashes, one
manifest, no orphans, a new UUID and a `derived_from` link to the original UUID
when available. `job.spec` payload bytes and scientific sections are preserved.

The report lists skipped members, including non-UTF-8 attachments and run text
for which no patterns were selected. Names in other fields, section IDs,
archive member paths, author names, arbitrary vendor payloads and binary
attachments are outside the scrub scope. Review them before sharing; a scrub
report does not certify that an archive contains no personal information.
This command does not alter the repository's publication-scanner policy.

## Run a corpus

`qvf conformance DIRECTORY` runs the positive `*.expected.json` cases there,
then `negative/` and `legacy/` siblings if present. `--negative DIRECTORY`
and `--legacy DIRECTORY` select explicit alternatives. `--no-negative` runs
only positives (also skipping automatic legacy cases); `--no-legacy` skips
just legacy. Corpus archives must resolve within their selected directory.
The corpus comes from the QVF source release and is not bundled into the
Python wheel. The installed command and `conformance/run_conformance.py`
share one runner and the same expected findings; neither needs NumPy.

## Text conversion (toolkit 0.6.0)

`qvf convert INPUT --from xyz --to qvf -o OUTPUT` imports geometry; use
`--from qvf --to extxyz --section ID` to export a selected structure or
trajectory. Scalar cube import/export supports explicit kind and periodic
endpoint choices. Toolkit 0.6.0 also supports static XSF geometry
and scalar grids, plus single-model PDB geometry with explicit element
columns and no alternate locations. Toolkit 0.7.0 adds bounded
molecular Molden conversion, where imports require an explicit `--spin`
interpretation, and PDB atom-site identity, alternate locations and model
ensembles.
Every operation returns a conversion report and requires a new output path.
See [supported subsets and
examples](converters.md).
