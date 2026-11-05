# @qvf/core

A reader, validator and writer for [QVF](../spec/qvf-format-spec.md) (Quantum
Visualization Format) archives, in pure TypeScript. It runs in Node and in
browsers, has one runtime dependency (`fflate`, for inflate and deflate), and is
certified against the conformance corpus in `../conformance/`.

The validator runs every semantic
check of spec § 6 and the JSON Schema, through a checker generated from the
schema (`src/schemaCheck.ts` interpreting `src/schemaTable.ts`, which
`python3 scripts/gen_kinds.py` writes), with no JSON Schema engine. The package is `private` for now; use it
from a checkout (`npm ci && npm run build` in `js/`) or install it from a path
(`npm install ../path/to/qvf/js`).

Requires Node 20 or later (or a current browser). SHA-256 uses WebCrypto
(`crypto.subtle`), which browsers expose on `https:`, `localhost` and `file:`
pages, so the read and write functions are `async`.

## Node example

```js
import { readFile } from "node:fs/promises";
import { QvfReader, validate } from "@qvf/core";

const bytes = new Uint8Array(await readFile("result.qvf"));
const report = await validate(bytes);           // never throws on a bad archive
for (const f of report.findings) console.log(f.severity, f.code, f.where, f.message);
if (!report.ok) process.exit(1);

const qvf = QvfReader.fromBytes(bytes);         // throws QvfError if the container is unreadable
console.log(qvf.manifest.source.program, qvf.sections.map((s) => `${s.id} [${s.kind}]`));
const structure = await qvf.readJson("structure", "structure");   // sha256 verified
console.log(structure.atoms.map((a) => a.symbol).join(" "));

for (const s of qvf.sectionsOf("volume.density")) {
  const grid = await qvf.readJson(s.id, "grid");
  const { dtype, shape, data } = await qvf.readBinary(s.id, "data");   // Float32Array | Float64Array
  let max = -Infinity;
  for (const v of data) max = Math.max(max, Number(v));
  console.log(s.id, dtype, shape, "max", max, "origin (bohr)", grid.origin);
}
```

## Browser example

`examples/browser.html` lets you pick a `.qvf` file and prints its sections and
the validation findings. It loads the built ESM (`dist/index.js`) directly, with
no bundler; an import map points `fflate` at its browser build.

```sh
npm ci && npm run build
npx --yes http-server .            # or: python3 -m http.server, run in js/
# open http://localhost:8080/examples/browser.html
```

In your own page:

```html
<script type="importmap">{ "imports": { "fflate": "/node_modules/fflate/esm/browser.js" } }</script>
<script type="module">
  import { QvfReader, validate } from "/node_modules/@qvf/core/dist/index.js";
  const bytes = new Uint8Array(await file.arrayBuffer());   // file: a File from <input type="file">
  const { ok, findings } = await validate(bytes);
  const qvf = QvfReader.fromBytes(bytes);
</script>
```

The built package contains no Node-only imports (a test enforces it); with a
bundler, `fflate` resolves to its browser build on its own.

## API

| | |
|---|---|
| `QvfReader.fromBytes(bytes, limits?)` | Open an archive (`Uint8Array` or `ArrayBuffer`). Throws `QvfError` when the container is unreadable. |
| `reader.manifest` | The parsed `manifest.json` (typed `Manifest`). |
| `reader.sections`, `sectionsOf(kind)`, `section(id)` | Well-formed sections in manifest order, those of one kind, one by id. |
| `reader.readJson(section, role, verify = true)` | A JSON member, parsed. |
| `reader.readBinary(section, role, verify = true)` | `{ dtype, shape, data }`; `data` is the matching typed array (`Float32Array`, `Float64Array`, `Int8Array` ... `BigInt64Array` / `BigUint64Array` for 64-bit integers), flat and C-ordered, decoded from the little-endian bytes. |
| `reader.readBytes(section, role, verify = true)` | The raw bytes of a member (opaque members such as `run.record` `log`). |
| `reader.versionReadable` | `qvf_version` is one a consumer opens (`1`, and the withdrawn `2`). |
| `validate(bytes, { level?, limits? })` | `{ ok, findings, nSections }`. |
| `DEFAULT_LIMITS`, `resolveLimits()` | The resource caps of spec § 6.1. |

`verify` checks the member's `sha256` before the bytes are used and throws
`QvfError` with code `E_SHA256_MISMATCH` on a mismatch. Reading a member throws
a `QvfError` whose `code` is a finding code of the specification's Appendix B;
looking up a section or role that does not exist throws a `RangeError`.

### Validation

`validate` returns `findings: { severity, code, where, message }[]`. The codes
are exactly those of spec Appendix B (`E_SHA256_MISMATCH`, `W_ORPHAN_MEMBER`,
...), `where` is `archive`, `entry:<name>`, `manifest.json` or a JSON Pointer
into the manifest. The Python and C++ validators report the same code at the
same place. `level` is cumulative: `container` (ZIP, paths, size caps, from the
central directory only), `manifest`, `payloads` (reads every member), `semantics`
(structure periodicity, the `meta` entry of every `atom_properties` custom
role, the per-frame members, `coords` precision, binary `energies` and
`lattice` / `pbc` of `trajectory` and `reaction.path` sections, the `wavefunction.gto` rules of spec § 5.3 including the Appendix A.1
check of an `overlap` member, `E_NORMALIZATION`, a structure's `charge` and
`multiplicity`, `E_CHARGE_MULTIPLICITY`, the step ids of a `reaction.scheme`, a
grid `units.data` against the unit spec § 4.3 fixes for its volume kind,
`E_UNIT_CONTRADICTED`),
`references` (operands, `trajectory_ref`, `job_ref`, `files` index, the
participants of a `reaction.scheme`, the atom indices of a `bonds` section
against the structure its `structure_ref` names), `all` (the default). The bytes of a
`preview` member are checked against its `media_type` at the `payloads` level.

The JSON Schema is checked at the `manifest` level (`E_SCHEMA`, `where` the
JSON Pointer of the failing value) and at the `payloads` level for every JSON
member with a payload definition (`E_PAYLOAD_SCHEMA`, `where` the member, the
payload's own pointer in the message), at the same locations as the Python
validator. `validate(bytes, { schemaCheck: false })` skips it and reports
`W_SCHEMA_SKIPPED` once instead, the behaviour of a validator without a schema
engine.

### Limits

`Limits` has the field names and defaults of the Python and C++ `Limits`:
`max_member_bytes` (8 GiB), `max_total_bytes` (64 GiB), `max_ratio` (0: disabled; a positive value opts into a stricter application policy),
`ratio_exempt_bytes` (1 MiB), `max_json_bytes` (256 MiB), `max_json_depth` (64).
Caps are checked from the ZIP central directory before anything is
decompressed, and an entry is never inflated past the size the directory
declares for it. Pass overrides as `fromBytes(bytes, { max_json_depth: 32 })`.

### Container

STORE and DEFLATE entries, ZIP64 (archives over 4 GiB, more than 65535
entries; sizes are JavaScript numbers, exact up to 2^53), CRC-32 verification,
duplicate entry names refused (`E_ZIP_DUPLICATE_ENTRY`) except `manifest.json`,
whose last entry is the manifest of an archive updated in place (spec § 2.3),
a member's declared `size` checked against its entry, encrypted entries and
multi-volume archives unsupported. Archives with data prepended to the ZIP
(self-extracting stubs) are not supported.

## Writing archives

`Writer` builds an archive in memory or streams it to a sink. Members go into
the ZIP in the order they are added, each DEFLATE-compressed and checksummed
once, with its `sha256` and `size` in the manifest. `manifest.json` is written
last and stored uncompressed (spec § 2.1). Every JSON payload with a schema
definition is checked against the JSON Schema before it is written, and the
manifest before `finish()` writes it, so a finished archive passes the schema
stage of the validator. Binary members are C-ordered and little-endian on every
host (spec § 2.2).

```js
import { writeFile } from "node:fs/promises";
import { Writer } from "@qvf/core";

const w = new Writer({ program: "my-sketcher", version: "0.3.0", calculation: "h2o" });
w.setMetadata({ title: "Water", license: "CC-BY-4.0" });
await w.addStructure(
  [
    { symbol: "O", position: [0.0, 0.0, 0.1173] },
    { symbol: "H", position: [0.0, 0.7572, -0.4692] },
    { symbol: "H", position: [0.0, -0.7572, -0.4692] },
  ],
  { charge: 0, multiplicity: 1, name: "water" },
);
const rho = new Float32Array(40 * 40 * 40);                 // flat, C-ordered
await w.addDensity({ data: rho, shape: [40, 40, 40] }, [-5, -5, -5],
                   [[0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.25]],   // bohr
                   { sectionId: "rho", structureRef: "structure" });
await writeFile("water.qvf", await w.toBytes());
```

In a browser, save the archive as a `Blob`:

```js
const blob = await w.toBlob();
const a = Object.assign(document.createElement("a"), { href: URL.createObjectURL(blob), download: "water.qvf" });
a.click();
```

### Streaming output

Pass `sink` and every chunk goes out as soon as it is written. Nothing is kept,
so memory is bounded by the largest member rather than the archive:

```js
const handle = await showSaveFilePicker({ suggestedName: "run.qvf" });   // File System Access API
const writer = (await handle.createWritable()).getWriter();
const w = new Writer(source, { sink: (chunk) => writer.write(chunk) });
// ... add sections ...
await w.finish();
await writer.close();
```

In Node, `sink: (chunk) => new Promise((ok, err) => fileStream.write(chunk, (e) => (e ? err(e) : ok())))`
writes to a file stream. A returned promise is awaited before the next chunk.
Each member is held whole while it is compressed. Slab-wise input of one large
member (the Python `begin_binary`) is not available yet, because WebCrypto has
no incremental SHA-256.

### Deterministic output

`created_utc` and a random `uuid` are added to the root `metadata` unless
`setMetadata` gives them (spec § 3.5). The ZIP entries carry the same time
stamp. Fix both and the output is byte-identical from run to run:

```js
new Writer(source, { timestamp: "2026-01-01T00:00:00Z", uuid: "0f1e2d3c-4b5a-4968-8776-655443322110" });
new Writer(source, { autoMetadata: false });   // no time stamp, no uuid; entries dated 1980-01-01
```

### Writer API

| | |
|---|---|
| `new Writer(source, options?)` | `source` is `{ program, version, calculation? }`. Options: `autoMetadata` (default `true`), `timestamp`, `uuid`, `compression` (`"deflate"` or `"store"`), `level` (0 to 9, default 6), `forceZip64`, `sink`, `checkSchema` (default `true`). |
| `setMetadata(fields)`, `setProvenance(fields)` | Merge into the root `metadata` (spec § 3.5) and `provenance` (spec § 3.2) blocks. `null` removes a key. |
| `setRoot(block, value)`, `setEnergies`, `setViewerDefaults` | Any other root block: `thermochemistry`, `dipole_moment`, `constraints`, `energies`, `extensions`, `viewer_defaults`. |
| `addSection({ id, kind, ...fields })` | A section without members (any canonical or `x_<vendor>.*` kind). Returns the id. |
| `addJson(id, role, payload, { path?, mediaType? })` | A JSON member. |
| `addBinary(id, role, data, { dtype?, shape?, path?, mediaType? })` | A binary array member from a typed array or a nested array. `dtype` defaults to the typed array's own, `float64` for plain arrays, and values are converted (64-bit integers from `BigInt64Array` / `BigUint64Array`, bigints or safe integers). |
| `addOpaque(id, role, bytes \| string, { path?, mediaType? })` | An opaque member without dtype and shape (a log, BibTeX, a PNG). PNG and JPEG are stored uncompressed. |
| `annotateSection(id, { structureRef?, runRef?, group?, label? })` | The generic fields of spec § 4.1. |
| `buildManifest()` | The manifest as it would be written now. |
| `finish()` | Check and write the manifest and the ZIP directory. A schema refusal leaves the writer open, so the problem can be fixed and `finish()` called again. |
| `toBytes()`, `toBlob(type?)`, `bytes()`, `blob(type?)` | The archive (`finish()` first for the `to*` forms), not available with a `sink`. A `Blob` can exceed the largest single buffer. |
| `ZipWriter` | The container layer alone: whole entries, STORE or DEFLATE, ZIP64 when needed. |

The typed builders mirror the Python `qvf.Writer` methods. Keyword arguments
become camelCase options (`lattice_vectors` is `latticeVectors`, `section_id`
is `sectionId`). Payload fields keep the snake_case names of the specification.
Each returns the section id.

| Python | JavaScript | Kind |
|---|---|---|
| `add_structure(atoms, ...)` | `addStructure(atoms, options)` | `structure` (+ `bonds` member) |
| `add_trajectory(metadata, coords, ...)` | `addTrajectory(metadata, coords, options)` | `trajectory` with `lattice`, `velocities`, `gradient`, `times`, `energies`, `cellVolumes`, `temperatures` |
| `add_wavefunction_gto(shells, ...)` | `addWavefunctionGto(shells, options)` | `wavefunction.gto`, restricted (`moCoefficients`) or unrestricted (`moCoefficientsAlpha`, `moCoefficientsBeta`), with `overlap` and density matrices |
| `add_volume(kind, data, origin, voxel_vectors, ...)` | `addVolume(kind, data, origin, voxelVectors, options)` | every `volume.*` kind, `volume.series` with `axis` |
| `add_density`, `add_orbital` | `addDensity`, `addOrbital` | `volume.density`, `volume.orbital` |
| `add_excited_states(states, ...)` | `addExcitedStates(states, options)` | `excited_states` |
| `add_spectrum(kind, ...)` | `addSpectrum(kind, options)` | every `spectra.*` kind |
| `add_reaction_scheme(steps, ...)` | `addReactionScheme(steps, options)` | `reaction.scheme` |

Array arguments of the builders take nested arrays, a typed array whose shape
follows from the other arguments (`coords` from the number of atoms,
coefficients from the basis), or `{ data, shape }`, which is also what
`readBinary` returns. Volumes need the shape, so a flat typed array goes in as
`{ data, shape }`. Other kinds are written with `addSection` and the member
methods.

### Errors

A refused call writes nothing and leaves the writer usable. Refusals that
correspond to a finding of spec Appendix B throw `QvfError` with that code:
`E_PAYLOAD_SCHEMA` (a payload fails its schema definition, and the message names
the definition and the JSON pointers), `E_SCHEMA` (the manifest, in
`finish()`), `E_SECTION_ID_DUPLICATE`, `E_RESERVED_ID`, `E_KIND_UNKNOWN`,
`E_MEMBER_PATH`, `E_ZIP_DUPLICATE_ENTRY`, `E_DTYPE`, `E_BINARY_SIZE`,
`E_REF_UNRESOLVED` (a reference to a section not yet added, or of the wrong
kind), `E_STRUCTURE_PERIODICITY`, `E_CHARGE_MULTIPLICITY` and `E_SOURCE`. Other
bad arguments throw `TypeError` or `RangeError`. A failure while a member is
being written (a sink that throws) ends the writer: `finish()` then refuses.

## Tests and conformance

```sh
npm ci
npm test                      # build, unit tests, the conformance corpus
npm run conformance           # the positive corpus only, with a per-archive report
npm run conformance -- --negative ../conformance/negative --exact
npm run conformance -- --legacy ../conformance/legacy --exact
npm run conformance -- --negative ../conformance/negative --exact --no-schema-check
npm run typecheck             # sources, scripts and tests
npm run gen-kinds             # regenerate src/kinds.ts from ../registry.json
```

`npm test` runs `../conformance/corpus` and, when `../conformance/negative`
exists (or `QVF_NEGATIVE_CORPUS` names a directory), the must-fail corpus too,
requiring exactly the expected finding codes; the legacy corpus
`../conformance/legacy` (archives stamped `qvf_version: 2`, which this
validator checks as version 1 and reports with the warning
`W_VERSION_WITHDRAWN`) is run the same way, or from
`QVF_LEGACY_CORPUS`. `scripts/gen-kinds.mjs` generates
`src/kinds.ts`, the list of canonical kinds, from `registry.json`; a test fails
when the committed file is stale. `src/schemaTable.ts`, the compiled JSON
Schema, is generated by `python3 scripts/gen_kinds.py` at the repository root
(`--check` fails when it is stale), and `scripts/schema-check.mjs` serves the
parity test against the Python reference
(`python/tests/test_schema_checker_parity.py`).

`test/writer.test.mjs` round-trips archives from every builder through the
reader and the validator and checks deterministic output, forced ZIP64, more
than 65535 entries, streaming and refusals. `test/writer_python.test.mjs` runs
the Python reference validator (`python -m qvf validate --json`) on archives
the JS writer produced and compares tensors decoded by the Python reader bit
for bit. It is skipped when Python with `qvf`, `jsonschema` and `numpy` is not
importable (`pip install -e "python/[numpy,validate]"` from the repository
root). `QVF_PYTHON` names another interpreter.

`node scripts/validate-cli.mjs [--json] FILE.qvf ...` validates files from the
command line.

## License

Apache-2.0.
