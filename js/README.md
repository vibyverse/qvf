# @qvf/core

A reader and validator for [QVF](../spec/qvf-format-spec.md) (Quantum
Visualization Format) archives, in pure TypeScript. It runs in Node and in
browsers, has one runtime dependency (`fflate`, for inflate), and is certified
against the conformance corpus in `../conformance/`.

This package **reads and validates**; it has no writer. It runs every semantic
check of spec § 6 and the JSON Schema, through a checker generated from the
schema (`src/schemaCheck.ts` interpreting `src/schemaTable.ts`, which
`python3 scripts/gen_kinds.py` writes), with no JSON Schema engine. The package is `private` for now; use it
from a checkout (`npm ci && npm run build` in `js/`) or install it from a path
(`npm install ../path/to/qvf/js`).

Requires Node 20 or later (or a current browser). SHA-256 uses WebCrypto
(`crypto.subtle`), which browsers expose on `https:`, `localhost` and `file:`
pages, so the read functions are `async`.

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

`node scripts/validate-cli.mjs [--json] FILE.qvf ...` validates files from the
command line.

## License

Apache-2.0.
