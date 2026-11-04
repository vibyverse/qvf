# libqvf — design of the library other codes link

> Status: design document (informative). It describes the target shape of the
> toolkit after the work packages in `work_packages.md` land. The original
> `qvf::QvfWriter`, `qvf::QvfReader`, `qvf_c.h`, and Python modules were the
> compatibility seed, not the current implementation limit. As of 2026-10-19,
> toolkit 0.5.0 includes the C++/C/Python facades and an initial Fortran slice;
> see [the roadmap](../ROADMAP.md) for delivered and remaining scope. Public
> legacy names remain available; the target APIs below are informative and
> do not replace the released [library guide](library_guide.md).

## 1. Goals and non-goals

**Goals**

1. **One library, every direction.** Write, read, validate, and update `.qvf`
   archives from C++, C, Fortran, Python, and (reader-first) JavaScript and
   Rust, with identical semantics everywhere, proven by the conformance corpus.
2. **Linkable by any code.** Zero external dependencies in the core, C++17,
   static or shared, single-header option, stable C ABI with symbol
   versioning, CMake/pkg-config/vcpkg/conan/PyPI distribution.
3. **Scales to the data.** Stream a multi-gigabyte grid straight from a
   producer's buffer into the archive with one copy and no intermediate file;
   read a member lazily; update an archive in place.
4. **Typed, not stringly.** Every canonical kind has a typed representation
   generated from the schema, so a producer cannot misspell a role or get a
   dtype wrong, and a consumer gets arrays with shapes it can trust.
5. **Safe on untrusted input.** The reader and validator are fuzzed and enforce
   size, ratio, depth, and path limits.

**Non-goals**

- Not a visualisation library; no grid evaluation of orbitals beyond the
  normalization helpers that already exist (consumers like vibe-view own that).
- Not a workflow engine; the editor implements the job-spec lifecycle rules,
  nothing schedules jobs.
- Not a general ZIP library; only what QVF needs (STORE, DEFLATE, ZIP64,
  data descriptors).

## 2. Layering

```
 ┌───────────────────────────────────────────────────────────────┐
 │ L7  CLI `qvf` (Python)      L6  convert (cube, xyz, molden…)   │
 ├───────────────────────────────────────────────────────────────┤
 │ L5  bindings: C ABI · Fortran module · Python (pure) · JS · Rust│
 ├───────────────────────────────────────────────────────────────┤
 │ L4  facades: Writer (streaming) · Reader (lazy) · Editor       │
 ├───────────────────────────────────────────────────────────────┤
 │ L3  kinds: typed section structs ⇄ manifest + members          │
 │     (generated tables from the schema + hand-written semantics)│
 ├───────────────────────────────────────────────────────────────┤
 │ L2  validate: container → manifest → payload → semantic → refs │
 ├───────────────────────────────────────────────────────────────┤
 │ L1  model: Manifest, Section, MemberSpec, Json (parse+emit),   │
 │     Kind registry, units vocabulary, Lattice/Grid helpers      │
 ├───────────────────────────────────────────────────────────────┤
 │ L0  container: ZipWriter/ZipReader (ZIP64, streaming), SHA-256,│
 │     CRC-32, DEFLATE/INFLATE, byte-order conversion             │
 └───────────────────────────────────────────────────────────────┘
```

The C++ library implements L0–L4 and is the reference. The Python package is a
second full implementation of L0–L4 on the standard library (it is what
adopters prototype with, and two independent implementations are what make
the corpus meaningful). Bindings in L5 wrap the C++ library through the C ABI
(Fortran, optional Python accelerator, Julia) or are independent
implementations certified by the corpus (JavaScript, Rust).

## 3. L0 — container

```cpp
namespace qvf::container {

enum class Method { Store, Deflate };

class ZipWriter {               // sequential, streaming, ZIP64-capable
public:
    explicit ZipWriter(std::unique_ptr<Sink> out);   // file, memory, or callback sink
    // Begin an entry; bytes are pushed through `EntrySink`; CRC-32, sizes and
    // the sha256 are accumulated incrementally. On seekable sinks the local
    // header is patched afterwards; on non-seekable sinks a data descriptor
    // (general-purpose bit 3) is written. Both are valid per APPNOTE.
    EntrySink begin(std::string path, Method m, bool force_zip64 = false);
    void      write(std::string path, std::span<const uint8_t> bytes, Method m);
    void      finish();         // central directory (+ ZIP64 records when needed)
};

class ZipReader {               // random access over a seekable source or memory
public:
    static ZipReader open(const std::string& path);
    static ZipReader from_bytes(std::vector<uint8_t>);
    const std::vector<Entry>& entries() const;          // name, sizes, method, offset
    bool contains(std::string_view name) const;         // duplicate names → error at open,
                                                         // except manifest.json: last wins (spec § 2.3)
    std::vector<uint8_t> read(std::string_view name, const Limits&) const;
    void read_into(std::string_view name, std::function<void(std::span<const uint8_t>)>) const; // streaming inflate
};

struct Limits {                  // enforced by ZipReader and the validator
    uint64_t max_member_bytes   = 8ull << 30;   // spec § 2.1
    uint64_t max_total_bytes    = 64ull << 30;
    double   max_ratio          = 1000.0;       // compressed → uncompressed
    size_t   max_json_bytes     = 256u << 20;
    size_t   max_json_depth     = 64;
};
}
```

Also in L0: `Sha256` (incremental), `crc32`, `deflate`/`inflate` (existing
codecs, made streaming), and `endian::to_le(span, itemsize)` used by the
writer when the host is big-endian.

## 4. L1 — model

As shipped (WP-41): `qvf/model.hpp` (C++, header-only) and `qvf.model`
(Python dataclasses; `qvf_model.py` remains a compatibility import). The model is **lossless**: `from_json` types a field
only when it has the JSON type the schema gives it and keeps everything else
verbatim in the `extra` object of the enclosing record (kind-specific section
fields, vendor keys, wrongly typed fields), and `to_json` writes the typed
fields followed by `extra`, so `Manifest::from_json(j).to_json() == j` for any
manifest object (round-trip tests over every conformance manifest, and over
mutated ones in Python). Absent and present are told apart, so every optional
field is a `std::optional`.

```cpp
namespace qvf {

struct Source { std::string program, version, calculation; };   // moved here from qvf.hpp

struct MemberSpec {              // one entry of a section's `members`
    enum class Format { Json, Binary, Opaque };   // Opaque: format "binary" without dtype
    std::optional<std::string> path;
    std::optional<Format> format;
    std::optional<std::string> dtype;              // Binary
    std::optional<std::vector<int64_t>> shape;     // Binary
    std::optional<std::string> sha256;
    std::optional<uint64_t> size;                  // P19
    Json extra;                                    // every other key
    static MemberSpec from_json(const Json&);  Json to_json() const;
};

struct Section {
    std::optional<std::string> id, kind, label, component, structure_ref, run_ref, group;
    std::optional<bool> critical;
    // role → spec, in manifest order
    std::optional<std::vector<std::pair<std::string, MemberSpec>>> members;
    Json extra;                  // kind-specific peer fields (ao_metadata, program, …)
    const MemberSpec* member(std::string_view role) const;
};

struct Manifest {
    std::optional<int64_t> qvf_version;
    std::optional<std::string> schema_uri, byte_order;
    std::optional<Source> source;
    std::optional<Json> provenance, metadata, viewer_defaults, thermochemistry,
        dipole_moment, constraints, energies, extensions;   // root blocks (P03, P05)
    std::optional<std::vector<Section>> sections;
    Json extra;
    const Section* section(std::string_view id) const;
    static Manifest from_json(const Json&);  Json to_json() const;
};
// QvfReader::typed_manifest() exposes it beside the raw QvfReader::manifest().

// Generated from the schema (scripts/gen_kinds.py, qvf/kinds_generated.hpp) —
// never hand-edited. C++17 has no std::span, hence pointer + count:
namespace gen {
enum class Format { Json, Binary, Opaque };
struct RoleInfo { std::string_view role; bool required; Format format;
                  const char* const* dtypes; size_t n_dtypes; int min_rank, max_rank;
                  std::string_view payload_def; bool pattern; };  // pattern: role is a regex (custom.<name>)
struct KindInfo { std::string_view kind; const RoleInfo* roles; size_t n_roles;
                  bool allow_extra_members; };
const KindInfo* kind_info(std::string_view kind);     // nullptr for unknown and vendor kinds
inline constexpr std::string_view kCanonicalKinds[] = { … };   // registry.json order
std::string_view payload_def(std::string_view kind, std::string_view role);

// Units vocabulary (P02), also generated; qvf/units.hpp re-exports it as qvf::Unit
// and keeps the conversions:
enum class Unit { Angstrom, Bohr, Hartree, ElectronVolt, … };
std::string_view unit_name(Unit);  std::optional<Unit> parse_unit(std::string_view);
// plus the finding codes of spec Appendix B and the schema $id.
}

// Geometry helpers that encode the two conventions (D06) once:
struct Lattice {                 // always rows, Å internally
    static Lattice from_rows_angstrom(const double (*rows)[3]);
    static Lattice from_columns_bohr(const double (*cols)[3]);   // reaction.path.lattice
    std::array<std::array<double,3>,3> rows_angstrom() const;
    std::array<std::array<double,3>,3> columns_bohr() const;
    double volume_angstrom3() const;
};
struct GridSpec { Vec3 origin_bohr; Mat3 voxel_bohr; std::array<int64_t,3> shape;
                  std::array<bool,3> periodic{}; bool includes_endpoint = true; };
}
```

`qvf::Json` keeps its current emitter and parser and gains a JSON-pointer
accessor (`at_pointer("/atoms/0/symbol")`) used by the validator and the
conformance runner.

## 5. L2 — validation

```cpp
namespace qvf::validate {
enum class Level { Container, Manifest, Payloads, Semantics, References, All };
struct Finding { enum class Severity { Error, Warning } severity;
                 std::string code;        // stable id, e.g. "E_SHA256_MISMATCH", "W_ORPHAN_MEMBER"
                 std::string where;       // "sections/3/members/data", "member:volumes/rho.dat"
                 std::string message; };
struct Report { bool ok; std::vector<Finding> findings; size_t n_sections; };
Report run(const ZipReader&, Level, const Limits& = {});
}
```

The structural check (what `jsonschema` does in Python) is **generated**
(WP-41): `scripts/gen_kinds.py` compiles the schema into a compact
instruction table — one node per schema object, one op per keyword, pooled
strings, numbers and JSON literals, every `pattern` compiled to a Pike-VM
program — emitted as `cpp/src/qvf_schema_table.hpp` and
`js/src/schemaTable.ts`, beside the `KindInfo`/`RoleInfo` tables. A small
interpreter (`cpp/src/qvf_schema_check.hpp`, `js/src/schemaCheck.ts`) knows
nothing about QVF and evaluates `manifest.json` against the root schema
(`E_SCHEMA`) and each JSON member against its payload definition
(`E_PAYLOAD_SCHEMA`), reporting errors where `jsonschema`'s `iter_errors`
does (applicators descend; `oneOf`, `anyOf`, `not`, `contains` report at
their instance), so `where` is the same JSON pointer in every validator. The
supported subset is what the schema uses — `type`, `const`, `enum`,
`required`, `properties`, `patternProperties`, `additionalProperties`,
`dependentRequired`, `minProperties`/`maxProperties`, `prefixItems`, `items`,
`contains`/`minContains`/`maxContains`, `minItems`/`maxItems`,
`minimum`/`maximum`/`exclusiveMinimum`/`exclusiveMaximum`,
`minLength`/`maxLength` (code points), `pattern`, `allOf`/`anyOf`/`oneOf`/
`not`, `if`/`then`/`else`, `$ref` into `$defs`, boolean schemas — plus
`uniqueItems`, `multipleOf`, `propertyNames` and `dependentSchemas`; `format`
is an annotation, as in the Python reference. The generator refuses anything
else, a regular expression outside its ECMA-262 subset, and a `*Payload`
definition no `(kind, role)` row names, so a schema change cannot leave the
checkers silently behind. `python/tests/test_schema_checker_parity.py` proves
the C++ and JavaScript checkers agree with `jsonschema` on every corpus,
negative and legacy archive (same codes, same schema locations) and on
thousands of documents mutated from the corpus. No JSON Schema engine ships in
C++ or JavaScript. `ValidationOptions::schema_check = false` (JS
`schemaCheck: false`) skips the check and reports `W_SCHEMA_SKIPPED`.

Finding codes are the contract for bindings and the CLI; the Python validator
emits the same codes. Semantic checks are the ones in spec § 6 plus D15.

## 6. L3 — kinds (typed sections)

One struct per canonical kind, with `to_section(Writer&)` and
`from_section(const Reader&, const Section&)`. The plain ones (volumes, DOS,
EOS, phonons, citations, spectra 1-D) are generated from the schema; the ones
with semantics (`Structure`, `WavefunctionGTO`, `Trajectory`, `Bands`,
`RunRecord`, `JobSpec`, `ExcitedStates`) are hand-written and tested against
the corpus. Sketch:

```cpp
namespace qvf::kinds {
struct Structure {
    std::vector<Atom> atoms;                 // Å
    std::array<bool,3> pbc{};  std::optional<Lattice> lattice;
    std::optional<Biomolecule> bio;  std::vector<Bond> bonds;
    std::string id = "structure";
};
struct Volume {                               // any volume.* kind
    std::string kind;  GridSpec grid;  TensorView data;   // float32/64, 3-D (or 4-D for series/vector)
    std::optional<std::string> label, component, structure_ref;
};
struct WavefunctionGTO { … existing fields …;  std::optional<Ecp> ecp;
    std::optional<TensorView> overlap, density_matrix; };
struct Trajectory { Json metadata; TensorView coords; std::optional<TensorView> lattice,
    velocities, forces, times, energies; };
// … Bands, Dos, Spectrum, ExcitedStates, Forces, Hessian, AtomProperties, …
}
```

`TensorView` is non-owning (`pointer + dtype + shape + host byte order`); the
writer converts endianness and dtype in chunks while streaming, so a 1 GB
grid is never copied whole.

## 7. L4 — facades

### Writer (streaming)

```cpp
class Writer {
public:
    Writer(Source, std::unique_ptr<container::Sink>);        // file / memory / callback
    explicit Writer(Source, const std::string& path);
    void set_root(RootBlock, Json);                           // provenance, metadata, energies, …
    std::string add(const kinds::Structure&);                 // typed overloads for every kind
    std::string add(const kinds::Volume&);
    // Streaming: declare the member, then push slabs; sha256/CRC/size computed on the fly.
    MemberSink begin_binary(std::string section_id, std::string role, DType, std::vector<int64_t> shape);
    std::string add_section(Section);                         // escape hatch (vendor, custom)
    void finish();                                            // writes manifest.json last
};
```

Members are written to the ZIP as they are added (manifest last, per spec
§ 2.1), so memory is bounded by the largest single slab the producer pushes.
`QvfWriter` (today's API) becomes a thin compatibility layer over `Writer`.

### Reader (lazy)

```cpp
class Reader {
public:
    static Reader open(const std::string& path, Limits = {});
    static Reader from_bytes(std::vector<uint8_t>, Limits = {});
    const Manifest& manifest() const;
    std::vector<const Section*> sections_of(std::string_view kind) const;
    const Section* section(std::string_view id) const;
    Json      json(const Section&, std::string_view role, bool verify = true) const;
    Tensor    tensor(const Section&, std::string_view role, bool verify = true) const;   // owning
    void      tensor_into(const Section&, std::string_view role, TensorSink) const;       // streaming
    template <class K> K get(std::string_view id) const;      // kinds::Structure etc.
    validate::Report validate(validate::Level = validate::Level::All) const;
};
```

### Editor (in-place update)

```cpp
class Editor {                       // implements P19 and the job.spec lifecycle (§ 5.9)
public:
    static Editor open(const std::string& path);              // reads central directory only
    Writer& writer();                                         // appends new members to the same file
    void remove_section(std::string_view id);                 // orphans the members
    void set_run_status(RunStatus);
    void commit();                                            // rewrites manifest + central directory
    static void repack(const std::string& in, const std::string& out,
                       Method = Method::Deflate);              // drops orphans, recompresses
};
```

Appending to an existing ZIP is standard: truncate at the old central
directory, append entries, write a new central directory. `commit` is atomic
on POSIX via write-to-temp + rename when the caller asks for it.

## 8. L5 — bindings

### C ABI (complete)

WP-44 implements this surface as additive ABI 1, documented in
[c_api.md](c_api.md). Canonical writers use schema-checked section descriptors
(JSON payloads plus typed pointer/shape arrays) instead of a distinct long
argument list per kind. New constructors return status with output handles;
existing signatures and per-writer error access remain compatible.

Principles: one opaque handle per object (`qvf_writer`, `qvf_reader`,
`qvf_editor`, `qvf_report`), every function returns `qvf_status` (0 = ok),
`qvf_last_error(handle)`; arrays are pointers + explicit dtype/shape; nested
payloads are JSON strings; `qvf_abi_version()` returns a monotonically
increasing integer and `qvf_toolkit_version()` the semver string; all symbols
exported through a `QVF_API` macro (shared builds on Windows and ELF
visibility). Surface (names indicative):

```
writer:   qvf_writer_create / destroy / finish / set_root_json
          qvf_writer_add_structure(w, n_atoms, symbols, Z, positions, pbc[3], lattice_rows[9] | NULL, bonds…)
          qvf_writer_add_volume(w, kind, grid_json, dtype, shape[3], data, label)
          qvf_writer_begin_binary(w, section_id, role, dtype, rank, shape) / qvf_writer_write(w, bytes, n) / qvf_writer_end_binary(w)
          qvf_writer_add_wavefunction_gto(w, …restricted…) / _unrestricted(…) / _set_ecp_json / _set_overlap(...)
          qvf_writer_add_trajectory(w, metadata_json, n_frames, n_atoms, coords, lattice | NULL, velocities | NULL, forces | NULL, times | NULL)
          qvf_writer_add_vibrations / add_spectrum_xy / add_spectrum_json / add_excited_states_json
          qvf_writer_add_bands(w, kpath_json, n_spin, n_k, n_bands, eig, n_channels, proj | NULL)
          qvf_writer_add_dos_total / add_dos_projected / add_dos_coop_cohp / add_fermi_surface / add_phonon_bands / add_phonon_dos / add_eos
          qvf_writer_add_atom_properties(w, role, n_atoms, values) / add_forces / add_hessian / add_bond_orders_json / add_scf_history_json
          qvf_writer_add_citations / add_run_record / add_job_spec_json / add_vendor_json / add_vendor_binary / set_extensions_json
reader:   qvf_reader_open / from_bytes / destroy
          qvf_reader_manifest_json(r, &str, &len)   qvf_reader_section_count / section_id / section_kind
          qvf_reader_member_info(r, section_id, role, &dtype, &rank, shape[], &nbytes)
          qvf_reader_read_json(r, section_id, role, &str, &len)   qvf_reader_read_binary(r, section_id, role, buffer, nbytes)
          qvf_reader_read_binary_chunked(r, section_id, role, callback, user)
validate: qvf_validate_path(path, level, &report) / qvf_report_count / qvf_report_get(i, &severity, &code, &where, &message) / qvf_report_destroy
editor:   qvf_editor_open / writer / set_run_status / commit / destroy    qvf_repack(in, out)
```

### Fortran module

`fortran/qvf.f90` provides `module qvf` with `bind(C)` interfaces for every
C ABI entry plus Fortran-friendly wrappers that take `character(len=*)`,
`real(c_double), dimension(:,:)` arrays (transposed to C row-major by the
wrapper where the spec's layout requires it), and `integer` shapes. Built as
part of the CMake project when a Fortran compiler is found (`QVF_FORTRAN=ON`),
tested by a gfortran example that writes a structure, a grid, a trajectory,
and a wavefunction and is validated by the Python reader in CI.

### Python

Distribution `qvf-writer` (maintainer decision 2026-10-15), package `qvf`.
WP-46 supplies the package, streaming writer, lazy typed reader, editor,
validation, model and unit modules. WP-47 supplies the CLI and cumulative
Python validation stages; converters remain WP-48. See the
[command-line guide](cli.md) for the implemented commands. Kind views expose typed member metadata for every kind,
with named convenience accessors for structure, volume and GTO sections:


```
qvf/
  __init__.py        # Writer, Reader, Editor, validate, __version__
  _container.py      # zipfile wrappers, limits, ZIP64 streaming (stdlib)
  model.py           # dataclasses: Manifest, Section, MemberSpec, Lattice, GridSpec
  _kinds.py          # GENERATED tables from the schema
  kinds/             # typed sections (Structure, Volume, WavefunctionGTO, …) with to/from
  writer.py          # streaming writer (writes members as added; manifest last)
  reader.py          # lazy reader; typed accessors; numpy optional (memoryview path without it)
  editor.py          # in-place update, repack
  operations.py      # explicit-output extraction, diff, repack, scrub
  conformance.py     # packaged corpus runner, shared with the old script
  validate.py        # finding codes identical to C++; jsonschema optional
  units.py           # vocabulary + conversions
  convert/           # cube, xyz/extxyz, molden, pdb, xsf, chgcar/poscar, cif
  cli.py             # `qvf` entry point
  schema/            # bundled qvf_manifest.schema.json (package data, not data-files)
qvf_writer.py, qvf_reader.py   # shims: re-export, DeprecationWarning, kept through 1.x
```

The Python implementation stays pure (stdlib + optional numpy); an optional
`qvf-native` wheel exposing the C ABI through `cffi` can come later for
speed, behind the same API.

### JavaScript / TypeScript (reader first)

`@qvf/core`: pure TypeScript reader + validator over `fflate` (inflate) and
WebCrypto (sha256), ESM, browser and Node, typed accessors returning
`Float32Array`/`Float64Array` views with shape. Certified by running the
conformance corpus in CI. A `Writer` with the typed builders of the Python
writer (structure, trajectory, wavefunction.gto, volume.*, excited_states,
spectra.*, reaction.scheme) writes archives in memory or to a streaming sink
(see `js/README.md`).

### Rust

`qvf` crate: independent implementation (`zip` + `serde_json` + `sha2`),
reader, writer, validator, `qvf-cli` optional; certified by the corpus. A
`qvf-sys` binding to the C ABI is *not* planned — Rust users prefer a native
crate and the corpus makes that safe.

## 9. Build, packaging, versioning

- **C++:** CMake ≥ 3.16; options `QVF_BUILD_SHARED`, `QVF_BUILD_TESTS`,
  `QVF_BUILD_EXAMPLES`, `QVF_FORTRAN`, `QVF_FUZZ`; `qvfConfig.cmake`,
  `qvf.pc`; single-header amalgamation kept for writer+reader; vcpkg and
  conan recipes once the ABI is complete.
- **Versioning:** toolkit semver (`VERSION`); C ABI integer (`QVF_ABI_VERSION`)
  bumped only on incompatible signature change (additions do not bump);
  shared-library SONAME follows the ABI integer; the format version is
  independent (`qvf_version` in the manifest).
- **Python:** tagged-source installation is available; PyPI wheel publication
  and trusted publishing remain WP-49 work. Add the `numpy` extra for numeric
  writer conveniences and NumPy views. The distribution keeps its existing
  name and contains the `qvf` package plus compatibility top-level imports.
- **Docs:** API reference generated (Doxygen for C++/C, mkdocstrings for
  Python), a "kind catalogue" page generated from the schema, and the three
  existing guides updated to the new entry points.

## 10. Testing strategy

| Layer | Tests |
|---|---|
| L0 | known-answer (CRC, SHA, inflate vectors), ZIP64 synthetic archives, data-descriptor archives, fuzz targets (zip, inflate, json, manifest) with libFuzzer |
| L1–L2 | generated checker ≡ `jsonschema` on corpus + negative corpus; finding-code parity C++ ≡ Python |
| L3 | round-trip property tests per kind (write → read → equal), normalization self-check via `overlap` member |
| L4 | streaming vs in-memory produce byte-identical members; editor append → validate → repack |
| L5 | each binding runs the conformance corpus; Fortran and C examples validated by the Python reader |
| CI | Linux gcc/clang (+ASan/UBSan), macOS, Windows MSVC; Python 3.9–3.14; optional big-endian job via QEMU (s390x) proving the writers byte-swap |

## 11. Migration of the current API

| Today | Target | Transition |
|---|---|---|
| `qvf::QvfWriter` | `qvf::Writer` | `QvfWriter` kept as subclass/alias through 1.x |
| `qvf::QvfReader` | `qvf::Reader` | same |
| `qvf_c.h` (16 functions) | complete C ABI | existing symbols unchanged; new ones added; `QVF_ABI_VERSION` 2 |
| `qvf_writer.py`, `qvf_reader.py` | `qvf` package | shims with `DeprecationWarning` |
| `qvf-validate` | `qvf validate` | alias kept |
| `cpp/qvf_single.hpp` | regenerated, now includes L0–L4 | unchanged usage |
