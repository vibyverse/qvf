// qvf/qvf_reader.hpp — reference C++ reader for QVF archives.
//
// Completes the round-trip: open a .qvf, parse manifest.json, verify member
// sha256s, and decode JSON / binary members. Zero external dependencies — the
// ZIP reader, INFLATE decompressor, JSON parser, and SHA-256 are all in the
// toolkit. License: Apache-2.0.
//
// The validator reports findings with the stable codes of the specification's
// Appendix B ("Validation finding codes") — the same codes the Python reference
// validator (python/qvf_reader.py) emits, including the JSON Schema findings:
// a checker generated from the schema (scripts/gen_kinds.py) checks
// manifest.json and every JSON member payload without a schema engine.
#ifndef QVF_QVF_READER_HPP
#define QVF_QVF_READER_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "qvf/container.hpp"
#include "qvf/json.hpp"
#include "qvf/model.hpp"

namespace qvf {

// A decoded binary member: little-endian raw bytes + declared dtype/shape.
struct BinaryMember {
    std::string dtype;            // e.g. "float64", "float32", "int32"
    std::vector<int64_t> shape;
    std::vector<uint8_t> bytes;   // C-contiguous, little-endian

    size_t element_count() const {
        size_t n = 1;
        for (int64_t d : shape) n *= static_cast<size_t>(d);
        return n;
    }
    // Convenience: reinterpret the bytes as a vector<double>. Only valid for
    // float64 members (throws otherwise).
    std::vector<double> as_doubles() const;
};

// Resource caps applied to untrusted archives (spec § 6): max_member_bytes,
// max_total_bytes, max_ratio, ratio_exempt_bytes, max_json_bytes and
// max_json_depth, with the reference values as defaults. The type lives in
// the container layer (qvf/container.hpp), whose ZipReader enforces it.
using Limits = container::Limits;

// Options of a validation run. `schema_check` runs the JSON Schema check of
// manifest.json (E_SCHEMA) and of every JSON member with a payload definition
// (E_PAYLOAD_SCHEMA) through the checker generated from the schema; false
// skips it and reports W_SCHEMA_SKIPPED instead, the behaviour of a validator
// without a schema engine.
struct ValidationOptions {
    bool schema_check = true;
};

// One validation finding. `code` is a stable identifier from spec Appendix B
// (E_* for errors, W_* for warnings). `where` is "archive", "entry:<zip name>",
// "manifest.json", or a JSON Pointer into manifest.json such as
// "/sections/2/members/data".
struct Finding {
    enum class Severity { Error, Warning } severity;
    std::string code, where, message;
};

// True if any finding is an error (warnings leave an archive valid).
bool has_errors(const std::vector<Finding>& findings);

class QvfReader {
public:
    // Open a .qvf file from disk, or from an in-memory archive. Throws on an
    // unreadable ZIP, on duplicate ZIP entry names (an ambiguous archive, spec
    // § 2.1; several manifest.json entries are an archive updated in place,
    // whose last manifest.json is read, § 2.3), when manifest.json is missing,
    // over `limits.max_json_bytes`, nested deeper than `limits.max_json_depth`,
    // or not valid JSON, and on a `byte_order` other than "little" (§ 2.2).
    // Member reads enforce every cap of `limits` (member and total size,
    // ratio; JSON members also the JSON caps), and refuse a member whose
    // declared `size` disagrees with its ZIP entry (E_BINARY_SIZE, § 4.2).
    // open() keeps the file open and reads entries from it on demand, so the
    // archive need not fit in memory.
    static QvfReader open(const std::string& path, const Limits& limits = Limits{});
    static QvfReader from_bytes(std::vector<uint8_t> archive,
                                const Limits& limits = Limits{});

    QvfReader(QvfReader&&) noexcept;
    QvfReader& operator=(QvfReader&&) noexcept;
    ~QvfReader();

    // The parsed manifest, and its `sections` array.
    const Json& manifest() const;
    const Json& sections() const;

    // The manifest as the typed, lossless model of qvf/model.hpp (built when
    // the archive is opened). Throws std::runtime_error when manifest.json is
    // not a JSON object.
    const Manifest& typed_manifest() const;

    bool has_section(const std::string& id) const;
    const Json& section(const std::string& id) const;   // throws if missing

    // Read a member (resolved through the manifest by section id + role),
    // verifying its sha256 by default.
    std::vector<uint8_t> read_member_bytes(const std::string& section_id,
                                           const std::string& role,
                                           bool verify = true) const;
    Json read_json_member(const std::string& section_id, const std::string& role,
                          bool verify = true) const;
    BinaryMember read_binary_member(const std::string& section_id,
                                    const std::string& role,
                                    bool verify = true) const;

    // Full validation (spec § 6) with stable finding codes: container caps and
    // path/duplicate/orphan rules (superseded manifest.json copies are
    // orphans, spec § 2.3), manifest and section well-formedness, member
    // existence, declared member size, sha256, JSON size/depth/parse, dtype
    // and binary size, UTF-8,
    // reference resolution, the structure periodicity contract, bands shapes,
    // the trajectory and reaction.path members (per-frame extents, coords
    // precision, binary energies, lattice and pbc), forces and hessian
    // shapes, excited_states state indices and orbital rows, the
    // wavefunction.gto rules of spec § 5.3 (including the Appendix A.1
    // normalization check of an `overlap` member, E_NORMALIZATION), the data
    // shapes of the vector, complex, series and Laplacian volumes, a series
    // axis against its slabs, the span of a periodic grid against its
    // structure's cell (spec § 5.2), a grid units.data against the unit spec
    // § 4.3 fixes for its volume kind (E_UNIT_CONTRADICTED), wavefunction.bloch shapes, k-point
    // weights and periodic structure, the reaction.scheme participants, step
    // ids and repeated structures (spec § 5.14), the bytes of preview members
    // against their media_type (spec § 5.15), a run.record job_ref (spec
    // § 5.8), a structure section's charge and multiplicity against its atoms
    // (E_CHARGE_MULTIPLICITY) and, with one structure section, against
    // provenance (spec § 5.1), the root energies consistency rule, unit names
    // (W_UNIT_UNKNOWN, spec § 4.3), and the JSON Schema (E_SCHEMA,
    // E_PAYLOAD_SCHEMA) unless `options.schema_check` is false (then
    // W_SCHEMA_SKIPPED). The withdrawn qvf_version 2 is checked as 1 and
    // reported as the warning W_VERSION_WITHDRAWN (spec § 3.1).
    std::vector<Finding> validate_findings(const Limits& limits = Limits{},
                                           const ValidationOptions& options = {}) const;

    // Validate an archive that may not even open as a QvfReader (unreadable
    // ZIP, duplicate entries, broken manifest): such failures are returned as
    // findings rather than thrown.
    static std::vector<Finding> validate_file(const std::string& path,
                                              const Limits& limits = Limits{},
                                              const ValidationOptions& options = {});
    static std::vector<Finding> validate_bytes(std::vector<uint8_t> archive,
                                               const Limits& limits = Limits{},
                                               const ValidationOptions& options = {});

    // Legacy form of validate_findings() with the default limits. Returns true
    // if there is no error; otherwise appends "CODE: message" strings for the
    // errors (warnings are not reported here).
    bool validate(std::vector<std::string>& errors) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    std::shared_ptr<const Manifest> typed_;   // null when manifest.json is not an object
    explicit QvfReader(std::unique_ptr<Impl> impl);
};

}  // namespace qvf

#endif  // QVF_QVF_READER_HPP
