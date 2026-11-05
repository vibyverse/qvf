// qvf_reader.cpp — implementation of the QVF C++ reader. License: Apache-2.0.
#include "qvf/qvf_reader.hpp"
#include "qvf/validation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>

#include "qvf/kinds_generated.hpp"
#include "qvf/units.hpp"
#include "qvf_schema_check.hpp"
#include "qvf_sha256.hpp"
#include "qvf_unzip.hpp"

namespace qvf {

static size_t dtype_itemsize(const std::string& d) {
    if (d == "int8" || d == "uint8") return 1;
    if (d == "int16" || d == "uint16") return 2;
    if (d == "int32" || d == "uint32" || d == "float32") return 4;
    if (d == "int64" || d == "uint64" || d == "float64") return 8;
    return 0;
}

// --- validation helpers (spec § 6, Appendix B) ------------------------------
namespace detail {
namespace vld {

// Canonical kinds (the schema's Section.oneOf branches), from the tables
// generated from the schema (qvf/kinds_generated.hpp); a vendor kind
// (x_<...>) is always accepted.
inline bool is_canonical_kind(const std::string& k) { return gen::is_canonical_kind(k); }

// A one-dimensional spectrum kind: its `spectrum` member is a SpectrumPayload.
inline bool is_spectrum_1d(const std::string& kind) {
    return gen::payload_def(kind, "spectrum") == "SpectrumPayload";
}

inline bool starts_with(const std::string& s, const char* prefix) {
    return s.compare(0, std::strlen(prefix), prefix) == 0;
}

inline std::string quote_str(const std::string& s) { return "'" + s + "'"; }

inline std::string format_number(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%g", v);
    return buf;
}

// Compact rendering of a JSON value for messages.
inline std::string show(const Json& j) {
    std::string s;
    try {
        s = j.dump(0);
    } catch (const std::exception&) {
        return "<value>";
    }
    if (s.size() > 80) s = s.substr(0, 77) + "...";
    return s;
}

// The value of `key` when `obj` is an object that has it (null included).
inline const Json* get(const Json& obj, const std::string& key) {
    if (!obj.is_object() || !obj.contains(key)) return nullptr;
    return &obj.at(key);
}

// Python-style truthiness, so both reference validators agree on loosely
// typed flags (`critical`) and non-empty checks (`source.program`).
inline bool truthy(const Json& j) {
    switch (j.type()) {
        case Json::Type::Null: return false;
        case Json::Type::Bool: return j.as_bool();
        case Json::Type::Int:
        case Json::Type::Double: return j.as_double() != 0.0;
        case Json::Type::String:
        case Json::Type::Raw: return !j.as_string().empty();
        case Json::Type::Array:
        case Json::Type::Object: return j.size() != 0;
    }
    return false;
}

// JSON Pointer (RFC 6901) pieces into manifest.json.
inline std::string escape_pointer(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '~') out += "~0";
        else if (c == '/') out += "~1";
        else out.push_back(c);
    }
    return out;
}
inline std::string section_pointer(size_t i) {
    return "/sections/" + std::to_string(i);
}
inline std::string member_pointer(size_t i, const std::string& role) {
    return section_pointer(i) + "/members/" + escape_pointer(role);
}

// Why a member path / ZIP entry name is illegal (spec § 2.1); empty if legal.
inline std::string path_problem(const std::string& p) {
    if (p.find('\\') != std::string::npos) return "contains a backslash";
    if (p.find('\0') != std::string::npos) return "contains a NUL character";
    if (!p.empty() && p[0] == '/') return "begins with '/'";
    if (p.size() >= 2 && p[1] == ':' &&
        ((p[0] >= 'A' && p[0] <= 'Z') || (p[0] >= 'a' && p[0] <= 'z')))
        return "begins with a drive letter";
    size_t start = 0;
    while (true) {
        size_t slash = p.find('/', start);
        size_t len = slash == std::string::npos ? std::string::npos : slash - start;
        if (p.compare(start, len, "..") == 0) return "contains a '..' segment";
        if (slash == std::string::npos) break;
        start = slash + 1;
    }
    return "";
}

// Strict UTF-8 (no overlongs, no surrogates, nothing above U+10FFFF), the
// same acceptance as Python's bytes.decode("utf-8").
inline bool valid_utf8(const std::vector<uint8_t>& s) {
    size_t i = 0, n = s.size();
    while (i < n) {
        uint8_t c = s[i];
        if (c < 0x80) { ++i; continue; }
        size_t len;
        uint32_t cp, min;
        if ((c & 0xE0) == 0xC0) { len = 2; cp = c & 0x1Fu; min = 0x80; }
        else if ((c & 0xF0) == 0xE0) { len = 3; cp = c & 0x0Fu; min = 0x800; }
        else if ((c & 0xF8) == 0xF0) { len = 4; cp = c & 0x07u; min = 0x10000; }
        else return false;
        if (len > n - i) return false;
        for (size_t k = 1; k < len; ++k) {
            if ((s[i + k] & 0xC0) != 0x80) return false;
            cp = (cp << 6) | (s[i + k] & 0x3Fu);
        }
        if (cp < min || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF))
            return false;
        i += len;
    }
    return true;
}

class Sink {
public:
    explicit Sink(std::vector<Finding>& out, bool references = true)
        : out_(out), references_(references) {}
    void error(const char* code, std::string where, std::string message) {
        // A kind checker can report both local shape and reference errors.
        // Resolve-reference findings belong to the last cumulative stage.
        if (!references_ && std::strcmp(code, "E_REF_UNRESOLVED") == 0) return;
        out_.push_back(Finding{Finding::Severity::Error, code, std::move(where),
                               std::move(message)});
    }
    void warning(const char* code, std::string where, std::string message) {
        out_.push_back(Finding{Finding::Severity::Warning, code, std::move(where),
                               std::move(message)});
    }

private:
    std::vector<Finding>& out_;
    bool references_;
};

// Semantic checks on a `structure` payload (spec § 5.1); mirrors
// _structure_payload_errors() in python/qvf_reader.py.
inline std::vector<std::string> structure_payload_errors(const std::string& sid,
                                                         const Json& payload) {
    std::vector<std::string> errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    if (!payload.is_object()) {
        errs.push_back(sec + "structure payload must be a JSON object");
        return errs;
    }
    auto non_null = [&](const char* key) -> const Json* {
        const Json* v = get(payload, key);
        return v && !v->is_null() ? v : nullptr;
    };
    const Json* pbc = non_null("pbc");
    const Json* lattice = non_null("lattice_vectors");
    const Json* dim = non_null("dimensionality");

    if (dim && !dim->is_number()) {
        errs.push_back(sec + "dimensionality must be an integer, got " + show(*dim) +
                       " (spec § 5.1)");
        dim = nullptr;
    }
    if (!pbc) {
        if (lattice)
            errs.push_back(sec + "structure payload has lattice_vectors but no "
                                 "'pbc'; pbc is REQUIRED for periodic structures "
                                 "(spec § 5.1)");
        if (dim && dim->as_double() != 0.0)
            errs.push_back(sec + "dimensionality " + show(*dim) +
                           " without 'pbc'; absent pbc means a molecule, so "
                           "dimensionality must be 0 (spec § 5.1)");
        return errs;
    }
    bool three_bools = pbc->is_array() && pbc->size() == 3;
    int n_periodic = 0;
    if (three_bools) {
        for (const Json& x : pbc->elements()) {
            if (!x.is_bool()) { three_bools = false; break; }
            n_periodic += x.as_bool() ? 1 : 0;
        }
    }
    if (!three_bools) {
        errs.push_back(sec + "pbc must be three booleans, got " + show(*pbc));
        return errs;
    }
    if (n_periodic && !(lattice && lattice->is_array()))
        errs.push_back(sec + "pbc marks " + std::to_string(n_periodic) +
                       " periodic axis/axes but lattice_vectors is " +
                       (lattice ? "not an array" : "null or absent"));
    if (dim && dim->as_double() != static_cast<double>(n_periodic))
        errs.push_back(sec + "dimensionality " + show(*dim) + " != sum(pbc) " +
                       std::to_string(n_periodic) + " (spec § 5.1 invariant)");
    return errs;
}

// Every declared unit of a member payload (spec § 4.3, "Declared units"),
// as (location inside the payload, value): each value of its `units` object,
// for a one-dimensional spectrum x_units, y_units and broadening.units, the
// `units` string of a volume.series `axis` member, and for excited_states
// each value of meta.units. Mirrors _declared_units() in python/qvf_reader.py.
inline std::vector<std::pair<std::string, const Json*>> declared_units(
    const std::string& kind, const std::string& role, const Json& payload) {
    std::vector<std::pair<std::string, const Json*>> out;
    if (!payload.is_object()) return out;
    if (const Json* units = get(payload, "units")) {
        if (units->is_object())
            for (const auto& kv : units->items())
                out.emplace_back("units/" + kv.first, &kv.second);
        else if (kind == "volume.series" && role == "axis")
            out.emplace_back("units", units);
    }
    if (kind == "excited_states")
        if (const Json* meta = get(payload, "meta"))
            if (const Json* units = get(*meta, "units"))
                if (units->is_object())
                    for (const auto& kv : units->items())
                        out.emplace_back("meta/units/" + kv.first, &kv.second);
    if (!is_spectrum_1d(kind)) return out;
    for (const char* key : {"x_units", "y_units"})
        if (const Json* v = get(payload, key)) out.emplace_back(key, v);
    if (const Json* b = get(payload, "broadening"))
        if (const Json* v = get(*b, "units")) out.emplace_back("broadening/units", v);
    return out;
}

// The equal-length rule of a 1-D `spectrum` payload (spec § 5.4):
// `frequencies`, `intensities` and, when present, `x` have one shared length.
// JSON Schema cannot say so, so this runs with or without the schema check.
// Mirrors _spectrum_payload_errors() in python/qvf_reader.py.
inline std::vector<std::string> spectrum_payload_errors(const std::string& sid,
                                                        const Json* payload) {
    std::vector<std::string> errs;
    if (!payload || !payload->is_object()) return errs;
    const Json* freqs = schema::find_last(*payload, "frequencies");
    if (!freqs || !freqs->is_array()) return errs;
    for (const char* key : {"intensities", "x"}) {
        const Json* other = schema::find_last(*payload, key);
        if (other && other->is_array() && other->size() != freqs->size())
            errs.push_back("section " + quote_str(sid) + ": spectrum " + key + " has " +
                           std::to_string(other->size()) + " entries but frequencies has " +
                           std::to_string(freqs->size()) +
                           "; they must have the same length (spec § 5.4)");
    }
    return errs;
}

// A member spec's `shape` when it is an array of non-negative integers.
inline bool shape_of(const Json* member, std::vector<int64_t>& out) {
    out.clear();
    const Json* shape = member ? get(*member, "shape") : nullptr;
    if (!shape || !shape->is_array()) return false;
    for (const Json& d : shape->elements()) {
        if (d.type() != Json::Type::Int || d.as_int() < 0) return false;
        out.push_back(d.as_int());
    }
    return true;
}

// "[a, b, c]" for dims[from, to).
inline std::string dims_str(const std::vector<int64_t>& dims, size_t from = 0,
                            size_t to = static_cast<size_t>(-1)) {
    std::string s = "[";
    for (size_t k = from; k < dims.size() && k < to; ++k) {
        if (k > from) s += ", ";
        s += std::to_string(dims[k]);
    }
    return s + "]";
}

// Spin-axis agreement of a `bands` section's binary members (spec § 5.6),
// checked from the manifest alone; mirrors _bands_shape_errors() in
// python/qvf_reader.py. Returns (role, message) pairs.
inline std::vector<std::pair<std::string, std::string>> bands_shape_errors(
    const std::string& sid, const Json& members) {
    std::vector<std::pair<std::string, std::string>> errs;
    std::vector<int64_t> eig, proj, occ;
    if (!shape_of(get(members, "eigenvalues"), eig) || eig.size() != 3) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    if (shape_of(get(members, "projections"), proj)) {
        if (proj.size() == 4 && !std::equal(eig.begin(), eig.end(), proj.begin()))
            errs.emplace_back("projections",
                              sec + "projections shape " + dims_str(proj) +
                                  " has [n_spin, n_k, n_bands] = " + dims_str(proj, 0, 3) +
                                  ", but eigenvalues is " + dims_str(eig) + " (spec § 5.6)");
        else if (proj.size() == 3 && (proj[0] != eig[1] || proj[1] != eig[2]))
            errs.emplace_back("projections",
                              sec + "rank-3 projections shape " + dims_str(proj) +
                                  " has [n_k, n_bands] = " + dims_str(proj, 0, 2) +
                                  ", but eigenvalues has " + dims_str(eig, 1) +
                                  " (spec § 5.6)");
    }
    if (shape_of(get(members, "occupations"), occ) && occ.size() == 3 && occ != eig)
        errs.emplace_back("occupations",
                          sec + "occupations shape " + dims_str(occ) +
                              " differs from eigenvalues " + dims_str(eig) +
                              "; both are [n_spin, n_k, n_bands] (spec § 5.6)");
    return errs;
}

// --- trajectory / reaction.path (spec § 5.5) ----------------------------------

// Largest allowed difference, in hartree, between the binary `energies` member
// of a trajectory or reaction path and its metadata.energies.
constexpr double kTrajectoryEnergyToleranceEh = 1e-10;

inline bool is_trajectory_kind(const std::string& kind) {
    return kind == "trajectory" || kind == "reaction.path";
}

// The values of a rank-1 float64 member whose byte length matches its shape,
// decoded little-endian; false for anything else.
inline bool float64_values(const Json* member, const std::vector<uint8_t>* bytes,
                           std::vector<double>& out) {
    std::vector<int64_t> shape;
    const Json* dtype = member ? get(*member, "dtype") : nullptr;
    if (!bytes || !dtype || !dtype->is_string() || dtype->as_string() != "float64" ||
        !shape_of(member, shape) || shape.size() != 1 ||
        bytes->size() != 8 * static_cast<uint64_t>(shape[0]))
        return false;
    out.assign(static_cast<size_t>(shape[0]), 0.0);
    for (size_t k = 0; k < out.size(); ++k) {
        uint64_t bits = 0;
        for (int b = 7; b >= 0; --b) bits = (bits << 8) | (*bytes)[8 * k + b];
        std::memcpy(&out[k], &bits, 8);
    }
    return true;
}

// Shape and consistency rules of a `trajectory` or `reaction.path` section:
// coords float32 only with metadata.precision "float32" (and precision, when
// present, naming the dtype); against a [n_frames, n_atoms, 3] coords (another
// rank is left to the schema), velocities and gradient [n_frames, n_atoms, 3]; times and energies (and, for a trajectory,
// cell_volumes and temperatures) [n_frames]; lattice [3, 3] or
// [n_frames, 3, 3], each kind's lattice read in its own convention (§ 4.3)
// but checked here by shape only; the binary energies equal to
// metadata.energies within 1e-10 hartree; a trajectory lattice needs
// metadata.pbc and a true pbc entry needs the lattice. Rules that need the
// metadata are skipped when it could not be read. Mirrors _trajectory_errors()
// in python/qvf_reader.py. Returns (role, message) pairs.
inline std::vector<std::pair<std::string, std::string>> trajectory_errors(
    const std::string& sid, const std::string& kind, const Json& members,
    const Json* metadata, const std::vector<uint8_t>* energies_bytes) {
    std::vector<std::pair<std::string, std::string>> errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    const Json* meta = metadata && metadata->is_object() ? metadata : nullptr;
    const Json* coords_spec = get(members, "coords");
    const Json* dtype = coords_spec ? get(*coords_spec, "dtype") : nullptr;
    if (meta && dtype && dtype->is_string() &&
        (dtype->as_string() == "float32" || dtype->as_string() == "float64")) {
        const Json* p = get(*meta, "precision");
        const bool f32 = p && p->is_string() && p->as_string() == "float32";
        if (dtype->as_string() == "float32" && !f32)
            errs.emplace_back("coords", sec + "coords are float32 but metadata.precision is " +
                                            (p ? show(*p) : std::string("absent")) +
                                            "; float32 coords require precision 'float32' "
                                            "(spec § 5.5)");
        else if (dtype->as_string() == "float64" && f32)
            errs.emplace_back("coords", sec + "metadata.precision is 'float32' but coords are "
                                            "float64; precision must name the dtype of "
                                            "coords (spec § 5.5)");
    }
    // A coords shape of another rank is left to the schema; the extents below
    // are only compared against a [n_frames, n_atoms, 3] coords.
    std::vector<int64_t> coords, shape;
    if (shape_of(coords_spec, coords) && coords.size() == 3 && coords[2] == 3) {
        const std::vector<int64_t> vectors = {coords[0], coords[1], 3};
        const std::vector<int64_t> scalars = {coords[0]};
        for (const char* role : {"velocities", "gradient"})
            if (shape_of(get(members, role), shape) && shape != vectors)
                errs.emplace_back(role, sec + role + " shape " + dims_str(shape) +
                                            " is not [n_frames, n_atoms, 3] = " +
                                            dims_str(vectors) +
                                            ", the extents of coords (spec § 5.5)");
        const bool trajectory = kind == "trajectory";
        for (const char* role : {"times", "energies", "cell_volumes", "temperatures"}) {
            const std::string r = role;
            if (!trajectory && (r == "cell_volumes" || r == "temperatures")) continue;
            if (shape_of(get(members, role), shape) && shape != scalars)
                errs.emplace_back(role, sec + role + " shape " + dims_str(shape) +
                                            " is not [n_frames] = " + dims_str(scalars) +
                                            ", the frames of coords (spec § 5.5)");
        }
        const std::vector<int64_t> fixed = {3, 3};
        const std::vector<int64_t> per_frame = {coords[0], 3, 3};
        if (shape_of(get(members, "lattice"), shape) && shape != fixed && shape != per_frame)
            errs.emplace_back("lattice", sec + "lattice shape " + dims_str(shape) +
                                             " is neither [3, 3] nor [n_frames, 3, 3] = " +
                                             dims_str(per_frame) + " (spec § 5.5)");
    }
    if (meta) {
        std::vector<double> binary;
        const Json* listed = get(*meta, "energies");
        bool numbers = listed && listed->is_array();
        if (numbers)
            for (const Json& x : listed->elements()) numbers = numbers && x.is_number();
        if (numbers && float64_values(get(members, "energies"), energies_bytes, binary)) {
            if (listed->size() != binary.size()) {
                errs.emplace_back("energies",
                                  sec + "the energies member has " +
                                      std::to_string(binary.size()) +
                                      " values but metadata.energies has " +
                                      std::to_string(listed->size()) +
                                      "; they must agree (spec § 5.5)");
            } else {
                size_t bad = 0, first = 0;
                for (size_t k = 0; k < binary.size(); ++k)
                    if (!(std::fabs(binary[k] - listed->at(k).as_double()) <=
                          kTrajectoryEnergyToleranceEh)) {
                        if (!bad) first = k;
                        ++bad;
                    }
                if (bad)
                    errs.emplace_back(
                        "energies",
                        sec + "the energies member disagrees with metadata.energies at " +
                            std::to_string(bad) + " frame(s), first frame " +
                            std::to_string(first) + ": " + format_number(binary[first]) +
                            " vs " + format_number(listed->at(first).as_double()) +
                            " hartree; they must agree within 1e-10 hartree (spec § 5.5)");
            }
        }
    }
    if (kind == "trajectory" && meta) {
        const Json* pbc = get(*meta, "pbc");
        const bool has_lattice = members.contains("lattice");
        bool periodic = false;
        if (pbc && pbc->is_array())
            for (const Json& x : pbc->elements())
                if (x.is_bool() && x.as_bool()) periodic = true;
        if (has_lattice && (!pbc || pbc->is_null()))
            errs.emplace_back("metadata", sec + "the section has a lattice member but its "
                                                "metadata has no pbc; pbc is REQUIRED when "
                                                "lattice is present (spec § 5.5)");
        else if (periodic && !has_lattice)
            errs.emplace_back("metadata", sec + "metadata.pbc " + show(*pbc) +
                                              " marks a periodic axis but the section has "
                                              "no lattice member (spec § 5.5)");
    }
    return errs;
}

// --- root energies (spec § 3.6) ---------------------------------------------

// Largest allowed difference, in hartree, between two root statements of one
// energy.
constexpr double kEnergyToleranceEh = 1e-8;

inline std::string ascii_lower(std::string s) {
    for (char& c : s)
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    return s;
}

// Hartree per unit for the energy units the root consistency check converts
// (CODATA 2018); false for any other unit. Mirrors _HARTREE_PER_UNIT in
// python/qvf_reader.py.
inline bool hartree_per_unit(const std::string& unit, double& factor) {
    const std::string u = ascii_lower(unit);
    if (u == "hartree" || u == "eh" || u == "ha") factor = 1.0;
    else if (u == "ev") factor = 1.0 / 27.211386245988;
    else if (u == "kcal_per_mol") factor = 1.0 / 627.5094740631;
    else if (u == "kj_per_mol") factor = 1.0 / 2625.4996394799;
    else return false;
    return true;
}

// A {value, units} energy in hartree; false when it is malformed or its unit
// is not one the check converts.
inline bool energy_eh(const Json* energy, double& out) {
    if (!energy || !energy->is_object()) return false;
    const Json* value = get(*energy, "value");
    const Json* units = get(*energy, "units");
    double factor = 0.0;
    if (!value || !value->is_number() || !units || !units->is_string() ||
        !hartree_per_unit(units->as_string(), factor))
        return false;
    out = value->as_double() * factor;
    return true;
}

inline std::string fixed10(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.10f", v);
    return buf;
}

// Root blocks that state one energy twice must agree: every `energies` total
// stated for the method provenance.method names against provenance.scf_energy,
// and energies.components.zpve against thermochemistry.zpve_eh. Mirrors
// _root_consistency_errors() in python/qvf_reader.py. Returns (where, message).
inline std::vector<std::pair<std::string, std::string>> root_consistency_errors(
    const Json& manifest) {
    std::vector<std::pair<std::string, std::string>> errs;
    const Json* energies = get(manifest, "energies");
    if (!energies || !energies->is_object()) return errs;
    auto compare = [&](const std::string& where, const std::string& label,
                       const Json* got_json, double ref, const std::string& ref_label) {
        double got = 0.0;
        if (energy_eh(got_json, got) && std::fabs(got - ref) > kEnergyToleranceEh)
            errs.emplace_back(where, label + " is " + fixed10(got) + " hartree but " +
                                         ref_label + " is " + fixed10(ref) +
                                         " hartree; they state one energy and must "
                                         "agree within 1e-8 hartree (spec § 3.6)");
    };
    const Json* prov = get(manifest, "provenance");
    const Json* method = prov ? get(*prov, "method") : nullptr;
    double scf = 0.0;
    if (method && method->is_string() && energy_eh(get(*prov, "scf_energy"), scf)) {
        const std::string key = ascii_lower(method->as_string());
        const std::string ref_label =
            "provenance.scf_energy (method " + quote_str(method->as_string()) + ")";
        const Json* em = get(*energies, "method");
        if (em && em->is_string() && ascii_lower(em->as_string()) == key)
            compare("/energies/total", "energies.total", get(*energies, "total"), scf,
                    ref_label);
        const Json* by = get(*energies, "by_method");
        if (by && by->is_array()) {
            size_t k = 0;
            for (const Json& entry : by->elements()) {
                const Json* m = get(entry, "method");
                if (m && m->is_string() && ascii_lower(m->as_string()) == key)
                    compare("/energies/by_method/" + std::to_string(k) + "/total",
                            "energies.by_method[" + std::to_string(k) + "].total",
                            get(entry, "total"), scf, ref_label);
                ++k;
            }
        }
    }
    const Json* thermo = get(manifest, "thermochemistry");
    const Json* zpve = thermo ? get(*thermo, "zpve_eh") : nullptr;
    const Json* components = get(*energies, "components");
    if (zpve && zpve->is_number() && components && components->is_object())
        compare("/energies/components/zpve", "energies.components.zpve",
                get(*components, "zpve"), zpve->as_double(), "thermochemistry.zpve_eh");
    return errs;
}

// --- forces / hessian (spec § 5.10, § 5.11) ------------------------------------

struct RoleFinding {
    const char* code;
    std::string role, message;
};

// Shape and reference rules of a `forces` or `hessian` section: gradient
// [n_atoms, 3], stress [3, 3], matrix square with an extent 3N divisible by 3,
// meta.atom_masses_amu with N entries, and n_atoms / N equal to the atom count
// of the section's structure when that resolves (the section-level
// structure_ref, else meta.structure_ref, else the only structure section). A
// meta.structure_ref must name a structure section. Mirrors
// _derivative_errors() in python/qvf_reader.py.
inline std::vector<RoleFinding> derivative_errors(
    const std::string& sid, const std::string& kind, const Json& section,
    const Json* members, const Json* meta,
    const std::map<std::string, int64_t>& n_atoms_of,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::string* only_structure) {
    std::vector<RoleFinding> errs;
    if (!members || !members->is_object()) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    const Json no_meta = Json::object();
    const Json& m = meta && meta->is_object() ? *meta : no_meta;
    const Json* meta_ref = get(m, "structure_ref");
    if (meta_ref) {
        auto it = meta_ref->is_string() ? kinds_by_id.find(meta_ref->as_string())
                                        : kinds_by_id.end();
        if (it == kinds_by_id.end() || it->second != "structure")
            errs.push_back({"E_REF_UNRESOLVED", "meta",
                            sec + "meta.structure_ref " + show(*meta_ref) +
                                " does not name a structure section (spec § 5.10)"});
    }
    bool have_structure = false;
    std::string structure;
    if (const Json* ref = get(section, "structure_ref")) {
        if (ref->is_string()) {
            have_structure = true;
            structure = ref->as_string();
        }
    } else if (meta_ref && meta_ref->is_string()) {
        have_structure = true;
        structure = meta_ref->as_string();
    } else if (only_structure) {
        have_structure = true;
        structure = *only_structure;
    }
    bool have_n = false;
    int64_t n_struct = 0;
    if (have_structure) {
        auto it = n_atoms_of.find(structure);
        if (it != n_atoms_of.end()) {
            have_n = true;
            n_struct = it->second;
        }
    }

    std::vector<int64_t> shape;
    if (kind == "forces") {
        if (shape_of(get(*members, "gradient"), shape)) {
            if (shape.size() != 2 || shape[1] != 3)
                errs.push_back({"E_MEMBER_MALFORMED", "gradient",
                                sec + "gradient shape " + dims_str(shape) +
                                    " is not [n_atoms, 3] (spec § 5.10)"});
            else if (have_n && shape[0] != n_struct)
                errs.push_back({"E_MEMBER_MALFORMED", "gradient",
                                sec + "gradient has " + std::to_string(shape[0]) +
                                    " rows but structure " + quote_str(structure) +
                                    " has " + std::to_string(n_struct) +
                                    " atoms (spec § 5.10)"});
        }
        if (shape_of(get(*members, "stress"), shape) &&
            (shape.size() != 2 || shape[0] != 3 || shape[1] != 3))
            errs.push_back({"E_MEMBER_MALFORMED", "stress",
                            sec + "stress shape " + dims_str(shape) +
                                " is not [3, 3] (spec § 5.10)"});
        return errs;
    }

    if (!shape_of(get(*members, "matrix"), shape)) return errs;
    if (shape.size() != 2 || shape[0] != shape[1]) {
        errs.push_back({"E_MEMBER_MALFORMED", "matrix",
                        sec + "matrix shape " + dims_str(shape) +
                            " is not square (spec § 5.11)"});
        return errs;
    }
    if (shape[0] % 3 != 0) {
        errs.push_back({"E_MEMBER_MALFORMED", "matrix",
                        sec + "matrix extent " + std::to_string(shape[0]) +
                            " is not divisible by 3, so the matrix is not [3N, 3N] "
                            "(spec § 5.11)"});
        return errs;
    }
    const int64_t n = shape[0] / 3;
    if (have_n && n != n_struct)
        errs.push_back({"E_MEMBER_MALFORMED", "matrix",
                        sec + "matrix " + dims_str(shape) + " is for " +
                            std::to_string(n) + " atoms but structure " +
                            quote_str(structure) + " has " + std::to_string(n_struct) +
                            " atoms (spec § 5.11)"});
    const Json* masses = get(m, "atom_masses_amu");
    if (masses && masses->is_array() && static_cast<int64_t>(masses->size()) != n)
        errs.push_back({"E_MEMBER_MALFORMED", "meta",
                        sec + "meta.atom_masses_amu has " +
                            std::to_string(masses->size()) +
                            " entries but the matrix is for " + std::to_string(n) +
                            " atoms (spec § 5.11)"});
    return errs;
}

// --- wavefunction.bloch (spec § 5.3.1) ----------------------------------------

// Largest allowed difference between 1 and the sum of the k-point weights.
constexpr double kBlochWeightTolerance = 1e-6;

// A JSON number with an integral value (2 or 2.0).
inline bool is_integral(const Json* j) {
    if (!j || !j->is_number()) return false;
    const double v = j->as_double();
    return std::isfinite(v) && std::floor(v) == v;
}

// Shape and structure rules of a `wavefunction.bloch` section that JSON Schema
// cannot express: mo_coefficients [n_k, n_spin, n_mo, n_ao, 2] with n_spin 1
// or 2; energies and occupations with its [n_k, n_spin, n_mo] (each other's
// when the coefficients are malformed); n_k k-points, labels and weights in
// the kpoints payload, the weights summing to 1 within 1e-6; mo_metadata
// n_spin (1 for restricted, 2 for unrestricted), n_mo, n_ao and n_mo_per_k,
// and basis.n_ao, agreeing with those axes; a structure (the section-level
// structure_ref, else basis.structure_ref) with at least one periodic axis.
// Mirrors _bloch_errors() in python/qvf_reader.py. The `role` of each finding
// is its location relative to the section: "members/<role>" or
// "structure_ref".
inline std::vector<RoleFinding> bloch_errors(
    const std::string& sid, const Json& section, const Json* members,
    const std::map<std::string, Json>& docs,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::map<std::string, bool>& periodic_of) {
    std::vector<RoleFinding> errs;
    if (!members || !members->is_object()) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    auto err = [&](const char* code, const std::string& where, const std::string& msg) {
        errs.push_back({code, where, sec + msg + " (spec § 5.3.1)"});
    };
    auto doc = [&](const char* role) -> const Json* {
        auto it = docs.find(role);
        return it == docs.end() || !it->second.is_object() ? nullptr : &it->second;
    };

    std::vector<int64_t> coeff, ref, shape;
    bool have_coeff = shape_of(get(*members, "mo_coefficients"), coeff);
    if (have_coeff && (coeff.size() != 5 || coeff[4] != 2 || coeff[1] < 1 || coeff[1] > 2)) {
        err("E_MEMBER_MALFORMED", "members/mo_coefficients",
            "mo_coefficients shape " + dims_str(coeff) +
                " is not [n_k, n_spin, n_mo, n_ao, 2] with n_spin 1 or 2");
        have_coeff = false;
    }
    bool have_ref = have_coeff;
    std::string ref_role = "mo_coefficients";
    if (have_coeff) ref.assign(coeff.begin(), coeff.begin() + 3);
    for (const char* role : {"energies", "occupations"}) {
        if (!shape_of(get(*members, role), shape)) continue;
        const std::string name = role;
        if (shape.size() != 3) {
            err("E_MEMBER_MALFORMED", "members/" + name,
                name + " shape " + dims_str(shape) + " is not [n_k, n_spin, n_mo]");
        } else if (!have_ref) {
            ref = shape;
            ref_role = name;
            have_ref = true;
        } else if (shape != ref) {
            err("E_MEMBER_MALFORMED", "members/" + name,
                name + " shape " + dims_str(shape) +
                    " differs from the [n_k, n_spin, n_mo] = " + dims_str(ref) + " of " +
                    ref_role);
        }
    }

    if (const Json* kp = doc("kpoints")) {
        const Json* points = get(*kp, "kpoints");
        const bool have_nk = points && points->is_array();
        const int64_t n_k = have_nk ? static_cast<int64_t>(points->size()) : 0;
        if (have_nk && have_ref && n_k != ref[0])
            err("E_MEMBER_MALFORMED", "members/kpoints",
                "kpoints lists " + std::to_string(n_k) + " k-points but axis 0 of " +
                    ref_role + " is " + std::to_string(ref[0]));
        for (const char* key : {"weights", "labels"}) {
            const Json* arr = get(*kp, key);
            if (arr && arr->is_array() && have_nk && static_cast<int64_t>(arr->size()) != n_k)
                err("E_MEMBER_MALFORMED", "members/kpoints",
                    std::string("kpoints.") + key + " has " + std::to_string(arr->size()) +
                        " entries for " + std::to_string(n_k) + " k-points");
        }
        const Json* weights = get(*kp, "weights");
        if (weights && weights->is_array()) {
            bool numbers = true;
            double sum = 0.0;
            for (const Json& x : weights->elements()) {
                if (!x.is_number()) {
                    numbers = false;
                    break;
                }
                sum += x.as_double();
            }
            if (numbers && std::fabs(sum - 1.0) > kBlochWeightTolerance)
                err("E_MEMBER_MALFORMED", "members/kpoints",
                    "kpoints.weights sum to " + format_number(sum) +
                        ", not to 1 within 1e-6");
        }
    }

    if (const Json* meta = doc("mo_metadata")) {
        const Json* spin = get(*meta, "spin");
        const Json* n_spin = get(*meta, "n_spin");
        int64_t want = 0;
        if (spin && spin->is_string()) {
            if (spin->as_string() == "restricted") want = 1;
            else if (spin->as_string() == "unrestricted") want = 2;
        }
        if (want && is_integral(n_spin) && n_spin->as_double() != static_cast<double>(want))
            err("E_MEMBER_MALFORMED", "members/mo_metadata",
                "mo_metadata.spin " + show(*spin) + " has n_spin " + std::to_string(want) +
                    ", not " + show(*n_spin));
        static const struct {
            const char* key;
            size_t axis;
            bool of_coeff;
        } kAxes[] = {{"n_spin", 1, false}, {"n_mo", 2, false}, {"n_ao", 3, true}};
        for (const auto& a : kAxes) {
            const Json* v = get(*meta, a.key);
            const bool have = a.of_coeff ? have_coeff : have_ref;
            const std::vector<int64_t>& dims = a.of_coeff ? coeff : ref;
            if (is_integral(v) && have &&
                v->as_double() != static_cast<double>(dims[a.axis]))
                err("E_MEMBER_MALFORMED", "members/mo_metadata",
                    std::string("mo_metadata.") + a.key + " is " + show(*v) + " but axis " +
                        std::to_string(a.axis) + " of " +
                        (a.of_coeff ? std::string("mo_coefficients") : ref_role) + " is " +
                        std::to_string(dims[a.axis]));
        }
        const Json* per_k = get(*meta, "n_mo_per_k");
        if (per_k && per_k->is_array()) {
            if (have_ref && static_cast<int64_t>(per_k->size()) != ref[0])
                err("E_MEMBER_MALFORMED", "members/mo_metadata",
                    "mo_metadata.n_mo_per_k has " + std::to_string(per_k->size()) +
                        " entries but axis 0 of " + ref_role + " is " +
                        std::to_string(ref[0]));
            bool have_n_mo = have_ref;
            double n_mo = have_ref ? static_cast<double>(ref[2]) : 0.0;
            if (!have_ref && is_integral(get(*meta, "n_mo"))) {
                have_n_mo = true;
                n_mo = get(*meta, "n_mo")->as_double();
            }
            if (have_n_mo)
                for (const Json& x : per_k->elements())
                    if (x.is_number() && x.as_double() > n_mo) {
                        err("E_MEMBER_MALFORMED", "members/mo_metadata",
                            "mo_metadata.n_mo_per_k has an entry above n_mo = " +
                                format_number(n_mo));
                        break;
                    }
        }
    }

    const Json* basis = doc("basis");
    if (basis && have_coeff) {
        const Json* n_ao = get(*basis, "n_ao");
        if (is_integral(n_ao) && n_ao->as_double() != static_cast<double>(coeff[3]))
            err("E_MEMBER_MALFORMED", "members/basis",
                "basis.n_ao is " + show(*n_ao) + " but axis 3 of mo_coefficients is " +
                    std::to_string(coeff[3]));
    }

    std::string structure, where;
    auto names_structure = [&](const Json* r) {
        if (!r || !r->is_string()) return false;
        auto it = kinds_by_id.find(r->as_string());
        return it != kinds_by_id.end() && it->second == "structure";
    };
    if (const Json* r = get(section, "structure_ref")) {
        if (!names_structure(r)) return errs;   // reported by the generic rule (§ 4.1)
        structure = r->as_string();
        where = "structure_ref";
    } else if (basis) {
        const Json* r = get(*basis, "structure_ref");
        where = "members/basis";
        if (!names_structure(r)) {
            err("E_REF_UNRESOLVED", where,
                "basis.structure_ref " + (r ? show(*r) : std::string("nothing")) +
                    " does not name a structure section");
            return errs;
        }
        structure = r->as_string();
    } else {
        return errs;   // no basis payload: its own finding says why
    }
    auto it = periodic_of.find(structure);
    if (it != periodic_of.end() && !it->second)
        err("E_MEMBER_MALFORMED", where,
            "structure " + quote_str(structure) +
                " has no periodic axis; Bloch orbitals need a periodic structure");
    return errs;
}

// A wavefunction.bloch section, checked after every structure has been read.
struct BlochSection {
    size_t index;
    std::string sid;
    const Json* section;
    const Json* members;
    std::map<std::string, Json> docs;
};

// A forces / hessian section, checked after every structure has been read.
struct DerivativeSection {
    size_t index;
    std::string sid, kind;
    const Json* section;
    const Json* members;
    bool has_meta;
    Json meta;
};
// --- excited states, the spectra that name them, transition densities ---------
// (spec § 5.2, § 5.4, § 5.12); mirrors _excited_states_errors(),
// _spectrum_states_errors() and _transition_density_errors() in
// python/qvf_reader.py.

inline bool is_int(const Json* j) { return j && j->type() == Json::Type::Int; }

// Orbital rows of a wavefunction.gto section by channel, from the manifest
// shapes: "" for mo_coefficients, "alpha" / "beta" for the unrestricted pair.
using MoRows = std::map<std::string, int64_t>;

inline MoRows mo_rows(const Json* members) {
    MoRows rows;
    if (!members || !members->is_object()) return rows;
    static const struct {
        const char* role;
        const char* key;
    } kRoles[] = {{"mo_coefficients", ""},
                  {"mo_coefficients_alpha", "alpha"},
                  {"mo_coefficients_beta", "beta"}};
    std::vector<int64_t> shape;
    for (const auto& r : kRoles)
        if (shape_of(get(*members, r.role), shape) && !shape.empty()) rows[r.key] = shape[0];
    return rows;
}

// The rows a contribution with `spin` may index: those of mo_coefficients,
// else of the channel spin selects, else (no spin) of both channels.
inline bool contribution_rows(const MoRows& rows, const Json* spin, int64_t& n) {
    auto it = rows.find("");
    if (it != rows.end()) {
        n = it->second;
        return true;
    }
    if (spin && spin->is_string() &&
        (spin->as_string() == "alpha" || spin->as_string() == "beta")) {
        auto c = rows.find(spin->as_string());
        if (c == rows.end()) return false;
        n = c->second;
        return true;
    }
    bool any = false;
    for (const char* ch : {"alpha", "beta"}) {
        auto c = rows.find(ch);
        if (c == rows.end()) continue;
        n = any ? std::min(n, c->second) : c->second;
        any = true;
    }
    return any;
}

// The `index` values of an excited_states payload's states; `known` is false
// when the payload has no `states` array.
struct StateIndices {
    bool known = false;
    std::set<int64_t> values;
};

inline StateIndices state_indices(const Json* payload) {
    StateIndices out;
    const Json* states = payload ? get(*payload, "states") : nullptr;
    if (!states || !states->is_array()) return out;
    out.known = true;
    for (const Json& s : states->elements()) {
        const Json* idx = get(s, "index");
        if (is_int(idx)) out.values.insert(idx->as_int());
    }
    return out;
}

// A finding about a section field (`field` non-empty) or, with an empty
// `field`, about the section's payload member.
struct FieldFinding {
    const char* code;
    std::string field, message;
};

inline std::vector<FieldFinding> excited_states_errors(
    const std::string& sid, const Json& section, const Json* payload,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::map<std::string, MoRows>& mo_rows_of) {
    std::vector<FieldFinding> errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    auto kind_of = [&](const Json* ref) -> std::string {
        if (!ref || !ref->is_string()) return "";
        auto it = kinds_by_id.find(ref->as_string());
        return it == kinds_by_id.end() ? "" : it->second;
    };
    bool have_wf = false;
    std::string wf_ref;
    if (const Json* ref = get(section, "wavefunction_ref")) {
        if (kind_of(ref) != "wavefunction.gto")
            errs.push_back({"E_REF_UNRESOLVED", "wavefunction_ref",
                            sec + "wavefunction_ref " + show(*ref) +
                                " does not name a wavefunction.gto section (spec § 5.12)"});
        else {
            have_wf = true;
            wf_ref = ref->as_string();
        }
    }
    const Json* states = payload ? get(*payload, "states") : nullptr;
    if (!states || !states->is_array()) return errs;
    MoRows rows;
    if (have_wf) {
        auto it = mo_rows_of.find(wf_ref);
        if (it != mo_rows_of.end()) rows = it->second;
    }
    std::set<int64_t> seen;
    bool have_prev = false;
    int64_t prev = 0;
    const std::vector<Json>& list = states->elements();
    for (size_t k = 0; k < list.size(); ++k) {
        const Json& state = list[k];
        if (!state.is_object()) continue;
        const std::string at = "states[" + std::to_string(k) + "]";
        const Json* idx = get(state, "index");
        if (is_int(idx)) {
            const int64_t v = idx->as_int();
            if (seen.count(v))
                errs.push_back({"E_MEMBER_MALFORMED", "",
                                sec + at + ".index " + std::to_string(v) +
                                    " repeats an earlier index; state indices are unique "
                                    "(spec § 5.12)"});
            else if (have_prev && v < prev)
                errs.push_back({"E_MEMBER_MALFORMED", "",
                                sec + at + ".index " + std::to_string(v) + " follows index " +
                                    std::to_string(prev) +
                                    "; states are listed in ascending order of index "
                                    "(spec § 5.12)"});
            seen.insert(v);
            prev = v;
            have_prev = true;
        }
        if (const Json* ref = get(state, "nto_ref"))
            if (kind_of(ref) != "wavefunction.gto")
                errs.push_back({"E_REF_UNRESOLVED", "",
                                sec + at + ".nto_ref " + show(*ref) +
                                    " does not name a wavefunction.gto section (spec § 5.12)"});
        if (const Json* ref = get(state, "density_ref")) {
            const std::string k2 = kind_of(ref);
            if (k2 != "volume.transition_density" && k2 != "volume.difference")
                errs.push_back({"E_REF_UNRESOLVED", "",
                                sec + at + ".density_ref " + show(*ref) +
                                    " does not name a volume.transition_density or "
                                    "volume.difference section (spec § 5.12)"});
        }
        const Json* contributions = get(state, "contributions");
        if (!have_wf || !contributions || !contributions->is_array()) continue;
        const std::vector<Json>& cs = contributions->elements();
        for (size_t c = 0; c < cs.size(); ++c) {
            if (!cs[c].is_object()) continue;
            int64_t n = 0;
            if (!contribution_rows(rows, get(cs[c], "spin"), n)) continue;
            for (const char* key : {"from_mo", "to_mo"}) {
                const Json* mo = get(cs[c], key);
                if (is_int(mo) && (mo->as_int() < 0 || mo->as_int() >= n))
                    errs.push_back({"E_MEMBER_MALFORMED", "",
                                    sec + at + ".contributions[" + std::to_string(c) +
                                        "]." + key + " " + std::to_string(mo->as_int()) +
                                        " is not one of the " + std::to_string(n) +
                                        " orbital rows of wavefunction " +
                                        quote_str(wf_ref) + " (spec § 5.12)"});
            }
        }
    }
    const Json* meta = get(*payload, "meta");
    const Json* n_states = meta ? get(*meta, "n_states") : nullptr;
    if (is_int(n_states) && n_states->as_int() != static_cast<int64_t>(list.size()))
        errs.push_back({"E_MEMBER_MALFORMED", "",
                        sec + "meta.n_states is " + std::to_string(n_states->as_int()) +
                            " but states has " + std::to_string(list.size()) +
                            " entries (spec § 5.12)"});
    return errs;
}

// A one-dimensional spectrum's states_ref names an excited_states section, and
// every peaks[].assignment.state_index is the index of one of its states.
inline std::vector<std::string> spectrum_states_errors(
    const std::string& sid, const Json& payload,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::map<std::string, StateIndices>& indices_of) {
    std::vector<std::string> errs;
    const Json* ref = get(payload, "states_ref");
    if (!ref) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    auto it = ref->is_string() ? kinds_by_id.find(ref->as_string()) : kinds_by_id.end();
    if (it == kinds_by_id.end() || it->second != "excited_states") {
        errs.push_back(sec + "spectrum states_ref " + show(*ref) +
                       " does not name an excited_states section (spec § 5.4)");
        return errs;
    }
    auto ix = indices_of.find(ref->as_string());
    const Json* peaks = get(payload, "peaks");
    if (ix == indices_of.end() || !ix->second.known || !peaks || !peaks->is_array())
        return errs;
    const std::vector<Json>& list = peaks->elements();
    for (size_t k = 0; k < list.size(); ++k) {
        const Json* assignment = get(list[k], "assignment");
        const Json* state = assignment ? get(*assignment, "state_index") : nullptr;
        if (is_int(state) && !ix->second.values.count(state->as_int()))
            errs.push_back(sec + "peaks[" + std::to_string(k) + "].assignment.state_index " +
                           std::to_string(state->as_int()) +
                           " names no state of excited_states section " +
                           quote_str(ref->as_string()) + " (spec § 5.4)");
    }
    return errs;
}


// --- property kinds (spec § 5.13) ----------------------------------------------

inline bool is_property_kind(const std::string& kind) {
    return kind == "properties.electric" || kind == "properties.elastic" ||
           kind == "properties.dielectric";
}

// Largest allowed disagreement between two statements of one property: matrix
// symmetry (relative to the largest entry), static against electronic + ionic,
// and the electric dipole against the root dipole (debye).
constexpr double kPropertyTolerance = 1e-6;

// `v` as a rows x cols matrix of finite numbers, row-major in `out`; false when
// it is anything else (the payload schema reports that). Mirrors _matrix_of()
// in python/qvf_reader.py.
inline bool matrix_of(const Json* v, size_t rows, size_t cols, std::vector<double>& out) {
    out.clear();
    if (!v || !v->is_array() || v->size() != rows) return false;
    for (const Json& row : v->elements()) {
        if (!row.is_array() || row.size() != cols) return false;
        for (const Json& x : row.elements()) {
            if (!x.is_number() || !std::isfinite(x.as_double())) return false;
            out.push_back(x.as_double());
        }
    }
    return true;
}

// `v` as three finite numbers in `out`; false when it is anything else.
inline bool vector3_of(const Json* v, std::vector<double>& out) {
    out.clear();
    if (!v || !v->is_array() || v->size() != 3) return false;
    for (const Json& x : v->elements()) {
        if (!x.is_number() || !std::isfinite(x.as_double())) return false;
        out.push_back(x.as_double());
    }
    return true;
}

// The first (i, j), i < j, whose two entries of the n x n matrix differ by more
// than 1e-6 of the largest absolute entry; false when it is symmetric.
inline bool first_asymmetry(const std::vector<double>& m, size_t n, size_t& bi, size_t& bj) {
    double scale = 0.0;
    for (double x : m) scale = std::max(scale, std::fabs(x));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            if (std::fabs(m[i * n + j] - m[j * n + i]) > kPropertyTolerance * scale) {
                bi = i;
                bj = j;
                return true;
            }
    return false;
}

struct CodeFinding {
    const char* code;
    std::string message;
};

// Rules of the property kinds that JSON Schema cannot express (spec § 5.13):
// symmetric elastic matrices, static = electronic + ionic, one Born tensor per
// atom of the section's structure (the section-level structure_ref, else the
// only structure section) when that resolves, and a dipole that agrees with the
// root dipole_moment within 1e-6 debye. Mirrors _property_errors() in
// python/qvf_reader.py.
inline std::vector<CodeFinding> property_errors(
    const std::string& sid, const std::string& kind, const Json& section,
    const Json& payload, const Json& manifest,
    const std::map<std::string, int64_t>& n_atoms_of, const std::string* only_structure) {
    std::vector<CodeFinding> errs;
    if (!payload.is_object()) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    std::vector<double> m;
    if (kind == "properties.elastic") {
        for (const char* name : {"stiffness", "compliance"}) {
            size_t i = 0, j = 0;
            if (matrix_of(get(payload, name), 6, 6, m) && first_asymmetry(m, 6, i, j))
                errs.push_back({"E_MEMBER_MALFORMED",
                                sec + name + " is not symmetric: entry (" + std::to_string(i) +
                                    ", " + std::to_string(j) + ") is " +
                                    format_number(m[i * 6 + j]) + " but entry (" +
                                    std::to_string(j) + ", " + std::to_string(i) + ") is " +
                                    format_number(m[j * 6 + i]) +
                                    "; they must agree within 1e-6 of the largest entry "
                                    "(spec § 5.13)"});
        }
    } else if (kind == "properties.dielectric") {
        std::vector<double> e, ion, s;
        if (matrix_of(get(payload, "electronic"), 3, 3, e) &&
            matrix_of(get(payload, "ionic"), 3, 3, ion) &&
            matrix_of(get(payload, "static"), 3, 3, s)) {
            for (size_t k = 0; k < 9; ++k) {
                if (std::fabs(s[k] - (e[k] + ion[k])) > kPropertyTolerance) {
                    errs.push_back({"E_MEMBER_MALFORMED",
                                    sec + "static[" + std::to_string(k / 3) + "][" +
                                        std::to_string(k % 3) + "] is " + format_number(s[k]) +
                                        " but electronic + ionic is " +
                                        format_number(e[k] + ion[k]) +
                                        "; static must equal their sum within 1e-6 "
                                        "(spec § 5.13)"});
                    break;
                }
            }
        }
        const Json* born = get(payload, "born_charges");
        if (born && born->is_array()) {
            bool have = false;
            std::string structure;
            if (const Json* ref = get(section, "structure_ref")) {
                if (ref->is_string()) {
                    have = true;
                    structure = ref->as_string();
                }
            } else if (only_structure) {
                have = true;
                structure = *only_structure;
            }
            auto it = have ? n_atoms_of.find(structure) : n_atoms_of.end();
            if (it != n_atoms_of.end() && static_cast<int64_t>(born->size()) != it->second)
                errs.push_back({"E_MEMBER_MALFORMED",
                                sec + "born_charges has " + std::to_string(born->size()) +
                                    " tensors but structure " + quote_str(structure) +
                                    " has " + std::to_string(it->second) +
                                    " atoms; there must be one per atom (spec § 5.13)"});
        }
    } else {  // properties.electric
        const Json* dipole = get(payload, "dipole");
        const Json* root = get(manifest, "dipole_moment");
        if (dipole && dipole->is_object() && root && root->is_object()) {
            std::vector<double> vec, root_vec;
            const Json* v = get(*dipole, "vector");
            const Json* rv = get(*root, "vector_debye");
            if (vector3_of(v, vec) && vector3_of(rv, root_vec)) {
                bool differs = false;
                for (size_t k = 0; k < 3; ++k)
                    if (std::fabs(vec[k] - root_vec[k]) > kPropertyTolerance) differs = true;
                if (differs)
                    errs.push_back({"E_ROOT_INCONSISTENT",
                                    sec + "dipole.vector " + show(*v) +
                                        " differs from the root dipole_moment.vector_debye " +
                                        show(*rv) +
                                        "; the two state one dipole and must agree within "
                                        "1e-6 debye (spec § 5.13)"});
            }
            const Json* total = get(*dipole, "total");
            const Json* root_total = get(*root, "total_debye");
            if (total && root_total && total->is_number() && root_total->is_number() &&
                std::isfinite(total->as_double()) && std::isfinite(root_total->as_double()) &&
                std::fabs(total->as_double() - root_total->as_double()) > kPropertyTolerance)
                errs.push_back({"E_ROOT_INCONSISTENT",
                                sec + "dipole.total " + show(*total) +
                                    " differs from the root dipole_moment.total_debye " +
                                    show(*root_total) +
                                    "; the two state one dipole and must agree within 1e-6 "
                                    "debye (spec § 5.13)"});
        }
    }
    return errs;
}

// A volume.transition_density section's excited_states_ref names an
// excited_states section, and its state_index the index of one of its states.
// Returns (field, message) pairs.
inline std::vector<std::pair<std::string, std::string>> transition_density_errors(
    const std::string& sid, const Json& section,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::map<std::string, StateIndices>& indices_of) {
    std::vector<std::pair<std::string, std::string>> errs;
    const Json* ref = get(section, "excited_states_ref");
    if (!ref) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    auto it = ref->is_string() ? kinds_by_id.find(ref->as_string()) : kinds_by_id.end();
    if (it == kinds_by_id.end() || it->second != "excited_states") {
        errs.emplace_back("excited_states_ref",
                          sec + "excited_states_ref " + show(*ref) +
                              " does not name an excited_states section (spec § 5.2)");
        return errs;
    }
    const Json* state = get(section, "state_index");
    auto ix = indices_of.find(ref->as_string());
    if (is_int(state) && ix != indices_of.end() && ix->second.known &&
        !ix->second.values.count(state->as_int()))
        errs.emplace_back("state_index",
                          sec + "state_index " + std::to_string(state->as_int()) +
                              " names no state of excited_states section " +
                              quote_str(ref->as_string()) + " (spec § 5.2)");
    return errs;
}

// Volume kinds that may name the wavefunction they were computed from
// (spec § 5.2, Volumes of a wavefunction).
inline bool is_wavefunction_volume_kind(const std::string& kind) {
    return kind == "volume.density" || kind == "volume.spin" || kind == "volume.orbital";
}

// A volume.density, volume.spin or volume.orbital section's wavefunction_ref
// names a wavefunction.gto section. For a volume.orbital, spin is absent when
// that section carries mo_coefficients and present, whenever mo_index is,
// when it carries the alpha and beta coefficients, and mo_index is a row of
// the member they select, from the manifest shapes. Mirrors
// _wavefunction_volume_errors() in python/qvf/validate.py. Returns (field,
// message) pairs, all E_REF_UNRESOLVED.
inline std::vector<std::pair<std::string, std::string>> wavefunction_volume_errors(
    const std::string& sid, const std::string& kind, const Json& section,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::map<std::string, MoRows>& mo_rows_of) {
    std::vector<std::pair<std::string, std::string>> errs;
    const Json* ref = get(section, "wavefunction_ref");
    if (!ref) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    auto it = ref->is_string() ? kinds_by_id.find(ref->as_string()) : kinds_by_id.end();
    if (it == kinds_by_id.end() || it->second != "wavefunction.gto") {
        errs.emplace_back("wavefunction_ref",
                          sec + "wavefunction_ref " + show(*ref) +
                              " does not name a wavefunction.gto section (spec § 5.2)");
        return errs;
    }
    if (kind != "volume.orbital") return errs;
    const std::string& wf = ref->as_string();
    MoRows rows;
    auto r = mo_rows_of.find(wf);
    if (r != mo_rows_of.end()) rows = r->second;
    const Json* mo = get(section, "mo_index");
    const Json* spin = get(section, "spin");
    bool known = false;
    int64_t n = 0;
    std::string member = "mo_coefficients";
    auto single = rows.find("");
    if (single != rows.end()) {
        if (spin)
            errs.emplace_back("spin", sec + "spin " + show(*spin) +
                                          " selects a spin channel, but wavefunction " +
                                          quote_str(wf) +
                                          " carries one set of orbitals in "
                                          "mo_coefficients (spec § 5.2)");
        known = true;
        n = single->second;
    } else if (rows.count("alpha") || rows.count("beta")) {
        if (!spin) {
            if (mo)
                errs.emplace_back("mo_index", sec + "mo_index " + show(*mo) +
                                                  " needs spin: wavefunction " +
                                                  quote_str(wf) +
                                                  " carries alpha and beta coefficients "
                                                  "(spec § 5.2)");
            return errs;
        }
        if (spin->is_string() &&
            (spin->as_string() == "alpha" || spin->as_string() == "beta")) {
            auto c = rows.find(spin->as_string());
            if (c != rows.end()) {
                known = true;
                n = c->second;
                member = "mo_coefficients_" + spin->as_string();
            }
        }
    }
    if (is_int(mo) && known && (mo->as_int() < 0 || mo->as_int() >= n))
        errs.emplace_back("mo_index", sec + "mo_index " + std::to_string(mo->as_int()) +
                                          " is not one of the " + std::to_string(n) +
                                          " rows of " + member + " of wavefunction " +
                                          quote_str(wf) + " (spec § 5.2)");
    return errs;
}

// An excited_states section, a spectrum with states_ref, a transition
// density, or a standalone bonds section, checked after every section has
// been read.
struct DeferredSection {
    size_t index;
    std::string sid;
    const Json* section;
    bool has_doc;
    Json doc;
};

// Atom indices of a standalone `bonds` section against the structure its
// structure_ref names (spec § 5.1). Only a section whose structure_ref names
// a structure section with a known atom count is checked; a reference that
// does not resolve is the generic reference check's finding, and a section
// without one is not checked. Every integer i and j of the payload's `pairs`
// must be an atom of that structure. Returns messages, located at the `bonds`
// member. Mirrors _bonds_errors() in python/qvf_reader.py.
inline std::vector<std::string> bonds_errors(
    const std::string& sid, const Json& section, const Json* payload,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::map<std::string, int64_t>& n_atoms_of) {
    std::vector<std::string> errs;
    const Json* ref = get(section, "structure_ref");
    if (!ref || !ref->is_string()) return errs;
    const std::string structure = ref->as_string();
    auto kind = kinds_by_id.find(structure);
    if (kind == kinds_by_id.end() || kind->second != "structure") return errs;
    auto count = n_atoms_of.find(structure);
    const Json* pairs = payload && payload->is_object() ? get(*payload, "pairs") : nullptr;
    if (count == n_atoms_of.end() || !pairs || !pairs->is_array()) return errs;
    const int64_t n = count->second;
    for (size_t k = 0; k < pairs->size(); ++k) {
        const Json& pair = pairs->at(k);
        if (!pair.is_object()) continue;
        for (const char* key : {"i", "j"}) {
            const Json* atom = get(pair, key);
            if (is_int(atom) && (atom->as_int() < 0 || atom->as_int() >= n))
                errs.push_back("section " + quote_str(sid) + ": pairs[" +
                               std::to_string(k) + "]." + key + " " +
                               std::to_string(atom->as_int()) + " is not one of the " +
                               std::to_string(n) + " atoms of structure " +
                               quote_str(structure) +
                               ", which structure_ref names (spec § 5.1)");
        }
    }
    return errs;
}

// A properties.* section, checked after every structure has been read.
struct PropertySection {
    size_t index;
    std::string sid, kind;
    const Json* section;
    Json payload;
};

// Custom roles of an `atom_properties` section need a `meta` entry (spec § 5.7);
// mirrors _atom_properties_meta_errors() in python/qvf_reader.py. Every
// custom.<name> member needs an object entry of that name in the section's
// `meta` member; a section without one has none. A `meta` member that is
// declared but could not be read (its own finding says why) is not judged.
// Returns (role, message) pairs.
inline std::vector<std::pair<std::string, std::string>> atom_properties_meta_errors(
    const std::string& sid, const Json& members, const std::map<std::string, Json>& docs) {
    std::vector<std::pair<std::string, std::string>> errs;
    std::vector<std::string> custom;
    for (const auto& kv : members.items())
        if (starts_with(kv.first, "custom.")) custom.push_back(kv.first);
    if (custom.empty()) return errs;
    const std::string sec = "section " + quote_str(sid) + ": custom role ";
    if (!members.contains("meta")) {
        for (const auto& r : custom)
            errs.emplace_back(r, sec + quote_str(r) +
                                     " has no meta member to describe it (spec § 5.7)");
        return errs;
    }
    const auto it = docs.find("meta");
    if (it == docs.end()) return errs;
    if (!it->second.is_object()) {
        for (const auto& r : custom)
            errs.emplace_back(r, sec + quote_str(r) +
                                     ": the meta member is not a JSON object (spec § 5.7)");
        return errs;
    }
    for (const auto& r : custom) {
        const Json* entry = get(it->second, r);
        if (!entry || !entry->is_object())
            errs.emplace_back(r, sec + quote_str(r) +
                                     " has no entry in the meta member (spec § 5.7)");
    }
    return errs;
}

// A member spec's `size` against the uncompressed size its ZIP entry records
// (spec § 4.2): nullptr when absent or equal; otherwise the finding code,
// E_MEMBER_MALFORMED (not a non-negative integer) or E_BINARY_SIZE, with
// `message` set. Reads the central directory only.
inline const char* size_problem(const Json& m, const ZipReader::Entry& e,
                                std::string& message) {
    const Json* sj = get(m, "size");
    if (!sj) return nullptr;
    if (sj->type() != Json::Type::Int || sj->as_int() < 0) {
        message = "member " + quote_str(e.name) +
                  ": size must be a non-negative integer, got " + show(*sj);
        return "E_MEMBER_MALFORMED";
    }
    if (static_cast<uint64_t>(sj->as_int()) != static_cast<uint64_t>(e.usize)) {
        message = "member " + quote_str(e.name) + ": declared size " +
                  std::to_string(sj->as_int()) + " != " + std::to_string(e.usize) +
                  " bytes recorded in the ZIP central directory";
        return "E_BINARY_SIZE";
    }
    return nullptr;
}

// --- wavefunction.gto (spec § 5.3, Appendix A.1) ------------------------------

// Tolerances of the electron count and of the Appendix A.1 normalization check.
constexpr double kNElectronsTolerance = 1e-6;
constexpr double kNormalizationRtol = 1e-6;

// One shell of a `basis` payload as the checks read it. `usable` is false when
// its primitive arrays are not numbers, an exponent is not above 0, or the
// lengths differ (the schema check reports those).
struct ShellInfo {
    int l = 0;
    bool pure = true;
    bool usable = false;
    std::vector<double> exponents, coefficients;
};

inline bool is_plain_number(const Json& j) {
    return j.type() == Json::Type::Int || j.type() == Json::Type::Double;
}

// The shells of a `basis` payload; false when the AO layout cannot be read
// (no shell list, or a shell whose `l` or `pure` is malformed). Mirrors
// _basis_shells() in python/qvf_reader.py.
inline bool basis_shells(const Json* basis, std::vector<ShellInfo>& out) {
    out.clear();
    const Json* shells = basis ? get(*basis, "shells") : nullptr;
    if (!shells || !shells->is_array()) return false;
    bool default_pure = true;
    if (const Json* p = get(*basis, "pure")) {
        if (!p->is_bool()) return false;
        default_pure = p->as_bool();
    }
    for (const Json& sh : shells->elements()) {
        if (!sh.is_object()) return false;
        const Json* l = get(sh, "l");
        if (!l || l->type() != Json::Type::Int || l->as_int() < 0 || l->as_int() > 7)
            return false;
        ShellInfo info;
        info.l = static_cast<int>(l->as_int());
        info.pure = default_pure;
        if (const Json* p = get(sh, "pure")) {
            if (!p->is_bool()) return false;
            info.pure = p->as_bool();
        }
        const Json* e = get(sh, "exponents");
        const Json* c = get(sh, "coefficients");
        info.usable = e && c && e->is_array() && c->is_array() &&
                      e->size() == c->size() && e->size() > 0;
        if (info.usable)
            for (size_t k = 0; k < e->size(); ++k) {
                const Json& a = e->at(k);
                const Json& b = c->at(k);
                if (!is_plain_number(a) || !(a.as_double() > 0) || !is_plain_number(b)) {
                    info.usable = false;
                    break;
                }
                info.exponents.push_back(a.as_double());
                info.coefficients.push_back(b.as_double());
            }
        if (!info.usable) {
            info.exponents.clear();
            info.coefficients.clear();
        }
        out.push_back(std::move(info));
    }
    return true;
}

// --- volumes (spec § 5.2) ------------------------------------------------------

// Data shape and series axis of a volume section: volume.laplacian data is
// [nx, ny, nz], volume.vector [nx, ny, nz, 3], volume.complex [nx, ny, nz, 2]
// and volume.series [n, nx, ny, nz]; a series axis has one `values` entry,
// and one `labels` entry when it has labels, per slab, and its
// trajectory_ref names a trajectory or reaction.path section. Mirrors
// _volume_errors() in python/qvf_reader.py.
inline std::vector<RoleFinding> volume_errors(
    const std::string& sid, const std::string& kind, const Json* members,
    const std::map<std::string, Json>& docs,
    const std::map<std::string, std::string>& kinds_by_id) {
    std::vector<RoleFinding> errs;
    size_t rank = 0;
    int64_t trailing = 0;
    const char* form = "";
    if (kind == "volume.laplacian") { rank = 3; form = "[nx, ny, nz]"; }
    else if (kind == "volume.vector") { rank = 4; trailing = 3; form = "[nx, ny, nz, 3]"; }
    else if (kind == "volume.complex") { rank = 4; trailing = 2; form = "[nx, ny, nz, 2]"; }
    else if (kind == "volume.series") { rank = 4; form = "[n, nx, ny, nz]"; }
    if (!rank || !members || !members->is_object()) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    std::vector<int64_t> shape;
    bool have_shape = shape_of(get(*members, "data"), shape);
    if (have_shape && (shape.size() != rank || (trailing && shape.back() != trailing))) {
        errs.push_back({"E_MEMBER_MALFORMED", "data",
                        sec + kind + " data shape " + dims_str(shape) + " is not " +
                            form + " (spec § 5.2)"});
        have_shape = false;
    }
    auto it = docs.find("axis");
    if (kind != "volume.series" || it == docs.end() || !it->second.is_object()) return errs;
    const Json& axis = it->second;
    for (const char* key : {"values", "labels"}) {
        const Json* seq = get(axis, key);
        if (have_shape && seq && seq->is_array() &&
            static_cast<int64_t>(seq->size()) != shape[0])
            errs.push_back({"E_MEMBER_MALFORMED", "axis",
                            sec + "axis." + key + " has " + std::to_string(seq->size()) +
                                " entries but data " + dims_str(shape) + " has " +
                                std::to_string(shape[0]) +
                                " slabs; it needs one per slab (spec § 5.2)"});
    }
    if (const Json* ref = get(axis, "trajectory_ref")) {
        auto k = ref->is_string() ? kinds_by_id.find(ref->as_string()) : kinds_by_id.end();
        if (k == kinds_by_id.end() ||
            (k->second != "trajectory" && k->second != "reaction.path"))
            errs.push_back({"E_REF_UNRESOLVED", "axis",
                            sec + "axis.trajectory_ref " + show(*ref) +
                                " does not name a trajectory or reaction.path section "
                                "(spec § 5.2)"});
    }
    return errs;
}

// True for a member whose payload is a volume grid ($defs/VolumeGridPayload):
// the `grid` of a canonical volume.* or basis.ao section, and the
// `volume_grid` of a reaction.path section.
inline bool is_grid_payload(const std::string& kind, const std::string& role) {
    if (role == "volume_grid") return kind == "reaction.path";
    return role == "grid" && is_canonical_kind(kind) &&
           (starts_with(kind, "volume.") || kind == "basis.ao");
}

// Three numbers, or false.
inline bool vec3_of(const Json* j, std::array<double, 3>& out) {
    if (!j || !j->is_array() || j->size() != 3) return false;
    for (size_t i = 0; i < 3; ++i) {
        if (!j->at(i).is_number()) return false;
        out[i] = j->at(i).as_double();
    }
    return true;
}

inline int64_t shell_size(int l, bool pure) {
    return pure ? 2 * l + 1 : (l + 1) * (l + 2) / 2;
}

// n_ao of a `basis` payload: counted from its shells, or its `n_ao` when the
// shells cannot be read (spec § 5.3); false when neither works.
inline bool basis_n_ao(const Json* basis, int64_t& n_ao) {
    std::vector<ShellInfo> shells;
    if (basis_shells(basis, shells)) {
        n_ao = 0;
        for (const auto& s : shells) n_ao += shell_size(s.l, s.pure);
        return true;
    }
    const Json* n = basis ? get(*basis, "n_ao") : nullptr;
    if (n && n->type() == Json::Type::Int && n->as_int() >= 0) {
        n_ao = n->as_int();
        return true;
    }
    return false;
}

// Self-overlap of the axial component (l, 0, 0) of a contracted shell on the
// N_i-normalized primitives of Appendix A.1:
// sum_ij c_i c_j (2 sqrt(a_i a_j) / (a_i + a_j))^(l + 3/2).
inline double axial_self_overlap(const ShellInfo& s) {
    double total = 0.0;
    for (size_t i = 0; i < s.exponents.size(); ++i)
        for (size_t j = 0; j < s.exponents.size(); ++j) {
            const double a = s.exponents[i], b = s.exponents[j];
            total += s.coefficients[i] * s.coefficients[j] *
                     std::pow(2.0 * std::sqrt(a * b) / (a + b), s.l + 1.5);
        }
    return total;
}

inline double odd_double_factorial(int n) {  // (2n - 1)!!, (-1)!! = 1
    double out = 1.0;
    for (int k = 1; k < 2 * n; k += 2) out *= k;
    return out;
}

// (name, factor) of each component of a shell in Appendix A.2 order: its
// diagonal overlap is factor * S_axial.
inline std::vector<std::pair<std::string, double>> shell_components(int l, bool pure) {
    std::vector<std::pair<std::string, double>> out;
    if (pure) {
        for (int m = -l; m <= l; ++m)
            out.emplace_back(m > 0 ? "m = +" + std::to_string(m)
                                   : "m = " + std::to_string(m),
                             1.0);
        return out;
    }
    for (int lx = l; lx >= 0; --lx)
        for (int ly = l - lx; ly >= 0; --ly) {
            const int lz = l - lx - ly;
            std::string name = std::string(static_cast<size_t>(lx), 'x') +
                               std::string(static_cast<size_t>(ly), 'y') +
                               std::string(static_cast<size_t>(lz), 'z');
            out.emplace_back("Cartesian " + (name.empty() ? std::string("s") : name) +
                                 " [" + std::to_string(lx) + ", " + std::to_string(ly) +
                                 ", " + std::to_string(lz) + "]",
                             odd_double_factorial(lx) * odd_double_factorial(ly) *
                                 odd_double_factorial(lz) / odd_double_factorial(l));
        }
    return out;
}

inline std::string g10(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.10g", v);
    return buf;
}

// The Appendix A.1 check of a wavefunction.gto `overlap` member: every
// diagonal entry agrees with the self-overlap the basis gives for that AO, to
// a relative tolerance of 1e-6. Runs only on a well-formed float64
// [n_ao, n_ao] overlap of a basis whose AO layout can be read. Returns at most
// one message. Mirrors _normalization_errors() in python/qvf_reader.py.
inline std::vector<std::string> normalization_errors(const std::string& sid,
                                                     const Json* basis,
                                                     const Json* overlap,
                                                     const std::vector<uint8_t>* bytes) {
    std::vector<std::string> errs;
    std::vector<ShellInfo> shells;
    std::vector<int64_t> shape;
    if (!bytes || !basis_shells(basis, shells) || !shape_of(overlap, shape)) return errs;
    int64_t n_ao = 0;
    for (const auto& s : shells) n_ao += shell_size(s.l, s.pure);
    const Json* dtype = get(*overlap, "dtype");
    if (shape.size() != 2 || shape[0] != n_ao || shape[1] != n_ao || !dtype ||
        !dtype->is_string() || dtype->as_string() != "float64" ||
        bytes->size() != static_cast<size_t>(8 * n_ao * n_ao))
        return errs;
    size_t n_bad = 0;
    std::string first;
    int64_t ao = 0;
    for (size_t k = 0; k < shells.size(); ++k) {
        const ShellInfo& s = shells[k];
        const auto components = shell_components(s.l, s.pure);
        if (s.usable) {
            const double s_axial = axial_self_overlap(s);
            for (size_t j = 0; j < components.size(); ++j) {
                const int64_t i = ao + static_cast<int64_t>(j);
                const double want = s_axial * components[j].second;
                uint64_t bits = 0;
                const size_t off = static_cast<size_t>(8 * (i * n_ao + i));
                for (int b = 0; b < 8; ++b)
                    bits |= static_cast<uint64_t>((*bytes)[off + static_cast<size_t>(b)])
                            << (8 * b);
                double got;
                std::memcpy(&got, &bits, 8);
                if (!(std::fabs(got - want) <= kNormalizationRtol * std::fabs(want))) {
                    if (n_bad++ == 0)
                        first = "AO " + std::to_string(i) + " (shell " + std::to_string(k) +
                                ", l = " + std::to_string(s.l) + ", " +
                                components[j].first + ") is " + g10(got) +
                                " but the basis gives " + g10(want);
                }
            }
        }
        ao += static_cast<int64_t>(components.size());
    }
    if (n_bad)
        errs.push_back("section " + quote_str(sid) + ": " + std::to_string(n_bad) + " of " +
                       std::to_string(n_ao) +
                       " overlap diagonal entries disagree with the self-overlaps "
                       "Appendix A.1 gives for the basis; first: " +
                       first +
                       ". The basis.json coefficients are not on N_i-normalized "
                       "primitives, or the AO order disagrees with Appendix A.2 "
                       "(spec § 5.3)");
    return errs;
}

// Whether mo_metadata declares complex coefficients (spec § 5.3): `complex`,
// else `coefficient_encoding`, else real. -1 when the declaration is
// malformed (the schema check reports it), else 0 or 1.
inline int coefficients_complex(const Json* meta) {
    if (!meta || !meta->is_object()) return -1;
    if (const Json* c = get(*meta, "complex")) return c->is_bool() ? (c->as_bool() ? 1 : 0) : -1;
    if (const Json* e = get(*meta, "coefficient_encoding")) {
        if (e->is_string() && e->as_string() == "complex_split_last_axis") return 1;
        if (e->is_string() && e->as_string() == "real") return 0;
        return -1;
    }
    return 0;
}

// Unrestricted NTOs (spec § 5.3): local channel rules, then forward pair
// references. Unit-sum normalization is conditional and is not imposed.
inline bool nto_text(const Json& meta, const char* key, const char* value) {
    const Json* p = get(meta, key);
    return p && p->is_string() && p->as_string() == value;
}
inline bool unrestricted_nto(const Json& meta) {
    return get(meta, "nto_role") && nto_text(meta, "orbital_kind", "nto") &&
           nto_text(meta, "spin", "unrestricted");
}
inline std::vector<std::string> unrestricted_nto_errors(const Json& meta, const Json* members) {
    std::vector<std::string> errors;
    if (!unrestricted_nto(meta)) return errors;
    const MoRows rows = mo_rows(members);
    for (const char* ch : {"alpha", "beta"}) {
        const Json* block = get(meta, ch);
        if (!block || !block->is_object()) continue;
        const auto n = rows.find(ch);
        if (n == rows.end()) {
            errors.push_back(std::string("unrestricted NTOs need mo_coefficients_") + ch + " (spec § 5.3)");
        } else {
            bool mismatch = false;
            for (const char* key : {"energies", "occupations"}) {
                const Json* a = get(*block, key);
                if (a && a->is_array() && static_cast<int64_t>(a->size()) != n->second) mismatch = true;
            }
            const Json* nm = get(*block, "n_mo");
            if (is_integral(nm) && nm->as_double() != n->second) mismatch = true;
            if (mismatch) errors.push_back(std::string("NTO ") + ch + " metadata must have one entry per coefficient row (spec § 5.3)");
        }
        const Json* occ = get(*block, "occupations");
        if (occ && occ->is_array()) {
            bool numbers = true, bad = false;
            double previous = 0;
            for (size_t i = 0; i < occ->size(); ++i) {
                const Json& x = occ->at(i);
                if (!x.is_number()) { numbers = false; break; }
                const double v = x.as_double();
                if (v < 0 || (i && previous < v)) bad = true;
                previous = v;
            }
            if (numbers && bad) errors.push_back(std::string("NTO ") + ch + " weights must be nonnegative and non-increasing (spec § 5.3)");
        }
    }
    return errors;
}

inline std::vector<std::pair<std::string, std::string>> nto_pair_errors(
    const DeferredSection& d, const std::vector<DeferredSection>& ntos,
    const std::map<std::string, std::string>& kinds,
    const std::map<std::string, StateIndices>& indices) {
    std::vector<std::pair<std::string, std::string>> errors;
    const Json& meta = d.doc;
    if (!unrestricted_nto(meta)) return errors;
    const Json* ref = get(meta, "state_ref");
    const Json* index = get(meta, "state_index");
    const Json* role = get(meta, "nto_role");
    if (!ref || !ref->is_string() || !is_integral(index) || !role || !role->is_string() ||
        (role->as_string() != "hole" && role->as_string() != "particle")) return errors;
    auto kind = kinds.find(ref->as_string());
    auto states = indices.find(ref->as_string());
    if (kind == kinds.end() || kind->second != "excited_states" ||
        (states != indices.end() && states->second.known && std::none_of(states->second.values.begin(), states->second.values.end(),
            [index](int64_t v) { return static_cast<double>(v) == index->as_double(); })))
        errors.emplace_back("E_REF_UNRESOLVED", "NTO state_ref/state_index must name an excited state (spec § 5.3)");
    const char* opposite = role->as_string() == "hole" ? "particle" : "hole";
    std::vector<const DeferredSection*> partners;
    for (const auto& p : ntos) {
        const Json* pr = get(p.doc, "state_ref");
        const Json* pi = get(p.doc, "state_index");
        if (nto_text(p.doc, "orbital_kind", "nto") && pr && pr->is_string() &&
            pr->as_string() == ref->as_string() && is_integral(pi) && pi->as_double() == index->as_double() &&
            nto_text(p.doc, "nto_role", opposite)) partners.push_back(&p);
    }
    if (partners.size() != 1) {
        errors.emplace_back("E_REF_UNRESOLVED", "NTOs must have one opposite-role section for the same state (spec § 5.3)");
        return errors;
    }
    const auto& partner = *partners.front();
    const Json& other = partner.doc;
    if (!nto_text(other, "spin", "unrestricted")) {
        errors.emplace_back("E_MEMBER_MALFORMED", "NTO partners must have the same spin layout (spec § 5.3)");
        return errors;
    }
    const MoRows rows = mo_rows(get(*d.section, "members"));
    const MoRows other_rows = mo_rows(get(*partner.section, "members"));
    for (const char* ch : {"alpha", "beta"}) {
        auto a = rows.find(ch), b = other_rows.find(ch);
        if (a != rows.end() && b != other_rows.end() && a->second != b->second)
            errors.emplace_back("E_MEMBER_MALFORMED", std::string("NTO partners must have the same ") + ch + " row count (spec § 5.3)");
        const Json* ab = get(meta, ch); const Json* bb = get(other, ch);
        const Json* av = ab ? get(*ab, "occupations") : nullptr;
        const Json* bv = bb ? get(*bb, "occupations") : nullptr;
        if (av && bv && av->is_array() && bv->is_array()) {
            const auto numbers = [](const Json& a) {
                return std::all_of(a.elements().begin(), a.elements().end(),
                                   [](const Json& x) { return x.is_number(); });
            };
            if (!numbers(*av) || !numbers(*bv)) continue;
            bool equal = av->size() == bv->size();
            if (equal) for (size_t i = 0; i < av->size(); ++i) {
                const Json& x = av->at(i); const Json& y = bv->at(i);
                if (x.is_number() && y.is_number()) equal = equal && x.as_double() == y.as_double();
                else equal = equal && x.dump() == y.dump();
            }
            if (!equal) errors.emplace_back("E_MEMBER_MALFORMED", std::string("NTO partners must have identical ") + ch + " weights (spec § 5.3)");
        }
    }
    return errors;
}

// Cross-member rules of a wavefunction.gto section (spec § 5.3): coefficient
// rank against `complex`, symmetry_labels lengths, [n_ao, n_ao] shapes and the
// pairing of the AO-basis matrices, basis.n_electrons against the occupations
// (E_MEMBER_MALFORMED) and against provenance.n_electrons
// (E_ROOT_INCONSISTENT). Mirrors _wavefunction_errors() in python/qvf_reader.py.
inline std::vector<RoleFinding> wavefunction_errors(const std::string& sid,
                                                    const Json& members,
                                                    const std::map<std::string, Json>& docs,
                                                    const Json* provenance) {
    std::vector<RoleFinding> errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    auto doc = [&](const char* role) -> const Json* {
        auto it = docs.find(role);
        return it == docs.end() ? nullptr : &it->second;
    };
    const Json* basis = doc("basis");
    const Json* meta = doc("mo_metadata");
    std::vector<int64_t> shape;

    const int cplx = coefficients_complex(meta);
    if (cplx >= 0)
        for (const char* role : {"mo_coefficients", "mo_coefficients_alpha",
                                 "mo_coefficients_beta"}) {
            if (!shape_of(get(members, role), shape)) continue;
            if (cplx == 1 && (shape.size() != 3 || shape[2] != 2))
                errs.push_back({"E_MEMBER_MALFORMED", role,
                                sec + role + " shape " + dims_str(shape) +
                                    " has no trailing [re, im] axis of 2, but "
                                    "mo_metadata declares complex coefficients "
                                    "(spec § 5.3)"});
            else if (cplx == 0 && shape.size() != 2)
                errs.push_back({"E_MEMBER_MALFORMED", role,
                                sec + role + " shape " + dims_str(shape) +
                                    " is not [n_mo, n_ao], and mo_metadata does not "
                                    "declare complex coefficients (spec § 5.3)"});
        }

    if (meta && meta->is_object()) {
        std::vector<std::pair<std::string, const Json*>> sets{{"", meta}};
        for (const char* ch : {"alpha", "beta"})
            if (const Json* c = get(*meta, ch))
                if (c->is_object()) sets.emplace_back(std::string(ch) + ".", c);
        for (const auto& set : sets) {
            const Json* labels = get(*set.second, "symmetry_labels");
            const Json* energies = get(*set.second, "energies");
            if (labels && energies && labels->is_array() && energies->is_array() &&
                labels->size() != energies->size())
                errs.push_back({"E_MEMBER_MALFORMED", "mo_metadata",
                                sec + "mo_metadata " + set.first + "symmetry_labels has " +
                                    std::to_string(labels->size()) + " entries but " +
                                    set.first + "energies has " +
                                    std::to_string(energies->size()) +
                                    "; there is one label per MO (spec § 5.3)"});
        }
    }

    int64_t n_ao = 0;
    const bool have_n_ao = basis_n_ao(basis, n_ao);
    for (const char* role : {"overlap", "density_matrix", "density_matrix_alpha",
                             "density_matrix_beta"})
        if (have_n_ao && shape_of(get(members, role), shape) &&
            (shape.size() != 2 || shape[0] != n_ao || shape[1] != n_ao))
            errs.push_back({"E_MEMBER_MALFORMED", role,
                            sec + role + " shape " + dims_str(shape) +
                                " is not [n_ao, n_ao] = [" + std::to_string(n_ao) + ", " +
                                std::to_string(n_ao) + "] for the basis (spec § 5.3)"});
    const bool unrestricted =
        members.contains("mo_coefficients_alpha") && members.contains("mo_coefficients_beta");
    if (members.contains("density_matrix") && !members.contains("mo_coefficients"))
        errs.push_back({"E_MEMBER_MALFORMED", "density_matrix",
                        sec + "density_matrix accompanies mo_coefficients only; an "
                              "unrestricted wavefunction carries density_matrix_alpha "
                              "and density_matrix_beta (spec § 5.3)"});
    for (const auto& pair : {std::make_pair("density_matrix_alpha", "density_matrix_beta"),
                             std::make_pair("density_matrix_beta", "density_matrix_alpha")}) {
        const std::string role = pair.first, other = pair.second;
        if (!members.contains(role)) continue;
        if (!unrestricted)
            errs.push_back({"E_MEMBER_MALFORMED", role,
                            sec + role + " accompanies mo_coefficients_alpha and "
                                         "mo_coefficients_beta only; a single-set "
                                         "wavefunction carries density_matrix (spec § 5.3)"});
        else if (!members.contains(other))
            errs.push_back({"E_MEMBER_MALFORMED", role,
                            sec + role + " without " + other +
                                "; the two come together (spec § 5.3)"});
    }

    const Json* n_el = basis ? get(*basis, "n_electrons") : nullptr;
    if (!n_el || !is_plain_number(*n_el)) return errs;
    const double n = n_el->as_double();
    const Json* semantics = meta ? get(*meta, "occupation_semantics") : nullptr;
    if (meta && meta->is_object() &&
        !(semantics && semantics->is_string() &&
          semantics->as_string() == "transition_weight")) {
        std::vector<const Json*> lists;
        if (const Json* occ = get(*meta, "occupations")) lists.push_back(occ);
        for (const char* ch : {"alpha", "beta"})
            if (const Json* c = get(*meta, ch))
                if (const Json* occ = get(*c, "occupations")) lists.push_back(occ);
        bool numbers = !lists.empty();
        double total = 0.0;
        for (const Json* occ : lists) {
            if (!occ->is_array()) { numbers = false; break; }
            for (const Json& x : occ->elements()) {
                if (!is_plain_number(x)) { numbers = false; break; }
                total += x.as_double();
            }
        }
        if (numbers && !(std::fabs(total - n) <= kNElectronsTolerance))
            errs.push_back({"E_MEMBER_MALFORMED", "basis",
                            sec + "basis.n_electrons is " + show(*n_el) +
                                " but the mo_metadata occupations sum to " + g10(total) +
                                "; they must agree within 1e-6 (spec § 5.3)"});
    }
    const Json* prov_n = provenance ? get(*provenance, "n_electrons") : nullptr;
    if (prov_n && is_plain_number(*prov_n) && prov_n->as_double() != n)
        errs.push_back({"E_ROOT_INCONSISTENT", "basis",
                        sec + "basis.n_electrons is " + show(*n_el) +
                            " but provenance.n_electrons is " + show(*prov_n) +
                            "; they state one electron count and must be equal "
                            "(spec § 3.2, § 5.3)"});
    return errs;
}

// Three booleans, or false.
inline bool bool3_of(const Json* j, std::array<bool, 3>& out) {
    if (!j || !j->is_array() || j->size() != 3) return false;
    for (size_t i = 0; i < 3; ++i) {
        if (!j->at(i).is_bool()) return false;
        out[i] = j->at(i).as_bool();
    }
    return true;
}

inline std::string vec_str(const std::array<double, 3>& v) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), "[%.6g, %.6g, %.6g]", v[0], v[1], v[2]);
    return buf;
}

// A grid payload, checked after every structure has been read.
struct GridMember {
    std::string where, sid;
    const Json* section;
    Json grid;
};

// The unit name spec § 4.3 fixes for the grid units.data of a volume kind, or
// nullptr. One grid payload serves every volume.* kind, so the schema cannot
// fix it per kind.
inline const char* volume_data_unit(const std::string& kind) {
    if (kind == "volume.density" || kind == "volume.transition_density") return "e_per_bohr3";
    if (kind == "volume.potential") return "hartree_per_e";
    if (kind == "volume.rdg") return "dimensionless";
    if (kind == "volume.laplacian") return "e_per_bohr5";
    return nullptr;
}

// A grid units.data that is a unit name other than the one § 4.3 fixes for the
// section's volume kind (E_UNIT_CONTRADICTED), as a message, or empty. A string
// that is not a unit name is W_UNIT_UNKNOWN alone. Mirrors
// _grid_data_unit_error() in python/qvf/validate.py.
inline std::string grid_data_unit_error(const std::string& sid, const Json& section,
                                        const Json& grid) {
    const Json* kind = get(section, "kind");
    const char* fixed = kind && kind->is_string() ? volume_data_unit(kind->as_string()) : nullptr;
    const Json* units = grid.is_object() ? get(grid, "units") : nullptr;
    const Json* unit = units && units->is_object() ? get(*units, "data") : nullptr;
    if (!fixed || !unit || !unit->is_string()) return {};
    const std::string& name = unit->as_string();
    if (name == fixed || !parse_unit(name)) return {};
    return "section " + quote_str(sid) + ": grid units.data " + quote_str(name) +
           " contradicts the unit " + quote_str(fixed) + " that spec § 4.3 fixes for " +
           kind->as_string() + " data";
}

// Reference and periodic-span rules of a grid payload: a grid structure_ref
// names a structure section; the grid's structure is the section-level
// structure_ref, else the grid's, else the archive's only structure. When it
// resolves to a structure whose payload was read (`cells`: id -> payload),
// every axis a with periodic[a] true has pbc[a] true and a span
// voxel_vectors[a] * (shape[a] - 1) (includes_endpoint true, the default) or
// voxel_vectors[a] * shape[a] (false) equal to lattice row a in bohr within
// 1e-6 relative. Malformed fields are left to the schema. Mirrors
// _grid_errors() in python/qvf_reader.py. Returns (code, message) pairs.
inline std::vector<std::pair<const char*, std::string>> grid_errors(
    const std::string& sid, const Json& section, const Json& grid,
    const std::map<std::string, std::string>& kinds_by_id,
    const std::map<std::string, Json>& cells, const std::string* only_structure) {
    std::vector<std::pair<const char*, std::string>> errs;
    if (!grid.is_object()) return errs;
    const std::string sec = "section " + quote_str(sid) + ": ";
    const Json* ref = get(grid, "structure_ref");
    if (ref) {
        auto it = ref->is_string() ? kinds_by_id.find(ref->as_string()) : kinds_by_id.end();
        if (it == kinds_by_id.end() || it->second != "structure")
            errs.emplace_back("E_REF_UNRESOLVED", sec + "grid structure_ref " + show(*ref) +
                                                      " does not name a structure section "
                                                      "(spec § 5.2)");
    }
    std::array<bool, 3> periodic{};
    if (!bool3_of(get(grid, "periodic"), periodic) ||
        !(periodic[0] || periodic[1] || periodic[2]))
        return errs;
    std::string structure;
    if (const Json* sref = get(section, "structure_ref")) {
        if (!sref->is_string()) return errs;
        structure = sref->as_string();
    } else if (ref && ref->is_string()) {
        structure = ref->as_string();
    } else if (only_structure) {
        structure = *only_structure;
    } else {
        return errs;
    }
    auto cell = cells.find(structure);
    if (cell == cells.end()) return errs;
    std::array<bool, 3> pbc{{false, false, false}};
    const Json* pj = get(cell->second, "pbc");
    if (pj && !pj->is_null() && !bool3_of(pj, pbc)) return errs;
    bool endpoint = true;
    if (const Json* e = get(grid, "includes_endpoint")) {
        if (!e->is_bool()) return errs;
        endpoint = e->as_bool();
    }
    const Json* shape = get(grid, "shape");
    const Json* vox = get(grid, "voxel_vectors");
    if (!shape || !shape->is_array() || shape->size() != 3 || !vox || !vox->is_array() ||
        vox->size() != 3)
        return errs;
    for (size_t a = 0; a < 3; ++a)
        if (shape->at(a).type() != Json::Type::Int) return errs;
    const Json* lattice = get(cell->second, "lattice_vectors");
    for (size_t a = 0; a < 3; ++a) {
        if (!periodic[a]) continue;
        const std::string axis = std::to_string(a);
        if (!pbc[a]) {
            errs.emplace_back("E_MEMBER_MALFORMED",
                              sec + "grid axis " + axis + " is periodic but structure " +
                                  quote_str(structure) + " has pbc[" + axis +
                                  "] false; a grid axis can be periodic only along a "
                                  "periodic axis of its structure (spec § 5.2)");
            continue;
        }
        std::array<double, 3> step{}, row{};
        if (!vec3_of(&vox->at(a), step) || !lattice || !lattice->is_array() ||
            lattice->size() != 3 || !vec3_of(&lattice->at(a), row))
            continue;
        const int64_t n = shape->at(a).as_int() - (endpoint ? 1 : 0);
        std::array<double, 3> span{}, target{};
        double diff = 0.0, norm = 0.0;
        for (size_t c = 0; c < 3; ++c) {
            span[c] = step[c] * static_cast<double>(n);
            target[c] = row[c] * kBohrPerAngstrom;   // spec § 4.3, Lattice conventions
            diff += (span[c] - target[c]) * (span[c] - target[c]);
            norm += target[c] * target[c];
        }
        if (std::sqrt(diff) > 1e-6 * std::sqrt(norm))
            errs.emplace_back("E_MEMBER_MALFORMED",
                              sec + "periodic grid axis " + axis + " spans " +
                                  vec_str(span) + " bohr (voxel_vectors[" + axis +
                                  "] x " + std::to_string(n) + ", includes_endpoint " +
                                  (endpoint ? "true" : "false") + ") but lattice vector " +
                                  axis + " of structure " + quote_str(structure) + " is " +
                                  vec_str(target) +
                                  " bohr; a periodic grid spans its cell within 1e-6 "
                                  "relative (spec § 5.2)");
    }
    return errs;
}

// --- preview, reaction.scheme, chemical identity (spec § 5.15, § 5.14, § 5.1) --

// The bytes of a `preview` member are what its `media_type` says (spec
// § 5.15): the PNG signature, FF D8 FF for JPEG, or a UTF-8 document with an
// `<svg` element. Mirrors _preview_member_errors() in python/qvf_reader.py.
inline std::vector<CodeFinding> preview_member_errors(const std::string& sid,
                                                      const std::string& role,
                                                      const Json& member,
                                                      const std::vector<uint8_t>& raw) {
    std::vector<CodeFinding> out;
    const Json* mt = get(member, "media_type");
    if (!mt || !mt->is_string()) return out;
    const std::string& media_type = mt->as_string();
    const std::string where = "section " + quote_str(sid) + " member " + quote_str(role);
    auto begins_with = [&raw](const char* signature, size_t n) {
        return raw.size() >= n && std::memcmp(raw.data(), signature, n) == 0;
    };
    if (media_type == "image/png" && !begins_with("\x89PNG\r\n\x1a\n", 8)) {
        out.push_back({"E_MEMBER_MALFORMED", where + ": image/png bytes do not begin with "
                                                     "the PNG signature (spec § 5.15)"});
    } else if (media_type == "image/jpeg" && !begins_with("\xff\xd8\xff", 3)) {
        out.push_back({"E_MEMBER_MALFORMED",
                       where + ": image/jpeg bytes do not begin with FF D8 FF (spec § 5.15)"});
    } else if (media_type == "image/svg+xml") {
        static const char kSvg[] = "<svg";
        if (!valid_utf8(raw))
            out.push_back({"E_UTF8", where + ": image/svg+xml bytes do not decode as "
                                             "UTF-8 (spec § 5.15)"});
        else if (std::search(raw.begin(), raw.end(), kSvg, kSvg + 4) == raw.end())
            out.push_back({"E_MEMBER_MALFORMED", where + ": image/svg+xml bytes contain no "
                                                         "<svg element (spec § 5.15)"});
    }
    return out;
}

// Reaction energy quantities and the term field of each absolute value
// (spec § 5.14, *Reaction energies*).
constexpr std::array<std::array<const char*, 2>, 3> kReactionQuantities{{
    {{"energy", "energy_eh"}},
    {{"enthalpy", "enthalpy_eh"}},
    {{"gibbs_free_energy", "gibbs_free_energy_eh"}},
}};
constexpr double kHartreeToKcalMol = 627.509474;  // CODATA 2018
constexpr double kReactionDeltaTolEh = 1e-8;
constexpr double kReactionKcalRelTol = 1e-6;
constexpr double kReactionCoefficientTol = 1e-9;

inline bool finite_number(const Json* j) {
    return j && j->is_number() && std::isfinite(j->as_double());
}

// Sum with one final rounding, as Python's math.fsum (Shewchuk partials).
inline double exact_sum(const std::vector<double>& values) {
    std::vector<double> partials;
    for (double x : values) {
        size_t i = 0;
        for (double y : partials) {
            if (std::fabs(x) < std::fabs(y)) std::swap(x, y);
            const double hi = x + y;
            const double lo = y - (hi - x);
            if (lo != 0.0) partials[i++] = lo;
            x = hi;
        }
        partials.resize(i);
        partials.push_back(x);
    }
    double total = 0.0;
    for (double p : partials) total += p;
    return total;
}

// A step's `energies` (spec § 5.14, *Reaction energies*): no energies on a
// resonance arrow; one term per reactant and product with the signed
// coefficient; each delta equals the sum of its terms and its kcal/mol value
// the conversion; one method (E_MEMBER_MALFORMED); every `run_ref` names a
// `run.record` (E_REF_UNRESOLVED). Mirrors _reaction_energies_errors() in
// python/qvf/validate.py.
inline std::vector<CodeFinding> reaction_energies_errors(
    const std::string& where, const Json& step, const Json& energies,
    const std::map<std::string, std::string>& kinds_by_id) {
    std::vector<CodeFinding> out;
    auto bad = [&](const std::string& msg) {
        out.push_back({"E_MEMBER_MALFORMED", where + " energies: " + msg + " (spec § 5.14)"});
    };
    const Json* arrow = get(step, "arrow");
    if (arrow && arrow->is_string() && arrow->as_string() == "resonance")
        bad("a resonance arrow has no reaction energy");
    const Json* list = get(energies, "terms");
    if (!list || !list->is_array()) return out;
    std::map<std::string, double> expected;
    for (const auto& role : {std::make_pair("reactants", -1.0), std::make_pair("products", 1.0)}) {
        const Json* items = get(step, role.first);
        if (!items || !items->is_array()) continue;
        for (const Json& item : items->elements()) {
            const Json* ref = get(item, "structure_ref");
            if (!ref || !ref->is_string()) continue;
            const Json* c = get(item, "coefficient");
            if (c && !finite_number(c)) continue;
            expected[ref->as_string()] += role.second * (c ? c->as_double() : 1.0);
        }
    }
    for (auto it = expected.begin(); it != expected.end();)
        it = it->second == 0.0 ? expected.erase(it) : std::next(it);
    std::vector<const Json*> terms;
    for (const Json& t : list->elements())
        if (t.is_object()) terms.push_back(&t);
    std::set<std::string> seen;
    for (size_t t = 0; t < terms.size(); ++t) {
        const Json* ref = get(*terms[t], "structure_ref");
        const Json* c = get(*terms[t], "coefficient");
        if (!ref || !ref->is_string()) continue;
        const std::string r = ref->as_string();
        if (!seen.insert(r).second) {
            bad("structure " + quote_str(r) + " has more than one term");
            continue;
        }
        auto want = expected.find(r);
        if (want == expected.end()) {
            bad("term " + std::to_string(t) + " names " + quote_str(r) +
                ", which is not a reactant or product");
        } else if (finite_number(c) &&
                   std::fabs(c->as_double() - want->second) >
                       kReactionCoefficientTol * std::max(1.0, std::fabs(c->as_double()))) {
            bad("term " + std::to_string(t) + " coefficient " + format_number(c->as_double()) +
                " for " + quote_str(r) + ", but the step gives " + format_number(want->second));
        }
        const Json* run_ref = get(*terms[t], "run_ref");
        if (run_ref && run_ref->is_string()) {
            auto kind = kinds_by_id.find(run_ref->as_string());
            if (kind == kinds_by_id.end() || kind->second != "run.record")
                out.push_back({"E_REF_UNRESOLVED",
                               where + " energies: term " + std::to_string(t) + " run_ref " +
                                   quote_str(run_ref->as_string()) +
                                   " does not name a run.record section (spec § 5.14)"});
        }
    }
    for (const auto& e : expected)
        if (!seen.count(e.first)) bad("no term for participant " + quote_str(e.first));
    const Json* first_method = nullptr;
    bool methods_differ = false;
    for (const Json* term : terms) {
        const Json* m = get(*term, "method");
        if (!m || !m->is_object()) continue;
        if (!first_method) first_method = m;
        else if (!schema::json_equal(*m, *first_method)) methods_differ = true;
    }
    if (methods_differ) bad("the terms have different methods");
    const Json* deltas = get(energies, "deltas");
    if (!deltas || !deltas->is_object()) return out;
    for (const auto& q : kReactionQuantities) {
        const Json* delta = get(*deltas, q[0]);
        const Json* value = delta ? get(*delta, "value_eh") : nullptr;
        if (!delta || !delta->is_object() || !finite_number(value)) continue;
        const double v = value->as_double();
        std::vector<double> products;
        long missing = -1;
        bool coefficients = true;
        for (size_t t = 0; t < terms.size(); ++t) {
            const Json* a = get(*terms[t], q[1]);
            const Json* c = get(*terms[t], "coefficient");
            if (!finite_number(a)) {
                if (missing < 0) missing = static_cast<long>(t);
                continue;
            }
            if (!finite_number(c)) coefficients = false;
            else products.push_back(c->as_double() * a->as_double());
        }
        if (missing >= 0) {
            bad(std::string(q[0]) + " is given but term " + std::to_string(missing) +
                " has no " + q[1]);
        } else if (coefficients) {
            const double total = exact_sum(products);
            if (std::fabs(total - v) > kReactionDeltaTolEh)
                bad(std::string(q[0]) + " value_eh " + format_number(v) +
                    " is not the sum of its terms, " + format_number(total));
        }
        const Json* kcal = get(*delta, "kcal_mol");
        if (finite_number(kcal)) {
            const double converted = v * kHartreeToKcalMol;
            if (std::fabs(kcal->as_double() - converted) >
                kReactionKcalRelTol * std::max(1.0, std::fabs(converted)))
                bad(std::string(q[0]) + " kcal_mol " + format_number(kcal->as_double()) +
                    " is not value_eh x 627.509474 = " + format_number(converted));
        }
    }
    return out;
}

// A `reaction.scheme` payload (spec § 5.14): every participant names a
// `structure` section (E_REF_UNRESOLVED), the step ids are unique and no
// structure repeats within one role of a step (E_MEMBER_MALFORMED). Mirrors
// _reaction_scheme_errors() in python/qvf_reader.py.
inline std::vector<CodeFinding> reaction_scheme_errors(
    const std::string& sid, const Json& payload,
    const std::map<std::string, std::string>& kinds_by_id) {
    std::vector<CodeFinding> out;
    const Json* steps = get(payload, "steps");
    if (!steps || !steps->is_array()) return out;
    const std::string sec = "section " + quote_str(sid) + ": ";
    std::set<std::string> seen;
    const std::vector<Json>& list = steps->elements();
    for (size_t k = 0; k < list.size(); ++k) {
        const Json& step = list[k];
        if (!step.is_object()) continue;
        const Json* id = get(step, "id");
        if (id && id->is_string()) {
            if (seen.count(id->as_string()))
                out.push_back({"E_MEMBER_MALFORMED",
                               sec + "step id " + quote_str(id->as_string()) +
                                   " is used by more than one step (spec § 5.14)"});
            seen.insert(id->as_string());
        }
        for (const char* role : {"reactants", "products", "reagents"}) {
            const Json* items = get(step, role);
            if (!items || !items->is_array()) continue;
            std::map<std::string, size_t> uses;   // sorted, as Python's sorted()
            for (const Json& item : items->elements()) {
                const Json* ref = get(item, "structure_ref");
                if (!ref || !ref->is_string()) continue;
                ++uses[ref->as_string()];
                auto it = kinds_by_id.find(ref->as_string());
                if (it == kinds_by_id.end() || it->second != "structure")
                    out.push_back({"E_REF_UNRESOLVED",
                                   sec + "step " + std::to_string(k) + " " + role +
                                       " structure_ref " + quote_str(ref->as_string()) +
                                       " does not name a structure section (spec § 5.14)"});
            }
            for (const auto& u : uses)
                if (u.second > 1)
                    out.push_back({"E_MEMBER_MALFORMED",
                                   sec + "step " + std::to_string(k) + " names " +
                                       quote_str(u.first) + " more than once among its " +
                                       role + " (spec § 5.14)"});
        }
        if (const Json* energies = get(step, "energies"); energies && energies->is_object())
            for (auto& e : reaction_energies_errors(sec + "step " + std::to_string(k), step,
                                                    *energies, kinds_by_id))
                out.push_back(std::move(e));
    }
    return out;
}

// The chemical identity of a `structure` section (spec § 5.1): an integer
// `charge` leaves a non-negative electron count N = sum_z - charge for the
// atoms of its payload (`sum_z`, nullptr when an atomic number is missing or
// not a non-negative integer), and an integer `multiplicity` >= 1 has
// multiplicity - 1 <= N with the parity of N (E_CHARGE_MULTIPLICITY); in an
// archive with exactly one structure section both equal the integer
// provenance.charge / provenance.multiplicity (E_ROOT_INCONSISTENT). The
// findings name the section field. Mirrors _identity_errors() in
// python/qvf_reader.py; the arithmetic is exact for every int64 value.
inline std::vector<FieldFinding> identity_errors(const std::string& sid, const Json& section,
                                                 const int64_t* sum_z, const Json* provenance,
                                                 bool only_structure) {
    std::vector<FieldFinding> out;
    const Json* cj = get(section, "charge");
    const Json* mj = get(section, "multiplicity");
    const bool has_charge = is_int(cj);
    const bool has_mult = is_int(mj) && mj->as_int() >= 1;
    const int64_t charge = has_charge ? cj->as_int() : 0;
    const int64_t mult = has_mult ? mj->as_int() : 0;
    const std::string sec = "structure " + quote_str(sid) + ": ";
    if (has_charge && sum_z) {
        const int64_t z = *sum_z;   // >= 0
        if (charge > z) {
            out.push_back({"E_CHARGE_MULTIPLICITY", "charge",
                           sec + "charge " + std::to_string(charge) + " leaves " +
                               std::to_string(z - charge) +
                               " electrons for atoms whose atomic numbers sum to " +
                               std::to_string(z) + " (spec § 5.1)"});
        } else if (has_mult) {
            // N = z - charge >= 0, exact in uint64 even for a very negative charge.
            const uint64_t electrons = static_cast<uint64_t>(z) - static_cast<uint64_t>(charge);
            const uint64_t unpaired = static_cast<uint64_t>(mult - 1);
            if (unpaired > electrons || ((electrons ^ unpaired) & 1u))
                out.push_back({"E_CHARGE_MULTIPLICITY", "multiplicity",
                               sec + "multiplicity " + std::to_string(mult) +
                                   " is impossible with " + std::to_string(electrons) +
                                   " electrons (charge " + std::to_string(charge) +
                                   "): the unpaired electrons must not exceed the electron "
                                   "count and must have its parity (spec § 5.1)"});
        }
    }
    if (only_structure && provenance && provenance->is_object()) {
        const struct {
            const char* field;
            bool has;
            int64_t value;
        } fields[] = {{"charge", has_charge, charge}, {"multiplicity", has_mult, mult}};
        for (const auto& fv : fields) {
            const Json* stated = get(*provenance, fv.field);
            if (fv.has && is_int(stated) && stated->as_int() != fv.value)
                out.push_back({"E_ROOT_INCONSISTENT", fv.field,
                               sec + fv.field + " " + std::to_string(fv.value) +
                                   " differs from provenance." + fv.field + " " +
                                   std::to_string(stated->as_int()) +
                                   ", which describes the archive's only structure "
                                   "(spec § 5.1)"});
        }
    }
    return out;
}

// The whole validator. Never throws for malformed input: every problem is a
// finding. Mirrors _validate_open_zip() in python/qvf_reader.py check for check.
inline void run_validation(const ZipReader& zip, const Limits& lim,
                           const ValidationOptions& opts, std::vector<Finding>& out,
                           validate::Level level = validate::Level::All,
                           size_t* n_sections = nullptr) {
    const bool payloads = level >= validate::Level::Payloads;
    const bool semantics = level >= validate::Level::Semantics;
    Sink f(out, level >= validate::Level::References);
    const std::vector<ZipReader::Entry>& entries = zip.entries();
    const std::vector<std::string> names = zip.names();
    std::map<std::string, size_t> counts;
    for (const auto& e : entries) ++counts[e.name];

    // --- container: central-directory checks, before decompressing anything
    for (const auto& name : names) {
        const size_t n = counts[name];
        if (n > 1 && name == "manifest.json") {
            // An archive updated in place: the last manifest.json is the
            // manifest, every earlier one is superseded (spec § 2.3).
            for (size_t copy = 1; copy < n; ++copy)
                f.warning("W_ORPHAN_MEMBER", "entry:" + name,
                          "ZIP entry " + quote_str(name) + " (" + std::to_string(copy) +
                              " of " + std::to_string(n) + ") is superseded by the last " +
                              quote_str(name) +
                              " entry (archive updated in place; consumers ignore it)");
        } else if (n > 1) {
            f.error("E_ZIP_DUPLICATE_ENTRY", "entry:" + name,
                    "ZIP entry name " + quote_str(name) + " occurs " + std::to_string(n) +
                        " times; readers disagree on which entry wins");
        }
    }
    for (const auto& name : names) {
        std::string why = path_problem(name);
        if (!why.empty())
            f.error("E_MEMBER_PATH", "entry:" + name,
                    "ZIP entry " + quote_str(name) + ": illegal path (" + why + ")");
    }
    std::set<std::string> unsafe;  // entries that must not be decompressed
    uint64_t total = 0;
    for (const auto& e : entries) {
        total += e.usize;
        if (e.usize > lim.max_member_bytes) {
            f.error("E_MEMBER_TOO_LARGE", "entry:" + e.name,
                    "ZIP entry " + quote_str(e.name) +
                        " exceeds the per-member size cap (" +
                        std::to_string(e.usize) + " > " +
                        std::to_string(lim.max_member_bytes) +
                        " bytes uncompressed)");
            unsafe.insert(e.name);
        }
        if (lim.max_ratio > 0 && e.method != 0 && e.usize >= lim.ratio_exempt_bytes &&
            (e.csize == 0 ||
             static_cast<double>(e.usize) > lim.max_ratio * static_cast<double>(e.csize))) {
            std::string ratio =
                e.csize ? format_number(static_cast<double>(e.usize) / e.csize) : "inf";
            f.error("E_COMPRESSION_RATIO", "entry:" + e.name,
                    "ZIP entry " + quote_str(e.name) + " expands " +
                        std::to_string(e.usize) + " bytes from " +
                        std::to_string(e.csize) + " (ratio " + ratio + " > " +
                        format_number(lim.max_ratio) + ")");
            unsafe.insert(e.name);
        }
    }
    const bool payloads_blocked = total > lim.max_total_bytes;
    if (payloads_blocked)
        f.error("E_TOTAL_TOO_LARGE", "archive",
                "archive expands to " + std::to_string(total) +
                    " bytes uncompressed, over the " +
                    std::to_string(lim.max_total_bytes) +
                    "-byte total cap; member payloads were not read");

    // --- manifest (find() resolves a duplicated name to its last central-
    // directory record: the current manifest of an updated archive, § 2.3)
    const ZipReader::Entry* me = zip.find("manifest.json");
    if (!me) {
        f.error("E_MANIFEST_MISSING", "manifest.json", "missing manifest.json");
        return;
    }
    if (level == validate::Level::Container) return;
    if (me->usize > lim.max_json_bytes) {
        f.error("E_JSON_TOO_LARGE", "manifest.json",
                "manifest.json is " + std::to_string(me->usize) +
                    " bytes, over the " + std::to_string(lim.max_json_bytes) +
                    "-byte JSON cap");
        return;
    }
    if (unsafe.count("manifest.json")) return;
    std::vector<uint8_t> raw;
    try {
        raw = zip.read("manifest.json", lim.max_member_bytes);
    } catch (const std::exception& ex) {
        f.error("E_ZIP_UNREADABLE", "entry:manifest.json",
                std::string("ZIP entry 'manifest.json' cannot be read: ") + ex.what());
        return;
    }
    Json manifest;
    try {
        manifest = Json::parse(std::string(raw.begin(), raw.end()), lim.max_json_depth);
    } catch (const JsonDepthError&) {
        f.error("E_JSON_TOO_DEEP", "manifest.json",
                "manifest.json nests deeper than " +
                    std::to_string(lim.max_json_depth) + " levels");
        return;
    } catch (const std::exception& ex) {
        f.error("E_MANIFEST_JSON", "manifest.json",
                std::string("manifest.json is not valid JSON: ") + ex.what());
        return;
    }
    if (!valid_utf8(raw)) {
        f.error("E_MANIFEST_JSON", "manifest.json",
                "manifest.json is not valid JSON: invalid UTF-8");
        return;
    }
    if (!manifest.is_object()) {
        f.error("E_MANIFEST_JSON", "manifest.json", "manifest.json must be a JSON object");
        return;
    }

    // Root sanity. The withdrawn version 2 is read as 1 (spec § 3.1): every
    // check sees version 1, the JSON Schema check included.
    const Json* version = get(manifest, "qvf_version");
    const Json* checked_manifest = &manifest;
    Json v1_manifest;
    if (version && version->is_number() && version->as_double() == 2.0) {
        f.warning("W_VERSION_WITHDRAWN", "/qvf_version",
                  "qvf_version 2 was withdrawn: the archive is checked as version 1, "
                  "which producers must write (spec § 3.1)");
        v1_manifest = manifest;
        v1_manifest.set("qvf_version", 1);
        checked_manifest = &v1_manifest;
    } else if (!version || !version->is_number() || version->as_double() != 1.0) {
        f.error("E_VERSION", "/qvf_version",
                "qvf_version must be 1, got " + (version ? show(*version) : "nothing"));
    }
    const Json* src = get(manifest, "source");
    if (src && !src->is_object()) {
        f.error("E_SOURCE", "/source", "source must be a JSON object");
        src = nullptr;
    }
    for (const char* key : {"program", "version", "calculation"}) {
        if (!src || !src->contains(key))
            f.error("E_SOURCE", "/source", std::string("source.") + key + " is missing");
    }
    auto non_empty = [&](const char* key) {
        const Json* v = src ? get(*src, key) : nullptr;
        return v && truthy(*v);
    };
    if (!non_empty("program") || !non_empty("version"))
        f.error("E_SOURCE", "/source", "source.program and source.version must be non-empty");
    if (semantics) for (const auto& e : root_consistency_errors(manifest))
        f.error("E_ROOT_INCONSISTENT", e.first, e.second);

    // The JSON Schema, through the checker generated from it.
    if (opts.schema_check) {
        for (const schema::Error& e : schema::check_manifest(*checked_manifest))
            f.error("E_SCHEMA", e.pointer.empty() ? "manifest.json" : e.pointer,
                    (e.location.empty() ? "(root)" : e.location) + ": " + e.message);
    } else {
        f.warning("W_SCHEMA_SKIPPED", "manifest.json",
                  "the JSON Schema check was disabled; structural schema check skipped");
    }

    // Declared units of {value, units} objects in provenance (spec § 4.3).
    if (semantics) if (const Json* prov = get(manifest, "provenance")) {
        if (prov->is_object())
            for (const auto& kv : prov->items()) {
                const Json& v = kv.second;
                if (!v.is_object() || !v.contains("value") || !v.contains("units"))
                    continue;
                const Json& u = v.at("units");
                if (!(u.is_string() && canonical_unit(u.as_string())))
                    f.warning("W_UNIT_UNKNOWN",
                              "/provenance/" + escape_pointer(kv.first) + "/units",
                              "provenance." + kv.first + ".units " + show(u) +
                                  " is not a unit name of spec § 4.3 or a legacy "
                                  "spelling");
            }
    }

    std::vector<const Json*> sections;
    if (const Json* secs = get(manifest, "sections")) {
        if (!secs->is_array())
            f.error("E_SECTION_MALFORMED", "/sections", "sections must be an array");
        else
            for (const Json& s : secs->elements()) sections.push_back(&s);
    }
    if (n_sections) *n_sections = sections.size();
    std::set<std::string> section_ids;
    std::map<std::string, std::string> kinds_by_id;  // first occurrence wins
    std::vector<std::string> structure_ids;
    // Section id -> kind ("" when not a string), for reference resolution
    // (spec § 4.1, § 5.5). With duplicate ids the first section wins.
    std::map<std::string, std::string> section_kinds;
    for (const Json* s : sections) {
        const Json* id = get(*s, "id");
        if (id && id->is_string()) {
            section_ids.insert(id->as_string());
            const Json* k = get(*s, "kind");
            if (k && k->is_string()) {
                kinds_by_id.emplace(id->as_string(), k->as_string());
                if (k->as_string() == "structure") structure_ids.push_back(id->as_string());
            }
            section_kinds.emplace(id->as_string(),
                                  k && k->is_string() ? k->as_string() : std::string());
        }
    }
    const std::string* only_structure =
        structure_ids.size() == 1 ? &structure_ids.front() : nullptr;
    std::map<std::string, int64_t> n_atoms_of;  // structure id -> atom count
    std::map<std::string, int64_t> sum_z_of;    // structure id -> sum of atomic numbers
    std::vector<DeferredSection> identity_sections;   // structure sections, checked last
    std::map<std::string, bool> periodic_of;    // structure id -> any(pbc)
    std::vector<DerivativeSection> derivative_sections;  // checked last
    std::vector<DeferredSection> bonds_sections;         // standalone bonds, checked last
    // Excited states (spec § 5.12): orbital rows of every wavefunction.gto,
    // state indices of every excited_states payload, and the sections whose
    // references are checked once every payload has been read.
    std::map<std::string, MoRows> mo_rows_of;
    for (const Json* s : sections) {
        const Json* id = get(*s, "id");
        const Json* k = get(*s, "kind");
        if (id && id->is_string() && k && k->is_string() &&
            k->as_string() == "wavefunction.gto" && !mo_rows_of.count(id->as_string()))
            mo_rows_of.emplace(id->as_string(), mo_rows(get(*s, "members")));
    }
    std::map<std::string, StateIndices> state_indices_of;
    std::vector<DeferredSection> nto_sections;
    std::vector<DeferredSection> excited_sections, states_ref_spectra,
        transition_densities, wavefunction_volumes;
    std::map<std::string, Json> cells;          // structure id -> payload
    std::vector<GridMember> grid_members;       // grid payloads, checked last
    std::vector<PropertySection> property_sections;      // checked last
    std::vector<BlochSection> bloch_sections;            // checked last
    const Json* declared_ext = get(manifest, "extensions");
    if (declared_ext && !declared_ext->is_object()) declared_ext = nullptr;
    std::set<std::string> declared_paths;
    std::set<std::string> seen_ids;

    for (size_t i = 0; i < sections.size(); ++i) {
        const Json& s = *sections[i];
        const std::string here = section_pointer(i);
        if (!s.is_object()) {
            f.error("E_SECTION_MALFORMED", here,
                    "section " + std::to_string(i) + " must be a JSON object");
            continue;
        }
        std::string sid;
        const Json* idj = get(s, "id");
        if (!idj || !idj->is_string()) {
            f.error("E_SECTION_MALFORMED", here,
                    "section " + std::to_string(i) + ": id must be a string, got " +
                        (idj ? show(*idj) : "nothing"));
            sid = "<section " + std::to_string(i) + ">";
        } else {
            sid = idj->as_string();
            if (seen_ids.count(sid))
                f.error("E_SECTION_ID_DUPLICATE", here + "/id",
                        "duplicate section id " + quote_str(sid));
            seen_ids.insert(sid);
            if (sid == "auto_open" || sid == "bookmarks" || sid == "camera" ||
                sid == "crossfade")
                f.error("E_RESERVED_ID", here + "/id",
                        "section id " + quote_str(sid) +
                            " is reserved (a viewer_defaults key, spec § 3.1, § 3.4)");
        }

        std::string kind;
        const Json* kj = get(s, "kind");
        if (!kj || !kj->is_string()) {
            f.error("E_SECTION_MALFORMED", here,
                    "section " + quote_str(sid) + ": kind must be a string, got " +
                        (kj ? show(*kj) : "nothing"));
        } else {
            kind = kj->as_string();
            if (!is_canonical_kind(kind) && !starts_with(kind, "x_"))
                f.error("E_KIND_UNKNOWN", here + "/kind",
                        "section " + quote_str(sid) + ": unknown non-vendor kind " +
                            quote_str(kind));
        }

        // critical governance (spec § 7.4).
        const Json* cj = get(s, "critical");
        if (cj && truthy(*cj) && starts_with(kind, "x_")) {
            std::string vendor = kind.substr(0, kind.find('.'));
            if (!declared_ext || !declared_ext->contains(vendor))
                f.error("E_CRITICAL_UNDECLARED", here + "/critical",
                        "section " + quote_str(sid) + " is critical but its namespace " +
                            quote_str(vendor) + " is not declared in root extensions");
        }

        const Json* members = get(s, "members");
        if (!members || !members->is_object()) {
            f.error("E_SECTION_MALFORMED", here,
                    "section " + quote_str(sid) + ": members must be a JSON object");
            members = nullptr;
        }

        std::map<std::string, Json> docs;  // parsed JSON members, by role
        std::vector<uint8_t> overlap_bytes;  // a wavefunction.gto overlap member
        bool have_overlap = false;
        std::vector<uint8_t> energies_bytes;  // trajectory-family energies member
        bool have_energies_bytes = false;
        if (members) {
            for (const auto& kv : members->items()) {
                const std::string& role = kv.first;
                const Json& m = kv.second;
                const std::string mwhere = member_pointer(i, role);
                const Json* pj = get(m, "path");
                if (!pj || !pj->is_string()) {
                    f.error("E_MEMBER_MALFORMED", mwhere,
                            "section " + quote_str(sid) + " member " + quote_str(role) +
                                ": malformed spec");
                    continue;
                }
                const std::string& mpath = pj->as_string();
                declared_paths.insert(mpath);
                std::string why = path_problem(mpath);
                if (!why.empty()) {
                    f.error("E_MEMBER_PATH", mwhere,
                            "member " + quote_str(mpath) + ": illegal path (" + why + ")");
                    continue;
                }
                const ZipReader::Entry* e = zip.find(mpath);
                if (!e) {
                    f.error("E_MEMBER_MISSING", mwhere,
                            "member " + quote_str(mpath) + " declared but not in archive");
                    continue;
                }
                // Declared size against the central directory, before reading.
                std::string size_msg;
                if (const char* code = size_problem(m, *e, size_msg))
                    f.error(code, mwhere, size_msg);
                if (!payloads || payloads_blocked || unsafe.count(mpath)) continue;
                const Json* fj = get(m, "format");
                const std::string fmt = fj && fj->is_string() ? fj->as_string() : "";
                if (fmt == "json" && e->usize > lim.max_json_bytes) {
                    f.error("E_JSON_TOO_LARGE", mwhere,
                            "member " + quote_str(mpath) + ": " + std::to_string(e->usize) +
                                " bytes, over the " + std::to_string(lim.max_json_bytes) +
                                "-byte JSON cap");
                    continue;
                }
                std::vector<uint8_t> bytes;
                try {
                    bytes = zip.read(mpath, lim.max_member_bytes);
                } catch (const std::exception& ex) {
                    f.error("E_ZIP_UNREADABLE", "entry:" + mpath,
                            "ZIP entry " + quote_str(mpath) + " cannot be read: " + ex.what());
                    continue;
                }
                if (const Json* hj = get(m, "sha256")) {
                    if (!hj->is_string() || hj->as_string() != Sha256::hex_of(bytes))
                        f.error("E_SHA256_MISMATCH", mwhere,
                                "member " + quote_str(mpath) + ": sha256 mismatch");
                }
                if (fmt == "json") {
                    Json doc;
                    bool parsed = false;
                    try {
                        doc = Json::parse(std::string(bytes.begin(), bytes.end()),
                                          lim.max_json_depth);
                        parsed = true;
                    } catch (const JsonDepthError&) {
                        f.error("E_JSON_TOO_DEEP", mwhere,
                                "member " + quote_str(mpath) + ": JSON nests deeper than " +
                                    std::to_string(lim.max_json_depth) + " levels");
                        continue;
                    } catch (const std::exception& ex) {
                        f.error("E_JSON_INVALID", mwhere,
                                "member " + quote_str(mpath) + ": not valid JSON (" +
                                    ex.what() + ")");
                    }
                    if (parsed && !valid_utf8(bytes)) {
                        f.error("E_JSON_INVALID", mwhere,
                                "member " + quote_str(mpath) +
                                    ": not valid JSON (invalid UTF-8)");
                        parsed = false;
                    }
                    if (parsed) {
                        if (semantics && kind == "structure" && role == "structure") {
                            for (const auto& msg : structure_payload_errors(sid, doc))
                                f.error("E_STRUCTURE_PERIODICITY", mwhere, msg);
                            const Json* atoms = get(doc, "atoms");
                            if (atoms && atoms->is_array()) {
                                n_atoms_of.emplace(sid,
                                                   static_cast<int64_t>(atoms->size()));
                                // Σ Z, when every atomic_number is a
                                // non-negative integer (spec § 5.1).
                                int64_t sum_z = 0;
                                bool known = true;
                                for (const Json& a : atoms->elements()) {
                                    const Json* z = get(a, "atomic_number");
                                    if (!is_int(z) || z->as_int() < 0 ||
                                        sum_z > INT64_MAX - z->as_int()) {
                                        known = false;
                                        break;
                                    }
                                    sum_z += z->as_int();
                                }
                                if (known) sum_z_of.emplace(sid, sum_z);
                            }
                            if (doc.is_object()) cells.emplace(sid, doc);
                            if (doc.is_object()) {
                                const Json* pbc = get(doc, "pbc");
                                bool periodic = false;
                                if (pbc && pbc->is_array())
                                    for (const Json& x : pbc->elements())
                                        periodic = periodic || (x.is_bool() && x.as_bool());
                                periodic_of.emplace(sid, periodic);
                            }
                        }
                        if (semantics && is_grid_payload(kind, role) && doc.is_object())
                            grid_members.push_back(GridMember{mwhere, sid, &s, doc});
                        // The member's payload definition (spec § 5, § 6).
                        const std::string_view def = gen::payload_def(kind, role);
                        if (opts.schema_check && !def.empty())
                            for (const schema::Error& e : schema::check_definition(def, doc))
                                f.error("E_PAYLOAD_SCHEMA", mwhere,
                                        "schema: " + mpath + ": " + e.location + ": " +
                                            e.message);
                        if (semantics && is_canonical_kind(kind))
                            for (const auto& du : declared_units(kind, role, doc))
                                if (!(du.second->is_string() &&
                                      parse_unit(du.second->as_string())))
                                    f.warning("W_UNIT_UNKNOWN", mwhere,
                                              "member " + quote_str(mpath) +
                                                  ": declared unit " + show(*du.second) +
                                                  " at " + du.first +
                                                  " is not a unit name of spec § 4.3");
                        if (semantics) docs[role] = std::move(doc);
                    }
                } else if (fmt == "binary" && m.contains("dtype")) {
                    const Json& dj = m.at("dtype");
                    size_t itemsize = dj.is_string() ? dtype_itemsize(dj.as_string()) : 0;
                    const Json* shape = get(m, "shape");
                    bool shape_ok = true;
                    if (shape) {
                        shape_ok = shape->is_array();
                        if (shape_ok)
                            for (const Json& d : shape->elements())
                                if (d.type() != Json::Type::Int || d.as_int() < 0)
                                    shape_ok = false;
                    }
                    if (!itemsize) {
                        f.error("E_DTYPE", mwhere,
                                "member " + quote_str(mpath) + ": bad dtype " + show(dj));
                    } else if (!shape_ok) {
                        f.error("E_MEMBER_MALFORMED", mwhere,
                                "member " + quote_str(mpath) +
                                    ": shape must be a list of non-negative integers, "
                                    "got " + show(*shape));
                    } else {
                        uint64_t n = 1;
                        bool overflow = false;
                        if (shape)
                            for (const Json& d : shape->elements()) {
                                uint64_t k = static_cast<uint64_t>(d.as_int());
                                if (k && n > UINT64_MAX / k) overflow = true;
                                else n *= k;
                            }
                        if (!overflow && n > UINT64_MAX / itemsize) overflow = true;
                        uint64_t want = overflow ? 0 : n * itemsize;
                        if (overflow || bytes.size() != want)
                            f.error("E_BINARY_SIZE", mwhere,
                                    "member " + quote_str(mpath) + ": byte length " +
                                        std::to_string(bytes.size()) +
                                        " != itemsize*prod(shape) " +
                                        (overflow ? std::string("(overflow)")
                                                  : std::to_string(want)));
                        else if (semantics && kind == "wavefunction.gto" && role == "overlap") {
                            overlap_bytes = bytes;
                            have_overlap = true;
                        }
                    }
                }
                if (semantics && is_trajectory_kind(kind) && role == "energies" && fmt == "binary") {
                    energies_bytes = bytes;
                    have_energies_bytes = true;
                }
                if (kind == "run.record" && (role == "input" || role == "log") &&
                    !valid_utf8(bytes))
                    f.error("E_UTF8", mwhere,
                            "member " + quote_str(mpath) + ": run.record " + role +
                                " must decode as UTF-8");
                if (kind == "preview" && fmt == "binary")
                    for (const auto& e : preview_member_errors(sid, role, m, bytes))
                        f.error(e.code, mwhere, e.message);
            }
        }

        if (!semantics) continue;

        // 1-D spectra: one shared array length (spec § 5.4).
        if (is_spectrum_1d(kind)) {
            auto it = docs.find("spectrum");
            for (const auto& msg :
                 spectrum_payload_errors(sid, it != docs.end() ? &it->second : nullptr))
                f.error("E_PAYLOAD_SCHEMA", member_pointer(i, "spectrum"), msg);
        }

        // Cross-references.
        if (kind == "volume.difference") {
            for (const char* op : {"operand_a", "operand_b"}) {
                const Json* ref = get(s, op);
                if (ref && !ref->is_null() &&
                    (!ref->is_string() || !section_ids.count(ref->as_string())))
                    f.error("E_REF_UNRESOLVED", here + "/" + op,
                            "section " + quote_str(sid) + ": " + op + " " + show(*ref) +
                                " does not resolve");
            }
        }
        if (kind == "reaction.waypoints") {
            const Json* ref = get(s, "trajectory_ref");
            if (!ref || !ref->is_string() || !section_ids.count(ref->as_string()))
                f.error("E_REF_UNRESOLVED", here + "/trajectory_ref",
                        "section " + quote_str(sid) + ": trajectory_ref " +
                            (ref ? show(*ref) : "nothing") + " does not resolve");
            else if (section_kinds.at(ref->as_string()) != "trajectory")
                f.error("E_REF_UNRESOLVED", here + "/trajectory_ref",
                        "section " + quote_str(sid) + ": trajectory_ref " +
                            show(*ref) + " names a section of kind " +
                            quote_str(section_kinds.at(ref->as_string())) +
                            ", not 'trajectory' (spec § 5.5)");
        }
        // Generic section references, on any kind (spec § 4.1).
        static const struct {
            const char* field;
            const char* target;
        } kSectionRefs[] = {{"structure_ref", "structure"}, {"run_ref", "run.record"}};
        for (const auto& r : kSectionRefs) {
            const Json* ref = get(s, r.field);
            if (!ref || ref->is_null()) continue;
            const std::string where = here + "/" + r.field;
            if (!ref->is_string() || !section_kinds.count(ref->as_string()))
                f.error("E_REF_UNRESOLVED", where,
                        "section " + quote_str(sid) + ": " + r.field + " " + show(*ref) +
                            " does not resolve");
            else if (section_kinds.at(ref->as_string()) != r.target)
                f.error("E_REF_UNRESOLVED", where,
                        "section " + quote_str(sid) + ": " + r.field + " " + show(*ref) +
                            " names a section of kind " +
                            quote_str(section_kinds.at(ref->as_string())) + ", not " +
                            quote_str(r.target) + " (spec § 4.1)");
        }
        if (kind == "bands" && members)
            for (const auto& e : bands_shape_errors(sid, *members))
                f.error("E_MEMBER_MALFORMED", member_pointer(i, e.first), e.second);
        if (is_trajectory_kind(kind) && members) {
            auto it = docs.find("metadata");
            for (const auto& e : trajectory_errors(
                     sid, kind, *members, it != docs.end() ? &it->second : nullptr,
                     have_energies_bytes ? &energies_bytes : nullptr))
                f.error("E_MEMBER_MALFORMED", member_pointer(i, e.first), e.second);
        }
        if (kind == "vibrations") {
            if (const Json* ref = get(s, "hessian_ref")) {
                auto it = ref->is_string() ? kinds_by_id.find(ref->as_string())
                                           : kinds_by_id.end();
                if (it == kinds_by_id.end() || it->second != "hessian")
                    f.error("E_REF_UNRESOLVED", here + "/hessian_ref",
                            "section " + quote_str(sid) + ": hessian_ref " + show(*ref) +
                                " does not name a hessian section (spec § 5.5)");
            }
        }
        if (kind == "job.spec") {
            // The payload's structure_ref names a structure section (spec § 5.9).
            auto it = docs.find("spec");
            if (it != docs.end()) {
                if (const Json* ref = get(it->second, "structure_ref")) {
                    auto kind_it = ref->is_string() ? kinds_by_id.find(ref->as_string())
                                                    : kinds_by_id.end();
                    if (kind_it == kinds_by_id.end() || kind_it->second != "structure")
                        f.error("E_REF_UNRESOLVED", member_pointer(i, "spec"),
                                "section " + quote_str(sid) + ": spec.structure_ref " +
                                    show(*ref) +
                                    " does not name a structure section (spec § 5.9)");
                }
            }
        }
        if (kind == "forces" || kind == "hessian") {
            auto it = docs.find("meta");
            const bool has_meta = it != docs.end();
            derivative_sections.push_back(DerivativeSection{
                i, sid, kind, &s, members, has_meta, has_meta ? it->second : Json()});
        }
        if (kind == "bonds" && s.has("structure_ref")) {
            auto it = docs.find("bonds");
            const bool has_doc = it != docs.end();
            bonds_sections.push_back(
                DeferredSection{i, sid, &s, has_doc, has_doc ? it->second : Json()});
        }
        if (kind == "excited_states") {
            auto it = docs.find("states");
            const bool has_doc = it != docs.end();
            if (!state_indices_of.count(sid))
                state_indices_of.emplace(sid,
                                         state_indices(has_doc ? &it->second : nullptr));
            excited_sections.push_back(DeferredSection{
                i, sid, &s, has_doc, has_doc ? it->second : Json()});
        }
        if (is_spectrum_1d(kind)) {
            auto it = docs.find("spectrum");
            if (it != docs.end() && it->second.is_object() &&
                it->second.contains("states_ref"))
                states_ref_spectra.push_back(DeferredSection{i, sid, &s, true, it->second});
        }
        if (kind == "volume.transition_density")
            transition_densities.push_back(DeferredSection{i, sid, &s, false, Json()});
        if (is_wavefunction_volume_kind(kind) && s.contains("wavefunction_ref"))
            wavefunction_volumes.push_back(DeferredSection{i, sid, &s, false, Json()});
        for (const auto& e : volume_errors(sid, kind, members, docs, kinds_by_id))
            f.error(e.code, member_pointer(i, e.role), e.message);
        if (is_property_kind(kind)) {
            auto it = docs.find("properties");
            if (it != docs.end())
                property_sections.push_back(PropertySection{i, sid, kind, &s, it->second});
        }
        if (kind == "wavefunction.bloch")
            bloch_sections.push_back(BlochSection{i, sid, &s, members, docs});
        if (kind == "atom_properties" && members)
            for (const auto& e : atom_properties_meta_errors(sid, *members, docs))
                f.error("E_MEMBER_MALFORMED", member_pointer(i, e.first), e.second);
        if (kind == "wavefunction.gto" && members) {
            auto meta = docs.find("mo_metadata");
            if (meta != docs.end()) {
                nto_sections.push_back(DeferredSection{i, sid, &s, true, meta->second});
                for (const auto& msg : unrestricted_nto_errors(meta->second, members))
                    f.error("E_MEMBER_MALFORMED", member_pointer(i, "mo_metadata"), msg);
            }
            for (const auto& e :
                 wavefunction_errors(sid, *members, docs, get(manifest, "provenance")))
                f.error(e.code, member_pointer(i, e.role), e.message);
            auto basis = docs.find("basis");
            for (const auto& msg : normalization_errors(
                     sid, basis == docs.end() ? nullptr : &basis->second,
                     get(*members, "overlap"), have_overlap ? &overlap_bytes : nullptr))
                f.error("E_NORMALIZATION", member_pointer(i, "overlap"), msg);
        }
        if (kind == "reaction.scheme") {
            auto it = docs.find("scheme");
            if (it != docs.end() && it->second.is_object())
                for (const auto& e : reaction_scheme_errors(sid, it->second, kinds_by_id))
                    f.error(e.code, member_pointer(i, "scheme"), e.message);
        }
        if (kind == "run.record") {
            // job_ref names a job.spec section (spec § 5.8).
            if (const Json* ref = get(s, "job_ref")) {
                auto it = ref->is_string() ? kinds_by_id.find(ref->as_string())
                                           : kinds_by_id.end();
                if (it == kinds_by_id.end() || it->second != "job.spec")
                    f.error("E_REF_UNRESOLVED", here + "/job_ref",
                            "section " + quote_str(sid) + ": job_ref " + show(*ref) +
                                " does not name a job.spec section (spec § 5.8)");
            }
        }
        if (kind == "structure")
            identity_sections.push_back(DeferredSection{i, sid, &s, false, Json()});
        if (kind == "run.record") {
            auto it = docs.find("files");
            if (it != docs.end() && it->second.is_object())
                for (const auto& kv : it->second.items())
                    if (!members || !members->contains(kv.first))
                        f.error("E_FILES_INDEX", member_pointer(i, "files"),
                                "section " + quote_str(sid) + ": files index names role " +
                                    quote_str(kv.first) + " with no matching member");
        }
    }

    if (semantics) {
    // Chemical identity of the structures (spec § 5.1).
    for (const auto& d : identity_sections) {
        auto z = sum_z_of.find(d.sid);
        for (const auto& e : identity_errors(d.sid, *d.section,
                                             z != sum_z_of.end() ? &z->second : nullptr,
                                             get(manifest, "provenance"),
                                             structure_ids.size() == 1))
            f.error(e.code, section_pointer(d.index) + "/" + e.field, e.message);
    }

    // forces / hessian: shapes against the structure, which may come later.
    for (const auto& d : derivative_sections)
        for (const auto& e : derivative_errors(d.sid, d.kind, *d.section, d.members,
                                               d.has_meta ? &d.meta : nullptr,
                                               n_atoms_of, kinds_by_id, only_structure))
            f.error(e.code, member_pointer(d.index, e.role), e.message);

    // Standalone bonds: atom indices against the structure they name, which
    // may come later (spec § 5.1).
    for (const auto& d : bonds_sections)
        for (const auto& msg : bonds_errors(d.sid, *d.section, d.has_doc ? &d.doc : nullptr,
                                            kinds_by_id, n_atoms_of))
            f.error("E_MEMBER_MALFORMED", member_pointer(d.index, "bonds"), msg);

    for (const auto& d : nto_sections)
        for (const auto& e : nto_pair_errors(d, nto_sections, kinds_by_id, state_indices_of))
            f.error(e.first.c_str(), member_pointer(d.index, "mo_metadata"), e.second);

    // Excited states, the spectra that name them, and transition densities
    // (spec § 5.2, § 5.4, § 5.12): the sections they name may come later.
    for (const auto& d : excited_sections)
        for (const auto& e : excited_states_errors(d.sid, *d.section,
                                                   d.has_doc ? &d.doc : nullptr,
                                                   kinds_by_id, mo_rows_of))
            f.error(e.code,
                    e.field.empty() ? member_pointer(d.index, "states")
                                    : section_pointer(d.index) + "/" + e.field,
                    e.message);
    for (const auto& d : states_ref_spectra)
        for (const auto& msg : spectrum_states_errors(d.sid, d.doc, kinds_by_id,
                                                      state_indices_of))
            f.error("E_REF_UNRESOLVED", member_pointer(d.index, "spectrum"), msg);
    for (const auto& d : transition_densities)
        for (const auto& e : transition_density_errors(d.sid, *d.section, kinds_by_id,
                                                       state_indices_of))
            f.error("E_REF_UNRESOLVED", section_pointer(d.index) + "/" + e.first, e.second);
    // Volumes that name their wavefunction and orbital (spec § 5.2).
    for (const auto& d : wavefunction_volumes) {
        const Json* k = get(*d.section, "kind");
        for (const auto& e : wavefunction_volume_errors(d.sid, k ? k->as_string() : "",
                                                        *d.section, kinds_by_id, mo_rows_of))
            f.error("E_REF_UNRESOLVED", section_pointer(d.index) + "/" + e.first, e.second);
    }

    // Grids: references and periodic spans against the structure (spec § 5.2).
    for (const auto& g : grid_members) {
        for (const auto& e : grid_errors(g.sid, *g.section, g.grid, kinds_by_id, cells,
                                         only_structure))
            f.error(e.first, g.where, e.second);
        const std::string unit_msg = grid_data_unit_error(g.sid, *g.section, g.grid);
        if (!unit_msg.empty()) f.error("E_UNIT_CONTRADICTED", g.where, unit_msg);
    }

    // properties.*: rules JSON Schema cannot express (spec § 5.13); the Born
    // charges are counted against a structure that may come later.
    for (const auto& p : property_sections)
        for (const auto& e : property_errors(p.sid, p.kind, *p.section, p.payload, manifest,
                                             n_atoms_of, only_structure))
            f.error(e.code, member_pointer(p.index, "properties"), e.message);

    // wavefunction.bloch: shapes, and a periodic structure, which may come later.
    for (const auto& b : bloch_sections)
        for (const auto& e : bloch_errors(b.sid, *b.section, b.members, b.docs, kinds_by_id,
                                          periodic_of))
            f.error(e.code, section_pointer(b.index) + "/" + e.role, e.message);

    }  // semantics

    // --- orphans: entries no manifest member references (spec § 2.1)
    for (const auto& name : names) {
        if (name != "manifest.json" && !(name.size() && name.back() == '/') &&
            !declared_paths.count(name))
            f.warning("W_ORPHAN_MEMBER", "entry:" + name,
                      "ZIP entry " + quote_str(name) +
                          " is not referenced by any manifest member (orphan; "
                          "consumers ignore it)");
    }
}

}  // namespace vld
}  // namespace detail

namespace validate {
Report run(const container::ZipReader& archive, Level level, const Limits& limits,
           const ValidationOptions& options) {
    if (level < Level::Container || level > Level::All)
        throw std::invalid_argument("qvf::validate: invalid validation level");
    Report report;
    detail::vld::run_validation(archive, limits, options, report.findings,
                                level, &report.n_sections);
    report.ok = !has_errors(report.findings);
    return report;
}
}  // namespace validate

bool has_errors(const std::vector<Finding>& findings) {
    for (const auto& f : findings)
        if (f.severity == Finding::Severity::Error) return true;
    return false;
}

struct QvfReader::Impl {
    detail::ZipReader zip;     // container::ZipReader, over memory or a file
    Limits limits;
    Json manifest;
    Impl(detail::ZipReader z, const Limits& lim)
        : zip(std::move(z)), limits(lim) {
        // Superseded manifest.json copies of an archive updated in place are
        // not ambiguous: the last one is the manifest (spec § 2.3).
        for (const std::string& name : zip.duplicates())
            if (name != "manifest.json")
                throw std::runtime_error("E_ZIP_DUPLICATE_ENTRY: duplicate ZIP entry name '" +
                                         name + "'; refusing an ambiguous archive");
        manifest = parse_json_entry("manifest.json");
    }

    // Read and parse a JSON entry under the JSON size and depth caps.
    Json parse_json_entry(const std::string& path) const {
        std::vector<uint8_t> b = read_json_bytes(path);
        return Json::parse(std::string(b.begin(), b.end()), limits.max_json_depth);
    }

    std::vector<uint8_t> read_json_bytes(const std::string& path) const {
        const detail::ZipReader::Entry* e = zip.find(path);
        if (!e) throw std::runtime_error("zip: member not found: " + path);
        if (e->usize > limits.max_json_bytes)
            throw std::runtime_error("E_JSON_TOO_LARGE: " + path + " is " +
                                     std::to_string(e->usize) + " bytes, over the " +
                                     std::to_string(limits.max_json_bytes) +
                                     "-byte JSON cap");
        return zip.read(path, limits.max_json_bytes);
    }

    const Json& member_spec(const std::string& sid, const std::string& role) const {
        const Json& secs = manifest.at("sections");
        for (size_t i = 0; i < secs.size(); ++i) {
            const Json& s = secs.at(i);
            if (s.contains("id") && s.at("id").as_string() == sid) {
                const Json& members = s.at("members");
                if (!members.contains(role))
                    throw std::runtime_error("section '" + sid +
                                             "' has no member '" + role + "'");
                return members.at(role);
            }
        }
        throw std::runtime_error("section id not found: " + sid);
    }

    // A member's bytes under every cap of `limits` (member, total, ratio).
    std::vector<uint8_t> read_spec(const Json& spec, bool verify) const {
        std::vector<uint8_t> bytes = zip.read(spec.at("path").as_string(), limits);
        if (verify && spec.contains("sha256")) {
            std::string got = detail::Sha256::hex_of(bytes);
            if (got != spec.at("sha256").as_string())
                throw std::runtime_error("sha256 mismatch for " +
                                         spec.at("path").as_string());
        }
        return bytes;
    }
};

std::vector<double> BinaryMember::as_doubles() const {
    if (dtype != "float64")
        throw std::runtime_error("BinaryMember::as_doubles: dtype is " + dtype);
    std::vector<double> out(element_count());
    if (bytes.size() != out.size() * 8)
        throw std::runtime_error("BinaryMember::as_doubles: size mismatch");
    for (size_t i = 0; i < out.size(); ++i) {
        uint64_t bits = 0;
        for (int b = 0; b < 8; ++b)
            bits |= static_cast<uint64_t>(bytes[i * 8 + static_cast<size_t>(b)])
                    << (8 * b);
        std::memcpy(&out[i], &bits, 8);  // little-endian host assumed (spec § 2.2)
    }
    return out;
}

QvfReader::QvfReader(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
QvfReader::QvfReader(QvfReader&&) noexcept = default;
QvfReader& QvfReader::operator=(QvfReader&&) noexcept = default;
QvfReader::~QvfReader() = default;

// "little" is the only v1 byte order; refuse rather than guess (spec § 2.2).
// (A schema violation, hence the E_SCHEMA prefix, as in the Python reader.)
static QvfReader checked_byte_order(QvfReader reader) {
    const Json* order = detail::vld::get(reader.manifest(), "byte_order");
    if (order && !(order->is_string() && order->as_string() == "little"))
        throw std::runtime_error("E_SCHEMA: byte_order " + detail::vld::show(*order) +
                                 " is not 'little', the only byte order of QVF v1 "
                                 "(spec § 2.2); refusing to guess");
    return reader;
}

QvfReader QvfReader::from_bytes(std::vector<uint8_t> archive, const Limits& limits) {
    QvfReader reader = checked_byte_order(QvfReader(std::unique_ptr<Impl>(
        new Impl(detail::ZipReader(std::move(archive)), limits))));
    if (reader.manifest().is_object())
        reader.typed_ =
            std::make_shared<const Manifest>(Manifest::from_json(reader.manifest()));
    return reader;
}

// Refuse a member whose declared `size` disagrees with its ZIP entry (spec
// § 4.2), before reading it. Other problems of the spec are left to the read.
static void check_declared_size(const detail::ZipReader& zip, const Json& spec) {
    const Json* path = detail::vld::get(spec, "path");
    if (!path || !path->is_string()) return;
    const detail::ZipReader::Entry* e = zip.find(path->as_string());
    if (!e) return;
    std::string message;
    if (const char* code = detail::vld::size_problem(spec, *e, message))
        throw std::runtime_error(std::string(code) + ": " + message);
}

// The archive file stays open; entries are read from it on demand.
QvfReader QvfReader::open(const std::string& path, const Limits& limits) {
    QvfReader reader = checked_byte_order(
        QvfReader(std::unique_ptr<Impl>(new Impl(detail::ZipReader::open(path), limits))));
    if (reader.manifest().is_object())
        reader.typed_ =
            std::make_shared<const Manifest>(Manifest::from_json(reader.manifest()));
    return reader;
}

const Json& QvfReader::manifest() const { return impl_->manifest; }

const Manifest& QvfReader::typed_manifest() const {
    if (!typed_) throw std::runtime_error("manifest.json is not a JSON object");
    return *typed_;
}

const Json& QvfReader::sections() const { return impl_->manifest.at("sections"); }

bool QvfReader::has_section(const std::string& id) const {
    const Json& secs = sections();
    for (size_t i = 0; i < secs.size(); ++i)
        if (secs.at(i).contains("id") && secs.at(i).at("id").as_string() == id)
            return true;
    return false;
}

const Json& QvfReader::section(const std::string& id) const {
    const Json& secs = sections();
    for (size_t i = 0; i < secs.size(); ++i)
        if (secs.at(i).contains("id") && secs.at(i).at("id").as_string() == id)
            return secs.at(i);
    throw std::runtime_error("section id not found: " + id);
}

std::vector<uint8_t> QvfReader::read_member_bytes(const std::string& sid,
                                                  const std::string& role,
                                                  bool verify) const {
    const Json& spec = impl_->member_spec(sid, role);
    check_declared_size(impl_->zip, spec);
    return impl_->read_spec(spec, verify);
}

Json QvfReader::read_json_member(const std::string& sid, const std::string& role,
                                 bool verify) const {
    const Json& spec = impl_->member_spec(sid, role);
    check_declared_size(impl_->zip, spec);
    const detail::ZipReader::Entry* e = impl_->zip.find(spec.at("path").as_string());
    if (e && e->usize > impl_->limits.max_json_bytes)
        throw std::runtime_error("E_JSON_TOO_LARGE: " + spec.at("path").as_string() +
                                 " is over the JSON size cap");
    std::vector<uint8_t> b = impl_->read_spec(spec, verify);
    return Json::parse(std::string(b.begin(), b.end()), impl_->limits.max_json_depth);
}

BinaryMember QvfReader::read_binary_member(const std::string& sid,
                                           const std::string& role,
                                           bool verify) const {
    const Json& spec = impl_->member_spec(sid, role);
    BinaryMember bm;
    if (spec.contains("dtype")) bm.dtype = spec.at("dtype").as_string();
    if (spec.contains("shape")) {
        const Json& sh = spec.at("shape");
        for (size_t i = 0; i < sh.size(); ++i)
            bm.shape.push_back(sh.at(i).as_int());
    }
    check_declared_size(impl_->zip, spec);
    bm.bytes = impl_->read_spec(spec, verify);
    return bm;
}

std::vector<Finding> QvfReader::validate_findings(const Limits& limits,
                                                 const ValidationOptions& options) const {
    std::vector<Finding> out;
    detail::vld::run_validation(impl_->zip, limits, options, out);
    return out;
}

std::vector<Finding> QvfReader::validate_bytes(std::vector<uint8_t> archive,
                                               const Limits& limits,
                                               const ValidationOptions& options) {
    std::vector<Finding> out;
    std::unique_ptr<detail::ZipReader> zip;
    try {
        zip.reset(new detail::ZipReader(std::move(archive)));
    } catch (const std::exception& ex) {
        out.push_back(Finding{Finding::Severity::Error, "E_ZIP_UNREADABLE", "archive",
                              std::string("not a readable zip: ") + ex.what()});
        return out;
    }
    detail::vld::run_validation(*zip, limits, options, out);
    return out;
}

// Reads the archive from the file on demand rather than loading it whole.
std::vector<Finding> QvfReader::validate_file(const std::string& path,
                                              const Limits& limits,
                                              const ValidationOptions& options) {
    {
        std::ifstream probe(path, std::ios::binary);
        if (!probe)
            return {Finding{Finding::Severity::Error, "E_ZIP_UNREADABLE", "archive",
                            "cannot open file: " + path}};
    }
    std::vector<Finding> out;
    std::unique_ptr<detail::ZipReader> zip;
    try {
        zip.reset(new detail::ZipReader(detail::ZipReader::open(path)));
    } catch (const std::exception& ex) {
        out.push_back(Finding{Finding::Severity::Error, "E_ZIP_UNREADABLE", "archive",
                              std::string("not a readable zip: ") + ex.what()});
        return out;
    }
    detail::vld::run_validation(*zip, limits, options, out);
    return out;
}

bool QvfReader::validate(std::vector<std::string>& errors) const {
    bool ok = true;
    for (const auto& f : validate_findings()) {
        if (f.severity != Finding::Severity::Error) continue;
        errors.push_back(f.code + ": " + f.message);
        ok = false;
    }
    return ok;
}

}  // namespace qvf
