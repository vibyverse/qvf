// qvf.cpp — implementation of the QVF writer library. License: Apache-2.0.
#include "qvf/qvf.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <stdexcept>

#include "qvf/container.hpp"
#include "qvf_sha256.hpp"

namespace qvf {

// ---------------------------------------------------------------------------
// dtype helpers
// ---------------------------------------------------------------------------

static const char* dtype_name(DType d) {
    switch (d) {
        case DType::Int8: return "int8";
        case DType::Int16: return "int16";
        case DType::Int32: return "int32";
        case DType::Int64: return "int64";
        case DType::UInt8: return "uint8";
        case DType::UInt16: return "uint16";
        case DType::UInt32: return "uint32";
        case DType::UInt64: return "uint64";
        case DType::Float32: return "float32";
        case DType::Float64: return "float64";
    }
    return "float64";
}

static size_t dtype_itemsize(DType d) {
    switch (d) {
        case DType::Int8: case DType::UInt8: return 1;
        case DType::Int16: case DType::UInt16: return 2;
        case DType::Int32: case DType::UInt32: case DType::Float32: return 4;
        case DType::Int64: case DType::UInt64: case DType::Float64: return 8;
    }
    return 8;
}

static void put_le(std::vector<uint8_t>& out, uint64_t v, int nbytes) {
    for (int i = 0; i < nbytes; ++i)
        out.push_back(static_cast<uint8_t>((v >> (8 * i)) & 0xFF));
}

Tensor tensor_f64(const std::vector<double>& v, std::vector<int64_t> shape) {
    Tensor t;
    t.dtype = DType::Float64;
    t.shape = std::move(shape);
    t.bytes.reserve(v.size() * 8);
    for (double x : v) {
        uint64_t bits;
        std::memcpy(&bits, &x, 8);
        put_le(t.bytes, bits, 8);
    }
    return t;
}

Tensor tensor_f32(const std::vector<double>& v, std::vector<int64_t> shape) {
    Tensor t;
    t.dtype = DType::Float32;
    t.shape = std::move(shape);
    t.bytes.reserve(v.size() * 4);
    for (double x : v) {
        float f = static_cast<float>(x);
        uint32_t bits;
        std::memcpy(&bits, &f, 4);
        put_le(t.bytes, bits, 4);
    }
    return t;
}

Tensor tensor_i32(const std::vector<int32_t>& v, std::vector<int64_t> shape) {
    Tensor t;
    t.dtype = DType::Int32;
    t.shape = std::move(shape);
    t.bytes.reserve(v.size() * 4);
    for (int32_t x : v) put_le(t.bytes, static_cast<uint32_t>(x), 4);
    return t;
}

Tensor tensor_i64(const std::vector<int64_t>& v, std::vector<int64_t> shape) {
    Tensor t;
    t.dtype = DType::Int64;
    t.shape = std::move(shape);
    t.bytes.reserve(v.size() * 8);
    for (int64_t x : v) put_le(t.bytes, static_cast<uint64_t>(x), 8);
    return t;
}

namespace detail {
// Portable pi. The POSIX math-constant macros are not standard C++ (MSVC only
// defines them when _USE_MATH_DEFINES is set before <cmath>) and C++17 has no
// std::numbers, so keep our own constant.
constexpr double kPi = 3.14159265358979323846;
}  // namespace detail

double primitive_norm(double alpha, int l) {
    double radial = std::pow(2.0 * alpha / detail::kPi, 0.75);
    double angular = std::pow(4.0 * alpha, l / 2.0);
    double df = 1.0;
    for (int k = 1; k < 2 * l; k += 2) df *= k;  // (2l-1)!!
    return radial * angular / std::sqrt(df);
}

// ---------------------------------------------------------------------------
// Lattice (spec § 4.3, "Lattice conventions")
// ---------------------------------------------------------------------------

static void require_finite(const Lattice::Mat3& m, const char* what) {
    for (const auto& row : m)
        for (double x : row)
            if (!std::isfinite(x))
                throw std::runtime_error(std::string("Lattice ") + what +
                                         " must be finite (no NaN or inf)");
}

Lattice Lattice::from_rows_angstrom(const Mat3& rows) {
    require_finite(rows, "rows_angstrom");
    return Lattice(rows);
}

Lattice Lattice::from_columns_bohr(const Mat3& columns) {
    require_finite(columns, "columns_bohr");
    Mat3 rows{};
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            rows[i][j] = columns[j][i] / kBohrPerAngstrom;   // R = L^T / f
    return Lattice(rows);
}

Lattice::Mat3 Lattice::rows_angstrom() const { return rows_; }

Lattice::Mat3 Lattice::columns_bohr() const {
    Mat3 cols{};
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            cols[i][j] = rows_[j][i] * kBohrPerAngstrom;     // L = f * R^T
    return cols;
}

double Lattice::volume_angstrom3() const {
    const auto& a = rows_[0];
    const auto& b = rows_[1];
    const auto& c = rows_[2];
    const double triple = a[0] * (b[1] * c[2] - b[2] * c[1]) -
                          a[1] * (b[0] * c[2] - b[2] * c[0]) +
                          a[2] * (b[0] * c[1] - b[1] * c[0]);
    return std::fabs(triple);
}

double Lattice::volume_bohr3() const {
    return volume_angstrom3() * kBohrPerAngstrom * kBohrPerAngstrom * kBohrPerAngstrom;
}

// Minimal Z -> symbol table (0..118).
static std::string symbol_for(int z) {
    static const char* S[] = {
        "X","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P",
        "S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh",
        "Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
        "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re",
        "Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
        "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db",
        "Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"};
    if (z >= 0 && z < static_cast<int>(sizeof(S) / sizeof(S[0]))) return S[z];
    return "X";
}

static bool is_valid_id(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s)
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' ||
              c == '.' || c == '-'))
            return false;
    return true;
}

// ---------------------------------------------------------------------------
// Declared units (QVF spec § 4.3, "Declared units"): every declared unit is a
// unit name, and a key whose unit the spec fixes can only restate it. The
// writers declare and never convert.
// ---------------------------------------------------------------------------

struct FixedUnit {
    const char* key;
    Unit unit;
};

static void check_unit_name(const std::string& what, const Json& value) {
    if (!(value.is_string() && parse_unit(value.as_string())))
        throw std::runtime_error(what + ": " + value.dump(0) +
                                 " is not a QVF unit name (spec § 4.3)");
}

static void check_units_object(const std::string& what, const Json& units,
                               std::initializer_list<FixedUnit> fixed) {
    if (!units.is_object())
        throw std::runtime_error(what + ": units must be an object of unit names");
    for (const auto& kv : units.items()) {
        check_unit_name(what + ": declared unit for '" + kv.first + "'", kv.second);
        const Unit u = *parse_unit(kv.second.as_string());
        for (const FixedUnit& f : fixed)
            if (kv.first == f.key && u != f.unit)
                throw std::runtime_error(
                    what + ": units['" + kv.first + "'] must be '" +
                    std::string(unit_name(f.unit)) +
                    "'; spec § 4.3 fixes that unit and a declaration only restates it");
    }
}

// The declared units of a JSON payload: its `units` object and, for a
// one-dimensional spectrum, x_units, y_units and broadening.units. A
// Json::raw fragment is parsed for the check; one that does not parse is left
// alone (it is written verbatim, as before).
static void check_payload_units(const std::string& what, const Json& payload,
                                std::initializer_list<FixedUnit> fixed,
                                bool spectrum_1d = false) {
    Json parsed;
    const Json* doc = &payload;
    if (payload.type() == Json::Type::Raw) {
        try {
            parsed = Json::parse(payload.as_string());
        } catch (const std::exception&) {
            return;
        }
        doc = &parsed;
    }
    if (!doc->is_object()) return;
    if (doc->contains("units")) check_units_object(what, doc->at("units"), fixed);
    if (!spectrum_1d) return;
    for (const char* key : {"x_units", "y_units"})
        if (doc->contains(key)) check_unit_name(what + ": " + key, doc->at(key));
    if (doc->contains("broadening") && doc->at("broadening").is_object() &&
        doc->at("broadening").contains("units"))
        check_unit_name(what + ": broadening units", doc->at("broadening").at("units"));
}

static bool is_spectrum_1d(const std::string& kind) {
    return kind == "spectra.ir" || kind == "spectra.raman" || kind == "spectra.uvvis" ||
           kind == "spectra.ecd" || kind == "spectra.vcd" || kind == "spectra.generic";
}

// Automatic root metadata (QVF spec § 3.5)
// ---------------------------------------------------------------------------

// `value` in decimal, left-padded with zeros to `width` digits (value >= 0).
static std::string metadata_zero_pad(int64_t value, size_t width) {
    std::string s = std::to_string(value);
    return s.size() < width ? std::string(width - s.size(), '0') + s : s;
}

// The current time as "YYYY-MM-DDThh:mm:ssZ" (UTC). The civil date is computed
// from the day count since 1970-01-01 (proleptic Gregorian calendar), which
// avoids the non-portable thread-safe variants of gmtime.
static std::string metadata_utc_now() {
    using namespace std::chrono;
    const int64_t secs =
        duration_cast<seconds>(system_clock::now().time_since_epoch()).count();
    int64_t days = secs / 86400, rem = secs % 86400;
    if (rem < 0) { rem += 86400; --days; }
    days += 719468;                                   // shift the epoch to 0000-03-01
    const int64_t era = (days >= 0 ? days : days - 146096) / 146097;
    const int64_t doe = days - era * 146097;          // [0, 146096]
    const int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const int64_t mp = (5 * doy + 2) / 153;
    const int64_t day = doy - (153 * mp + 2) / 5 + 1;
    const int64_t month = mp < 10 ? mp + 3 : mp - 9;
    const int64_t year = yoe + era * 400 + (month <= 2 ? 1 : 0);
    return metadata_zero_pad(year, 4) + "-" + metadata_zero_pad(month, 2) + "-" +
           metadata_zero_pad(day, 2) + "T" + metadata_zero_pad(rem / 3600, 2) + ":" +
           metadata_zero_pad(rem % 3600 / 60, 2) + ":" + metadata_zero_pad(rem % 60, 2) +
           "Z";
}

// A random (version 4, RFC 4122 variant) UUID in lowercase 8-4-4-4-12 form.
static std::string metadata_random_uuid() {
    std::random_device device;
    const auto tick = static_cast<uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::seed_seq seed{device(), device(), device(), device(),
                       static_cast<unsigned>(tick & 0xFFFFFFFFu),
                       static_cast<unsigned>(tick >> 32)};
    std::mt19937_64 engine(seed);
    uint8_t bytes[16];
    for (int half = 0; half < 2; ++half) {
        const uint64_t word = engine();
        for (int i = 0; i < 8; ++i)
            bytes[half * 8 + i] = static_cast<uint8_t>(word >> (8 * i));
    }
    bytes[6] = static_cast<uint8_t>((bytes[6] & 0x0F) | 0x40);   // version 4
    bytes[8] = static_cast<uint8_t>((bytes[8] & 0x3F) | 0x80);   // variant 10
    static const char* const kHex = "0123456789abcdef";
    std::string out;
    for (int i = 0; i < 16; ++i) {
        if (i == 4 || i == 6 || i == 8 || i == 10) out += '-';
        out += kHex[bytes[i] >> 4];
        out += kHex[bytes[i] & 0x0F];
    }
    return out;
}

// ---------------------------------------------------------------------------
// QvfWriter
// ---------------------------------------------------------------------------

namespace detail {

// The archive file a streaming or updating writer is bound to.
struct WriterOutput {
    // Metadata stays available to legacy setters after payloads stream out.
    std::map<std::string, Json> json_cache;
    std::unique_ptr<container::ZipReader> original; // lazy reads during update
    std::string path;
    bool update = false;                         // the file holds an archive: append
    bool sink_only = false;                      // one-shot caller-owned destination
    std::unique_ptr<container::ZipWriter> zip;   // open while members stream into it
    Json preserved_root = Json::object();        // update: root keys kept verbatim
    Json source;                                 // update: the archive's source object
    bool has_source = false;
    // Entry names the archive held before this update: none may be reused
    // except manifest.json (spec § 2.3).
    std::set<std::string> archived;
    // metadata.modified_utc as the archive has it, and the time stamp this
    // update writes instead (made once, so to_bytes() previews what write()
    // writes).
    Json loaded_modified;
    std::string stamp;
};

}  // namespace detail

struct MemberSink::State {
    QvfWriter* writer = nullptr;
    std::string section_id, role, path;
    DType dtype = DType::Float64;
    std::vector<int64_t> shape;
    uint64_t expected = 0;
    uint64_t written = 0;
    std::unique_ptr<container::EntrySink> entry;   // bound writer: into the archive
    std::vector<uint8_t> buffer;                    // buffered writer
    detail::Sha256 sha;                             // buffered writer
    bool open = true;
};

QvfWriter::QvfWriter(Source source) : source_(std::move(source)) {
    if (source_.program.empty() || source_.version.empty())
        throw std::runtime_error("source.program and source.version are required");
}

void QvfWriter::check_new_id(const std::string& id) const {
    if (!is_valid_id(id))
        throw std::runtime_error("invalid section id: " + id);
    if (id == "auto_open" || id == "bookmarks" || id == "camera" ||
        id == "crossfade")   // viewer_defaults keys (§ 3.1, § 3.4)
        throw std::runtime_error("reserved section id: " + id);
    if (ids_.count(id)) throw std::runtime_error("duplicate section id: " + id);
}

std::string QvfWriter::reserve_id(const std::string& id) {
    check_new_id(id);
    ids_.insert(id);
    return id;
}

// Spec § 2.1: relative, '/'-separated, no '..' segment, no drive letter.
static void check_member_path(const std::string& path) {
    const bool drive = path.size() >= 2 && path[1] == ':' &&
                       std::isalpha(static_cast<unsigned char>(path[0]));
    if (path.empty() || path.front() == '/' || drive ||
        path.find('\\') != std::string::npos || path.find('\0') != std::string::npos ||
        path == ".." || path.find("/../") != std::string::npos ||
        path.rfind("../", 0) == 0 ||
        (path.size() >= 3 && path.compare(path.size() - 3, 3, "/..") == 0))
        throw std::runtime_error("illegal member path: " + path);
}

// A member file written in this session that has not become part of the
// archive yet: buffered, or streamed into a bound archive not yet finished.
bool QvfWriter::pending_file(const std::string& path) const {
    for (const auto& f : files_)
        if (f.first == path) return true;
    return out_ && out_->zip && !out_->archived.count(path) && out_->zip->contains(path);
}

// Spec § 2.3: an update never rewrites an entry the archive holds. A member
// whose default path is taken by the archive being updated (or named by its
// manifest) gets a fresh one, "<stem>.<n><suffix>" with n = 2, 3, ..., the
// stem ending at the first '.' of the last path segment (as the Python
// writer's _add_file names it). A path a member of this session already
// uses is an error, as in a new archive.
std::string QvfWriter::entry_name(const std::string& path) const {
    check_member_path(path);
    const auto taken = [this](const std::string& p) {
        return paths_.count(p) != 0 || (out_ && out_->archived.count(p) != 0);
    };
    if (!taken(path)) return path;
    if (!out_ || pending_file(path))
        throw std::runtime_error("duplicate member path: " + path);
    const size_t slash = path.rfind('/');
    const std::string head = slash == std::string::npos ? "" : path.substr(0, slash + 1);
    const std::string tail = slash == std::string::npos ? path : path.substr(slash + 1);
    const size_t dot = tail.find('.');
    const std::string stem = tail.substr(0, dot);
    const std::string suffix = dot == std::string::npos ? "" : tail.substr(dot);
    for (unsigned long long n = 2;; ++n) {
        const std::string fresh = head + stem + "." + std::to_string(n) + suffix;
        if (!taken(fresh)) return fresh;
    }
}

std::string QvfWriter::add_file(const std::string& path, std::vector<uint8_t> bytes,
                                bool store) {
    if (member_open_)
        throw std::runtime_error("a member sink is open; close it before adding members");
    const std::string name = entry_name(path);
    const container::Method method = store ? container::Method::Store
                                           : container::Method::Deflate;
    if (out_ && out_->zip) {
        // Streaming: straight into the archive. A path the writer dropped and
        // re-adds in this same session supersedes the entry it wrote.
        out_->zip->drop(name);
        out_->zip->write(name, bytes, method, force_zip64);
    } else {
        files_.emplace_back(name, std::move(bytes));
    }
    paths_.insert(name);
    if (store) stored_.insert(name);
    return name;
}

// Forget a member file of this session: a buffered file is discarded, and an
// entry already streamed into the bound archive is dropped from its directory
// (its bytes stay as unreferenced space). An entry the archive held before
// this session, or a path its manifest names, is never touched (spec § 2.3):
// it stays taken.
void QvfWriter::drop_file(const std::string& path) {
    if (!pending_file(path)) return;
    paths_.erase(path);
    stored_.erase(path);
    files_.erase(std::remove_if(files_.begin(), files_.end(),
                                [&path](const std::pair<std::string, std::vector<uint8_t>>& f) {
                                    return f.first == path;
                                }),
                 files_.end());
    if (out_) out_->json_cache.erase(path);
    if (out_ && out_->zip) out_->zip->drop(path);
}

Json QvfWriter::json_member(const std::string& path, const Json& obj) {
    std::string dumped = obj.dump(2);
    std::vector<uint8_t> bytes(dumped.begin(), dumped.end());
    std::string sha = detail::Sha256::hex_of(bytes);
    const long long size = static_cast<long long>(bytes.size());
    const std::string name = add_file(path, std::move(bytes));
    if (out_) out_->json_cache[name] = obj;
    Json m = Json::object();
    m.set("path", name).set("format", "json").set("sha256", sha).set("size", size);
    return m;
}

Json QvfWriter::binary_member(const std::string& path, const Tensor& t,
                              bool with_shape) {
    // Sanity: byte length must equal itemsize * prod(shape).
    size_t n = 1;
    for (int64_t d : t.shape) n *= static_cast<size_t>(d);
    if (with_shape && t.bytes.size() != dtype_itemsize(t.dtype) * n)
        throw std::runtime_error("binary member byte length != itemsize*prod(shape): "
                                 + path);
    std::string sha = detail::Sha256::hex_of(t.bytes);
    std::vector<uint8_t> bytes = t.bytes;
    const std::string name = add_file(path, std::move(bytes));
    Json m = Json::object();
    m.set("path", name).set("format", "binary");
    if (with_shape) {
        m.set("dtype", std::string(dtype_name(t.dtype)));
        Json shape = Json::array();
        for (int64_t d : t.shape) shape.push_back(Json(static_cast<long long>(d)));
        m.set("shape", shape);
    }
    m.set("sha256", sha);
    m.set("size", static_cast<long long>(t.bytes.size()));
    return m;
}

// "type/subtype" without parameters (RFC 6838): the pattern of the schema's
// $defs/MediaType, ^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]*/[A-Za-z0-9][...]*$.
static bool is_media_type(const std::string& s) {
    const auto token = [&s](size_t from, size_t to) {
        if (from >= to || !std::isalnum(static_cast<unsigned char>(s[from]))) return false;
        for (size_t i = from; i < to; ++i) {
            const unsigned char c = static_cast<unsigned char>(s[i]);
            if (!(std::isalnum(c) || (c != 0 && std::strchr("!#$&^_.+-", c)))) return false;
        }
        return true;
    };
    const size_t slash = s.find('/');
    return slash != std::string::npos && token(0, slash) && token(slash + 1, s.size());
}

// An opaque binary member (no dtype or shape, spec § 4.2) of `bytes`, with
// its media_type when not empty (checked before anything is written).
Json QvfWriter::opaque_member(const std::string& path, std::vector<uint8_t> bytes,
                              const std::string& media_type, bool store) {
    if (!media_type.empty() && !is_media_type(media_type))
        throw std::runtime_error("member " + path + ": media type '" + media_type +
                                 "' is not type/subtype (QVF spec § 4.2)");
    const std::string sha = detail::Sha256::hex_of(bytes);
    const long long size = static_cast<long long>(bytes.size());
    const std::string name = add_file(path, std::move(bytes), store);
    Json m = Json::object();
    m.set("path", name).set("format", "binary").set("sha256", sha).set("size", size);
    if (!media_type.empty()) m.set("media_type", media_type);
    return m;
}

// -- typed provenance and job.spec fields (spec § 3.2, § 5.9) ------------------

namespace {

[[noreturn]] void field_error(const std::string& what, const std::string& rule) {
    throw std::runtime_error(what + " " + rule);
}

// A pre-serialized fragment nested in a built object (the C ABI path) is
// parsed so that it can be checked; anything else is returned as it is.
Json resolved(const Json& v, const std::string& what) {
    if (v.type() != Json::Type::Raw) return v;
    try {
        return Json::parse(v.as_string());
    } catch (const std::exception& e) {
        throw std::runtime_error(what + " is not valid JSON: " + e.what());
    }
}

bool is_digit(char c) { return c >= '0' && c <= '9'; }

// 7 to 64 hexadecimal characters (provenance.git_commit).
bool is_git_commit(const std::string& s) {
    if (s.size() < 7 || s.size() > 64) return false;
    for (char c : s)
        if (!(is_digit(c) || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F')))
            return false;
    return true;
}

// ISO 8601 UTC, YYYY-MM-DDThh:mm:ssZ with an optional fractional second: the
// pattern of $defs/Metadata created_utc (spec § 3.5).
bool is_utc_timestamp(const std::string& s) {
    static const size_t kDigits[] = {0, 1, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18};
    if (s.size() < 20) return false;
    for (size_t i : kDigits)
        if (!is_digit(s[i])) return false;
    if (s[4] != '-' || s[7] != '-' || s[10] != 'T' || s[13] != ':' || s[16] != ':')
        return false;
    auto two = [&s](size_t i) { return (s[i] - '0') * 10 + (s[i + 1] - '0'); };
    if (two(5) < 1 || two(5) > 12 || two(8) < 1 || two(8) > 31) return false;
    if (two(11) > 23 || two(14) > 59 || two(17) > 60) return false;
    size_t pos = 19;
    if (s[pos] == '.') {
        const size_t start = ++pos;
        while (pos < s.size() && is_digit(s[pos])) ++pos;
        if (pos == start) return false;
    }
    return pos + 1 == s.size() && s[pos] == 'Z';
}

bool is_vendor_key(const std::string& key) { return key.rfind("x_", 0) == 0; }

const std::string& checked_text(const Json& v, const std::string& what) {
    if (!v.is_string() || v.as_string().empty())
        field_error(what, "must be a non-empty string");
    return v.as_string();
}

int64_t checked_int(const Json& v, const std::string& what, int64_t minimum) {
    if (v.type() != Json::Type::Int || v.as_int() < minimum)
        field_error(what, "must be an integer >= " + std::to_string(minimum));
    return v.as_int();
}

double checked_number(const Json& v, const std::string& what, bool positive) {
    if (!v.is_number() || !std::isfinite(v.as_double()))
        field_error(what, "must be a finite number");
    if (positive && !(v.as_double() > 0.0)) field_error(what, "must be a number > 0");
    return v.as_double();
}

void check_software(const Json& list) {
    if (!list.is_array())
        throw std::runtime_error(
            "provenance software must be an array of {name, version?, url?} "
            "objects (QVF spec § 3.2)");
    size_t index = 0;
    for (const Json& raw : list.elements()) {
        const std::string what = "provenance software[" + std::to_string(index++) + "]";
        const Json entry = resolved(raw, what);
        if (!entry.is_object())
            field_error(what, "must be an object with a 'name'");
        for (const auto& kv : entry.items())
            if (kv.first != "name" && kv.first != "version" && kv.first != "url" &&
                !is_vendor_key(kv.first))
                throw std::runtime_error(
                    what + ": key '" + kv.first +
                    "' is not one of name, url, version; vendor keys must begin "
                    "with 'x_' (QVF spec § 3.2)");
        if (!entry.has("name")) field_error(what, "needs a 'name'");
        checked_text(entry.at("name"), what + " name");
        for (const char* key : {"version", "url"})
            if (entry.has(key) && !entry.at(key).is_string())
                field_error(what + " " + key, "must be a string");
    }
}

// The typed fields of the root provenance block (spec § 3.2). The block is
// open: every other key passes unchecked.
void check_provenance_fields(const Json& block) {
    for (const auto& kv : block.items()) {
        const std::string& key = kv.first;
        const std::string what = "provenance " + key;
        if (key == "software") {
            check_software(resolved(kv.second, what));
        } else if (key == "git_commit") {
            const Json v = resolved(kv.second, what);
            if (!v.is_string() || !is_git_commit(v.as_string()))
                field_error(what, "must be a git commit hash of 7 to 64 "
                                  "hexadecimal characters (QVF spec § 3.2)");
        } else if (key == "n_threads" || key == "n_mpi_ranks") {
            checked_int(resolved(kv.second, what), what, 1);
        } else if (key == "n_gpus") {
            checked_int(resolved(kv.second, what), what, 0);
        } else if (key == "started_utc" || key == "finished_utc") {
            const Json v = resolved(kv.second, what);
            if (!v.is_string() || !is_utc_timestamp(v.as_string()))
                field_error(what, "must be an ISO 8601 UTC time stamp "
                                  "'YYYY-MM-DDThh:mm:ssZ' (QVF spec § 3.2)");
        }
    }
}

enum class Rule { Text, Number, PositiveNumber, Count, Coordinates };

struct FieldRule {
    const char* name;
    Rule rule;
};

// One typed block of a job.spec payload (spec § 5.9): every key is one of
// `rules` or begins with x_, and `required` is present.
void check_block(const std::string& block, const Json& value,
                 std::initializer_list<FieldRule> rules, const char* required) {
    const std::string whole = "job.spec " + block;
    const Json v = resolved(value, whole);
    if (!v.is_object())
        field_error(whole, "must be an object (QVF spec § 5.9)");
    std::string names;
    for (const FieldRule& r : rules) names += (names.empty() ? "" : ", ") + std::string(r.name);
    for (const auto& kv : v.items()) {
        if (is_vendor_key(kv.first)) continue;
        const FieldRule* found = nullptr;
        for (const FieldRule& r : rules)
            if (kv.first == r.name) found = &r;
        if (!found)
            throw std::runtime_error(
                whole + ": key '" + kv.first + "' is not one of " + names +
                "; vendor keys must begin with 'x_' and engine keywords go under "
                "options (QVF spec § 5.9)");
        const std::string what = whole + "." + kv.first;
        switch (found->rule) {
            case Rule::Text: checked_text(kv.second, what); break;
            case Rule::Number: checked_number(kv.second, what, false); break;
            case Rule::PositiveNumber: checked_number(kv.second, what, true); break;
            case Rule::Count: checked_int(kv.second, what, 1); break;
            case Rule::Coordinates: {
                const std::string& c = checked_text(kv.second, what);
                if (c != "cartesian" && c != "internal" && c != "redundant")
                    field_error(what, "must be one of cartesian, internal, redundant");
                break;
            }
        }
    }
    if (required && !v.has(required))
        field_error(whole, std::string("requires '") + required + "' (QVF spec § 5.9)");
}

}  // namespace

// -- root metadata ----------------------------------------------------------

void QvfWriter::set_provenance(Json v) { set_provenance(std::move(v), false); }
void QvfWriter::set_provenance(Json v, bool allow_hostname) {
    // A pre-serialized fragment (the C ABI path) is parsed to check the typed
    // fields and to look for `hostname`; without a key to drop, the caller's
    // text stays verbatim.
    if (v.type() == Json::Type::Raw) {
        Json parsed;
        try {
            parsed = Json::parse(v.as_string());
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("provenance is not valid JSON: ") +
                                     e.what());
        }
        if (!parsed.is_object())
            throw std::runtime_error(
                "provenance must be a JSON object (QVF spec § 3.2)");
        check_provenance_fields(parsed);
        // `hostname` is deprecated (spec § 3.2) and dropped unless allowed.
        if (!allow_hostname && parsed.erase("hostname")) v = std::move(parsed);
    } else {
        if (!v.is_object())
            throw std::runtime_error(
                "provenance must be a JSON object (QVF spec § 3.2)");
        check_provenance_fields(v);
        if (!allow_hostname) v.erase("hostname");
    }
    root_.set("provenance", std::move(v));
}
void QvfWriter::set_thermochemistry(Json v) { root_.set("thermochemistry", std::move(v)); }
void QvfWriter::set_dipole_moment(Json v) { root_.set("dipole_moment", std::move(v)); }
void QvfWriter::set_constraints(Json v) { root_.set("constraints", std::move(v)); }
void QvfWriter::set_energies(Json v) {
    // Root energies block (spec § 3.6). A Json::raw fragment (the C ABI path)
    // is emitted verbatim; a built object must carry `total`.
    if (v.type() == Json::Type::Object) {
        if (!v.has("total"))
            throw std::runtime_error("energies requires a total {value, units}");
    } else if (v.type() != Json::Type::Raw) {
        throw std::runtime_error("energies must be a JSON object");
    }
    root_.set("energies", std::move(v));
}
void QvfWriter::set_extensions(Json v) { root_.set("extensions", std::move(v)); }
void QvfWriter::set_viewer_defaults(Json v) { root_.set("viewer_defaults", std::move(v)); }

void QvfWriter::set_metadata(Json v) {
    if (v.type() == Json::Type::Raw) {
        // A pre-serialized fragment (the C ABI path) is parsed, so that the
        // automatic created_utc / uuid can be merged into it.
        try {
            v = Json::parse(v.as_string());
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("metadata is not valid JSON: ") +
                                     e.what());
        }
    }
    if (!v.is_object())
        throw std::runtime_error("metadata must be a JSON object (QVF spec § 3.5)");
    root_.set("metadata", std::move(v));
}

// -- section references and grouping (spec § 4.1) ---------------------------

size_t QvfWriter::section_index(const std::string& id) const {
    const std::vector<Json>& sections = sections_.elements();
    for (size_t i = 0; i < sections.size(); ++i)
        if (sections[i].at("id").as_string() == id) return i;
    throw std::runtime_error("no section with id '" + id + "'");
}

void QvfWriter::check_section_ref(const std::string& key, const Json& value) const {
    const char* target = key == "structure_ref" ? "structure"
                         : key == "run_ref"     ? "run.record"
                                                : nullptr;
    if (!target) return;
    if (!value.is_string())
        throw std::runtime_error(key + " must be a string naming a '" + target +
                                 "' section");
    const std::string& ref = value.as_string();
    const std::vector<Json>& sections = sections_.elements();
    for (const Json& s : sections)
        if (s.at("id").as_string() == ref) {
            if (s.at("kind").as_string() != target)
                throw std::runtime_error(key + " '" + ref + "' names a '" +
                                         s.at("kind").as_string() +
                                         "' section, not a '" + target +
                                         "' section (QVF spec § 4.1)");
            return;
        }
    throw std::runtime_error(key + " '" + ref + "' names no section");
}

void QvfWriter::set_section_field(const std::string& id, const std::string& key,
                                  Json value) {
    if (key.empty() || key == "id" || key == "kind" || key == "members")
        throw std::runtime_error("set_section_field: field '" + key +
                                 "' cannot be set");
    const size_t index = section_index(id);
    check_section_ref(key, value);
    // Json has no mutable element access: rebuild the array with the section
    // replaced.
    const std::vector<Json>& sections = sections_.elements();
    Json updated = sections[index];
    updated.set(key, std::move(value));
    Json rebuilt = Json::array();
    for (size_t i = 0; i < sections.size(); ++i)
        rebuilt.push_back(i == index ? updated : sections[i]);
    sections_ = std::move(rebuilt);
}

std::string QvfWriter::annotate_section(const std::string& id,
                                        const std::string& structure_ref,
                                        const std::string& run_ref,
                                        const std::string& group) {
    section_index(id);   // an unknown id throws before anything is checked
    // Check both references first, so a failure changes nothing.
    if (!structure_ref.empty()) check_section_ref("structure_ref", structure_ref);
    if (!run_ref.empty()) check_section_ref("run_ref", run_ref);
    if (!structure_ref.empty()) set_section_field(id, "structure_ref", structure_ref);
    if (!run_ref.empty()) set_section_field(id, "run_ref", run_ref);
    if (!group.empty()) set_section_field(id, "group", group);
    return id;
}

void QvfWriter::require_section_kind(const std::string& what, const std::string& ref,
                                     const char* kind) const {
    for (const Json& s : sections_.elements())
        if (s.at("id").as_string() == ref) {
            if (s.at("kind").as_string() != kind)
                throw std::runtime_error(what + " '" + ref + "' must name an existing '" +
                                         kind + "' section, found a '" +
                                         s.at("kind").as_string() + "' section");
            return;
        }
    throw std::runtime_error(what + " '" + ref + "' must name an existing '" + kind +
                             "' section, found no such section");
}

// -- removing and replacing sections (spec § 2.3) ---------------------------

Json QvfWriter::remove_section(const std::string& id) {
    if (member_open_)
        throw std::runtime_error("a member sink is open; close it before removing a section");
    const size_t index = section_index(id);
    const std::vector<Json>& sections = sections_.elements();
    Json removed = sections[index];
    Json rebuilt = Json::array();
    for (size_t i = 0; i < sections.size(); ++i)
        if (i != index) rebuilt.push_back(sections[i]);
    sections_ = std::move(rebuilt);
    ids_.erase(id);
    structure_sum_z_.erase(id);
    // Members not yet part of the archive are dropped; drop_file() leaves an
    // entry the archive already holds, which becomes a superseded entry.
    if (removed.has("members") && removed.at("members").is_object())
        for (const auto& kv : removed.at("members").items())
            if (kv.second.is_object() && kv.second.has("path") &&
                kv.second.at("path").is_string())
                drop_file(kv.second.at("path").as_string());
    return removed;
}

// -- structure --------------------------------------------------------------

std::string QvfWriter::add_structure(const std::vector<Atom>& atoms,
                                     const std::array<bool, 3>& pbc,
                                     const std::vector<std::array<double, 3>>& lattice,
                                     const std::vector<Bond>& bonds,
                                     const std::string& section_id) {
    check_new_id(section_id);
    const bool periodic = pbc[0] || pbc[1] || pbc[2];
    if (periodic && lattice.empty()) {
        throw std::runtime_error(
            "add_structure: pbc marks a periodic axis, so lattice_vectors is "
            "required (QVF spec § 5.1)");
    }
    if (!lattice.empty() && lattice.size() != 3) {
        throw std::runtime_error(
            "add_structure: lattice_vectors must be 3 row vectors");
    }
    Json atoms_json = Json::array();
    for (const auto& a : atoms) {
        Json aj = Json::object();
        std::string sym = a.symbol.empty() ? symbol_for(a.atomic_number) : a.symbol;
        Json pos = Json::array();
        for (double x : a.position) pos.push_back(Json(x));
        aj.set("symbol", sym).set("position", pos)
          .set("atomic_number", a.atomic_number);
        atoms_json.push_back(aj);
    }
    Json pbc_json = Json::array();
    for (bool p : pbc) pbc_json.push_back(Json(p));
    Json structure = Json::object();
    structure.set("atoms", atoms_json).set("pbc", pbc_json);
    if (!lattice.empty()) {
        Json lat = Json::array();
        for (const auto& row : lattice) {
            Json r = Json::array();
            for (double x : row) r.push_back(Json(x));
            lat.push_back(r);
        }
        structure.set("lattice_vectors", lat);
    } else {
        structure.set("lattice_vectors", Json(nullptr));
    }

    Json members = Json::object();
    members.set("structure", json_member(section_id + "/structure.json", structure));
    if (!bonds.empty()) {
        Json pairs = Json::array();
        for (const auto& b : bonds) {
            Json p = Json::object();
            p.set("i", b.i).set("j", b.j).set("order", b.order);
            pairs.push_back(p);
        }
        Json bj = Json::object();
        bj.set("pairs", pairs);
        members.set("bonds", json_member(section_id + "/bonds.json", bj));
    }
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "structure")
           .set("members", members);
    sections_.push_back(section);
    // Σ Z for set_structure_identity's electron count (spec § 5.1).
    int64_t sum_z = 0;
    bool known = true;
    for (const auto& a : atoms) {
        known = known && a.atomic_number >= 0;
        sum_z += a.atomic_number;
    }
    if (known) structure_sum_z_[section_id] = sum_z;
    return section_id;
}

std::string QvfWriter::add_structure(const std::vector<Atom>& atoms,
                                     const std::array<bool, 3>& pbc,
                                     const Lattice& lattice,
                                     const std::vector<Bond>& bonds,
                                     const std::string& section_id) {
    const Lattice::Mat3 rows = lattice.rows_angstrom();
    return add_structure(
        atoms, pbc,
        std::vector<std::array<double, 3>>(rows.begin(), rows.end()), bonds,
        section_id);
}

// -- chemical identity (spec § 5.1) ------------------------------------------

namespace {

// ECMAScript \s (the JSON Schema pattern dialect): a SMILES or InChI holds
// none of these code points. Bytes that are not UTF-8 count as non-space.
bool has_pattern_space(const std::string& s) {
    for (size_t i = 0; i < s.size();) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        uint32_t cp = c;
        size_t len = 1;
        if ((c & 0xE0) == 0xC0) { cp = c & 0x1Fu; len = 2; }
        else if ((c & 0xF0) == 0xE0) { cp = c & 0x0Fu; len = 3; }
        else if ((c & 0xF8) == 0xF0) { cp = c & 0x07u; len = 4; }
        else if (c >= 0x80) { ++i; continue; }   // not a lead byte
        if (len > s.size() - i) return false;
        for (size_t k = 1; k < len; ++k)
            cp = (cp << 6) | (static_cast<unsigned char>(s[i + k]) & 0x3Fu);
        i += len;
        if ((cp >= 0x09 && cp <= 0x0D) || cp == 0x20 || cp == 0xA0 || cp == 0x1680 ||
            (cp >= 0x2000 && cp <= 0x200A) || cp == 0x2028 || cp == 0x2029 ||
            cp == 0x202F || cp == 0x205F || cp == 0x3000 || cp == 0xFEFF)
            return true;
    }
    return false;
}

// An InChIKey: 14 uppercase letters, '-', 10 uppercase letters, '-', one.
bool is_inchikey(const std::string& s) {
    if (s.size() != 27 || s[14] != '-' || s[25] != '-') return false;
    for (size_t i = 0; i < s.size(); ++i)
        if (i != 14 && i != 25 && !(s[i] >= 'A' && s[i] <= 'Z')) return false;
    return true;
}

// Why `charge` and `multiplicity` (has_mult) are impossible for atoms whose
// atomic numbers sum to `sum_z` >= 0 (spec § 5.1), or "": the electron count
// N = sum_z - charge is not negative, and multiplicity - 1 neither exceeds N
// nor differs from it in parity. Exact for every int64 value.
std::string charge_multiplicity_problem(int64_t sum_z, int64_t charge, bool has_mult,
                                        int64_t mult) {
    if (charge > sum_z)
        return "charge " + std::to_string(charge) + " leaves " +
               std::to_string(sum_z - charge) +
               " electrons for atoms with atomic numbers summing to " + std::to_string(sum_z);
    if (!has_mult) return "";
    const uint64_t electrons = static_cast<uint64_t>(sum_z) - static_cast<uint64_t>(charge);
    const uint64_t unpaired = static_cast<uint64_t>(mult - 1);
    if (unpaired > electrons || ((electrons ^ unpaired) & 1u))
        return "multiplicity " + std::to_string(mult) + " is impossible with " +
               std::to_string(electrons) + " electrons (charge " + std::to_string(charge) +
               "): the unpaired electrons must not exceed the electron count and must "
               "have its parity";
    return "";
}

// Σ Z of a structure payload, when every atom has a non-negative integer
// atomic_number (the validator's rule, spec § 5.1).
bool payload_sum_z(const Json& payload, int64_t& out) {
    if (!payload.is_object() || !payload.has("atoms") || !payload.at("atoms").is_array())
        return false;
    int64_t sum = 0;
    for (const Json& a : payload.at("atoms").elements()) {
        if (!a.is_object() || !a.has("atomic_number") ||
            a.at("atomic_number").type() != Json::Type::Int)
            return false;
        const int64_t z = a.at("atomic_number").as_int();
        if (z < 0 || sum > INT64_MAX - z) return false;
        sum += z;
    }
    out = sum;
    return true;
}

}  // namespace

std::string QvfWriter::set_structure_identity(const std::string& id, const Json& identity_in) {
    const size_t index = section_index(id);
    Json section = sections_.elements()[index];
    if (!section.has("kind") || !section.at("kind").is_string() ||
        section.at("kind").as_string() != "structure")
        throw std::runtime_error("set_structure_identity: section '" + id +
                                 "' is not a structure section");
    const Json identity = resolved(identity_in, "structure identity");
    if (!identity.is_object())
        throw std::runtime_error("structure identity must be a JSON object (QVF spec § 5.1)");
    static const char* const kKeys[] = {"charge", "multiplicity", "name",
                                        "smiles", "inchi",        "inchikey"};
    for (const auto& kv : identity.items())
        if (std::find_if(std::begin(kKeys), std::end(kKeys),
                         [&kv](const char* k) { return kv.first == k; }) == std::end(kKeys))
            throw std::runtime_error("structure identity: key '" + kv.first +
                                     "' is not one of charge, multiplicity, name, smiles, "
                                     "inchi, inchikey (QVF spec § 5.1)");
    for (const char* key : kKeys) {
        if (!identity.has(key)) continue;
        const Json& v = identity.at(key);
        const std::string what = std::string("structure ") + key;
        if (v.is_null()) {
            section.erase(key);
            continue;
        }
        if (std::strcmp(key, "charge") == 0) {
            if (v.type() != Json::Type::Int) field_error(what, "must be an integer");
        } else if (std::strcmp(key, "multiplicity") == 0) {
            checked_int(v, what, 1);
        } else {
            const std::string& text = checked_text(v, what);
            const bool ok = std::strcmp(key, "name") == 0 ||
                            (std::strcmp(key, "smiles") == 0 && !has_pattern_space(text)) ||
                            (std::strcmp(key, "inchi") == 0 && text.size() > 6 &&
                             text.compare(0, 6, "InChI=") == 0 && !has_pattern_space(text)) ||
                            (std::strcmp(key, "inchikey") == 0 && is_inchikey(text));
            if (!ok)
                throw std::runtime_error(what + " '" + text +
                                         "' is not of the form QVF spec § 5.1 gives");
        }
        section.set(key, v);
    }
    // The charge and multiplicity the section will have, against its atoms.
    const auto sum_z = structure_sum_z_.find(id);
    if (sum_z != structure_sum_z_.end() && section.has("charge") &&
        section.at("charge").type() == Json::Type::Int) {
        const bool has_mult = section.has("multiplicity") &&
                              section.at("multiplicity").type() == Json::Type::Int &&
                              section.at("multiplicity").as_int() >= 1;
        const std::string problem = charge_multiplicity_problem(
            sum_z->second, section.at("charge").as_int(), has_mult,
            has_mult ? section.at("multiplicity").as_int() : 0);
        if (!problem.empty())
            throw std::runtime_error("structure '" + id + "': " + problem +
                                     " (QVF spec § 5.1)");
    }
    replace_section(index, std::move(section));
    return id;
}

std::string QvfWriter::add_bonds(const std::vector<Bond>& pairs,
                                 const std::string& section_id,
                                 const std::string& structure_ref) {
    // Check the reference and the id first, so a failure adds nothing.
    if (!structure_ref.empty()) check_section_ref("structure_ref", structure_ref);
    check_new_id(section_id);
    Json pj = Json::array();
    for (const auto& b : pairs) {
        Json p = Json::object();
        p.set("i", b.i).set("j", b.j).set("order", b.order);
        pj.push_back(p);
    }
    Json payload = Json::object();
    payload.set("pairs", pj);
    Json members = Json::object();
    members.set("bonds", json_member(section_id + "/connectivity.json", payload));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "bonds");
    if (!structure_ref.empty()) section.set("structure_ref", structure_ref);
    section.set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_bond_orders(const std::string& method, const Json& pairs,
                                       const std::string& section_id) {
    check_new_id(section_id);
    Json payload = Json::object();
    payload.set("method", method).set("pairs", pairs);
    Json members = Json::object();
    members.set("bond_orders", json_member(section_id + "/orders.json", payload));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "bond_orders")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

// -- wavefunction -----------------------------------------------------------

static Tensor matrix_to_tensor(const Matrix& m) {
    if (m.data.size() != static_cast<size_t>(m.rows * m.cols))
        throw std::runtime_error("Matrix data size does not match rows*cols");
    return tensor_f64(m.data, {m.rows, m.cols});
}

static bool matrix_empty(const Matrix& m) {
    return m.rows == 0 && m.cols == 0 && m.data.empty();
}

// One MO coefficient matrix (QVF spec § 5.3): [rows, cols] when real, or
// [rows, cols, 2] from (re, im)-interleaved data when complex; cols must be
// the AO extent (n_ao, or 2 n_ao for a generalized wavefunction).
static Tensor coefficient_tensor(const Matrix& m, bool complex_coefficients,
                                 int64_t ao_extent, const std::string& what) {
    if (m.rows < 0 || m.cols < 0)
        throw std::runtime_error(what + ": rows and cols must be >= 0");
    const size_t want = static_cast<size_t>(m.rows) * static_cast<size_t>(m.cols) *
                        (complex_coefficients ? 2u : 1u);
    if (m.data.size() != want)
        throw std::runtime_error(what + ": data has " + std::to_string(m.data.size()) +
                                 " doubles but rows * cols" +
                                 (complex_coefficients ? " * 2" : "") + " is " +
                                 std::to_string(want));
    if (m.cols != ao_extent)
        throw std::runtime_error(what + " has " + std::to_string(m.cols) +
                                 " AO columns but the basis gives " +
                                 std::to_string(ao_extent));
    if (complex_coefficients) return tensor_f64(m.data, {m.rows, m.cols, 2});
    return tensor_f64(m.data, {m.rows, m.cols});
}

// An overlap or density matrix: real [n_ao, n_ao].
static Tensor ao_matrix_tensor(const Matrix& m, int64_t n_ao, const std::string& what) {
    if (m.rows != n_ao || m.cols != n_ao)
        throw std::runtime_error(what + " must be [n_ao, n_ao] = [" + std::to_string(n_ao) +
                                 ", " + std::to_string(n_ao) + "], got [" +
                                 std::to_string(m.rows) + ", " + std::to_string(m.cols) +
                                 "]");
    return matrix_to_tensor(m);
}

static bool is_finite_number(const Json& j) {
    return (j.type() == Json::Type::Int || j.type() == Json::Type::Double) &&
           std::isfinite(j.as_double());
}

// The `ecp` entries of a basis payload (QVF spec § 5.3), checked and parsed:
// at most one per atom; parallel projector arrays of finite numbers, exponents
// above 0 and integer powers; at most one local part per centre, with the
// centre's highest l. Keys outside the definitions must begin with x_.
static Json checked_ecp(const Json& ecp_in) {
    Json ecp = ecp_in;
    if (ecp.type() == Json::Type::Raw) {
        try {
            ecp = Json::parse(ecp.as_string());
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("ecp is not valid JSON: ") + e.what());
        }
    }
    if (!ecp.is_array())
        throw std::runtime_error(
            "ecp must be a JSON array of {center, n_core_electrons, label?, projectors?}");
    auto is_x = [](const std::string& key) { return key.compare(0, 2, "x_") == 0; };
    auto non_negative_int = [](const Json& o, const char* key) {
        return o.has(key) && o.at(key).type() == Json::Type::Int && o.at(key).as_int() >= 0;
    };
    std::set<int64_t> centres;
    size_t k = 0;
    for (const Json& entry : ecp.elements()) {
        const std::string what = "ecp[" + std::to_string(k++) + "]";
        if (!entry.is_object()) throw std::runtime_error(what + " must be an object");
        for (const auto& kv : entry.items())
            if (kv.first != "center" && kv.first != "n_core_electrons" &&
                kv.first != "label" && kv.first != "projectors" && !is_x(kv.first))
                throw std::runtime_error(what + ": unknown key '" + kv.first + "'");
        if (!non_negative_int(entry, "center"))
            throw std::runtime_error(what + ": center must be an integer >= 0");
        if (!non_negative_int(entry, "n_core_electrons"))
            throw std::runtime_error(what + ": n_core_electrons must be an integer >= 0");
        if (!centres.insert(entry.at("center").as_int()).second)
            throw std::runtime_error(what + ": atom " +
                                     std::to_string(entry.at("center").as_int()) +
                                     " already has an ecp entry; there is at most one per atom");
        if (entry.has("label") && !entry.at("label").is_string())
            throw std::runtime_error(what + ": label must be a string");
        if (!entry.has("projectors")) continue;
        const Json& projectors = entry.at("projectors");
        if (!projectors.is_array())
            throw std::runtime_error(what + ": projectors must be an array");
        int n_local = 0;
        int64_t local_l = -1, max_semi_l = -1;
        size_t j = 0;
        for (const Json& p : projectors.elements()) {
            const std::string pw = what + ".projectors[" + std::to_string(j++) + "]";
            if (!p.is_object()) throw std::runtime_error(pw + " must be an object");
            for (const auto& kv : p.items())
                if (kv.first != "l" && kv.first != "exponents" && kv.first != "coefficients" &&
                    kv.first != "powers" && kv.first != "local" && !is_x(kv.first))
                    throw std::runtime_error(pw + ": unknown key '" + kv.first + "'");
            if (!non_negative_int(p, "l"))
                throw std::runtime_error(pw + ": l must be an integer >= 0");
            size_t n = 0;
            for (const char* key : {"exponents", "coefficients", "powers"}) {
                if (!p.has(key) || !p.at(key).is_array() || p.at(key).size() == 0)
                    throw std::runtime_error(pw + ": " + key + " must be a non-empty array");
                for (const Json& x : p.at(key).elements()) {
                    const bool ok = std::string(key) == "powers"
                                        ? x.type() == Json::Type::Int
                                        : is_finite_number(x) &&
                                              (std::string(key) != "exponents" ||
                                               x.as_double() > 0);
                    if (!ok)
                        throw std::runtime_error(
                            pw + ": " + key +
                            (std::string(key) == "powers"
                                 ? " must hold integers"
                                 : std::string(key) == "exponents"
                                       ? " must hold finite numbers above 0"
                                       : " must hold finite numbers"));
                }
                if (n && p.at(key).size() != n)
                    throw std::runtime_error(
                        pw + ": exponents, coefficients and powers must have the same length");
                n = p.at(key).size();
            }
            bool local = false;
            if (p.has("local")) {
                if (!p.at("local").is_bool())
                    throw std::runtime_error(pw + ": local must be a boolean");
                local = p.at("local").as_bool();
            }
            const int64_t l = p.at("l").as_int();
            if (local) {
                ++n_local;
                local_l = l;
            } else if (l > max_semi_l) {
                max_semi_l = l;
            }
        }
        if (n_local > 1) throw std::runtime_error(what + ": more than one local projector");
        if (n_local == 1 && max_semi_l >= local_l)
            throw std::runtime_error(what + ": the local projector's l (" +
                                     std::to_string(local_l) +
                                     ") must exceed the l of every semi-local channel");
    }
    return ecp;
}

// The `shells` array of a basis payload (spec § 5.3, Appendix A) and its AO
// count; with `libint_normalized` every contraction coefficient is divided by
// the primitive norm N_i. Shared by add_wavefunction_gto and
// add_wavefunction_bloch.
static Json basis_shells_json(const std::vector<Shell>& shells, bool libint_normalized,
                              int64_t& n_ao) {
    Json shells_json = Json::array();
    n_ao = 0;
    for (const auto& sh : shells) {
        if (sh.exponents.size() != sh.coefficients.size())
            throw std::runtime_error("shell exponents/coefficients length mismatch");
        Json exps = Json::array(), coeffs = Json::array();
        for (size_t i = 0; i < sh.exponents.size(); ++i) {
            double c = sh.coefficients[i];
            if (libint_normalized)
                c /= primitive_norm(sh.exponents[i], sh.l);
            exps.push_back(Json(sh.exponents[i]));
            coeffs.push_back(Json(c));
        }
        Json sj = Json::object();
        sj.set("center", sh.center).set("l", sh.l)
          .set("exponents", exps).set("coefficients", coeffs)
          .set("pure", sh.pure);
        shells_json.push_back(sj);
        n_ao += sh.pure ? (2 * sh.l + 1) : ((sh.l + 1) * (sh.l + 2) / 2);
    }
    return shells_json;
}

std::string QvfWriter::add_wavefunction_gto(const WavefunctionGTO& wf) {
    check_new_id(wf.section_id);
    int64_t n_ao = 0;
    Json shells_json = basis_shells_json(wf.shells, wf.coeffs_are_libint_normalized, n_ao);
    // Natural transition orbitals (spec § 5.3): checked before anything is
    // written.
    if (!wf.nto_role.empty() &&
        (wf.orbital_kind != "nto" || (wf.nto_role != "hole" && wf.nto_role != "particle")))
        throw std::runtime_error(
            "wavefunction.gto nto_role is 'hole' or 'particle', for orbital_kind 'nto' "
            "(QVF spec § 5.3)");
    if (wf.state_index < 0)
        throw std::runtime_error("wavefunction.gto state_index must be >= 1 (0 omits it)");
    if (wf.state_index > 0 && wf.state_ref.empty())
        throw std::runtime_error("wavefunction.gto state_index needs state_ref (QVF spec § 5.3)");
    auto add_state_fields = [&wf](Json& meta) {
        if (!wf.state_ref.empty()) meta.set("state_ref", wf.state_ref);
        if (wf.state_index > 0) meta.set("state_index", wf.state_index);
        if (!wf.nto_role.empty()) meta.set("nto_role", wf.nto_role);
        if (wf.orbital_kind == "nto") meta.set("occupation_semantics", "transition_weight");
    };

    // spin (QVF spec § 5.3)
    const std::string spin =
        !wf.spin.empty() ? wf.spin : (wf.unrestricted ? "unrestricted" : "restricted");
    if (spin != "restricted" && spin != "unrestricted" && spin != "restricted_open" &&
        spin != "generalized")
        throw std::runtime_error("spin '" + spin +
                                 "' is not restricted, unrestricted, restricted_open or "
                                 "generalized");
    const bool unrestricted = spin == "unrestricted";
    if (wf.orbital_kind == "nto" && unrestricted) {
        if (wf.state_ref.empty() || wf.state_index <= 0 || wf.nto_role.empty())
            throw std::runtime_error("unrestricted NTOs require state_ref, state_index and nto_role");
        auto check_nto = [](const Matrix& c, const std::vector<double>& e,
                            const std::vector<double>& o) {
            if (c.rows < 0 || static_cast<int64_t>(o.size()) != c.rows ||
                (!e.empty() && static_cast<int64_t>(e.size()) != c.rows))
                throw std::runtime_error("NTOs need one weight and energy per coefficient row");
            for (size_t i = 0; i < o.size(); ++i)
                if (!std::isfinite(o[i]) || o[i] < 0 || (i && o[i - 1] < o[i]))
                    throw std::runtime_error("NTO weights must be finite, nonnegative and non-increasing");
        };
        check_nto(wf.mo_coefficients_alpha, wf.energies_alpha, wf.occupations_alpha);
        check_nto(wf.mo_coefficients_beta, wf.energies_beta, wf.occupations_beta);
    }
    if (wf.unrestricted && !unrestricted)
        throw std::runtime_error("unrestricted = true disagrees with spin '" + spin + "'");
    if (spin == "generalized" && !wf.complex_coefficients)
        throw std::runtime_error(
            "a generalized wavefunction needs complex_coefficients (QVF spec § 5.3)");
    const int64_t ao_extent = spin == "generalized" ? 2 * n_ao : n_ao;

    auto padded = [](const std::vector<double>& v, int64_t n) {
        std::vector<double> out(static_cast<size_t>(n), 0.0);
        for (size_t i = 0; i < out.size() && i < v.size(); ++i) out[i] = v[i];
        return out;
    };
    auto labels = [](const std::vector<std::string>& l, int64_t n, const std::string& what) {
        if (static_cast<int64_t>(l.size()) != n)
            throw std::runtime_error(what + " has " + std::to_string(l.size()) +
                                     " entries but there are " + std::to_string(n) +
                                     " MOs; there is one label per MO");
        return Json::from_strings(l);
    };

    if (wf.localization_reference_basis && wf.localization_reference_basis->empty())
        throw std::runtime_error("localization_reference_basis must be a non-empty string");
    Json mo_meta = Json::object();
    if (wf.localization_reference_basis)
        mo_meta.set("localization_reference_basis", *wf.localization_reference_basis);
    std::vector<std::pair<std::string, Tensor>> coefficients;
    double occupation_sum = 0.0;
    if (!unrestricted) {
        if (!wf.symmetry_labels_alpha.empty() || !wf.symmetry_labels_beta.empty())
            throw std::runtime_error(
                "symmetry_labels_alpha / _beta belong to an unrestricted wavefunction");
        Tensor C = coefficient_tensor(wf.mo_coefficients, wf.complex_coefficients,
                                      ao_extent, "mo_coefficients");
        const int64_t n_mo = wf.mo_coefficients.rows;
        const std::vector<double> occ = padded(wf.occupations, n_mo);
        mo_meta.set("n_mo", static_cast<long long>(n_mo))
               .set("n_ao", static_cast<long long>(n_ao))
               .set("spin", spin).set("orbital_kind", wf.orbital_kind)
               .set("energies", Json::from_doubles(padded(wf.energies, n_mo)))
               .set("occupations", Json::from_doubles(occ));
        if (!wf.symmetry_labels.empty())
            mo_meta.set("symmetry_labels", labels(wf.symmetry_labels, n_mo, "symmetry_labels"));
        for (double x : occ) occupation_sum += x;
        coefficients.emplace_back("mo_coefficients", std::move(C));
    } else {
        if (!wf.symmetry_labels.empty())
            throw std::runtime_error("an unrestricted wavefunction takes symmetry_labels_alpha "
                                     "and symmetry_labels_beta");
        Tensor Ca = coefficient_tensor(wf.mo_coefficients_alpha, wf.complex_coefficients,
                                       ao_extent, "mo_coefficients_alpha");
        Tensor Cb = coefficient_tensor(wf.mo_coefficients_beta, wf.complex_coefficients,
                                       ao_extent, "mo_coefficients_beta");
        auto channel = [&](int64_t n, const std::vector<double>& e,
                           const std::vector<double>& o, const std::vector<std::string>& l,
                           const char* what) {
            const std::vector<double> occ = padded(o, n);
            for (double x : occ) occupation_sum += x;
            Json c = Json::object();
            c.set("n_mo", static_cast<long long>(n))
             .set("energies", Json::from_doubles(padded(e, n)))
             .set("occupations", Json::from_doubles(occ));
            if (!l.empty()) c.set("symmetry_labels", labels(l, n, what));
            return c;
        };
        Json alpha = channel(wf.mo_coefficients_alpha.rows, wf.energies_alpha,
                             wf.occupations_alpha, wf.symmetry_labels_alpha,
                             "symmetry_labels_alpha");
        Json beta = channel(wf.mo_coefficients_beta.rows, wf.energies_beta,
                            wf.occupations_beta, wf.symmetry_labels_beta,
                            "symmetry_labels_beta");
        mo_meta.set("n_ao", static_cast<long long>(n_ao))
               .set("spin", "unrestricted").set("orbital_kind", wf.orbital_kind)
               .set("alpha", alpha).set("beta", beta);
        coefficients.emplace_back("mo_coefficients_alpha", std::move(Ca));
        coefficients.emplace_back("mo_coefficients_beta", std::move(Cb));
    }
    if (wf.has_kpoint) {
        Json k = Json::array();
        for (double x : wf.kpoint) k.push_back(Json(x));
        mo_meta.set("k_point", k);
    }
    if (wf.complex_coefficients) mo_meta.set("complex", true);
    add_state_fields(mo_meta);

    Json basis = Json::object();
    basis.set("structure_ref", wf.structure_ref).set("pure", wf.pure)
         .set("n_ao", static_cast<long long>(n_ao)).set("shells", shells_json);
    if (!wf.ecp.is_null()) basis.set("ecp", checked_ecp(wf.ecp));
    if (wf.n_electrons >= 0) {
        if (wf.orbital_kind != "nto" && std::fabs(occupation_sum - wf.n_electrons) > 1e-6)
            throw std::runtime_error("n_electrons is " + std::to_string(wf.n_electrons) +
                                     " but the occupations sum to " +
                                     std::to_string(occupation_sum) +
                                     "; they must agree (QVF spec § 5.3)");
        basis.set("n_electrons", wf.n_electrons);
    }

    std::vector<std::pair<std::string, Tensor>> matrices;
    if (!matrix_empty(wf.overlap))
        matrices.emplace_back("overlap", ao_matrix_tensor(wf.overlap, n_ao, "overlap"));
    if (!matrix_empty(wf.density_matrix)) {
        if (unrestricted)
            throw std::runtime_error("density_matrix accompanies mo_coefficients; an "
                                     "unrestricted wavefunction takes density_matrix_alpha "
                                     "and density_matrix_beta");
        matrices.emplace_back("density_matrix",
                              ao_matrix_tensor(wf.density_matrix, n_ao, "density_matrix"));
    }
    const bool has_da = !matrix_empty(wf.density_matrix_alpha);
    const bool has_db = !matrix_empty(wf.density_matrix_beta);
    if (has_da || has_db) {
        if (!unrestricted)
            throw std::runtime_error("density_matrix_alpha and density_matrix_beta accompany "
                                     "the alpha and beta coefficients");
        if (!(has_da && has_db))
            throw std::runtime_error(
                "density_matrix_alpha and density_matrix_beta come together");
        matrices.emplace_back("density_matrix_alpha",
                              ao_matrix_tensor(wf.density_matrix_alpha, n_ao,
                                               "density_matrix_alpha"));
        matrices.emplace_back("density_matrix_beta",
                              ao_matrix_tensor(wf.density_matrix_beta, n_ao,
                                               "density_matrix_beta"));
    }

    // Everything is checked: write the section.
    const std::string sid = reserve_id(wf.section_id);
    Json members = Json::object();
    members.set("basis", json_member(sid + "/basis.json", basis));
    members.set("mo_metadata", json_member(sid + "/mo_metadata.json", mo_meta));
    for (const auto& c : coefficients)
        members.set(c.first, binary_member(sid + "/" + c.first + ".dat", c.second));
    for (const auto& m : matrices)
        members.set(m.first, binary_member(sid + "/" + m.first + ".dat", m.second));

    Json section = Json::object();
    section.set("id", sid).set("kind", "wavefunction.gto").set("members", members);
    sections_.push_back(section);
    return sid;
}

size_t QvfWriter::wavefunction_index(const std::string& id) const {
    const size_t index = section_index(id);
    if (sections_.elements()[index].at("kind").as_string() != "wavefunction.gto")
        throw std::runtime_error("section '" + id + "' is not a wavefunction.gto section");
    return index;
}

Json QvfWriter::member_json(const Json& member) const {
    const std::string& path = member.at("path").as_string();
    for (const auto& f : files_)
        if (f.first == path) return Json::parse(std::string(f.second.begin(), f.second.end()));
    if (out_) {
        const auto cached = out_->json_cache.find(path);
        if (cached != out_->json_cache.end()) return Json::parse(cached->second.dump());
        if (out_->original) {
            const auto bytes = out_->original->read(path, container::Limits{}.max_json_bytes);
            return Json::parse(std::string(bytes.begin(), bytes.end()));
        }
    }
    throw std::runtime_error("no member file " + path);
}

void QvfWriter::replace_section(size_t index, Json section) {
    // Json has no mutable element access: rebuild the array.
    const std::vector<Json>& sections = sections_.elements();
    Json rebuilt = Json::array();
    for (size_t i = 0; i < sections.size(); ++i)
        rebuilt.push_back(i == index ? section : sections[i]);
    sections_ = std::move(rebuilt);
}

std::string QvfWriter::set_wavefunction_overlap(const std::string& section_id,
                                                const Matrix& overlap) {
    const size_t index = wavefunction_index(section_id);
    Json section = sections_.elements()[index];
    Json members = section.at("members");
    if (members.has("overlap"))
        throw std::runtime_error("section '" + section_id + "' already has an overlap member");
    const int64_t n_ao = member_json(members.at("basis")).at("n_ao").as_int();
    const Tensor t = ao_matrix_tensor(overlap, n_ao, "overlap");
    members.set("overlap", binary_member(section_id + "/overlap.dat", t));
    section.set("members", members);
    replace_section(index, std::move(section));
    return section_id;
}

std::string QvfWriter::set_wavefunction_ecp(const std::string& section_id, const Json& ecp) {
    const size_t index = wavefunction_index(section_id);
    const Json checked = checked_ecp(ecp);
    Json section = sections_.elements()[index];
    Json members = section.at("members");
    Json basis = member_json(members.at("basis"));
    basis.set("ecp", checked);
    const std::string path = members.at("basis").at("path").as_string();
    drop_file(path);
    members.set("basis", json_member(path, basis));
    section.set("members", members);
    replace_section(index, std::move(section));
    return section_id;
}

std::string QvfWriter::add_wavefunction_bloch(const WavefunctionBloch& wf) {
    check_new_id(wf.section_id);
    int64_t n_spin = 0;
    if (wf.spin == "restricted") n_spin = 1;
    else if (wf.spin == "unrestricted") n_spin = 2;
    else throw std::runtime_error("wavefunction.bloch spin must be 'restricted' or "
                                  "'unrestricted', got '" + wf.spin + "'");
    const int64_t n_k = static_cast<int64_t>(wf.kpoints.size());
    if (n_k < 1) throw std::runtime_error("wavefunction.bloch needs at least one k-point");
    if (wf.n_mo < 0) throw std::runtime_error("wavefunction.bloch n_mo must be >= 0");
    int64_t n_ao = 0;
    Json shells_json = basis_shells_json(wf.shells, wf.coeffs_are_libint_normalized, n_ao);
    const size_t n_orb = static_cast<size_t>(n_k) * static_cast<size_t>(n_spin) *
                         static_cast<size_t>(wf.n_mo);
    if (wf.mo_coefficients.size() != n_orb * static_cast<size_t>(n_ao) * 2)
        throw std::runtime_error("wavefunction.bloch mo_coefficients must hold n_k * n_spin * "
                                 "n_mo * n_ao * 2 doubles");
    if (wf.energies.size() != n_orb || wf.occupations.size() != n_orb)
        throw std::runtime_error("wavefunction.bloch energies and occupations must hold "
                                 "n_k * n_spin * n_mo doubles");
    Json kp = Json::object();
    Json points = Json::array();
    for (const auto& k : wf.kpoints) {
        Json v = Json::array();
        for (double x : k) v.push_back(Json(x));
        points.push_back(v);
    }
    kp.set("kpoints", points);
    if (!wf.weights.empty()) {
        double sum = 0.0;
        bool negative = false;
        for (double x : wf.weights) {
            sum += x;
            negative = negative || !(x >= 0.0);
        }
        if (static_cast<int64_t>(wf.weights.size()) != n_k || negative ||
            !(std::fabs(sum - 1.0) <= 1e-6))
            throw std::runtime_error("wavefunction.bloch weights must be n_k non-negative "
                                     "numbers that sum to 1 within 1e-6");
        kp.set("weights", Json::from_doubles(wf.weights));
    }
    if (!wf.labels.empty()) {
        if (static_cast<int64_t>(wf.labels.size()) != n_k)
            throw std::runtime_error("wavefunction.bloch labels must have n_k entries");
        kp.set("labels", Json::from_strings(wf.labels));
    }
    if (!wf.mesh.empty()) {
        if (wf.mesh.size() != 3 || wf.mesh[0] < 1 || wf.mesh[1] < 1 || wf.mesh[2] < 1)
            throw std::runtime_error("wavefunction.bloch mesh must be three integers >= 1");
        kp.set("mesh", Json::from_ints(wf.mesh));
    }
    if (wf.has_shift) {
        Json v = Json::array();
        for (double x : wf.shift) v.push_back(Json(x));
        kp.set("shift", v);
    }
    if (wf.has_time_reversal) kp.set("time_reversal", wf.time_reversal);
    Json meta = Json::object();
    meta.set("spin", wf.spin).set("n_spin", static_cast<long long>(n_spin))
        .set("n_mo", static_cast<long long>(wf.n_mo))
        .set("n_ao", static_cast<long long>(n_ao))
        .set("orbital_kind", wf.orbital_kind).set("complex", true)
        .set("k_point_convention", "fractional");
    if (!wf.n_mo_per_k.empty()) {
        if (static_cast<int64_t>(wf.n_mo_per_k.size()) != n_k)
            throw std::runtime_error("wavefunction.bloch n_mo_per_k must have n_k entries");
        Json per_k = Json::array();
        for (int64_t n : wf.n_mo_per_k) {
            if (n < 0 || n > wf.n_mo)
                throw std::runtime_error("wavefunction.bloch n_mo_per_k entries must lie "
                                         "between 0 and n_mo");
            per_k.push_back(Json(static_cast<long long>(n)));
        }
        meta.set("n_mo_per_k", per_k);
    }
    if (wf.has_fermi_energy) {
        Json fermi = Json::object();
        fermi.set("value", wf.fermi_energy_ev).set("units", "ev");
        meta.set("fermi_energy", fermi);
    }
    Json basis = Json::object();
    basis.set("structure_ref", wf.structure_ref).set("pure", wf.pure)
         .set("n_ao", static_cast<long long>(n_ao)).set("shells", shells_json);

    const std::string sid = reserve_id(wf.section_id);
    Json members = Json::object();
    members.set("basis", json_member(sid + "/basis.json", basis));
    members.set("kpoints", json_member(sid + "/kpoints.json", kp));
    members.set("mo_metadata", json_member(sid + "/mo_metadata.json", meta));
    members.set("mo_coefficients",
                binary_member(sid + "/mo_coefficients.bin",
                              tensor_f64(wf.mo_coefficients, {n_k, n_spin, wf.n_mo, n_ao, 2})));
    members.set("energies", binary_member(sid + "/energies.bin",
                                          tensor_f64(wf.energies, {n_k, n_spin, wf.n_mo})));
    members.set("occupations",
                binary_member(sid + "/occupations.bin",
                              tensor_f64(wf.occupations, {n_k, n_spin, wf.n_mo})));
    Json section = Json::object();
    section.set("id", sid).set("kind", "wavefunction.bloch").set("members", members);
    sections_.push_back(section);
    return sid;
}

// -- spectra ----------------------------------------------------------------

std::string QvfWriter::add_spectrum(const std::string& kind, const Spectrum& spectrum,
                                    const std::string& label,
                                    const std::string& section_id) {
    Json payload = Json::object();
    payload.set("frequencies", Json::from_doubles(spectrum.frequencies));
    payload.set("intensities", Json::from_doubles(spectrum.intensities));
    if (!spectrum.x_units.empty()) payload.set("x_units", spectrum.x_units);
    if (!spectrum.y_units.empty()) payload.set("y_units", spectrum.y_units);
    return add_spectrum(kind, payload, label, section_id);
}

std::string QvfWriter::add_spectrum(const std::string& kind, const Json& payload,
                                    const std::string& label,
                                    const std::string& section_id) {
    if (is_spectrum_1d(kind)) {
        check_payload_units(kind, payload, {}, true);
        // states_ref names an added excited_states section (spec § 5.4). A
        // Json::raw fragment is parsed for the check; one that does not parse
        // is left alone (it is written verbatim, as before).
        Json parsed;
        const Json* doc = &payload;
        if (payload.type() == Json::Type::Raw) {
            try {
                parsed = Json::parse(payload.as_string());
            } catch (const std::exception&) {
                parsed = Json();
            }
            doc = &parsed;
        }
        if (doc->is_object() && doc->has("states_ref")) {
            const Json& ref = doc->at("states_ref");
            bool found = false;
            if (ref.is_string())
                for (const Json& s : sections_.elements())
                    if (s.at("id").as_string() == ref.as_string() &&
                        s.at("kind").as_string() == "excited_states")
                        found = true;
            if (!found)
                throw std::runtime_error(kind + ": states_ref " + ref.dump(0) +
                                         " must name an added excited_states section "
                                         "(QVF spec § 5.4)");
        }
    }
    std::string sid = section_id;
    if (sid.empty()) {
        sid = kind;
        for (char& c : sid) if (c == '.') c = '_';
    }
    check_new_id(sid);
    Json members = Json::object();
    members.set("spectrum", json_member("spectra/" + sid + ".json", payload));
    Json section = Json::object();
    section.set("id", reserve_id(sid)).set("kind", kind).set("members", members);
    if (!label.empty()) section.set("label", label);
    sections_.push_back(section);
    return sid;
}

// -- analysis / provenance --------------------------------------------------

// The atom_properties member roles, in schema order ($defs/SectionAtomProperties,
// QVF spec § 5.7). The schema, the spec, both reference writers and the C ABI
// share this one list; the scalar roles are float64 [n_atoms], the two vector
// roles float64 [n_atoms, 3].
static const char* const kAtomPropertyRoles[] = {
    "mulliken_charge", "loewdin_charge", "hirshfeld_charge", "iao_charge",
    "spin_population",
    "cm5_charge", "chelpg_charge", "resp_charge", "esp_charge", "nbo_charge",
    "mbis_charge", "bader_charge", "ddec_charge",
    "valence", "bonded_valence", "fukui_plus", "fukui_minus", "fukui_zero",
    "dual_descriptor", "atomic_energy_eh", "atomic_volume_bohr3",
    "nmr_shielding_iso_ppm",
    "magnetic_moment", "force_eh_bohr"};
// The vector roles of the list above.
static const char* const kAtomVectorRoles[] = {"magnetic_moment", "force_eh_bohr"};

static bool is_atom_vector_role(const std::string& role) {
    for (const char* r : kAtomVectorRoles)
        if (role == r) return true;
    return false;
}

static bool is_custom_role(const std::string& role) {
    return role.compare(0, 7, "custom.") == 0;
}

// custom.<name> with <name> one or more characters of [A-Za-z0-9_].
static bool is_valid_custom_role(const std::string& role) {
    if (!is_custom_role(role) || role.size() == 7) return false;
    for (size_t i = 7; i < role.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(role[i]);
        if (!(std::isalnum(c) || c == '_')) return false;
    }
    return true;
}

// Throws unless `role` can carry a scalar (vector == false) or a vector member:
// a listed role of that shape, or a custom.<name> role whose entry in
// `registry` (set_atom_properties_meta) has that kind.
static void check_atom_role(const std::string& role, bool vector,
                            const Json& registry) {
    const std::string shape = vector ? "vector" : "scalar";
    if (is_custom_role(role)) {
        if (!is_valid_custom_role(role))
            throw std::runtime_error("atom_properties custom role '" + role +
                                     "' must be custom.<name> with <name> in "
                                     "[A-Za-z0-9_]+");
        if (!registry.has(role))
            throw std::runtime_error(
                "atom_properties " + role + " has no meta entry; describe it with "
                "set_atom_properties_meta first (QVF spec § 5.7)");
        const std::string& kind = registry.at(role).at("kind").as_string();
        if (kind != shape)
            throw std::runtime_error("atom_properties " + role + " is described as a " +
                                     kind + " in its meta entry but is written as a " +
                                     shape);
        return;
    }
    bool known = false;
    std::string expected;
    for (const char* r : kAtomPropertyRoles) {
        known = known || role == r;
        expected += expected.empty() ? r : std::string(", ") + r;
    }
    if (!known)
        throw std::runtime_error("unknown atom_properties role '" + role +
                                 "' (expected one of " + expected +
                                 ", or custom.<name>)");
    if (is_atom_vector_role(role) != vector)
        throw std::runtime_error(
            "atom_properties role '" + role + "' is a " +
            (vector ? "scalar role; use add_atom_property" :
                      "vector role; use add_atom_vector_property"));
}

static Tensor atom_property_tensor(const std::vector<double>& flat, size_t n_atoms,
                                   bool vector) {
    std::vector<int64_t> shape = {static_cast<int64_t>(n_atoms)};
    if (vector) shape.push_back(3);
    return tensor_f64(flat, shape);
}

void QvfWriter::set_atom_properties_meta(const Json& meta) {
    Json doc = meta;
    if (doc.type() == Json::Type::Raw) {
        try {
            doc = Json::parse(meta.as_string());
        } catch (const std::exception& e) {
            throw std::runtime_error(
                std::string("atom_properties meta is not valid JSON: ") + e.what());
        }
    }
    if (!doc.is_object())
        throw std::runtime_error(
            "atom_properties meta must be a JSON object mapping custom.<name> to "
            "{units, description, kind}");
    Json merged = atom_properties_meta_;
    for (const auto& kv : doc.items()) {
        const std::string& role = kv.first;
        if (!is_valid_custom_role(role))
            throw std::runtime_error("atom_properties meta key '" + role +
                                     "' must be custom.<name> with <name> in "
                                     "[A-Za-z0-9_]+");
        const Json& e = kv.second;
        if (!e.is_object())
            throw std::runtime_error("atom_properties meta entry '" + role +
                                     "' must be a JSON object");
        Json entry = Json::object();
        for (const char* key : {"units", "description"}) {
            if (!e.has(key) || !e.at(key).is_string() || e.at(key).as_string().empty())
                throw std::runtime_error("atom_properties meta entry '" + role +
                                         "' needs a non-empty string '" + key + "'");
            entry.set(key, e.at(key));
        }
        if (!e.has("kind") || !e.at("kind").is_string() ||
            (e.at("kind").as_string() != "scalar" && e.at("kind").as_string() != "vector"))
            throw std::runtime_error("atom_properties meta entry '" + role +
                                     "' needs 'kind' \"scalar\" or \"vector\"");
        entry.set("kind", e.at("kind"));
        for (const auto& field : e.items()) {
            if (field.first == "units" || field.first == "description" ||
                field.first == "kind")
                continue;
            if (field.first.compare(0, 2, "x_") != 0)
                throw std::runtime_error("atom_properties meta entry '" + role +
                                         "' has unknown key '" + field.first +
                                         "' (allowed: units, description, kind, x_*)");
            entry.set(field.first, field.second);
        }
        merged.set(role, entry);
    }
    atom_properties_meta_ = std::move(merged);
}

Json QvfWriter::atom_properties_meta_payload(
    const std::vector<std::string>& custom_roles) const {
    Json payload = Json::object();
    for (const std::string& role : custom_roles)
        payload.set(role, atom_properties_meta_.at(role));
    return payload;
}

std::string QvfWriter::add_atom_properties(const std::vector<double>& mulliken,
                                           const std::vector<double>& loewdin,
                                           const std::vector<double>& spin,
                                           const std::string& section_id) {
    std::map<std::string, std::vector<double>> roles;
    if (!mulliken.empty()) roles["mulliken_charge"] = mulliken;
    if (!loewdin.empty()) roles["loewdin_charge"] = loewdin;
    if (!spin.empty()) roles["spin_population"] = spin;
    return add_atom_properties(roles, section_id);
}

std::string QvfWriter::add_atom_properties(std::initializer_list<double> mulliken,
                                           const std::vector<double>& loewdin,
                                           const std::vector<double>& spin,
                                           const std::string& section_id) {
    return add_atom_properties(std::vector<double>(mulliken), loewdin, spin,
                               section_id);
}

std::string QvfWriter::add_atom_properties(
    const std::map<std::string, std::vector<double>>& roles,
    const std::string& section_id) {
    check_new_id(section_id);
    if (roles.empty())
        throw std::runtime_error("atom_properties needs at least one array");
    const size_t n_atoms = roles.begin()->second.size();
    for (const auto& kv : roles) {
        check_atom_role(kv.first, false, atom_properties_meta_);
        if (kv.second.empty())
            throw std::runtime_error("atom_properties " + kv.first + " has no values");
        if (kv.second.size() != n_atoms)
            throw std::runtime_error("atom_properties roles differ in length");
    }
    // Listed roles in schema order, whatever order the map holds them in, then
    // the custom roles by name, then the `meta` member that describes them.
    Json members = Json::object();
    for (const char* role : kAtomPropertyRoles) {
        const auto it = roles.find(role);
        if (it != roles.end())
            members.set(role, binary_member("atom_properties/" + std::string(role) +
                                                ".bin",
                                            atom_property_tensor(it->second, n_atoms,
                                                                 false)));
    }
    std::vector<std::string> custom_roles;
    for (const auto& kv : roles) {
        if (!is_custom_role(kv.first)) continue;
        custom_roles.push_back(kv.first);
        members.set(kv.first, binary_member("atom_properties/" + kv.first + ".bin",
                                            atom_property_tensor(kv.second, n_atoms,
                                                                 false)));
    }
    if (!custom_roles.empty())
        members.set("meta", json_member("atom_properties/meta.json",
                                        atom_properties_meta_payload(custom_roles)));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "atom_properties")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_atom_property(const std::string& role,
                                         const std::vector<double>& values,
                                         const std::string& section_id) {
    return add_atom_role(role, values, values.size(), false, section_id);
}

std::string QvfWriter::add_atom_vector_property(
    const std::string& role, const std::vector<std::array<double, 3>>& values,
    const std::string& section_id) {
    std::vector<double> flat;
    flat.reserve(values.size() * 3);
    for (const auto& v : values) flat.insert(flat.end(), v.begin(), v.end());
    return add_atom_role(role, flat, values.size(), true, section_id);
}

std::string QvfWriter::add_atom_role(const std::string& role,
                                     const std::vector<double>& flat, size_t n_atoms,
                                     bool vector, const std::string& section_id) {
    check_atom_role(role, vector, atom_properties_meta_);
    if (n_atoms == 0)
        throw std::runtime_error("atom_properties " + role + " has no values");
    const bool custom = is_custom_role(role);
    const std::string path = "atom_properties/" + role + ".bin";
    const std::vector<Json>& sections = sections_.elements();
    size_t index = 0;
    while (index < sections.size() &&
           sections[index].at("id").as_string() != section_id)
        ++index;

    if (index == sections.size()) {
        check_new_id(section_id);
        Json members = Json::object();
        members.set(role, binary_member(path, atom_property_tensor(flat, n_atoms, vector)));
        if (custom)
            members.set("meta", json_member("atom_properties/meta.json",
                                            atom_properties_meta_payload({role})));
        Json section = Json::object();
        section.set("id", reserve_id(section_id)).set("kind", "atom_properties")
               .set("members", members);
        sections_.push_back(section);
        return section_id;
    }

    const Json& existing = sections[index];
    if (existing.at("kind").as_string() != "atom_properties")
        throw std::runtime_error("section '" + section_id +
                                 "' is not an atom_properties section");
    Json members = existing.at("members");
    if (members.has(role))
        throw std::runtime_error("atom_properties " + role +
                                 " is already in section '" + section_id + "'");
    std::vector<std::string> custom_roles;
    for (const auto& kv : members.items()) {
        if (kv.first == "meta") continue;
        const int64_t have = kv.second.at("shape").at(0).as_int();
        if (static_cast<int64_t>(n_atoms) != have)
            throw std::runtime_error(
                "atom_properties " + role + " has " + std::to_string(n_atoms) +
                " values but section '" + section_id + "' has " +
                std::to_string(have) + " atoms");
        if (is_custom_role(kv.first)) {
            // The meta entries are the writer's current ones: each must still
            // describe the role it was written for.
            check_atom_role(kv.first, kv.second.at("shape").size() == 2,
                            atom_properties_meta_);
            custom_roles.push_back(kv.first);
        }
    }
    // Everything is checked: change the section. `meta` stays the last member.
    const bool had_meta = members.has("meta");
    Json meta_spec = had_meta ? members.at("meta") : Json();
    members.erase("meta");
    members.set(role, binary_member(path, atom_property_tensor(flat, n_atoms, vector)));
    if (custom) {
        custom_roles.push_back(role);
        if (had_meta && meta_spec.has("path") && meta_spec.at("path").is_string())
            drop_file(meta_spec.at("path").as_string());
        meta_spec = json_member("atom_properties/meta.json",
                                atom_properties_meta_payload(custom_roles));
        members.set("meta", meta_spec);
    } else if (had_meta) {
        members.set("meta", meta_spec);
    }
    // Json has no mutable element access: rebuild the array with the section
    // replaced. A handful of roles makes this cheap.
    Json updated = existing;
    updated.set("members", members);
    Json rebuilt = Json::array();
    for (size_t i = 0; i < sections.size(); ++i)
        rebuilt.push_back(i == index ? updated : sections[i]);
    sections_ = std::move(rebuilt);
    return section_id;
}

std::string QvfWriter::add_scf_history(const Json& iterations,
                                       const std::string& section_id) {
    check_new_id(section_id);
    Json payload = Json::object();
    payload.set("iterations", iterations);
    Json members = Json::object();
    members.set("iterations", json_member(section_id + "/iterations.json", payload));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "scf_history")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_citations(const std::string& bibtex,
                                     const std::string& section_id) {
    check_new_id(section_id);
    Json members = Json::object();
    members.set("references",
                opaque_member(section_id + "/references.bib",
                              std::vector<uint8_t>(bibtex.begin(), bibtex.end())));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "citations")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_run_record(const RunRecord& record,
                                      const std::string& section_id) {
    check_new_id(section_id);
    if (record.program.empty())
        throw std::runtime_error("run.record requires a non-empty program");
    if (record.input_text.empty() && record.log_text.empty())
        throw std::runtime_error(
            "run.record requires input_text and/or log_text");
    // A job set (spec § 5.8, § 5.9): the job this run executed, and its outcome.
    if (!record.job_ref.empty())
        require_section_kind("run.record job_ref", record.job_ref, "job.spec");
    if (!record.run_status.empty() && record.run_status != "running" &&
        record.run_status != "converged" && record.run_status != "failed")
        throw std::runtime_error("run.record run_status must be one of running, converged, "
                                 "failed, got '" + record.run_status + "'");
    Json members = Json::object();
    auto opaque = [&](const std::string& role, const std::string& path,
                      std::vector<uint8_t> bytes) {
        members.set(role, opaque_member(path, std::move(bytes)));
    };
    if (!record.input_text.empty())
        opaque("input", section_id + "/input.txt",
               std::vector<uint8_t>(record.input_text.begin(),
                                    record.input_text.end()));
    if (!record.log_text.empty())
        opaque("log", section_id + "/log.txt",
               std::vector<uint8_t>(record.log_text.begin(),
                                    record.log_text.end()));
    for (const auto& att : record.attachments)
        opaque("attachment." + att.first,
               section_id + "/attachments/" + att.first, att.second);
    if (record.files.type() == Json::Type::Object)
        members.set("files",
                    json_member(section_id + "/files.json", record.files));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "run.record")
           .set("program", record.program).set("members", members);
    if (!record.program_version.empty())
        section.set("program_version", record.program_version);
    if (!record.command.empty()) section.set("command", record.command);
    if (record.has_exit_status) section.set("exit_status", record.exit_status);
    if (!record.started_utc.empty())
        section.set("started_utc", record.started_utc);
    if (!record.finished_utc.empty())
        section.set("finished_utc", record.finished_utc);
    if (record.sequence >= 0) section.set("sequence", record.sequence);
    if (!record.job_ref.empty()) section.set("job_ref", record.job_ref);
    if (!record.run_status.empty()) section.set("run_status", record.run_status);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_job_spec(const Json& spec,
                                    const std::string& section_id) {
    check_new_id(section_id);
    if (spec.type() != Json::Type::Object)
        throw std::runtime_error("job.spec spec must be a JSON object");
    if (!spec.has("job_type") || !spec.at("job_type").is_string())
        throw std::runtime_error("job.spec requires a job_type string");
    const std::string& jt = spec.at("job_type").as_string();
    if (jt != "molecular" && jt != "periodic")
        throw std::runtime_error(
            "job.spec job_type must be 'molecular' or 'periodic', got: " + jt);
    // v0.3: null is explicit "none", absence leaves the choice to the engine.
    for (const char* key : {"basis", "functional"})
        if (spec.has(key)) {
            const Json value = resolved(spec.at(key), std::string("job.spec ") + key);
            if (!value.is_null() && !value.is_string())
                throw std::runtime_error(std::string("job.spec ") + key + " must be a string or null");
        }
    // The typed fields added by spec § 5.9; everything is checked before the
    // id is reserved, so a refused spec adds nothing. The root stays open.
    if (spec.has("structure_ref")) {
        const std::string what = "job.spec structure_ref";
        const Json ref_json = resolved(spec.at("structure_ref"), what);
        const std::string& ref = checked_text(ref_json, what);
        if (!is_valid_id(ref))
            throw std::runtime_error(what + " '" + ref +
                                     "' is not a section id ([A-Za-z0-9_.-]+)");
        for (const Json& s : sections_.elements())
            if (s.at("id").as_string() == ref && s.at("kind").as_string() != "structure")
                throw std::runtime_error(what + " '" + ref + "' names a '" +
                                         s.at("kind").as_string() +
                                         "' section, not a 'structure' section "
                                         "(QVF spec § 5.9)");
    }
    for (const char* key : {"dispersion", "relativistic"})
        if (spec.has(key))
            checked_text(resolved(spec.at(key), std::string("job.spec ") + key),
                         std::string("job.spec ") + key);
    if (spec.has("solvation"))
        check_block("solvation", spec.at("solvation"),
                    {{"model", Rule::Text}, {"solvent", Rule::Text},
                     {"epsilon", Rule::PositiveNumber}},
                    "model");
    if (spec.has("scf"))
        check_block("scf", spec.at("scf"),
                    {{"convergence", Rule::PositiveNumber},
                     {"max_iterations", Rule::Count}, {"guess", Rule::Text},
                     {"damping", Rule::Number}},
                    nullptr);
    if (spec.has("geometry"))
        check_block("geometry", spec.at("geometry"),
                    {{"convergence", Rule::PositiveNumber},
                     {"max_steps", Rule::Count}, {"coordinates", Rule::Coordinates}},
                    nullptr);
    if (spec.has("resources"))
        check_block("resources", spec.at("resources"),
                    {{"n_cores", Rule::Count}, {"memory_mb", Rule::Count},
                     {"walltime_s", Rule::PositiveNumber}},
                    nullptr);
    Json members = Json::object();
    members.set("spec", json_member(section_id + "/spec.json", spec));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "job.spec");
    // The section-level structure_ref, which a reader finds in the manifest
    // (spec § 5.9, "Job sets"); the payload's stays as runners read it.
    if (spec.has("structure_ref"))
        section.set("structure_ref",
                    resolved(spec.at("structure_ref"), "job.spec structure_ref"));
    section.set("members", members);
    sections_.push_back(section);
    return section_id;
}

// -- vendor -----------------------------------------------------------------

std::string QvfWriter::add_vendor_section(const std::string& kind,
                                          const Json& json_members, bool critical,
                                          const std::string& schema_uri,
                                          const std::string& section_id) {
    return add_vendor_section(kind, json_members, std::vector<OpaqueMember>{}, critical,
                              schema_uri, section_id);
}

std::string QvfWriter::add_vendor_section(const std::string& kind,
                                          const Json& json_members,
                                          const std::vector<OpaqueMember>& opaque_members,
                                          bool critical, const std::string& schema_uri,
                                          const std::string& section_id) {
    if (kind.rfind("x_", 0) != 0)
        throw std::runtime_error("vendor kind must start with x_: " + kind);
    std::string sid = section_id;
    if (sid.empty()) {
        sid = kind;
        for (char& c : sid) if (c == '.') c = '_';
    }
    check_new_id(sid);
    // Every role once; the opaque members' media types before anything is
    // written.
    std::set<std::string> roles;
    if (!json_members.is_null())
        for (const auto& kv : json_members.items()) roles.insert(kv.first);
    for (const OpaqueMember& om : opaque_members) {
        if (om.role.empty())
            throw std::runtime_error("vendor section: an opaque member has no role");
        if (!roles.insert(om.role).second)
            throw std::runtime_error("vendor section: member role '" + om.role +
                                     "' is given twice");
        if (!om.media_type.empty() && !is_media_type(om.media_type))
            throw std::runtime_error("vendor section member '" + om.role + "': media type '" +
                                     om.media_type + "' is not type/subtype (QVF spec § 4.2)");
    }
    // json_members is an object mapping role -> JSON payload; each becomes a
    // JSON member under sid/<role>.json, each opaque member sid/<role>.bin.
    Json members = Json::object();
    if (!json_members.is_null())
        for (const auto& kv : json_members.items())
            members.set(kv.first, json_member(sid + "/" + kv.first + ".json", kv.second));
    for (const OpaqueMember& om : opaque_members)
        members.set(om.role,
                    opaque_member(sid + "/" + om.role + ".bin", om.bytes, om.media_type));
    if (members.size() == 0)
        throw std::runtime_error("vendor section needs at least one member");
    Json section = Json::object();
    section.set("id", reserve_id(sid)).set("kind", kind).set("members", members);
    if (critical) section.set("critical", true);
    if (!schema_uri.empty()) section.set("schema_uri", schema_uri);
    sections_.push_back(section);
    return sid;
}

std::string QvfWriter::add_opaque_member(const std::string& section_id, const std::string& role,
                                         const std::vector<uint8_t>& bytes,
                                         const std::string& media_type) {
    const size_t index = section_index(section_id);
    Json section = sections_.elements()[index];
    const Json& kind = section.at("kind");
    if (!kind.is_string() || kind.as_string().rfind("x_", 0) != 0)
        throw std::runtime_error("add_opaque_member: section '" + section_id +
                                 "' is not a vendor (x_<vendor>) section (QVF spec § 7.3)");
    if (role.empty()) throw std::runtime_error("add_opaque_member: empty role");
    Json members = section.at("members");
    if (members.has(role))
        throw std::runtime_error("section '" + section_id + "' already has a member '" +
                                 role + "'");
    members.set(role, opaque_member(section_id + "/" + role + ".bin", bytes, media_type));
    section.set("members", members);
    replace_section(index, std::move(section));
    return section_id;
}

// -- preview (spec § 5.15) ----------------------------------------------------

namespace {

const char kPngSignature[] = "\x89PNG\r\n\x1a\n";   // 8 bytes
const char kJpegSignature[] = "\xff\xd8\xff";       // 3 bytes

bool starts_with_bytes(const std::vector<uint8_t>& data, const char* signature, size_t n) {
    return data.size() >= n && std::memcmp(data.data(), signature, n) == 0;
}

bool contains_svg_tag(std::vector<uint8_t>::const_iterator first,
                      std::vector<uint8_t>::const_iterator last) {
    static const char kSvg[] = "<svg";
    return std::search(first, last, kSvg, kSvg + 4) != last;
}

// Strict UTF-8 (no overlongs, surrogates or code points above U+10FFFF), as
// the validator checks an image/svg+xml member.
bool is_strict_utf8(const std::vector<uint8_t>& s) {
    size_t i = 0;
    const size_t n = s.size();
    while (i < n) {
        const uint8_t c = s[i];
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
        if (cp < min || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF)) return false;
        i += len;
    }
    return true;
}

// Why `data` are not what `media_type` says (spec § 5.15), or "".
std::string preview_bytes_problem(const std::vector<uint8_t>& data,
                                  const std::string& media_type) {
    if (media_type == "image/png" && !starts_with_bytes(data, kPngSignature, 8))
        return "image/png bytes do not begin with the PNG signature";
    if (media_type == "image/jpeg" && !starts_with_bytes(data, kJpegSignature, 3))
        return "image/jpeg bytes do not begin with FF D8 FF";
    if (media_type == "image/svg+xml") {
        if (!is_strict_utf8(data)) return "image/svg+xml bytes do not decode as UTF-8";
        if (!contains_svg_tag(data.begin(), data.end()))
            return "image/svg+xml bytes contain no <svg element";
    }
    return "";
}

// The media type of preview bytes: `given` when not empty, checked against
// the bytes, else sniffed from them (spec § 5.15).
std::string preview_media_type(const std::vector<uint8_t>& data,
                               std::initializer_list<const char*> allowed,
                               const std::string& what, const std::string& given) {
    std::string media_type = given;
    bool svg_allowed = false;
    std::string names;
    for (const char* a : allowed) {
        svg_allowed = svg_allowed || std::strcmp(a, "image/svg+xml") == 0;
        names += (names.empty() ? "" : ", ") + std::string(a);
    }
    if (media_type.empty()) {
        if (starts_with_bytes(data, kPngSignature, 8))
            media_type = "image/png";
        else if (starts_with_bytes(data, kJpegSignature, 3))
            media_type = "image/jpeg";
        else if (svg_allowed &&
                 contains_svg_tag(data.begin(),
                                  data.begin() + static_cast<std::ptrdiff_t>(
                                                     std::min<size_t>(data.size(), 4096))))
            media_type = "image/svg+xml";
        else
            throw std::runtime_error(what + ": not PNG or JPEG bytes; pass " + what +
                                     "_media_type to say what they are");
    }
    bool listed = false;
    for (const char* a : allowed) listed = listed || media_type == a;
    if (!listed)
        throw std::runtime_error(what + ": media type '" + media_type + "' is not one of " +
                                 names + " (QVF spec § 5.15)");
    const std::string problem = preview_bytes_problem(data, media_type);
    if (!problem.empty()) throw std::runtime_error(what + ": " + problem);
    return media_type;
}

const char* image_extension(const std::string& media_type) {
    return media_type == "image/png" ? "png" : media_type == "image/jpeg" ? "jpg" : "svg";
}

}  // namespace

std::string QvfWriter::add_preview(const std::vector<uint8_t>& thumbnail,
                                   const std::vector<uint8_t>& image, const std::string& label,
                                   const std::string& section_id,
                                   const std::string& thumbnail_media_type,
                                   const std::string& image_media_type) {
    const std::string thumb_type = preview_media_type(
        thumbnail, {"image/png", "image/jpeg"}, "thumbnail", thumbnail_media_type);
    std::string image_type;
    if (!image.empty())
        image_type = preview_media_type(image, {"image/svg+xml", "image/png", "image/jpeg"},
                                        "image", image_media_type);
    check_new_id(section_id);
    // Already compressed images gain nothing from DEFLATE: stored as they are.
    Json members = Json::object();
    members.set("thumbnail",
                opaque_member(section_id + "/thumbnail." + image_extension(thumb_type),
                              thumbnail, thumb_type, /*store=*/true));
    if (!image.empty()) {
        try {
            members.set("image",
                        opaque_member(section_id + "/image." + image_extension(image_type),
                                      image, image_type, /*store=*/true));
        } catch (...) {
            drop_file(members.at("thumbnail").at("path").as_string());
            throw;
        }
    }
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "preview");
    if (!label.empty()) section.set("label", label);
    section.set("members", members);
    sections_.push_back(section);
    return section_id;
}

// -- volumes & basis AO -----------------------------------------------------

static Json grid_json(const Grid& g, const std::vector<int64_t>& shape) {
    Json grid = Json::object();
    Json origin = Json::array();
    for (double x : g.origin) origin.push_back(Json(x));
    Json vox = Json::array();
    for (const auto& v : g.voxel_vectors) {
        Json r = Json::array();
        for (double x : v) r.push_back(Json(x));
        vox.push_back(r);
    }
    Json sh = Json::array();
    for (int64_t d : shape) sh.push_back(Json(static_cast<long long>(d)));
    grid.set("origin", origin).set("voxel_vectors", vox).set("shape", sh);
    if (!g.units.empty()) {
        Json units = Json::object();
        for (const auto& kv : g.units) units.set(kv.first, kv.second);
        check_units_object("grid", units,
                           {{"origin", Unit::Bohr}, {"voxel_vectors", Unit::Bohr}});
        grid.set("units", units);
    }
    // Periodic grid semantics (spec § 5.2), written only when set.
    if (g.periodic) {
        Json flags = Json::array();
        for (bool p : *g.periodic) flags.push_back(Json(p));
        grid.set("periodic", flags);
    }
    if (g.includes_endpoint) grid.set("includes_endpoint", *g.includes_endpoint);
    if (!g.structure_ref.empty()) grid.set("structure_ref", g.structure_ref);
    return grid;
}

static void merge_peers(Json& section, const Json& extra) {
    if (extra.type() != Json::Type::Object) return;
    for (const auto& kv : extra.items()) section.set(kv.first, kv.second);
}

// The extent the trailing data axis of a 4-D volume kind must have (spec
// § 5.2): 3 for volume.vector, 2 for volume.complex; 0 for any other kind.
static int64_t volume_trailing_axis(const std::string& kind) {
    if (kind == "volume.vector") return 3;
    if (kind == "volume.complex") return 2;
    return 0;
}

static void check_volume_dtype(const std::string& kind, const Tensor& data) {
    if (data.dtype != DType::Float32 && data.dtype != DType::Float64)
        throw std::runtime_error(kind + ": data must be float32 or float64");
}

std::string QvfWriter::add_volume(const std::string& kind, const Grid& grid,
                                  const Tensor& data, const std::string& label,
                                  const std::string& component,
                                  const std::string& section_id) {
    if (kind == "volume.series")
        throw std::runtime_error("volume.series: use add_volume_series, which "
                                 "takes the axis member");
    const int64_t trailing = volume_trailing_axis(kind);
    std::vector<int64_t> spatial = data.shape;
    if (trailing) {
        if (data.shape.size() != 4 || data.shape[3] != trailing)
            throw std::runtime_error(kind + ": data must be [nx, ny, nz, " +
                                     std::to_string(trailing) + "]");
        spatial.pop_back();
    } else if (data.shape.size() != 3) {
        throw std::runtime_error(kind + ": data must be 3-D");
    }
    check_volume_dtype(kind, data);
    if (!grid.structure_ref.empty()) check_section_ref("structure_ref", grid.structure_ref);
    const Json grid_doc = grid_json(grid, spatial);   // checks the units first
    std::string sid = section_id;
    if (sid.empty()) {
        sid = kind;
        for (char& c : sid) if (c == '.') c = '_';
        sid += "_" + std::to_string(sections_.size());
    }
    check_new_id(sid);
    Json members = Json::object();
    members.set("grid", json_member("volumes/" + sid + "_grid.json", grid_doc));
    members.set("data", binary_member("volumes/" + sid + ".dat", data));
    Json section = Json::object();
    section.set("id", reserve_id(sid)).set("kind", kind).set("members", members);
    if (!label.empty()) section.set("label", label);
    if (!component.empty()) section.set("component", component);
    sections_.push_back(section);
    return sid;
}

// A volume.series `axis` payload (spec § 5.2, $defs/VolumeSeriesAxisPayload)
// checked against the n slabs of the data; `sections` are the sections added so
// far, for trajectory_ref. A Json::raw fragment is parsed. Throws otherwise.
static Json checked_series_axis(const Json& axis_in, int64_t n, const Json& sections) {
    const Json axis = axis_in.type() == Json::Type::Raw
                          ? Json::parse(axis_in.as_string()) : axis_in;
    if (!axis.is_object())
        throw std::runtime_error("volume.series axis must be a JSON object");
    for (const auto& kv : axis.items())
        if (kv.first != "kind" && kv.first != "values" && kv.first != "units" &&
            kv.first != "labels" && kv.first != "trajectory_ref" &&
            kv.first.rfind("x_", 0) != 0)
            throw std::runtime_error("volume.series axis: unknown key '" + kv.first + "'");
    const std::string kind =
        axis.has("kind") && axis.at("kind").is_string() ? axis.at("kind").as_string() : "";
    if (kind != "time" && kind != "frame" && kind != "state" && kind != "custom")
        throw std::runtime_error("volume.series axis kind must be time, frame, state "
                                 "or custom");
    if (!axis.has("values") || !axis.at("values").is_array())
        throw std::runtime_error("volume.series axis: values must be an array");
    const Json& values = axis.at("values");
    if (static_cast<int64_t>(values.size()) != n)
        throw std::runtime_error("volume.series axis has " +
                                 std::to_string(values.size()) + " values but data has " +
                                 std::to_string(n) + " slabs");
    for (const Json& v : values.elements()) {
        if (!v.is_number())
            throw std::runtime_error("volume.series axis values must be numbers");
        if (kind == "frame" && (v.as_double() < 0 || v.as_double() != std::floor(v.as_double())))
            throw std::runtime_error("volume.series frame axis values must be 0-based "
                                     "frame indices");
    }
    if (axis.has("units")) check_unit_name("volume.series axis units", axis.at("units"));
    if (axis.has("labels")) {
        const Json& labels = axis.at("labels");
        if (!labels.is_array() || static_cast<int64_t>(labels.size()) != n)
            throw std::runtime_error("volume.series axis labels must be one string per slab");
        for (const Json& l : labels.elements())
            if (!l.is_string())
                throw std::runtime_error("volume.series axis labels must be strings");
    }
    if (axis.has("trajectory_ref")) {
        const Json& ref = axis.at("trajectory_ref");
        bool found = false;
        if (ref.is_string())
            for (const Json& s : sections.elements()) {
                const std::string& k = s.at("kind").as_string();
                if (s.at("id").as_string() == ref.as_string() &&
                    (k == "trajectory" || k == "reaction.path"))
                    found = true;
            }
        if (!found)
            throw std::runtime_error("volume.series axis trajectory_ref names no "
                                     "trajectory or reaction.path section added so far");
    }
    return axis;
}

std::string QvfWriter::add_volume_series(const Grid& grid, const Tensor& data,
                                         const Json& axis, const std::string& label,
                                         const std::string& section_id) {
    if (data.shape.size() != 4)
        throw std::runtime_error("volume.series: data must be [n, nx, ny, nz]");
    check_volume_dtype("volume.series", data);
    const Json axis_doc = checked_series_axis(axis, data.shape[0], sections_);
    if (!grid.structure_ref.empty()) check_section_ref("structure_ref", grid.structure_ref);
    const Json grid_doc = grid_json(
        grid, std::vector<int64_t>(data.shape.begin() + 1, data.shape.end()));
    std::string sid = section_id;
    if (sid.empty()) sid = "volume_series_" + std::to_string(sections_.size());
    check_new_id(sid);
    Json members = Json::object();
    members.set("grid", json_member("volumes/" + sid + "_grid.json", grid_doc));
    members.set("data", binary_member("volumes/" + sid + ".dat", data));
    members.set("axis", json_member("volumes/" + sid + "_axis.json", axis_doc));
    Json section = Json::object();
    section.set("id", reserve_id(sid)).set("kind", "volume.series").set("members", members);
    if (!label.empty()) section.set("label", label);
    sections_.push_back(section);
    return sid;
}

std::string QvfWriter::add_volume_difference(const Grid& grid, const Tensor& data,
                                             const std::string& operand_a,
                                             const std::string& operand_b,
                                             const std::string& description,
                                             const std::string& label,
                                             const std::string& section_id) {
    check_new_id(section_id);
    if (data.shape.size() != 3)
        throw std::runtime_error("volume.difference: data must be 3-D");
    Json members = Json::object();
    members.set("grid", json_member("volumes/" + section_id + "_grid.json",
                                    grid_json(grid, data.shape)));
    members.set("data", binary_member("volumes/" + section_id + ".dat", data));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "volume.difference")
           .set("members", members)
           .set("operand_a", operand_a).set("operand_b", operand_b);
    if (!description.empty()) section.set("description", description);
    if (!label.empty()) section.set("label", label);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_volume_transition_density(const Grid& grid, const Tensor& data,
                                                     const std::string& excited_states_ref,
                                                     int state_index,
                                                     const std::string& label,
                                                     const std::string& section_id) {
    if (state_index < 0)
        throw std::runtime_error(
            "volume.transition_density: state_index must be >= 1 (0 omits it)");
    if (state_index > 0 && excited_states_ref.empty())
        throw std::runtime_error(
            "volume.transition_density: state_index needs excited_states_ref "
            "(QVF spec § 5.2)");
    const std::string sid =
        add_volume("volume.transition_density", grid, data, label, "", section_id);
    if (!excited_states_ref.empty())
        set_section_field(sid, "excited_states_ref", excited_states_ref);
    if (state_index > 0) set_section_field(sid, "state_index", state_index);
    return sid;
}

std::string QvfWriter::link_volume_wavefunction(const std::string& volume_id,
                                                const std::string& wavefunction_ref,
                                                int64_t mo_index,
                                                const std::string& spin) {
    const Json& volume = sections_.elements()[section_index(volume_id)];
    const std::string kind = volume.at("kind").as_string();
    if (kind != "volume.density" && kind != "volume.spin" && kind != "volume.orbital")
        throw std::runtime_error(kind + ": wavefunction_ref applies to volume.density, "
                                        "volume.spin and volume.orbital only");
    if (mo_index < -1)
        throw std::runtime_error(kind + ": mo_index must be >= 0 (-1 omits it)");
    if ((mo_index >= 0 || !spin.empty()) && kind != "volume.orbital")
        throw std::runtime_error(kind + ": mo_index and spin apply to volume.orbital only");
    if (!spin.empty() && mo_index < 0)
        throw std::runtime_error(kind + ": spin needs mo_index (QVF spec § 5.2)");
    const Json& wf = sections_.elements()[section_index(wavefunction_ref)];
    if (wf.at("kind").as_string() != "wavefunction.gto")
        throw std::runtime_error(kind + ": wavefunction_ref '" + wavefunction_ref +
                                 "' names a '" + wf.at("kind").as_string() +
                                 "' section, not a 'wavefunction.gto' section");
    if (mo_index >= 0) {
        const Json& members = wf.at("members");
        std::string role = "mo_coefficients";
        if (members.has("mo_coefficients")) {
            if (!spin.empty())
                throw std::runtime_error(kind + ": wavefunction '" + wavefunction_ref +
                                         "' has one set of orbitals; spin must be "
                                         "omitted (QVF spec § 5.2)");
        } else {
            if (spin != "alpha" && spin != "beta")
                throw std::runtime_error(kind + ": wavefunction '" + wavefunction_ref +
                                         "' is unrestricted; spin must be 'alpha' or "
                                         "'beta' (QVF spec § 5.2)");
            role += "_" + spin;
        }
        if (!members.has(role) || !members.at(role).has("shape") ||
            members.at(role).at("shape").size() == 0)
            throw std::runtime_error(kind + ": wavefunction '" + wavefunction_ref +
                                     "' has no " + role + " member with a shape");
        const int64_t rows = members.at(role).at("shape").elements()[0].as_int();
        if (mo_index >= rows)
            throw std::runtime_error(kind + ": mo_index " + std::to_string(mo_index) +
                                     " is outside the " + std::to_string(rows) +
                                     " rows of " + role + " of '" + wavefunction_ref + "'");
    }
    set_section_field(volume_id, "wavefunction_ref", wavefunction_ref);
    if (mo_index >= 0) set_section_field(volume_id, "mo_index", mo_index);
    if (!spin.empty()) set_section_field(volume_id, "spin", spin);
    return volume_id;
}

std::string QvfWriter::add_basis_ao(const Json& ao_metadata, const Grid& grid,
                                    const Tensor& data, const std::string& label,
                                    const std::string& section_id) {
    if (data.shape.size() != 3)
        throw std::runtime_error("basis.ao: data must be 3-D");
    std::string sid = section_id.empty()
                          ? "ao_" + std::to_string(sections_.size()) : section_id;
    check_new_id(sid);
    Json members = Json::object();
    members.set("grid", json_member("basis_ao/" + sid + "_grid.json",
                                    grid_json(grid, data.shape)));
    members.set("data", binary_member("basis_ao/" + sid + ".dat", data));
    Json section = Json::object();
    section.set("id", reserve_id(sid)).set("kind", "basis.ao").set("members", members)
           .set("ao_metadata", ao_metadata);
    if (!label.empty()) section.set("label", label);
    sections_.push_back(section);
    return sid;
}

// -- bands & DOS ------------------------------------------------------------

std::string QvfWriter::add_bands(const Json& kpath, const Tensor& eigenvalues,
                                 const Tensor* projections,
                                 const std::string& section_id) {
    return add_bands(kpath, eigenvalues, projections, nullptr, section_id);
}

std::string QvfWriter::add_bands(const Json& kpath, const Tensor& eigenvalues,
                                 const Tensor* projections, const Tensor* occupations,
                                 const std::string& section_id) {
    check_new_id(section_id);
    const auto& e = eigenvalues.shape;
    if (e.size() != 3)
        throw std::runtime_error("bands: eigenvalues must be [n_spin, n_k, n_bands]");
    if (projections) {
        const auto& p = projections->shape;
        if (projections->dtype != DType::Float64)
            throw std::runtime_error("bands: projections must be float64");
        if (p.size() == 4) {
            if (p[0] != e[0] || p[1] != e[1] || p[2] != e[2])
                throw std::runtime_error(
                    "bands: projections [n_spin, n_k, n_bands, n_channels] must "
                    "repeat the shape of eigenvalues");
        } else if (p.size() == 3) {
            if (p[0] != e[1] || p[1] != e[2])
                throw std::runtime_error(
                    "bands: projections [n_k, n_bands, n_channels] must repeat "
                    "n_k and n_bands of eigenvalues");
        } else {
            throw std::runtime_error(
                "bands: projections must be [n_spin, n_k, n_bands, n_channels] "
                "or [n_k, n_bands, n_channels]");
        }
    }
    if (occupations) {
        if (occupations->dtype != DType::Float64)
            throw std::runtime_error("bands: occupations must be float64");
        if (occupations->shape != e)
            throw std::runtime_error(
                "bands: occupations must have the shape of eigenvalues "
                "[n_spin, n_k, n_bands]");
    }
    Json members = Json::object();
    members.set("kpath", json_member(section_id + "/kpath.json", kpath));
    members.set("eigenvalues",
                binary_member(section_id + "/eigenvalues.bin", eigenvalues));
    if (projections)
        members.set("projections",
                    binary_member(section_id + "/projections.bin", *projections));
    if (occupations)
        members.set("occupations",
                    binary_member(section_id + "/occupations.bin", *occupations));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "bands")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_dos_total(const Tensor& energies, const Tensor& dos,
                                     const Json& section_meta,
                                     const std::string& section_id) {
    check_new_id(section_id);
    Json members = Json::object();
    members.set("energies",
                binary_member("dos/" + section_id + "_energies.bin", energies));
    members.set("dos", binary_member("dos/" + section_id + ".bin", dos));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "dos.total")
           .set("members", members);
    merge_peers(section, section_meta);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_dos_projected(const Tensor& energies,
                                         const Tensor& projections,
                                         const Json& section_meta,
                                         const std::string& section_id) {
    check_new_id(section_id);
    const auto& p = projections.shape;
    if (p.size() != 2 && p.size() != 3)
        throw std::runtime_error(
            "dos.projected: projections must be [n_spin, n_channels, n_points] "
            "or [n_channels, n_points]");
    if (p.size() == 3 && section_meta.is_object() && section_meta.contains("n_spin")) {
        const Json& ns = section_meta.at("n_spin");
        if (ns.type() == Json::Type::Int && ns.as_int() != p[0])
            throw std::runtime_error(
                "dos.projected: the spin axis of projections must equal n_spin");
    }
    Json members = Json::object();
    members.set("energies",
                binary_member("dos/" + section_id + "_energies.bin", energies));
    members.set("projections",
                binary_member("dos/" + section_id + "_projections.bin", projections));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "dos.projected")
           .set("members", members);
    merge_peers(section, section_meta);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_dos_coop_cohp(const std::string& kind,
                                         const Tensor& energies,
                                         const Tensor& projections,
                                         const Tensor& integrated, const Json& meta,
                                         const std::string& section_id) {
    if (kind != "dos.coop" && kind != "dos.cohp")
        throw std::runtime_error("kind must be dos.coop or dos.cohp");
    std::string sid = section_id;
    if (sid.empty()) { sid = kind; for (char& c : sid) if (c == '.') c = '_'; }
    check_new_id(sid);
    Json members = Json::object();
    members.set("energies", binary_member("dos/" + sid + "_energies.bin", energies));
    members.set("projections",
                binary_member("dos/" + sid + "_projections.bin", projections));
    members.set("integrated",
                binary_member("dos/" + sid + "_integrated.bin", integrated));
    members.set("meta", json_member("dos/" + sid + "_meta.json", meta));
    Json section = Json::object();
    section.set("id", reserve_id(sid)).set("kind", kind).set("members", members);
    sections_.push_back(section);
    return sid;
}

// -- trajectories / reactions / scans / vibrations --------------------------

// Throws unless `t` has `dtype` and itemsize * prod(shape) bytes.
static void require_tensor(const Tensor& t, DType dtype, const std::string& what) {
    if (t.dtype != dtype)
        throw std::runtime_error(what + " must be " + dtype_name(dtype));
    size_t n = 1;
    for (int64_t d : t.shape) {
        if (d < 0) throw std::runtime_error(what + " has a negative extent");
        n *= static_cast<size_t>(d);
    }
    if (t.bytes.size() != dtype_itemsize(t.dtype) * n)
        throw std::runtime_error(what + " byte length != itemsize * prod(shape)");
}

static std::string shape_str(const std::vector<int64_t>& shape) {
    std::string s = "[";
    for (size_t k = 0; k < shape.size(); ++k)
        s += (k ? ", " : "") + std::to_string(shape[k]);
    return s + "]";
}

// The little-endian float64 values of a float64 tensor.
static std::vector<double> f64_values(const Tensor& t) {
    std::vector<double> out(t.bytes.size() / 8);
    for (size_t k = 0; k < out.size(); ++k) {
        uint64_t bits = 0;
        for (int b = 7; b >= 0; --b) bits = (bits << 8) | t.bytes[8 * k + b];
        std::memcpy(&out[k], &bits, 8);
    }
    return out;
}

// Largest difference, in hartree, between a binary `energies` member and the
// metadata `energies` array of one section (spec § 5.5).
constexpr double kTrajectoryEnergyToleranceEh = 1e-10;

std::string QvfWriter::add_frames_section(const std::string& kind, const Json& metadata,
                                          const Tensor& coords, const Tensor* lattice,
                                          const TrajectoryExtras& extras,
                                          const std::string& section_id) {
    check_new_id(section_id);
    const bool trajectory = kind == "trajectory";
    if (coords.shape.size() != 3 || coords.shape[2] != 3)
        throw std::runtime_error(kind + " coords must be [n_frames, n_atoms, 3]");
    if (trajectory && coords.dtype == DType::Float32)
        require_tensor(coords, DType::Float32, kind + " coords");
    else
        require_tensor(coords, DType::Float64,
                       kind + (trajectory ? " coords (float64 or float32)" : " coords"));
    const int64_t n_frames = coords.shape[0];
    const int64_t n_atoms = coords.shape[1];

    // The metadata, parsed when it is a raw fragment, for the checks below.
    Json meta = metadata;
    if (metadata.type() == Json::Type::Raw) {
        try {
            meta = Json::parse(metadata.as_string());
        } catch (const std::exception& e) {
            throw std::runtime_error(kind + " metadata is not valid JSON: " + e.what());
        }
    }
    if (!meta.is_object()) throw std::runtime_error(kind + " metadata must be a JSON object");
    if (trajectory)
        check_payload_units(kind, meta,
                            {{"energies", Unit::Hartree}, {"coords", Unit::Angstrom},
                             {"times", Unit::Femtosecond}, {"timestep_fs", Unit::Femtosecond},
                             {"lattice", Unit::Angstrom}, {"cell_volumes", Unit::Angstrom3},
                             {"temperatures", Unit::Kelvin}});
    else
        check_payload_units(kind, meta,
                            {{"energies", Unit::Hartree}, {"coords", Unit::Angstrom},
                             {"times", Unit::Femtosecond}, {"timestep_fs", Unit::Femtosecond},
                             {"lattice", Unit::Bohr}});

    // precision names the dtype of coords; float32 coords add it when absent.
    const std::string dtype = dtype_name(coords.dtype);
    bool add_precision = false;
    if (meta.contains("precision")) {
        const Json& p = meta.at("precision");
        if (!p.is_string() || (p.as_string() != "float64" && p.as_string() != "float32"))
            throw std::runtime_error(kind + " metadata precision must be \"float64\" or "
                                            "\"float32\"");
        if (p.as_string() != dtype)
            throw std::runtime_error(kind + " metadata precision \"" + p.as_string() +
                                     "\" contradicts " + dtype + " coords");
    } else if (coords.dtype == DType::Float32) {
        add_precision = true;
    }

    if (lattice) {
        const auto& s = lattice->shape;
        const bool fixed = (s.size() == 2 && s[0] == 3 && s[1] == 3);
        const bool per_frame = (s.size() == 3 && s[1] == 3 && s[2] == 3);
        if (!fixed && !per_frame)
            throw std::runtime_error(kind + " lattice must be [3,3] or [n_frames,3,3]");
        if (per_frame && s[0] != n_frames)
            throw std::runtime_error(kind + " per-frame lattice frame count != coords frames");
        require_tensor(*lattice, DType::Float64, kind + " lattice");
    }
    if (trajectory) {
        const Json* pbc = meta.contains("pbc") ? &meta.at("pbc") : nullptr;
        bool periodic = false;
        if (pbc) {
            if (!pbc->is_array() || pbc->size() != 3)
                throw std::runtime_error("trajectory metadata pbc must be three booleans");
            for (const Json& x : pbc->elements()) {
                if (!x.is_bool())
                    throw std::runtime_error("trajectory metadata pbc must be three booleans");
                periodic = periodic || x.as_bool();
            }
        }
        if (lattice && !pbc)
            throw std::runtime_error("trajectory lattice needs metadata \"pbc\": pbc is "
                                     "REQUIRED when lattice is present (spec § 5.5)");
        if (!lattice && periodic)
            throw std::runtime_error("trajectory metadata pbc marks a periodic axis, so "
                                     "a lattice is required (spec § 5.5)");
    }

    struct Frame {
        const char* role;
        const Tensor* t;
        bool vectors;
    };
    const Frame frames[] = {
        {"velocities", extras.velocities, true},  {"gradient", extras.gradient, true},
        {"times", extras.times, false},           {"energies", extras.energies, false},
        {"cell_volumes", extras.cell_volumes, false},
        {"temperatures", extras.temperatures, false},
    };
    for (const Frame& f : frames) {
        if (!f.t) continue;
        const std::string what = kind + " " + f.role;
        const std::vector<int64_t> want =
            f.vectors ? std::vector<int64_t>{n_frames, n_atoms, 3}
                      : std::vector<int64_t>{n_frames};
        if (f.t->shape != want)
            throw std::runtime_error(what + " must be " + shape_str(want) +
                                     " to match coords, got " + shape_str(f.t->shape));
        require_tensor(*f.t, DType::Float64, what);
    }
    if (extras.energies && meta.contains("energies")) {
        const Json& listed = meta.at("energies");
        const std::vector<double> binary = f64_values(*extras.energies);
        bool agree = listed.is_array() && listed.size() == binary.size();
        for (size_t k = 0; agree && k < binary.size(); ++k)
            agree = listed.at(k).is_number() &&
                    std::fabs(listed.at(k).as_double() - binary[k]) <=
                        kTrajectoryEnergyToleranceEh;
        if (!agree)
            throw std::runtime_error(kind + " energies and metadata energies must agree "
                                            "within 1e-10 hartree (spec § 5.5)");
    }

    const std::string sid = reserve_id(section_id);
    Json members = Json::object();
    if (add_precision) {
        meta.set("precision", "float32");
        members.set("metadata", json_member(sid + "/metadata.json", meta));
    } else {
        members.set("metadata", json_member(sid + "/metadata.json", metadata));
    }
    members.set("coords", binary_member(sid + "/coords.bin", coords));
    if (lattice) members.set("lattice", binary_member(sid + "/lattice.bin", *lattice));
    for (const Frame& f : frames)
        if (f.t)
            members.set(f.role, binary_member(sid + "/" + f.role + ".bin", *f.t));
    Json section = Json::object();
    section.set("id", sid).set("kind", kind).set("members", members);
    sections_.push_back(section);
    return sid;
}

std::string QvfWriter::add_trajectory(const Json& metadata, const Tensor& coords,
                                      const std::string& section_id) {
    return add_trajectory(metadata, coords, TrajectoryExtras(), section_id);
}

std::string QvfWriter::add_trajectory(const Json& metadata, const Tensor& coords,
                                      const TrajectoryExtras& extras,
                                      const std::string& section_id) {
    return add_frames_section("trajectory", metadata, coords, extras.lattice, extras,
                              section_id);
}

// The [3, 3] rows-in-Å matrices of `lattices`, flattened frame by frame.
static Tensor lattice_rows_tensor(const std::vector<Lattice>& lattices, bool per_frame) {
    std::vector<double> flat;
    flat.reserve(9 * lattices.size());
    for (const Lattice& cell : lattices)
        for (const auto& row : cell.rows_angstrom())
            flat.insert(flat.end(), row.begin(), row.end());
    if (!per_frame) return tensor_f64(flat, {3, 3});
    return tensor_f64(flat, {static_cast<int64_t>(lattices.size()), 3, 3});
}

std::string QvfWriter::add_trajectory(const Json& metadata, const Tensor& coords,
                                      const Lattice& lattice,
                                      const TrajectoryExtras& extras,
                                      const std::string& section_id) {
    if (extras.lattice)
        throw std::runtime_error("trajectory: pass the lattice once, as a Lattice or "
                                 "as extras.lattice");
    const Tensor t = lattice_rows_tensor({lattice}, false);
    return add_frames_section("trajectory", metadata, coords, &t, extras, section_id);
}

std::string QvfWriter::add_trajectory(const Json& metadata, const Tensor& coords,
                                      const std::vector<Lattice>& lattices,
                                      const TrajectoryExtras& extras,
                                      const std::string& section_id) {
    if (extras.lattice)
        throw std::runtime_error("trajectory: pass the lattice once, as Lattices or "
                                 "as extras.lattice");
    const Tensor t = lattice_rows_tensor(lattices, true);
    return add_frames_section("trajectory", metadata, coords, &t, extras, section_id);
}

std::string QvfWriter::add_reaction_path(const Json& metadata, const Tensor& coords,
                                         const Tensor* lattice,
                                         const std::string& section_id) {
    return add_reaction_path(metadata, coords, lattice, TrajectoryExtras(), section_id);
}

std::string QvfWriter::add_reaction_path(const Json& metadata, const Tensor& coords,
                                         const Tensor* lattice,
                                         const TrajectoryExtras& extras,
                                         const std::string& section_id) {
    if (extras.lattice || extras.cell_volumes || extras.temperatures)
        throw std::runtime_error(
            "reaction.path carries velocities, gradient, times and energies; its cell "
            "is the lattice argument (columns, bohr), and it has no cell_volumes or "
            "temperatures");
    return add_frames_section("reaction.path", metadata, coords, lattice, extras,
                              section_id);
}

std::string QvfWriter::add_reaction_path(const Json& metadata, const Tensor& coords,
                                         const Lattice* lattice,
                                         const std::string& section_id) {
    if (!lattice)
        return add_reaction_path(metadata, coords,
                                 static_cast<const Tensor*>(nullptr), section_id);
    const Lattice::Mat3 cols = lattice->columns_bohr();
    std::vector<double> flat;
    flat.reserve(9);
    for (const auto& row : cols) flat.insert(flat.end(), row.begin(), row.end());
    const Tensor t = tensor_f64(flat, {3, 3});
    return add_reaction_path(metadata, coords, &t, section_id);
}

std::string QvfWriter::add_reaction_path(const Json& metadata, const Tensor& coords,
                                         std::nullptr_t,
                                         const std::string& section_id) {
    return add_reaction_path(metadata, coords,
                             static_cast<const Tensor*>(nullptr), section_id);
}

std::string QvfWriter::add_reaction_waypoints(const std::string& trajectory_ref,
                                              const Json& waypoints,
                                              const std::string& section_id) {
    check_new_id(section_id);
    if (!ids_.count(trajectory_ref))
        throw std::runtime_error("trajectory_ref not yet added: " + trajectory_ref);
    Json members = Json::object();
    members.set("waypoints", json_member(section_id + "/waypoints.json", waypoints));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "reaction.waypoints")
           .set("trajectory_ref", trajectory_ref).set("members", members);
    sections_.push_back(section);
    return section_id;
}

// -- reaction.scheme (spec § 5.14) ------------------------------------------

namespace {

// `key` is present in the object `obj` and not null (null means absent).
bool given_field(const Json& obj, const char* key) {
    return obj.has(key) && !obj.at(key).is_null();
}

// A strict payload object: every key is in `allowed` or begins with x_.
void check_strict_keys(const Json& obj, std::initializer_list<const char*> allowed,
                       const std::string& what) {
    if (!obj.is_object()) field_error(what, "must be a JSON object");
    std::string names;
    for (const char* a : allowed) names += (names.empty() ? "" : ", ") + std::string(a);
    for (const auto& kv : obj.items()) {
        bool known = is_vendor_key(kv.first);
        for (const char* a : allowed) known = known || kv.first == a;
        if (!known)
            throw std::runtime_error(what + ": unknown key '" + kv.first + "'; allowed are " +
                                     names + " and keys beginning x_");
    }
}

// One role of a reaction.scheme step as participant objects: each item is a
// structure section id or {structure_ref, coefficient?} (further keys x_).
Json scheme_participants(const Json* items, const std::string& what, bool required) {
    if (!items || items->is_null()) {
        if (required) field_error(what, "is required");
        return Json::array();
    }
    const Json list = resolved(*items, what);
    if (!list.is_array()) field_error(what, "must be an array of participants");
    Json out = Json::array();
    std::map<std::string, size_t> uses;   // sorted, for the message
    size_t k = 0;
    for (const Json& raw : list.elements()) {
        const std::string where = what + "[" + std::to_string(k++) + "]";
        Json entry = Json::object();
        if (raw.is_string()) {
            entry.set("structure_ref", raw);
        } else {
            entry = resolved(raw, where);
            if (!entry.is_object())
                throw std::runtime_error(where + " must be a structure section id or an "
                                                 "object, got " + entry.dump(0));
            check_strict_keys(entry, {"structure_ref", "coefficient"}, where);
        }
        if (!entry.has("structure_ref") || !entry.at("structure_ref").is_string() ||
            entry.at("structure_ref").as_string().empty())
            field_error(where + ": structure_ref", "must be a non-empty string");
        if (entry.has("coefficient"))
            checked_number(entry.at("coefficient"), where + " coefficient", true);
        ++uses[entry.at("structure_ref").as_string()];
        out.push_back(entry);
    }
    if (required && out.size() == 0) field_error(what, "needs at least one participant");
    std::string repeated;
    for (const auto& u : uses)
        if (u.second > 1) repeated += (repeated.empty() ? "'" : ", '") + u.first + "'";
    if (!repeated.empty())
        throw std::runtime_error(what + " names " + repeated +
                                 " more than once; a coefficient says how many "
                                 "(QVF spec § 5.14)");
    return out;
}

}  // namespace

std::string QvfWriter::add_reaction_scheme(const Json& steps_in, const std::string& label,
                                           const std::string& section_id) {
    const Json steps = resolved(steps_in, "reaction.scheme steps");
    if (!steps.is_array())
        throw std::runtime_error("reaction.scheme steps must be a JSON array of step "
                                 "objects (QVF spec § 5.14)");
    Json out = Json::array();
    std::set<std::string> seen;
    size_t k = 0;
    for (const Json& raw : steps.elements()) {
        const std::string what = "reaction.scheme step " + std::to_string(k++);
        const Json step = resolved(raw, what);
        check_strict_keys(step,
                          {"id", "label", "reactants", "products", "reagents", "arrow",
                           "conditions", "yield_percent", "yield_text", "energies"},
                          what);
        const Json* id = step.has("id") ? &step.at("id") : nullptr;
        if (!id || !id->is_string() || !is_valid_id(id->as_string()))
            throw std::runtime_error(what + ": id must match [A-Za-z0-9_.-]+, got " +
                                     (id ? id->dump(0) : std::string("nothing")));
        if (!seen.insert(id->as_string()).second)
            throw std::runtime_error(what + ": id '" + id->as_string() +
                                     "' is used by another step");
        Json entry = Json::object();
        entry.set("id", *id);
        if (given_field(step, "label"))
            entry.set("label", checked_text(step.at("label"), what + " label"));
        for (const char* role : {"reactants", "products", "reagents"}) {
            const bool reagents = std::strcmp(role, "reagents") == 0;
            if (reagents && !given_field(step, role)) continue;
            const Json parts = scheme_participants(step.has(role) ? &step.at(role) : nullptr,
                                                   what + " " + role, !reagents);
            for (const Json& p : parts.elements())
                require_section_kind(what + " " + role + " structure_ref",
                                     p.at("structure_ref").as_string(), "structure");
            entry.set(role, parts);
        }
        if (given_field(step, "arrow")) {
            const Json& arrow = step.at("arrow");
            const std::string a = arrow.is_string() ? arrow.as_string() : std::string();
            if (a != "forward" && a != "equilibrium" && a != "retrosynthetic" &&
                a != "resonance")
                throw std::runtime_error(what + " arrow must be one of forward, equilibrium, "
                                                "retrosynthetic, resonance, got " +
                                         arrow.dump(0));
            entry.set("arrow", arrow);
        }
        if (given_field(step, "conditions")) {
            const Json lines = resolved(step.at("conditions"), what + " conditions");
            bool ok = lines.is_array();
            if (ok)
                for (const Json& line : lines.elements()) ok = ok && line.is_string();
            if (!ok)
                field_error(what + ": conditions", "must be an array of strings, one per line");
            entry.set("conditions", lines);
        }
        if (given_field(step, "yield_percent")) {
            const Json& y = step.at("yield_percent");
            const double v = checked_number(y, what + " yield_percent", false);
            if (v < 0.0 || v > 100.0)
                throw std::runtime_error(what + ": yield_percent must be from 0 to 100, got " +
                                         y.dump(0));
            entry.set("yield_percent", y);
        }
        if (given_field(step, "yield_text"))
            entry.set("yield_text", checked_text(step.at("yield_text"), what + " yield_text"));
        // Reaction energies (spec § 5.14) are written as given; the validator
        // checks their terms and sums against this step.
        if (given_field(step, "energies")) {
            const Json energies = resolved(step.at("energies"), what + " energies");
            if (!energies.is_object())
                field_error(what + ": energies", "must be an object ($defs/ReactionStepEnergies)");
            for (const char* key : {"version", "deltas", "terms"})
                if (!energies.has(key))
                    field_error(what + ": energies", std::string("needs ") + key +
                                                        " (QVF spec § 5.14)");
            entry.set("energies", energies);
        }
        for (const auto& kv : step.items())
            if (is_vendor_key(kv.first)) entry.set(kv.first, kv.second);
        out.push_back(entry);
    }
    if (out.size() == 0) throw std::runtime_error("reaction.scheme needs at least one step");
    check_new_id(section_id);
    Json payload = Json::object();
    payload.set("steps", out);
    Json members = Json::object();
    members.set("scheme", json_member(section_id + "/scheme.json", payload));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "reaction.scheme");
    if (!label.empty()) section.set("label", label);
    section.set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_scan_surface(const Json& metadata, const Tensor& axis_a,
                                        const Tensor& axis_b, const Tensor& energies,
                                        const Tensor* geometries,
                                        const std::string& section_id) {
    check_new_id(section_id);
    check_payload_units("scan.surface", metadata, {});
    Json members = Json::object();
    members.set("metadata", json_member(section_id + "/metadata.json", metadata));
    members.set("axis_a", binary_member(section_id + "/axis_a.bin", axis_a));
    members.set("axis_b", binary_member(section_id + "/axis_b.bin", axis_b));
    members.set("energies", binary_member(section_id + "/energies.bin", energies));
    if (geometries)
        members.set("geometries",
                    binary_member(section_id + "/geometries.bin", *geometries));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "scan.surface")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_vibrations(const Json& metadata,
                                      const Tensor& displacements,
                                      const std::string& section_id,
                                      const std::string& hessian_ref) {
    check_new_id(section_id);
    if (displacements.shape.size() != 3 || displacements.shape[2] != 3)
        throw std::runtime_error("displacements must be [n_modes, n_atoms, 3]");
    if (!hessian_ref.empty()) {  // spec § 5.5: names an existing hessian section
        bool found = false;
        for (const Json& s : sections_.elements())
            if (s.at("id").as_string() == hessian_ref &&
                s.at("kind").as_string() == "hessian")
                found = true;
        if (!found)
            throw std::runtime_error("vibrations hessian_ref names no hessian section "
                                     "added so far: " + hessian_ref);
    }
    Json members = Json::object();
    members.set("metadata", json_member(section_id + "/metadata.json", metadata));
    members.set("displacements",
                binary_member(section_id + "/displacements.bin", displacements));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "vibrations");
    if (!hessian_ref.empty()) section.set("hessian_ref", hessian_ref);
    section.set("members", members);
    sections_.push_back(section);
    return section_id;
}

// -- energy derivatives (spec § 5.10, § 5.11) -------------------------------

// Throws unless `t` is float64 with itemsize * prod(shape) bytes, so a refused
// call is refused before anything is written.
static void require_f64_tensor(const Tensor& t, const char* what) {
    if (t.dtype != DType::Float64)
        throw std::runtime_error(std::string(what) + " must be float64");
    size_t n = 1;
    for (int64_t d : t.shape) {
        if (d < 0) throw std::runtime_error(std::string(what) + " has a negative extent");
        n *= static_cast<size_t>(d);
    }
    if (t.bytes.size() != 8 * n)
        throw std::runtime_error(std::string(what) +
                                 " byte length != 8 * prod(shape)");
}

static void require_meta(const Json& meta, const char* what) {
    if (!meta.is_null() && meta.type() != Json::Type::Object &&
        meta.type() != Json::Type::Raw)
        throw std::runtime_error(std::string(what) + " meta must be a JSON object");
}

std::string QvfWriter::add_forces(const Tensor& gradient, const Tensor* stress,
                                  const Json& meta, const std::string& section_id) {
    check_new_id(section_id);
    if (gradient.shape.size() != 2 || gradient.shape[1] != 3)
        throw std::runtime_error("forces gradient must be [n_atoms, 3]");
    require_f64_tensor(gradient, "forces gradient");
    if (stress) {
        if (stress->shape.size() != 2 || stress->shape[0] != 3 || stress->shape[1] != 3)
            throw std::runtime_error("forces stress must be [3, 3]");
        require_f64_tensor(*stress, "forces stress");
    }
    require_meta(meta, "forces");
    const std::string sid = reserve_id(section_id);
    Json members = Json::object();
    members.set("gradient", binary_member(sid + "/gradient.bin", gradient));
    if (stress) members.set("stress", binary_member(sid + "/stress.bin", *stress));
    if (!meta.is_null()) members.set("meta", json_member(sid + "/meta.json", meta));
    Json section = Json::object();
    section.set("id", sid).set("kind", "forces").set("members", members);
    sections_.push_back(section);
    return sid;
}

std::string QvfWriter::add_hessian(const Tensor& matrix, const Json& meta,
                                   const std::string& section_id) {
    check_new_id(section_id);
    const auto& s = matrix.shape;
    if (s.size() != 2 || s[0] != s[1])
        throw std::runtime_error("hessian matrix must be square [3N, 3N]");
    if (s[0] % 3 != 0)
        throw std::runtime_error("hessian matrix extent is not divisible by 3");
    require_f64_tensor(matrix, "hessian matrix");
    require_meta(meta, "hessian");
    if (meta.type() == Json::Type::Object) {
        const bool weighted = meta.has("mass_weighted") &&
                              meta.at("mass_weighted").is_bool() &&
                              meta.at("mass_weighted").as_bool();
        if (weighted && !meta.has("atom_masses_amu"))
            throw std::runtime_error(
                "hessian meta: a mass-weighted matrix needs atom_masses_amu");
        if (meta.has("atom_masses_amu")) {
            const Json& masses = meta.at("atom_masses_amu");
            if (!masses.is_array() || static_cast<int64_t>(masses.size()) != s[0] / 3)
                throw std::runtime_error(
                    "hessian meta atom_masses_amu must hold one mass per atom "
                    "(N = matrix extent / 3)");
        }
    }
    const std::string sid = reserve_id(section_id);
    Json members = Json::object();
    members.set("matrix", binary_member(sid + "/matrix.bin", matrix));
    if (!meta.is_null()) members.set("meta", json_member(sid + "/meta.json", meta));
    Json section = Json::object();
    section.set("id", sid).set("kind", "hessian").set("members", members);
    sections_.push_back(section);
    return sid;
}

// -- excited states (spec § 5.12) ---------------------------------------------

// The fields of a state ($defs/ExcitedState) and of a contribution
// ($defs/ExcitedStateContribution); other keys must begin with x_.
static const char* const kStateFields[] = {
    "index", "excitation_energy", "total_energy", "multiplicity", "symmetry",
    "oscillator_strength", "rotatory_strength_cgs", "transition_dipole",
    "transition_magnetic_dipole", "contributions", "nto_ref", "density_ref"};
static const char* const kContributionFields[] = {"from_mo", "to_mo", "spin", "weight"};

template <size_t N>
static bool is_listed(const std::string& key, const char* const (&list)[N]) {
    for (const char* k : list)
        if (key == k) return true;
    return key.compare(0, 2, "x_") == 0;
}

static bool is_json_int(const Json& j) { return j.type() == Json::Type::Int; }

// A {value, units} energy whose unit is a unit name; a number is written in
// `default_unit`.
static Json state_energy(const Json& v, const char* default_unit, const std::string& what) {
    if (is_finite_number(v)) {
        Json e = Json::object();
        e.set("value", v.as_double()).set("units", default_unit);
        return e;
    }
    if (!v.is_object() || v.size() != 2 || !v.has("value") || !v.has("units") ||
        !is_finite_number(v.at("value")))
        throw std::runtime_error(what + " must be a number or a {value, units} object");
    check_unit_name(what + " units", v.at("units"));
    return v;
}

// Orbital rows of a wavefunction.gto section by channel: "" for
// mo_coefficients, "alpha" and "beta" for the unrestricted pair.
static std::map<std::string, int64_t> section_mo_rows(const Json& section) {
    std::map<std::string, int64_t> rows;
    if (!section.has("members")) return rows;
    const Json& members = section.at("members");
    static const struct {
        const char* role;
        const char* key;
    } kRoles[] = {{"mo_coefficients", ""},
                  {"mo_coefficients_alpha", "alpha"},
                  {"mo_coefficients_beta", "beta"}};
    for (const auto& r : kRoles)
        if (members.has(r.role) && members.at(r.role).has("shape") &&
            members.at(r.role).at("shape").size() > 0)
            rows[r.key] = members.at(r.role).at("shape").at(size_t{0}).as_int();
    return rows;
}

std::string QvfWriter::add_excited_states(const Json& states_in, const Json& meta_in,
                                          const std::string& wavefunction_ref,
                                          const std::string& section_id) {
    check_new_id(section_id);
    const Json states = states_in.type() == Json::Type::Raw
                            ? Json::parse(states_in.as_string())
                            : states_in;
    const Json meta = meta_in.type() == Json::Type::Raw ? Json::parse(meta_in.as_string())
                                                        : meta_in;
    if (!states.is_array())
        throw std::runtime_error("excited_states: states must be a JSON array (QVF spec § 5.12)");
    if (!meta.is_null() && !meta.is_object())
        throw std::runtime_error("excited_states: meta must be a JSON object");

    auto kind_of = [this](const std::string& id) -> std::string {
        for (const Json& s : sections_.elements())
            if (s.at("id").as_string() == id) return s.at("kind").as_string();
        return "";
    };
    auto need = [&kind_of](const std::string& what, const Json& ref,
                           std::initializer_list<const char*> kinds) {
        const std::string found = ref.is_string() ? kind_of(ref.as_string()) : "";
        std::string names;
        for (const char* k : kinds) {
            if (found == k) return;
            names += names.empty() ? std::string(k) : std::string(" or ") + k;
        }
        throw std::runtime_error("excited_states: " + what + " " + ref.dump(0) +
                                 " must name an added " + names +
                                 " section (QVF spec § 5.12)");
    };

    bool have_rows = false;
    std::map<std::string, int64_t> rows;
    if (!wavefunction_ref.empty()) {
        need("wavefunction_ref", Json(wavefunction_ref), {"wavefunction.gto"});
        rows = section_mo_rows(sections_.at(section_index(wavefunction_ref)));
        have_rows = true;
    }

    Json out_states = Json::array();
    int64_t prev = 0;
    for (size_t k = 0; k < states.size(); ++k) {
        const Json& st = states.at(k);
        const std::string what = "excited_states: states[" + std::to_string(k) + "]";
        if (!st.is_object()) throw std::runtime_error(what + " must be a JSON object");
        if (!st.has("index") || !st.has("excitation_energy"))
            throw std::runtime_error(what + " needs index and excitation_energy "
                                            "(QVF spec § 5.12)");
        Json out = Json::object();
        for (const auto& kv : st.items()) {
            const std::string& key = kv.first;
            const Json& v = kv.second;
            const std::string field = what + "." + key;
            if (!is_listed(key, kStateFields))
                throw std::runtime_error(what + ": unknown field '" + key +
                                         "' (QVF spec § 5.12)");
            if (key == "index") {
                if (!is_json_int(v) || v.as_int() < 1)
                    throw std::runtime_error(field + " must be an integer >= 1");
                if (k > 0 && v.as_int() <= prev)
                    throw std::runtime_error(field + " " + std::to_string(v.as_int()) +
                                             " follows index " + std::to_string(prev) +
                                             "; indices are unique and ascending "
                                             "(QVF spec § 5.12)");
                prev = v.as_int();
                out.set(key, v);
            } else if (key == "excitation_energy") {
                out.set(key, state_energy(v, "ev", field));
            } else if (key == "total_energy") {
                out.set(key, state_energy(v, "hartree", field));
            } else if (key == "multiplicity") {
                if (!is_json_int(v) || v.as_int() < 1)
                    throw std::runtime_error(field + " must be an integer >= 1");
                out.set(key, v);
            } else if (key == "symmetry") {
                if (!v.is_string()) throw std::runtime_error(field + " must be a string");
                out.set(key, v);
            } else if (key == "oscillator_strength" || key == "rotatory_strength_cgs") {
                if (!is_finite_number(v) ||
                    (key == "oscillator_strength" && v.as_double() < 0))
                    throw std::runtime_error(field + " must be a finite number" +
                                             (key == "oscillator_strength" ? " >= 0" : ""));
                out.set(key, v);
            } else if (key == "transition_dipole" || key == "transition_magnetic_dipole") {
                bool ok = v.is_array() && v.size() == 3;
                for (size_t c = 0; ok && c < 3; ++c) ok = is_finite_number(v.at(c));
                if (!ok) throw std::runtime_error(field + " must be three numbers");
                out.set(key, v);
            } else if (key == "contributions") {
                if (!v.is_array()) throw std::runtime_error(field + " must be an array");
                for (size_t c = 0; c < v.size(); ++c) {
                    const Json& con = v.at(c);
                    const std::string cw = field + "[" + std::to_string(c) + "]";
                    if (!con.is_object() || !con.has("from_mo") || !con.has("to_mo") ||
                        !con.has("weight"))
                        throw std::runtime_error(cw + " needs from_mo, to_mo and weight "
                                                      "(QVF spec § 5.12)");
                    for (const auto& ckv : con.items())
                        if (!is_listed(ckv.first, kContributionFields))
                            throw std::runtime_error(cw + ": unknown field '" + ckv.first +
                                                     "' (QVF spec § 5.12)");
                    if (!is_finite_number(con.at("weight")))
                        throw std::runtime_error(cw + ".weight must be a number");
                    std::string spin;
                    if (con.has("spin")) {
                        const Json& sj = con.at("spin");
                        if (!sj.is_string() ||
                            (sj.as_string() != "alpha" && sj.as_string() != "beta"))
                            throw std::runtime_error(cw + ".spin must be 'alpha' or 'beta'");
                        spin = sj.as_string();
                    }
                    bool have_n = false;
                    int64_t n = 0;
                    if (have_rows) {
                        if (rows.count("")) {
                            have_n = true;
                            n = rows.at("");
                        } else if (!spin.empty()) {
                            if (rows.count(spin)) {
                                have_n = true;
                                n = rows.at(spin);
                            }
                        } else {
                            for (const char* ch : {"alpha", "beta"})
                                if (rows.count(ch)) {
                                    n = have_n ? std::min(n, rows.at(ch)) : rows.at(ch);
                                    have_n = true;
                                }
                        }
                    }
                    for (const char* mo_key : {"from_mo", "to_mo"}) {
                        const Json& mo = con.at(mo_key);
                        if (!is_json_int(mo) || mo.as_int() < 0)
                            throw std::runtime_error(cw + "." + mo_key +
                                                     " must be an integer >= 0");
                        if (have_n && mo.as_int() >= n)
                            throw std::runtime_error(
                                cw + "." + mo_key + " " + std::to_string(mo.as_int()) +
                                " is not one of the " + std::to_string(n) +
                                " orbital rows of the reference wavefunction "
                                "(QVF spec § 5.12)");
                    }
                }
                out.set(key, v);
            } else if (key == "nto_ref") {
                need("states[" + std::to_string(k) + "].nto_ref", v, {"wavefunction.gto"});
                out.set(key, v);
            } else if (key == "density_ref") {
                need("states[" + std::to_string(k) + "].density_ref", v,
                     {"volume.transition_density", "volume.difference"});
                out.set(key, v);
            } else {
                out.set(key, v);  // an x_ key
            }
        }
        out_states.push_back(out);
    }

    Json payload = Json::object();
    payload.set("states", out_states);
    if (!meta.is_null()) {
        if (meta.has("units"))
            check_units_object("excited_states meta", meta.at("units"),
                               {{"oscillator_strength", Unit::Dimensionless}});
        if (meta.has("n_states") &&
            (!is_json_int(meta.at("n_states")) ||
             meta.at("n_states").as_int() != static_cast<int64_t>(states.size())))
            throw std::runtime_error("excited_states: meta.n_states must equal the number "
                                     "of states (QVF spec § 5.12)");
        payload.set("meta", meta);
    }
    const std::string sid = reserve_id(section_id);
    Json members = Json::object();
    members.set("states", json_member(sid + "/states.json", payload));
    Json section = Json::object();
    section.set("id", sid).set("kind", "excited_states");
    if (!wavefunction_ref.empty()) section.set("wavefunction_ref", wavefunction_ref);
    section.set("members", members);
    sections_.push_back(section);
    return sid;
}

// -- symmetry & periodic kinds ----------------------------------------------

std::string QvfWriter::add_structure_symmetry(const Json& data,
                                              const std::string& section_id) {
    check_new_id(section_id);
    Json members = Json::object();
    members.set("data", json_member("structure/symmetry.json", data));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "structure.symmetry")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_fermi_surface(const Json& mesh, const Tensor& energies,
                                         const std::string& section_id) {
    check_new_id(section_id);
    if (energies.shape.size() != 4)
        throw std::runtime_error("fermi_surface energies must be 4-D");
    Json members = Json::object();
    members.set("mesh", json_member("fermi/" + section_id + "_mesh.json", mesh));
    members.set("energies", binary_member("fermi/" + section_id + ".bin", energies));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "fermi_surface")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_phonon_bands(const Json& qpath, const Tensor& frequencies,
                                        const Tensor* eigenvectors,
                                        const std::string& section_id) {
    check_new_id(section_id);
    if (eigenvectors) {  // checked before any member is written (spec § 5.6)
        const auto& s = eigenvectors->shape;
        const bool is_real = s.size() == 4 && s[3] == 3;
        const bool is_complex = s.size() == 5 && s[3] == 3 && s[4] == 2;
        if (!is_real && !is_complex)
            throw std::runtime_error(
                "phonon_bands eigenvectors must be [n_q, n_modes, n_atoms, 3] (real) "
                "or [n_q, n_modes, n_atoms, 3, 2] (complex, last axis [re, im])");
        if (eigenvectors->dtype != DType::Float64)
            throw std::runtime_error("phonon_bands eigenvectors must be float64");
        const auto& f = frequencies.shape;
        if (f.size() != 2 || s[0] != f[0] || s[1] != f[1])
            throw std::runtime_error(
                "phonon_bands eigenvectors [n_q, n_modes] must match frequencies "
                "[n_q, n_modes]");
    }
    Json members = Json::object();
    members.set("qpath", json_member("phonons/" + section_id + "_qpath.json", qpath));
    members.set("frequencies",
                binary_member("phonons/" + section_id + "_freq.bin", frequencies));
    if (eigenvectors)
        members.set("eigenvectors",
                    binary_member("phonons/" + section_id + "_eig.bin", *eigenvectors));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "phonon_bands")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_phonon_dos(const Json& meta, const Tensor& frequencies,
                                      const Tensor& dos, const Tensor* projected,
                                      const std::string& section_id) {
    check_new_id(section_id);
    if (projected) {  // checked before any member is written (spec § 5.6)
        const auto& s = projected->shape;
        if (s.size() != 2)
            throw std::runtime_error(
                "phonon_dos projected must be [n_channels, n_points]");
        if (projected->dtype != DType::Float64)
            throw std::runtime_error("phonon_dos projected must be float64");
        if (frequencies.shape.size() != 1 || s[1] != frequencies.shape[0])
            throw std::runtime_error(
                "phonon_dos projected [n_channels, n_points] must match frequencies "
                "[n_points]");
        if (meta.type() == Json::Type::Object && meta.has("channels")) {
            const Json& channels = meta.at("channels");
            if (channels.type() == Json::Type::Array &&
                static_cast<int64_t>(channels.size()) != s[0])
                throw std::runtime_error(
                    "phonon_dos meta channels must have one entry per projected row");
        }
    }
    Json members = Json::object();
    members.set("meta", json_member("phonons/" + section_id + "_meta.json", meta));
    members.set("frequencies",
                binary_member("phonons/" + section_id + "_freq.bin", frequencies));
    members.set("dos", binary_member("phonons/" + section_id + "_dos.bin", dos));
    if (projected)
        members.set("projected",
                    binary_member("phonons/" + section_id + "_proj.bin", *projected));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "phonon_dos")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_equation_of_state(const Tensor& volumes,
                                             const Tensor& energies, const Json& fit,
                                             const std::string& section_id) {
    check_new_id(section_id);
    Json members = Json::object();
    members.set("volumes", binary_member("eos/" + section_id + "_volumes.bin", volumes));
    members.set("energies",
                binary_member("eos/" + section_id + "_energies.bin", energies));
    members.set("fit", json_member("eos/" + section_id + "_fit.json", fit));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "equation_of_state")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

std::string QvfWriter::add_topology_qtaim(const Json& critical_points,
                                          const std::string& section_id) {
    check_new_id(section_id);
    check_payload_units("topology.qtaim", critical_points, {});
    Json members = Json::object();
    members.set("critical_points",
                json_member("topology/" + section_id + ".json", critical_points));
    Json section = Json::object();
    section.set("id", reserve_id(section_id)).set("kind", "topology.qtaim")
           .set("members", members);
    sections_.push_back(section);
    return section_id;
}

// -- property kinds (spec § 5.13) --------------------------------------------
//
// The writer has no JSON Schema engine, so these checks repeat what the payload
// definitions say that can be judged from the object alone (shapes, required
// fields, enumerations, strictness) plus the two rules of the validators that
// involve only the payload (symmetric elastic matrices, static = electronic +
// ionic). The rules that need other parts of the archive (one Born tensor per
// atom, the dipole against the root block) are left to the validators.

// Largest disagreement the validators accept between two statements of one
// property: matrix symmetry (relative to the largest entry) and static against
// electronic + ionic.
static const double kPropertyTolerance = 1e-6;

// The payload as a parsed JSON object (a Json::raw fragment is parsed; malformed
// JSON throws).
static Json property_object(const Json& properties, const std::string& what) {
    Json doc = properties.type() == Json::Type::Raw ? Json::parse(properties.as_string())
                                                    : properties;
    if (!doc.is_object())
        throw std::runtime_error(what + ": the payload must be a JSON object");
    return doc;
}

// A strict payload object: only `allowed` keys and keys beginning x_.
static void check_property_keys(const Json& obj, std::initializer_list<const char*> allowed,
                                const std::string& what) {
    if (!obj.is_object()) throw std::runtime_error(what + " must be a JSON object");
    for (const auto& kv : obj.items()) {
        bool ok = kv.first.compare(0, 2, "x_") == 0;
        for (const char* a : allowed)
            if (kv.first == a) ok = true;
        if (!ok)
            throw std::runtime_error(what + ": unknown key '" + kv.first + "'");
    }
}

static const Json& require_property_member(const Json& obj, const char* key,
                                           const std::string& what) {
    if (!obj.has(key)) throw std::runtime_error(what + " requires " + key);
    return obj.at(key);
}

// `v` as an array of n finite numbers.
static std::vector<double> property_numbers(const Json& v, size_t n, const std::string& what) {
    if (!v.is_array() || v.size() != n)
        throw std::runtime_error(what + " must have " + std::to_string(n) + " numbers");
    std::vector<double> out;
    for (const Json& x : v.elements()) {
        if (!is_finite_number(x)) throw std::runtime_error(what + " must hold finite numbers");
        out.push_back(x.as_double());
    }
    return out;
}

// `v` as rows x cols finite numbers, row-major.
static std::vector<double> property_matrix(const Json& v, size_t rows, size_t cols,
                                           const std::string& what) {
    if (!v.is_array() || v.size() != rows)
        throw std::runtime_error(what + " must be " + std::to_string(rows) + " x " +
                                 std::to_string(cols));
    std::vector<double> out;
    for (const Json& row : v.elements()) {
        if (!row.is_array() || row.size() != cols)
            throw std::runtime_error(what + " must be " + std::to_string(rows) + " x " +
                                     std::to_string(cols));
        for (const double x : property_numbers(row, cols, what)) out.push_back(x);
    }
    return out;
}

// Throws unless the n x n matrix is symmetric to within 1e-6 of its largest
// absolute entry, the rule the validators apply.
static void require_symmetric(const std::vector<double>& m, size_t n, const std::string& what) {
    double scale = 0.0;
    for (double x : m) scale = std::max(scale, std::fabs(x));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            if (std::fabs(m[i * n + j] - m[j * n + i]) > kPropertyTolerance * scale)
                throw std::runtime_error(what + " must be symmetric: entry (" +
                                         std::to_string(i) + ", " + std::to_string(j) +
                                         ") differs from entry (" + std::to_string(j) + ", " +
                                         std::to_string(i) + ") (spec § 5.13)");
}

static void require_one_of(const Json& v, std::initializer_list<const char*> options,
                           const std::string& what) {
    if (v.is_string())
        for (const char* o : options)
            if (v.as_string() == o) return;
    std::string list;
    for (const char* o : options) list += (list.empty() ? "" : ", ") + std::string(o);
    throw std::runtime_error(what + " must be one of " + list);
}

static void validate_electric_properties(const Json& doc) {
    const std::string what = "properties.electric";
    check_property_keys(doc, {"origin", "dipole", "quadrupole", "octupole", "polarizability",
                              "hyperpolarizability", "units"}, what);
    if (doc.has("origin")) {
        const Json& o = doc.at("origin");
        if (o.is_string()) {
            if (o.as_string().empty()) throw std::runtime_error(what + " origin must not be empty");
        } else {
            property_numbers(o, 3, what + " origin");
        }
    }
    bool any_block = false;
    if (doc.has("dipole")) {
        any_block = true;
        const Json& d = doc.at("dipole");
        const std::string w = what + " dipole";
        check_property_keys(d, {"vector", "total", "units"}, w);
        property_numbers(require_property_member(d, "vector", w), 3, w + " vector");
        if (d.has("total") && (!is_finite_number(d.at("total")) || d.at("total").as_double() < 0))
            throw std::runtime_error(w + " total must be a number >= 0");
        if (d.has("units") &&
            !(d.at("units").is_string() && d.at("units").as_string() == "debye"))
            throw std::runtime_error(w + " units must be 'debye'; spec § 5.13 fixes it");
    }
    if (doc.has("quadrupole")) {
        any_block = true;
        const Json& q = doc.at("quadrupole");
        const std::string w = what + " quadrupole";
        check_property_keys(q, {"matrix", "traceless", "convention"}, w);
        const std::vector<double> m =
            property_matrix(require_property_member(q, "matrix", w), 3, 3, w + " matrix");
        require_symmetric(m, 3, w + " matrix");
        const Json& traceless = require_property_member(q, "traceless", w);
        if (!traceless.is_bool()) throw std::runtime_error(w + " traceless must be a boolean");
        const Json& convention = require_property_member(q, "convention", w);
        require_one_of(convention, {"buckingham", "cartesian"}, w + " convention");
        if (convention.as_string() == "buckingham" && !traceless.as_bool())
            throw std::runtime_error(w + ": the buckingham convention is traceless");
    }
    if (doc.has("octupole")) {
        any_block = true;
        const Json& o = doc.at("octupole");
        const std::string w = what + " octupole";
        check_property_keys(o, {"tensor", "convention"}, w);
        property_numbers(require_property_member(o, "tensor", w), 27, w + " tensor");
        if (o.has("convention"))
            require_one_of(o.at("convention"), {"buckingham", "cartesian"}, w + " convention");
    }
    if (doc.has("polarizability")) {
        any_block = true;
        const Json& p = doc.at("polarizability");
        const std::string w = what + " polarizability";
        check_property_keys(p, {"matrix", "frequency_ev", "isotropic"}, w);
        property_matrix(require_property_member(p, "matrix", w), 3, 3, w + " matrix");
        if (p.has("frequency_ev") &&
            (!is_finite_number(p.at("frequency_ev")) || p.at("frequency_ev").as_double() < 0))
            throw std::runtime_error(w + " frequency_ev must be a number >= 0");
        if (p.has("isotropic") && !is_finite_number(p.at("isotropic")))
            throw std::runtime_error(w + " isotropic must be a number");
    }
    if (doc.has("hyperpolarizability")) {
        any_block = true;
        const Json& b = doc.at("hyperpolarizability");
        const std::string w = what + " hyperpolarizability";
        check_property_keys(b, {"tensor", "frequencies_ev"}, w);
        property_numbers(require_property_member(b, "tensor", w), 27, w + " tensor");
        if (b.has("frequencies_ev")) property_numbers(b.at("frequencies_ev"), 2, w + " frequencies_ev");
    }
    if (!any_block)
        throw std::runtime_error(what + " needs at least one of dipole, quadrupole, octupole, "
                                        "polarizability, hyperpolarizability");
    if (doc.has("units"))
        check_units_object(what, doc.at("units"),
                           {{"dipole", Unit::Debye}, {"origin", Unit::Bohr}});
}

static void validate_elastic_properties(const Json& doc) {
    const std::string what = "properties.elastic";
    check_property_keys(doc, {"stiffness", "compliance", "bulk_modulus_gpa", "shear_modulus_gpa",
                              "youngs_modulus_gpa", "poisson_ratio", "average", "units"}, what);
    const std::vector<double> c =
        property_matrix(require_property_member(doc, "stiffness", what), 6, 6, what + " stiffness");
    require_symmetric(c, 6, what + " stiffness");
    if (doc.has("compliance")) {
        const std::vector<double> s =
            property_matrix(doc.at("compliance"), 6, 6, what + " compliance");
        require_symmetric(s, 6, what + " compliance");
    }
    for (const char* key : {"bulk_modulus_gpa", "shear_modulus_gpa", "youngs_modulus_gpa",
                            "poisson_ratio"})
        if (doc.has(key) && !is_finite_number(doc.at(key)))
            throw std::runtime_error(what + " " + key + " must be a number");
    if (doc.has("average"))
        require_one_of(doc.at("average"), {"voigt", "reuss", "hill"}, what + " average");
    if (doc.has("units"))
        check_units_object(what, doc.at("units"),
                           {{"stiffness", Unit::Gigapascal},
                            {"bulk_modulus_gpa", Unit::Gigapascal},
                            {"shear_modulus_gpa", Unit::Gigapascal},
                            {"youngs_modulus_gpa", Unit::Gigapascal},
                            {"poisson_ratio", Unit::Dimensionless}});
}

static void validate_dielectric_properties(const Json& doc) {
    const std::string what = "properties.dielectric";
    check_property_keys(doc, {"electronic", "ionic", "static", "born_charges", "piezoelectric",
                              "units"}, what);
    const std::vector<double> e = property_matrix(
        require_property_member(doc, "electronic", what), 3, 3, what + " electronic");
    std::vector<double> ion, st;
    if (doc.has("ionic")) ion = property_matrix(doc.at("ionic"), 3, 3, what + " ionic");
    if (doc.has("static")) st = property_matrix(doc.at("static"), 3, 3, what + " static");
    if (!ion.empty() && !st.empty())
        for (size_t k = 0; k < 9; ++k)
            if (std::fabs(st[k] - (e[k] + ion[k])) > kPropertyTolerance)
                throw std::runtime_error(what + " static[" + std::to_string(k / 3) + "][" +
                                         std::to_string(k % 3) +
                                         "] is not electronic + ionic; static must equal their "
                                         "sum within 1e-6 (spec § 5.13)");
    if (doc.has("born_charges")) {
        const Json& b = doc.at("born_charges");
        if (!b.is_array() || b.size() == 0)
            throw std::runtime_error(what + " born_charges must be a non-empty array of 3 x 3 "
                                            "tensors, one per atom");
        for (const Json& tensor : b.elements())
            property_matrix(tensor, 3, 3, what + " born_charges tensor");
    }
    if (doc.has("piezoelectric"))
        property_matrix(doc.at("piezoelectric"), 3, 6, what + " piezoelectric");
    if (doc.has("units"))
        check_units_object(what, doc.at("units"),
                           {{"electronic", Unit::Dimensionless},
                            {"ionic", Unit::Dimensionless},
                            {"static", Unit::Dimensionless},
                            {"born_charges", Unit::ElementaryCharge}});
}

std::string QvfWriter::add_properties_section(const std::string& kind,
                                              const std::string& section_id,
                                              const Json& payload) {
    check_new_id(section_id);
    const std::string sid = reserve_id(section_id);
    Json members = Json::object();
    members.set("properties", json_member(sid + "/properties.json", payload));
    Json section = Json::object();
    section.set("id", sid).set("kind", kind).set("members", members);
    sections_.push_back(section);
    return sid;
}

std::string QvfWriter::add_electric_properties(const Json& properties,
                                               const std::string& section_id) {
    validate_electric_properties(property_object(properties, "properties.electric"));
    return add_properties_section("properties.electric", section_id, properties);
}

std::string QvfWriter::add_elastic_properties(const Json& properties,
                                              const std::string& section_id) {
    validate_elastic_properties(property_object(properties, "properties.elastic"));
    return add_properties_section("properties.elastic", section_id, properties);
}

std::string QvfWriter::add_dielectric_properties(const Json& properties,
                                                 const std::string& section_id) {
    validate_dielectric_properties(property_object(properties, "properties.dielectric"));
    return add_properties_section("properties.dielectric", section_id, properties);
}

// -- output -----------------------------------------------------------------

Json QvfWriter::build_manifest() const {
    Json m = Json::object();
    m.set("qvf_version", 1);
    m.set("schema_uri", "https://vibe-qc.org/spec/qvf/1/manifest.schema.json");
    m.set("byte_order", "little");  // the only v1 byte order (spec § 2.2)
    Json src = Json::object();
    src.set("program", source_.program).set("version", source_.version)
       .set("calculation", source_.calculation);
    m.set("source", src);
    // Root metadata (spec § 3.5), with the automatic created_utc / uuid unless
    // set_metadata() supplied them; made once, then reused.
    Json metadata = root_.has("metadata") ? root_.at("metadata") : Json::object();
    if (auto_metadata) {
        if (auto_uuid_.empty()) {
            auto_created_utc_ = metadata_utc_now();
            auto_uuid_ = metadata_random_uuid();
        }
        if (!metadata.has("created_utc")) metadata.set("created_utc", auto_created_utc_);
        if (!metadata.has("uuid")) metadata.set("uuid", auto_uuid_);
    }
    if (metadata.size() > 0) m.set("metadata", metadata);
    // Merge optional root blocks in a deterministic order.
    for (const char* key : {"provenance", "thermochemistry", "dipole_moment",
                            "constraints", "energies", "extensions",
                            "viewer_defaults"}) {
        for (const auto& kv : root_.items())
            if (kv.first == key) m.set(key, kv.second);
    }
    m.set("sections", sections_);
    return m;
}

// -- streaming and in-place update (spec § 2.1; container layer) -------------

static std::string qvf_path(const std::string& path) {
    std::string out = path;
    if (out.size() < 4 || out.substr(out.size() - 4) != ".qvf") out += ".qvf";
    return out;
}

static bool same_file(const std::string& a, const std::string& b) {
    std::error_code ec;
    if (std::filesystem::equivalent(a, b, ec)) return true;
    return a == b;
}

static size_t itemsize_of(DType d) { return dtype_itemsize(d); }

QvfWriter::QvfWriter(Source source, const std::string& path) : QvfWriter(std::move(source)) {
    out_ = std::make_shared<detail::WriterOutput>();
    out_->path = qvf_path(path);
    out_->zip.reset(new container::ZipWriter(
        std::unique_ptr<container::Sink>(new container::FileSink(out_->path))));
}

QvfWriter::QvfWriter(Source source, std::unique_ptr<container::Sink> sink)
    : QvfWriter(std::move(source)) {
    if (!sink) throw std::invalid_argument("QvfWriter: null output sink");
    out_ = std::make_shared<detail::WriterOutput>();
    out_->sink_only = true;
    out_->zip.reset(new container::ZipWriter(std::move(sink)));
}

// Root blocks build_manifest() writes from root_; any other root key of an
// updated archive is kept verbatim.
static bool is_modelled_root_block(const std::string& key) {
    for (const char* k : {"metadata", "provenance", "thermochemistry", "dipole_moment",
                          "constraints", "energies", "extensions", "viewer_defaults"})
        if (key == k) return true;
    return false;
}

QvfWriter QvfWriter::open_for_update(const std::string& path) {
    container::ZipReader zip = container::ZipReader::open(path);
    // Superseded manifest.json copies of an earlier update are expected; any
    // other repeated name makes the archive ambiguous (spec § 2.1, § 2.3).
    for (const std::string& name : zip.duplicates())
        if (name != "manifest.json")
            throw std::runtime_error("E_ZIP_DUPLICATE_ENTRY: duplicate ZIP entry name '" +
                                     name + "'; refusing to update");
    const container::Limits lim{};
    if (!zip.has("manifest.json"))
        throw std::runtime_error("not a QVF archive (no manifest.json): " + path);
    const std::vector<uint8_t> raw = zip.read("manifest.json", lim.max_json_bytes);
    const Json m = Json::parse(std::string(raw.begin(), raw.end()), lim.max_json_depth);
    if (!m.is_object()) throw std::runtime_error("manifest.json is not a JSON object");
    // The withdrawn qvf_version 2 is read as 1 (spec § 3.1); the update
    // writes 1 and the v1 schema_uri.
    const double version = m.has("qvf_version") && m.at("qvf_version").is_number()
                               ? m.at("qvf_version").as_double() : 0.0;
    if (version != 1.0 && version != 2.0)
        throw std::runtime_error("only qvf_version 1 archives (or the withdrawn 2, read as 1) "
                                 "can be updated");
    const bool withdrawn_version = version == 2.0;
    const Json* src = m.has("source") ? &m.at("source") : nullptr;
    if (!src || !src->is_object() || !src->has("program") || !src->has("version") ||
        !src->at("program").is_string() || !src->at("version").is_string())
        throw std::runtime_error("manifest.json has no usable source block");
    Source source{src->at("program").as_string(), src->at("version").as_string(),
                  src->has("calculation") && src->at("calculation").is_string()
                      ? src->at("calculation").as_string() : std::string()};
    if (!m.has("sections") || !m.at("sections").is_array())
        throw std::runtime_error("manifest.json has no sections array");
    if (m.has("byte_order") &&
        !(m.at("byte_order").is_string() && m.at("byte_order").as_string() == "little"))
        throw std::runtime_error("byte_order is not 'little' (spec § 2.2); refusing to update");

    QvfWriter w(std::move(source));
    w.auto_metadata = false;   // created_utc / uuid stay as the archive has them
    std::shared_ptr<detail::WriterOutput> out = std::make_shared<detail::WriterOutput>();
    out->path = path;
    out->update = true;
    out->source = *src;
    out->has_source = true;
    for (const auto& kv : m.items()) {
        if (kv.first == "qvf_version" || kv.first == "source" || kv.first == "sections")
            continue;
        if (kv.first == "schema_uri" && withdrawn_version)
            continue;   // build_manifest() writes the v1 URI
        if (is_modelled_root_block(kv.first))
            w.root_.set(kv.first, kv.second);
        else
            out->preserved_root.set(kv.first, kv.second);
    }
    w.sections_ = m.at("sections");
    for (const Json& s : w.sections_.elements()) {
        if (!s.is_object() || !s.has("id") || !s.at("id").is_string() || !s.has("members") ||
            !s.at("members").is_object())
            throw std::runtime_error("manifest.json has a malformed section");
        w.ids_.insert(s.at("id").as_string());
        for (const auto& mv : s.at("members").items())
            if (mv.second.is_object() && mv.second.has("path") &&
                mv.second.at("path").is_string())
                w.paths_.insert(mv.second.at("path").as_string());
    }
    // Every existing entry name is taken: a new member may not reuse one.
    for (const std::string& name : zip.names()) {
        w.paths_.insert(name);
        out->archived.insert(name);
    }
    if (m.has("metadata") && m.at("metadata").is_object() &&
        m.at("metadata").has("modified_utc"))
        out->loaded_modified = m.at("metadata").at("modified_utc");
    // The custom atom_properties roles described so far, so that more can
    // be added (best effort: a malformed meta member is left as it is).
    for (const Json& s : w.sections_.elements()) {
        if (!s.has("kind") || !s.at("kind").is_string() ||
            s.at("kind").as_string() != "atom_properties" || !s.at("members").has("meta"))
            continue;
        try {
            const Json& spec = s.at("members").at("meta");
            const std::vector<uint8_t> b = zip.read(spec.at("path").as_string(), lim);
            w.set_atom_properties_meta(Json::parse(std::string(b.begin(), b.end())));
        } catch (const std::exception&) {
        }
    }
    // Σ Z of each structure, for set_structure_identity (best effort, as above).
    for (const Json& s : w.sections_.elements()) {
        if (!s.has("kind") || !s.at("kind").is_string() ||
            s.at("kind").as_string() != "structure" || !s.at("members").has("structure"))
            continue;
        try {
            const Json& spec = s.at("members").at("structure");
            const std::vector<uint8_t> b = zip.read(spec.at("path").as_string(), lim);
            int64_t sum_z = 0;
            if (payload_sum_z(Json::parse(std::string(b.begin(), b.end()), lim.max_json_depth),
                              sum_z))
                w.structure_sum_z_.emplace(s.at("id").as_string(), sum_z);
        } catch (const std::exception&) {
        }
    }
    out->original = std::make_unique<container::ZipReader>(std::move(zip));
    w.out_ = std::move(out);
    return w;
}

// Start appending to the bound archive: reopen it after its last entry and
// flush the members buffered so far.
void QvfWriter::open_bound_archive() {
    if (out_ && out_->sink_only && out_->zip->finished())
        throw std::runtime_error("QvfWriter: output sink is already finished");
    if (!out_ || out_->zip) return;
    out_->zip.reset(new container::ZipWriter(container::ZipWriter::append(out_->path)));
    for (const auto& f : files_)
        out_->zip->write(f.first, f.second,
                         stored_.count(f.first) ? container::Method::Store
                                                : container::Method::Deflate,
                         force_zip64);
    files_.clear();
}

// The manifest to write. An update keeps the archive's source and unmodelled
// root keys, and stamps metadata.modified_utc with the time of the update
// (spec § 2.3) unless the caller set a value of its own.
Json QvfWriter::final_manifest() const {
    Json m = build_manifest();
    if (out_ && out_->update) {
        if (out_->has_source) m.set("source", out_->source);
        for (const auto& kv : out_->preserved_root.items()) m.set(kv.first, kv.second);
        Json metadata = m.has("metadata") ? m.at("metadata") : Json::object();
        const bool caller_set =
            metadata.has("modified_utc") &&
            (out_->loaded_modified.is_null() ||
             metadata.at("modified_utc").dump(0) != out_->loaded_modified.dump(0));
        if (!caller_set && metadata.is_object()) {
            if (out_->stamp.empty()) out_->stamp = metadata_utc_now();
            metadata.set("modified_utc", out_->stamp);
            m.set("metadata", metadata);
        }
    }
    return m;
}

// The buffered members, then manifest.json (STOREd, last), then the central
// directory. When `zip` appends to an archive (spec § 2.3), the previous
// manifest stays listed as a superseded entry: readers use the last one.
void QvfWriter::emit_archive(container::ZipWriter& zip) {
    if (member_open_)
        throw std::runtime_error("a member sink is open; close it before writing");
    if (force_zip64) zip.set_force_zip64(true);
    const std::string manifest = final_manifest().dump(2);
    for (const auto& f : files_)
        zip.write(f.first, f.second,
                  stored_.count(f.first) ? container::Method::Store : container::Method::Deflate,
                  force_zip64);
    if (out_ && out_->update) zip.allow_repeated_name("manifest.json");
    zip.write("manifest.json", manifest, container::Method::Store, force_zip64);
    zip.finish();
}

std::vector<uint8_t> QvfWriter::to_bytes() {
    if (out_) {
        if (out_->zip || !out_->update)
            throw std::runtime_error("to_bytes() is not available while members stream "
                                     "to a file; call write()");
        // An update not yet started: the updated archive, in memory.
        std::ifstream fh(out_->path, std::ios::binary);
        if (!fh) throw std::runtime_error("cannot open: " + out_->path);
        std::vector<uint8_t> current((std::istreambuf_iterator<char>(fh)),
                                     std::istreambuf_iterator<char>());
        container::ZipWriter zip = container::ZipWriter::append(std::move(current));
        emit_archive(zip);
        return static_cast<container::MemorySink&>(zip.sink()).take();
    }
    container::MemorySink* mem = new container::MemorySink;
    container::ZipWriter zip{std::unique_ptr<container::Sink>(mem)};
    emit_archive(zip);
    return mem->take();
}

std::string QvfWriter::write(const std::string& path) {
    const std::string out = qvf_path(path);
    if (out_) {
        if (!same_file(out, out_->path) && !same_file(path, out_->path))
            throw std::runtime_error("this writer is bound to " + out_->path +
                                     "; call write() to finish it");
        return write();
    }
    if (member_open_)
        throw std::runtime_error("a member sink is open; close it before writing");
    container::ZipWriter zip{std::unique_ptr<container::Sink>(new container::FileSink(out))};
    emit_archive(zip);
    return out;
}

std::string QvfWriter::write() {
    if (!out_)
        throw std::runtime_error("write() needs a path unless the writer is bound to a "
                                 "file (QvfWriter(source, path) or open_for_update)");
    if (member_open_)
        throw std::runtime_error("a member sink is open; close it before writing");
    open_bound_archive();
    const Json written = final_manifest();
    emit_archive(*out_->zip);
    if (out_->sink_only) return {}; // retain the sink for a memory caller
    // A later write() appends to the finished archive: every name in it is
    // taken, and its modified_utc is the one just written.
    for (const container::EntryInfo& e : out_->zip->entries()) out_->archived.insert(e.name);
    out_->zip.reset();
    out_->update = true;
    if (written.has("metadata")) {
        root_.set("metadata", written.at("metadata"));
        out_->loaded_modified = written.at("metadata").has("modified_utc")
                                    ? written.at("metadata").at("modified_utc") : Json();
    }
    out_->stamp.clear();
    return out_->path;
}

namespace {
void generic_member_role(const std::string& role) {
    if (role.empty() || !((role[0] >= 'A' && role[0] <= 'Z') ||
        (role[0] >= 'a' && role[0] <= 'z') || (role[0] >= '0' && role[0] <= '9') || role[0] == '_'))
        throw std::invalid_argument("invalid member role: " + role);
    for (char c : role)
        if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
              (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.'))
            throw std::invalid_argument("invalid member role: " + role);
}
}
std::string QvfWriter::add_section(Section section) {
    if (member_open_) throw std::runtime_error("a member sink is open");
    Json s = section.to_json();
    if (!s.has("id") || !s.at("id").is_string() ||
        !s.has("kind") || !s.at("kind").is_string() || s.at("kind").as_string().empty())
        throw std::invalid_argument("add_section requires id and kind");
    const std::string id = s.at("id").as_string();
    check_new_id(id);
    if (!s.has("members")) s.set("members", Json::object());
    if (!s.at("members").is_object()) throw std::invalid_argument("members must be an object");
    reserve_id(id); sections_.push_back(std::move(s)); return id;
}
void QvfWriter::add_json_member(const std::string& id, const std::string& role, Json value) {
    if (member_open_) throw std::runtime_error("a member sink is open");
    generic_member_role(role);
    const size_t index = section_index(id);
    Json s = sections_.elements()[index];
    Json members = s.at("members");
    if (members.has(role)) throw std::invalid_argument("duplicate member role: " + role);
    members.set(role, json_member(id + "/" + role + ".json", value));
    s.set("members", std::move(members)); replace_section(index, std::move(s));
}
void QvfWriter::add_section_opaque(const std::string& id, const std::string& role,
                                 const std::vector<uint8_t>& bytes,
                                 const std::string& media_type) {
    if (member_open_) throw std::runtime_error("a member sink is open");
    generic_member_role(role);
    const size_t index = section_index(id);
    Json s = sections_.elements()[index];
    Json members = s.at("members");
    if (members.has(role)) throw std::invalid_argument("duplicate member role: " + role);
    members.set(role, opaque_member(id + "/" + role + ".bin", bytes, media_type));
    s.set("members", std::move(members)); replace_section(index, std::move(s));
}

MemberSink QvfWriter::begin_binary(const std::string& section_id, const std::string& role,
                                   DType dtype, std::vector<int64_t> shape) {
    if (member_open_)
        throw std::runtime_error("a member sink is already open");
    const Json& section = sections_.elements()[section_index(section_id)];
    if (role.empty()) throw std::runtime_error("begin_binary: empty role");
    if (section.at("members").has(role))
        throw std::runtime_error("section '" + section_id + "' already has a member '" +
                                 role + "'");
    uint64_t expected = itemsize_of(dtype);
    for (int64_t d : shape) {
        if (d < 0) throw std::runtime_error("begin_binary: negative extent in shape");
        const uint64_t u = static_cast<uint64_t>(d);
        if (u != 0 && expected > UINT64_MAX / u)
            throw std::runtime_error("begin_binary: shape too large");
        expected *= u;
    }
    const std::string path = entry_name(section_id + "/" + role + ".bin");

    std::unique_ptr<MemberSink::State> st(new MemberSink::State);
    st->writer = this;
    st->section_id = section_id;
    st->role = role;
    st->path = path;
    st->dtype = dtype;
    st->shape = std::move(shape);
    st->expected = expected;
    if (out_) {
        open_bound_archive();
        // DEFLATE in blocks adds at most 5 bytes per 32 KiB; begin with ZIP64
        // when the entry may reach 4 GiB.
        const bool zip64 = force_zip64 || expected >= 0xFFFFFFFFull - expected / 4096 - 65536;
        st->entry.reset(new container::EntrySink(
            out_->zip->begin(path, container::Method::Deflate, zip64)));
    }
    paths_.insert(path);
    member_open_ = true;
    return MemberSink(std::move(st));
}

// -- MemberSink ---------------------------------------------------------------

MemberSink::MemberSink(std::unique_ptr<State> state) : s_(std::move(state)) {}
MemberSink::MemberSink(MemberSink&&) noexcept = default;
MemberSink& MemberSink::operator=(MemberSink&& other) noexcept {
    if (this != &other) {
        abandon();
        s_ = std::move(other.s_);
    }
    return *this;
}
MemberSink::~MemberSink() { abandon(); }

void MemberSink::abandon() noexcept {
    if (!s_ || !s_->open) return;
    s_->open = false;
    s_->writer->member_open_ = false;
    s_->writer->paths_.erase(s_->path);
    s_->entry.reset();          // an unfinished ZIP entry is left out of the directory
    s_->buffer.clear();
}

uint64_t MemberSink::bytes_written() const { return s_ ? s_->written : 0; }
uint64_t MemberSink::bytes_expected() const { return s_ ? s_->expected : 0; }
const std::string& MemberSink::path() const {
    static const std::string empty;
    return s_ ? s_->path : empty;
}

void MemberSink::write(const void* data, size_t n) {
    if (!s_ || !s_->open) throw std::runtime_error("write to a closed member sink");
    State& st = *s_;
    if (n > st.expected - st.written)
        throw std::runtime_error("member " + st.path + ": more bytes than dtype x shape (" +
                                 std::to_string(st.expected) + ")");
    const uint8_t* p = static_cast<const uint8_t*>(data);
    if (st.entry) {
        st.entry->write(p, n);
    } else {
        st.buffer.insert(st.buffer.end(), p, p + n);
        st.sha.update(p, n);
    }
    st.written += n;
}

void MemberSink::close() {
    if (!s_ || !s_->open) throw std::runtime_error("member sink already closed");
    State& st = *s_;
    if (st.written != st.expected) {
        const std::string msg = "member " + st.path + ": " + std::to_string(st.written) +
                                " bytes written, dtype x shape needs " +
                                std::to_string(st.expected);
        abandon();
        throw std::runtime_error(msg);
    }
    std::string sha;
    if (st.entry) {
        try {
            sha = st.entry->close().sha256;
        } catch (...) {
            abandon();
            throw;
        }
    } else {
        sha = st.sha.hex();
    }
    QvfWriter& w = *st.writer;
    st.open = false;
    w.member_open_ = false;
    Json spec = Json::object();
    spec.set("path", st.path).set("format", "binary");
    spec.set("dtype", std::string(dtype_name(st.dtype)));
    Json shape = Json::array();
    for (int64_t d : st.shape) shape.push_back(Json(static_cast<long long>(d)));
    spec.set("shape", shape).set("sha256", sha);
    spec.set("size", static_cast<long long>(st.written));
    if (!st.entry) w.files_.emplace_back(st.path, std::move(st.buffer));
    // Json has no mutable element access: rebuild the sections array.
    Json rebuilt = Json::array();
    for (const Json& s : w.sections_.elements()) {
        if (s.has("id") && s.at("id").is_string() && s.at("id").as_string() == st.section_id) {
            Json members = s.at("members");
            members.set(st.role, spec);
            Json updated = s;
            updated.set("members", members);
            rebuilt.push_back(updated);
        } else {
            rebuilt.push_back(s);
        }
    }
    w.sections_ = std::move(rebuilt);
}

// -- repack -------------------------------------------------------------------

void repack(const std::string& in, const std::string& out) {
    const bool in_place = same_file(in, out);
    const std::string target = in_place ? out + ".repack-tmp" : out;
    {
        container::ZipReader zip = container::ZipReader::open(in);
        // Superseded manifest.json copies are what repack removes (spec
        // § 2.3); any other repeated name makes the archive ambiguous.
        for (const std::string& name : zip.duplicates())
            if (name != "manifest.json")
                throw std::runtime_error("E_ZIP_DUPLICATE_ENTRY: duplicate ZIP entry name '" +
                                         name + "'; refusing to repack");
        const container::Limits lim{};
        if (!zip.has("manifest.json"))
            throw std::runtime_error("not a QVF archive (no manifest.json): " + in);
        const std::vector<uint8_t> manifest = zip.read("manifest.json", lim.max_json_bytes);
        const Json m =
            Json::parse(std::string(manifest.begin(), manifest.end()), lim.max_json_depth);
        std::set<std::string> referenced;
        if (m.is_object() && m.has("sections") && m.at("sections").is_array())
            for (const Json& s : m.at("sections").elements())
                if (s.is_object() && s.has("members") && s.at("members").is_object())
                    for (const auto& mv : s.at("members").items())
                        if (mv.second.is_object() && mv.second.has("path") &&
                            mv.second.at("path").is_string())
                            referenced.insert(mv.second.at("path").as_string());
        container::ZipWriter w{
            std::unique_ptr<container::Sink>(new container::FileSink(target))};
        const uint64_t kWhole = uint64_t(16) << 20;   // larger members are streamed
        for (const container::ZipEntry& e : zip.entries()) {
            if (e.name == "manifest.json" || !referenced.count(e.name)) continue;
            const container::Method method =
                e.method == 0 ? container::Method::Store : container::Method::Deflate;
            if (e.usize <= kWhole) {
                w.write(e.name, zip.read(e.name, lim), method, e.zip64);
            } else {
                const bool zip64 = e.zip64 || e.usize >= 0xFFFFFFFFull - e.usize / 4096 - 65536;
                container::EntrySink sink = w.begin(e.name, method, zip64);
                zip.read_into(e.name, [&sink](const uint8_t* p, size_t n) { sink.write(p, n); },
                              lim);
                sink.close();
            }
        }
        w.write("manifest.json", manifest, container::Method::Store);
        w.finish();
    }
    if (in_place) {
        std::filesystem::remove(out);
        std::filesystem::rename(target, out);
    }
}

}  // namespace qvf
