// test_reader.cpp — validator hardening tests for the QVF C++ reader.
//
// Two modes:
//   qvf_reader_tests                    in-memory unit tests: malformed archives
//                                       are built with detail::ZipWriter (or by
//                                       patching its bytes) and each must yield
//                                       its stable finding code (spec App. B).
//   qvf_reader_tests --negative DIR     runs the negative conformance corpus
//                                       (conformance/negative): every archive
//                                       must produce exactly the codes its
//                                       expected.json names — the same codes the
//                                       Python reference validator reports —
//                                       once with the generated JSON Schema
//                                       checker and once without it.
//                                       Exits 77 (skipped) if DIR is absent.
//   qvf_reader_tests --model-roundtrip DIR...
//                                       the typed manifest model of every
//                                       archive in the directories that opens
//                                       round-trips to its manifest.json
//                                       exactly (qvf/model.hpp).
//
// The unit tests of the semantic checks validate with the schema check
// disabled (semantic_findings(), which reports W_SCHEMA_SKIPPED), so each
// archive isolates one rule; test_schema_checker() covers the checker.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "qvf/qvf.hpp"
#include "qvf/qvf_reader.hpp"
#include "qvf_inflate.hpp"
#include "qvf_deflate.hpp"
#include "qvf_schema_check.hpp"
#include "qvf_sha256.hpp"
#include "qvf_zip.hpp"

using namespace qvf;

static int g_failures = 0;

#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::cerr << "FAIL: " << (msg) << " (" << #cond << ")\n";      \
            ++g_failures;                                                  \
        }                                                                  \
    } while (0)

// --- archive construction helpers -------------------------------------------

struct Entry {
    std::string name;
    std::string data;
    bool compress = true;
};

static std::vector<uint8_t> make_zip(const std::vector<Entry>& entries) {
    detail::ZipWriter zw;
    for (const auto& e : entries) zw.add(e.name, e.data, e.compress);
    return zw.finish();
}

static std::string sha(const std::string& s) { return detail::Sha256::hex_of(s); }

template <class F>
static bool throws_runtime(F fn) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

static Json json_member(const std::string& path, const std::string& data) {
    Json m = Json::object();
    m.set("path", path).set("format", "json").set("sha256", sha(data));
    return m;
}

static Json binary_member(const std::string& path, const std::string& data,
                          const std::string& dtype, int64_t n) {
    Json shape = Json::array();
    shape.push_back(Json(static_cast<long long>(n)));
    Json m = Json::object();
    m.set("path", path).set("format", "binary").set("dtype", dtype)
        .set("shape", shape).set("sha256", sha(data));
    return m;
}

static Json section(const std::string& id, const std::string& kind, Json members) {
    Json s = Json::object();
    s.set("id", id).set("kind", kind).set("members", std::move(members));
    return s;
}

static const std::string kWater =
    R"({"atoms":[{"symbol":"O","position":[0.0,0.0,0.1173],"atomic_number":8},)"
    R"({"symbol":"H","position":[0.0,0.7572,-0.4692],"atomic_number":1}],)"
    R"("pbc":[false,false,false],"lattice_vectors":null})";

static Json structure_section(const std::string& id = "structure",
                              const std::string& path = "structure/structure.json",
                              const std::string& data = kWater) {
    Json members = Json::object();
    members.set("structure", json_member(path, data));
    return section(id, "structure", members);
}

static Json manifest_of(const std::vector<Json>& sections, Json extra = Json::object()) {
    Json src = Json::object();
    src.set("program", "t").set("version", "1").set("calculation", "c");
    Json arr = Json::array();
    for (const auto& s : sections) arr.push_back(s);
    Json m = Json::object();
    m.set("qvf_version", 1).set("source", src).set("sections", arr);
    for (const auto& kv : extra.items()) m.set(kv.first, kv.second);
    return m;
}

// The archive: payload entries followed by a STOREd manifest.json.
static std::vector<uint8_t> archive(const Json& manifest, std::vector<Entry> files) {
    files.push_back({"manifest.json", manifest.dump(), false});
    return make_zip(files);
}

static std::vector<uint8_t> base_archive() {
    return archive(manifest_of({structure_section()}),
                   {{"structure/structure.json", kWater}});
}

// Validation without the JSON Schema check: the semantic checks alone, as a
// validator without a schema engine runs them (W_SCHEMA_SKIPPED).
static ValidationOptions no_schema() {
    ValidationOptions o;
    o.schema_check = false;
    return o;
}

static std::vector<Finding> semantic_findings(std::vector<uint8_t> bytes,
                                              const Limits& lim = Limits{}) {
    return QvfReader::validate_bytes(std::move(bytes), lim, no_schema());
}

static std::set<std::string> codes(const std::vector<Finding>& fs) {
    std::set<std::string> out;
    for (const auto& f : fs) out.insert(f.code);
    return out;
}

static std::string describe(const std::vector<Finding>& fs) {
    std::string s;
    for (const auto& f : fs) s += "\n    " + f.code + " [" + f.where + "] " + f.message;
    return s;
}

static bool has_code(const std::vector<Finding>& fs, const std::string& code,
                     Finding::Severity sev = Finding::Severity::Error) {
    for (const auto& f : fs)
        if (f.code == code && f.severity == sev) return true;
    return false;
}

// Expect `code` (and nothing else but W_SCHEMA_SKIPPED) from validating bytes.
static void expect_only(const std::string& label, std::vector<uint8_t> bytes,
                        const std::string& code, const Limits& lim = Limits{},
                        Finding::Severity sev = Finding::Severity::Error) {
    std::vector<Finding> fs = semantic_findings(std::move(bytes), lim);
    std::set<std::string> got = codes(fs);
    got.erase("W_SCHEMA_SKIPPED");
    bool ok = has_code(fs, code, sev) && got == std::set<std::string>{code} &&
              has_errors(fs) == (sev == Finding::Severity::Error);
    if (!ok) {
        std::cerr << "FAIL: " << label << ": expected only " << code << ", got"
                  << describe(fs) << "\n";
        ++g_failures;
    }
}

// Expect a valid archive: no finding at all besides W_SCHEMA_SKIPPED.
static void expect_clean(const std::string& label, std::vector<uint8_t> bytes,
                         const Limits& lim = Limits{}) {
    std::vector<Finding> fs = semantic_findings(std::move(bytes), lim);
    if (codes(fs) != std::set<std::string>{"W_SCHEMA_SKIPPED"}) {
        std::cerr << "FAIL: " << label << ": expected a clean archive, got"
                  << describe(fs) << "\n";
        ++g_failures;
    }
}

// Overwrite one entry's declared uncompressed size in both ZIP headers.
static void patch_usize(std::vector<uint8_t>& z, const std::string& name, uint32_t size) {
    int patched = 0;
    for (size_t p = 0; p + 46 <= z.size(); ++p) {
        uint32_t sig = z[p] | (z[p + 1] << 8) | (z[p + 2] << 16) |
                       (static_cast<uint32_t>(z[p + 3]) << 24);
        size_t size_off, nlen_off, name_off;
        if (sig == 0x04034b50u) { size_off = 22; nlen_off = 26; name_off = 30; }
        else if (sig == 0x02014b50u) { size_off = 24; nlen_off = 28; name_off = 46; }
        else continue;
        size_t n = z[p + nlen_off] | (z[p + nlen_off + 1] << 8);
        if (p + name_off + n > z.size()) continue;
        if (std::string(z.begin() + p + name_off, z.begin() + p + name_off + n) != name)
            continue;
        for (int b = 0; b < 4; ++b)
            z[p + size_off + b] = static_cast<uint8_t>((size >> (8 * b)) & 0xFF);
        ++patched;
    }
    CHECK(patched == 2, "patch_usize patched local + central headers");
}

static std::string nested(int depth) {
    return std::string(depth, '[') + std::string(depth, ']');
}

// --- unit tests ---------------------------------------------------------------

static void test_json_depth() {
    CHECK(Json::parse(nested(64)).is_array(), "depth 64 parses with the default limit");
    bool threw = false;
    try {
        Json::parse(nested(65));
    } catch (const JsonDepthError&) {
        threw = true;
    }
    CHECK(threw, "depth 65 throws JsonDepthError");
    threw = false;
    try {
        Json::parse(R"({"a":[{"b":1}]})", 2);
    } catch (const JsonDepthError&) {
        threw = true;
    }
    CHECK(threw, "explicit max_depth 2 rejects depth 3");
    CHECK(Json::parse(R"({"s":"[[[[{{{{"})", 1).is_object(),
          "brackets inside strings do not count toward depth");
    threw = false;
    try {
        Json::parse("[1,");
    } catch (const JsonDepthError&) {
        threw = false;
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw, "malformed JSON throws a plain runtime_error, not JsonDepthError");
}

static void test_inflate_bound() {
    std::vector<uint8_t> zeros(3000, 0);
    std::vector<uint8_t> comp = detail::deflate(zeros);
    CHECK(detail::inflate(comp.data(), comp.size(), zeros.size(), zeros.size()) == zeros,
          "inflate within bound");
    bool threw = false;
    try {
        detail::inflate(comp.data(), comp.size(), 100, 100);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw, "inflate past max_out throws");
}

static void test_good_archive() {
    QvfWriter w({"t", "1", "good"});
    w.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}}});
    w.add_atom_properties({-0.68, 0.34});
    std::vector<uint8_t> bytes = w.to_bytes();
    std::vector<Finding> fs = semantic_findings(bytes);
    CHECK(!has_errors(fs), "writer output has no errors" + describe(fs));
    CHECK(has_code(fs, "W_SCHEMA_SKIPPED", Finding::Severity::Warning),
          "without the schema check the validator reports W_SCHEMA_SKIPPED");
    const std::vector<Finding> full = QvfReader::validate_bytes(bytes);
    CHECK(full.empty(), "writer output is schema-clean" + describe(full));
    QvfReader r = QvfReader::from_bytes(bytes);
    CHECK(!has_errors(r.validate_findings()), "validate_findings() clean");
    std::vector<std::string> errs;
    CHECK(r.validate(errs) && errs.empty(), "legacy validate() clean");
    expect_clean("base archive", base_archive());

    // The writer refuses what the validator rejects.
    bool refused = false;
    try {
        QvfWriter w2({"t", "1", "reserved"});
        w2.add_structure({{"H", 1, {{0, 0, 0}}}}, {{false, false, false}}, {}, {},
                         "bookmarks");
    } catch (const std::exception&) {
        refused = true;
    }
    CHECK(refused, "writer refuses a reserved section id");
    refused = false;
    try {
        QvfWriter w3({"t", "1", "bad-role"});
        Json m = Json::object();
        m.set("a\\b", Json::object());
        w3.add_vendor_section("x_test.thing", m);
    } catch (const std::exception&) {
        refused = true;
    }
    CHECK(refused, "writer refuses a member path with a backslash");
}

static void test_container_rules() {
    // Duplicate entry names: reported by the validator, refused by the reader.
    std::vector<uint8_t> dup = archive(manifest_of({structure_section()}),
                                       {{"structure/structure.json", "{\"atoms\":[]}"},
                                        {"structure/structure.json", kWater}});
    expect_only("duplicate entry", dup, "E_ZIP_DUPLICATE_ENTRY");
    bool threw = false;
    try {
        QvfReader::from_bytes(dup);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("E_ZIP_DUPLICATE_ENTRY") != std::string::npos;
    }
    CHECK(threw, "QvfReader refuses duplicate entry names");

    for (const std::string path : {"structure\\structure.json", "C:/structure.json",
                                   "c:structure.json", "/structure.json",
                                   "../structure.json", "a/../b.json", "a/.."}) {
        expect_only("path " + path,
                    archive(manifest_of({structure_section("structure", path)}),
                            {{path, kWater}}),
                    "E_MEMBER_PATH");
    }
    // Legal look-alikes stay legal.
    for (const std::string path : {"a..b/structure.json", "dir/..hidden.json",
                                   "dir/.../structure.json"}) {
        std::vector<Finding> fs = semantic_findings(
            archive(manifest_of({structure_section("structure", path)}), {{path, kWater}}));
        CHECK(!has_errors(fs), "legal path " + path + describe(fs));
    }

    // Caps, with deliberately lowered limits.
    std::string bin(4096, '\x5a');
    Json big = Json::object();
    big.set("data", binary_member("x/big.bin", bin, "uint8", 4096));
    std::vector<uint8_t> with_big =
        archive(manifest_of({structure_section(), section("big", "x_test.big", big)}),
                {{"structure/structure.json", kWater}, {"x/big.bin", bin}});
    Limits member_cap;
    member_cap.max_member_bytes = 2048;
    expect_only("member cap", with_big, "E_MEMBER_TOO_LARGE", member_cap);
    Limits total_cap;
    total_cap.max_total_bytes = 4096;
    expect_only("total cap", with_big, "E_TOTAL_TOO_LARGE", total_cap);
    Limits ratio_cap;
    ratio_cap.max_ratio = 10;
    ratio_cap.ratio_exempt_bytes = 1024;
    expect_only("ratio cap", with_big, "E_COMPRESSION_RATIO", ratio_cap);
    ratio_cap.ratio_exempt_bytes = 8192;   // entry is under the exemption floor
    expect_clean("ratio exemption", with_big, ratio_cap);
    expect_clean("default caps accept it", with_big);

    // A STOREd entry is exempt from the ratio check whatever its size (a cap
    // below 1 would otherwise flag every entry, compressed or not).
    std::vector<uint8_t> stored =
        archive(manifest_of({structure_section(), section("big", "x_test.big", big)}),
                {{"structure/structure.json", kWater}, {"x/big.bin", bin, false}});
    Limits tiny_ratio;
    tiny_ratio.max_ratio = 0.5;
    tiny_ratio.ratio_exempt_bytes = 1024;
    expect_clean("stored exempt", stored, tiny_ratio);

    // Declared size smaller than the real inflated size: never inflate past it.
    std::vector<uint8_t> liar = with_big;
    patch_usize(liar, "x/big.bin", 1024);
    expect_only("understated size", liar, "E_ZIP_UNREADABLE");

    // Orphan entry: a warning, the archive stays valid.
    std::vector<Finding> fs = semantic_findings(
        archive(manifest_of({structure_section()}),
                {{"structure/structure.json", kWater}, {"stale/old.json", kWater},
                 {"stale/", ""}}));
    CHECK(has_code(fs, "W_ORPHAN_MEMBER", Finding::Severity::Warning) && !has_errors(fs),
          "orphan is a warning" + describe(fs));
    size_t orphans = 0;
    for (const auto& f : fs) orphans += f.code == "W_ORPHAN_MEMBER";
    CHECK(orphans == 1, "directory entries are not orphans");

    // Broken containers never throw from validate_bytes().
    std::vector<uint8_t> trunc = base_archive();
    trunc.resize(trunc.size() - 30);
    fs = semantic_findings(trunc);
    CHECK(has_code(fs, "E_ZIP_UNREADABLE"), "truncated archive is unreadable");
    fs = semantic_findings(std::vector<uint8_t>{'n', 'o', 't', ' ', 'z', 'i', 'p'});
    CHECK(has_code(fs, "E_ZIP_UNREADABLE"), "garbage is unreadable");
    fs = QvfReader::validate_file("no/such/file.qvf");
    CHECK(has_code(fs, "E_ZIP_UNREADABLE"), "missing file is unreadable");
}

// Member size, superseded manifests and byte_order (spec § 2.2, § 2.3, § 4.2).
static void test_size_updates_and_byte_order() {
    // The writer declares `size` on every member and `byte_order` at the root.
    QvfWriter w({"t", "1", "sizes"});
    w.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}}});
    w.add_atom_properties({-0.68, 0.34});
    w.add_citations("@misc{x, title={y}}\n");
    RunRecord rec;
    rec.program = "democode";
    rec.input_text = "input\n";
    rec.log_text = "log\n";
    rec.attachments.push_back({"geom", std::vector<uint8_t>{'x', 'y', 'z'}});
    rec.files = Json::object();
    w.add_run_record(rec);
    std::vector<uint8_t> written = w.to_bytes();
    expect_clean("writer output with sizes", written);
    QvfReader wr = QvfReader::from_bytes(written);
    CHECK(wr.manifest().contains("byte_order") &&
              wr.manifest().at("byte_order").as_string() == "little",
          "writer emits byte_order little");
    size_t n_members = 0;
    for (const Json& s : wr.sections().elements())
        for (const auto& kv : s.at("members").items()) {
            ++n_members;
            const std::string sid = s.at("id").as_string();
            CHECK(kv.second.contains("size") &&
                      static_cast<size_t>(kv.second.at("size").as_int()) ==
                          wr.read_member_bytes(sid, kv.first).size(),
                  "writer size of " + sid + "/" + kv.first);
        }
    // structure; mulliken_charge; references; input, log, attachment.geom, files
    CHECK(n_members == 7, "every member checked (" + std::to_string(n_members) + ")");

    // A declared size is checked against the central directory.
    auto sized = [](Json size) {
        Json members = Json::object();
        Json m = json_member("structure/structure.json", kWater);
        m.set("size", std::move(size));
        members.set("structure", m);
        return archive(manifest_of({section("structure", "structure", members)}),
                       {{"structure/structure.json", kWater}});
    };
    expect_clean("matching size", sized(Json(static_cast<long long>(kWater.size()))));
    expect_only("size mismatch", sized(Json(static_cast<long long>(kWater.size() + 8))),
                "E_BINARY_SIZE");
    expect_only("negative size", sized(Json(-1)), "E_MEMBER_MALFORMED");
    expect_only("string size", sized(Json("12")), "E_MEMBER_MALFORMED");
    expect_only("fractional size", sized(Json(1.5)), "E_MEMBER_MALFORMED");
    bool threw = false;
    try {
        QvfReader::from_bytes(sized(Json(static_cast<long long>(kWater.size() + 1))))
            .read_json_member("structure", "structure");
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("E_BINARY_SIZE") != std::string::npos;
    }
    CHECK(threw, "reader refuses a member whose size disagrees with its entry");

    // An archive updated in place: the last manifest.json is the manifest,
    // earlier ones are orphans and not duplicates.
    Json stale_prov = Json::object();
    stale_prov.set("run_status", "pending");
    Json prov = Json::object();
    prov.set("run_status", "converged");
    Json pending = Json::object();
    pending.set("provenance", stale_prov);
    Json converged = Json::object();
    converged.set("provenance", prov);
    const std::string stale = manifest_of({structure_section()}, pending).dump();
    std::vector<uint8_t> updated =
        make_zip({{"structure/structure.json", kWater},
                  {"manifest.json", stale, false},
                  {"old.json", "{}"},
                  {"manifest.json", stale, false},
                  {"manifest.json", manifest_of({structure_section()}, converged).dump(),
                   false}});
    std::vector<Finding> fs = semantic_findings(updated);
    size_t superseded = 0, orphans = 0;
    for (const auto& f : fs) {
        orphans += f.code == "W_ORPHAN_MEMBER";
        superseded += f.code == "W_ORPHAN_MEMBER" && f.where == "entry:manifest.json";
    }
    CHECK(!has_errors(fs) && orphans == 3 && superseded == 2 &&
              codes(fs) == (std::set<std::string>{"W_ORPHAN_MEMBER", "W_SCHEMA_SKIPPED"}),
          "superseded manifests are orphans" + describe(fs));
    QvfReader ur = QvfReader::from_bytes(updated);
    CHECK(ur.manifest().at("provenance").at("run_status").as_string() == "converged",
          "the reader uses the last manifest.json");

    // byte_order: only "little"; the schema rejects others (E_SCHEMA, see
    // test_schema_checker; without the schema check the archive is clean), and
    // the reader refuses to guess.
    Json little = Json::object();
    little.set("byte_order", "little");
    Json big = Json::object();
    big.set("byte_order", "big");
    std::vector<uint8_t> le = archive(manifest_of({structure_section()}, little),
                                      {{"structure/structure.json", kWater}});
    std::vector<uint8_t> be = archive(manifest_of({structure_section()}, big),
                                      {{"structure/structure.json", kWater}});
    expect_clean("byte_order little", le);
    expect_clean("byte_order big is schema-only here", be);
    QvfReader::from_bytes(le);
    threw = false;
    try {
        QvfReader::from_bytes(be);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("byte_order") != std::string::npos;
    }
    CHECK(threw, "reader refuses byte_order big");
}

static void test_json_caps() {
    std::string deep = nested(100);
    Json members = Json::object();
    members.set("data", json_member("x/deep.json", deep));
    expect_only("deep member",
                archive(manifest_of({structure_section(), section("deep", "x_test.deep", members)}),
                        {{"structure/structure.json", kWater}, {"x/deep.json", deep}}),
                "E_JSON_TOO_DEEP");

    std::string deep_manifest = manifest_of({structure_section()}).dump(0);
    deep_manifest.insert(deep_manifest.size() - 1, ",\"provenance\":{\"x\":" + nested(70) + "}");
    std::vector<uint8_t> dm = make_zip({{"structure/structure.json", kWater},
                                        {"manifest.json", deep_manifest, false}});
    expect_only("deep manifest", dm, "E_JSON_TOO_DEEP");
    bool threw = false;
    try {
        QvfReader::from_bytes(dm);
    } catch (const JsonDepthError&) {
        threw = true;
    }
    CHECK(threw, "QvfReader refuses a manifest over the depth cap");

    Limits json_cap;
    json_cap.max_json_bytes = 64;
    std::vector<Finding> fs = semantic_findings(base_archive(), json_cap);
    CHECK(has_code(fs, "E_JSON_TOO_LARGE") && fs.size() == 1,
          "manifest over the JSON cap stops validation" + describe(fs));
    json_cap.max_json_bytes = 1024;
    std::string padded = kWater + std::string(2000, ' ');
    expect_only("large member",
                archive(manifest_of({structure_section("structure", "s.json", padded)}),
                        {{"s.json", padded}}),
                "E_JSON_TOO_LARGE", json_cap);
}

static void test_manifest_and_sections() {
    expect_only("missing manifest", make_zip({{"structure/structure.json", kWater}}),
                "E_MANIFEST_MISSING");
    expect_only("manifest not JSON",
                make_zip({{"manifest.json", "{\"qvf_version\": 1,", false}}),
                "E_MANIFEST_JSON");
    expect_only("manifest not an object", make_zip({{"manifest.json", "[1, 2]", false}}),
                "E_MANIFEST_JSON");
    expect_only("manifest bad UTF-8",
                make_zip({{"manifest.json", "{\"x\": \"\xff\"}", false}}),
                "E_MANIFEST_JSON");

    Json v3 = Json::object();
    v3.set("qvf_version", 3);
    std::vector<uint8_t> bad_version =
        archive(manifest_of({structure_section()}, v3), {{"structure/structure.json", kWater}});
    expect_only("version", bad_version, "E_VERSION");

    Json no_src = Json::object();
    no_src.set("source", Json::object());
    expect_only("source",
                archive(manifest_of({structure_section()}, no_src),
                        {{"structure/structure.json", kWater}}),
                "E_SOURCE");

    Json not_array = Json::object();
    not_array.set("sections", Json::object());
    expect_only("sections not array", archive(manifest_of({}, not_array), {}),
                "E_SECTION_MALFORMED");
    expect_only("section not object",
                archive(manifest_of({structure_section(), Json(42)}),
                        {{"structure/structure.json", kWater}}),
                "E_SECTION_MALFORMED");

    expect_only("duplicate id",
                archive(manifest_of({structure_section(),
                                     structure_section("structure", "s2.json")}),
                        {{"structure/structure.json", kWater}, {"s2.json", kWater}}),
                "E_SECTION_ID_DUPLICATE");
    for (const char* id : {"auto_open", "bookmarks", "camera", "crossfade"})
        expect_only(std::string("reserved ") + id,
                    archive(manifest_of({structure_section(id)}),
                            {{"structure/structure.json", kWater}}),
                    "E_RESERVED_ID");

    Json m = Json::object();
    m.set("d", json_member("d.json", "{}"));
    expect_only("unknown kind",
                archive(manifest_of({structure_section(), section("v", "volume.bogus", m)}),
                        {{"structure/structure.json", kWater}, {"d.json", "{}"}}),
                "E_KIND_UNKNOWN");

    Json crit = section("ecp", "x_test.ecp", m);
    crit.set("critical", true);
    std::vector<uint8_t> undeclared =
        archive(manifest_of({structure_section(), crit}),
                {{"structure/structure.json", kWater}, {"d.json", "{}"}});
    expect_only("critical undeclared", undeclared, "E_CRITICAL_UNDECLARED");
    Json ext = Json::object(), ns = Json::object(), root = Json::object();
    ns.set("version", "1");
    ext.set("x_test", ns);
    root.set("extensions", ext);
    std::vector<Finding> fs = semantic_findings(
        archive(manifest_of({structure_section(), crit}, root),
                {{"structure/structure.json", kWater}, {"d.json", "{}"}}));
    CHECK(!has_errors(fs), "declared critical namespace is fine" + describe(fs));
}

static void test_members() {
    Json bad_spec = Json::object();
    bad_spec.set("format", "json");
    Json members = Json::object();
    members.set("structure", json_member("structure/structure.json", kWater))
        .set("bonds", bad_spec);
    expect_only("malformed member",
                archive(manifest_of({section("structure", "structure", members)}),
                        {{"structure/structure.json", kWater}}),
                "E_MEMBER_MALFORMED");

    expect_only("missing member", archive(manifest_of({structure_section()}), {}),
                "E_MEMBER_MISSING");

    Json wrong = Json::object();
    wrong.set("path", "structure/structure.json").set("format", "json")
        .set("sha256", sha("other"));
    Json wm = Json::object();
    wm.set("structure", wrong);
    expect_only("sha256",
                archive(manifest_of({section("structure", "structure", wm)}),
                        {{"structure/structure.json", kWater}}),
                "E_SHA256_MISMATCH");

    std::string broken = "{\"atoms\": [";
    expect_only("invalid JSON",
                archive(manifest_of({structure_section("structure", "s.json", broken)}),
                        {{"s.json", broken}}),
                "E_JSON_INVALID");

    std::string charges(16, '\0');
    Json ap = Json::object();
    ap.set("mulliken_charge", binary_member("q.bin", charges, "float64", 3));
    expect_only("binary size",
                archive(manifest_of({structure_section(), section("ap", "atom_properties", ap)}),
                        {{"structure/structure.json", kWater}, {"q.bin", charges}}),
                "E_BINARY_SIZE");
    // Fuzzing regression (WP-50): an integer literal beyond int64 in a shape
    // parses as a double. Json::as_int() converted it unchecked (undefined
    // behaviour); read_binary_member now refuses it and validation reports it.
    {
        Json huge = binary_member("q.bin", charges, "float64", 2);
        Json shape = Json::array();
        shape.push_back(Json::raw("1" + std::string(98, '0')));
        huge.set("shape", shape);
        Json aph = Json::object();
        aph.set("mulliken_charge", huge);
        const std::vector<uint8_t> a =
            archive(manifest_of({structure_section(), section("ap", "atom_properties", aph)}),
                    {{"structure/structure.json", kWater}, {"q.bin", charges}});
        CHECK(throws_runtime([&] {
                  QvfReader::from_bytes(a).read_binary_member("ap", "mulliken_charge");
              }),
              "shape beyond int64: read_binary_member throws");
        CHECK(throws_runtime([] { (void)Json::parse("1e300").as_int(); }) &&
                  throws_runtime([] { (void)Json::parse("-1e19").as_int(); }) &&
                  Json::parse("-9.2e18").as_int() == -9200000000000000000LL,
              "as_int refuses doubles outside the int64 range");
        bool err = false;
        try {
            err = has_errors(QvfReader::validate_bytes(a));
        } catch (const std::exception& e) {
            std::cerr << "    validate_bytes threw: " << e.what() << "\n";
        }
        CHECK(err, "shape beyond int64: validation reports an error and does not throw");
    }
    Json ap16 = Json::object();
    ap16.set("mulliken_charge", binary_member("q.bin", charges, "float16", 8));
    expect_only("dtype",
                archive(manifest_of({structure_section(), section("ap", "atom_properties", ap16)}),
                        {{"structure/structure.json", kWater}, {"q.bin", charges}}),
                "E_DTYPE");

    // Custom atom_properties roles need a `meta` entry (spec § 5.7).
    {
        std::string two(16, '\0');
        std::string meta_ok =
            R"({"custom.a": {"units": "e", "description": "d", "kind": "scalar"}})";
        std::string meta_other =
            R"({"custom.b": {"units": "e", "description": "d", "kind": "scalar"}})";
        std::string meta_scalar = R"({"custom.a": "scalar"})";
        auto custom_archive = [&](const std::string* meta) {
            Json cm = Json::object();
            cm.set("custom.a", binary_member("ap/custom.a.bin", two, "float64", 2));
            std::vector<Entry> files = {{"structure/structure.json", kWater},
                                        {"ap/custom.a.bin", two}};
            if (meta) {
                cm.set("meta", json_member("ap/meta.json", *meta));
                files.push_back({"ap/meta.json", *meta});
            }
            return archive(manifest_of({structure_section(),
                                        section("ap", "atom_properties", cm)}),
                           files);
        };
        expect_clean("custom role with a meta entry", custom_archive(&meta_ok));
        expect_only("custom role without a meta member", custom_archive(nullptr),
                    "E_MEMBER_MALFORMED");
        expect_only("meta without the role's entry", custom_archive(&meta_other),
                    "E_MEMBER_MALFORMED");
        expect_only("meta entry that is not an object", custom_archive(&meta_scalar),
                    "E_MEMBER_MALFORMED");
        std::string meta_array = "[]";
        expect_only("meta that is not an object", custom_archive(&meta_array),
                    "E_MEMBER_MALFORMED");
        std::string meta_broken = "{not json";
        std::vector<Finding> mf = semantic_findings(custom_archive(&meta_broken));
        CHECK(has_code(mf, "E_JSON_INVALID") && !has_code(mf, "E_MEMBER_MALFORMED"),
              "an unreadable meta member is not judged twice" + describe(mf));
        Json listed = Json::object();
        listed.set("mulliken_charge", binary_member("ap/q.bin", two, "float64", 2));
        expect_clean("a listed role needs no meta",
                     archive(manifest_of({structure_section(),
                                          section("ap", "atom_properties", listed)}),
                             {{"structure/structure.json", kWater}, {"ap/q.bin", two}}));
        // The finding points at the custom member.
        std::vector<Finding> pf = semantic_findings(custom_archive(nullptr));
        bool pointed = false;
        for (const auto& f : pf)
            pointed = pointed || (f.code == "E_MEMBER_MALFORMED" &&
                                  f.where == "/sections/1/members/custom.a");
        CHECK(pointed, "finding locates the custom member" + describe(pf));
    }

    std::string log = "\xff\xfe not UTF-8\n";
    Json rr_members = Json::object();
    Json opaque = Json::object();
    opaque.set("path", "log.txt").set("format", "binary").set("sha256", sha(log));
    rr_members.set("log", opaque);
    Json rr = section("rr", "run.record", rr_members);
    rr.set("program", "t");
    expect_only("run.record UTF-8",
                archive(manifest_of({structure_section(), rr}),
                        {{"structure/structure.json", kWater}, {"log.txt", log}}),
                "E_UTF8");

    std::string ok_log = "done\n";
    std::string files = R"({"attachment.restart": {"filename": "r.chk"}})";
    Json fm = Json::object(), op2 = Json::object();
    op2.set("path", "log.txt").set("format", "binary").set("sha256", sha(ok_log));
    fm.set("log", op2).set("files", json_member("files.json", files));
    Json rr2 = section("rr", "run.record", fm);
    rr2.set("program", "t");
    expect_only("files index",
                archive(manifest_of({structure_section(), rr2}),
                        {{"structure/structure.json", kWater}, {"log.txt", ok_log},
                         {"files.json", files}}),
                "E_FILES_INDEX");

    std::string wp = R"({"waypoints": []})";
    Json wpm = Json::object();
    wpm.set("waypoints", json_member("wp.json", wp));
    Json wps = section("wp", "reaction.waypoints", wpm);
    wps.set("trajectory_ref", "nope");
    expect_only("trajectory_ref",
                archive(manifest_of({structure_section(), wps}),
                        {{"structure/structure.json", kWater}, {"wp.json", wp}}),
                "E_REF_UNRESOLVED");

    std::string grid = "{}";
    std::string vol(8, '\0');
    Json vm = Json::object(), vshape = Json::array();
    vshape.push_back(1).push_back(1).push_back(2);
    Json vdata = Json::object();
    vdata.set("path", "v.bin").set("format", "binary").set("dtype", "float32")
        .set("shape", vshape).set("sha256", sha(vol));
    vm.set("grid", json_member("g.json", grid)).set("data", vdata);
    Json diff = section("diff", "volume.difference", vm);
    diff.set("operand_a", "structure").set("operand_b", "missing");
    expect_only("operand",
                archive(manifest_of({structure_section(), diff}),
                        {{"structure/structure.json", kWater}, {"g.json", grid},
                         {"v.bin", vol}}),
                "E_REF_UNRESOLVED");

    std::string periodic =
        R"({"atoms":[],"lattice_vectors":[[5,0,0],[0,5,0],[0,0,5]]})";
    expect_only("periodicity: lattice without pbc",
                archive(manifest_of({structure_section("structure", "s.json", periodic)}),
                        {{"s.json", periodic}}),
                "E_STRUCTURE_PERIODICITY");
    std::string dim_mismatch =
        R"({"atoms":[],"pbc":[true,true,false],"lattice_vectors":[[5,0,0],[0,5,0],[0,0,5]],"dimensionality":3})";
    expect_only("periodicity: dimensionality",
                archive(manifest_of({structure_section("structure", "s.json", dim_mismatch)}),
                        {{"s.json", dim_mismatch}}),
                "E_STRUCTURE_PERIODICITY");
}

// A `bands` section whose eigenvalues / projections / occupations have the
// given shapes (zero-filled float64 payloads of the right byte length).
static std::vector<uint8_t> bands_archive(const std::vector<int64_t>& eig,
                                          const std::vector<int64_t>& proj,
                                          const std::vector<int64_t>& occ = {}) {
    auto member = [](const std::string& path, const std::vector<int64_t>& shape,
                     std::string& data) {
        size_t n = 1;
        Json sh = Json::array();
        for (int64_t d : shape) {
            n *= static_cast<size_t>(d);
            sh.push_back(Json(static_cast<long long>(d)));
        }
        data.assign(8 * n, '\0');
        Json m = Json::object();
        m.set("path", path).set("format", "binary").set("dtype", "float64")
            .set("shape", sh).set("sha256", sha(data));
        return m;
    };
    std::string kpath = R"({"n_bands": 2})", e, p, o;
    Json members = Json::object();
    members.set("kpath", json_member("b/kpath.json", kpath));
    members.set("eigenvalues", member("b/e.bin", eig, e));
    std::vector<Entry> files = {{"structure/structure.json", kWater},
                                {"b/kpath.json", kpath}, {"b/e.bin", e}};
    if (!proj.empty()) {
        members.set("projections", member("b/p.bin", proj, p));
        files.push_back({"b/p.bin", p});
    }
    if (!occ.empty()) {
        members.set("occupations", member("b/o.bin", occ, o));
        files.push_back({"b/o.bin", o});
    }
    return archive(manifest_of({structure_section(), section("bands", "bands", members)}),
                   files);
}

static void test_bands_spin_axis() {
    expect_clean("bands: rank-4 projections + occupations agree",
                 bands_archive({2, 3, 2}, {2, 3, 2, 4}, {2, 3, 2}));
    expect_clean("bands: legacy rank-3 projections", bands_archive({1, 3, 2}, {3, 2, 4}));
    expect_only("bands: rank-4 projections disagree on n_spin",
                bands_archive({1, 3, 2}, {2, 3, 2, 4}), "E_MEMBER_MALFORMED");
    expect_only("bands: rank-4 projections disagree on n_bands",
                bands_archive({2, 3, 2}, {2, 3, 5, 4}), "E_MEMBER_MALFORMED");
    expect_only("bands: rank-3 projections disagree on n_k",
                bands_archive({1, 3, 2}, {4, 2, 4}), "E_MEMBER_MALFORMED");
    expect_only("bands: occupations disagree with eigenvalues",
                bands_archive({2, 3, 2}, {}, {1, 3, 2}), "E_MEMBER_MALFORMED");
}

// --- unit names (spec § 4.3, "Declared units") --------------------------------

// The base archive plus `provenance` and, when `kind` is set, one section of
// that kind whose `role` member is the JSON text `payload`.
static std::vector<uint8_t> units_archive(const std::string& provenance,
                                          const std::string& kind = "",
                                          const std::string& role = "",
                                          const std::string& payload = "") {
    Json extra = Json::object();
    if (!provenance.empty()) extra.set("provenance", Json::parse(provenance));
    std::vector<Json> sections{structure_section()};
    std::vector<Entry> files{{"structure/structure.json", kWater}};
    if (!kind.empty()) {
        Json members = Json::object();
        members.set(role, json_member("p/payload.json", payload));
        sections.push_back(section("p", kind, members));
        files.push_back({"p/payload.json", payload});
    }
    return archive(manifest_of(sections, extra), files);
}

static void test_unit_names() {
    const auto warn = Finding::Severity::Warning;
    expect_clean("provenance legacy spelling",
                 units_archive(R"({"scf_energy":{"value":-1,"units":"Eh"},)"
                               R"("fermi_energy":{"value":0,"units":"eV"}})"));
    expect_clean("provenance unit name",
                 units_archive(R"({"scf_energy":{"value":-1,"units":"hartree"}})"));
    expect_only("provenance unknown unit",
                units_archive(R"({"scf_energy":{"value":-1,"units":"kJ/mol"}})"),
                "W_UNIT_UNKNOWN", Limits{}, warn);
    expect_only("grid units: unknown name",
                units_archive("", "volume.density", "grid",
                              R"({"origin":[0,0,0],"voxel_vectors":[[1,0,0],[0,1,0],[0,0,1]],)"
                              R"("shape":[1,1,1],"units":{"data":"e_per_bohr3"}})"),
                "W_UNIT_UNKNOWN", Limits{}, warn);
    expect_only("spectrum x_units", units_archive("", "spectra.uvvis", "spectrum",
                                                  R"({"frequencies":[1],"intensities":[1],)"
                                                  R"("x_units":"1/cm"})"),
                "W_UNIT_UNKNOWN", Limits{}, warn);
    expect_clean("spectrum unit names",
                 units_archive("", "spectra.uvvis", "spectrum",
                               R"({"frequencies":[1],"intensities":[1],"x_units":"nm",)"
                               R"("broadening":{"type":"none","units":"ev"}})"));
    expect_clean("vendor units are free",
                 units_archive("", "x_test.p", "m", R"({"units":{"a":"furlong"}})"));
    std::vector<Finding> fs = semantic_findings(
        units_archive("", "spectra.ir", "spectrum",
                      R"({"frequencies":[1],"intensities":[1],"x_units":"a","y_units":"b",)"
                      R"("broadening":{"type":"none","units":"c"},"units":{"x":7}})"));
    size_t n = 0;
    for (const auto& f : fs)
        if (f.code == "W_UNIT_UNKNOWN") {
            ++n;
            CHECK(f.where == "/sections/1/members/spectrum", "unit finding points at the member");
        }
    CHECK(n == 4, "one W_UNIT_UNKNOWN per unknown declared unit");
}

// --- forces / hessian / root energies (spec § 3.6, § 5.10, § 5.11) ----------

// A zero-filled float64 member of the given shape; its bytes go to `files`.
static Json f64_member(const std::string& path, const std::vector<int64_t>& shape,
                       std::vector<Entry>& files) {
    size_t n = 1;
    Json sh = Json::array();
    for (int64_t d : shape) {
        n *= static_cast<size_t>(d);
        sh.push_back(Json(static_cast<long long>(d)));
    }
    std::string data(8 * n, '\0');
    files.push_back({path, data});
    Json m = Json::object();
    m.set("path", path).set("format", "binary").set("dtype", "float64")
        .set("shape", sh).set("sha256", sha(data));
    return m;
}

// The base water structure (2 atoms in kWater) plus one forces or hessian
// section with the given member shape and optional meta / section fields.
static std::vector<uint8_t> derivative_archive(const std::string& kind,
                                               const std::vector<int64_t>& shape,
                                               const std::string& meta = "",
                                               Json fields = Json::object(),
                                               const std::vector<int64_t>& stress = {}) {
    std::vector<Entry> files = {{"structure/structure.json", kWater}};
    Json members = Json::object();
    members.set(kind == "forces" ? "gradient" : "matrix",
                f64_member("d/main.bin", shape, files));
    if (!stress.empty()) members.set("stress", f64_member("d/stress.bin", stress, files));
    if (!meta.empty()) {
        members.set("meta", json_member("d/meta.json", meta));
        files.push_back({"d/meta.json", meta});
    }
    Json sec = section(kind, kind, members);
    for (const auto& kv : fields.items()) sec.set(kv.first, kv.second);
    return archive(manifest_of({structure_section(), sec}), files);
}

// --- volumes (spec § 5.2) --------------------------------------------------------

static const std::string kCubic =
    R"({"atoms":[{"symbol":"Na","position":[0.0,0.0,0.0],"atomic_number":11}],)"
    R"("pbc":[true,true,true],)"
    R"("lattice_vectors":[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]})";

// A structure (`structure`) plus one volume section of `kind` whose grid is
// `grid` and whose float32 data (zeros) has `shape`, with an optional axis
// member and section fields.
static std::vector<uint8_t> volume_archive(const std::string& kind, const std::string& grid,
                                           const std::vector<int64_t>& shape,
                                           const std::string& axis = "",
                                           const std::string& structure = kWater,
                                           Json fields = Json::object()) {
    std::vector<Entry> files = {{"structure/structure.json", structure},
                                {"v/grid.json", grid}};
    size_t n = 1;
    Json sh = Json::array();
    for (int64_t d : shape) {
        n *= static_cast<size_t>(d);
        sh.push_back(Json(static_cast<long long>(d)));
    }
    const std::string data(4 * n, '\0');
    files.push_back({"v/data.bin", data});
    Json dm = Json::object();
    dm.set("path", "v/data.bin").set("format", "binary").set("dtype", "float32")
        .set("shape", sh).set("sha256", sha(data));
    Json members = Json::object();
    members.set("grid", json_member("v/grid.json", grid)).set("data", dm);
    if (!axis.empty()) {
        members.set("axis", json_member("v/axis.json", axis));
        files.push_back({"v/axis.json", axis});
    }
    Json sec = section("v", kind, members);
    for (const auto& kv : fields.items()) sec.set(kv.first, kv.second);
    return archive(manifest_of({structure_section("structure", "structure/structure.json",
                                                  structure),
                                sec}),
                   files);
}

static std::string unit_grid(const std::string& shape) {
    return R"({"origin":[0,0,0],"voxel_vectors":[[1,0,0],[0,1,0],[0,0,1]],"shape":)" +
           shape + "}";
}

// A grid periodic on all three axes of the 4 Å cubic cell, shape [4, 4, 4],
// with `step` bohr per point on every axis.
static std::string periodic_grid(const std::string& step, const std::string& extra) {
    return R"({"origin":[0,0,0],"voxel_vectors":[[)" + step + ",0,0],[0," + step +
           ",0],[0,0," + step + R"(]],"shape":[4,4,4],"periodic":[true,true,true])" +
           extra + "}";
}

static void test_volumes() {
    const std::string g3 = unit_grid("[2,2,2]");
    expect_clean("volume.vector [2, 2, 2, 3]", volume_archive("volume.vector", g3, {2, 2, 2, 3}));
    expect_clean("volume.complex [2, 2, 2, 2]",
                 volume_archive("volume.complex", g3, {2, 2, 2, 2}));
    expect_clean("volume.laplacian [2, 2, 2]", volume_archive("volume.laplacian", g3, {2, 2, 2}));
    expect_clean("volume.series [3, 2, 2, 2]",
                 volume_archive("volume.series", g3, {3, 2, 2, 2},
                                R"({"kind":"time","values":[0,1,2],"units":"fs"})"));
    expect_only("volume.vector trailing axis 2",
                volume_archive("volume.vector", g3, {2, 2, 2, 2}), "E_MEMBER_MALFORMED");
    expect_only("volume.complex trailing axis 3",
                volume_archive("volume.complex", g3, {2, 2, 2, 3}), "E_MEMBER_MALFORMED");
    expect_only("volume.laplacian rank 4",
                volume_archive("volume.laplacian", g3, {2, 2, 2, 1}), "E_MEMBER_MALFORMED");
    expect_only("volume.series rank 3",
                volume_archive("volume.series", g3, {2, 2, 2}, R"({"kind":"custom","values":[0,1]})"),
                "E_MEMBER_MALFORMED");
    expect_only("volume.series: two values for three slabs",
                volume_archive("volume.series", g3, {3, 2, 2, 2},
                               R"({"kind":"time","values":[0,1]})"),
                "E_MEMBER_MALFORMED");
    expect_only("volume.series: two labels for three slabs",
                volume_archive("volume.series", g3, {3, 2, 2, 2},
                               R"({"kind":"state","values":[1,2,3],"labels":["S1","S2"]})"),
                "E_MEMBER_MALFORMED");
    expect_only("volume.series: trajectory_ref names a structure",
                volume_archive("volume.series", g3, {3, 2, 2, 2},
                               R"({"kind":"frame","values":[0,1,2],"trajectory_ref":"structure"})"),
                "E_REF_UNRESOLVED");
    expect_only("volume.series: an unknown axis unit is a warning",
                volume_archive("volume.series", g3, {3, 2, 2, 2},
                               R"({"kind":"time","values":[0,1,2],"units":"sec"})"),
                "W_UNIT_UNKNOWN", Limits{}, Finding::Severity::Warning);

    // Periodic grids against the 4 Å cubic cell, 7.5589045 bohr.
    const std::string chgcar = "1.8897261246";          // 4 points per period
    const std::string cube = "2.5196348328";            // 3 steps per period
    expect_clean("periodic: CHGCAR style spans the cell",
                 volume_archive("volume.density",
                                periodic_grid(chgcar, R"(,"includes_endpoint":false)"),
                                {4, 4, 4}, "", kCubic));
    expect_clean("periodic: cube style spans the cell",
                 volume_archive("volume.density",
                                periodic_grid(cube, R"(,"includes_endpoint":true)"),
                                {4, 4, 4}, "", kCubic));
    expect_clean("periodic: includes_endpoint defaults to true",
                 volume_archive("volume.density", periodic_grid(cube, ""), {4, 4, 4}, "",
                                kCubic));
    expect_only("periodic: CHGCAR numbers declared cube style",
                volume_archive("volume.density",
                               periodic_grid(chgcar, R"(,"includes_endpoint":true)"),
                               {4, 4, 4}, "", kCubic),
                "E_MEMBER_MALFORMED");
    expect_only("periodic: a grid that misses the cell",
                volume_archive("volume.density",
                               periodic_grid("1.9", R"(,"includes_endpoint":false)"),
                               {4, 4, 4}, "", kCubic),
                "E_MEMBER_MALFORMED");
    expect_only("periodic: a periodic grid on a molecule",
                volume_archive("volume.density", periodic_grid(chgcar, ""), {4, 4, 4}),
                "E_MEMBER_MALFORMED");
    expect_only("grid structure_ref names no structure",
                volume_archive("volume.density",
                               periodic_grid(chgcar, R"(,"structure_ref":"v")"), {4, 4, 4},
                               "", kCubic),
                "E_REF_UNRESOLVED");
    {   // The section-level structure_ref wins over the grid's.
        std::vector<Entry> files = {{"structure/structure.json", kCubic},
                                    {"mol/structure.json", kWater}};
        const std::string grid =
            periodic_grid(chgcar, R"(,"includes_endpoint":false,"structure_ref":"mol")");
        files.push_back({"v/grid.json", grid});
        const std::string data(4 * 64, '\0');
        files.push_back({"v/data.bin", data});
        Json dm = Json::object();
        Json sh = Json::array();
        for (int k = 0; k < 3; ++k) sh.push_back(Json(4));
        dm.set("path", "v/data.bin").set("format", "binary").set("dtype", "float32")
            .set("shape", sh).set("sha256", sha(data));
        Json members = Json::object();
        members.set("grid", json_member("v/grid.json", grid)).set("data", dm);
        Json vol = section("v", "volume.density", members);
        const auto bytes = [&](const char* ref) {
            Json sec = vol;
            if (ref) sec.set("structure_ref", ref);
            return archive(manifest_of({structure_section("structure",
                                                          "structure/structure.json", kCubic),
                                        structure_section("mol", "mol/structure.json", kWater),
                                        sec}),
                           files);
        };
        expect_only("grid structure_ref names a molecule", bytes(nullptr),
                    "E_MEMBER_MALFORMED");
        expect_clean("the section-level structure_ref takes precedence", bytes("structure"));
    }
}

static std::vector<uint8_t> energies_archive(const std::string& energies,
                                             const std::string& provenance,
                                             const std::string& thermo = "") {
    Json extra = Json::object();
    extra.set("energies", Json::parse(energies));
    if (!provenance.empty()) extra.set("provenance", Json::parse(provenance));
    if (!thermo.empty()) extra.set("thermochemistry", Json::parse(thermo));
    return archive(manifest_of({structure_section()}, extra),
                   {{"structure/structure.json", kWater}});
}

static void test_derivatives_and_energies() {
    expect_clean("forces: [2, 3] gradient + stress for a 2-atom structure",
                 derivative_archive("forces", {2, 3}, R"({"structure_ref":"structure"})",
                                    Json::object(), {3, 3}));
    expect_only("forces: gradient rows differ from the structure's atoms",
                derivative_archive("forces", {3, 3}, R"({"structure_ref":"structure"})"),
                "E_MEMBER_MALFORMED");
    expect_only("forces: the only structure is used when none is named",
                derivative_archive("forces", {3, 3}), "E_MEMBER_MALFORMED");
    expect_only("forces: [2, 2] gradient", derivative_archive("forces", {2, 2}),
                "E_MEMBER_MALFORMED");
    expect_only("forces: [3, 2] stress",
                derivative_archive("forces", {2, 3}, "", Json::object(), {3, 2}),
                "E_MEMBER_MALFORMED");
    expect_only("forces: meta.structure_ref names no structure",
                derivative_archive("forces", {2, 3}, R"({"structure_ref":"forces"})"),
                "E_REF_UNRESOLVED");
    Json own_ref = Json::object();
    own_ref.set("structure_ref", "elsewhere");
    // Since the generic section references landed (spec § 4.1), the unresolved
    // section-level structure_ref is itself reported; the atom count is still
    // not checked against it.
    expect_only("forces: an unresolvable section-level structure_ref is reported, "
                "the count skipped",
                derivative_archive("forces", {5, 3}, R"({"structure_ref":"structure"})",
                                   own_ref),
                "E_REF_UNRESOLVED");

    expect_clean("hessian: [6, 6] for a 2-atom structure",
                 derivative_archive("hessian", {6, 6}, R"({"atom_masses_amu":[16.0,1.0]})"));
    expect_only("hessian: non-square matrix", derivative_archive("hessian", {6, 3}),
                "E_MEMBER_MALFORMED");
    expect_only("hessian: extent not divisible by 3", derivative_archive("hessian", {4, 4}),
                "E_MEMBER_MALFORMED");
    expect_only("hessian: [9, 9] for a 2-atom structure",
                derivative_archive("hessian", {9, 9}), "E_MEMBER_MALFORMED");
    expect_only("hessian: masses not one per atom",
                derivative_archive("hessian", {6, 6}, R"({"atom_masses_amu":[16.0]})"),
                "E_MEMBER_MALFORMED");

    {   // vibrations.hessian_ref must name a hessian section
        std::string meta = R"({"frequencies":[1.0]})";
        std::vector<Entry> files = {{"structure/structure.json", kWater},
                                    {"v/metadata.json", meta}};
        Json members = Json::object();
        members.set("metadata", json_member("v/metadata.json", meta));
        members.set("displacements", f64_member("v/d.bin", {1, 2, 3}, files));
        Json vib = section("vibrations", "vibrations", members);
        vib.set("hessian_ref", "structure");
        expect_only("vibrations: hessian_ref names a structure",
                    archive(manifest_of({structure_section(), vib}), files),
                    "E_REF_UNRESOLVED");
    }

    const std::string rhf = R"({"method":"RHF","scf_energy":{"value":-74.9659011,"units":"hartree"}})";
    expect_only("energies: total for the provenance method disagrees",
                energies_archive(R"({"total":{"value":-74.96,"units":"hartree"},"method":"rhf"})",
                                 rhf),
                "E_ROOT_INCONSISTENT");
    expect_only("energies: a by_method entry disagrees",
                energies_archive(R"({"total":{"value":-75.2,"units":"hartree"},"method":"mp2",)"
                                 R"("by_method":[{"method":"rhf","total":{"value":-74.95,)"
                                 R"("units":"hartree"}}]})",
                                 rhf),
                "E_ROOT_INCONSISTENT");
    expect_clean("energies: the same energy in eV agrees",
                 energies_archive(R"({"total":{"value":-2039.9260901106366,"units":"ev"},)"
                                  R"("method":"rhf"})",
                                  rhf));
    expect_clean("energies: a total of another method is not compared",
                 energies_archive(R"({"total":{"value":-75.2,"units":"hartree"},"method":"mp2"})",
                                  rhf));
    expect_clean("energies: an unrecognised unit is not compared",
                 energies_archive(R"({"total":{"value":-75.2,"units":"furlong"},"method":"rhf"})",
                                  rhf));
    expect_only("energies: zpve disagrees with thermochemistry",
                energies_archive(R"({"total":{"value":-75.0,"units":"hartree"},)"
                                 R"("components":{"zpve":{"value":0.0213,"units":"hartree"}}})",
                                 "", R"({"zpve_eh":0.0214})"),
                "E_ROOT_INCONSISTENT");
    expect_clean("energies: zpve agrees with thermochemistry",
                 energies_archive(R"({"total":{"value":-75.0,"units":"hartree"},)"
                                  R"("components":{"zpve":{"value":0.0213,"units":"hartree"}}})",
                                  "", R"({"zpve_eh":0.0213})"));

}

// --- section references (spec § 4.1, § 5.5) ----------------------------------

// structure + run.record "run" + trajectory "traj" + scf_history "scf" (with
// `scf_fields` merged onto its section object) + reaction.waypoints "wp".
static std::vector<uint8_t> refs_archive(const Json& scf_fields,
                                         const std::string& trajectory_ref = "traj") {
    const std::string log = "done\n";
    Json op = Json::object();
    op.set("path", "run/log.txt").set("format", "binary").set("sha256", sha(log));
    Json rm = Json::object();
    rm.set("log", op);
    Json run = section("run", "run.record", rm);
    run.set("program", "t");

    const std::string meta = R"({"symbols":["O","H"]})";
    const std::string coords(48, '\0');                       // float64 [1, 2, 3]
    Json tm = Json::object();
    tm.set("metadata", json_member("traj/metadata.json", meta))
        .set("coords", binary_member("traj/coords.bin", coords, "float64", 6));
    Json traj = section("traj", "trajectory", tm);

    const std::string iters = R"({"iterations":[{"iter":1,"energy_eh":-74.9}]})";
    Json sm = Json::object();
    sm.set("iterations", json_member("scf.json", iters));
    Json scf = section("scf", "scf_history", sm);
    for (const auto& kv : scf_fields.items()) scf.set(kv.first, kv.second);

    const std::string wp =
        R"({"waypoints":[{"frame_index":0,"label":"R","kind":"reactant"}]})";
    Json wpm = Json::object();
    wpm.set("waypoints", json_member("wp.json", wp));
    Json wps = section("wp", "reaction.waypoints", wpm);
    wps.set("trajectory_ref", trajectory_ref);

    return archive(manifest_of({structure_section(), run, traj, scf, wps}),
                   {{"structure/structure.json", kWater}, {"run/log.txt", log},
                    {"traj/metadata.json", meta}, {"traj/coords.bin", coords},
                    {"scf.json", iters}, {"wp.json", wp}});
}

static Json one_field(const std::string& key, Json value) {
    Json o = Json::object();
    o.set(key, std::move(value));
    return o;
}

static void test_section_references() {
    Json all = Json::object();
    all.set("structure_ref", "structure").set("run_ref", "run").set("group", "step 1");
    expect_clean("refs: structure_ref, run_ref and group resolve", refs_archive(all));
    expect_clean("refs: a null structure_ref is absent",
                 refs_archive(one_field("structure_ref", Json(nullptr))));
    expect_only("refs: structure_ref names a run.record",
                refs_archive(one_field("structure_ref", "run")), "E_REF_UNRESOLVED");
    expect_only("refs: structure_ref names a trajectory",
                refs_archive(one_field("structure_ref", "traj")), "E_REF_UNRESOLVED");
    expect_only("refs: structure_ref names no section",
                refs_archive(one_field("structure_ref", "nowhere")), "E_REF_UNRESOLVED");
    expect_only("refs: run_ref names a structure",
                refs_archive(one_field("run_ref", "structure")), "E_REF_UNRESOLVED");
    expect_only("refs: run_ref is not a string",
                refs_archive(one_field("run_ref", Json(5))), "E_REF_UNRESOLVED");
    expect_only("refs: trajectory_ref names a structure",
                refs_archive(Json::object(), "structure"), "E_REF_UNRESOLVED");
    expect_only("refs: trajectory_ref names an scf_history",
                refs_archive(Json::object(), "scf"), "E_REF_UNRESOLVED");
}

// --- wavefunction.gto (spec § 5.3, Appendix A.1) -----------------------------

// A float64 member holding `values` (little-endian); its bytes go to `files`.
static Json f64_values(const std::string& path, const std::vector<int64_t>& shape,
                       const std::vector<double>& values, std::vector<Entry>& files) {
    const Tensor t = tensor_f64(values, shape);
    const std::string data(t.bytes.begin(), t.bytes.end());
    files.push_back({path, data});
    Json sh = Json::array();
    for (int64_t d : shape) sh.push_back(Json(static_cast<long long>(d)));
    Json m = Json::object();
    m.set("path", path).set("format", "binary").set("dtype", "float64")
        .set("shape", sh).set("sha256", sha(data));
    return m;
}

struct WfMember {
    std::string role;
    std::vector<int64_t> shape;
    std::vector<double> values;   // empty: zeros
};

// structure + one wavefunction.gto section "wf" with these payloads and
// binary members; `extra` is merged into the manifest root.
static std::vector<uint8_t> wf_archive(const std::string& basis, const std::string& meta,
                                       const std::vector<WfMember>& binaries,
                                       Json extra = Json::object()) {
    std::vector<Entry> files = {{"structure/structure.json", kWater},
                                {"wf/basis.json", basis},
                                {"wf/mo_metadata.json", meta}};
    Json members = Json::object();
    members.set("basis", json_member("wf/basis.json", basis));
    members.set("mo_metadata", json_member("wf/mo_metadata.json", meta));
    for (const auto& b : binaries) {
        size_t n = 1;
        for (int64_t d : b.shape) n *= static_cast<size_t>(d);
        members.set(b.role, f64_values("wf/" + b.role + ".dat", b.shape,
                                       b.values.empty() ? std::vector<double>(n, 0.0)
                                                        : b.values,
                                       files));
    }
    return archive(manifest_of({structure_section(), section("wf", "wavefunction.gto",
                                                             members)},
                               extra),
                   files);
}

static std::vector<double> diagonal(const std::vector<double>& d) {
    std::vector<double> m(d.size() * d.size(), 0.0);
    for (size_t i = 0; i < d.size(); ++i) m[i * d.size() + i] = d[i];
    return m;
}

static void test_wavefunction_rules() {
    // An s and a Cartesian d shell, one primitive each with coefficient 1: the
    // Appendix A.1 diagonal is 1, 1, 1/3, 1/3, 1, 1/3, 1.
    const std::string cart =
        R"({"structure_ref":"structure","pure":false,"n_ao":7,"shells":[)"
        R"({"center":0,"l":0,"exponents":[1.2],"coefficients":[1.0]},)"
        R"({"center":0,"l":2,"exponents":[0.8],"coefficients":[1.0]}]})";
    const std::string meta7 =
        R"({"spin":"restricted","energies":[0,0,0,0,0,0,0],"occupations":[2,0,0,0,0,0,0]})";
    const double t = 1.0 / 3.0;
    const std::vector<double> diag7 = {1, 1, t, t, 1, t, 1};
    const WfMember c7{"mo_coefficients", {7, 7}, diagonal(std::vector<double>(7, 1.0))};
    expect_clean("wavefunction: the Appendix A.1 diagonal of a Cartesian d shell",
                 wf_archive(cart, meta7, {c7, {"overlap", {7, 7}, diagonal(diag7)}}));
    std::vector<double> bad = diag7;
    bad[2] = 1.0;   // d_xy normalized on its own: the per-component bug
    expect_only("wavefunction: d_xy at 1 instead of 1/3",
                wf_archive(cart, meta7, {c7, {"overlap", {7, 7}, diagonal(bad)}}),
                "E_NORMALIZATION");
    {
        std::vector<Finding> fs = QvfReader::validate_bytes(
            wf_archive(cart, meta7, {c7, {"overlap", {7, 7}, diagonal(bad)}}));
        for (const auto& f : fs)
            if (f.code == "E_NORMALIZATION")
                CHECK(f.where == "/sections/1/members/overlap" &&
                          f.message.find("AO 2") != std::string::npos,
                      "E_NORMALIZATION points at the overlap and names AO 2");
    }
    // A spherical d shell: every component has the axial self-overlap.
    const std::string sph =
        R"({"structure_ref":"structure","shells":[)"
        R"({"center":0,"l":2,"exponents":[0.8],"coefficients":[1.0]}]})";
    const std::string meta5 = R"({"spin":"restricted","energies":[0,0,0,0,0],)"
                              R"("occupations":[0,0,0,0,0]})";
    const WfMember c5{"mo_coefficients", {5, 5}, {}};
    expect_clean("wavefunction: a spherical d shell has a diagonal of 1",
                 wf_archive(sph, meta5, {c5, {"overlap", {5, 5}, diagonal({1, 1, 1, 1, 1})}}));
    expect_only("wavefunction: the Cartesian factor applied to a spherical shell",
                wf_archive(sph, meta5, {c5, {"overlap", {5, 5}, diagonal({t, 1, 1, 1, 1})}}),
                "E_NORMALIZATION");
    // A contracted s shell that is not normalized: the diagonal is S_axial.
    const std::string s2 =
        R"({"structure_ref":"structure","shells":[)"
        R"({"center":0,"l":0,"exponents":[2.0,0.5],"coefficients":[0.6,0.5]}]})";
    const std::string meta1 = R"({"spin":"restricted","energies":[0],"occupations":[2]})";
    const double s01 = std::pow(2.0 * std::sqrt(2.0 * 0.5) / 2.5, 1.5);
    const double s_axial = 0.36 + 0.25 + 2 * 0.6 * 0.5 * s01;
    const WfMember c1{"mo_coefficients", {1, 1}, {1.0}};
    expect_clean("wavefunction: an unnormalized contraction",
                 wf_archive(s2, meta1, {c1, {"overlap", {1, 1}, {s_axial}}}));
    expect_clean("wavefunction: within the relative tolerance",
                 wf_archive(s2, meta1, {c1, {"overlap", {1, 1}, {s_axial * (1 + 5e-7)}}}));
    expect_only("wavefunction: beyond the relative tolerance",
                wf_archive(s2, meta1, {c1, {"overlap", {1, 1}, {s_axial * (1 + 2e-6)}}}),
                "E_NORMALIZATION");

    // s + spherical p: n_ao 4.
    const std::string sp =
        R"({"structure_ref":"structure","shells":[)"
        R"({"center":0,"l":0,"exponents":[1.0],"coefficients":[1.0]},)"
        R"({"center":0,"l":1,"exponents":[0.5],"coefficients":[1.0]}]})";
    const std::string sp_ne =
        R"({"structure_ref":"structure","n_electrons":2,"shells":[)"
        R"({"center":0,"l":0,"exponents":[1.0],"coefficients":[1.0]},)"
        R"({"center":0,"l":1,"exponents":[0.5],"coefficients":[1.0]}]})";
    auto meta4 = [](const std::string& extra) {
        return R"({"spin":"restricted","energies":[0,0,0,0],"occupations":[2,0,0,0])" +
               extra + "}";
    };
    const std::string uhf =
        R"({"spin":"unrestricted","alpha":{"energies":[0,0,0,0],"occupations":[1,0,0,0]},)"
        R"("beta":{"energies":[0,0,0,0],"occupations":[1,0,0,0]}})";
    const WfMember c4{"mo_coefficients", {4, 4}, {}};
    const WfMember ca{"mo_coefficients_alpha", {4, 4}, {}};
    const WfMember cb{"mo_coefficients_beta", {4, 4}, {}};
    expect_only("wavefunction: complex true with rank-2 coefficients",
                wf_archive(sp, meta4(R"(,"complex":true)"), {c4}), "E_MEMBER_MALFORMED");
    expect_only("wavefunction: real with rank-3 coefficients",
                wf_archive(sp, meta4(""), {{"mo_coefficients", {4, 4, 2}, {}}}),
                "E_MEMBER_MALFORMED");
    expect_clean("wavefunction: complex rank-3 coefficients",
                 wf_archive(sp, meta4(R"(,"complex":true)"),
                            {{"mo_coefficients", {4, 4, 2}, {}}}));
    expect_clean("wavefunction: coefficient_encoding read when complex is absent",
                 wf_archive(sp, meta4(R"(,"coefficient_encoding":"complex_split_last_axis")"),
                            {{"mo_coefficients", {4, 4, 2}, {}}}));
    expect_only("wavefunction: symmetry_labels one short",
                wf_archive(sp, meta4(R"(,"symmetry_labels":["a","b","c"])"), {c4}),
                "E_MEMBER_MALFORMED");
    expect_only("wavefunction: overlap of the wrong shape",
                wf_archive(sp, meta4(""), {c4, {"overlap", {3, 3}, {}}}), "E_MEMBER_MALFORMED");
    expect_only("wavefunction: density_matrix with alpha and beta coefficients",
                wf_archive(sp, uhf, {ca, cb, {"density_matrix", {4, 4}, {}}}),
                "E_MEMBER_MALFORMED");
    expect_only("wavefunction: density_matrix_alpha without beta",
                wf_archive(sp, uhf, {ca, cb, {"density_matrix_alpha", {4, 4}, {}}}),
                "E_MEMBER_MALFORMED");
    expect_only("wavefunction: density_matrix_alpha with a single set",
                wf_archive(sp, meta4(""), {c4, {"density_matrix_alpha", {4, 4}, {}},
                                           {"density_matrix_beta", {4, 4}, {}}}),
                "E_MEMBER_MALFORMED");
    expect_clean("wavefunction: the unrestricted density pair",
                 wf_archive(sp_ne, uhf, {ca, cb, {"density_matrix_alpha", {4, 4}, {}},
                                         {"density_matrix_beta", {4, 4}, {}}}));
    expect_only("wavefunction: n_electrons disagrees with the occupations",
                wf_archive(sp_ne, R"({"spin":"restricted","energies":[0,0,0,0],)"
                                  R"("occupations":[1,0,0,0]})",
                           {c4}),
                "E_MEMBER_MALFORMED");
    Json prov = Json::object();
    prov.set("provenance", Json::parse(R"({"n_electrons":2})"));
    expect_clean("wavefunction: provenance.n_electrons agrees",
                 wf_archive(sp_ne, meta4(""), {c4}, prov));
    prov.set("provenance", Json::parse(R"({"n_electrons":3})"));
    expect_only("wavefunction: provenance.n_electrons disagrees",
                wf_archive(sp_ne, meta4(""), {c4}, prov), "E_ROOT_INCONSISTENT");
    expect_clean("wavefunction: NTO weights are not electrons",
                 wf_archive(sp_ne, R"({"spin":"restricted","occupation_semantics":)"
                                   R"("transition_weight","energies":[0,0,0,0],)"
                                   R"("occupations":[0.9,0.1,0,0]})",
                            {c4}));
}

// --- negative corpus ----------------------------------------------------------

// --- trajectory / reaction.path members (spec § 5.5) ---------------------------

struct FrameMember {
    std::string role;
    std::vector<int64_t> shape;
    std::vector<double> values;  // empty: zeros
};

static std::string f64_bytes(const std::vector<double>& v) {
    std::string out(8 * v.size(), '\0');
    for (size_t k = 0; k < v.size(); ++k) {
        uint64_t bits = 0;
        std::memcpy(&bits, &v[k], 8);
        for (int b = 0; b < 8; ++b)
            out[8 * k + static_cast<size_t>(b)] = static_cast<char>((bits >> (8 * b)) & 0xFF);
    }
    return out;
}

// One `kind` section with metadata `meta`, coords [4, 3, 3] of `coords_dtype`,
// and the float64 `extra` members.
static std::vector<uint8_t> trajectory_archive(const std::string& kind, const std::string& meta,
                                               const std::string& coords_dtype,
                                               const std::vector<FrameMember>& extra) {
    std::vector<Entry> files = {{"t/metadata.json", meta}};
    Json members = Json::object();
    members.set("metadata", json_member("t/metadata.json", meta));
    const std::string coords(36 * (coords_dtype == "float32" ? 4 : 8), '\0');
    files.push_back({"t/coords.bin", coords});
    Json c = Json::object();
    c.set("path", "t/coords.bin").set("format", "binary").set("dtype", coords_dtype)
        .set("shape", Json::from_ints({4, 3, 3})).set("sha256", sha(coords));
    members.set("coords", c);
    for (const FrameMember& m : extra) {
        size_t n = 1;
        Json sh = Json::array();
        for (int64_t d : m.shape) {
            n *= static_cast<size_t>(d);
            sh.push_back(Json(static_cast<long long>(d)));
        }
        const std::string data =
            f64_bytes(m.values.empty() ? std::vector<double>(n, 0.0) : m.values);
        const std::string path = "t/" + m.role + ".bin";
        files.push_back({path, data});
        Json spec = Json::object();
        spec.set("path", path).set("format", "binary").set("dtype", "float64")
            .set("shape", sh).set("sha256", sha(data));
        members.set(m.role, spec);
    }
    return archive(manifest_of({section("t", kind, members)}), files);
}

// Expect E_MEMBER_MALFORMED (and nothing else) located at members/`role`.
static void expect_member_malformed(const std::string& label, std::vector<uint8_t> bytes,
                                    const std::string& role) {
    std::vector<Finding> fs = semantic_findings(std::move(bytes));
    std::set<std::string> got = codes(fs);
    got.erase("W_SCHEMA_SKIPPED");
    bool at_role = false;
    for (const auto& f : fs)
        if (f.code == "E_MEMBER_MALFORMED" && f.where == "/sections/0/members/" + role)
            at_role = true;
    if (got != std::set<std::string>{"E_MEMBER_MALFORMED"} || !at_role) {
        std::cerr << "FAIL: " << label << ": expected only E_MEMBER_MALFORMED at " << role
                  << ", got" << describe(fs) << "\n";
        ++g_failures;
    }
}

static void test_trajectory_members() {
    const std::string sym = R"({"symbols":["O","H","H"])";
    const std::string pbc = sym + R"(,"pbc":[true,true,false])";
    const std::vector<double> e4 = {-76.0, -76.1, -76.2, -76.3};
    expect_clean("trajectory: every member with matching extents",
                 trajectory_archive(
                     "trajectory", pbc + R"(,"energies":[-76.0,-76.1,-76.2,-76.3]})",
                     "float64",
                     {{"lattice", {4, 3, 3}, {}}, {"velocities", {4, 3, 3}, {}},
                      {"gradient", {4, 3, 3}, {}}, {"times", {4}, {}},
                      {"energies", {4}, e4}, {"cell_volumes", {4}, {}},
                      {"temperatures", {4}, {}}}));
    expect_clean("trajectory: float32 coords with precision float32",
                 trajectory_archive("trajectory", sym + R"(,"precision":"float32"})",
                                    "float32", {}));
    expect_clean("trajectory: a fixed [3, 3] cell and all-false pbc",
                 trajectory_archive("trajectory", sym + R"(,"pbc":[false,false,false]})",
                                    "float64", {{"lattice", {3, 3}, {}}}));
    expect_clean("reaction.path: a lattice needs no pbc; shared members",
                 trajectory_archive("reaction.path", sym + "}", "float64",
                                    {{"lattice", {4, 3, 3}, {}}, {"velocities", {4, 3, 3}, {}},
                                     {"gradient", {4, 3, 3}, {}}, {"times", {4}, {}},
                                     {"energies", {4}, {}}}));
    expect_clean("trajectory: energies within 1e-10 of metadata.energies",
                 trajectory_archive("trajectory",
                                    sym + R"(,"energies":[-76.0,-76.1,-76.2,-76.3]})",
                                    "float64",
                                    {{"energies", {4}, {-76.0, -76.1 + 5e-11, -76.2, -76.3}}}));

    expect_member_malformed("trajectory: float32 coords without precision",
                            trajectory_archive("trajectory", sym + "}", "float32", {}),
                            "coords");
    expect_member_malformed(
        "trajectory: precision float32 with float64 coords",
        trajectory_archive("trajectory", sym + R"(,"precision":"float32"})", "float64", {}),
        "coords");
    expect_member_malformed("reaction.path: float32 coords without precision",
                            trajectory_archive("reaction.path", sym + "}", "float32", {}),
                            "coords");
    expect_member_malformed(
        "trajectory: velocities with the wrong frame count",
        trajectory_archive("trajectory", sym + "}", "float64", {{"velocities", {5, 3, 3}, {}}}),
        "velocities");
    expect_member_malformed(
        "trajectory: gradient with the wrong atom count",
        trajectory_archive("trajectory", sym + "}", "float64", {{"gradient", {4, 2, 3}, {}}}),
        "gradient");
    expect_member_malformed(
        "trajectory: times with the wrong frame count",
        trajectory_archive("trajectory", sym + "}", "float64", {{"times", {3}, {}}}), "times");
    expect_member_malformed(
        "trajectory: cell_volumes of rank 2",
        trajectory_archive("trajectory", sym + "}", "float64", {{"cell_volumes", {4, 1}, {}}}),
        "cell_volumes");
    expect_member_malformed(
        "trajectory: temperatures with the wrong frame count",
        trajectory_archive("trajectory", sym + "}", "float64", {{"temperatures", {1}, {}}}),
        "temperatures");
    expect_member_malformed(
        "reaction.path: gradient with the wrong frame count",
        trajectory_archive("reaction.path", sym + "}", "float64", {{"gradient", {3, 3, 3}, {}}}),
        "gradient");
    expect_member_malformed(
        "trajectory: per-frame lattice with the wrong frame count",
        trajectory_archive("trajectory", pbc + "}", "float64", {{"lattice", {3, 3, 3}, {}}}),
        "lattice");
    expect_member_malformed(
        "reaction.path: per-frame lattice with the wrong frame count",
        trajectory_archive("reaction.path", sym + "}", "float64", {{"lattice", {5, 3, 3}, {}}}),
        "lattice");
    expect_member_malformed(
        "trajectory: binary energies disagree with metadata.energies",
        trajectory_archive("trajectory", sym + R"(,"energies":[-76.0,-76.1,-76.2,-76.3]})",
                           "float64", {{"energies", {4}, {-76.0, -76.1, -76.2 + 2e-10, -76.3}}}),
        "energies");
    expect_member_malformed(
        "reaction.path: binary energies longer than metadata.energies",
        trajectory_archive("reaction.path", sym + R"(,"energies":[-76.0,-76.1,-76.2]})",
                           "float64", {{"energies", {4}, e4}}),
        "energies");
    expect_member_malformed(
        "trajectory: a lattice without pbc",
        trajectory_archive("trajectory", sym + "}", "float64", {{"lattice", {3, 3}, {}}}),
        "metadata");
    expect_member_malformed(
        "trajectory: a periodic pbc without a lattice",
        trajectory_archive("trajectory", sym + R"(,"pbc":[false,true,false]})", "float64", {}),
        "metadata");
}

// --- the generated JSON Schema checker ---------------------------------------------

// The schema errors of a JSON text against a $defs entry ("" = the root schema).
static std::vector<detail::schema::Error> schema_errors(const std::string& def,
                                                        const std::string& text) {
    const Json doc = Json::parse(text);
    return def.empty() ? detail::schema::check_manifest(doc)
                       : detail::schema::check_definition(def, doc);
}

// The pointers of those errors, in order.
static std::vector<std::string> schema_pointers(const std::string& def, const std::string& text) {
    std::vector<std::string> out;
    for (const auto& e : schema_errors(def, text)) out.push_back(e.pointer);
    return out;
}

static std::string joined(const std::vector<std::string>& v) {
    std::string s;
    for (const auto& x : v) s += " '" + x + "'";
    return s;
}

#define CHECK_POINTERS(def, text, ...)                                                 \
    do {                                                                               \
        const std::vector<std::string> got_ = schema_pointers(def, text);              \
        const std::vector<std::string> want_ = __VA_ARGS__;                            \
        CHECK(got_ == want_, std::string(def) + " " + (text) + ": got" + joined(got_)); \
    } while (0)

using Ptrs = std::vector<std::string>;

static void test_schema_checker() {
    // Patterns: the hand-free Pike VM over code points.
    const std::string hex64(64, 'a');
    CHECK_POINTERS("Sha256Hex", "\"" + hex64 + "\"", Ptrs{});
    CHECK_POINTERS("Sha256Hex", "\"" + std::string(64, 'A') + "\"", Ptrs{""});
    CHECK_POINTERS("Sha256Hex", "\"" + hex64 + "a\"", Ptrs{""});
    CHECK_POINTERS("Sha256Hex", "\"" + hex64 + "\\n\"", Ptrs{""});   // ECMA-262 `$`
    CHECK_POINTERS("ZipPath", R"("dir/file.json")", Ptrs{});
    CHECK_POINTERS("ZipPath", R"("a/..b/c")", Ptrs{});
    CHECK_POINTERS("ZipPath", R"("a/../b")", Ptrs{""});
    CHECK_POINTERS("ZipPath", R"("..")", Ptrs{""});
    CHECK_POINTERS("ZipPath", R"("/abs")", Ptrs{""});
    CHECK_POINTERS("ZipPath", R"("")", (Ptrs{"", ""}));   // minLength and pattern
    CHECK_POINTERS("ZipPath", R"("caf\u00e9/\u00fc.json")", Ptrs{});
    CHECK_POINTERS("Metadata", R"({"uuid":"3f2504e0-4f89-41d3-9a0c-0305e82c3301"})", Ptrs{});
    CHECK_POINTERS("Metadata", R"({"uuid":"3f2504e04f8941d39a0c0305e82c3301"})", Ptrs{"/uuid"});
    CHECK_POINTERS("Metadata", R"({"created_utc":"2026-06-30T23:59:60.25Z"})", Ptrs{});
    CHECK_POINTERS("Metadata", R"({"created_utc":"2026-13-01T00:00:00Z"})",
                   Ptrs{"/created_utc"});
    CHECK_POINTERS("Metadata", R"({"x_tool":{"anything":[1]},"license":"CC-BY-4.0"})", Ptrs{});
    CHECK_POINTERS("Metadata", R"({"licence":"x"})", Ptrs{""});   // additionalProperties
    // Types and numbers: integers are numbers with a zero fraction, booleans
    // are not numbers, 1 == 1.0 for const.
    CHECK_POINTERS("MemberSize", "2.0", Ptrs{});
    CHECK_POINTERS("MemberSize", "2.5", Ptrs{""});
    CHECK_POINTERS("MemberSize", "true", Ptrs{""});
    CHECK_POINTERS("MemberSize", "-1", Ptrs{""});
    CHECK_POINTERS("Camera",
                   R"({"position":[0,0,1],"focal_point":[0,0,0],"view_up":[0,1,0],"view_angle":180})",
                   Ptrs{"/view_angle"});   // exclusiveMaximum
    // Code points, not bytes: a three-letter symbol of two-byte characters.
    CHECK_POINTERS("StructureAtom",
                   R"({"symbol":"\u00c5\u00c5\u00c5","position":[0,0,0],"atomic_number":0})",
                   Ptrs{});
    CHECK_POINTERS("StructureAtom",
                   R"({"symbol":"\u00c5\u00c5\u00c5\u00c5","position":[0,0],"atomic_number":-1})",
                   (Ptrs{"/symbol", "/position", "/atomic_number"}));
    // An escaped surrogate pair is one code point (as in Python and JS).
    CHECK(Json::parse(R"("😀")").as_string() == "\xF0\x9F\x98\x80",
          "an escaped surrogate pair decodes to one UTF-8 code point");
    CHECK_POINTERS("StructureAtom",
                   R"({"symbol":"😀😀😀","position":[0,0,0],"atomic_number":0})",
                   Ptrs{});
    // A repeated key keeps its last value, as in Python.
    CHECK_POINTERS("ValueWithUnits", R"({"value":"x","value":1,"units":"ev"})", Ptrs{});
    CHECK_POINTERS("ValueWithUnits", R"({"value":1,"units":"Eh"})", Ptrs{});   // oneOf
    CHECK_POINTERS("ValueWithUnits", R"({"value":1,"units":"EV"})", Ptrs{"/units"});
    // Pointers escape '~' and '/'; patternProperties descend.
    {
        const auto e = schema_errors("EnergyComponents",
                                     R"({"x_a/b~c":{"value":"no","units":"ev"}})");
        CHECK(e.size() == 1 && e[0].pointer == "/x_a~1b~0c/value" &&
                  e[0].location == "x_a/b~c/value",
              "pointer escaping in patternProperties");
    }
    // prefixItems and items.
    CHECK_POINTERS("AtomVectorBinary",
                   R"({"path":"p","format":"binary","dtype":"float64","shape":[5,4],"sha256":")" +
                       hex64 + "\"}",
                   Ptrs{"/shape/1"});
    CHECK_POINTERS("Vec3", R"([0, "y", 2])", Ptrs{"/1"});
    // contains / minContains / maxContains, if / then.
    CHECK_POINTERS("StructurePayload",
                   R"({"atoms":[],"pbc":[true,true,false],"dimensionality":2,)"
                   R"("lattice_vectors":[[1,0,0],[0,1,0],[0,0,1]]})",
                   Ptrs{});
    CHECK_POINTERS("StructurePayload",
                   R"({"atoms":[],"pbc":[true,false,false],"dimensionality":2,)"
                   R"("lattice_vectors":[[1,0,0],[0,1,0],[0,0,1]]})",
                   Ptrs{"/pbc"});
    CHECK_POINTERS("StructurePayload", R"({"atoms":[],"lattice_vectors":[[1,0,0],[0,1,0],[0,0,1]]})",
                   Ptrs{""});
    // if / then / else with not + anyOf.
    CHECK_POINTERS("MoMetadataPayload", R"({"spin":"restricted","energies":[],"occupations":[]})",
                   Ptrs{});
    CHECK_POINTERS("MoMetadataPayload",
                   R"({"spin":"unrestricted","energies":[1],"occupations":[1]})",
                   (Ptrs{"", "", ""}));
    // dependentRequired.
    CHECK_POINTERS("SectionVolumeDifference", R"({"operand_a":"x"})", Ptrs{""});
    // oneOf over the section kinds: one error at the section, whatever is
    // wrong inside it; a non-object section matches many branches.
    CHECK_POINTERS("Section", R"({"id":"s","kind":"bonds","members":{}})", Ptrs{""});
    CHECK_POINTERS("Section", "42", (Ptrs{"", ""}));
    CHECK_POINTERS("Section", R"({"id":"a b","kind":"x_v.k","members":{}})", Ptrs{"/id"});

    // Through the validator: findings, codes and locations.
    const std::vector<Finding> clean = QvfReader::validate_bytes(base_archive());
    CHECK(clean.empty(), "the base archive is schema-clean" + describe(clean));
    Json big = Json::object();
    big.set("byte_order", "big");
    std::vector<Finding> fs = QvfReader::validate_bytes(
        archive(manifest_of({structure_section()}, big), {{"structure/structure.json", kWater}}));
    CHECK(fs.size() == 1 && fs[0].code == "E_SCHEMA" && fs[0].where == "/byte_order",
          "byte_order big is E_SCHEMA at /byte_order" + describe(fs));
    fs = semantic_findings(
        archive(manifest_of({structure_section()}, big), {{"structure/structure.json", kWater}}));
    CHECK(codes(fs) == std::set<std::string>{"W_SCHEMA_SKIPPED"} && !has_errors(fs),
          "without the schema check byte_order big is clean" + describe(fs));
    const std::string bad_atom =
        R"({"atoms":[{"symbol":"O","position":[0,0,0],"atomic_number":-8}]})";
    fs = QvfReader::validate_bytes(
        archive(manifest_of({structure_section("structure", "structure/structure.json", bad_atom)}),
                {{"structure/structure.json", bad_atom}}));
    CHECK(fs.size() == 1 && fs[0].code == "E_PAYLOAD_SCHEMA" &&
              fs[0].where == "/sections/0/members/structure" &&
              fs[0].message.find("atoms/0/atomic_number") != std::string::npos,
          "a payload error is E_PAYLOAD_SCHEMA at the member" + describe(fs));
    Json version = Json::object();
    version.set("qvf_version", 1.0);
    fs = QvfReader::validate_bytes(archive(manifest_of({structure_section()}, version),
                                           {{"structure/structure.json", kWater}}));
    CHECK(fs.empty(), "qvf_version 1.0 equals the const 1" + describe(fs));
    // The spectrum equal-length rule runs with and without the schema check.
    const std::string spec = R"({"frequencies":[1,2,3],"intensities":[1,2]})";
    Json sm = Json::object();
    sm.set("spectrum", json_member("s/spectrum.json", spec));
    const std::vector<uint8_t> sa =
        archive(manifest_of({structure_section(), section("s", "spectra.ir", sm)}),
                {{"structure/structure.json", kWater}, {"s/spectrum.json", spec}});
    fs = QvfReader::validate_bytes(sa);
    CHECK(codes(fs) == std::set<std::string>{"E_PAYLOAD_SCHEMA"} &&
              fs[0].where == "/sections/1/members/spectrum",
          "spectrum lengths" + describe(fs));
    fs = semantic_findings(sa);
    CHECK(codes(fs) == (std::set<std::string>{"E_PAYLOAD_SCHEMA", "W_SCHEMA_SKIPPED"}),
          "spectrum lengths without the schema check" + describe(fs));

    // The generated tables.
    const gen::KindInfo* ap = gen::kind_info("atom_properties");
    CHECK(ap && !ap->allow_extra_members && gen::role_info(*ap, "magnetic_moment") &&
              gen::role_info(*ap, "magnetic_moment")->min_rank == 2 &&
              gen::role_info(*ap, "magnetic_moment")->format == gen::Format::Binary,
          "KindInfo / RoleInfo of atom_properties");
    CHECK(gen::kind_info("x_vendor.thing") == nullptr && gen::kind_info("bogus") == nullptr,
          "kind_info of non-canonical kinds");
    CHECK(gen::payload_def("structure", "structure") == "StructurePayload" &&
              gen::payload_def("structure", "nope").empty(),
          "payload_def");
    CHECK(gen::finding_code("E_SCHEMA") && gen::finding_code("E_SCHEMA")->error &&
              !gen::finding_code("W_ORPHAN_MEMBER")->error && !gen::finding_code("E_NOPE"),
          "finding codes");
    CHECK(gen::kSchemaId == "https://vibe-qc.org/spec/qvf/1/manifest.schema.json",
          "schema $id");

    // The typed model of an opened archive.
    QvfReader r = QvfReader::from_bytes(base_archive());
    const Manifest& m = r.typed_manifest();
    CHECK(m.qvf_version == 1 && m.source && m.source->program == "t" && m.sections &&
              m.sections->size() == 1 && m.section("structure") &&
              m.section("structure")->member("structure") &&
              m.section("structure")->member("structure")->format == MemberSpec::Format::Json,
          "typed_manifest()");
    CHECK(detail::schema::json_equal(m.to_json(), r.manifest()), "typed manifest round trip");
}

// --- the typed manifest model (qvf/model.hpp) ---------------------------------------

// Exact JSON equality for the model's round trip: an integer is not a double,
// objects compare as Python sees them (each key once, last value wins, order
// ignored).
static bool same_json(const Json& a, const Json& b) {
    if (a.type() != b.type()) return false;
    switch (a.type()) {
        case Json::Type::Null: return true;
        case Json::Type::Bool: return a.as_bool() == b.as_bool();
        case Json::Type::Int: return a.as_int() == b.as_int();
        case Json::Type::Double: return a.as_double() == b.as_double();
        case Json::Type::String: return a.as_string() == b.as_string();
        case Json::Type::Array:
            if (a.size() != b.size()) return false;
            for (size_t i = 0; i < a.size(); ++i)
                if (!same_json(a.at(i), b.at(i))) return false;
            return true;
        case Json::Type::Object: {
            const auto ua = detail::schema::unique_members(a);
            if (ua.size() != detail::schema::unique_members(b).size()) return false;
            for (const auto& kv : ua) {
                const Json* other = detail::schema::find_last(b, *kv.first);
                if (!other || !same_json(*kv.second, *other)) return false;
            }
            return true;
        }
        default: return false;
    }
}

static void test_model() {
    // Every field with the wrong type is kept verbatim; nothing is lost.
    const Json doc = Json::parse(R"({
        "qvf_version": 1.0, "schema_uri": 7, "byte_order": "big",
        "source": {"program": "p", "version": "1"},
        "provenance": [1], "metadata": {"title": "x"}, "x_vendor": {"a": 1},
        "sections": [
          {"id": 3, "kind": "structure", "label": null, "critical": "yes",
           "structure_ref": "s", "operand_a": "a",
           "members": {"structure": {"path": "s.json", "format": "json", "sha256": "0",
                                     "size": -1, "shape": [1.5], "x": true,
                                     "media_type": 5}}},
          {"id": "b", "kind": "bonds", "members": {"bonds": "not an object"}},
          {"id": "c", "kind": "x_v.k", "members": {"m": {"format": "text", "dtype": 4}}},
          {"id": "d", "kind": "x_v.k", "critical": false,
           "members": {"o": {"path": "o.bin", "format": "binary", "sha256": "1", "size": 3,
                             "media_type": "image/png"},
                       "t": {"path": "t.bin", "format": "binary", "dtype": "int8",
                             "shape": [3], "sha256": "2"}}}
        ]})");
    const Manifest m = Manifest::from_json(doc);
    CHECK(!m.qvf_version && m.extra.has("qvf_version") && m.byte_order == std::string("big") &&
              !m.source && m.extra.has("source") && !m.provenance && m.metadata,
          "model: root fields typed or kept");
    CHECK(m.sections && m.sections->size() == 4, "model: sections typed");
    const Section& s0 = (*m.sections)[0];
    CHECK(!s0.id && s0.extra.has("id") && !s0.critical && s0.extra.has("operand_a") &&
              s0.structure_ref == std::string("s") && s0.member("structure") &&
              !s0.member("structure")->size && !s0.member("structure")->shape &&
              !s0.member("structure")->media_type &&
              s0.member("structure")->extra.size() == 4,
          "model: section and member fields typed or kept");
    CHECK(!(*m.sections)[1].members && (*m.sections)[1].extra.has("members"),
          "model: a members object with a non-object spec is kept whole");
    const Section& s3 = (*m.sections)[3];
    CHECK(s3.critical == false && s3.member("o")->format == MemberSpec::Format::Opaque &&
              s3.member("o")->size == 3u &&
              s3.member("o")->media_type == std::string("image/png") &&
              !s3.member("o")->extra.has("media_type") &&
              s3.member("t")->format == MemberSpec::Format::Binary &&
              s3.member("t")->shape == std::vector<int64_t>{3},
          "model: opaque and typed binary members");
    CHECK(same_json(m.to_json(), doc), "model: lossless round trip of a malformed manifest");
    bool threw = false;
    try {
        Manifest::from_json(Json::array());
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw, "model: a non-object manifest throws std::invalid_argument");
}

// Every archive in the directories that opens: typed_manifest() round-trips
// to manifest() exactly, and to_json() is a fixed point.
static int run_model_roundtrip(const std::vector<std::string>& dirs) {
    size_t tested = 0, skipped = 0, failed = 0;
    bool any_dir = false;
    for (const std::string& dir : dirs) {
        if (!std::filesystem::is_directory(dir)) continue;
        any_dir = true;
        std::vector<std::filesystem::path> archives;
        for (const auto& e : std::filesystem::directory_iterator(dir))
            if (e.path().extension() == ".qvf") archives.push_back(e.path());
        std::sort(archives.begin(), archives.end());
        for (const std::filesystem::path& path : archives) {
            Json manifest;
            Json back;
            try {
                QvfReader r = QvfReader::open(path.string());
                manifest = r.manifest();
                back = r.typed_manifest().to_json();
            } catch (const std::exception&) {
                ++skipped;   // does not open (the negative corpus), or not an object
                continue;
            }
            ++tested;
            const Json again = Manifest::from_json(back).to_json();
            if (!same_json(back, manifest) || again.dump() != back.dump()) {
                ++failed;
                std::cout << "FAIL " << path.filename().string() << "\n";
            }
        }
    }
    if (!any_dir) {
        std::cout << "no corpus directory found; skipped\n";
        return 77;
    }
    std::cout << "model round trip (C++): " << tested - failed << "/" << tested
              << " manifests lossless (" << skipped << " archives do not open)\n";
    return failed ? 1 : 0;
}

static bool set_limit(Limits& lim, const std::string& key, const Json& v) {
    if (!v.is_number()) return false;
    if (key == "max_member_bytes") lim.max_member_bytes = static_cast<uint64_t>(v.as_int());
    else if (key == "max_total_bytes") lim.max_total_bytes = static_cast<uint64_t>(v.as_int());
    else if (key == "max_ratio") lim.max_ratio = v.as_double();
    else if (key == "ratio_exempt_bytes") lim.ratio_exempt_bytes = static_cast<uint64_t>(v.as_int());
    else if (key == "max_json_bytes") lim.max_json_bytes = static_cast<uint64_t>(v.as_int());
    else if (key == "max_json_depth") lim.max_json_depth = static_cast<size_t>(v.as_int());
    else return false;
    return true;
}

static int run_negative_corpus(const std::string& dir) {
    namespace fs = std::filesystem;
    if (!fs::is_directory(dir)) {
        std::cout << "negative corpus not found at " << dir << "; skipped\n";
        return 77;
    }
    std::vector<fs::path> expected_files;
    for (const auto& e : fs::directory_iterator(dir)) {
        const std::string name = e.path().filename().string();
        const std::string suffix = ".expected.json";
        if (name.size() > suffix.size() &&
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
            expected_files.push_back(e.path());
    }
    std::sort(expected_files.begin(), expected_files.end());
    if (expected_files.empty()) {
        std::cerr << "no *.expected.json in " << dir << "\n";
        return 1;
    }

    int failed = 0;
    size_t schema_only = 0;
    for (const bool with_schema : {true, false}) {
        ValidationOptions opts;
        opts.schema_check = with_schema;
        const char* mode = with_schema ? "schema check" : "no schema check";
        for (const auto& ef : expected_files) {
            std::ifstream fh(ef, std::ios::binary);
            std::stringstream text;
            text << fh.rdbuf();
            std::vector<std::string> problems;
            std::string archive_name = ef.filename().string();
            bool skipped_all = false;
            try {
                Json exp = Json::parse(text.str());
                archive_name = exp.at("archive").as_string();
                Limits lim;
                if (exp.contains("limits"))
                    for (const auto& kv : exp.at("limits").items())
                        if (!set_limit(lim, kv.first, kv.second))
                            problems.push_back("unknown limit '" + kv.first + "'");
                std::vector<Finding> findings = QvfReader::validate_file(
                    (fs::path(dir) / archive_name).string(), lim, opts);
                // Exact code sets apart from the schema codes, which the
                // negative archives do not all name (a schema engine also
                // reports E_SCHEMA for several of them); every expected code
                // must be present. Without the schema check (W_SCHEMA_SKIPPED)
                // the schema-only expectations are dropped and the verdict is
                // judged on the other codes, as conformance/run_conformance.py
                // does.
                const auto schema_code = [](const std::string& c) {
                    return c == "E_SCHEMA" || c == "E_PAYLOAD_SCHEMA" || c == "W_SCHEMA_SKIPPED";
                };
                std::set<std::string> got, want;
                bool other_errors = false;
                for (const auto& f : findings) {
                    if (!schema_code(f.code)) got.insert(f.code);
                    if (f.severity == Finding::Severity::Error &&
                        (with_schema || !schema_code(f.code)))
                        other_errors = true;
                }
                bool expects_error = false;
                for (const Json& w : exp.at("expect").elements()) {
                    const std::string code = w.at("code").as_string();
                    if (!with_schema && schema_code(code)) continue;
                    const bool is_error = w.at("severity").as_string() == "error";
                    expects_error = expects_error || is_error;
                    if (!schema_code(code)) want.insert(code);
                    if (!has_code(findings, code,
                                  is_error ? Finding::Severity::Error
                                           : Finding::Severity::Warning))
                        problems.push_back("missing " + code);
                }
                skipped_all = !with_schema && want.empty() && !expects_error &&
                              !exp.at("expect").elements().empty();
                if (with_schema && has_code(findings, "W_SCHEMA_SKIPPED",
                                            Finding::Severity::Warning))
                    problems.push_back("W_SCHEMA_SKIPPED with the schema check enabled");
                if (other_errors != expects_error)
                    problems.push_back(expects_error ? "validator reported no error"
                                                     : "validator reported an unexpected error");
                if (got != want)
                    problems.push_back("codes differ from the expected set:" +
                                       describe(findings));
            } catch (const std::exception& e) {
                problems.push_back(std::string("runner error: ") + e.what());
            }
            if (problems.empty() && skipped_all) {
                ++schema_only;
                std::cout << "PASS " << archive_name << "  [" << mode
                          << "; expects only JSON-Schema codes: accepted without errors]\n";
            } else if (problems.empty()) {
                std::cout << "PASS " << archive_name << "  [" << mode << "]\n";
            } else {
                ++failed;
                std::cout << "FAIL " << archive_name << "  [" << mode << "]\n";
                for (const auto& p : problems) std::cout << "     - " << p << "\n";
            }
        }
    }
    std::cout << "\nnegative conformance (C++): " << 2 * expected_files.size() - failed << "/"
              << 2 * expected_files.size()
              << " runs as expected (each archive with and without the schema check; "
              << schema_only << " schema-only archives accepted without the check)\n";
    return failed ? 1 : 0;
}

// --- excited states, spectra that name them, transition densities -----------
// (spec § 5.2, § 5.4, § 5.12)

static const std::string kTwoS =
    R"({"structure_ref":"structure","pure":true,"n_ao":2,"shells":[)"
    R"({"center":0,"l":0,"exponents":[1.0],"coefficients":[1.0]},)"
    R"({"center":1,"l":0,"exponents":[1.0],"coefficients":[1.0]}]})";

// The base structure, a wavefunction.gto "wf" (restricted with 2 orbital
// rows, or unrestricted with 2 alpha and 3 beta rows), an excited_states
// section "es" with `payload` and `es_fields`, and the caller's `extra`
// sections and files.
static std::vector<uint8_t> excited_archive(const std::string& payload,
                                            Json es_fields = Json::object(),
                                            bool unrestricted = false,
                                            const std::vector<Json>& extra = {},
                                            const std::vector<Entry>& extra_files = {}) {
    std::vector<Entry> files = {{"structure/structure.json", kWater}};
    const std::string mo =
        unrestricted
            ? R"({"n_ao":2,"spin":"unrestricted",)"
              R"("alpha":{"energies":[0,0],"occupations":[1,0]},)"
              R"("beta":{"energies":[0,0,0],"occupations":[1,0,0]}})"
            : R"({"n_mo":2,"n_ao":2,"spin":"restricted","energies":[-0.5,0.5],)"
              R"("occupations":[2,0]})";
    files.push_back({"wf/basis.json", kTwoS});
    files.push_back({"wf/mo.json", mo});
    Json wm = Json::object();
    wm.set("basis", json_member("wf/basis.json", kTwoS))
        .set("mo_metadata", json_member("wf/mo.json", mo));
    if (unrestricted) {
        wm.set("mo_coefficients_alpha", f64_member("wf/ca.bin", {2, 2}, files));
        wm.set("mo_coefficients_beta", f64_member("wf/cb.bin", {3, 2}, files));
    } else {
        wm.set("mo_coefficients", f64_member("wf/c.bin", {2, 2}, files));
    }
    Json es_members = Json::object();
    es_members.set("states", json_member("es/states.json", payload));
    files.push_back({"es/states.json", payload});
    Json es = section("es", "excited_states", es_members);
    for (const auto& kv : es_fields.items()) es.set(kv.first, kv.second);
    std::vector<Json> secs = {structure_section(), section("wf", "wavefunction.gto", wm), es};
    for (const auto& s : extra) secs.push_back(s);
    for (const auto& e : extra_files) files.push_back(e);
    return archive(manifest_of(secs), files);
}

static Json with_field(Json obj, const std::string& key, Json value) {
    obj.set(key, std::move(value));
    return obj;
}

// A volume.transition_density section "td" (a 2x2x2 float64 grid) with the
// given section fields; its files go to `files`.
static Json transition_density(Json fields, std::vector<Entry>& files) {
    const std::string grid =
        R"({"origin":[0,0,0],"voxel_vectors":[[1,0,0],[0,1,0],[0,0,1]],"shape":[2,2,2]})";
    files.push_back({"td/grid.json", grid});
    Json members = Json::object();
    members.set("grid", json_member("td/grid.json", grid))
        .set("data", f64_member("td/data.bin", {2, 2, 2}, files));
    Json s = section("td", "volume.transition_density", members);
    for (const auto& kv : fields.items()) s.set(kv.first, kv.second);
    return s;
}

// A spectra.uvvis section "uv" with `payload`; its file goes to `files`.
static Json uvvis(const std::string& payload, std::vector<Entry>& files) {
    files.push_back({"uv/spectrum.json", payload});
    Json members = Json::object();
    members.set("spectrum", json_member("uv/spectrum.json", payload));
    return section("uv", "spectra.uvvis", members);
}

static void test_excited_states() {
    const std::string two =
        R"({"states":[{"index":1,"excitation_energy":{"value":7.4,"units":"ev"},)"
        R"("contributions":[{"from_mo":0,"to_mo":1,"weight":0.9}]},)"
        R"({"index":2,"excitation_energy":{"value":9.2,"units":"ev"}}],)"
        R"("meta":{"method":"tddft","n_states":2}})";
    Json wf_ref = Json::object();
    wf_ref.set("wavefunction_ref", "wf");
    expect_clean("excited_states: two states, contributions within the 2 rows",
                 excited_archive(two, wf_ref));

    auto one = [](const std::string& state) {
        return R"({"states":[)" + state + "]}";
    };
    const std::string e1 = R"("excitation_energy":{"value":7.4,"units":"ev"})";
    expect_only("excited_states: a repeated index",
                excited_archive(R"({"states":[{"index":1,)" + e1 + R"(},{"index":1,)" + e1 +
                                "}]}"),
                "E_MEMBER_MALFORMED");
    expect_only("excited_states: indices out of order",
                excited_archive(R"({"states":[{"index":2,)" + e1 + R"(},{"index":1,)" + e1 +
                                "}]}"),
                "E_MEMBER_MALFORMED");
    expect_only("excited_states: meta.n_states differs from the states",
                excited_archive(R"({"states":[{"index":1,)" + e1 +
                                R"(}],"meta":{"n_states":2}})"),
                "E_MEMBER_MALFORMED");
    expect_only("excited_states: wavefunction_ref names a structure",
                excited_archive(two, with_field(Json::object(), "wavefunction_ref",
                                                "structure")),
                "E_REF_UNRESOLVED");
    expect_only("excited_states: wavefunction_ref names no section",
                excited_archive(two, with_field(Json::object(), "wavefunction_ref", "nope")),
                "E_REF_UNRESOLVED");
    expect_clean("excited_states: nto_ref names a wavefunction.gto",
                 excited_archive(one(R"({"index":1,)" + e1 + R"(,"nto_ref":"wf"})")));
    expect_only("excited_states: nto_ref names the excited_states section itself",
                excited_archive(one(R"({"index":1,)" + e1 + R"(,"nto_ref":"es"})")),
                "E_REF_UNRESOLVED");
    expect_only("excited_states: density_ref names a wavefunction.gto",
                excited_archive(one(R"({"index":1,)" + e1 + R"(,"density_ref":"wf"})")),
                "E_REF_UNRESOLVED");
    const std::string to2 =
        one(R"({"index":1,)" + e1 + R"(,"contributions":[{"from_mo":0,"to_mo":2,"weight":1}]})");
    expect_only("excited_states: to_mo beyond the 2 orbital rows",
                excited_archive(to2, wf_ref), "E_MEMBER_MALFORMED");
    expect_clean("excited_states: without wavefunction_ref the rows are not checked",
                 excited_archive(to2));
    auto spin_to2 = [&](const std::string& spin) {
        return one(R"({"index":1,)" + e1 + R"(,"contributions":[{"from_mo":0,"to_mo":2,)" +
                   spin + R"("weight":1}]})");
    };
    expect_clean("excited_states: beta to_mo 2 is one of the 3 beta rows",
                 excited_archive(spin_to2(R"("spin":"beta",)"), wf_ref, true));
    expect_only("excited_states: alpha to_mo 2 is beyond the 2 alpha rows",
                excited_archive(spin_to2(R"("spin":"alpha",)"), wf_ref, true),
                "E_MEMBER_MALFORMED");
    expect_only("excited_states: without spin, to_mo must be a row of both channels",
                excited_archive(spin_to2(""), wf_ref, true), "E_MEMBER_MALFORMED");
    expect_only("excited_states: an unknown meta.units name is a warning",
                excited_archive(R"({"states":[{"index":1,)" + e1 +
                                R"(}],"meta":{"units":{"excitation_energy":"electronvolt"}}})"),
                "W_UNIT_UNKNOWN", Limits{}, Finding::Severity::Warning);

    {   // volume.transition_density names the table and a state of it
        std::vector<Entry> files;
        Json fields = Json::object();
        fields.set("excited_states_ref", "es").set("state_index", 2);
        const std::string dens = R"({"states":[{"index":1,)" + e1 +
                                 R"(,"density_ref":"td"},{"index":2,)" + e1 + "}]}";
        Json td = transition_density(fields, files);
        expect_clean("transition density: names state 2 of es",
                     excited_archive(dens, Json::object(), false, {td}, files));
    }
    {
        std::vector<Entry> files;
        Json fields = Json::object();
        fields.set("excited_states_ref", "es").set("state_index", 5);
        Json td = transition_density(fields, files);
        expect_only("transition density: state_index names no state",
                    excited_archive(two, Json::object(), false, {td}, files),
                    "E_REF_UNRESOLVED");
    }
    {
        std::vector<Entry> files;
        Json td = transition_density(with_field(Json::object(), "excited_states_ref", "wf"),
                                     files);
        expect_only("transition density: excited_states_ref names a wavefunction",
                    excited_archive(two, Json::object(), false, {td}, files),
                    "E_REF_UNRESOLVED");
    }
    {   // spectra name states through states_ref
        const std::string peaks =
            R"({"frequencies":[7.4,9.2],"intensities":[0.05,0.0],)"
            R"("peaks":[{"x":7.4,"assignment":{"state_index":%d}}],"states_ref":"%s"})";
        auto spectrum = [&](int state, const std::string& ref) {
            char buf[256];
            std::snprintf(buf, sizeof(buf), peaks.c_str(), state, ref.c_str());
            return std::string(buf);
        };
        std::vector<Entry> f1;
        Json ok = uvvis(spectrum(2, "es"), f1);
        expect_clean("spectrum: states_ref and state_index resolve",
                     excited_archive(two, Json::object(), false, {ok}, f1));
        std::vector<Entry> f2;
        Json bad_state = uvvis(spectrum(7, "es"), f2);
        expect_only("spectrum: state_index names no state of es",
                    excited_archive(two, Json::object(), false, {bad_state}, f2),
                    "E_REF_UNRESOLVED");
        std::vector<Entry> f3;
        Json bad_ref = uvvis(spectrum(1, "wf"), f3);
        expect_only("spectrum: states_ref names a wavefunction",
                    excited_archive(two, Json::object(), false, {bad_ref}, f3),
                    "E_REF_UNRESOLVED");
    }
}

// --- wavefunction.bloch (spec § 5.3.1) ---------------------------------------

static const std::string kChain =
    R"({"atoms":[{"symbol":"H","position":[0.0,0.0,0.0],"atomic_number":1},)"
    R"({"symbol":"H","position":[0.8,0.0,0.0],"atomic_number":1}],)"
    R"("pbc":[true,false,false],)"
    R"("lattice_vectors":[[2.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,10.0]]})";
static const std::string kBlochBasis =
    R"({"structure_ref":"structure","pure":true,"n_ao":2,"shells":[)"
    R"({"center":0,"l":0,"exponents":[0.4166],"coefficients":[1.0]},)"
    R"({"center":1,"l":0,"exponents":[0.4166],"coefficients":[1.0]}]})";
static const std::string kBlochKpoints =
    R"({"kpoints":[[0.0,0.0,0.0],[0.5,0.0,0.0]],"weights":[0.5,0.5]})";
static const std::string kBlochMeta =
    R"({"spin":"unrestricted","n_spin":2,"n_mo":2,"n_ao":2,"complex":true,)"
    R"("k_point_convention":"fractional"})";

struct BlochSpec {
    std::vector<int64_t> coeff{2, 2, 2, 2, 2}, energies{2, 2, 2}, occupations{2, 2, 2};
    std::string kpoints = kBlochKpoints, meta = kBlochMeta, basis = kBlochBasis;
    std::string structure = kChain;
    Json fields = Json::object();
    bool structure_last = false;
    std::vector<Json> extra_sections;   // placed before the Bloch section
    std::vector<Entry> extra_files;
};

static std::vector<uint8_t> bloch_archive(const BlochSpec& b) {
    std::vector<Entry> files = {{"structure/structure.json", b.structure}};
    for (const auto& e : b.extra_files) files.push_back(e);
    Json members = Json::object();
    for (const auto& kv : {std::make_pair(std::string("basis"), b.basis),
                           std::make_pair(std::string("kpoints"), b.kpoints),
                           std::make_pair(std::string("mo_metadata"), b.meta)}) {
        members.set(kv.first, json_member("wf/" + kv.first + ".json", kv.second));
        files.push_back({"wf/" + kv.first + ".json", kv.second});
    }
    members.set("mo_coefficients", f64_member("wf/c.bin", b.coeff, files));
    members.set("energies", f64_member("wf/e.bin", b.energies, files));
    members.set("occupations", f64_member("wf/o.bin", b.occupations, files));
    Json wf = section("wf", "wavefunction.bloch", members);
    for (const auto& kv : b.fields.items()) wf.set(kv.first, kv.second);
    std::vector<Json> secs;
    if (!b.structure_last) secs.push_back(structure_section("structure",
                                                            "structure/structure.json",
                                                            b.structure));
    for (const auto& e : b.extra_sections) secs.push_back(e);
    secs.push_back(wf);
    if (b.structure_last) secs.push_back(structure_section("structure",
                                                           "structure/structure.json",
                                                           b.structure));
    return archive(manifest_of(secs), files);
}

// Expect only E_MEMBER_MALFORMED (or `code`), with one finding at `where`.
static void expect_bloch(const std::string& label, const BlochSpec& b, const std::string& where,
                         const std::string& code = "E_MEMBER_MALFORMED") {
    const std::vector<Finding> fs = semantic_findings(bloch_archive(b));
    std::set<std::string> got = codes(fs);
    got.erase("W_SCHEMA_SKIPPED");
    bool at = false;
    for (const auto& f : fs) at = at || (f.code == code && f.where == where);
    if (got != std::set<std::string>{code} || !at) {
        std::cerr << "FAIL: " << label << ": expected only " << code << " at " << where
                  << ", got" << describe(fs) << "\n";
        ++g_failures;
    }
}

static void test_wavefunction_bloch() {
    expect_clean("bloch: a valid section", bloch_archive(BlochSpec{}));
    {
        BlochSpec b;
        b.structure_last = true;
        expect_clean("bloch: the structure may come after the section", bloch_archive(b));
    }
    {
        BlochSpec b;
        b.kpoints = R"({"kpoints":[[0,0,0],[0.5,0,0]],"weights":[0.5,0.5000005]})";
        expect_clean("bloch: weights within 1e-6 of 1", bloch_archive(b));
        b.kpoints = R"({"kpoints":[[0,0,0],[0.5,0,0]]})";
        expect_clean("bloch: weights are optional", bloch_archive(b));
    }
    const std::string m = "/sections/1/members/";
    auto with = [](auto fn) {
        BlochSpec b;
        fn(b);
        return b;
    };
    expect_bloch("bloch: trailing axis of 1",
                 with([](BlochSpec& b) { b.coeff = {2, 2, 2, 2, 1}; }), m + "mo_coefficients");
    expect_bloch("bloch: n_spin 3", with([](BlochSpec& b) {
                     b.coeff = {2, 3, 2, 2, 2};
                     b.energies = b.occupations = {2, 3, 2};
                     b.meta = R"({"spin":"unrestricted","n_spin":2,"n_mo":2,"n_ao":2,)"
                              R"("complex":true,"k_point_convention":"fractional"})";
                 }),
                 m + "mo_coefficients");
    expect_bloch("bloch: rank-4 coefficients",
                 with([](BlochSpec& b) { b.coeff = {2, 2, 2, 2}; }), m + "mo_coefficients");
    expect_bloch("bloch: energies [2, 2, 3]",
                 with([](BlochSpec& b) { b.energies = {2, 2, 3}; }), m + "energies");
    expect_bloch("bloch: rank-2 energies",
                 with([](BlochSpec& b) { b.energies = {2, 2}; }), m + "energies");
    expect_bloch("bloch: occupations [1, 2, 2]",
                 with([](BlochSpec& b) { b.occupations = {1, 2, 2}; }), m + "occupations");
    expect_bloch("bloch: three k-points for two",
                 with([](BlochSpec& b) {
                     b.kpoints = R"({"kpoints":[[0,0,0],[0.5,0,0],[0.25,0,0]],)"
                                 R"("weights":[0.5,0.25,0.25]})";
                 }),
                 m + "kpoints");
    expect_bloch("bloch: one weight for two k-points",
                 with([](BlochSpec& b) {
                     b.kpoints = R"({"kpoints":[[0,0,0],[0.5,0,0]],"weights":[1.0]})";
                 }),
                 m + "kpoints");
    expect_bloch("bloch: weights sum to 0.75",
                 with([](BlochSpec& b) {
                     b.kpoints = R"({"kpoints":[[0,0,0],[0.5,0,0]],"weights":[0.5,0.25]})";
                 }),
                 m + "kpoints");
    expect_bloch("bloch: one label for two k-points",
                 with([](BlochSpec& b) {
                     b.kpoints = R"({"kpoints":[[0,0,0],[0.5,0,0]],"labels":["G"]})";
                 }),
                 m + "kpoints");
    const std::string head = R"({"complex":true,"k_point_convention":"fractional",)";
    expect_bloch("bloch: mo_metadata n_mo 3",
                 with([&](BlochSpec& b) {
                     b.meta = head + R"("spin":"unrestricted","n_spin":2,"n_mo":3,"n_ao":2})";
                 }),
                 m + "mo_metadata");
    expect_bloch("bloch: mo_metadata n_ao 3",
                 with([&](BlochSpec& b) {
                     b.meta = head + R"("spin":"unrestricted","n_spin":2,"n_mo":2,"n_ao":3})";
                 }),
                 m + "mo_metadata");
    expect_bloch("bloch: restricted with n_spin 2",
                 with([&](BlochSpec& b) {
                     b.meta = head + R"("spin":"restricted","n_spin":2,"n_mo":2,"n_ao":2})";
                 }),
                 m + "mo_metadata");
    expect_bloch("bloch: n_mo_per_k with one entry",
                 with([&](BlochSpec& b) {
                     b.meta = head + R"("spin":"unrestricted","n_spin":2,"n_mo":2,"n_ao":2,)"
                                     R"("n_mo_per_k":[2]})";
                 }),
                 m + "mo_metadata");
    expect_bloch("bloch: n_mo_per_k above n_mo",
                 with([&](BlochSpec& b) {
                     b.meta = head + R"("spin":"unrestricted","n_spin":2,"n_mo":2,"n_ao":2,)"
                                     R"("n_mo_per_k":[2,3]})";
                 }),
                 m + "mo_metadata");
    expect_bloch("bloch: basis n_ao 3",
                 with([](BlochSpec& b) {
                     b.basis = R"({"structure_ref":"structure","n_ao":3,"shells":[]})";
                 }),
                 m + "basis");
    expect_bloch("bloch: molecular structure",
                 with([](BlochSpec& b) { b.structure = kWater; }), m + "basis");
    expect_bloch("bloch: basis.structure_ref names nothing",
                 with([](BlochSpec& b) {
                     b.basis = R"({"structure_ref":"nowhere","n_ao":2,"shells":[]})";
                 }),
                 m + "basis", "E_REF_UNRESOLVED");
    {   // The section-level structure_ref takes precedence over the basis one.
        BlochSpec b;
        b.extra_sections.push_back(structure_section("mol", "mol/structure.json", kWater));
        b.extra_files.push_back({"mol/structure.json", kWater});
        b.basis = R"({"structure_ref":"mol","n_ao":2,"shells":[]})";
        b.fields.set("structure_ref", "structure");
        expect_clean("bloch: section-level periodic structure wins", bloch_archive(b));
        b.basis = kBlochBasis;
        b.fields.set("structure_ref", "mol");
        expect_bloch("bloch: section-level molecular structure", b, "/sections/2/structure_ref");
    }
}

// --- withdrawn version, chemical identity, job sets, reaction.scheme, preview --
// (spec § 3.1, § 5.1, § 5.8, § 5.14, § 5.15)

// The `where` of every finding with `code`.
static std::set<std::string> where_of(const std::vector<Finding>& fs, const std::string& code) {
    std::set<std::string> out;
    for (const auto& f : fs)
        if (f.code == code) out.insert(f.where);
    return out;
}

// Expect only `code` (as expect_only), every finding of it at `where`.
static void expect_at(const std::string& label, std::vector<uint8_t> bytes,
                      const std::string& code, const std::string& where) {
    const std::vector<Finding> fs = semantic_findings(bytes);
    CHECK(where_of(fs, code) == std::set<std::string>{where},
          label + ": " + code + " at " + where + describe(fs));
    expect_only(label, std::move(bytes), code);
}

static Json opaque_member(const std::string& path, const std::string& data,
                          const std::string& media_type = "") {
    Json m = Json::object();
    m.set("path", path).set("format", "binary").set("sha256", sha(data));
    if (!media_type.empty()) m.set("media_type", media_type);
    return m;
}

static void test_version_withdrawn() {
    Json v2 = Json::object();
    v2.set("qvf_version", 2);
    const std::vector<uint8_t> withdrawn = archive(manifest_of({structure_section()}, v2),
                                                   {{"structure/structure.json", kWater}});
    expect_only("version 2 is withdrawn", withdrawn, "W_VERSION_WITHDRAWN", Limits{},
                Finding::Severity::Warning);
    // The JSON Schema check sees version 1, so it reports nothing for it.
    const std::vector<Finding> with_schema = QvfReader::validate_bytes(withdrawn);
    CHECK(where_of(with_schema, "W_VERSION_WITHDRAWN") == std::set<std::string>{"/qvf_version"} &&
              where_of(with_schema, "E_SCHEMA").count("/qvf_version") == 0 &&
              !has_code(with_schema, "E_VERSION"),
          "version 2 is checked as 1, the schema check included" + describe(with_schema));
    Json v2f = Json::object();
    v2f.set("qvf_version", 2.0);
    expect_only("version 2.0 is withdrawn too",
                archive(manifest_of({structure_section()}, v2f),
                        {{"structure/structure.json", kWater}}),
                "W_VERSION_WITHDRAWN", Limits{}, Finding::Severity::Warning);
    for (const Json& bad : {Json(0), Json(3), Json(1.5), Json(true), Json("2")}) {
        Json v = Json::object();
        v.set("qvf_version", bad);
        expect_at("version " + bad.dump(0) + " is E_VERSION",
                  archive(manifest_of({structure_section()}, v),
                          {{"structure/structure.json", kWater}}),
                  "E_VERSION", "/qvf_version");
    }
}

// structure "structure" (kWater: O and H, sum Z = 9) with identity `fields`;
// `root` is merged into the manifest root, `second` adds another structure.
static std::vector<uint8_t> identity_archive(Json fields, Json root = Json::object(),
                                             bool second = false,
                                             const std::string& data = kWater) {
    Json s = structure_section("structure", "structure/structure.json", data);
    for (const auto& kv : fields.items()) s.set(kv.first, kv.second);
    std::vector<Json> sections{s};
    std::vector<Entry> files{{"structure/structure.json", data}};
    if (second) {
        sections.push_back(structure_section("other", "other/structure.json"));
        files.push_back({"other/structure.json", kWater});
    }
    return archive(manifest_of(sections, root), files);
}

static Json charge_mult(Json charge, Json mult) {
    Json f = Json::object();
    if (!charge.is_null()) f.set("charge", charge);
    if (!mult.is_null()) f.set("multiplicity", mult);
    return f;
}

static void test_chemical_identity() {
    Json named = charge_mult(-1, 1);
    named.set("name", "hydroxide").set("smiles", "[OH-]")
        .set("inchi", "InChI=1S/H2O/h1H2/p-1").set("inchikey", "XLYOFNOQVPJJNP-UHFFFAOYSA-M");
    expect_clean("identity: OH- with ten paired electrons", identity_archive(named));
    expect_clean("identity: OH radical, one unpaired electron",
                 identity_archive(charge_mult(0, 2)));
    expect_clean("identity: no electrons left, singlet", identity_archive(charge_mult(9, 1)));
    expect_clean("identity: multiplicity alone states no electron count",
                 identity_archive(charge_mult(Json(), 1)));
    expect_clean("identity: a charge alone is possible", identity_archive(charge_mult(0, Json())));
    expect_clean("identity: a non-integer charge is left to the schema",
                 identity_archive(charge_mult(0.0, 1)));
    expect_clean("identity: a multiplicity below 1 is left to the schema",
                 identity_archive(charge_mult(0, 0)));
    expect_at("identity: nine electrons cannot be a singlet",
              identity_archive(charge_mult(0, 1)), "E_CHARGE_MULTIPLICITY",
              "/sections/0/multiplicity");
    expect_at("identity: unpaired electrons beyond the count",
              identity_archive(charge_mult(-1, 13)), "E_CHARGE_MULTIPLICITY",
              "/sections/0/multiplicity");
    expect_at("identity: a negative electron count", identity_archive(charge_mult(10, 1)),
              "E_CHARGE_MULTIPLICITY", "/sections/0/charge");
    // Exact for a charge near the int64 limit: N = 9 + 2^63 - 1 is even.
    const long long huge = -9223372036854775807LL;
    expect_clean("identity: exact electron count for a huge negative charge",
                 identity_archive(charge_mult(Json(huge), 1)));
    expect_at("identity: parity of a huge electron count",
              identity_archive(charge_mult(Json(huge), 2)), "E_CHARGE_MULTIPLICITY",
              "/sections/0/multiplicity");
    // Atoms without a usable atomic number: the electron count is unknown.
    const std::string no_z =
        R"({"atoms":[{"symbol":"O","position":[0,0,0]}],"pbc":[false,false,false]})";
    expect_clean("identity: missing atomic_number skips the rule",
                 identity_archive(charge_mult(5, 1), Json::object(), false, no_z));
    const std::string negative_z =
        R"({"atoms":[{"symbol":"X","position":[0,0,0],"atomic_number":-1}],)"
        R"("pbc":[false,false,false]})";
    expect_clean("identity: a negative atomic_number skips the rule",
                 identity_archive(charge_mult(5, 1), Json::object(), false, negative_z));

    // Agreement with provenance in an archive with one structure section.
    Json prov = Json::object(), root = Json::object();
    prov.set("charge", 1).set("multiplicity", 3);
    root.set("provenance", prov);
    const std::vector<Finding> both = semantic_findings(identity_archive(charge_mult(0, 2), root));
    CHECK(codes(both) == (std::set<std::string>{"E_ROOT_INCONSISTENT", "W_SCHEMA_SKIPPED"}) &&
              where_of(both, "E_ROOT_INCONSISTENT") ==
                  (std::set<std::string>{"/sections/0/charge", "/sections/0/multiplicity"}),
          "identity: charge and multiplicity differ from provenance" + describe(both));
    expect_clean("identity: provenance agrees", identity_archive(charge_mult(1, 3), root));
    expect_clean("identity: provenance describes no one structure of two",
                 identity_archive(charge_mult(0, 2), root, true));
    Json loose = Json::object(), loose_root = Json::object();
    loose.set("charge", 1.0).set("multiplicity", "3");
    loose_root.set("provenance", loose);
    expect_clean("identity: a non-integer provenance value is not compared",
                 identity_archive(charge_mult(0, 2), loose_root));
}

// structure + job.spec "job" + run.record "run" with the section fields `run`.
static std::vector<uint8_t> job_set_archive(const Json& run_fields) {
    const std::string spec = R"({"job_type":"molecular","structure_ref":"structure"})";
    Json job_members = Json::object();
    job_members.set("spec", json_member("job/spec.json", spec));
    Json job = section("job", "job.spec", job_members);
    job.set("structure_ref", "structure");
    Json run_members = Json::object();
    run_members.set("log", opaque_member("run/log.txt", "converged\n"));
    Json run = section("run", "run.record", run_members);
    run.set("program", "t");
    for (const auto& kv : run_fields.items()) run.set(kv.first, kv.second);
    return archive(manifest_of({structure_section(), job, run}),
                   {{"structure/structure.json", kWater}, {"job/spec.json", spec},
                    {"run/log.txt", "converged\n"}});
}

static void test_job_sets() {
    Json ok = Json::object();
    ok.set("job_ref", "job").set("run_status", "converged");
    expect_clean("job set: job_ref names the job.spec", job_set_archive(ok));
    for (const Json& bad : {Json("structure"), Json("nowhere"), Json(5), Json(nullptr)}) {
        Json fields = Json::object();
        fields.set("job_ref", bad);
        expect_at("job set: job_ref " + bad.dump(0), job_set_archive(fields),
                  "E_REF_UNRESOLVED", "/sections/2/job_ref");
    }
}

// Structures "a" and "b" and a reaction.scheme "scheme" with payload `scheme`.
static std::vector<uint8_t> scheme_archive(const std::string& scheme) {
    Json members = Json::object();
    members.set("scheme", json_member("scheme/scheme.json", scheme));
    return archive(manifest_of({structure_section("a", "a/structure.json"),
                                structure_section("b", "b/structure.json"),
                                section("scheme", "reaction.scheme", members)}),
                   {{"a/structure.json", kWater},
                    {"b/structure.json", kWater},
                    {"scheme/scheme.json", scheme}});
}

static void test_reaction_scheme() {
    const std::string m = "/sections/2/members/scheme";
    expect_clean("scheme: two steps, a structure in two roles of a step",
                 scheme_archive(R"({"steps":[
                     {"id":"r1","reactants":[{"structure_ref":"a"}],
                      "products":[{"structure_ref":"b"}],
                      "reagents":[{"structure_ref":"a"}],"arrow":"equilibrium"},
                     {"id":"r2","reactants":[{"structure_ref":"b"}],
                      "products":[{"structure_ref":"a","coefficient":2}]}]})"));
    expect_at("scheme: a participant names no section",
              scheme_archive(R"({"steps":[{"id":"r1","reactants":[{"structure_ref":"a"}],
                                 "products":[{"structure_ref":"nowhere"}]}]})"),
              "E_REF_UNRESOLVED", m);
    expect_at("scheme: a participant names the scheme itself",
              scheme_archive(R"({"steps":[{"id":"r1","reactants":[{"structure_ref":"scheme"}],
                                 "products":[{"structure_ref":"b"}]}]})"),
              "E_REF_UNRESOLVED", m);
    expect_at("scheme: duplicate step ids",
              scheme_archive(R"({"steps":[
                  {"id":"r1","reactants":[{"structure_ref":"a"}],
                   "products":[{"structure_ref":"b"}]},
                  {"id":"r1","reactants":[{"structure_ref":"b"}],
                   "products":[{"structure_ref":"a"}]}]})"),
              "E_MEMBER_MALFORMED", m);
    expect_at("scheme: a structure twice among the reactants",
              scheme_archive(R"({"steps":[{"id":"r1",
                  "reactants":[{"structure_ref":"a"},{"structure_ref":"a"}],
                  "products":[{"structure_ref":"b"}]}]})"),
              "E_MEMBER_MALFORMED", m);
    const std::vector<Finding> fs = semantic_findings(
        scheme_archive(R"({"steps":[{"id":"r1",
            "reactants":[{"structure_ref":"x"},{"structure_ref":"x"}],
            "products":[{"structure_ref":"b"}]}]})"));
    size_t unresolved = 0;
    for (const auto& f : fs) unresolved += f.code == "E_REF_UNRESOLVED";
    CHECK(unresolved == 2 && has_code(fs, "E_MEMBER_MALFORMED") &&
              codes(fs) == (std::set<std::string>{"E_REF_UNRESOLVED", "E_MEMBER_MALFORMED",
                                                  "W_SCHEMA_SKIPPED"}),
          "scheme: each unresolved occurrence and the repeat are reported" + describe(fs));
    expect_clean("scheme: malformed steps are left to the schema",
                 scheme_archive(R"({"steps":[7,{"id":5,"reactants":"a"}]})"));
}

static const std::string kPng = std::string("\x89PNG\r\n\x1a\n", 8) + "IHDR";
static const std::string kJpeg = "\xff\xd8\xff\xe0JFIF";
static const std::string kSvg =
    R"(<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg"/>)";

// structure + a preview "preview" with the given members (role -> bytes, type).
struct PreviewPart {
    std::string role, data, media_type;
};
static std::vector<uint8_t> preview_archive(const std::vector<PreviewPart>& parts,
                                            const std::string& kind = "preview") {
    Json members = Json::object();
    std::vector<Entry> files{{"structure/structure.json", kWater}};
    for (const auto& p : parts) {
        members.set(p.role, opaque_member("preview/" + p.role, p.data, p.media_type));
        files.push_back({"preview/" + p.role, p.data, false});
    }
    return archive(manifest_of({structure_section(), section("preview", kind, members)}),
                   files);
}

static void test_preview_members() {
    expect_clean("preview: PNG thumbnail and SVG image",
                 preview_archive({{"thumbnail", kPng, "image/png"},
                                  {"image", kSvg, "image/svg+xml"}}));
    expect_clean("preview: JPEG thumbnail and PNG image",
                 preview_archive({{"thumbnail", kJpeg, "image/jpeg"},
                                  {"image", kPng, "image/png"}}));
    expect_at("preview: image/png with JPEG bytes",
              preview_archive({{"thumbnail", kJpeg, "image/png"}}), "E_MEMBER_MALFORMED",
              "/sections/1/members/thumbnail");
    expect_at("preview: image/jpeg with PNG bytes",
              preview_archive({{"thumbnail", kPng, "image/jpeg"}}), "E_MEMBER_MALFORMED",
              "/sections/1/members/thumbnail");
    expect_at("preview: a truncated PNG signature",
              preview_archive({{"thumbnail", kPng.substr(0, 7), "image/png"}}),
              "E_MEMBER_MALFORMED", "/sections/1/members/thumbnail");
    expect_at("preview: SVG that is not UTF-8",
              preview_archive({{"thumbnail", kPng, "image/png"},
                               {"image", "<svg>\xff</svg>", "image/svg+xml"}}),
              "E_UTF8", "/sections/1/members/image");
    expect_at("preview: UTF-8 without an <svg element",
              preview_archive({{"thumbnail", kPng, "image/png"},
                               {"image", "<html/>", "image/svg+xml"}}),
              "E_MEMBER_MALFORMED", "/sections/1/members/image");
    // The media_type itself is the schema's business.
    expect_clean("preview: a missing media_type is left to the schema",
                 preview_archive({{"thumbnail", kPng, ""}}));
    expect_clean("preview: another media_type is left to the schema",
                 preview_archive({{"thumbnail", kJpeg, "image/gif"}}));
    const std::vector<Finding> fs =
        QvfReader::validate_bytes(preview_archive({{"thumbnail", kPng, ""}}));
    CHECK(where_of(fs, "E_SCHEMA").count("/sections/1"),
          "preview: the schema check requires the media_type" + describe(fs));
    // Opaque members of a vendor section are legal, with or without a
    // media_type, and their bytes are not judged.
    expect_clean("preview: a vendor section's opaque members",
                 preview_archive({{"drawing", kJpeg, "image/png"}, {"blob", "\x01\x02", ""}},
                                 "x_test.document"));
}

// --- toolkit 0.3.0: what the family's programs write or read ------------------

static const std::string kH2 =
    R"({"atoms":[{"symbol":"H","position":[0.0,0.0,0.0],"atomic_number":1},)"
    R"({"symbol":"H","position":[0.0,0.0,0.74],"atomic_number":1}],)"
    R"("pbc":[false,false,false],"lattice_vectors":null})";
// Water with all three atoms (kWater holds O and one H).
static const std::string kWater3 =
    R"({"atoms":[{"symbol":"O","position":[0.0,0.0,0.1173],"atomic_number":8},)"
    R"({"symbol":"H","position":[0.0,0.7572,-0.4692],"atomic_number":1},)"
    R"({"symbol":"H","position":[0.0,-0.7572,-0.4692],"atomic_number":1}],)"
    R"("pbc":[false,false,false],"lattice_vectors":null})";

// Water, H2, a job.spec "job" and a bonds section with payload `bonds`
// (structure_ref `ref` when not empty); with `bonds_first` the bonds section
// comes first and there is no job.
static std::vector<uint8_t> bonds_archive(const std::string& bonds, const std::string& ref,
                                          bool bonds_first = false) {
    Json bm = Json::object();
    bm.set("bonds", json_member("bonds/connectivity.json", bonds));
    Json b = section("bonds", "bonds", bm);
    if (!ref.empty()) b.set("structure_ref", ref);
    Json jm = Json::object();
    jm.set("spec", json_member("job/spec.json", R"({"job_type":"molecular"})"));
    const Json water = structure_section("structure_water", "water/structure.json", kWater3);
    const Json h2 = structure_section("structure_h2", "h2/structure.json", kH2);
    std::vector<Entry> files{{"water/structure.json", kWater3},
                             {"h2/structure.json", kH2},
                             {"bonds/connectivity.json", bonds}};
    if (bonds_first) return archive(manifest_of({b, water, h2}), files);
    files.push_back({"job/spec.json", R"({"job_type":"molecular"})"});
    return archive(manifest_of({water, h2, section("job", "job.spec", jm), b}), files);
}

// A water structure and an H2 wavefunction.gto with mo_metadata `meta`.
static std::vector<uint8_t> localized_archive(const std::string& meta) {
    const std::string basis =
        R"({"structure_ref":"structure","pure":true,"n_ao":2,"shells":[)"
        R"({"center":1,"l":0,"exponents":[1.24],"coefficients":[1.0]},)"
        R"({"center":2,"l":0,"exponents":[1.24],"coefficients":[1.0]}]})";
    std::vector<Entry> files{{"structure/structure.json", kWater}};
    Json m = Json::object();
    m.set("basis", json_member("wf/basis.json", basis));
    m.set("mo_metadata", json_member("wf/mo_metadata.json", meta));
    m.set("mo_coefficients", f64_values("wf/c.dat", {2, 2}, {1.0, 0.0, 0.0, 1.0}, files));
    files.push_back({"wf/basis.json", basis});
    files.push_back({"wf/mo_metadata.json", meta});
    return archive(manifest_of({structure_section(), section("wf", "wavefunction.gto", m)}),
                   files);
}

// A water structure and a job.spec with payload `spec`.
static std::vector<uint8_t> job_spec_archive(const std::string& spec) {
    Json m = Json::object();
    m.set("spec", json_member("job/spec.json", spec));
    return archive(manifest_of({structure_section(), section("job", "job.spec", m)}),
                   {{"structure/structure.json", kWater}, {"job/spec.json", spec}});
}

// Every finding of a validation with the generated JSON Schema check.
static std::vector<Finding> schema_findings(std::vector<uint8_t> bytes) {
    return QvfReader::validate_bytes(std::move(bytes));
}

static void test_family_drift() {
    // structure_ref on a standalone bonds section (spec § 5.1).
    const std::string two = R"({"pairs":[{"i":0,"j":1},{"i":0,"j":2}]})";
    expect_clean("bonds: indices within the referenced water", bonds_archive(two, "structure_water"));
    expect_at("bonds: index 2 of a two-atom structure", bonds_archive(two, "structure_h2"),
              "E_MEMBER_MALFORMED", "/sections/3/members/bonds");
    expect_at("bonds: the referenced structure comes later",
              bonds_archive(R"({"pairs":[{"i":1,"j":2}]})", "structure_h2", true),
              "E_MEMBER_MALFORMED", "/sections/0/members/bonds");
    expect_clean("bonds: without structure_ref no atom count is checked",
                 bonds_archive(R"({"pairs":[{"i":0,"j":7}]})", ""));
    for (const char* ref : {"missing", "job"})
        expect_at(std::string("bonds: structure_ref ") + ref,
                  bonds_archive(R"({"pairs":[{"i":0,"j":7}]})", ref), "E_REF_UNRESOLVED",
                  "/sections/3/structure_ref");
    {
        const std::vector<Finding> fs = semantic_findings(bonds_archive(two, "structure_h2"));
        bool named = false;
        for (const auto& f : fs)
            named = named || (f.message.find("pairs[1].j 2") != std::string::npos &&
                              f.message.find("'structure_h2'") != std::string::npos);
        CHECK(named, "bonds: the message names the pair and the structure" + describe(fs));
    }

    // mo_metadata: the IAO reference basis and the deprecated keys (spec § 5.3).
    const std::string mo = R"({"n_mo":2,"n_ao":2,"spin":"restricted",)"
                           R"("energies":[-0.58,0.67],"occupations":[2.0,0.0],)";
    for (const std::string& extra :
         {std::string(R"("orbital_kind":"ibo","localization_method":"ibo",)"
                      R"("localization_reference_basis":"minao"})"),
          std::string(R"("orbital_kind":"localized","localization_reference_basis":"mini",)"
                      R"("iao_charges":[0.0,0.1,-0.1]})"),
          std::string(R"("orbital_kind":"natural","occupation_semantics":"transition_weight",)"
                      R"("excitation_index":1,"excitation_energy_ev":10.2})")}) {
        const std::vector<Finding> fs = schema_findings(localized_archive(mo + extra));
        CHECK(fs.empty(), "mo_metadata " + extra + " is schema-clean" + describe(fs));
    }
    for (const std::string& extra :
         {std::string(R"("localization_reference_basis":""})"),
          std::string(R"("iao_charges":["0.1"]})"), std::string(R"("excitation_index":-1})"),
          std::string(R"("excitation_energy_ev":"10.2"})")}) {
        const std::vector<Finding> fs = schema_findings(localized_archive(mo + extra));
        CHECK(where_of(fs, "E_PAYLOAD_SCHEMA") ==
                      std::set<std::string>{"/sections/1/members/mo_metadata"} &&
                  codes(fs) == std::set<std::string>{"E_PAYLOAD_SCHEMA"},
              "mo_metadata " + extra + " is E_PAYLOAD_SCHEMA" + describe(fs));
    }

    // job.spec functional and basis: null, or the deprecated "" (spec § 5.9).
    for (const char* spec :
         {R"({"job_type":"molecular","method":"rhf","functional":null,"basis":"sto-3g"})",
          R"({"job_type":"molecular","method":"gfn2_xtb","functional":null,"basis":null})",
          R"({"job_type":"molecular","method":"gfn2_xtb","basis":""})"}) {
        const std::vector<Finding> fs = schema_findings(job_spec_archive(spec));
        CHECK(fs.empty(), std::string("job.spec ") + spec + " is schema-clean" + describe(fs));
    }
    for (const char* spec : {R"({"job_type":"molecular","functional":3})",
                             R"({"job_type":"molecular","basis":["sto-3g"]})"}) {
        const std::vector<Finding> fs = schema_findings(job_spec_archive(spec));
        CHECK(where_of(fs, "E_PAYLOAD_SCHEMA") ==
                  std::set<std::string>{"/sections/1/members/spec"},
              std::string("job.spec ") + spec + " is E_PAYLOAD_SCHEMA" + describe(fs));
    }

    // viewer_defaults camera, crossfade and energy_window never fail (spec § 3.4).
    for (const char* hints :
         {R"({"camera":{"position":[0,-12,4],"focal_point":[0,0,0],"view_up":[0,0,1],)"
          R"("view_angle":30},"crossfade":{"a":"a","b":"b","blend":0.25},)"
          R"("structure":{"energy_window":[-15,10]}})",
          R"({"camera":"front","crossfade":5,"structure":{"energy_window":"valence"}})",
          R"({"camera":null,"crossfade":[],"structure":{"energy_window":[10,-15]}})"}) {
        Json root = Json::object();
        root.set("viewer_defaults", Json::parse(hints));
        const std::vector<Finding> fs = schema_findings(
            archive(manifest_of({structure_section()}, root), {{"structure/structure.json", kWater}}));
        CHECK(fs.empty(), std::string("viewer hints ") + hints + " are clean" + describe(fs));
    }
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--negative") return run_negative_corpus(argv[2]);
    if (argc >= 3 && std::string(argv[1]) == "--model-roundtrip")
        return run_model_roundtrip(std::vector<std::string>(argv + 2, argv + argc));
    if (argc != 1) {
        std::cerr << "usage: qvf_reader_tests [--negative DIR | --model-roundtrip DIR...]\n";
        return 2;
    }
    test_json_depth();
    test_inflate_bound();
    test_good_archive();
    test_container_rules();
    test_size_updates_and_byte_order();
    test_json_caps();
    test_manifest_and_sections();
    test_members();
    test_bands_spin_axis();
    test_unit_names();
    test_derivatives_and_energies();
    test_excited_states();
    test_trajectory_members();
    test_wavefunction_rules();
    test_volumes();
    test_wavefunction_bloch();
    test_version_withdrawn();
    test_chemical_identity();
    test_job_sets();
    test_reaction_scheme();
    test_preview_members();
    test_family_drift();

    test_section_references();
    test_schema_checker();
    test_model();
    if (g_failures == 0) {
        std::cout << "all C++ reader tests passed\n";
        return 0;
    }
    std::cerr << g_failures << " C++ reader test(s) failed\n";
    return 1;
}
