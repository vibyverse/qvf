// test_container.cpp — tests for the container layer (qvf/container.hpp) and
// the streaming / updating QvfWriter entry points built on it (WP-40).
//
//   * the codecs: deflate output unchanged (golden digests of the encoder
//     before streaming), chunk-independence, the block framing, streaming
//     inflate, inflate hardening;
//   * 32-bit archives byte-identical to the original writer (golden archive);
//   * ZIP64: forced on small members, a 70 000-entry archive, offsets past
//     4 GiB (through a sink that starts there), all read back;
//   * data descriptors on non-seekable sinks, read back;
//   * read_into() and the caps enforced while inflating; CRC-32 checks;
//   * QvfWriter: to_bytes() == write(), streaming vs buffered members,
//     begin_binary(), open_for_update() then validate then repack(), the
//     orphan warning gone after repack, superseded entries.
//
// Scratch files go to the current directory (CTest runs in the build tree).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "qvf/container.hpp"
#include "qvf/qvf.hpp"
#include "qvf/qvf_reader.hpp"
#include "qvf_crc32.hpp"
#include "qvf_deflate.hpp"
#include "qvf_inflate.hpp"
#include "qvf_sha256.hpp"
#include "qvf_zip.hpp"

using namespace qvf;
namespace fs = std::filesystem;
using container::Method;

static int g_failures = 0;

#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::cerr << "FAIL: " << (msg) << " (" << #cond << ")\n";      \
            ++g_failures;                                                  \
        }                                                                  \
    } while (0)

static bool throws(const std::function<void()>& fn, const std::string& needle = "") {
    try {
        fn();
    } catch (const std::exception& e) {
        return needle.empty() || std::string(e.what()).find(needle) != std::string::npos;
    }
    return false;
}

static std::string sha(const std::vector<uint8_t>& v) { return detail::Sha256::hex_of(v); }

static std::vector<uint8_t> bytes_of(const std::string& s) {
    return std::vector<uint8_t>(s.begin(), s.end());
}

static std::vector<uint8_t> read_file(const std::string& path) {
    std::ifstream fh(path, std::ios::binary);
    return std::vector<uint8_t>((std::istreambuf_iterator<char>(fh)),
                                std::istreambuf_iterator<char>());
}

static void write_file(const std::string& path, const std::vector<uint8_t>& b) {
    std::ofstream fh(path, std::ios::binary);
    fh.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
}

static uint32_t rd32(const std::vector<uint8_t>& z, size_t o) {
    return static_cast<uint32_t>(z[o]) | (static_cast<uint32_t>(z[o + 1]) << 8) |
           (static_cast<uint32_t>(z[o + 2]) << 16) | (static_cast<uint32_t>(z[o + 3]) << 24);
}
static uint16_t rd16(const std::vector<uint8_t>& z, size_t o) {
    return static_cast<uint16_t>(z[o] | (z[o + 1] << 8));
}
static uint64_t rd64(const std::vector<uint8_t>& z, size_t o) {
    return static_cast<uint64_t>(rd32(z, o)) | (static_cast<uint64_t>(rd32(z, o + 4)) << 32);
}

static size_t count_sig(const std::vector<uint8_t>& z, uint32_t sig) {
    size_t n = 0;
    for (size_t p = 0; p + 4 <= z.size(); ++p) n += rd32(z, p) == sig;
    return n;
}

static std::string describe(const std::vector<Finding>& found) {
    std::string s;
    for (const auto& f : found) s += "\n    " + f.code + ": " + f.message + " [" + f.where + "]";
    return s;
}

// A clean archive: no error and no finding but W_SCHEMA_SKIPPED.
static bool clean(const std::vector<Finding>& found) {
    for (const auto& f : found)
        if (f.code != "W_SCHEMA_SKIPPED") return false;
    return true;
}

// --- deterministic inputs (the same as the golden generator used) -----------

static uint32_t g_lcg;
static uint32_t lcg() {
    g_lcg = g_lcg * 1664525u + 1013904223u;
    return g_lcg;
}

static std::vector<uint8_t> input(int k) {
    std::vector<uint8_t> v;
    switch (k) {
        case 0: break;
        case 1: v = {'a'}; break;
        case 2: v = {'a', 'b', 'c'}; break;
        case 3:
            for (int i = 0; i < 1000; ++i) {
                const char* s = "the quick brown fox ";
                v.insert(v.end(), s, s + 20);
            }
            break;
        case 4: v.assign(4096, 0); break;
        case 5:
            g_lcg = 1;
            for (int i = 0; i < 70000; ++i) v.push_back(static_cast<uint8_t>(lcg() >> 24));
            break;
        case 6: {
            g_lcg = 7;
            for (int i = 0; i < 40000; ++i)
                v.push_back(static_cast<uint8_t>('a' + (lcg() >> 24) % 7));
            const size_t dists[] = {32768, 32769, 32767, 1, 258, 30000};
            for (int r = 0; r < 60; ++r) {
                const size_t d = dists[r % 6];
                const size_t len = 100 + (lcg() >> 20) % 900;
                for (size_t i = 0; i < len; ++i) v.push_back(v[v.size() - d]);
                for (int i = 0; i < 500; ++i) v.push_back(static_cast<uint8_t>(lcg() >> 24));
            }
            break;
        }
        case 7:
            // Doubles in [-1, 1) from the generator, scaled by a power of two: exact
            // on every platform (a libm sine differs in its last bits between
            // implementations, which changed the input and so the digest).
            for (int i = 0; i < 100000; ++i) {
                const uint64_t hi = lcg(), lo = lcg();   // 53 random mantissa bits
                const uint64_t n = ((hi << 32) | lo) >> 11;
                const double x = static_cast<double>(n) * (1.0 / 4503599627370496.0) - 1.0;
                uint64_t b;
                std::memcpy(&b, &x, 8);
                for (int j = 0; j < 8; ++j) v.push_back(static_cast<uint8_t>(b >> (8 * j)));
            }
            break;
    }
    return v;
}

// Digests of deflate(input(k)) from the encoder as it was before it streamed
// (input 7 re-recorded when its generation stopped depending on libm).
struct Golden {
    int k;
    size_t in_size, out_size;
    const char* sha256;
};
static const Golden kGolden[] = {
    {0, 0, 2, "9b4fb24edd6d1d8830e272398263cdbf026b97392cc35387b991dc0248a628f9"},
    {1, 1, 3, "59fced9a2ff83dabdd7055cf3abd6c12b802e6397b6fba5b8251cdf95d2bc991"},
    {2, 3, 5, "7ce573d21753aedbe7d6d1f5893e3f1265a35b9b51e42a154174a4ce9602ae39"},
    {3, 20000, 178, "40411c2d66f03f3e27e64118a427666f608271ef06d3fefd4316a312fb9d8d7f"},
    {4, 4096, 29, "e383e14e1349cb4f3ca96f558f50de4ab2b71ad3bfe055d855faa2855f633779"},
    {5, 70000, 73812, "dcdaa4067c01f327e13a509d4434b4f2d8b89031c4118a7cf6dbbd0ca986617f"},
    {6, 102566, 58787, "e9944a42b3b91bd21cc4165ab523976eddbbc2e969460c932e077dcc1fc11666"},
    {7, 800000, 838784, "be2f2837bb48fbebd69a7d7070741e74c16a424d165862f044de2d5153e495ec"},
};

// --- codecs -------------------------------------------------------------------

static void test_deflate_golden_and_chunking() {
    for (const Golden& g : kGolden) {
        const std::vector<uint8_t> in = input(g.k);
        CHECK(in.size() == g.in_size, "golden input size " + std::to_string(g.k));
        const std::vector<uint8_t> out = detail::deflate(in);
        CHECK(out.size() == g.out_size && sha(out) == g.sha256,
              "deflate output unchanged, input " + std::to_string(g.k) + " (got " +
                  std::to_string(out.size()) + " bytes, sha256 " + sha(out) + ")");
        // Any chunking gives the same bytes.
        for (size_t chunk : {size_t(1), size_t(7), size_t(4096), size_t(65537)}) {
            if (chunk == 1 && in.size() > 30000) continue;   // keep the test quick
            detail::Deflater d;
            std::vector<uint8_t> s;
            for (size_t p = 0; p < in.size(); p += chunk)
                d.write(in.data() + p, std::min(chunk, in.size() - p), s);
            d.finish(s);
            CHECK(s == out, "chunked deflate == one-shot, input " + std::to_string(g.k) +
                                " chunk " + std::to_string(chunk));
        }
        // The block framing round-trips and never grows data much.
        detail::Deflater b(detail::Deflater::Framing::Blocks);
        std::vector<uint8_t> bs;
        for (size_t p = 0; p < in.size(); p += 10000)
            b.write(in.data() + p, std::min<size_t>(10000, in.size() - p), bs);
        b.finish(bs);
        CHECK(detail::inflate(bs, in.size()) == in, "block framing round-trips " +
                                                        std::to_string(g.k));
        CHECK(bs.size() <= in.size() + 5 * (in.size() / 32768 + 2),
              "block framing bounded growth " + std::to_string(g.k));
    }
}

static void test_inflate_streaming() {
    const std::vector<uint8_t> in = input(6);
    const std::vector<uint8_t> comp = detail::deflate(in);
    // Input pulled in odd-sized pieces, output pushed through a small window.
    size_t pos = 0;
    const size_t piece = 333;
    detail::InflateBitReader br(nullptr, 0, [&](const uint8_t*& q) -> size_t {
        if (pos >= comp.size()) return 0;
        const size_t n = std::min(piece, comp.size() - pos);
        q = comp.data() + pos;
        pos += n;
        return n;
    });
    std::vector<uint8_t> out;
    size_t chunks = 0;
    detail::InflateWindow win(in.size(), [&](const uint8_t* p, size_t n) {
        out.insert(out.end(), p, p + n);
        ++chunks;
    });
    detail::inflate_stream(br, win);
    CHECK(out == in, "streaming inflate == input");
    // A large output is delivered in several chunks.
    const std::vector<uint8_t> zeros(3u << 20, 0);
    const std::vector<uint8_t> zc = detail::deflate(zeros);
    detail::InflateBitReader br2(zc.data(), zc.size());
    size_t got = 0;
    chunks = 0;
    detail::InflateWindow win2(zeros.size(), [&](const uint8_t*, size_t n) {
        got += n;
        ++chunks;
    });
    detail::inflate_stream(br2, win2);
    CHECK(got == zeros.size() && chunks >= 2, "large output arrives in chunks");
    // Output past the bound is refused.
    CHECK(throws([&] { detail::inflate(zc.data(), zc.size(), 0, 1000); }, "exceeds"),
          "inflate past max_out throws");
    // Stored block whose NLEN is not the complement of LEN.
    const std::vector<uint8_t> bad_stored = {0x01, 0x03, 0x00, 0x00, 0x00, 'a', 'b', 'c'};
    CHECK(throws([&] { detail::inflate(bad_stored); }, "stored block"),
          "stored block length check");
    const std::vector<uint8_t> good_stored = {0x01, 0x03, 0x00, 0xFC, 0xFF, 'a', 'b', 'c'};
    CHECK(detail::inflate(good_stored) == bytes_of("abc"), "stored block decodes");
    // CRC-32 known answer through the incremental class.
    detail::Crc32 c;
    c.update(reinterpret_cast<const uint8_t*>("1234"), 4);
    c.update(reinterpret_cast<const uint8_t*>("56789"), 5);
    CHECK(c.value() == 0xCBF43926u, "incremental crc32");
}

// --- archives -------------------------------------------------------------------

static std::vector<uint8_t> golden_entries_rand() {
    g_lcg = 3;
    std::vector<uint8_t> r;
    for (int i = 0; i < 1000; ++i) r.push_back(static_cast<uint8_t>(lcg() >> 24));
    return r;
}

static void test_golden_archive() {
    const std::string manifest = "{\"qvf_version\": 1}";
    const std::string json = "{\"atoms\": [1, 2, 3], \"name\": \"water water water water\"}";
    const std::vector<uint8_t> zeros(4096, 0);
    const std::vector<uint8_t> rnd = golden_entries_rand();
    // The original interface.
    detail::ZipWriter old;
    old.add("manifest.json", manifest, false);
    old.add("a/b.json", json);
    old.add("x/zeros.bin", zeros);
    old.add("x/rand.bin", rnd);
    old.add("empty", std::string());
    const std::vector<uint8_t> a = old.finish();
    CHECK(a.size() == 1583 &&
              sha(a) == "3ea65fa5d3f93bb43d49f950d597daf339e5d74bb1fd325c428df83f31229d59",
          "32-bit archive byte-identical to the original writer");
    // The container writer, whole entries.
    container::MemorySink* mem = new container::MemorySink;
    container::ZipWriter zw{std::unique_ptr<container::Sink>(mem)};
    zw.write("manifest.json", manifest, Method::Store);
    zw.write("a/b.json", json);
    zw.write("x/zeros.bin", zeros);
    zw.write("x/rand.bin", rnd);
    zw.write("empty", std::string());
    zw.finish();
    CHECK(mem->bytes() == a, "container::ZipWriter::write == original writer");
    // Streamed STORE entries on a seekable sink (patched local headers) give
    // the same bytes as whole STORE entries.
    container::MemorySink* m1 = new container::MemorySink;
    container::MemorySink* m2 = new container::MemorySink;
    container::ZipWriter w1{std::unique_ptr<container::Sink>(m1)};
    container::ZipWriter w2{std::unique_ptr<container::Sink>(m2)};
    w1.write("x/rand.bin", rnd, Method::Store);
    container::EntrySink es = w2.begin("x/rand.bin", Method::Store);
    es.write(rnd.data(), 10);
    es.write(rnd.data() + 10, rnd.size() - 10);
    container::EntryInfo info = es.close();
    w1.finish();
    w2.finish();
    CHECK(m1->bytes() == m2->bytes(), "streamed STORE == whole STORE, byte for byte");
    CHECK(info.sha256 == sha(rnd) && info.size == rnd.size() &&
              info.crc32 == detail::crc32_of(rnd.data(), rnd.size()),
          "EntrySink computes sha256, size, crc32");
}

static void test_zip64_forced() {
    container::MemorySink* mem = new container::MemorySink;
    container::ZipWriter zw{std::unique_ptr<container::Sink>(mem), /*force_zip64=*/true};
    const std::vector<uint8_t> six = input(6);
    zw.write("a.txt", std::string("hello zip64"), Method::Store);
    zw.write("b.bin", six);
    container::EntrySink s = zw.begin("c.bin", Method::Deflate);
    s.write(six);
    s.close();
    zw.finish();
    const std::vector<uint8_t>& z = mem->bytes();
    CHECK(count_sig(z, 0x06064b50u) == 1 && count_sig(z, 0x07064b50u) == 1,
          "forced ZIP64: end record and locator");
    const size_t eocd = z.size() - 22;
    CHECK(rd32(z, eocd) == 0x06054b50u && rd16(z, eocd + 10) == 0xFFFF &&
              rd32(z, eocd + 12) == 0xFFFFFFFFu && rd32(z, eocd + 16) == 0xFFFFFFFFu,
          "forced ZIP64: end record saturated");
    container::ZipReader r(z);
    CHECK(r.entries().size() == 3, "forced ZIP64 read: 3 entries");
    for (const auto& e : r.entries()) {
        CHECK(e.zip64, "central record carries the ZIP64 extra: " + e.name);
        CHECK(r.local_header(e).zip64, "local header carries the ZIP64 extra: " + e.name);
    }
    CHECK(r.read("a.txt") == bytes_of("hello zip64"), "forced ZIP64 read a.txt");
    CHECK(r.read("b.bin") == six && r.read("c.bin") == six, "forced ZIP64 read b, c");
    CHECK(r.local_header(*r.find("b.bin")).usize == six.size(),
          "local ZIP64 extra holds the size");

    // A per-entry force on an otherwise 32-bit archive: no end record needed.
    container::MemorySink* m2 = new container::MemorySink;
    container::ZipWriter w2{std::unique_ptr<container::Sink>(m2)};
    w2.write("plain", std::string("p"), Method::Store);
    container::EntrySink f = w2.begin("forced", Method::Deflate, true);
    f.write(std::string("forced entry"));
    f.close();
    w2.finish();
    CHECK(count_sig(m2->bytes(), 0x06064b50u) == 0, "per-entry force needs no end record");
    container::ZipReader r2(m2->bytes());
    CHECK(!r2.find("plain")->zip64 && r2.find("forced")->zip64,
          "only the forced entry has ZIP64 fields");
    CHECK(r2.read("forced") == bytes_of("forced entry"), "per-entry forced read");

    // A QVF archive with every ZIP64 structure validates and reads.
    QvfWriter w({"zip64-test", "1.0", "forced"});
    w.auto_metadata = false;
    w.force_zip64 = true;
    w.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}}});
    w.add_atom_properties({-0.68, 0.34});
    const std::vector<uint8_t> q = w.to_bytes();
    CHECK(count_sig(q, 0x06064b50u) == 1, "QvfWriter force_zip64 writes the end record");
    const std::vector<Finding> found = QvfReader::validate_bytes(q);
    CHECK(clean(found), "forced ZIP64 QVF archive validates" + describe(found));
    QvfReader qr = QvfReader::from_bytes(q);
    CHECK(qr.read_binary_member("atom_properties", "mulliken_charge").as_doubles().at(1) == 0.34,
          "forced ZIP64 QVF member reads");
}

static void test_zip64_entry_count() {
    container::MemorySink* mem = new container::MemorySink;
    container::ZipWriter zw{std::unique_ptr<container::Sink>(mem)};
    const size_t n = 70000;
    char name[32];
    for (size_t i = 0; i < n; ++i) {
        std::snprintf(name, sizeof name, "e/%05zu", i);
        zw.write(name, std::string(name), Method::Store);
    }
    zw.finish();
    const std::vector<uint8_t>& z = mem->bytes();
    CHECK(count_sig(z, 0x07064b50u) >= 1, "70 000 entries: ZIP64 locator written");
    const size_t eocd = z.size() - 22;
    CHECK(rd16(z, eocd + 10) == 0xFFFF, "70 000 entries: 16-bit count saturated");
    CHECK(rd32(z, eocd + 16) != 0xFFFFFFFFu, "70 000 entries: offset not saturated");
    container::ZipReader r(z);
    CHECK(r.entries().size() == n, "70 000 entries read back");
    CHECK(r.read("e/69999") == bytes_of("e/69999") && r.read("e/00000") == bytes_of("e/00000"),
          "70 000 entries: first and last read");
    // 65535 entries exactly: the count is saturated, so ZIP64 records too.
    container::MemorySink* m2 = new container::MemorySink;
    container::ZipWriter w2{std::unique_ptr<container::Sink>(m2)};
    for (size_t i = 0; i < 65535; ++i) {
        std::snprintf(name, sizeof name, "%05zu", i);
        w2.write(name, std::string(), Method::Store);
    }
    w2.finish();
    CHECK(container::ZipReader(m2->bytes()).entries().size() == 65535,
          "65 535 entries read back");
}

// A seekable sink whose first byte is at `base` (only the bytes from there
// on are kept): offsets past 4 GiB without writing 4 GiB.
class OffsetSink : public container::Sink {
public:
    explicit OffsetSink(uint64_t base) : base_(base) {}
    void write(const uint8_t* d, size_t n) override { bytes.insert(bytes.end(), d, d + n); }
    uint64_t position() const override { return base_ + bytes.size(); }
    bool seekable() const override { return true; }
    void write_at(uint64_t off, const uint8_t* d, size_t n) override {
        std::memcpy(bytes.data() + (off - base_), d, n);
    }
    std::vector<uint8_t> bytes;

private:
    uint64_t base_;
};

static void test_zip64_offsets() {
    const uint64_t base = 5ull << 30;   // 5 GiB
    OffsetSink* sink = new OffsetSink(base);
    container::ZipWriter zw{std::unique_ptr<container::Sink>(sink)};
    zw.write("far.txt", std::string("beyond 4 GiB"), Method::Store);
    container::EntrySink s = zw.begin("far2.bin", Method::Deflate);
    s.write(input(3));
    s.close();
    zw.finish();
    const std::vector<uint8_t>& z = sink->bytes;
    // Central records: only the offset moved to the ZIP64 extra field.
    size_t cd = 0;
    for (size_t p = 0; p + 4 <= z.size(); ++p)
        if (rd32(z, p) == 0x02014b50u) { cd = p; break; }
    CHECK(cd != 0, "central directory found");
    CHECK(rd32(z, cd + 42) == 0xFFFFFFFFu && rd32(z, cd + 20) != 0xFFFFFFFFu,
          "offset saturated, sizes not");
    const uint16_t namelen = rd16(z, cd + 28);
    const size_t extra = cd + 46 + namelen;
    CHECK(rd16(z, cd + 30) == 12 && rd16(z, extra) == 1 && rd16(z, extra + 2) == 8 &&
              rd64(z, extra + 4) == base,
          "ZIP64 extra holds exactly the local-header offset");
    CHECK(rd16(z, cd + 6) == 45, "version needed 45 with a ZIP64 extra");
    // End records: the 32-bit offset is saturated and the ZIP64 record has it.
    const size_t eocd = z.size() - 22;
    CHECK(rd32(z, eocd + 16) == 0xFFFFFFFFu && rd16(z, eocd + 10) == 2,
          "end record: offset saturated, count kept");
    const size_t loc = eocd - 20;
    CHECK(rd32(z, loc) == 0x07064b50u, "locator before the end record");
    const size_t rec = static_cast<size_t>(rd64(z, loc + 8) - base);
    CHECK(rd32(z, rec) == 0x06064b50u && rd64(z, rec + 32) == 2 &&
              rd64(z, rec + 48) == base + cd,
          "ZIP64 end record: count and central-directory offset");
}

static void test_data_descriptors() {
    std::vector<uint8_t> out;
    container::ZipWriter zw{std::unique_ptr<container::Sink>(
        new container::CallbackSink([&](const uint8_t* p, size_t n) {
            out.insert(out.end(), p, p + n);
        }))};
    const std::vector<uint8_t> six = input(6);
    container::EntrySink a = zw.begin("deflated.bin", Method::Deflate);
    for (size_t p = 0; p < six.size(); p += 5000)
        a.write(six.data() + p, std::min<size_t>(5000, six.size() - p));
    a.close();
    container::EntrySink b = zw.begin("stored.txt", Method::Store);
    b.write(std::string("stored with a descriptor"));
    b.close();
    container::EntrySink c = zw.begin("zip64.bin", Method::Deflate, true);
    c.write(six);
    c.close();
    // Whole entries need no descriptor even on a non-seekable sink.
    zw.write("whole.txt", std::string("whole entry"), Method::Store);
    zw.finish();
    CHECK(count_sig(out, 0x08074b50u) == 3, "three data descriptors");
    container::ZipReader r(out);
    for (const char* n : {"deflated.bin", "stored.txt", "zip64.bin"}) {
        const container::ZipEntry* e = r.find(n);
        CHECK(e && (e->flags & 0x08), std::string("bit 3 set: ") + n);
        const container::LocalHeader lh = r.local_header(*e);
        CHECK(lh.crc32 == 0 && lh.csize == 0 && lh.usize == 0,
              std::string("local header zeros: ") + n);
    }
    CHECK(!(r.find("whole.txt")->flags & 0x08), "whole entry has no descriptor");
    CHECK(r.read("deflated.bin") == six && r.read("zip64.bin") == six,
          "descriptor entries read from the central directory");
    CHECK(r.read("stored.txt") == bytes_of("stored with a descriptor"), "stored descriptor");
    // The ZIP64 descriptor carries 8-byte sizes.
    const container::ZipEntry* z64 = r.find("zip64.bin");
    const size_t desc = static_cast<size_t>(r.local_header(*z64).data_offset + z64->csize);
    CHECK(rd32(out, desc) == 0x08074b50u && rd64(out, desc + 8) == z64->csize &&
              rd64(out, desc + 16) == six.size(),
          "ZIP64 data descriptor layout");

    // A whole QVF archive streamed through a callback validates.
    QvfWriter w({"descriptor-test", "1.0", "streamed"});
    w.auto_metadata = false;
    w.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}}});
    w.add_atom_properties({-0.68, 0.34});
    container::ZipReader plain(w.to_bytes());
    std::vector<uint8_t> streamed;
    container::ZipWriter sw{std::unique_ptr<container::Sink>(
        new container::CallbackSink([&](const uint8_t* p, size_t n) {
            streamed.insert(streamed.end(), p, p + n);
        }))};
    for (const auto& e : plain.entries()) {
        container::EntrySink es = sw.begin(e.name, e.name == "manifest.json" ? Method::Store
                                                                             : Method::Deflate);
        es.write(plain.read(e.name));
        es.close();
    }
    sw.finish();
    const std::vector<Finding> found = QvfReader::validate_bytes(streamed);
    CHECK(clean(found), "data-descriptor QVF archive validates" + describe(found));
}

static void test_reader_streaming_and_caps() {
    // A 3 MiB compressible entry, a highly compressible one, a stored one.
    std::vector<uint8_t> big;
    char line[64];
    for (int i = 0; big.size() < (3u << 20); ++i) {
        const int n = std::snprintf(line, sizeof line, "row %d value %d\n", i, (i * 7) % 1000);
        big.insert(big.end(), line, line + n);
    }
    const std::vector<uint8_t> zeros(4u << 20, 0);   // DEFLATEs about 159:1
    container::MemorySink* mem = new container::MemorySink;
    container::ZipWriter zw{std::unique_ptr<container::Sink>(mem)};
    zw.write("big.bin", big);
    zw.write("zeros.bin", zeros);
    zw.write("small.txt", std::string("small"), Method::Store);
    zw.finish();
    const std::vector<uint8_t> z = mem->bytes();
    container::ZipReader r(z);

    std::vector<uint8_t> got;
    size_t chunks = 0;
    r.read_into("big.bin", [&](const uint8_t* p, size_t n) {
        got.insert(got.end(), p, p + n);
        ++chunks;
    }, container::Limits{});
    CHECK(got == big && chunks >= 2, "read_into streams a large entry in chunks");

    // Caps, checked before anything is delivered.
    container::Limits lim;
    lim.max_member_bytes = 1000;
    bool delivered = false;
    CHECK(throws([&] {
              r.read_into("big.bin", [&](const uint8_t*, size_t) { delivered = true; }, lim);
          }, "member cap"),
          "read_into: member cap");
    CHECK(!delivered, "nothing delivered over the member cap");
    container::Limits total;
    total.max_total_bytes = 1u << 20;
    CHECK(throws([&] { r.read_into("small.txt", [](const uint8_t*, size_t) {}, total); },
                 "total cap"),
          "read_into: total cap");
    CHECK(r.find("big.bin")->method == 8 && r.find("zeros.bin")->method == 8, "deflated");
    container::Limits tight;
    tight.max_ratio = 100;
    CHECK(throws([&] { r.read_into("zeros.bin", [](const uint8_t*, size_t) {}, tight); },
                 "ratio"),
          "read_into: compression-ratio cap (4 MiB of zeros, cap 100)");
    const container::Limits loose{};
    size_t n_zeros = 0;
    r.read_into("zeros.bin", [&](const uint8_t*, size_t n) { n_zeros += n; }, loose);
    CHECK(n_zeros == zeros.size(), "read_into under the reference ratio cap");
    CHECK(throws([&] { r.read("zeros.bin", tight); }, "ratio"),
          "read(name, Limits): ratio cap");
    CHECK(r.read("zeros.bin") == zeros, "read(name) without limits");

    // An understated size is never inflated past, while streaming either.
    std::vector<uint8_t> liar = z;
    const container::ZipEntry* e = r.find("big.bin");
    for (size_t p = 0; p + 46 <= liar.size(); ++p)
        if (rd32(liar, p) == 0x02014b50u && rd16(liar, p + 28) == 7 &&
            std::memcmp(&liar[p + 46], "big.bin", 7) == 0) {
            const uint32_t small = 1000;
            for (int b = 0; b < 4; ++b) liar[p + 24 + b] = static_cast<uint8_t>(small >> (8 * b));
        }
    container::ZipReader lr(liar);
    CHECK(lr.find("big.bin")->usize == 1000 && e->usize == big.size(), "size patched");
    CHECK(throws([&] { lr.read_into("big.bin", [](const uint8_t*, size_t) {}, loose); },
                 "exceeds the declared size"),
          "streaming inflate stops at the declared size");

    // CRC-32 mismatch in a stored entry.
    std::vector<uint8_t> flipped = z;
    const container::LocalHeader lh = r.local_header(*r.find("small.txt"));
    flipped[static_cast<size_t>(lh.data_offset)] ^= 0x20;
    CHECK(throws([&] { container::ZipReader(flipped).read("small.txt"); }, "CRC-32"),
          "CRC-32 mismatch detected");

    // A local header naming another entry.
    std::vector<uint8_t> renamed = z;
    const size_t name_at = static_cast<size_t>(r.find("small.txt")->local_offset + 30);
    renamed[name_at] = 'S';
    CHECK(throws([&] { container::ZipReader(renamed).read("small.txt"); }, "local header name"),
          "local-header name mismatch detected");

    // File-backed reading.
    write_file("wp40_reader.zip", z);
    {   // closed before the file is removed (Windows cannot remove an open file)
        container::ZipReader fr = container::ZipReader::open("wp40_reader.zip");
        std::vector<uint8_t> fbig;
        fr.read_into("big.bin",
                     [&](const uint8_t* p, size_t n) { fbig.insert(fbig.end(), p, p + n); },
                     container::Limits{});
        CHECK(fbig == big && fr.read("small.txt") == bytes_of("small"), "file-backed reads");
    }
    fs::remove("wp40_reader.zip");
}

static void test_writer_rules() {
    container::MemorySink* mem = new container::MemorySink;
    container::ZipWriter zw{std::unique_ptr<container::Sink>(mem)};
    zw.write("a", std::string("1"));
    CHECK(throws([&] { zw.write("a", std::string("2")); }, "duplicate"), "duplicate refused");
    container::EntrySink s = zw.begin("b");
    CHECK(throws([&] { zw.write("c", std::string("3")); }, "streamed"),
          "no second entry while one streams");
    {
        container::EntrySink abandoned = std::move(s);
        abandoned.write(std::string("never listed"));
    }   // destroyed open: abandoned
    zw.write("c", std::string("3"));
    CHECK(zw.contains("c") && !zw.contains("b"), "abandoned entry not listed");
    CHECK(zw.drop("a") && !zw.contains("a"), "drop removes an entry");
    zw.write("\xc3\xa9t\xc3\xa9.txt", std::string("utf-8"));
    zw.finish();
    CHECK(throws([&] { zw.write("d", std::string("4")); }, "finished"), "no write after finish");
    container::ZipReader r(mem->bytes());
    CHECK(r.names() == std::vector<std::string>({"c", "\xc3\xa9t\xc3\xa9.txt"}),
          "only the kept entries are listed");
    CHECK(r.find("\xc3\xa9t\xc3\xa9.txt")->flags & 0x0800, "UTF-8 name flag (bit 11)");
    CHECK(!(r.find("c")->flags & 0x0800), "ASCII name has no UTF-8 flag");
}

// --- QvfWriter: streaming, update, repack ------------------------------------------

static QvfWriter water(const std::string& calc = "water") {
    QvfWriter w({"container-test", "1.0", calc});
    w.auto_metadata = false;
    w.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}},
                     {"H", 1, {{0, -0.757, -0.469}}}});
    w.add_atom_properties({-0.70, 0.35, 0.35});
    Json payload = Json::object();
    payload.set("note", "streamed grid follows");
    Json members = Json::object();
    members.set("meta", payload);
    w.add_vendor_section("x_wp40.grid", members, false, "", "grid");
    return w;
}

static std::vector<double> grid_values(size_t n) {
    std::vector<double> v(n);
    for (size_t i = 0; i < n; ++i) v[i] = std::sin(static_cast<double>(i) * 0.37) * 1e-3;
    return v;
}

static void push_doubles(MemberSink& s, const std::vector<double>& v, size_t slab) {
    std::vector<uint8_t> buf;
    for (size_t p = 0; p < v.size(); p += slab) {
        const size_t n = std::min(slab, v.size() - p);
        buf.assign(n * 8, 0);
        for (size_t i = 0; i < n; ++i) {
            uint64_t b;
            std::memcpy(&b, &v[p + i], 8);
            for (int j = 0; j < 8; ++j) buf[i * 8 + j] = static_cast<uint8_t>(b >> (8 * j));
        }
        s.write(buf);
    }
}

static void test_writer_to_bytes_equals_write() {
    QvfWriter w = water();
    const std::vector<uint8_t> mem = w.to_bytes();
    const std::string path = w.write("wp40_plain");
    CHECK(path == "wp40_plain.qvf", "write() appends .qvf");
    CHECK(read_file(path) == mem, "write() streams the same bytes to_bytes() returns");
    container::ZipReader r(mem);
    CHECK(r.entries().back().name == "manifest.json" && r.entries().back().method == 0,
          "manifest.json last and STOREd");
    fs::remove(path);
}

static void test_streaming_vs_buffered() {
    const std::vector<double> v = grid_values(50000);
    const std::vector<int64_t> shape = {50, 1000};
    // Buffered writer: begin_binary keeps the bytes until to_bytes().
    QvfWriter b = water();
    MemberSink bs = b.begin_binary("grid", "data", DType::Float64, shape);
    CHECK(throws([&] { b.add_atom_properties({1.0}, {}, {}, "other"); }, "member sink"),
          "no add_* while a member sink is open");
    push_doubles(bs, v, 777);
    bs.close();
    const std::vector<uint8_t> buffered = b.to_bytes();
    // The same through a tensor, the old way.
    QvfReader br = QvfReader::from_bytes(buffered);
    const std::vector<uint8_t> raw = br.read_member_bytes("grid", "data");
    CHECK(raw.size() == v.size() * 8, "buffered streamed member size");
    // Streaming writer: members go to the file as they are added.
    QvfWriter s({"container-test", "1.0", "water"}, "wp40_stream.qvf");
    s.auto_metadata = false;
    s.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}},
                     {"H", 1, {{0, -0.757, -0.469}}}});
    s.add_atom_properties({-0.70, 0.35, 0.35});
    Json payload = Json::object();
    payload.set("note", "streamed grid follows");
    Json members = Json::object();
    members.set("meta", payload);
    s.add_vendor_section("x_wp40.grid", members, false, "", "grid");
    CHECK(throws([&] { s.to_bytes(); }, "write()"), "to_bytes() refused while streaming");
    MemberSink ss = s.begin_binary("grid", "data", DType::Float64, shape);
    push_doubles(ss, v, 4096);
    CHECK(ss.bytes_written() == ss.bytes_expected(), "all bytes pushed");
    ss.close();
    CHECK(s.write() == "wp40_stream.qvf", "streaming write() returns the bound path");
    const std::vector<Finding> found = QvfReader::validate_file("wp40_stream.qvf");
    CHECK(clean(found), "streamed archive validates" + describe(found));
    {   // the reader keeps the file open; close it before removing the file
        QvfReader sr = QvfReader::open("wp40_stream.qvf");
        CHECK(sr.read_member_bytes("grid", "data") == raw,
              "streaming and buffered members byte-identical");
        for (const char* sec : {"structure", "atom_properties", "grid"}) {
            const Json& a = br.section(sec).at("members");
            const Json& c = sr.section(sec).at("members");
            for (const auto& kv : a.items())
                CHECK(c.at(kv.first).at("sha256").as_string() ==
                          kv.second.at("sha256").as_string(),
                      std::string("same sha256: ") + sec + "/" + kv.first);
        }
        CHECK(sr.manifest().dump(2) == br.manifest().dump(2), "same manifest either way");
        CHECK(sr.read_binary_member("grid", "data").as_doubles() == v, "streamed values");
    }

    // Short and long writes are refused; an abandoned member is not listed.
    QvfWriter x = water();
    MemberSink shortw = x.begin_binary("grid", "data", DType::Float64, {4});
    shortw.write(std::vector<uint8_t>(16, 0));
    CHECK(throws([&] { shortw.close(); }, "needs 32"), "short member refused");
    MemberSink longw = x.begin_binary("grid", "data", DType::Int32, {2});
    CHECK(throws([&] { longw.write(std::vector<uint8_t>(12, 0)); }, "more bytes"),
          "long member refused");
    { MemberSink dropped = std::move(longw); }
    CHECK(!x.build_manifest().at("sections").at(2).at("members").has("data"),
          "abandoned members are not listed");
    CHECK(throws([&] { x.begin_binary("nope", "data", DType::Int8, {1}); }, "no section"),
          "begin_binary needs an existing section");
    CHECK(throws([&] { x.begin_binary("grid", "meta", DType::Int8, {1}); }, "already has"),
          "begin_binary refuses an existing role");
    fs::remove("wp40_stream.qvf");
}

static size_t count_code(const std::vector<Finding>& found, const std::string& code) {
    size_t n = 0;
    for (const auto& f : found) n += f.code == code;
    return n;
}

static void test_update_validate_repack() {
    const std::string path = "wp40_update.qvf";
    {
        QvfWriter w = water();
        w.set_atom_properties_meta(Json::parse(
            R"({"custom.flux": {"units": "e", "description": "Net flux.", "kind": "scalar"}})"));
        w.add_atom_property("custom.flux", {0.1, -0.05, -0.05});
        w.write(path);
    }
    const std::vector<uint8_t> before = read_file(path);

    // An update never rewrites an existing entry (spec § 2.3): a new custom
    // role rewrites the section's meta member, which the archive holds as
    // atom_properties/meta.json, so the new one gets a fresh path and the old
    // one stays as a superseded entry.
    {
        QvfWriter fresh = QvfWriter::open_for_update(path);
        fresh.set_atom_properties_meta(Json::parse(
            R"({"custom.flow": {"units": "e", "description": "Flow.", "kind": "scalar"}})"));
        fresh.add_atom_property("custom.flow", {0.2, -0.1, -0.1});
        const Json manifest = fresh.build_manifest();
        std::string meta_path;
        for (const Json& s : manifest.at("sections").elements())
            if (s.at("id").as_string() == "atom_properties")
                meta_path = s.at("members").at("meta").at("path").as_string();
        CHECK(meta_path == "atom_properties/meta.2.json",
              "a member whose path the archive holds gets a fresh one: " + meta_path);
        const std::vector<Finding> found = QvfReader::validate_bytes(fresh.to_bytes());
        CHECK(!has_errors(found) && count_code(found, "W_ORPHAN_MEMBER") == 2,
              "the old meta.json and the superseded manifest are orphans" + describe(found));
        CHECK(read_file(path) == before, "to_bytes() of an update leaves the archive untouched");
    }

    // 1. Append a section and a vendor section (spec § 2.3).
    QvfWriter u = QvfWriter::open_for_update(path);
    Spectrum ir;
    ir.frequencies = {1595.0, 3657.0};
    ir.intensities = {67.0, 3.0};
    u.add_spectrum("spectra.ir", ir, "", "ir");
    Json note = Json::object();
    note.set("step", 1);
    Json note_members = Json::object();
    note_members.set("note", note);
    u.add_vendor_section("x_wp40.note", note_members);
    CHECK(throws([&] { u.add_spectrum("spectra.ir", ir, "", "ir"); }, "duplicate"),
          "existing section ids stay taken");
    const std::vector<uint8_t> preview = u.to_bytes();
    CHECK(read_file(path) == before, "to_bytes() leaves the archive untouched");
    u.write();
    CHECK(read_file(path) == preview, "write() appends exactly what to_bytes() previewed");
    std::vector<Finding> found = QvfReader::validate_file(path);
    CHECK(!has_errors(found) && count_code(found, "W_ORPHAN_MEMBER") == 1,
          "updated archive validates; the superseded manifest is an orphan" + describe(found));
    // Readers are scoped: they keep the file open, and Windows cannot replace
    // or remove an open file.
    {
        container::ZipReader zr = container::ZipReader::open(path);
        CHECK(zr.entries().back().name == "manifest.json" && zr.entries().back().method == 0,
              "updated manifest.json last and STOREd");
        CHECK(zr.duplicates() == std::vector<std::string>{"manifest.json"},
              "the superseded manifest stays listed; no other name repeats");
        CHECK(zr.find("manifest.json") == &zr.entries().back(),
              "a lookup resolves manifest.json to the last record");
    }
    {
        const std::vector<uint8_t> after = read_file(path);
        const size_t old_cd = static_cast<size_t>(
            container::ZipReader(before).central_directory_offset());
        CHECK(after.size() > old_cd &&
                  std::equal(before.begin(), before.begin() + static_cast<std::ptrdiff_t>(old_cd),
                             after.begin()),
              "the old entries' bytes are untouched");
    }
    {
        QvfReader r = QvfReader::open(path);
        CHECK(r.has_section("ir") && r.has_section("structure") && r.has_section("x_wp40_note"),
              "old and new sections");
        CHECK(r.read_binary_member("atom_properties", "custom.flux").as_doubles().at(0) == 0.1 &&
                  r.read_json_member("ir", "spectrum").at("frequencies").at(0).as_double() ==
                      1595.0,
              "old and new members read, sha256-verified");
        CHECK(r.manifest().at("source").at("calculation").as_string() == "water",
              "source kept");
        const Json& meta = r.manifest().at("metadata");
        CHECK(meta.has("modified_utc") && meta.at("modified_utc").as_string().size() == 20,
              "the update stamps metadata.modified_utc");
    }

    // 2. Stream a member into the archive in place.
    const std::vector<double> v = grid_values(20000);
    {
        QvfWriter u2 = QvfWriter::open_for_update(path);
        MemberSink s = u2.begin_binary("grid", "data", DType::Float64, {20000});
        push_doubles(s, v, 3000);
        s.close();
        u2.write(path);   // its own path is accepted
    }
    found = QvfReader::validate_file(path);
    CHECK(!has_errors(found) && count_code(found, "W_ORPHAN_MEMBER") == 2,
          "a member streamed in place: valid, two superseded manifests" + describe(found));
    CHECK(QvfReader::open(path).read_binary_member("grid", "data").as_doubles() == v,
          "member streamed in place reads back");

    // 3. An orphan (an entry appended without a manifest reference), then repack.
    {
        container::ZipWriter raw = container::ZipWriter::append(path);
        raw.write("stale/old.bin", std::vector<uint8_t>(5000, 7));
        raw.finish();
    }
    found = QvfReader::validate_file(path);
    CHECK(count_code(found, "W_ORPHAN_MEMBER") == 3 && !has_errors(found),
          "orphans (the stale entry, two superseded manifests) are warnings" + describe(found));
    const uint64_t with_dead = fs::file_size(path);
    repack(path, "wp40_repacked.qvf");
    found = QvfReader::validate_file("wp40_repacked.qvf");
    CHECK(clean(found), "repacked archive has no orphan" + describe(found));
    CHECK(fs::file_size("wp40_repacked.qvf") < with_dead, "repack drops dead space");
    {
        QvfReader rp = QvfReader::open("wp40_repacked.qvf");
        CHECK(rp.read_binary_member("grid", "data").as_doubles() == v &&
                  rp.read_json_member("ir", "spectrum").at("frequencies").at(1).as_double() ==
                      3657.0,
              "repacked members intact");
        container::ZipReader pz = container::ZipReader::open("wp40_repacked.qvf");
        CHECK(pz.entries().back().name == "manifest.json" && pz.entries().back().method == 0 &&
                  !pz.has("stale/old.bin"),
              "repack: manifest last, orphan gone");
    }
    repack(path, path);   // in place
    CHECK(clean(QvfReader::validate_file(path)), "in-place repack");
    CHECK(read_file(path) == read_file("wp40_repacked.qvf"), "repack is deterministic");

    // 4. Refusals.
    CHECK(throws([&] { QvfWriter::open_for_update("wp40_missing.qvf"); }, "cannot open"),
          "open_for_update on a missing file");
    QvfWriter bound = QvfWriter::open_for_update(path);
    CHECK(throws([&] { bound.write("elsewhere.qvf"); }, "bound"),
          "a bound writer refuses another path");
    QvfWriter plain = water();
    CHECK(throws([&] { plain.write(); }, "path"), "write() without a path needs a bound writer");
    fs::remove(path);
    fs::remove("wp40_repacked.qvf");
}

// One empty STORE entry, hand-built: local header, central directory and end
// record, with the local and central names given separately.
static std::vector<uint8_t> one_entry_zip(const std::string& local_name,
                                          const std::string& cd_name) {
    std::vector<uint8_t> z;
    auto p16 = [&z](uint32_t v) { z.push_back(uint8_t(v)); z.push_back(uint8_t(v >> 8)); };
    auto p32 = [&](uint32_t v) { p16(v & 0xFFFF); p16(v >> 16); };
    p32(0x04034b50u); p16(20); p16(0); p16(0); p16(0); p16(0);   // local header
    p32(0); p32(0); p32(0); p16(uint32_t(local_name.size())); p16(0);
    z.insert(z.end(), local_name.begin(), local_name.end());
    const uint32_t cd = uint32_t(z.size());
    p32(0x02014b50u); p16(20); p16(20); p16(0); p16(0); p16(0); p16(0);   // central
    p32(0); p32(0); p32(0); p16(uint32_t(cd_name.size())); p16(0); p16(0);
    p16(0); p16(0); p32(0); p32(0);
    z.insert(z.end(), cd_name.begin(), cd_name.end());
    const uint32_t cd_size = uint32_t(z.size()) - cd;
    p32(0x06054b50u); p16(0); p16(0); p16(1); p16(1); p32(cd_size); p32(cd); p16(0);
    return z;
}

// Fuzzing regression (WP-50): an entry with an empty name and a local header
// without name or extra field once passed a null pointer to memcmp (UBSan).
static void test_empty_entry_name() {
    const container::ZipReader r(one_entry_zip("", ""));
    CHECK(r.entries().size() == 1, "empty-name entry is listed");
    bool ok = false;
    try {
        ok = r.local_header(r.entries()[0]).data_offset == 30;
    } catch (const std::exception& e) {
        std::cerr << "    " << e.what() << "\n";
    }
    CHECK(ok, "empty-name local header reads (data at offset 30)");
    const container::ZipReader r2(one_entry_zip("x", ""));
    CHECK(throws([&] { (void)r2.local_header(r2.entries()[0]); }, "differs"),
          "local name 'x' against an empty central name is refused");
}

static void test_zip_append_container() {
    container::MemorySink* mem = new container::MemorySink;
    container::ZipWriter zw{std::unique_ptr<container::Sink>(mem)};
    zw.write("keep.txt", std::string("keep"));
    zw.write("replace.txt", std::string("old"));
    zw.finish();
    container::ZipWriter ap = container::ZipWriter::append(mem->take());
    CHECK(ap.contains("keep.txt") && ap.contains("replace.txt"), "append keeps entries");
    CHECK(throws([&] { ap.write("keep.txt", std::string("x")); }, "duplicate"),
          "append refuses a duplicate name");
    CHECK(ap.drop("replace.txt"), "drop a kept entry");
    ap.write("replace.txt", std::string("new"));
    ap.finish();
    container::ZipReader r(static_cast<container::MemorySink&>(ap.sink()).bytes());
    // An explicitly allowed repeat: both records listed, the last one wins,
    // and a further append keeps both.
    container::ZipWriter rep = container::ZipWriter::append(
        static_cast<container::MemorySink&>(ap.sink()).bytes());
    rep.allow_repeated_name("keep.txt");
    rep.write("keep.txt", std::string("kept twice"));
    rep.finish();
    container::ZipWriter rep2 = container::ZipWriter::append(
        static_cast<container::MemorySink&>(rep.sink()).bytes());
    rep2.finish();
    container::ZipReader rr(static_cast<container::MemorySink&>(rep2.sink()).bytes());
    CHECK(rr.duplicates() == std::vector<std::string>{"keep.txt"} &&
              rr.read("keep.txt") == bytes_of("kept twice") && rr.entries().size() == 3,
          "allow_repeated_name: both records kept, the last one read");
    CHECK(r.duplicates().empty() && r.read("replace.txt") == bytes_of("new") &&
              r.read("keep.txt") == bytes_of("keep"),
          "appended archive: kept entry and replacement read");
    // Appending to a forced-ZIP64 archive keeps its ZIP64 central records.
    container::MemorySink* m64 = new container::MemorySink;
    container::ZipWriter w64{std::unique_ptr<container::Sink>(m64), true};
    w64.write("a", std::string("aaa"));
    w64.finish();
    container::ZipWriter ap64 = container::ZipWriter::append(m64->take());
    ap64.write("b", std::string("bbb"));
    ap64.finish();
    container::ZipReader r64(static_cast<container::MemorySink&>(ap64.sink()).bytes());
    CHECK(r64.find("a")->zip64 && !r64.find("b")->zip64 && r64.read("a") == bytes_of("aaa"),
          "kept ZIP64 records survive an append");
}

// The corpus' ZIP64 and data-descriptor archives read, update and repack.
static void test_corpus_containers(const std::string& corpus) {
    for (const char* name : {"zip64.qvf", "data_descriptor.qvf", "updated_in_place.qvf"}) {
        const std::string src = corpus + "/" + name;
        if (!fs::exists(src)) {
            CHECK(std::string(name) == "updated_in_place.qvf",
                  std::string("corpus archive missing: ") + name);
            continue;
        }
        const std::vector<Finding> original = QvfReader::validate_file(src);
        CHECK(!has_errors(original), std::string("corpus reads: ") + name + describe(original));
        const size_t orphans = count_code(original, "W_ORPHAN_MEMBER");
        size_t n_sections = 0;
        {
            QvfReader r = QvfReader::open(src);
            n_sections = r.sections().size();
        }
        const std::string copy = std::string("wp40_corpus_") + name;
        fs::copy_file(src, copy, fs::copy_options::overwrite_existing);
        QvfWriter u = QvfWriter::open_for_update(copy);
        Json payload = Json::object();
        payload.set("container", name);
        Json members = Json::object();
        members.set("note", payload);
        u.add_vendor_section("x_wp40.note", members);
        u.write();
        std::vector<Finding> found = QvfReader::validate_file(copy);
        CHECK(!has_errors(found) && count_code(found, "W_ORPHAN_MEMBER") == orphans + 1,
              std::string("updated corpus archive validates: ") + name + describe(found));
        repack(copy, copy);
        found = QvfReader::validate_file(copy);
        CHECK(clean(found), std::string("repacked corpus archive validates: ") + name + describe(found));
        {
            QvfReader r = QvfReader::open(copy);
            CHECK(r.has_section("x_wp40_note") && r.sections().size() == n_sections + 1,
                  std::string("old and new sections after update: ") + name);
        }
        fs::remove(copy);
    }
}

// Dynamic-Huffman reference-writer data exceeds 1000:1. Exercise both the
// buffered and streamed public container APIs under their default limits.
static void test_zero_grid(const std::string& path) {
    auto zip = container::ZipReader::open(path);
    const uint64_t size = 64ull * 64 * 64 * 8;
    const container::ZipEntry* entry = nullptr;
    for (const auto& e : zip.entries()) if (e.usize == size) entry = &e;
    CHECK(entry && entry->csize * 1000 < size, "reference grid exceeds ratio 1000");
    if (!entry) return;
    const auto bytes = zip.read(entry->name, container::Limits{});
    CHECK(bytes.size() == size && std::all_of(bytes.begin(), bytes.end(),
          [](uint8_t b) { return b == 0; }), "default buffered read accepts zero grid");
    uint64_t received = 0;
    zip.read_into(entry->name, [&](const uint8_t* data, size_t n) {
        received += n;
        CHECK(std::all_of(data, data + n, [](uint8_t b) { return b == 0; }),
              "streamed grid contains only zeros");
    });
    CHECK(received == size, "default streaming read accepts complete zero grid");
    container::Limits explicit_ratio;
    explicit_ratio.max_ratio = 1000;
    bool delivered = false;
    CHECK(throws([&] { zip.read_into(entry->name, [&](const uint8_t*, size_t) {
        delivered = true;
    }, explicit_ratio); }, "ratio"), "explicit ratio policy still rejects");
    CHECK(!delivered, "explicit cap rejects before delivering data");
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--zero-grid") {
        test_zero_grid(argv[2]);
        return g_failures ? 1 : 0;
    }
    test_deflate_golden_and_chunking();
    test_inflate_streaming();
    test_golden_archive();
    test_zip64_forced();
    test_zip64_entry_count();
    test_zip64_offsets();
    test_data_descriptors();
    test_reader_streaming_and_caps();
    test_writer_rules();
    test_writer_to_bytes_equals_write();
    test_streaming_vs_buffered();
    test_update_validate_repack();
    test_zip_append_container();
    test_empty_entry_name();
    if (argc > 1 && fs::is_directory(argv[1])) {
        test_corpus_containers(argv[1]);
        test_zero_grid((fs::path(argv[1]) / "volume_spin_zero.qvf").string());
    }
    if (g_failures) {
        std::cerr << g_failures << " container test(s) failed\n";
        return 1;
    }
    std::cout << "container tests passed\n";
    return 0;
}
