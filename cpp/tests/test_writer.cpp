// test_writer.cpp — unit tests for the QVF C++ writer.
//
// Known-answer tests for the vendored-free SHA-256 / CRC-32 / primitive_norm
// primitives, plus an archive smoke test. Full schema + integrity validation of
// the emitted archive is done cross-language by python/qvf_reader.py (run in the
// build verification step), which is the authoritative check.
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "qvf/qvf.hpp"
#include "qvf/qvf_c.h"
#include "qvf/qvf_reader.hpp"
#include "qvf_deflate.hpp"
#include "qvf_inflate.hpp"
#include "qvf_sha256.hpp"
#include "qvf_zip.hpp"

#include <cstdio>
#include <fstream>

static int g_failures = 0;

#define CHECK(cond, msg)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::cerr << "FAIL: " << (msg) << " (" << #cond << ")\n";      \
            ++g_failures;                                                  \
        }                                                                  \
    } while (0)

static bool contains(const std::vector<uint8_t>& hay, const std::string& needle) {
    std::string s(hay.begin(), hay.end());
    return s.find(needle) != std::string::npos;
}

int main() {
    using namespace qvf;

    // --- SHA-256 known-answer (FIPS 180-4 "abc") ---
    CHECK(detail::Sha256::hex_of(std::string("abc")) ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
          "sha256(abc)");
    CHECK(detail::Sha256::hex_of(std::string("")) ==
          "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
          "sha256(empty)");
    // A > 64-byte message to exercise multi-block processing.
    CHECK(detail::Sha256::hex_of(std::string(
              "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")) ==
          "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
          "sha256(56-byte)");

    // --- CRC-32 known-answer ("123456789" -> 0xCBF43926) ---
    {
        std::string m = "123456789";
        uint32_t crc = detail::crc32_of(reinterpret_cast<const uint8_t*>(m.data()),
                                        m.size());
        CHECK(crc == 0xCBF43926u, "crc32(123456789)");
    }

    // --- primitive_norm ---
    CHECK(std::abs(primitive_norm(1.0, 0) - 0.7127054703549902) < 1e-12,
          "primitive_norm(1,0)");
    // The constexpr pi inside primitive_norm must agree with libm's pi.
    CHECK(std::abs(primitive_norm(1.0, 0) -
                   std::pow(2.0 / std::acos(-1.0), 0.75)) < 1e-15,
          "primitive_norm(1,0) matches libm pi");
    CHECK(std::abs(primitive_norm(1.0, 1) - 1.4254109407099804) < 1e-12,
          "primitive_norm(1,1)");

    // --- Json integer constructors: size_t and uint64_t must be unambiguous
    //     on every ABI (size_t is `unsigned long long` on 64-bit Windows) ---
    {
        const size_t sz = 3;
        const uint64_t u64 = 4;
        const unsigned long long ull = 5;
        CHECK(Json(sz).as_int() == 3, "Json(size_t)");
        CHECK(Json(u64).as_int() == 4, "Json(uint64_t)");
        CHECK(Json(ull).as_int() == 5, "Json(unsigned long long)");
    }

    // --- tensor byte layout: little-endian, correct size ---
    {
        Tensor t = tensor_f64({1.0, 2.0}, {2});
        CHECK(t.bytes.size() == 16, "tensor_f64 byte size");
        // 1.0 == 0x3FF0000000000000; LE first byte is 0x00, 8th is 0x3F.
        CHECK(t.bytes[7] == 0x3F, "tensor_f64 little-endian");
    }

    // --- DEFLATE compresses repetitive data (round-trip verified externally
    //     by unzip / Python zipfile / the sha256 check in validate_qvf) ---
    {
        std::vector<uint8_t> repetitive(4096, 0);  // all zeros — highly compressible
        auto d0 = detail::deflate(repetitive);
        CHECK(d0.size() < repetitive.size() / 4, "deflate shrinks zeros >4x");
        std::vector<uint8_t> patterned;
        for (int i = 0; i < 1000; ++i) {
            const char* s = "the quick brown fox ";
            patterned.insert(patterned.end(), s, s + 20);
        }
        auto d1 = detail::deflate(patterned);
        CHECK(d1.size() < patterned.size() / 2, "deflate shrinks repeated text >2x");
        // Empty input must not crash.
        auto d2 = detail::deflate(std::vector<uint8_t>{});
        CHECK(!d2.empty(), "deflate of empty emits a valid block");
    }

    // --- archive smoke test ---
    QvfWriter w({"test-code", "9.9", "smoke"});
    w.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}}});
    Spectrum ir;
    ir.frequencies = {1595.0};
    ir.intensities = {67.0};
    w.add_spectrum("spectra.ir", ir);
    w.add_atom_properties({-0.68, 0.34});
    WavefunctionGTO wf;
    wf.shells = {{0, 0, {1.0}, {1.0}, true}};
    wf.mo_coefficients = {1, 1, {1.0}};
    wf.energies = {-0.5};
    wf.occupations = {2};
    w.add_wavefunction_gto(wf);
    Json epr = Json::object();
    Json gt = Json::array();
    gt.push_back(Json(2.0023));
    epr.set("g_iso", Json(2.0023)).set("g_tensor_diag", gt);
    Json ext = Json::object();
    Json orca_ext = Json::object();
    orca_ext.set("version", "1.0").set("critical", true);
    ext.set("x_orca", orca_ext);
    w.set_extensions(ext);
    Json epr_members = Json::object();
    epr_members.set("epr", epr);
    w.add_vendor_section("x_orca.epr", epr_members, /*critical=*/true);

    std::vector<uint8_t> bytes = w.to_bytes();
    CHECK(bytes.size() > 200, "archive non-trivial size");
    CHECK(bytes[0] == 'P' && bytes[1] == 'K' && bytes[2] == 3 && bytes[3] == 4,
          "zip local-file magic");
    CHECK(contains(bytes, "manifest.json"), "contains manifest.json entry");
    CHECK(contains(bytes, "wavefunction.gto"), "contains wavefunction kind");
    CHECK(contains(bytes, "x_orca.epr"), "contains vendor section");

    // --- duplicate id is rejected ---
    {
        bool threw = false;
        try {
            QvfWriter w2({"t", "1", ""});
            w2.add_structure({{"H", 1, {{0, 0, 0}}}}, {{false, false, false}}, {}, {},
                             "dup");
            w2.add_bonds({{0, 0, 1.0}}, "dup");
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(threw, "duplicate id rejected");
    }

    // --- volumes that name their wavefunction and orbital (spec § 5.2) ---
    {
        auto throws = [](auto&& fn) {
            try {
                fn();
            } catch (const std::exception&) {
                return true;
            }
            return false;
        };
        QvfWriter vw({"t", "1", "wavefunction volumes"});
        vw.auto_metadata = false;
        vw.add_structure({{"H", 1, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.74}}}});
        WavefunctionGTO rhf;
        rhf.shells = {{0, 0, {1.0}, {1.0}, true}, {1, 0, {1.0}, {1.0}, true}};
        rhf.mo_coefficients = {2, 2, {1, 0, 0, 1}};
        rhf.energies = {-0.5, 0.5};
        rhf.occupations = {2, 0};
        vw.add_wavefunction_gto(rhf);
        WavefunctionGTO uhf;
        uhf.shells = rhf.shells;
        uhf.unrestricted = true;
        uhf.mo_coefficients_alpha = {2, 2, {1, 0, 0, 1}};
        uhf.mo_coefficients_beta = {1, 2, {1, 0}};
        uhf.energies_alpha = {-0.5, 0.5};
        uhf.occupations_alpha = {1, 0};
        uhf.energies_beta = {-0.4};
        uhf.occupations_beta = {1};
        uhf.section_id = "wf_u";
        vw.add_wavefunction_gto(uhf);
        Grid g;
        const Tensor field = tensor_f32(std::vector<double>(8, 0.125), {2, 2, 2});
        vw.add_volume("volume.density", g, field, "rho", "", "rho");
        vw.add_volume("volume.orbital", g, field, "MO 1", "", "mo_1");
        vw.add_volume("volume.orbital", g, field, "MO 0 beta", "", "mo_b0");
        vw.add_volume("volume.elf", g, field, "ELF", "", "elf");
        CHECK(vw.link_volume_wavefunction("rho", "wf") == "rho" &&
                  vw.link_volume_wavefunction("mo_1", "wf", 1) == "mo_1" &&
                  vw.link_volume_wavefunction("mo_b0", "wf_u", 0, "beta") == "mo_b0",
              "link_volume_wavefunction returns the volume id");
        CHECK(throws([&] { vw.link_volume_wavefunction("elf", "wf"); }),
              "link_volume_wavefunction refuses a volume.elf");
        CHECK(throws([&] { vw.link_volume_wavefunction("rho", "wf", 0); }),
              "link_volume_wavefunction refuses mo_index on a density");
        CHECK(throws([&] { vw.link_volume_wavefunction("mo_1", "structure", 0); }),
              "link_volume_wavefunction refuses a reference to a structure");
        CHECK(throws([&] { vw.link_volume_wavefunction("mo_1", "wf", 2); }),
              "link_volume_wavefunction refuses a row outside mo_coefficients");
        CHECK(throws([&] { vw.link_volume_wavefunction("mo_1", "wf", 0, "alpha"); }),
              "link_volume_wavefunction refuses spin with mo_coefficients");
        CHECK(throws([&] { vw.link_volume_wavefunction("mo_b0", "wf_u", 0); }),
              "link_volume_wavefunction needs spin for an unrestricted wavefunction");
        CHECK(throws([&] { vw.link_volume_wavefunction("mo_b0", "wf_u", 1, "beta"); }),
              "link_volume_wavefunction refuses a row outside mo_coefficients_beta");
        QvfReader vr = QvfReader::from_bytes(vw.to_bytes());
        const Json& secs = vr.manifest().at("sections");
        const Json& mo_b0 = secs.at(5);
        CHECK(secs.at(3).at("wavefunction_ref").as_string() == "wf" &&
                  !secs.at(3).has("mo_index"),
              "a linked density carries wavefunction_ref only");
        CHECK(mo_b0.at("wavefunction_ref").as_string() == "wf_u" &&
                  mo_b0.at("mo_index").as_int() == 0 &&
                  mo_b0.at("spin").as_string() == "beta",
              "a linked unrestricted orbital carries wavefunction_ref, mo_index and spin");
        CHECK(!secs.at(4).has("spin"), "a restricted orbital carries no spin");
        CHECK(!has_errors(vr.validate_findings()), "linked volumes validate");
    }

    // --- reserved ids and a bonds structure_ref (spec § 3.4, § 5.1) ---
    {
        auto throws = [](auto&& fn) {
            try {
                fn();
            } catch (const std::exception&) {
                return true;
            }
            return false;
        };
        QvfWriter bw({"t", "1", "bonds"});
        for (const char* id : {"auto_open", "bookmarks", "camera", "crossfade"})
            CHECK(throws([&] { bw.add_bonds({{0, 1, 1.0}}, id); }),
                  std::string("reserved section id refused: ") + id);
        bw.add_structure({{"O", 8, {{0, 0, 0.1173}}},
                          {"H", 1, {{0, 0.7572, -0.4692}}},
                          {"H", 1, {{0, -0.7572, -0.4692}}}},
                         {{false, false, false}}, {}, {}, "structure_water");
        bw.add_structure({{"H", 1, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.74}}}},
                         {{false, false, false}}, {}, {}, "structure_h2");
        bw.add_job_spec(Json::parse(R"({"job_type": "molecular"})"), "job");
        CHECK(throws([&] { bw.add_bonds({{0, 1, 1.0}}, "bonds_h2", "missing"); }) &&
                  throws([&] { bw.add_bonds({{0, 1, 1.0}}, "bonds_h2", "job"); }),
              "add_bonds refuses a structure_ref that names no structure");
        CHECK(bw.add_bonds({{0, 1, 1.0}}, "bonds_h2", "structure_h2") == "bonds_h2",
              "add_bonds with structure_ref (the id was not consumed by the refusals)");
        bw.add_bonds({{0, 1, 1.0}, {0, 2, 1.0}}, "bonds_water");
        QvfReader br = QvfReader::from_bytes(bw.to_bytes());
        const Json& secs = br.manifest().at("sections");
        CHECK(secs.size() == 5, "bonds archive has five sections");
        const Json& bonds_h2 = secs.at(3);
        CHECK(bonds_h2.at("id").as_string() == "bonds_h2" &&
                  bonds_h2.at("structure_ref").as_string() == "structure_h2",
              "add_bonds writes the section-level structure_ref");
        CHECK(!secs.at(4).has("structure_ref"), "add_bonds without structure_ref writes none");
        const std::vector<Finding> fs = br.validate_findings();
        CHECK(!has_errors(fs), "bonds with structure_ref validate");
    }

    // --- C ABI smoke: create → add → write via the extern "C" surface ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "smoke");
        CHECK(cw != nullptr, "qvf_create returns a handle");
        const char* syms[2] = {"O", "H"};
        int Z[2] = {8, 1};
        double pos[6] = {0, 0, 0, 0, 0, 0.96};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 2) == 0, "C: add_structure ok");
        double f[1] = {1595.0}, in[1] = {67.0};
        CHECK(qvf_add_spectrum_xy(cw, "spectra.ir", f, in, 1) == 0, "C: spectrum ok");
        CHECK(qvf_add_spectrum_json(
                  cw, "spectra.epr",
                  "{\"g_tensor\":{\"principal\":[2.0,2.0,2.0]}}") == 0,
              "C: spectrum_json ok");
        CHECK(qvf_set_provenance_json(cw, "{\"method\":\"RHF\"}") == 0,
              "C: provenance_json ok");
        const char* path = "c_abi_smoke.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        std::ifstream f2(path, std::ios::binary);
        char magic[4] = {0};
        f2.read(magic, 4);
        CHECK(magic[0] == 'P' && magic[1] == 'K', "C: output is a zip");
        f2.close();
        std::remove(path);
        qvf_destroy(cw);
    }

    // --- Json::erase ---
    {
        Json o = Json::object();
        o.set("a", 1).set("b", 2);
        CHECK(o.erase("a") && !o.has("a") && o.has("b") && o.size() == 1,
              "Json::erase removes the key");
        CHECK(!o.erase("missing"), "Json::erase of an absent key is false");
        Json n(1);
        CHECK(!n.erase("a"), "Json::erase on a non-object is a no-op");
        Json dup = Json::parse("{\"k\":1,\"x\":2,\"k\":3}");
        CHECK(dup.erase("k") && dup.size() == 1 && dup.has("x"),
              "Json::erase removes a repeated key entirely");
    }

    // --- provenance.hostname is deprecated: dropped unless allowed (§ 3.2) ---
    {
        Json prov = Json::object();
        prov.set("method", "RHF").set("hostname", "example-host");
        QvfWriter hw({"t", "1", "provenance"});

        hw.set_provenance(prov);
        std::string text = hw.build_manifest().dump(2);
        CHECK(text.find("\"hostname\"") == std::string::npos,
              "provenance hostname dropped by default");
        CHECK(text.find("\"RHF\"") != std::string::npos,
              "other provenance fields kept");

        hw.set_provenance(prov, /*allow_hostname=*/true);
        text = hw.build_manifest().dump(2);
        CHECK(text.find("\"hostname\"") != std::string::npos &&
                  text.find("example-host") != std::string::npos,
              "provenance hostname kept with allow_hostname");

        // A pre-serialized fragment is filtered too ...
        hw.set_provenance(Json::raw(
            "{\"method\":\"RHF\",\"hostname\":\"example-host\"}"));
        text = hw.build_manifest().dump(2);
        CHECK(text.find("\"hostname\"") == std::string::npos &&
                  text.find("\"RHF\"") != std::string::npos,
              "raw provenance hostname dropped by default");
        // ... stays verbatim when there is nothing to drop ...
        hw.set_provenance(Json::raw("{\"method\":\"RHF\"}"));
        text = hw.build_manifest().dump(2);
        CHECK(text.find("{\"method\":\"RHF\"}") != std::string::npos,
              "raw provenance without hostname emitted verbatim");
        // ... and is rejected when it is not JSON at all.
        bool threw = false;
        try {
            hw.set_provenance(Json::raw("{\"method\":"));
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(threw, "malformed raw provenance rejected");
    }

    // --- C ABI: qvf_set_provenance_json inherits the hostname drop ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "provenance");
        const char* syms[1] = {"H"};
        int Z[1] = {1};
        double pos[3] = {0, 0, 0};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 1) == 0, "C: add_structure ok");
        CHECK(qvf_set_provenance_json(cw, "{\"method\":\"RHF\",\"hostname\":"
                                          "\"example-host\"}") == 0,
              "C: provenance_json with hostname ok");
        CHECK(qvf_set_provenance_json(cw, "{\"method\":") != 0,
              "C: malformed provenance_json is an error");
        CHECK(std::string(qvf_last_error(cw)).find("provenance") !=
                  std::string::npos,
              "C: malformed provenance_json error message");
        const char* path = "c_abi_hostname.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        qvf_destroy(cw);
        try {
            QvfReader r = QvfReader::open(path);
            std::string text = r.manifest().dump(2);
            CHECK(text.find("\"hostname\"") == std::string::npos,
                  "C: provenance hostname dropped");
            CHECK(r.manifest().at("provenance").at("method").as_string() == "RHF",
                  "C: other provenance fields kept");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading hostname archive: ") + e.what());
        }
        std::remove(path);
    }

    // --- INFLATE round-trips DEFLATE (in-memory) ---
    {
        std::vector<uint8_t> src;
        for (int i = 0; i < 2000; ++i)
            src.push_back(static_cast<uint8_t>("QVF round-trip "[i % 15]));
        auto comp = detail::deflate(src);
        auto back = detail::inflate(comp, src.size());
        CHECK(back == src, "inflate(deflate(x)) == x for repeated text");
        std::vector<uint8_t> zeros(4096, 0);
        CHECK(detail::inflate(detail::deflate(zeros), zeros.size()) == zeros,
              "inflate(deflate(zeros)) == zeros");
    }

    // --- INFLATE of a hand-built stored block. LEN = 0x0102 has different low
    //     and high bytes, so reading them in the wrong order fails. ---
    {
        std::vector<uint8_t> payload(0x0102);
        for (size_t i = 0; i < payload.size(); ++i)
            payload[i] = static_cast<uint8_t>(i * 7 + 1);
        std::vector<uint8_t> stored = {0x01,        // BFINAL=1, BTYPE=00
                                       0x02, 0x01,  // LEN  = 0x0102 (LE)
                                       0xFD, 0xFE}; // NLEN = ~LEN   (LE)
        stored.insert(stored.end(), payload.begin(), payload.end());
        CHECK(detail::inflate(stored, payload.size()) == payload,
              "inflate of a stored block (LEN byte order)");
    }

    // --- QvfReader round-trip: write then read back, verify values ---
    {
        QvfWriter rw({"rt", "1.0", "roundtrip"});
        rw.add_structure({{"O", 8, {{0.0, 0.0, 0.1173}}},
                          {"H", 1, {{0.0, 0.7572, -0.4692}}}});
        Spectrum s;
        s.frequencies = {1595.0, 3657.0};
        s.intensities = {67.0, 5.0};
        rw.add_spectrum("spectra.ir", s);
        rw.add_atom_properties({-0.68, 0.34});
        std::vector<uint8_t> archive = rw.to_bytes();

        QvfReader reader = QvfReader::from_bytes(archive);
        CHECK(reader.manifest().at("qvf_version").as_int() == 1,
              "reader: qvf_version == 1");
        CHECK(reader.has_section("structure"), "reader: has structure section");
        // JSON member: decode the structure and check an atom.
        Json st = reader.read_json_member("structure", "structure");
        CHECK(st.at("atoms").at(0).at("symbol").as_string() == "O",
              "reader: first atom is O");
        CHECK(std::abs(st.at("atoms").at(1).at("position").at(1).as_double() -
                       0.7572) < 1e-9, "reader: decoded H y-position");
        // JSON member: the spectrum.
        Json spec = reader.read_json_member("spectra_ir", "spectrum");
        CHECK(std::abs(spec.at("frequencies").at(0).as_double() - 1595.0) < 1e-9,
              "reader: decoded IR frequency");
        // Binary member: Mulliken charges (float64 [2]).
        BinaryMember bm = reader.read_binary_member("atom_properties",
                                                    "mulliken_charge");
        CHECK(bm.dtype == "float64" && bm.shape.size() == 1 && bm.shape[0] == 2,
              "reader: binary dtype/shape");
        auto charges = bm.as_doubles();
        CHECK(charges.size() == 2 && std::abs(charges[0] + 0.68) < 1e-9,
              "reader: decoded Mulliken charge");
        // sha256 verification: corrupting a byte must be caught. Tamper with the
        // stored bytes and expect a mismatch on verified read.
        std::vector<std::string> errs;
        CHECK(reader.validate(errs), "reader: validate() clean on good archive");
    }

    // --- Lattice: the two lattice conventions (spec § 4.3 worked example) ---
    //
    // A hexagonal cell with a = 2.46 Å and a 15.875 Å normal, as
    // structure.lattice_vectors (rows, Å) and as the equivalent
    // reaction.path.lattice (columns, bohr, 4 decimals as printed in the spec).
    {
        const Lattice::Mat3 rows = {{{{2.46, 0.0, 0.0}},
                                     {{1.23, 2.1304, 0.0}},
                                     {{0.0, 0.0, 15.875}}}};
        const Lattice::Mat3 cols_spec = {{{{4.6487, 2.3244, 0.0}},
                                          {{0.0, 4.0259, 0.0}},
                                          {{0.0, 0.0, 29.9994}}}};
        auto near = [](double a, double b, double tol) {
            return std::abs(a - b) <= tol;
        };
        auto mat_near = [&](const Lattice::Mat3& a, const Lattice::Mat3& b,
                            double tol) {
            for (size_t i = 0; i < 3; ++i)
                for (size_t j = 0; j < 3; ++j)
                    if (!near(a[i][j], b[i][j], tol)) return false;
            return true;
        };

        CHECK(kBohrPerAngstrom == 1.8897261246, "kBohrPerAngstrom value");
        CHECK(near(kBohrPerAngstrom, 1.0 / 0.529177210903, 1e-10),
              "kBohrPerAngstrom is the CODATA 2018 reciprocal");

        const Lattice cell = Lattice::from_rows_angstrom(rows);
        CHECK(cell.rows_angstrom() == rows, "Lattice rows_angstrom round-trips exactly");
        CHECK(mat_near(cell.columns_bohr(), cols_spec, 5e-5),
              "Lattice columns_bohr matches the spec worked example");
        // Column j is row j scaled: the bohr matrix is the Å matrix transposed.
        {
            const auto cols = cell.columns_bohr();
            bool scaled = true;
            for (size_t i = 0; i < 3; ++i)
                for (size_t j = 0; j < 3; ++j)
                    scaled = scaled && near(cols[i][j], kBohrPerAngstrom * rows[j][i],
                                            1e-12);
            CHECK(scaled, "Lattice columns are the scaled rows");
            CHECK(!near(cols[0][1], cols[1][0], 1e-3),
                  "worked example is not symmetric (transpose is visible)");
        }

        // Round trips, both directions, 1e-9.
        CHECK(mat_near(Lattice::from_columns_bohr(cell.columns_bohr()).rows_angstrom(),
                       rows, 1e-9),
              "Lattice rows -> columns -> rows");
        CHECK(mat_near(Lattice::from_columns_bohr(cols_spec).rows_angstrom(), rows,
                       1e-4),
              "Lattice from_columns_bohr(spec) gives the Å rows");
        {
            const Lattice::Mat3 cols = {{{{3.5, -1.25, 0.5}},
                                         {{0.75, 6.0, -2.0}},
                                         {{-0.5, 1.5, 9.25}}}};
            CHECK(mat_near(Lattice::from_columns_bohr(cols).columns_bohr(), cols,
                           1e-9),
                  "Lattice columns -> rows -> columns");
        }

        // Volume: |det| of a triangular cell is the product of its diagonal.
        const double v = 2.46 * 2.1304 * 15.875;
        CHECK(near(cell.volume_angstrom3(), v, 1e-9), "Lattice volume_angstrom3");
        CHECK(near(cell.volume_bohr3(),
                   v * kBohrPerAngstrom * kBohrPerAngstrom * kBohrPerAngstrom, 1e-8),
              "Lattice volume_bohr3");
        CHECK(near(Lattice::from_rows_angstrom({{{{0, 2, 0}}, {{1, 0, 0}}, {{0, 0, 3}}}})
                       .volume_angstrom3(),
                   6.0, 1e-12),
              "Lattice volume is positive for a left-handed cell");

        // Non-finite entries are rejected by both constructors.
        for (double bad : {std::nan(""), HUGE_VAL, -HUGE_VAL}) {
            Lattice::Mat3 m = rows;
            m[1][1] = bad;
            bool threw_rows = false, threw_cols = false;
            try { Lattice::from_rows_angstrom(m); } catch (const std::exception&) {
                threw_rows = true;
            }
            try { Lattice::from_columns_bohr(m); } catch (const std::exception&) {
                threw_cols = true;
            }
            CHECK(threw_rows, "Lattice::from_rows_angstrom rejects non-finite");
            CHECK(threw_cols, "Lattice::from_columns_bohr rejects non-finite");
        }

        // Writer overloads, checked by reading the archive back.
        QvfWriter lw({"lattice-test", "1.0", "lattice"});
        const std::vector<Atom> c_atom = {{"C", 6, {{0, 0, 0}}}};
        const std::array<bool, 3> slab_pbc = {{true, true, false}};
        const std::vector<std::array<double, 3>> row_vectors(rows.begin(), rows.end());
        lw.add_structure(c_atom, slab_pbc, row_vectors, {}, "slab_rows");
        lw.add_structure(c_atom, slab_pbc, cell, {}, "slab_cell");

        // The metadata payload must name the atoms (ReactionPathMetadataPayload).
        Json meta = Json::object();
        meta.set("symbols", Json::from_strings({"C"}));
        const Tensor coords = tensor_f64(std::vector<double>(9, 0.0), {3, 1, 3});
        lw.add_reaction_path(meta, coords, &cell, "path_cell");
        lw.add_reaction_path(meta, coords, nullptr, "path_nullptr");       // compiles, molecular
        lw.add_reaction_path(meta, coords, static_cast<const Lattice*>(nullptr),
                             "path_null_lattice");                         // molecular
        std::vector<double> raw;                                           // Tensor overload
        for (const auto& r : cols_spec) raw.insert(raw.end(), r.begin(), r.end());
        const Tensor raw_cols = tensor_f64(raw, {3, 3});
        lw.add_reaction_path(meta, coords, &raw_cols, "path_tensor");

        QvfReader lr = QvfReader::from_bytes(lw.to_bytes());
        std::vector<std::string> lerrs;
        CHECK(lr.validate(lerrs), "Lattice archive: validate() clean");
        for (const auto& e : lerrs) std::cerr << "  " << e << "\n";

        // add_structure(Lattice) writes exactly the payload of the row vectors.
        CHECK(lr.read_member_bytes("slab_rows", "structure") ==
              lr.read_member_bytes("slab_cell", "structure"),
              "add_structure(Lattice) payload equals the row-vector payload");
        {
            Json st = lr.read_json_member("slab_cell", "structure");
            bool ok = st.at("lattice_vectors").size() == 3;
            for (size_t i = 0; ok && i < 3; ++i)
                for (size_t j = 0; j < 3; ++j)
                    ok = ok && near(st.at("lattice_vectors").at(i).at(j).as_double(),
                                    rows[i][j], 1e-12);
            CHECK(ok, "structure.lattice_vectors read back as the Å rows");
            CHECK(st.at("pbc").at(0).as_bool() && st.at("pbc").at(1).as_bool() &&
                  !st.at("pbc").at(2).as_bool(), "slab pbc read back");
        }

        // add_reaction_path(Lattice*) writes the bohr column matrix.
        {
            BinaryMember bm = lr.read_binary_member("path_cell", "lattice");
            CHECK(bm.dtype == "float64" && bm.shape.size() == 2 && bm.shape[0] == 3 &&
                  bm.shape[1] == 3, "reaction.path lattice is float64 [3, 3]");
            const auto d = bm.as_doubles();
            bool spec_ok = d.size() == 9, exact_ok = d.size() == 9;
            for (size_t i = 0; spec_ok && i < 3; ++i)
                for (size_t j = 0; j < 3; ++j) {
                    spec_ok = spec_ok && near(d[3 * i + j], cols_spec[i][j], 5e-5);
                    exact_ok = exact_ok &&
                               near(d[3 * i + j], kBohrPerAngstrom * rows[j][i], 1e-12);
                }
            CHECK(spec_ok, "reaction.path lattice matches the spec bohr columns");
            CHECK(exact_ok, "reaction.path lattice is the transposed Å rows times f");
        }
        // The existing Tensor overload still writes its argument verbatim, and a
        // null pointer (either type, or a literal nullptr) means a molecular path.
        {
            const auto d = lr.read_binary_member("path_tensor", "lattice").as_doubles();
            CHECK(d == raw, "add_reaction_path(Tensor*) writes the tensor unchanged");
            CHECK(!lr.section("path_nullptr").at("members").has("lattice"),
                  "add_reaction_path(nullptr) writes no lattice");
            CHECK(!lr.section("path_null_lattice").at("members").has("lattice"),
                  "add_reaction_path((const Lattice*)nullptr) writes no lattice");
        }
    }

    // Returns true when `fn` throws a std::exception.
    auto throws = [](auto&& fn) {
        try {
            fn();
        } catch (const std::exception&) {
            return true;
        }
        return false;
    };

    // --- phonon eigenvectors and projected phonon DOS (spec § 5.6) ---
    {
        const int64_t nq = 3, nm = 6, na = 2;
        std::vector<double> freq(static_cast<size_t>(nq * nm));
        for (size_t i = 0; i < freq.size(); ++i) freq[i] = 10.0 * static_cast<double>(i);
        std::vector<double> eig(static_cast<size_t>(nq * nm * na * 3 * 2));
        for (size_t i = 0; i < eig.size(); ++i) eig[i] = 0.001 * static_cast<double>(i);
        const Tensor freq_t = tensor_f64(freq, {nq, nm});
        const Tensor eig_complex = tensor_f64(eig, {nq, nm, na, 3, 2});
        const std::vector<double> eig_gamma(static_cast<size_t>(nm * na * 3), 0.5);
        const Tensor freq_gamma =
            tensor_f64(std::vector<double>(static_cast<size_t>(nm), 1.0), {1, nm});
        const Tensor eig_real = tensor_f64(eig_gamma, {1, nm, na, 3});
        Json qpath = Json::object();
        qpath.set("n_atoms", na).set("n_modes", nm).set("has_eigenvectors", true);

        Json meta = Json::object();
        Json channels = Json::array();
        {
            Json c = Json::object();
            c.set("atom", 0).set("label", "Mg");
            channels.push_back(c);
        }
        {
            Json c = Json::object();
            c.set("atom", 1).set("label", "O");
            channels.push_back(c);
        }
        meta.set("channels", channels);
        const Tensor grid = tensor_f64({0, 100, 200, 300}, {4});
        const Tensor dos = tensor_f64({0.0, 1.0, 2.0, 0.5}, {4});
        const std::vector<double> proj = {0.0, 0.4, 0.8, 0.2, 0.0, 0.6, 1.2, 0.3};
        const Tensor proj_t = tensor_f64(proj, {2, 4});

        QvfWriter pw({"phonon-test", "1.0", "phonons"});
        pw.add_structure({{"Mg", 12, {{0, 0, 0}}}, {"O", 8, {{2.1, 2.1, 2.1}}}},
                         {{true, true, true}},
                         {{{4.2, 0, 0}}, {{0, 4.2, 0}}, {{0, 0, 4.2}}});
        pw.add_phonon_bands(qpath, freq_t, &eig_complex);
        pw.add_phonon_bands(qpath, freq_gamma, &eig_real, "phonon_gamma");
        pw.add_phonon_dos(meta, grid, dos, &proj_t);

        QvfReader pr = QvfReader::from_bytes(pw.to_bytes());
        std::vector<std::string> perrs;
        CHECK(pr.validate(perrs), "phonon archive: validate() clean");
        BinaryMember ev = pr.read_binary_member("phonon_bands", "eigenvectors");
        CHECK(ev.dtype == "float64" &&
                  ev.shape == (std::vector<int64_t>{nq, nm, na, 3, 2}),
              "complex phonon eigenvectors are float64 [n_q, n_modes, n_atoms, 3, 2]");
        CHECK(ev.as_doubles() == eig, "complex phonon eigenvectors read back exactly");
        BinaryMember ev_gamma = pr.read_binary_member("phonon_gamma", "eigenvectors");
        CHECK(ev_gamma.dtype == "float64" &&
                  ev_gamma.shape == (std::vector<int64_t>{1, nm, na, 3}),
              "real phonon eigenvectors are float64 [n_q, n_modes, n_atoms, 3]");
        CHECK(ev_gamma.as_doubles() == eig_gamma, "real phonon eigenvectors read back");
        BinaryMember pp = pr.read_binary_member("phonon_dos", "projected");
        CHECK(pp.dtype == "float64" && pp.shape == (std::vector<int64_t>{2, 4}),
              "projected phonon DOS is float64 [n_channels, n_points]");
        CHECK(pp.as_doubles() == proj, "projected phonon DOS read back exactly");
        CHECK(pr.read_json_member("phonon_dos", "meta")
                      .at("channels").at(1).at("label").as_string() == "O",
              "phonon DOS channel labels read back");

        // Every bad input throws before anything is written.
        QvfWriter bw({"phonon-test", "1.0", "rejects"});
        const Tensor rank3 = tensor_f64(std::vector<double>(36, 0.0), {nq, nm, na});
        const Tensor axis2 = tensor_f64(std::vector<double>(72, 0.0), {nq, nm, na, 2});
        const Tensor axis33 =
            tensor_f64(std::vector<double>(324, 0.0), {nq, nm, na, 3, 3});
        const Tensor eig_f32 = tensor_f32(eig, {nq, nm, na, 3, 2});
        const Tensor other_q =
            tensor_f64(std::vector<double>(144, 0.0), {nq + 1, nm, na, 3});
        const Tensor freq_flat = tensor_f64(freq, {nq * nm});
        CHECK(throws([&] { bw.add_phonon_bands(qpath, freq_t, &rank3); }),
              "phonon eigenvectors of rank 3 rejected");
        CHECK(throws([&] { bw.add_phonon_bands(qpath, freq_t, &axis2); }),
              "phonon eigenvectors with a trailing axis of 2 rejected");
        CHECK(throws([&] { bw.add_phonon_bands(qpath, freq_t, &axis33); }),
              "phonon eigenvectors with a complex axis of 3 rejected");
        CHECK(throws([&] { bw.add_phonon_bands(qpath, freq_t, &eig_f32); }),
              "float32 phonon eigenvectors rejected");
        CHECK(throws([&] { bw.add_phonon_bands(qpath, freq_t, &other_q); }),
              "phonon eigenvectors with another n_q rejected");
        CHECK(throws([&] { bw.add_phonon_bands(qpath, freq_flat, &eig_complex); }),
              "phonon eigenvectors with 1-D frequencies rejected");
        const Tensor proj_1d = tensor_f64({0, 1, 2, 3}, {4});
        const Tensor proj_f32 = tensor_f32(proj, {2, 4});
        const Tensor proj_3pts = tensor_f64({0, 1, 2, 3, 4, 5}, {2, 3});
        const Tensor proj_3rows = tensor_f64(std::vector<double>(12, 0.0), {3, 4});
        CHECK(throws([&] { bw.add_phonon_dos(meta, grid, dos, &proj_1d); }),
              "projected phonon DOS of rank 1 rejected");
        CHECK(throws([&] { bw.add_phonon_dos(meta, grid, dos, &proj_f32); }),
              "float32 projected phonon DOS rejected");
        CHECK(throws([&] { bw.add_phonon_dos(meta, grid, dos, &proj_3pts); }),
              "projected phonon DOS on another grid rejected");
        CHECK(throws([&] { bw.add_phonon_dos(meta, grid, dos, &proj_3rows); }),
              "projected phonon DOS with more rows than channels rejected");
        CHECK(bw.build_manifest().at("sections").size() == 0,
              "rejected phonon sections leave no section behind");
        // ... and no member file either: the same ids and paths still work.
        CHECK(!throws([&] { bw.add_phonon_bands(qpath, freq_t, &eig_complex); }) &&
                  !throws([&] { bw.add_phonon_dos(meta, grid, dos, &proj_t); }),
              "rejected phonon sections leave no member behind");
    }

    // --- atom_properties: one role list across every writer form (spec § 5.7) ---
    {
        const std::vector<double> q_m = {-0.7, 0.35, 0.35};
        const std::vector<double> q_l = {-0.5, 0.25, 0.25};
        const std::vector<double> q_h = {-0.3, 0.15, 0.15};
        const std::vector<double> q_i = {-0.6, 0.3, 0.3};
        const std::vector<double> s_p = {0.0, 0.0, 0.0};
        const std::vector<std::string> schema_order = {
            "mulliken_charge", "loewdin_charge", "hirshfeld_charge", "iao_charge",
            "spin_population"};
        auto roles_of = [](const Json& section) {
            std::vector<std::string> roles;
            for (const auto& kv : section.at("members").items()) roles.push_back(kv.first);
            return roles;
        };

        // The map form takes all five roles and writes them in schema order.
        QvfWriter aw({"atom-props", "1.0", "all roles"});
        aw.add_structure({{"O", 8, {{0, 0, 0.117}}}, {"H", 1, {{0, 0.757, -0.469}}},
                          {"H", 1, {{0, -0.757, -0.469}}}});
        aw.add_atom_properties({{"spin_population", s_p}, {"iao_charge", q_i},
                                {"hirshfeld_charge", q_h}, {"loewdin_charge", q_l},
                                {"mulliken_charge", q_m}});
        QvfReader ar = QvfReader::from_bytes(aw.to_bytes());
        std::vector<std::string> aerrs;
        CHECK(ar.validate(aerrs), "atom_properties archive: validate() clean");
        CHECK(roles_of(ar.section("atom_properties")) == schema_order,
              "map add_atom_properties writes the five roles in schema order");
        CHECK(ar.read_binary_member("atom_properties", "hirshfeld_charge").as_doubles() ==
                  q_h,
              "hirshfeld_charge read back");
        CHECK(ar.read_binary_member("atom_properties", "iao_charge").as_doubles() == q_i,
              "iao_charge read back");

        // The three-vector form keeps its member order; add_atom_property
        // appends to the section it created.
        QvfWriter vw({"atom-props", "1.0", "append"});
        vw.add_atom_properties(q_m, {}, s_p);
        vw.add_atom_property("iao_charge", q_i);
        vw.add_atom_property("hirshfeld_charge", q_h);
        QvfReader vr = QvfReader::from_bytes(vw.to_bytes());
        std::vector<std::string> verrs;
        CHECK(vr.validate(verrs), "appended atom_properties archive: validate() clean");
        CHECK(vr.sections().size() == 1, "add_atom_property adds no second section");
        CHECK(roles_of(vr.section("atom_properties")) ==
                  (std::vector<std::string>{"mulliken_charge", "spin_population",
                                            "iao_charge", "hirshfeld_charge"}),
              "add_atom_property appends to the existing section");
        CHECK(vr.read_binary_member("atom_properties", "iao_charge").as_doubles() == q_i,
              "appended iao_charge read back");

        // Every form rejects what the schema would.
        QvfWriter ew({"atom-props", "1.0", "rejects"});
        CHECK(throws([&] { ew.add_atom_properties(); }),
              "atom_properties without a role rejected");
        CHECK(throws([&] {
                  ew.add_atom_properties(std::map<std::string, std::vector<double>>{});
              }),
              "atom_properties with an empty role map rejected");
        CHECK(throws([&] { ew.add_atom_properties({{"bogus_charge", q_m}}); }),
              "atom_properties role outside the schema list rejected");
        CHECK(throws([&] { ew.add_atom_properties({{"mulliken_charge", {}}}); }),
              "atom_properties role without values rejected");
        CHECK(throws([&] {
                  ew.add_atom_properties({{"mulliken_charge", q_m},
                                          {"iao_charge", {0.1, -0.1}}});
              }),
              "atom_properties roles of different lengths rejected");
        CHECK(throws([&] { ew.add_atom_properties(q_m, {0.1, -0.1}); }),
              "three-vector atom_properties of different lengths rejected");
        CHECK(throws([&] { ew.add_atom_property("charge", q_m); }),
              "add_atom_property role outside the schema list rejected");
        CHECK(ew.build_manifest().at("sections").size() == 0,
              "rejected atom_properties leave no section behind");
        ew.add_atom_property("mulliken_charge", q_m);
        CHECK(throws([&] { ew.add_atom_property("mulliken_charge", q_m); }),
              "add_atom_property of a role already present rejected");
        CHECK(throws([&] { ew.add_atom_property("iao_charge", {0.1, -0.1}); }),
              "add_atom_property with another atom count rejected");
        ew.add_bonds({{0, 1, 1.0}}, "bonds");
        CHECK(throws([&] { ew.add_atom_property("iao_charge", q_i, "bonds"); }),
              "add_atom_property into a section of another kind rejected");
    }

    // --- C ABI: qvf_add_atom_property creates, then appends (spec § 5.7) ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "atom properties");
        const char* syms[2] = {"O", "H"};
        int Z[2] = {8, 1};
        double pos[6] = {0, 0, 0, 0, 0, 0.96};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 2) == 0, "C: add_structure ok");
        double mul[2] = {-0.34, 0.34}, hir[2] = {-0.2, 0.2}, iao[2] = {-0.4, 0.4};
        CHECK(qvf_add_atom_properties(cw, mul, nullptr, nullptr, 2) == 0,
              "C: add_atom_properties ok");
        CHECK(qvf_add_atom_property(cw, "hirshfeld_charge", hir, 2) == 0,
              "C: add_atom_property hirshfeld_charge ok");
        CHECK(qvf_add_atom_property(cw, "iao_charge", iao, 2) == 0,
              "C: add_atom_property iao_charge ok");
        CHECK(qvf_add_atom_property(cw, "iao_charge", iao, 2) != 0,
              "C: add_atom_property of a role already present is an error");
        CHECK(qvf_add_atom_property(cw, "bogus_charge", iao, 2) != 0 &&
                  std::string(qvf_last_error(cw)).find("bogus_charge") !=
                      std::string::npos,
              "C: add_atom_property of an unknown role is an error");
        CHECK(qvf_add_atom_property(cw, "loewdin_charge", iao, 1) != 0,
              "C: add_atom_property with another n_atoms is an error");
        CHECK(qvf_add_atom_property(cw, nullptr, iao, 2) != 0,
              "C: add_atom_property with a NULL role is an error");
        CHECK(qvf_add_atom_property(cw, "loewdin_charge", nullptr, 2) != 0,
              "C: add_atom_property with NULL values is an error");
        const char* path = "c_abi_atom_property.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        qvf_destroy(cw);
        try {
            QvfReader r = QvfReader::open(path);
            std::vector<std::string> errs;
            CHECK(r.validate(errs), "C: atom_properties archive validates");
            CHECK(r.section("atom_properties").at("members").size() == 3,
                  "C: one atom_properties section with three roles");
            CHECK(r.read_binary_member("atom_properties", "iao_charge").as_doubles() ==
                      (std::vector<double>{-0.4, 0.4}),
                  "C: iao_charge read back");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading atom_properties archive: ") + e.what());
        }
        std::remove(path);

    }

    // --- atom_properties: every scalar and vector role, custom roles (spec § 5.7) ---
    {
        const std::vector<std::string> scalar_roles = {
            "mulliken_charge", "loewdin_charge", "hirshfeld_charge", "iao_charge",
            "spin_population", "cm5_charge", "chelpg_charge", "resp_charge",
            "esp_charge", "nbo_charge", "mbis_charge", "bader_charge", "ddec_charge",
            "valence", "bonded_valence", "fukui_plus", "fukui_minus", "fukui_zero",
            "dual_descriptor", "atomic_energy_eh", "atomic_volume_bohr3",
            "nmr_shielding_iso_ppm"};
        const std::vector<std::string> vector_roles = {"magnetic_moment",
                                                       "force_eh_bohr"};
        std::vector<std::string> schema_order = scalar_roles;
        schema_order.insert(schema_order.end(), vector_roles.begin(), vector_roles.end());
        auto roles_of = [](const Json& section) {
            std::vector<std::string> roles;
            for (const auto& kv : section.at("members").items()) roles.push_back(kv.first);
            return roles;
        };
        const std::vector<std::array<double, 3>> moments = {{{0, 0, 1.0}}, {{0, 0, -0.5}}};

        // Every scalar role through the map form, vector roles appended: the
        // section lists them in schema order or in append order.
        QvfWriter fw({"atom-props", "1.0", "every role"});
        fw.add_structure({{"O", 8, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.96}}}});
        std::map<std::string, std::vector<double>> all;
        for (size_t i = 0; i < scalar_roles.size(); ++i)
            all[scalar_roles[i]] = {0.1 * static_cast<double>(i), -0.1 * static_cast<double>(i)};
        fw.add_atom_properties(all);
        fw.add_atom_vector_property("magnetic_moment", moments);
        fw.add_atom_vector_property("force_eh_bohr", {{{0, 0, 0.01}}, {{0, 0, -0.01}}});
        QvfReader fr = QvfReader::from_bytes(fw.to_bytes());
        std::vector<std::string> ferrs;
        CHECK(fr.validate(ferrs), "every-role atom_properties archive: validate() clean");
        CHECK(roles_of(fr.section("atom_properties")) == schema_order,
              "every role written in schema order");
        for (const auto& role : vector_roles) {
            const Json& spec = fr.section("atom_properties").at("members").at(role);
            CHECK(spec.at("dtype").as_string() == "float64" &&
                      spec.at("shape").size() == 2 &&
                      spec.at("shape").at(0).as_int() == 2 &&
                      spec.at("shape").at(1).as_int() == 3,
                  role + " is float64 [n_atoms, 3]");
        }
        CHECK(fr.read_binary_member("atom_properties", "magnetic_moment").as_doubles() ==
                  (std::vector<double>{0, 0, 1.0, 0, 0, -0.5}),
              "magnetic_moment read back row-major");
        CHECK(fr.read_binary_member("atom_properties", "nmr_shielding_iso_ppm")
                      .as_doubles() == all["nmr_shielding_iso_ppm"],
              "nmr_shielding_iso_ppm read back");

        // Custom roles: described first, written with a meta member.
        const std::string meta_json =
            R"({"custom.flux": {"units": "e", "description": "Net flux.", "kind": "scalar"},)"
            R"( "custom.axis": {"units": "dimensionless", "description": "Axis.", "kind": "vector",)"
            R"( "x_note": 1}})";
        QvfWriter cw2({"atom-props", "1.0", "custom"});
        cw2.add_structure({{"O", 8, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.96}}}});
        CHECK(throws([&] { cw2.add_atom_property("custom.flux", {0.1, -0.1}); }),
              "a custom role without a meta entry is rejected");
        CHECK(cw2.build_manifest().at("sections").size() == 1,
              "the rejected custom role left no section behind");
        CHECK(throws([&] { cw2.set_atom_properties_meta(Json::raw("{not json")); }),
              "set_atom_properties_meta rejects invalid JSON");
        CHECK(throws([&] { cw2.set_atom_properties_meta(Json::raw("[]")); }),
              "set_atom_properties_meta rejects a non-object");
        CHECK(throws([&] {
                  cw2.set_atom_properties_meta(Json::raw(
                      R"({"flux": {"units": "e", "description": "d", "kind": "scalar"}})"));
              }),
              "set_atom_properties_meta rejects a name without custom.");
        CHECK(throws([&] {
                  cw2.set_atom_properties_meta(Json::raw(
                      R"({"custom.a b": {"units": "e", "description": "d", "kind": "scalar"}})"));
              }),
              "set_atom_properties_meta rejects a bad custom name");
        for (const char* bad : {
                 R"({"custom.a": {"description": "d", "kind": "scalar"}})",
                 R"({"custom.a": {"units": "", "description": "d", "kind": "scalar"}})",
                 R"({"custom.a": {"units": "e", "kind": "scalar"}})",
                 R"({"custom.a": {"units": "e", "description": "d"}})",
                 R"({"custom.a": {"units": "e", "description": "d", "kind": "tensor"}})",
                 R"({"custom.a": {"units": "e", "description": "d", "kind": "scalar", "extra": 1}})",
                 R"({"custom.a": "scalar"})"})
            CHECK(throws([&] { cw2.set_atom_properties_meta(Json::raw(bad)); }),
                  std::string("set_atom_properties_meta rejects ") + bad);
        CHECK(throws([&] { cw2.add_atom_property("custom.a", {0.1, -0.1}); }),
              "a rejected meta changes nothing");
        cw2.set_atom_properties_meta(Json::raw(meta_json));
        cw2.add_atom_properties({{"mulliken_charge", {0.3, -0.3}},
                                 {"custom.flux", {0.9, -0.9}}});
        CHECK(throws([&] { cw2.add_atom_property("custom.axis", {0.1, -0.1}); }),
              "a vector custom role is rejected as a scalar");
        CHECK(throws([&] {
                  cw2.add_atom_vector_property("custom.flux", moments);
              }),
              "a scalar custom role is rejected as a vector");
        cw2.add_atom_vector_property("custom.axis", moments);
        QvfReader cr = QvfReader::from_bytes(cw2.to_bytes());
        std::vector<std::string> cerrs;
        CHECK(cr.validate(cerrs), "custom-role archive: validate() clean");
        CHECK(roles_of(cr.section("atom_properties")) ==
                  (std::vector<std::string>{"mulliken_charge", "custom.flux",
                                            "custom.axis", "meta"}),
              "custom roles follow the listed roles, then meta");
        Json meta = cr.read_json_member("atom_properties", "meta");
        CHECK(meta.size() == 2 && meta.at("custom.flux").at("units").as_string() == "e" &&
                  meta.at("custom.axis").at("kind").as_string() == "vector" &&
                  meta.at("custom.axis").at("x_note").as_int() == 1,
              "meta member holds both entries");
        CHECK(cr.read_binary_member("atom_properties", "custom.flux").as_doubles() ==
                  (std::vector<double>{0.9, -0.9}),
              "custom.flux read back");

        // A custom role added to a section that has no meta yet creates it; a
        // second one replaces it, and meta stays the last member.
        QvfWriter iw({"atom-props", "1.0", "incremental"});
        iw.add_structure({{"O", 8, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.96}}}});
        iw.set_atom_properties_meta(Json::raw(meta_json));
        iw.add_atom_property("mulliken_charge", {0.3, -0.3});
        iw.add_atom_property("custom.flux", {0.9, -0.9});
        iw.add_atom_property("hirshfeld_charge", {0.2, -0.2});
        iw.add_atom_vector_property("custom.axis", moments);
        QvfReader ir = QvfReader::from_bytes(iw.to_bytes());
        std::vector<std::string> ierrs;
        CHECK(ir.validate(ierrs), "incremental custom-role archive: validate() clean");
        CHECK(roles_of(ir.section("atom_properties")) ==
                  (std::vector<std::string>{"mulliken_charge", "custom.flux",
                                            "hirshfeld_charge", "custom.axis", "meta"}),
              "meta is rewritten and stays the last member");
        Json meta2 = ir.read_json_member("atom_properties", "meta");
        CHECK(meta2.size() == 2 && meta2.has("custom.flux") && meta2.has("custom.axis"),
              "the rewritten meta describes both custom roles");
        // Only the final meta file is in the archive (no orphan from the rewrite).
        bool orphan = false;
        for (const auto& f : ir.validate_findings())
            orphan = orphan || f.code == "W_ORPHAN_MEMBER";
        CHECK(!orphan, "no orphan meta.json entry after the rewrite");

        // Listed roles still need no meta; the meta of a section only lists its roles.
        QvfWriter lw({"atom-props", "1.0", "unused meta"});
        lw.add_structure({{"O", 8, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.96}}}});
        lw.set_atom_properties_meta(Json::raw(meta_json));
        lw.add_atom_properties({{"custom.flux", {0.9, -0.9}}});
        QvfReader lr = QvfReader::from_bytes(lw.to_bytes());
        CHECK(lr.read_json_member("atom_properties", "meta").size() == 1,
              "meta lists only the custom roles the section carries");

        // Role / shape mismatches and bad vectors.
        QvfWriter sw({"atom-props", "1.0", "shapes"});
        CHECK(throws([&] { sw.add_atom_property("magnetic_moment", {0.1, 0.2}); }),
              "a vector role is rejected as a scalar");
        CHECK(throws([&] { sw.add_atom_properties({{"force_eh_bohr", {0.1, 0.2}}}); }),
              "a vector role is rejected in the scalar map");
        CHECK(throws([&] { sw.add_atom_vector_property("mulliken_charge", moments); }),
              "a scalar role is rejected as a vector");
        CHECK(throws([&] { sw.add_atom_vector_property("bogus_vector", moments); }),
              "an unknown vector role is rejected");
        CHECK(throws([&] { sw.add_atom_vector_property("magnetic_moment", {}); }),
              "a vector role without values is rejected");
        sw.add_atom_property("mulliken_charge", {0.1, -0.1});
        CHECK(throws([&] {
                  sw.add_atom_vector_property("magnetic_moment", {{{0, 0, 1.0}}});
              }),
              "a vector role of another atom count is rejected");
        CHECK(sw.build_manifest().at("sections").at(0).at("members").size() == 1,
              "rejected roles leave the section unchanged");
        sw.add_atom_vector_property("magnetic_moment", moments);
        CHECK(throws([&] { sw.add_atom_vector_property("magnetic_moment", moments); }),
              "a vector role already present is rejected");
    }

    // --- C ABI: vector roles and custom roles (spec § 5.7) ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "atom vectors");
        const char* syms[2] = {"O", "H"};
        int Z[2] = {8, 1};
        double pos[6] = {0, 0, 0, 0, 0, 0.96};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 2) == 0, "C: add_structure ok");
        double moments[6] = {0, 0, 1.0, 0, 0, -0.5};
        double flux[2] = {0.9, -0.9};
        double cm5[2] = {-0.3, 0.3};
        CHECK(qvf_add_atom_property(cw, "cm5_charge", cm5, 2) == 0,
              "C: add_atom_property cm5_charge ok");
        CHECK(qvf_add_atom_vector_property(cw, "magnetic_moment", moments, 2) == 0,
              "C: add_atom_vector_property magnetic_moment ok");
        CHECK(qvf_add_atom_vector_property(cw, "magnetic_moment", moments, 2) != 0,
              "C: a vector role already present is an error");
        CHECK(qvf_add_atom_vector_property(cw, "cm5_charge", moments, 2) != 0,
              "C: a scalar role as a vector is an error");
        CHECK(qvf_add_atom_property(cw, "force_eh_bohr", flux, 2) != 0,
              "C: a vector role as a scalar is an error");
        CHECK(qvf_add_atom_vector_property(cw, "force_eh_bohr", moments, 1) != 0,
              "C: a vector role with another n_atoms is an error");
        CHECK(qvf_add_atom_vector_property(cw, nullptr, moments, 2) != 0,
              "C: a vector role with a NULL role is an error");
        CHECK(qvf_add_atom_vector_property(cw, "force_eh_bohr", nullptr, 2) != 0,
              "C: a vector role with NULL values is an error");
        CHECK(qvf_add_atom_vector_property(cw, "force_eh_bohr", moments, 0) != 0,
              "C: a vector role with n_atoms 0 is an error");
        CHECK(qvf_add_atom_property(cw, "custom.flux", flux, 2) != 0 &&
                  std::string(qvf_last_error(cw)).find("custom.flux") != std::string::npos,
              "C: a custom role without meta is an error");
        CHECK(qvf_set_atom_properties_meta_json(cw, "{not json") != 0,
              "C: invalid meta JSON is an error");
        CHECK(qvf_set_atom_properties_meta_json(cw, nullptr) != 0,
              "C: NULL meta JSON is an error");
        CHECK(qvf_set_atom_properties_meta_json(
                  cw, "{\"custom.flux\": {\"units\": \"e\", \"description\": \"Flux\", "
                      "\"kind\": \"scalar\"}, \"custom.axis\": {\"units\": \"debye\", "
                      "\"description\": \"Axis\", \"kind\": \"vector\"}}") == 0,
              "C: set_atom_properties_meta_json ok");
        CHECK(qvf_add_atom_property(cw, "custom.flux", flux, 2) == 0,
              "C: custom scalar role ok");
        CHECK(qvf_add_atom_vector_property(cw, "custom.axis", moments, 2) == 0,
              "C: custom vector role ok");
        const char* path = "c_abi_atom_vector.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        qvf_destroy(cw);
        try {
            QvfReader r = QvfReader::open(path);
            std::vector<std::string> errs;
            CHECK(r.validate(errs), "C: vector/custom archive validates");
            CHECK(r.section("atom_properties").at("members").size() == 5,
                  "C: cm5_charge, magnetic_moment, two custom roles and meta");
            CHECK(r.read_binary_member("atom_properties", "custom.axis").as_doubles() ==
                      (std::vector<double>{0, 0, 1.0, 0, 0, -0.5}),
                  "C: custom.axis read back");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading vector archive: ") + e.what());
        }
        std::remove(path);
    }

    // --- Spin axis on bands / dos.projected projections (spec § 5.6) ---
    {
        auto throws = [](auto fn) {
            try {
                fn();
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        QvfWriter sw({"spin-test", "1.0", "spin"});
        Json kp = Json::object();
        kp.set("n_spin", 2).set("n_kpoints", 3).set("n_bands", 2);
        const Tensor eig = tensor_f64(std::vector<double>(12, 0.0), {2, 3, 2});
        const Tensor proj4 = tensor_f64(std::vector<double>(24, 0.5), {2, 3, 2, 2});
        const Tensor occ = tensor_f64(std::vector<double>(12, 1.0), {2, 3, 2});
        sw.add_bands(kp, eig, &proj4, &occ, "bands_spin");
        Json kp1 = Json::object();
        kp1.set("n_spin", 1).set("n_kpoints", 3).set("n_bands", 2);
        const Tensor eig1 = tensor_f64(std::vector<double>(6, 0.0), {1, 3, 2});
        const Tensor proj3 = tensor_f64(std::vector<double>(12, 0.5), {3, 2, 2});
        sw.add_bands(kp1, eig1, &proj3, "bands_legacy");       // rank 3: n_spin = 1
        sw.add_bands(kp, eig, nullptr, nullptr, "bands_plain");

        const Tensor bad_spin = tensor_f64(std::vector<double>(12, 0.5), {1, 3, 2, 2});
        CHECK(throws([&] { sw.add_bands(kp, eig, &bad_spin, nullptr, "bad1"); }),
              "add_bands: rank-4 projections with another n_spin throws");
        const Tensor bad_k = tensor_f64(std::vector<double>(8, 0.5), {2, 2, 2});
        CHECK(throws([&] { sw.add_bands(kp1, eig1, &bad_k, "bad2"); }),
              "add_bands: rank-3 projections with another n_k throws");
        const Tensor bad_occ = tensor_f64(std::vector<double>(6, 1.0), {1, 3, 2});
        CHECK(throws([&] { sw.add_bands(kp, eig, nullptr, &bad_occ, "bad3"); }),
              "add_bands: occupations with another shape throws");
        const Tensor f32 = tensor_f32(std::vector<double>(24, 0.5), {2, 3, 2, 2});
        CHECK(throws([&] { sw.add_bands(kp, eig, &f32, nullptr, "bad4"); }),
              "add_bands: float32 projections throw");

        Json dm2 = Json::object();
        dm2.set("n_spin", 2);
        const Tensor de = tensor_f64({-1, 0, 1}, {3});
        const Tensor dp3 = tensor_f64(std::vector<double>(12, 0.1), {2, 2, 3});
        sw.add_dos_projected(de, dp3, dm2, "pdos_spin");
        Json dm1 = Json::object();
        dm1.set("n_spin", 1);
        CHECK(throws([&] { sw.add_dos_projected(de, dp3, dm1, "bad5"); }),
              "add_dos_projected: spin axis differing from n_spin throws");
        CHECK(throws([&] {
                  sw.add_dos_projected(de, tensor_f64({0, 0, 0}, {3}), dm1, "bad6");
              }),
              "add_dos_projected: rank-1 projections throw");

        const std::vector<uint8_t> bytes = sw.to_bytes();
        CHECK(!has_errors(QvfReader::validate_bytes(bytes)),
              "spin-axis archive validates without errors");
        QvfReader sr = QvfReader::from_bytes(bytes);
        const BinaryMember bp = sr.read_binary_member("bands_spin", "projections");
        CHECK(bp.dtype == "float64" && bp.shape == std::vector<int64_t>({2, 3, 2, 2}),
              "rank-4 projections written as float64 [2, 3, 2, 2]");
        const BinaryMember bo = sr.read_binary_member("bands_spin", "occupations");
        CHECK(bo.shape == std::vector<int64_t>({2, 3, 2}) && bo.as_doubles()[5] == 1.0,
              "occupations written as [n_spin, n_k, n_bands]");
        CHECK(!sr.section("bands_plain").at("members").has("occupations") &&
                  !sr.section("bands_plain").at("members").has("projections"),
              "add_bands with null pointers writes neither member");
        CHECK(sr.read_binary_member("pdos_spin", "projections").shape.size() == 3,
              "rank-3 dos.projected projections written");
        CHECK(!sr.has_section("bad1") && !sr.has_section("bad5"),
              "a refused section leaves nothing behind");
    }

    // --- unit vocabulary and declared units (QVF spec § 4.3) ---
    {
        auto near = [](double a, double b, double rel) {
            return std::fabs(a - b) <= rel * std::fabs(b);
        };
        auto throws = [](auto&& fn) {
            try {
                fn();
            } catch (const std::exception&) {
                return true;
            }
            return false;
        };
        CHECK(kUnitCount == 45 && unit_name(Unit::Wavenumber) == "cm-1" &&
                  unit_name(Unit::Dimensionless) == "dimensionless" &&
                  static_cast<std::size_t>(Unit::Dimensionless) == 24,
              "unit_name");
        CHECK(unit_name(Unit::HartreePerBohr) == "hartree_per_bohr" &&
                  unit_name(Unit::ReciprocalAngstrom) == "angstrom-1" &&
                  parse_unit("e_per_bohr3") == Unit::ElementaryChargePerBohr3 &&
                  !parse_unit("hartree/bohr") && !parse_unit("e_per_angstrom3"),
              "compound unit names are closed names");
        bool round_trip = true;
        for (std::size_t i = 0; i < kUnitCount; ++i) {
            const Unit u = static_cast<Unit>(i);
            round_trip = round_trip && parse_unit(unit_name(u)) == u;
        }
        CHECK(round_trip, "parse_unit(unit_name(u)) == u for every unit");
        CHECK(!parse_unit("eV") && !parse_unit("Bohr") && !parse_unit(""),
              "parse_unit is case-sensitive and rejects legacy spellings");
        CHECK(canonical_unit("Eh") == Unit::Hartree &&
                  canonical_unit("Hartree") == Unit::Hartree &&
                  canonical_unit("eV") == Unit::ElectronVolt &&
                  canonical_unit("ev") == Unit::ElectronVolt && !canonical_unit("rad"),
              "canonical_unit maps the legacy spellings");
        CHECK(convert(1.0, Unit::Angstrom, Unit::Bohr) == kBohrPerAngstrom,
              "convert angstrom -> bohr");
        CHECK(near(convert(1.0, Unit::Hartree, Unit::ElectronVolt), 27.211386245988, 1e-15),
              "convert hartree -> ev");
        CHECK(near(convert(1.0, Unit::Hartree, Unit::KcalPerMol), 627.5094740631, 1e-12),
              "convert hartree -> kcal_per_mol");
        CHECK(near(convert(1.0, Unit::KcalPerMol, Unit::KjPerMol), 4.184, 1e-14),
              "convert kcal_per_mol -> kj_per_mol");
        CHECK(near(convert(1.0, Unit::Hartree, Unit::Wavenumber), 219474.6313632, 1e-15),
              "convert hartree -> cm-1");
        CHECK(near(convert(300.0, Unit::Nanometer, Unit::ElectronVolt), 4.13280661444, 1e-11),
              "convert nm -> ev");
        CHECK(near(convert(20000.0, Unit::Wavenumber, Unit::Nanometer), 500.0, 1e-11),
              "convert cm-1 -> nm");
        CHECK(convert(5.0, Unit::Kelvin, Unit::Kelvin) == 5.0, "convert to itself");
        CHECK(!can_convert(Unit::Angstrom, Unit::ElectronVolt) &&
                  !can_convert(Unit::Bohr, Unit::Nanometer) &&
                  throws([] { convert(1.0, Unit::Kelvin, Unit::Atmosphere); }),
              "unsupported conversions are refused");

        // Writers declare units and refuse names outside the vocabulary.
        QvfWriter uw(Source{"t", "1", "units"});
        uw.add_structure({Atom{"O", 8, {{0, 0, 0}}}});
        Grid g;
        g.units = {{"origin", "bohr"}, {"voxel_vectors", "bohr"}};
        uw.add_volume("volume.density", g, tensor_f32(std::vector<double>(8, 0.1), {2, 2, 2}),
                      "", "", "rho");
        Json meta = Json::object();
        Json tu = Json::object();
        tu.set("energies", "hartree").set("coords", "angstrom");
        meta.set("symbols", Json::from_strings({"O"})).set("units", tu);
        uw.add_trajectory(meta, tensor_f64({0, 0, 0}, {1, 1, 3}), "traj");
        Spectrum sp;
        sp.frequencies = {200.0};
        sp.intensities = {0.1};
        sp.x_units = "nm";
        sp.y_units = "dimensionless";
        uw.add_spectrum("spectra.uvvis", sp, "", "uv");
        Grid bad = g;
        bad.units["origin"] = "angstrom";
        CHECK(throws([&] {
                  uw.add_volume("volume.density", bad,
                                tensor_f32(std::vector<double>(8, 0.1), {2, 2, 2}), "", "",
                                "bad1");
              }),
              "grid units.origin must restate bohr");
        Grid unknown = g;
        unknown.units["data"] = "e_per_angstrom3";
        CHECK(throws([&] {
                  uw.add_volume("volume.density", unknown,
                                tensor_f32(std::vector<double>(8, 0.1), {2, 2, 2}), "", "",
                                "bad2");
              }),
              "an unknown grid unit name throws");
        Json bad_meta = Json::object();
        Json bu = Json::object();
        bu.set("energies", "Eh");
        bad_meta.set("symbols", Json::from_strings({"O"})).set("units", bu);
        CHECK(throws([&] { uw.add_trajectory(bad_meta, tensor_f64({0, 0, 0}, {1, 1, 3}), "bad3"); }),
              "a legacy spelling is not a declared unit");
        Spectrum bad_sp = sp;
        bad_sp.x_units = "1/cm";
        CHECK(throws([&] { uw.add_spectrum("spectra.ir", bad_sp, "", "bad4"); }),
              "an unknown spectrum x_units throws");
        CHECK(throws([&] {
                  uw.add_spectrum("spectra.ir",
                                  Json::raw(R"({"frequencies":[1],"intensities":[1],)"
                                            R"("broadening":{"type":"none","units":"x"}})"),
                                  "", "bad5");
              }),
              "an unknown broadening unit in a raw payload throws");
        Json qt = Json::object();
        Json qu = Json::object();
        qu.set("position", "au");
        qt.set("points", Json::array()).set("units", qu);
        CHECK(throws([&] { uw.add_topology_qtaim(qt, "bad6"); }),
              "an unknown QTAIM unit throws");

        const std::vector<uint8_t> bytes = uw.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        bool unit_warning = false;
        for (const auto& f : fs) unit_warning = unit_warning || f.code == "W_UNIT_UNKNOWN";
        CHECK(!has_errors(fs) && !unit_warning, "declared-units archive validates cleanly");
        QvfReader ur = QvfReader::from_bytes(bytes);
        CHECK(ur.read_json_member("rho", "grid").at("units").at("origin").as_string() == "bohr",
              "grid units written");
        CHECK(ur.read_json_member("uv", "spectrum").at("x_units").as_string() == "nm",
              "spectrum x_units written");
        for (const char* sid : {"bad1", "bad2", "bad3", "bad4", "bad5", "bad6"})
            CHECK(!ur.has_section(sid), std::string("refused section left nothing: ") + sid);
    }

    // --- forces, hessian and the root energies block (spec § 3.6, 5.10, 5.11) ---
    {
        auto throws = [](auto fn) {
            try {
                fn();
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        QvfWriter dw({"deriv-test", "1.0", "derivatives"});
        dw.add_structure({{"O", 8, {{0.0, 0.0, 0.117}}},
                          {"H", 1, {{0.0, 0.757, -0.469}}},
                          {"H", 1, {{0.0, -0.757, -0.469}}}});
        const Tensor grad = tensor_f64({0, 0, 0.02, 0, 0.01, -0.01, 0, -0.01, -0.01}, {3, 3});
        const Tensor stress = tensor_f64({-1e-4, 0, 0, 0, -1e-4, 0, 0, 0, -2e-4}, {3, 3});
        Json units = Json::object();
        units.set("energy", "hartree").set("length", "bohr");
        Json fmeta = Json::object();
        fmeta.set("structure_ref", "structure").set("method", "rhf")
             .set("convention", "gradient").set("units", units);
        CHECK(dw.add_forces(grad, &stress, fmeta) == "forces", "add_forces default id");
        std::vector<double> h(81, 0.0);
        for (int k = 0; k < 9; ++k) h[static_cast<size_t>(k * 9 + k)] = 0.5;
        const Tensor hess = tensor_f64(h, {9, 9});
        Json hmeta = Json::object();
        hmeta.set("structure_ref", "structure").set("mass_weighted", false)
             .set("projected", false)
             .set("atom_masses_amu", Json::from_doubles({15.999, 1.008, 1.008}));
        CHECK(dw.add_hessian(hess, hmeta) == "hessian", "add_hessian default id");
        Json vmeta = Json::object();
        vmeta.set("frequencies", Json::from_doubles({1650.0, 3800.0, 3900.0}));
        dw.add_vibrations(vmeta, tensor_f64(std::vector<double>(27, 0.0), {3, 3, 3}),
                          "modes", "hessian");
        Json e_total = Json::object();
        e_total.set("value", -74.9659011).set("units", "hartree");
        Json prov = Json::object();
        prov.set("method", "RHF").set("scf_energy", e_total);
        dw.set_provenance(prov);
        Json energies = Json::object();
        energies.set("total", e_total).set("method", "rhf");
        dw.set_energies(energies);

        CHECK(throws([&] { dw.add_forces(tensor_f64({0, 0, 0}, {3}), nullptr, Json(), "b1"); }),
              "add_forces: rank-1 gradient throws");
        CHECK(throws([&] {
                  dw.add_forces(tensor_f64(std::vector<double>(6, 0.0), {3, 2}), nullptr,
                                Json(), "b2");
              }),
              "add_forces: [3, 2] gradient throws");
        CHECK(throws([&] {
                  dw.add_forces(tensor_f32(std::vector<double>(9, 0.0), {3, 3}), nullptr,
                                Json(), "b3");
              }),
              "add_forces: float32 gradient throws");
        const Tensor bad_stress = tensor_f64(std::vector<double>(6, 0.0), {2, 3});
        CHECK(throws([&] { dw.add_forces(grad, &bad_stress, Json(), "b4"); }),
              "add_forces: [2, 3] stress throws");
        CHECK(throws([&] { dw.add_forces(grad, nullptr, Json("x"), "b5"); }),
              "add_forces: non-object meta throws");
        CHECK(throws([&] {
                  dw.add_hessian(tensor_f64(std::vector<double>(54, 0.0), {9, 6}), Json(),
                                 "b6");
              }),
              "add_hessian: non-square matrix throws");
        CHECK(throws([&] {
                  dw.add_hessian(tensor_f64(std::vector<double>(64, 0.0), {8, 8}), Json(),
                                 "b7");
              }),
              "add_hessian: extent not divisible by 3 throws");
        Json weighted = Json::object();
        weighted.set("mass_weighted", true);
        CHECK(throws([&] { dw.add_hessian(hess, weighted, "b8"); }),
              "add_hessian: mass-weighted without masses throws");
        Json two_masses = Json::object();
        two_masses.set("atom_masses_amu", Json::from_doubles({1.0, 2.0}));
        CHECK(throws([&] { dw.add_hessian(hess, two_masses, "b9"); }),
              "add_hessian: masses not one per atom throws");
        CHECK(throws([&] {
                  dw.add_vibrations(vmeta,
                                    tensor_f64(std::vector<double>(27, 0.0), {3, 3, 3}),
                                    "b10", "structure");
              }),
              "add_vibrations: hessian_ref naming a non-hessian throws");
        Json no_total = Json::object();
        no_total.set("method", "rhf");
        CHECK(throws([&] { dw.set_energies(no_total); }), "set_energies without total throws");
        CHECK(throws([&] { dw.set_energies(Json(1.0)); }), "set_energies of a number throws");

        const std::vector<uint8_t> bytes = dw.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        CHECK(!has_errors(fs), "derivatives archive validates without errors");
        QvfReader dr = QvfReader::from_bytes(bytes);
        CHECK(dr.read_binary_member("forces", "gradient").as_doubles()[2] == 0.02,
              "gradient read back");
        CHECK(dr.read_binary_member("forces", "stress").shape ==
                  std::vector<int64_t>({3, 3}),
              "stress written as [3, 3]");
        CHECK(dr.read_binary_member("hessian", "matrix").shape ==
                  std::vector<int64_t>({9, 9}),
              "matrix written as [9, 9]");
        CHECK(dr.read_json_member("hessian", "meta").at("atom_masses_amu").size() == 3,
              "hessian meta read back");
        CHECK(dr.section("modes").at("hessian_ref").as_string() == "hessian",
              "vibrations carries hessian_ref");
        CHECK(dr.manifest().at("energies").at("method").as_string() == "rhf",
              "root energies block written");
        for (const char* id : {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10"})
            CHECK(!dr.has_section(id), std::string("refused section left nothing: ") + id);

        // An energies block that disagrees with provenance.scf_energy.
        QvfWriter bw({"deriv-test", "1.0", "inconsistent"});
        bw.add_structure({{"H", 1, {{0, 0, 0}}}});
        bw.set_provenance(prov);
        Json off = Json::object();
        off.set("value", -74.96).set("units", "hartree");
        Json bad_energies = Json::object();
        bad_energies.set("total", off).set("method", "Rhf");
        bw.set_energies(bad_energies);
        bool flagged = false;
        for (const auto& f : QvfReader::validate_bytes(bw.to_bytes()))
            if (f.code == "E_ROOT_INCONSISTENT" && f.where == "/energies/total") flagged = true;
        CHECK(flagged, "energies.total disagreeing with scf_energy is E_ROOT_INCONSISTENT");
    }

    // --- C ABI: forces, hessian, energies ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "derivatives");
        const char* syms[2] = {"H", "H"};
        int Z[2] = {1, 1};
        double pos[6] = {0, 0, 0, 0, 0, 0.74};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 2) == 0, "C: add_structure ok");
        double grad[6] = {0, 0, -0.01, 0, 0, 0.01};
        double stress[9] = {-1e-4, 0, 0, 0, -1e-4, 0, 0, 0, -1e-4};
        CHECK(qvf_add_forces(cw, 2, grad, stress) == 0, "C: add_forces ok");
        CHECK(qvf_add_forces(cw, 0, grad, nullptr) != 0, "C: add_forces with 0 atoms fails");
        CHECK(qvf_add_forces(cw, 2, nullptr, nullptr) != 0, "C: add_forces NULL fails");
        double hess[36] = {0};
        for (int k = 0; k < 6; ++k) hess[k * 6 + k] = 0.37;
        CHECK(qvf_add_hessian(cw, 2, hess) == 0, "C: add_hessian ok");
        CHECK(qvf_add_hessian(cw, -1, hess) != 0, "C: add_hessian with n_atoms < 1 fails");
        CHECK(qvf_set_provenance_json(
                  cw, "{\"method\":\"rhf\",\"scf_energy\":{\"value\":-1.1167,"
                      "\"units\":\"hartree\"}}") == 0,
              "C: provenance ok");
        CHECK(qvf_set_energies_json(
                  cw, "{\"total\":{\"value\":-1.1167,\"units\":\"hartree\"},"
                      "\"method\":\"rhf\"}") == 0,
              "C: set_energies_json ok");
        CHECK(qvf_set_energies_json(cw, "") != 0, "C: empty energies_json fails");
        const char* path = "c_abi_derivatives.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        qvf_destroy(cw);
        try {
            QvfReader r = QvfReader::open(path);
            std::vector<std::string> errs;
            CHECK(r.validate(errs), "C: derivatives archive validates");
            CHECK(r.read_binary_member("forces", "gradient").shape ==
                      std::vector<int64_t>({2, 3}),
                  "C: gradient is [2, 3]");
            CHECK(r.read_binary_member("forces", "stress").as_doubles()[8] == -1e-4,
                  "C: stress read back");
            CHECK(r.read_binary_member("hessian", "matrix").shape ==
                      std::vector<int64_t>({6, 6}),
                  "C: matrix is [6, 6]");
            CHECK(r.manifest().at("energies").at("total").at("value").as_double() == -1.1167,
                  "C: energies written verbatim");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading derivatives archive: ") + e.what());
        }

    }

    // --- excited states, NTOs and a transition density (spec § 5.2-5.4, 5.12) ---
    {
        auto throws = [](auto fn) {
            try {
                fn();
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        QvfWriter ew({"es-test", "1.0", "excited states"});
        ew.add_structure({{"O", 8, {{0.0, 0.0, 0.117}}},
                          {"H", 1, {{0.0, 0.757, -0.469}}},
                          {"H", 1, {{0.0, -0.757, -0.469}}}});
        WavefunctionGTO wf;
        wf.shells = {{0, 0, {1.0}, {1.0}, true}, {1, 0, {1.0}, {1.0}, true}};
        wf.mo_coefficients = {2, 2, {1, 0, 0, 1}};
        wf.energies = {-0.5, 0.5};
        wf.occupations = {2, 0};
        CHECK(ew.add_wavefunction_gto(wf) == "wf", "reference wavefunction added");
        WavefunctionGTO hole = wf;
        hole.orbital_kind = "nto";
        hole.mo_coefficients = {1, 2, {1, 0}};
        hole.energies = {};
        hole.occupations = {1.0};
        hole.state_ref = "es";
        hole.state_index = 1;
        hole.nto_role = "hole";
        hole.section_id = "nto_hole";
        ew.add_wavefunction_gto(hole);
        WavefunctionGTO particle = hole;
        particle.mo_coefficients = {1, 2, {0, 1}};
        particle.nto_role = "particle";
        particle.section_id = "nto_particle";
        ew.add_wavefunction_gto(particle);
        Grid g;
        const Tensor dens = tensor_f32(std::vector<double>(8, 0.125), {2, 2, 2});
        CHECK(ew.add_volume_transition_density(g, dens, "es", 1, "S1", "tden") == "tden",
              "add_volume_transition_density returns its id");

        const Json states = Json::parse(
            R"([{"index":1,"excitation_energy":7.4,"total_energy":-74.69,)"
            R"("multiplicity":1,"symmetry":"B1","oscillator_strength":0.05,)"
            R"("rotatory_strength_cgs":0.0,"transition_dipole":[0.5,0,0],)"
            R"("transition_magnetic_dipole":[0,0.2,0],)"
            R"("contributions":[{"from_mo":0,"to_mo":1,"weight":0.99}],)"
            R"("nto_ref":"nto_hole","density_ref":"tden","x_note":"demo"},)"
            R"({"index":2,"excitation_energy":{"value":9.2,"units":"ev"}}])");
        const Json meta = Json::parse(
            R"({"method":"tddft","n_states":2,"tda":false,)"
            R"("units":{"excitation_energy":"ev","oscillator_strength":"dimensionless"}})");
        auto bad = [&](const std::string& states_json, const char* id,
                       const std::string& wf_ref = "wf", const Json& m = Json()) {
            return throws([&] { ew.add_excited_states(Json::parse(states_json), m, wf_ref, id); });
        };
        const std::string e1 = R"("excitation_energy":7.4)";
        CHECK(bad(R"([{"index":1,)" + e1 + R"(},{"index":1,)" + e1 + "}]", "x1"),
              "add_excited_states: a repeated index throws");
        CHECK(bad(R"([{"index":2,)" + e1 + R"(},{"index":1,)" + e1 + "}]", "x2"),
              "add_excited_states: descending indices throw");
        CHECK(bad(R"([{"index":0,)" + e1 + "}]", "x3"), "add_excited_states: index 0 throws");
        CHECK(bad(R"([{"index":1}])", "x4"),
              "add_excited_states: a state without excitation_energy throws");
        CHECK(bad(R"([{"index":1,)" + e1 + R"(,"energy_ev":7.4}])", "x5"),
              "add_excited_states: an unknown state field throws");
        CHECK(bad(R"([{"index":1,)" + e1 + R"(,"nto_ref":"structure"}])", "x6"),
              "add_excited_states: nto_ref naming a structure throws");
        CHECK(bad(R"([{"index":1,)" + e1 + R"(,"density_ref":"wf"}])", "x7"),
              "add_excited_states: density_ref naming a wavefunction throws");
        CHECK(bad(R"([{"index":1,)" + e1 +
                      R"(,"contributions":[{"from_mo":0,"to_mo":2,"weight":1}]}])",
                  "x8"),
              "add_excited_states: to_mo beyond the reference's rows throws");
        CHECK(bad(R"([{"index":1,)" + e1 +
                      R"(,"contributions":[{"from_mo":0,"to_mo":1,"spin":"up","weight":1}]}])",
                  "x9"),
              "add_excited_states: a spin other than alpha / beta throws");
        CHECK(bad(R"([{"index":1,)" + e1 + "}]", "x10", "structure"),
              "add_excited_states: wavefunction_ref naming a structure throws");
        CHECK(bad(R"([{"index":1,)" + e1 + "}]", "x11", "wf",
                  Json::parse(R"({"n_states":3})")),
              "add_excited_states: meta.n_states disagreeing throws");
        CHECK(bad(R"([{"index":1,)" + e1 + "}]", "x12", "wf",
                  Json::parse(R"({"units":{"oscillator_strength":"ev"}})")),
              "add_excited_states: a contradicted fixed unit throws");
        CHECK(bad(R"([{"index":1,"excitation_energy":{"value":7.4,"units":"eV"}}])", "x13"),
              "add_excited_states: an energy unit that is not a unit name throws");
        CHECK(throws([&] { ew.add_excited_states(Json::object(), Json(), "", "x14"); }),
              "add_excited_states: states that are not an array throw");
        CHECK(ew.add_excited_states(states, meta, "wf", "es") == "es",
              "add_excited_states returns its id");
        const Json uv = Json::parse(
            R"({"frequencies":[7.4,9.2],"intensities":[0.05,0.0],)"
            R"("peaks":[{"x":7.4,"assignment":{"state_index":1}}],"states_ref":"es"})");
        ew.add_spectrum("spectra.uvvis", uv);
        CHECK(throws([&] {
                  ew.add_spectrum("spectra.uvvis",
                                  Json::parse(R"({"frequencies":[1],"intensities":[1],)"
                                              R"("states_ref":"wf"})"),
                                  "", "x15");
              }),
              "add_spectrum: states_ref naming a wavefunction throws");
        WavefunctionGTO bad_role = hole;
        bad_role.orbital_kind = "canonical";
        bad_role.section_id = "x16";
        CHECK(throws([&] { ew.add_wavefunction_gto(bad_role); }),
              "add_wavefunction_gto: nto_role without orbital_kind nto throws");
        WavefunctionGTO no_ref = hole;
        no_ref.state_ref = "";
        no_ref.section_id = "x17";
        CHECK(throws([&] { ew.add_wavefunction_gto(no_ref); }),
              "add_wavefunction_gto: state_index without state_ref throws");
        CHECK(throws([&] { ew.add_volume_transition_density(g, dens, "", 1, "", "x18"); }),
              "add_volume_transition_density: state_index without a reference throws");

        const std::vector<uint8_t> bytes = ew.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        CHECK(!has_errors(fs), "excited-states archive validates without errors");
        QvfReader er = QvfReader::from_bytes(bytes);
        const Json payload = er.read_json_member("es", "states");
        CHECK(payload.at("states").at(size_t{0}).at("excitation_energy").at("units")
                      .as_string() == "ev",
              "a number excitation_energy is written in ev");
        CHECK(payload.at("states").at(size_t{0}).at("total_energy").at("units")
                      .as_string() == "hartree",
              "a number total_energy is written in hartree");
        CHECK(payload.at("states").at(size_t{0}).at("x_note").as_string() == "demo",
              "x_ keys of a state are kept");
        CHECK(payload.at("meta").at("n_states").as_int() == 2, "meta written");
        CHECK(er.section("es").at("wavefunction_ref").as_string() == "wf",
              "section carries wavefunction_ref");
        CHECK(er.section("tden").at("kind").as_string() == "volume.transition_density" &&
                  er.section("tden").at("excited_states_ref").as_string() == "es" &&
                  er.section("tden").at("state_index").as_int() == 1,
              "transition density section fields written");
        const Json hm = er.read_json_member("nto_hole", "mo_metadata");
        CHECK(hm.at("orbital_kind").as_string() == "nto" &&
                  hm.at("state_ref").as_string() == "es" &&
                  hm.at("state_index").as_int() == 1 &&
                  hm.at("nto_role").as_string() == "hole" &&
                  hm.at("occupation_semantics").as_string() == "transition_weight",
              "NTO mo_metadata fields written");
        CHECK(!er.read_json_member("wf", "mo_metadata").has("occupation_semantics"),
              "a canonical wavefunction gets no occupation_semantics");
        for (const char* id : {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10",
                               "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18"})
            CHECK(!er.has_section(id), std::string("refused section left nothing: ") + id);
    }

    // --- C ABI: excited states ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "excited states");
        const char* syms[2] = {"H", "H"};
        int Z[2] = {1, 1};
        double pos[6] = {0, 0, 0, 0, 0, 0.74};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 2) == 0, "C: add_structure ok");
        int center[2] = {0, 1}, l[2] = {0, 0}, nprim[2] = {1, 1};
        double exps[2] = {1.24, 1.24}, coeffs[2] = {1.0, 1.0};
        double C[4] = {0.5, 0.5, 0.5, -0.5};
        double eps[2] = {-0.58, 0.67}, occ[2] = {2, 0};
        CHECK(qvf_add_wavefunction_gto(cw, 2, center, l, nprim, exps, coeffs, 1, 0, 2, 2, C,
                                       eps, occ) == 0,
              "C: add_wavefunction_gto ok");
        const char* states =
            "[{\"index\":1,\"excitation_energy\":{\"value\":12.9,\"units\":\"ev\"},"
            "\"oscillator_strength\":0.55,"
            "\"contributions\":[{\"from_mo\":0,\"to_mo\":1,\"weight\":1.0}]}]";
        CHECK(qvf_add_excited_states_json(cw, nullptr, nullptr, nullptr) != 0,
              "C: NULL states_json fails");
        CHECK(qvf_add_excited_states_json(cw, "[", nullptr, nullptr) != 0,
              "C: malformed states_json fails");
        CHECK(qvf_add_excited_states_json(
                  cw,
                  "[{\"index\":1,\"excitation_energy\":12.9,"
                  "\"contributions\":[{\"from_mo\":0,\"to_mo\":5,\"weight\":1}]}]",
                  nullptr, "wf") != 0,
              "C: a contribution beyond the orbital rows fails");
        CHECK(qvf_add_excited_states_json(cw, states, "{\"method\":\"cis\",\"n_states\":1}",
                                          "wf") == 0,
              "C: add_excited_states_json ok");
        const char* path = "c_abi_excited_states.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        qvf_destroy(cw);
        try {
            QvfReader r = QvfReader::open(path);
            std::vector<std::string> errs;
            CHECK(r.validate(errs), "C: excited-states archive validates");
            CHECK(r.section("excited_states").at("wavefunction_ref").as_string() == "wf",
                  "C: wavefunction_ref written");
            const Json p = r.read_json_member("excited_states", "states");
            CHECK(p.at("states").at(size_t{0}).at("oscillator_strength").as_double() == 0.55,
                  "C: states read back");
            CHECK(p.at("meta").at("method").as_string() == "cis", "C: meta read back");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading excited-states archive: ") + e.what());
        }
        std::remove(path);
    }

    // --- root metadata and section references (spec § 3.5, § 4.1) ---
    {
        // The form of an automatic value: version-4 uuid, "YYYY-MM-DDThh:mm:ssZ".
        auto is_uuid_v4 = [](const std::string& u) {
            if (u.size() != 36) return false;
            for (size_t i = 0; i < u.size(); ++i) {
                const char c = u[i];
                if (i == 8 || i == 13 || i == 18 || i == 23) {
                    if (c != '-') return false;
                } else if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) {
                    return false;
                }
            }
            return u[14] == '4' && std::string("89ab").find(u[19]) != std::string::npos;
        };
        auto is_utc = [](const std::string& t) {
            if (t.size() != 20 || t[4] != '-' || t[7] != '-' || t[10] != 'T' ||
                t[13] != ':' || t[16] != ':' || t[19] != 'Z')
                return false;
            for (size_t i = 0; i < 19; ++i)
                if (i != 4 && i != 7 && i != 10 && i != 13 && i != 16 &&
                    (t[i] < '0' || t[i] > '9'))
                    return false;
            return t.compare(0, 2, "20") == 0;
        };
        const std::vector<Atom> water = {{"O", 8, {{0, 0, 0.117}}},
                                         {"H", 1, {{0, 0.757, -0.469}}},
                                         {"H", 1, {{0, -0.757, -0.469}}}};

        QvfWriter mw({"meta", "1.0", "auto"});
        mw.add_structure(water);
        const Json first = mw.build_manifest();
        CHECK(first.has("metadata"), "auto metadata block present by default");
        const std::string uuid = first.at("metadata").at("uuid").as_string();
        const std::string created = first.at("metadata").at("created_utc").as_string();
        CHECK(is_uuid_v4(uuid), "auto uuid is a lowercase version-4 UUID: " + uuid);
        CHECK(is_utc(created), "auto created_utc is YYYY-MM-DDThh:mm:ssZ: " + created);
        CHECK(mw.build_manifest().dump() == first.dump(),
              "automatic values are made once and reused");
        QvfWriter other({"meta", "1.0", "auto"});
        CHECK(other.build_manifest().at("metadata").at("uuid").as_string() != uuid,
              "two writers get different uuids");

        QvfWriter off({"meta", "1.0", "off"});
        off.auto_metadata = false;
        CHECK(!off.build_manifest().has("metadata"), "auto_metadata = false: no block");

        Json meta = Json::object();
        meta.set("title", "explicit").set("uuid", "6f1c2a4e-8b3d-4c5e-9f70-1a2b3c4d5e6f");
        mw.set_metadata(meta);
        const Json m2 = mw.build_manifest().at("metadata");
        CHECK(m2.at("uuid").as_string() == "6f1c2a4e-8b3d-4c5e-9f70-1a2b3c4d5e6f" &&
                  m2.at("title").as_string() == "explicit" &&
                  m2.at("created_utc").as_string() == created,
              "explicit metadata wins; the missing created_utc is still added");
        CHECK(throws([&] { mw.set_metadata(Json::array()); }),
              "set_metadata: a non-object throws");
        CHECK(throws([&] { mw.set_metadata(Json::raw("{\"title\": ")); }),
              "set_metadata: a malformed raw fragment throws");

        // Section references.
        const std::vector<Atom> ammonia_n = {{"N", 7, {{0, 0, 0.116}}}};
        mw.add_structure(ammonia_n, {{false, false, false}},
                         std::vector<std::array<double, 3>>{}, std::vector<Bond>{}, "nh3");
        RunRecord rr;
        rr.program = "meta";
        rr.log_text = "done\n";
        mw.add_run_record(rr, "run");
        Json iters = Json::array();
        Json it1 = Json::object();
        it1.set("iter", 1).set("energy_eh", -74.9);
        iters.push_back(it1);
        mw.add_scf_history(iters, "scf");
        CHECK(mw.annotate_section("scf", "structure", "run", "step 1") == "scf",
              "annotate_section returns the id");
        mw.annotate_section("scf", "", "", "step 2");     // "" leaves the refs
        mw.set_section_field("nh3", "label", "ammonia");
        CHECK(throws([&] { mw.annotate_section("scf", "run"); }),
              "structure_ref naming a run.record throws");
        CHECK(throws([&] { mw.annotate_section("scf", "nh3", "structure"); }),
              "run_ref naming a structure throws");
        CHECK(throws([&] { mw.annotate_section("scf", "nowhere"); }),
              "structure_ref naming no section throws");
        CHECK(throws([&] { mw.annotate_section("missing", "structure"); }),
              "annotating an unknown section throws");
        CHECK(throws([&] { mw.set_section_field("scf", "kind", "bonds"); }),
              "set_section_field refuses kind");
        CHECK(throws([&] { mw.set_section_field("scf", "run_ref", Json(5)); }),
              "set_section_field: a non-string run_ref throws");
        const Json sec = mw.build_manifest().at("sections").at(3);
        CHECK(sec.at("id").as_string() == "scf" &&
                  sec.at("structure_ref").as_string() == "structure" &&
                  sec.at("run_ref").as_string() == "run" &&
                  sec.at("group").as_string() == "step 2",
              "a refused annotation changes nothing (atomic)");
        const std::vector<uint8_t> bytes = mw.to_bytes();
        CHECK(!has_errors(QvfReader::validate_bytes(bytes)),
              "annotated archive with metadata validates without errors");
        QvfReader mr = QvfReader::from_bytes(bytes);
        CHECK(mr.manifest().at("metadata").at("uuid").as_string() ==
                  "6f1c2a4e-8b3d-4c5e-9f70-1a2b3c4d5e6f",
              "metadata round-trips");
        CHECK(mr.section("nh3").at("label").as_string() == "ammonia",
              "set_section_field label round-trips");

        // C ABI.
        qvf_writer* cw = qvf_create("c-abi", "1.0", "refs");
        const char* syms[1] = {"He"};
        int Z[1] = {2};
        double pos[3] = {0, 0, 0};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 1) == 0, "C: structure");
        CHECK(qvf_add_run_record(cw, "c-abi", nullptr, "ok\n", nullptr, nullptr, 0, 0,
                                 nullptr, nullptr, -1) == 0,
              "C: run record");
        CHECK(qvf_add_citations(cw, "@misc{x, title={T}}") == 0, "C: citations");
        CHECK(qvf_annotate_section(cw, "citations", "structure", "run_record", nullptr) == 0,
              "C: annotate_section");
        CHECK(qvf_annotate_section(cw, "citations", "run_record", nullptr, nullptr) != 0 &&
                  std::string(qvf_last_error(cw)).find("structure_ref") != std::string::npos,
              "C: a structure_ref of the wrong kind is an error");
        CHECK(qvf_annotate_section(cw, nullptr, nullptr, nullptr, "g") != 0,
              "C: a NULL id is an error");
        CHECK(qvf_set_metadata_json(cw, "[1]") != 0, "C: non-object metadata is an error");
        CHECK(qvf_set_metadata_json(cw, nullptr) != 0, "C: NULL metadata is an error");
        CHECK(qvf_set_metadata_json(cw, "{\"title\":\"C\",\"license\":\"MIT\"}") == 0,
              "C: set_metadata_json");
        const char* path = "c_abi_refs.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write");
        qvf_destroy(cw);
        std::vector<Finding> cfs = QvfReader::validate_file(path);
        CHECK(!has_errors(cfs), "C: annotated archive validates");
        QvfReader cr = QvfReader::open(path);
        const Json cm = cr.manifest().at("metadata");
        CHECK(cm.at("title").as_string() == "C" && is_uuid_v4(cm.at("uuid").as_string()),
              "C: metadata written with the automatic uuid added");
        CHECK(cr.section("citations").at("structure_ref").as_string() == "structure" &&
                  cr.section("citations").at("run_ref").as_string() == "run_record" &&
                  !cr.section("citations").has("group"),
              "C: references written, NULL group left unset");
        std::remove(path);
    }

    // --- typed provenance fields (spec § 3.2) and the grown job.spec payload
    //     (spec § 5.9): the writers check them, the validator resolves
    //     structure_ref ---
    {
        const std::vector<Atom> hydrogen = {{"H", 1, {{0, 0, 0}}}};
        auto software_entry = [](const std::string& name, const std::string& version) {
            Json s = Json::object();
            s.set("name", name).set("version", version);
            return s;
        };
        auto base = [&](const std::string& calc) {
            QvfWriter bw({"wp27", "1.0", calc});
            bw.auto_metadata = false;
            bw.add_structure(hydrogen);
            return bw;
        };
        auto provenance_with = [&](const std::string& key, Json value) {
            Json prov = Json::object();
            prov.set("method", "rhf").set(key, std::move(value));
            return prov;
        };
        auto error_codes = [](const std::vector<Finding>& findings) {
            std::vector<std::string> codes;
            for (const Finding& f : findings)
                if (f.severity == Finding::Severity::Error) codes.push_back(f.code);
            return codes;
        };

        // A full, valid block is written and the archive validates.
        {
            QvfWriter pw = base("prov_full");
            Json software = Json::array();
            software.push_back(software_entry("my-code", "1.2.3"));
            Json lib = software_entry("integral-lib", "2.9.0");
            lib.set("url", "https://example.org/integral-lib").set("x_build", "release");
            software.push_back(lib);
            Json prov = Json::object();
            prov.set("method", "rks").set("software", software)
                .set("git_commit", "9f2b6c1d4e7a").set("n_threads", 8)
                .set("n_mpi_ranks", 4).set("n_gpus", 0)
                .set("started_utc", "2026-10-07T12:00:00Z")
                .set("finished_utc", "2026-10-07T12:07:31.25Z")
                .set("wall_seconds", 451.25).set("x_other", "kept");
            pw.set_provenance(prov);
            const Json written = pw.build_manifest().at("provenance");
            CHECK(written.at("software").at(1).at("url").as_string() ==
                      "https://example.org/integral-lib" &&
                      written.at("n_mpi_ranks").as_int() == 4 &&
                      written.at("finished_utc").as_string() == "2026-10-07T12:07:31.25Z" &&
                      written.has("x_other"),
                  "typed provenance fields are written, other keys pass through");
            CHECK(!has_errors(QvfReader::validate_bytes(pw.to_bytes())),
                  "an archive with the typed provenance fields validates");
        }

        // A wrong field throws and the earlier block stays as it was.
        {
            QvfWriter pw = base("prov_bad");
            Json good = Json::object();
            good.set("method", "rhf");
            pw.set_provenance(good);
            const std::string before = pw.build_manifest().at("provenance").dump();
            auto sw = [&](Json entries) { return provenance_with("software", entries); };
            Json no_name = Json::array();
            no_name.push_back(Json::object().set("version", "1"));
            Json empty_name = Json::array();
            empty_name.push_back(software_entry("", "1"));
            Json bad_version = Json::array();
            bad_version.push_back(Json::object().set("name", "a").set("version", 2));
            Json extra_key = Json::array();
            extra_key.push_back(Json::object().set("name", "a").set("license", "MIT"));
            Json not_object = Json::array();
            not_object.push_back(Json("a"));
            const std::vector<std::pair<std::string, Json>> bad = {
                {"software not an array", sw(Json::object())},
                {"software item without name", sw(no_name)},
                {"software item with empty name", sw(empty_name)},
                {"software version not a string", sw(bad_version)},
                {"software item with an unknown key", sw(extra_key)},
                {"software item not an object", sw(not_object)},
                {"git_commit too short", provenance_with("git_commit", Json("9f2b6c"))},
                {"git_commit too long", provenance_with("git_commit", Json(std::string(65, 'a')))},
                {"git_commit not hex", provenance_with("git_commit", Json("9f2b6c1d4e7g"))},
                {"git_commit not a string", provenance_with("git_commit", Json(1234567))},
                {"n_threads zero", provenance_with("n_threads", Json(0))},
                {"n_threads a float", provenance_with("n_threads", Json(2.0))},
                {"n_threads a string", provenance_with("n_threads", Json("4"))},
                {"n_threads a bool", provenance_with("n_threads", Json(true))},
                {"n_mpi_ranks zero", provenance_with("n_mpi_ranks", Json(0))},
                {"n_gpus negative", provenance_with("n_gpus", Json(-1))},
                {"started_utc without Z", provenance_with("started_utc", Json("2026-10-07T12:00:00"))},
                {"started_utc with an offset",
                 provenance_with("started_utc", Json("2026-10-07T14:00:00+02:00"))},
                {"finished_utc with a space", provenance_with("finished_utc", Json("2026-10-07 12:00:00Z"))},
                {"finished_utc month 13", provenance_with("finished_utc", Json("2026-13-07T12:00:00Z"))},
                {"finished_utc hour 24", provenance_with("finished_utc", Json("2026-10-07T24:00:00Z"))},
                {"finished_utc empty fraction", provenance_with("finished_utc", Json("2026-10-07T12:00:00.Z"))},
                {"finished_utc a number", provenance_with("finished_utc", Json(1790000000))},
            };
            for (const auto& entry : bad) {
                CHECK(throws([&] { pw.set_provenance(entry.second); }),
                      "set_provenance refuses: " + entry.first);
                CHECK(pw.build_manifest().at("provenance").dump() == before,
                      "a refused provenance changes nothing: " + entry.first);
            }
            for (const char* ts : {"2016-12-31T23:59:60Z", "2026-10-07T12:00:00.25Z"})
                CHECK(!throws([&] {
                          pw.set_provenance(provenance_with("started_utc", Json(ts)));
                      }),
                      std::string("set_provenance accepts ") + ts);
            for (const char* hash : {"9f2b6c1", "9F2B6C1D4E7A", "0123456789abcdefABCDEF"})
                CHECK(!throws([&] {
                          pw.set_provenance(provenance_with("git_commit", Json(hash)));
                      }),
                      std::string("set_provenance accepts git_commit ") + hash);
            CHECK(throws([&] { pw.set_provenance(Json::array()); }),
                  "set_provenance: a non-object throws");
        }

        // A raw fragment (the C ABI path) is checked and stays verbatim.
        {
            QvfWriter pw = base("prov_raw");
            const std::string text =
                "{\"method\":\"rhf\",\"git_commit\":\"9f2b6c1d4e7a\",\"n_threads\":2}";
            pw.set_provenance(Json::raw(text));
            CHECK(pw.build_manifest().dump().find(text) != std::string::npos,
                  "a valid raw provenance is emitted verbatim");
            CHECK(throws([&] {
                      pw.set_provenance(Json::raw("{\"git_commit\":\"nothex1\"}"));
                  }),
                  "a raw provenance with a wrong typed field throws");
            CHECK(throws([&] { pw.set_provenance(Json::raw("[1]")); }),
                  "a raw provenance that is not an object throws");
            CHECK(throws([&] {
                      pw.set_provenance(Json::raw("{\"n_gpus\":-1}"), true);
                  }),
                  "a raw provenance is checked also with allow_hostname");
            // A nested raw fragment of a built object is checked too.
            Json prov = Json::object();
            prov.set("software", Json::raw("[{\"version\":\"1\"}]"));
            CHECK(throws([&] { pw.set_provenance(prov); }),
                  "a nested raw software array is checked");
            CHECK(!has_errors(QvfReader::validate_bytes(pw.to_bytes())),
                  "the archive of a verbatim raw provenance validates");
        }

        // job.spec: the new typed fields.
        auto job_full = [] {
            Json spec = Json::object();
            Json solvation = Json::object();
            solvation.set("model", "cpcm").set("solvent", "water").set("epsilon", 78.4);
            Json scf = Json::object();
            scf.set("convergence", 1e-8).set("max_iterations", 100)
               .set("guess", "sad").set("damping", 0.3);
            Json geometry = Json::object();
            geometry.set("convergence", 4.5e-4).set("max_steps", 50)
                    .set("coordinates", "redundant");
            Json resources = Json::object();
            resources.set("n_cores", 16).set("memory_mb", 32000)
                     .set("walltime_s", 86400.0);
            spec.set("job_type", "periodic").set("method", "rks")
                .set("structure_ref", "structure").set("solvation", solvation)
                .set("dispersion", "d3bj").set("relativistic", "x2c")
                .set("scf", scf).set("geometry", geometry)
                .set("resources", resources)
                .set("options", Json::object().set("cutoff_ha", 300.0));
            return spec;
        };
        {
            QvfWriter jw = base("job_full");
            jw.add_job_spec(job_full());
            Json prov = Json::object();
            prov.set("run_status", "pending");
            jw.set_provenance(prov);
            const std::vector<uint8_t> bytes = jw.to_bytes();
            CHECK(!has_errors(QvfReader::validate_bytes(bytes)),
                  "an archive with every new job.spec field validates");
            QvfReader jr = QvfReader::from_bytes(bytes);
            const Json spec = jr.read_json_member("job_spec", "spec");
            CHECK(spec.at("structure_ref").as_string() == "structure" &&
                      spec.at("scf").at("max_iterations").as_int() == 100 &&
                      spec.at("geometry").at("coordinates").as_string() == "redundant" &&
                      spec.at("resources").at("n_cores").as_int() == 16 &&
                      spec.at("solvation").at("epsilon").as_double() == 78.4,
                  "the new job.spec fields round-trip");
            CHECK(jr.section("job_spec").at("structure_ref").as_string() == "structure",
                  "the structure_ref is written on the job.spec section object too");
        }
        {
            // Each wrong value throws and adds nothing; the id stays free.
            QvfWriter jw = base("job_bad");
            const size_t sections_before = jw.build_manifest().at("sections").size();
            auto with = [&](const std::string& key, Json value) {
                Json spec = Json::object();
                spec.set("job_type", "molecular").set(key, std::move(value));
                return spec;
            };
            auto block = [](const std::string& key, Json value) {
                return Json::object().set(key, std::move(value));
            };
            const std::vector<std::pair<std::string, Json>> bad = {
                {"structure_ref empty", with("structure_ref", Json(""))},
                {"structure_ref a number", with("structure_ref", Json(5))},
                {"structure_ref not an id", with("structure_ref", Json("a/b"))},
                {"dispersion empty", with("dispersion", Json(""))},
                {"relativistic a number", with("relativistic", Json(1))},
                {"solvation a string", with("solvation", Json("pcm"))},
                {"solvation without model", with("solvation", block("solvent", Json("water")))},
                {"solvation epsilon zero",
                 with("solvation", Json::object().set("model", "pcm").set("epsilon", 0))},
                {"solvation epsilon negative",
                 with("solvation", Json::object().set("model", "pcm").set("epsilon", -1.5))},
                {"solvation unknown key",
                 with("solvation", Json::object().set("model", "pcm").set("cavity", 1))},
                {"scf convergence zero", with("scf", block("convergence", Json(0.0)))},
                {"scf max_iterations zero", with("scf", block("max_iterations", Json(0)))},
                {"scf max_iterations a float", with("scf", block("max_iterations", Json(2.5)))},
                {"scf guess empty", with("scf", block("guess", Json("")))},
                {"scf damping a string", with("scf", block("damping", Json("high")))},
                {"scf unknown key", with("scf", block("diis", Json(true)))},
                {"geometry convergence negative", with("geometry", block("convergence", Json(-1.0)))},
                {"geometry max_steps zero", with("geometry", block("max_steps", Json(0)))},
                {"geometry coordinates zmatrix", with("geometry", block("coordinates", Json("zmatrix")))},
                {"geometry an array", with("geometry", Json::array())},
                {"resources n_cores zero", with("resources", block("n_cores", Json(0)))},
                {"resources memory_mb zero", with("resources", block("memory_mb", Json(0)))},
                {"resources walltime_s zero", with("resources", block("walltime_s", Json(0.0)))},
                {"resources unknown key", with("resources", block("gpus", Json(1)))},
            };
            for (const auto& entry : bad) {
                CHECK(throws([&] { jw.add_job_spec(entry.second); }),
                      "add_job_spec refuses: " + entry.first);
                CHECK(jw.build_manifest().at("sections").size() == sections_before,
                      "a refused job.spec adds nothing: " + entry.first);
            }
            // Vendor keys pass through the blocks; the root stays open.
            Json spec = Json::object();
            spec.set("job_type", "molecular").set("x_future", 1)
                .set("scf", block("x_engine_mode", Json("direct")))
                .set("resources", block("x_queue", Json("short")));
            CHECK(!throws([&] { jw.add_job_spec(spec); }),
                  "add_job_spec accepts vendor keys in the blocks and an open root");
        }
        {
            // structure_ref: a section of another kind is refused, a later one allowed.
            QvfWriter jw = base("job_ref");
            jw.add_citations("@misc{x, title={T}}", "refs");
            Json spec = Json::object();
            spec.set("job_type", "molecular").set("structure_ref", "refs");
            CHECK(throws([&] { jw.add_job_spec(spec, "job_a"); }),
                  "structure_ref naming a citations section throws");
            spec.set("structure_ref", "structure");
            CHECK(!throws([&] { jw.add_job_spec(spec, "job_a"); }),
                  "structure_ref naming a structure is accepted");
            spec.set("structure_ref", "later");
            CHECK(!throws([&] { jw.add_job_spec(spec, "job_b"); }),
                  "a structure_ref to a section added later is accepted");
        }

        // The validator resolves structure_ref (E_REF_UNRESOLVED, spec § 5.9).
        {
            // The job.spec comes first and names a section added after it: the
            // writer accepts a reference it cannot yet check, the validator
            // resolves it against the whole archive.
            enum class Later { Nothing, Structure, Citations };
            auto job_archive = [&](const std::string& ref, Later later) {
                QvfWriter jw({"wp27", "1.0", "ref"});
                jw.auto_metadata = false;
                Json spec = Json::object();
                spec.set("job_type", "molecular").set("structure_ref", ref);
                jw.add_job_spec(spec);
                if (later == Later::Structure)
                    jw.add_structure(hydrogen, {{false, false, false}},
                                     std::vector<std::array<double, 3>>{},
                                     std::vector<Bond>{}, "later");
                if (later == Later::Citations)
                    jw.add_citations("@misc{x, title={T}}", "later");
                return QvfReader::validate_bytes(jw.to_bytes());
            };
            CHECK(!has_errors(job_archive("later", Later::Structure)),
                  "a structure_ref to a structure added later resolves");
            // The structure_ref is written on the section object too (spec
            // § 5.9, "Job sets"), so a dangling one is reported twice: at the
            // section field and at the spec member.
            const std::vector<std::string> twice{"E_REF_UNRESOLVED", "E_REF_UNRESOLVED"};
            std::vector<Finding> wrong = job_archive("later", Later::Citations);
            std::set<std::string> where_wrong;
            for (const Finding& f : wrong)
                if (f.code == "E_REF_UNRESOLVED") where_wrong.insert(f.where);
            const std::set<std::string> both{"/sections/0/members/spec",
                                             "/sections/0/structure_ref"};
            CHECK(error_codes(wrong) == twice && where_wrong == both,
                  "a structure_ref naming a citations section is E_REF_UNRESOLVED at the "
                  "spec member and the section field");
            std::vector<Finding> missing = job_archive("nowhere", Later::Nothing);
            CHECK(error_codes(missing) == twice,
                  "a structure_ref naming no section is E_REF_UNRESOLVED");
        }

        // C ABI.
        {
            qvf_writer* cw = qvf_create("c-abi", "1.0", "wp27");
            const char* syms[1] = {"He"};
            int Z[1] = {2};
            double pos[3] = {0, 0, 0};
            CHECK(qvf_add_structure(cw, syms, Z, pos, 1) == 0, "C: structure");
            CHECK(qvf_add_job_spec_json(cw, nullptr) != 0, "C: NULL spec_json is an error");
            CHECK(qvf_add_job_spec_json(cw, "") != 0, "C: empty spec_json is an error");
            CHECK(qvf_add_job_spec_json(cw, "{\"job_type\":") != 0 &&
                      std::string(qvf_last_error(cw)).find("not valid JSON") != std::string::npos,
                  "C: malformed spec_json is an error");
            CHECK(qvf_add_job_spec_json(cw, "[1]") != 0, "C: a non-object spec_json is an error");
            CHECK(qvf_add_job_spec_json(cw, "{\"method\":\"rks\"}") != 0,
                  "C: a spec_json without job_type is an error");
            CHECK(qvf_add_job_spec_json(cw, "{\"job_type\":\"cluster\"}") != 0,
                  "C: a wrong job_type is an error");
            CHECK(qvf_add_job_spec_json(
                      cw, "{\"job_type\":\"molecular\",\"scf\":{\"max_iterations\":0}}") != 0 &&
                      std::string(qvf_last_error(cw)).find("max_iterations") != std::string::npos,
                  "C: a wrong typed field names the field");
            CHECK(qvf_add_job_spec_json(
                      cw,
                      "{\"job_type\":\"periodic\",\"method\":\"rks\","
                      "\"structure_ref\":\"structure\",\"dispersion\":\"d3bj\","
                      "\"solvation\":{\"model\":\"pcm\",\"epsilon\":78.4},"
                      "\"scf\":{\"convergence\":1e-8,\"max_iterations\":100},"
                      "\"geometry\":{\"coordinates\":\"internal\"},"
                      "\"resources\":{\"n_cores\":16,\"walltime_s\":86400}}") == 0,
                  "C: qvf_add_job_spec_json");
            CHECK(qvf_add_job_spec_json(cw, "{\"job_type\":\"molecular\"}") != 0 &&
                      std::string(qvf_last_error(cw)).find("duplicate") != std::string::npos,
                  "C: a second job_spec section is a duplicate id, an error");
            CHECK(qvf_set_provenance_json(
                      cw,
                      "{\"method\":\"rks\",\"software\":[{\"name\":\"c-abi\",\"version\":\"1.0\"}],"
                      "\"git_commit\":\"9f2b6c1d4e7a\",\"n_threads\":8,\"n_mpi_ranks\":4,"
                      "\"n_gpus\":0,\"started_utc\":\"2026-10-07T12:00:00Z\","
                      "\"finished_utc\":\"2026-10-07T12:07:31Z\"}") == 0,
                  "C: provenance with the typed fields");
            CHECK(qvf_set_provenance_json(cw, "{\"git_commit\":\"xyz\"}") != 0 &&
                      std::string(qvf_last_error(cw)).find("git_commit") != std::string::npos,
                  "C: a wrong git_commit is an error");
            CHECK(qvf_set_provenance_json(cw, "{\"started_utc\":\"2026-10-07T12:00:00+02:00\"}") != 0,
                  "C: a started_utc with an offset is an error");
            CHECK(qvf_set_provenance_json(cw, "{\"n_threads\":0}") != 0,
                  "C: n_threads 0 is an error");
            CHECK(qvf_set_provenance_json(cw, nullptr) != 0,
                  "C: NULL provenance_json is an error");
            const char* path = "c_abi_job_spec_full.qvf";
            CHECK(qvf_write(cw, path) == 0, "C: write");
            qvf_destroy(cw);
            CHECK(!has_errors(QvfReader::validate_file(path)),
                  "C: the archive with the typed fields validates");
            QvfReader cr = QvfReader::open(path);
            CHECK(cr.manifest().at("provenance").at("n_mpi_ranks").as_int() == 4 &&
                      cr.read_json_member("job_spec", "spec").at("geometry").at("coordinates")
                              .as_string() == "internal",
                  "C: the typed fields round-trip");
            std::remove(path);
        }

    }

    // --- trajectory extensions: cell, velocities, gradient, times, energies,
    //     float32 MD coords (spec § 4.3, § 5.5) ---
    {
        auto throws = [](auto fn) {
            try {
                fn();
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        const int64_t nf = 3, na = 2;
        std::vector<double> xyz(static_cast<size_t>(nf * na * 3));
        for (size_t k = 0; k < xyz.size(); ++k) xyz[k] = 0.25 * static_cast<double>(k);
        const Tensor coords64 = tensor_f64(xyz, {nf, na, 3});
        const Tensor coords32 = tensor_f32(xyz, {nf, na, 3});
        const Tensor vel = tensor_f64(std::vector<double>(18, 0.001), {nf, na, 3});
        const Tensor grad = tensor_f64(std::vector<double>(18, -0.01), {nf, na, 3});
        const Tensor times = tensor_f64({0.0, 0.5, 1.0}, {nf});
        const Tensor energies = tensor_f64({-1.0, -1.1, -1.2}, {nf});
        const Tensor volumes = tensor_f64({216.0, 217.0, 218.0}, {nf});
        const Tensor temps = tensor_f64({300.0, 301.0, 299.0}, {nf});
        // A non-symmetric cell, so a transposed write would show.
        const Lattice::Mat3 rows = {{{6.0, 0.0, 0.0}, {0.5, 6.0, 0.0}, {0.25, 0.75, 6.0}}};
        std::vector<Lattice> cells;
        for (int k = 0; k < 3; ++k) {
            Lattice::Mat3 r = rows;
            r[0][0] += 0.1 * k;
            cells.push_back(Lattice::from_rows_angstrom(r));
        }
        auto meta_with = [](const std::string& extra) {
            return Json::parse(R"({"symbols":["H","H"],"kind":"md")" + extra + "}");
        };
        const Json md_meta =
            meta_with(R"(,"ensemble":"npt","timestep_fs":0.5,"pbc":[true,true,true],)"
                      R"("energies":[-1.0,-1.1,-1.2])");

        QvfWriter tw({"traj-test", "1.0", "md"});
        TrajectoryExtras extras;
        extras.velocities = &vel;
        extras.gradient = &grad;
        extras.times = &times;
        extras.energies = &energies;
        extras.cell_volumes = &volumes;
        extras.temperatures = &temps;
        tw.add_trajectory(md_meta, coords32, cells, extras, "md");
        tw.add_trajectory(meta_with(R"(,"pbc":[true,true,false])"), coords64,
                          Lattice::from_rows_angstrom(rows), TrajectoryExtras(), "fixed");
        tw.add_trajectory(meta_with(""), coords64, "plain");
        const Tensor cols = tensor_f64({1, 0, 0, 0, 1, 0, 0, 0, 1}, {3, 3});
        TrajectoryExtras shared;
        shared.velocities = &vel;
        shared.gradient = &grad;
        shared.times = &times;
        shared.energies = &energies;
        tw.add_reaction_path(meta_with(R"(,"energies":[-1.0,-1.1,-1.2])"), coords64, &cols,
                             shared, "path");
        tw.add_reaction_path(meta_with(""), coords64, nullptr, shared, "path_molecular");

        // Refusals, each before anything is written.
        CHECK(throws([&] {
                  TrajectoryExtras e;
                  const Tensor lat = tensor_f64(std::vector<double>(9, 1.0), {3, 3});
                  e.lattice = &lat;
                  tw.add_trajectory(meta_with(""), coords64, e, "b1");
              }),
              "a trajectory lattice needs metadata pbc");
        CHECK(throws([&] {
                  tw.add_trajectory(meta_with(R"(,"pbc":[false,true,false])"), coords64, "b2");
              }),
              "a periodic pbc needs a lattice");
        CHECK(throws([&] { tw.add_trajectory(meta_with(""), coords32, cells, extras, "b3"); }),
              "a lattice without pbc is refused for Lattice lists too");
        CHECK(throws([&] {
                  tw.add_trajectory(md_meta, coords32,
                                    std::vector<Lattice>(cells.begin(), cells.begin() + 2),
                                    TrajectoryExtras(), "b4");
              }),
              "a per-frame lattice must have n_frames cells");
        CHECK(throws([&] {
                  TrajectoryExtras e;
                  const Tensor bad = tensor_f64(std::vector<double>(24, 0.0), {4, na, 3});
                  e.velocities = &bad;
                  tw.add_trajectory(meta_with(""), coords64, e, "b5");
              }),
              "velocities with the wrong frame count");
        CHECK(throws([&] {
                  TrajectoryExtras e;
                  const Tensor bad = tensor_f32({0.0, 0.5, 1.0}, {nf});
                  e.times = &bad;
                  tw.add_trajectory(meta_with(""), coords64, e, "b6");
              }),
              "float32 times are refused");
        CHECK(throws([&] {
                  TrajectoryExtras e;
                  const Tensor off = tensor_f64({-1.0, -1.1, -1.2 + 1e-9}, {nf});
                  e.energies = &off;
                  tw.add_trajectory(meta_with(R"(,"energies":[-1.0,-1.1,-1.2])"), coords64,
                                    e, "b7");
              }),
              "binary energies must match metadata energies within 1e-10");
        CHECK(throws([&] {
                  tw.add_trajectory(meta_with(R"(,"precision":"float64")"), coords32, "b8");
              }),
              "precision float64 contradicts float32 coords");
        CHECK(throws([&] {
                  tw.add_trajectory(meta_with(R"(,"precision":"float32")"), coords64, "b9");
              }),
              "precision float32 contradicts float64 coords");
        CHECK(throws([&] {
                  tw.add_trajectory(meta_with(R"(,"units":{"lattice":"bohr"})"), coords64,
                                    "b10");
              }),
              "a trajectory lattice is declared in angstrom only");
        CHECK(throws([&] {
                  TrajectoryExtras e;
                  e.temperatures = &temps;
                  tw.add_reaction_path(meta_with(""), coords64, nullptr, e, "b11");
              }),
              "a reaction path has no temperatures");
        CHECK(throws([&] { tw.add_reaction_path(meta_with(""), coords32, nullptr, shared, "b12"); }),
              "reaction path coords stay float64");
        CHECK(throws([&] {
                  TrajectoryExtras e;
                  e.lattice = &cols;
                  tw.add_trajectory(md_meta, coords64, Lattice::from_rows_angstrom(rows), e,
                                    "b13");
              }),
              "a lattice passed twice is refused");
        CHECK(throws([&] {
                  tw.add_trajectory(Json::raw("{not json"), coords64, "b14");
              }),
              "raw metadata that does not parse is refused");

        std::vector<Finding> fs = QvfReader::validate_bytes(tw.to_bytes());
        CHECK(!has_errors(fs), "trajectory extensions archive validates");
        QvfReader r = QvfReader::from_bytes(tw.to_bytes());
        const Json md = r.section("md").at("members");
        CHECK(md.at("coords").at("dtype").as_string() == "float32",
              "float32 coords are written as float32");
        CHECK(r.read_json_member("md", "metadata").at("precision").as_string() == "float32",
              "float32 coords add metadata precision");
        CHECK(r.read_binary_member("md", "lattice").shape == std::vector<int64_t>({3, 3, 3}),
              "one cell per frame");
        CHECK(r.read_binary_member("md", "lattice").as_doubles()[9 * 2 + 0] == 6.2 &&
                  r.read_binary_member("md", "lattice").as_doubles()[3] == 0.5,
              "trajectory lattice is written as rows in angstrom");
        CHECK(r.read_binary_member("fixed", "lattice").as_doubles()[7] == 0.75,
              "fixed cell rows in angstrom (row 2, column 1)");
        for (const char* role : {"velocities", "gradient", "times", "energies",
                                 "cell_volumes", "temperatures"})
            CHECK(md.contains(role), std::string("md section carries ") + role);
        CHECK(!r.section("plain").at("members").contains("lattice"),
              "a plain trajectory has no lattice");
        CHECK(r.read_binary_member("path", "energies").as_doubles()[2] == -1.2 &&
                  r.section("path").at("members").contains("lattice"),
              "reaction path shares the per-frame members");
        for (const char* id : {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10",
                               "b11", "b12", "b13", "b14"})
            CHECK(!r.has_section(id), std::string("refused trajectory left nothing: ") + id);
    }

    // --- C ABI: qvf_add_trajectory (spec § 5.5) ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "md");
        double xyz[12] = {0, 0, 0, 0, 0, 0.74, 0, 0, 0.01, 0, 0, 0.75};   // 2 frames, 2 atoms
        double rows[18] = {5, 0, 0, 0.5, 5, 0, 0, 0, 5, 5.1, 0, 0, 0.5, 5, 0, 0, 0, 5};
        double vel[12] = {0};
        double grad[12] = {0, 0, -0.01, 0, 0, 0.01, 0, 0, -0.02, 0, 0, 0.02};
        double times[2] = {0.0, 1.0};
        double energies[2] = {-1.10, -1.11};
        const char* meta =
            "{\"symbols\":[\"H\",\"H\"],\"kind\":\"md\",\"pbc\":[true,true,true],"
            "\"precision\":\"float32\",\"energies\":[-1.10,-1.11],\"timestep_fs\":1.0}";
        CHECK(qvf_add_trajectory(cw, "{\"symbols\":[\"H\",\"H\"]}", 2, 2, xyz, rows, 0, nullptr,
                                 nullptr, nullptr, nullptr) != 0 &&
                  std::string(qvf_last_error(cw)).find("pbc") != std::string::npos,
              "C: a lattice without pbc is an error");
        CHECK(qvf_add_trajectory(cw, nullptr, 2, 2, xyz, nullptr, 0, nullptr, nullptr, nullptr,
                                 nullptr) != 0,
              "C: NULL metadata is an error");
        CHECK(qvf_add_trajectory(cw, "{\"symbols\":[\"H\",\"H\"]}", 0, 2, xyz, nullptr, 0,
                                 nullptr, nullptr, nullptr, nullptr) != 0,
              "C: zero frames is an error");
        CHECK(qvf_add_trajectory(cw, "{not json", 2, 2, xyz, nullptr, 0, nullptr, nullptr,
                                 nullptr, nullptr) != 0,
              "C: malformed metadata is an error");
        CHECK(qvf_add_trajectory(cw, meta, 2, 2, xyz, rows, 1, vel, grad, times, energies) == 0,
              "C: add_trajectory with every member");
        const char* path = "c_abi_trajectory.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write");
        qvf_destroy(cw);
        std::vector<Finding> cfs = QvfReader::validate_file(path);
        CHECK(!has_errors(cfs), "C: trajectory archive validates");
        QvfReader cr = QvfReader::open(path);
        CHECK(cr.section("trajectory").at("members").at("coords").at("dtype").as_string() ==
                  "float32",
              "C: precision float32 stores float32 coords");
        CHECK(cr.read_binary_member("trajectory", "lattice").shape ==
                      std::vector<int64_t>({2, 3, 3}) &&
                  cr.read_binary_member("trajectory", "lattice").as_doubles()[9] == 5.1 &&
                  cr.read_binary_member("trajectory", "lattice").as_doubles()[3] == 0.5,
              "C: per-frame lattice rows written as given");
        CHECK(cr.read_binary_member("trajectory", "gradient").as_doubles()[11] == 0.02,
              "C: gradient written");
        CHECK(cr.read_binary_member("trajectory", "energies").as_doubles()[1] == -1.11,
              "C: binary energies written");
        std::remove(path);
    }

    // --- wavefunction.gto extensions (spec § 5.3, Appendix A.1) ---
    {
        auto throws = [](auto fn) {
            try {
                fn();
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        auto eye = [](int64_t n) {
            Matrix m;
            m.rows = m.cols = n;
            m.data.assign(static_cast<size_t>(n * n), 0.0);
            for (int64_t i = 0; i < n; ++i) m.data[static_cast<size_t>(i * n + i)] = 1.0;
            return m;
        };
        // s + spherical p on the oxygen, contraction-normalized: n_ao 4, the
        // overlap is the identity (one centre, unit-normalized single primitives).
        auto base = [&]() {
            WavefunctionGTO wf;
            wf.shells = {{0, 0, {1.0}, {1.0}, true}, {0, 1, {0.5}, {1.0}, true}};
            wf.mo_coefficients = eye(4);
            wf.occupations = {2, 0, 0, 0};
            return wf;
        };
        QvfWriter ww({"wf-test", "1.0", "wavefunction"});
        ww.add_structure({{"O", 8, {{0.0, 0.0, 0.117}}}, {"H", 1, {{0.0, 0.757, -0.469}}},
                          {"H", 1, {{0.0, -0.757, -0.469}}}});
        WavefunctionGTO wf = base();
        wf.ecp = Json::raw(
            R"([{"center":0,"n_core_electrons":2,"label":"toy","projectors":[)"
            R"({"l":1,"exponents":[2.0],"coefficients":[-1.0],"powers":[2],"local":true},)"
            R"({"l":0,"exponents":[3.0,1.0],"coefficients":[5.0,-0.5],"powers":[0,2]}]}])");
        wf.n_electrons = 2;
        wf.symmetry_labels = {"a1", "b2", "a1", "b1"};
        wf.overlap = eye(4);
        wf.density_matrix = eye(4);
        CHECK(ww.add_wavefunction_gto(wf) == "wf", "wavefunction with every option");

        WavefunctionGTO gw = base();     // generalized: 2 n_ao complex columns
        gw.section_id = "wf_ghf";
        gw.spin = "generalized";
        gw.complex_coefficients = true;
        gw.mo_coefficients.rows = 8;
        gw.mo_coefficients.cols = 8;
        gw.mo_coefficients.data.assign(128, 0.0);
        for (int i = 0; i < 8; ++i) gw.mo_coefficients.data[static_cast<size_t>(16 * i + 2 * i)] = 1.0;
        gw.occupations = {1, 1, 0, 0, 0, 0, 0, 0};
        CHECK(ww.add_wavefunction_gto(gw) == "wf_ghf", "generalized complex wavefunction");

        WavefunctionGTO uw = base();     // unrestricted with labels and densities
        uw.section_id = "wf_uhf";
        uw.unrestricted = true;
        uw.mo_coefficients = Matrix();
        uw.mo_coefficients_alpha = eye(4);
        uw.mo_coefficients_beta = eye(4);
        uw.occupations_alpha = {1, 0, 0, 0};
        uw.occupations_beta = {1, 0, 0, 0};
        uw.symmetry_labels_beta = {"a", "b", "c", "d"};
        uw.density_matrix_alpha = eye(4);
        uw.density_matrix_beta = eye(4);
        uw.n_electrons = 2;
        CHECK(ww.add_wavefunction_gto(uw) == "wf_uhf", "unrestricted wavefunction");

        // Refusals leave nothing behind.
        auto refused = [&](const char* label, const std::function<void(WavefunctionGTO&)>& edit) {
            WavefunctionGTO bad = base();
            bad.section_id = "refused";
            edit(bad);
            CHECK(throws([&] { ww.add_wavefunction_gto(bad); }), label);
        };
        refused("overlap of the wrong shape", [&](WavefunctionGTO& b) { b.overlap = eye(3); });
        refused("density_matrix with unrestricted", [&](WavefunctionGTO& b) {
            b.unrestricted = true;
            b.mo_coefficients_alpha = b.mo_coefficients_beta = eye(4);
            b.density_matrix = eye(4);
        });
        refused("density_matrix_alpha without beta", [&](WavefunctionGTO& b) {
            b.unrestricted = true;
            b.mo_coefficients_alpha = b.mo_coefficients_beta = eye(4);
            b.density_matrix_alpha = eye(4);
        });
        refused("density_matrix_alpha with a single set", [&](WavefunctionGTO& b) {
            b.density_matrix_alpha = b.density_matrix_beta = eye(4);
        });
        refused("n_electrons against the occupations",
                [&](WavefunctionGTO& b) { b.n_electrons = 3; });
        refused("generalized without complex", [&](WavefunctionGTO& b) {
            b.spin = "generalized";
            b.mo_coefficients = eye(8);
        });
        refused("unknown spin", [&](WavefunctionGTO& b) { b.spin = "ghf"; });
        refused("unrestricted flag with a single-set spin", [&](WavefunctionGTO& b) {
            b.unrestricted = true;
            b.spin = "restricted_open";
        });
        refused("complex flag with real-sized data",
                [&](WavefunctionGTO& b) { b.complex_coefficients = true; });
        refused("AO columns", [&](WavefunctionGTO& b) { b.mo_coefficients = eye(5); });
        refused("symmetry_labels length",
                [&](WavefunctionGTO& b) { b.symmetry_labels = {"a", "b"}; });
        refused("ecp not an array", [&](WavefunctionGTO& b) { b.ecp = Json::object(); });
        refused("ecp invalid JSON", [&](WavefunctionGTO& b) { b.ecp = Json::raw("[{"); });
        refused("ecp duplicate centre", [&](WavefunctionGTO& b) {
            b.ecp = Json::raw(R"([{"center":0,"n_core_electrons":2},)"
                              R"({"center":0,"n_core_electrons":2}])");
        });
        refused("ecp ragged projector", [&](WavefunctionGTO& b) {
            b.ecp = Json::raw(R"([{"center":0,"n_core_electrons":2,"projectors":[)"
                              R"({"l":0,"exponents":[1.0,2.0],"coefficients":[1.0],)"
                              R"("powers":[2,2]}]}])");
        });
        refused("ecp local not the highest l", [&](WavefunctionGTO& b) {
            b.ecp = Json::raw(R"([{"center":0,"n_core_electrons":2,"projectors":[)"
                              R"({"l":1,"exponents":[1.0],"coefficients":[1.0],"powers":[2],)"
                              R"("local":true},{"l":2,"exponents":[1.0],)"
                              R"("coefficients":[1.0],"powers":[2]}]}])");
        });
        refused("ecp unknown key", [&](WavefunctionGTO& b) {
            b.ecp = Json::raw(R"([{"center":0,"n_core_electrons":2,"core":2}])");
        });

        // set_wavefunction_overlap / set_wavefunction_ecp on an added section.
        WavefunctionGTO later = base();
        later.section_id = "wf_later";
        ww.add_wavefunction_gto(later);
        CHECK(throws([&] { ww.set_wavefunction_overlap("wf_later", eye(5)); }),
              "set_wavefunction_overlap: wrong n_ao throws");
        CHECK(throws([&] { ww.set_wavefunction_overlap("structure", eye(4)); }),
              "set_wavefunction_overlap: another kind throws");
        CHECK(ww.set_wavefunction_overlap("wf_later", eye(4)) == "wf_later",
              "set_wavefunction_overlap");
        CHECK(throws([&] { ww.set_wavefunction_overlap("wf_later", eye(4)); }),
              "set_wavefunction_overlap: a second overlap throws");
        CHECK(throws([&] { ww.set_wavefunction_ecp("wf_later", Json::raw("{}")); }),
              "set_wavefunction_ecp: a non-array throws");
        CHECK(ww.set_wavefunction_ecp("wf_later",
                                      Json::raw(R"([{"center":1,"n_core_electrons":0}])")) ==
                  "wf_later",
              "set_wavefunction_ecp");

        const std::vector<uint8_t> bytes = ww.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        CHECK(!has_errors(fs), "wavefunction extensions archive validates");
        for (const auto& f : fs)
            if (f.severity == Finding::Severity::Error)
                std::cerr << "  " << f.code << " " << f.where << " " << f.message << "\n";
        QvfReader wr = QvfReader::from_bytes(bytes);
        CHECK(!wr.has_section("refused"), "refused wavefunctions left nothing");
        const Json basis = wr.read_json_member("wf", "basis");
        CHECK(basis.at("n_electrons").as_int() == 2 &&
                  basis.at("ecp").at(0).at("projectors").at(0).at("local").as_bool(),
              "basis carries n_electrons and the ecp");
        CHECK(wr.read_json_member("wf", "mo_metadata").at("symmetry_labels").size() == 4,
              "symmetry_labels written");
        CHECK(wr.read_binary_member("wf", "overlap").shape == std::vector<int64_t>({4, 4}),
              "overlap written as [n_ao, n_ao]");
        CHECK(wr.read_binary_member("wf", "density_matrix").shape ==
                  std::vector<int64_t>({4, 4}),
              "density_matrix written");
        const Json gmeta = wr.read_json_member("wf_ghf", "mo_metadata");
        CHECK(gmeta.at("spin").as_string() == "generalized" && gmeta.at("complex").as_bool() &&
                  gmeta.at("n_ao").as_int() == 4,
              "generalized mo_metadata");
        CHECK(wr.read_binary_member("wf_ghf", "mo_coefficients").shape ==
                  std::vector<int64_t>({8, 8, 2}),
              "generalized coefficients are [n_mo, 2 n_ao, 2]");
        CHECK(wr.read_json_member("wf_uhf", "mo_metadata")
                      .at("beta").at("symmetry_labels").size() == 4,
              "per-channel symmetry_labels written");
        CHECK(wr.read_binary_member("wf_later", "overlap").shape ==
                  std::vector<int64_t>({4, 4}),
              "set_wavefunction_overlap wrote the member");
        CHECK(wr.read_json_member("wf_later", "basis").at("ecp").at(0).at("center").as_int() == 1,
              "set_wavefunction_ecp rewrote the basis");

        // An overlap that breaks Appendix A.1: d_xy normalized on its own.
        QvfWriter nw({"wf-test", "1.0", "normalization"});
        nw.add_structure({{"Ne", 10, {{0.0, 0.0, 0.0}}}});
        WavefunctionGTO cw;
        cw.pure = false;
        cw.shells = {{0, 2, {0.8}, {1.0}, false}};
        cw.mo_coefficients = eye(6);
        cw.overlap = eye(6);          // the diagonal should be 1, 1/3, 1/3, 1, 1/3, 1
        nw.add_wavefunction_gto(cw);
        bool flagged = false;
        for (const auto& f : QvfReader::validate_bytes(nw.to_bytes()))
            if (f.code == "E_NORMALIZATION" && f.where == "/sections/1/members/overlap")
                flagged = true;
        CHECK(flagged, "a unit d_xy self-overlap is E_NORMALIZATION");
    }

    // --- C ABI: unrestricted wavefunction, overlap, ECP ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "uhf");
        const char* syms[1] = {"Li"};
        int Z[1] = {3};
        double pos[3] = {0, 0, 0};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 1) == 0, "C: add_structure ok");
        int center[2] = {0, 0}, l[2] = {0, 1}, nprim[2] = {1, 1};
        double exps[2] = {1.0, 0.5}, coefs[2] = {1.0, 1.0};
        double c[16] = {0};
        for (int i = 0; i < 4; ++i) c[i * 4 + i] = 1.0;
        double occ_a[4] = {1, 1, 0, 0}, occ_b[4] = {1, 0, 0, 0};
        CHECK(qvf_add_wavefunction_gto_unrestricted(cw, 2, center, l, nprim, exps, coefs, 1,
                                                    0, 4, 4, c, nullptr, occ_a, 4, c,
                                                    nullptr, occ_b) == 0,
              "C: add_wavefunction_gto_unrestricted ok");
        CHECK(qvf_add_wavefunction_gto_unrestricted(cw, 2, center, l, nprim, exps, coefs, 1,
                                                    0, 4, 4, c, nullptr, occ_a, 4, nullptr,
                                                    nullptr, occ_b) != 0,
              "C: NULL beta coefficients fail");
        CHECK(qvf_set_wavefunction_overlap(cw, nullptr, 3, c) != 0,
              "C: overlap of the wrong n_ao fails");
        CHECK(qvf_set_wavefunction_overlap(cw, "wf", 4, c) == 0, "C: set_wavefunction_overlap");
        CHECK(qvf_set_wavefunction_overlap(cw, "wf", 4, c) != 0, "C: a second overlap fails");
        CHECK(qvf_set_wavefunction_ecp_json(cw, nullptr, "[{") != 0, "C: invalid ECP JSON fails");
        CHECK(std::string(qvf_last_error(cw)).find("ecp") != std::string::npos,
              "C: the ECP error is reported");
        CHECK(qvf_set_wavefunction_ecp_json(
                  cw, "wf", R"([{"center":0,"n_core_electrons":2,"label":"toy"}])") == 0,
              "C: set_wavefunction_ecp_json");
        CHECK(qvf_set_wavefunction_ecp_json(cw, "missing", "[]") != 0,
              "C: ECP on a missing section fails");
        const char* path = "c_abi_uhf.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        qvf_destroy(cw);
        try {
            QvfReader r = QvfReader::open(path);
            CHECK(!has_errors(r.validate_findings()), "C: unrestricted archive validates");
            const Json meta = r.read_json_member("wf", "mo_metadata");
            CHECK(meta.at("spin").as_string() == "unrestricted" &&
                      meta.at("alpha").at("occupations").at(1).as_double() == 1.0,
                  "C: unrestricted mo_metadata");
            CHECK(r.read_binary_member("wf", "mo_coefficients_beta").shape ==
                      std::vector<int64_t>({4, 4}),
                  "C: beta coefficients written");
            CHECK(r.read_binary_member("wf", "overlap").as_doubles()[5] == 1.0,
                  "C: overlap written");
            CHECK(r.read_json_member("wf", "basis").at("ecp").at(0).at("label").as_string() ==
                      "toy",
                  "C: ecp written");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading the unrestricted archive: ") + e.what());
        }
        std::remove(path);
    }

    // --- volumes: vector, complex, series, Laplacian; periodic grids (spec § 5.2) ---
    {
        const Lattice::Mat3 cubic = {{{{4.0, 0, 0}}, {{0, 4.0, 0}}, {{0, 0, 4.0}}}};
        const double step = kBohrPerAngstrom;          // 4 Å over 4 points
        auto periodic_writer = [&](bool endpoint) {
            QvfWriter w({"vol-test", "1.0", "volumes"});
            w.auto_metadata = false;
            w.add_structure({Atom{"Na", 11, {{0, 0, 0}}}}, {{true, true, true}},
                            Lattice::from_rows_angstrom(cubic));
            Grid pg;
            pg.voxel_vectors = {{{{step, 0, 0}}, {{0, step, 0}}, {{0, 0, step}}}};
            pg.periodic = std::array<bool, 3>{{true, true, true}};
            pg.includes_endpoint = endpoint;
            pg.structure_ref = "structure";
            w.add_volume("volume.density", pg,
                         tensor_f32(std::vector<double>(64, 0.5), {4, 4, 4}), "", "", "rho");
            return w;
        };
        QvfWriter vw = periodic_writer(false);
        const Grid g;
        vw.add_volume("volume.vector", g, tensor_f32(std::vector<double>(24, 0.1), {2, 2, 2, 3}),
                      "gradient", "", "vec");
        vw.add_volume("volume.complex", g,
                      tensor_f64(std::vector<double>(16, 0.2), {2, 2, 2, 2}), "", "", "psi");
        vw.add_volume("volume.laplacian", g,
                      tensor_f32(std::vector<double>(8, -1.0), {2, 2, 2}), "", "", "lap");
        Json tmeta = Json::object();
        tmeta.set("symbols", Json::from_strings({"Na"}));
        vw.add_trajectory(tmeta, tensor_f64(std::vector<double>(9, 0.0), {3, 1, 3}), "traj");
        Json axis = Json::object();
        axis.set("kind", "frame").set("values", Json::from_ints({0, 2}))
            .set("labels", Json::from_strings({"first", "last"})).set("trajectory_ref", "traj");
        CHECK(vw.add_volume_series(g, tensor_f32(std::vector<double>(16, 0.3), {2, 2, 2, 2}),
                                   axis, "frames", "ser") == "ser",
              "add_volume_series returns its id");
        CHECK(vw.add_volume_series(g, tensor_f64(std::vector<double>(8, 0.3), {1, 2, 2, 2}),
                                   Json::raw(R"({"kind":"time","values":[0.5],"units":"fs"})"))
                  .rfind("volume_series_", 0) == 0,
              "add_volume_series: a raw axis and the default id");

        const Tensor t4 = tensor_f32(std::vector<double>(16, 0.0), {2, 2, 2, 2});
        CHECK(throws([&] { vw.add_volume("volume.vector", g, t4, "", "", "b1"); }),
              "add_volume: volume.vector with a trailing axis of 2 throws");
        CHECK(throws([&] {
                  vw.add_volume("volume.complex", g,
                                tensor_f32(std::vector<double>(24, 0.0), {2, 2, 2, 3}), "", "",
                                "b2");
              }),
              "add_volume: volume.complex with a trailing axis of 3 throws");
        CHECK(throws([&] { vw.add_volume("volume.laplacian", g, t4, "", "", "b3"); }),
              "add_volume: a rank-4 Laplacian throws");
        CHECK(throws([&] { vw.add_volume("volume.series", g, t4, "", "", "b4"); }),
              "add_volume: volume.series needs add_volume_series");
        CHECK(throws([&] {
                  vw.add_volume("volume.density", g, tensor_i32({0, 0, 0, 0, 0, 0, 0, 0}, {2, 2, 2}),
                                "", "", "b5");
              }),
              "add_volume: integer data throws");
        Grid badref = g;
        badref.structure_ref = "traj";
        CHECK(throws([&] {
                  vw.add_volume("volume.density", badref,
                                tensor_f32(std::vector<double>(8, 0.0), {2, 2, 2}), "", "", "b6");
              }),
              "add_volume: a grid structure_ref naming a trajectory throws");
        auto series_throws = [&](const std::string& axis_json, const char* id) {
            return throws([&] { vw.add_volume_series(g, t4, Json::raw(axis_json), "", id); });
        };
        CHECK(series_throws(R"({"kind":"time","values":[0]})", "b7"),
              "add_volume_series: one value for two slabs throws");
        CHECK(series_throws(R"({"kind":"epoch","values":[0,1]})", "b8"),
              "add_volume_series: an unknown axis kind throws");
        CHECK(series_throws(R"({"kind":"frame","values":[0.5,1]})", "b9"),
              "add_volume_series: a fractional frame index throws");
        CHECK(series_throws(R"({"kind":"frame","values":[0,1],"trajectory_ref":"structure"})",
                            "b10"),
              "add_volume_series: a trajectory_ref naming a structure throws");
        CHECK(series_throws(R"({"kind":"time","values":[0,1],"units":"sec"})", "b11"),
              "add_volume_series: an unknown unit throws");
        CHECK(series_throws(R"({"kind":"time","values":[0,1],"step":1})", "b12"),
              "add_volume_series: an unknown axis key throws");
        CHECK(throws([&] {
                  vw.add_volume_series(g, tensor_f32(std::vector<double>(8, 0.0), {2, 2, 2}),
                                       axis, "", "b13");
              }),
              "add_volume_series: rank-3 data throws");

        const std::vector<uint8_t> bytes = vw.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        CHECK(!has_errors(fs), "volume archive validates without errors");
        QvfReader vr = QvfReader::from_bytes(bytes);
        const Json grid = vr.read_json_member("rho", "grid");
        CHECK(grid.at("periodic").at(2).as_bool() && !grid.at("includes_endpoint").as_bool() &&
                  grid.at("structure_ref").as_string() == "structure",
              "periodic grid fields written");
        CHECK(!vr.read_json_member("vec", "grid").has("periodic") &&
                  !vr.read_json_member("vec", "grid").has("includes_endpoint") &&
                  !vr.read_json_member("vec", "grid").has("structure_ref"),
              "unset grid fields are omitted");
        CHECK(vr.read_binary_member("vec", "data").shape == std::vector<int64_t>({2, 2, 2, 3}),
              "vector data written as [2, 2, 2, 3]");
        CHECK(vr.read_json_member("vec", "grid").at("shape").size() == 3 &&
                  vr.read_json_member("ser", "grid").at("shape").at(0).as_int() == 2,
              "grid shape is the spatial extents");
        CHECK(vr.read_json_member("ser", "axis").at("trajectory_ref").as_string() == "traj",
              "series axis read back");
        for (const char* id : {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10",
                               "b11", "b12", "b13"})
            CHECK(!vr.has_section(id), std::string("refused volume left nothing: ") + id);

        // The same CHGCAR-style numbers declared cube style fail validation.
        QvfWriter ew = periodic_writer(true);
        int flagged = 0;
        for (const auto& f : QvfReader::validate_bytes(ew.to_bytes()))
            if (f.code == "E_MEMBER_MALFORMED" && f.where == "/sections/1/members/grid")
                ++flagged;
        CHECK(flagged == 3, "includes_endpoint true with CHGCAR numbers: one finding per axis");
    }

    // --- C ABI: volumes ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "volumes");
        const char* syms[1] = {"Na"};
        int Z[1] = {11};
        double pos[3] = {0, 0, 0};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 1) == 0, "C: structure");
        const char* grid =
            R"({"origin":[0,0,0],"voxel_vectors":[[1,0,0],[0,1,0],[0,0,1]],)"
            R"("units":{"origin":"bohr"}})";
        float vec[24];
        for (int k = 0; k < 24; ++k) vec[k] = 0.25f * static_cast<float>(k);
        const int64_t vshape[4] = {2, 2, 2, 3};
        CHECK(qvf_add_volume(cw, "volume.vector", grid, "float32", 4, vshape, vec, "gradient") == 0,
              "C: add_volume volume.vector");
        double psi[16];
        for (int k = 0; k < 16; ++k) psi[k] = -0.5 * k;
        const int64_t cshape[4] = {2, 2, 2, 2};
        CHECK(qvf_add_volume(cw, "volume.complex", grid, "float64", 4, cshape, psi, nullptr) == 0,
              "C: add_volume volume.complex");
        float lap[8] = {-1, -2, -3, -4, 1, 2, 3, 4};
        const int64_t lshape[3] = {2, 2, 2};
        CHECK(qvf_add_volume(cw, "volume.laplacian",
                             R"({"origin":[0,0,0],"voxel_vectors":[[1,0,0],[0,1,0],[0,0,1]],)"
                             R"("shape":[2,2,2],"periodic":[false,false,false],)"
                             R"("includes_endpoint":true,"structure_ref":"structure"})",
                             nullptr, 3, lshape, lap, nullptr) == 0,
              "C: add_volume with every grid key and the default dtype");
        CHECK(qvf_add_volume(cw, "volume.series", grid, "float32", 4, cshape, vec, nullptr) != 0 &&
                  std::string(qvf_last_error(cw)).find("qvf_add_volume_series") !=
                      std::string::npos,
              "C: qvf_add_volume refuses volume.series");
        CHECK(qvf_add_volume(cw, "volume.density", grid, "int32", 3, lshape, lap, nullptr) != 0,
              "C: an integer dtype is an error");
        CHECK(qvf_add_volume(cw, "volume.vector", grid, "float32", 4, cshape, vec, nullptr) != 0,
              "C: a vector with a trailing axis of 2 is an error");
        CHECK(qvf_add_volume(cw, "volume.density",
                             R"({"origin":[0,0,0],"voxel_vectors":[[1,0,0],[0,1,0],[0,0,1]],)"
                             R"("shape":[2,2,3]})",
                             "float32", 3, lshape, lap, nullptr) != 0,
              "C: a grid shape that differs from the data is an error");
        CHECK(qvf_add_volume(cw, "volume.density",
                             R"({"origin":[0,0,0],"voxel_vectors":[[1,0,0],[0,1,0],[0,0,1]],)"
                             R"("x_note":"hi"})",
                             "float32", 3, lshape, lap, nullptr) != 0,
              "C: an unknown grid key is an error");
        CHECK(qvf_add_volume(cw, "volume.density", R"({"origin":[0,0,0]})", "float32", 3,
                             lshape, lap, nullptr) != 0,
              "C: a grid without voxel_vectors is an error");
        CHECK(qvf_add_volume(cw, "volume.density", nullptr, "float32", 3, lshape, lap, nullptr) != 0,
              "C: a NULL grid is an error");
        float sdata[16] = {0};
        CHECK(qvf_add_volume_series(cw, grid, R"({"kind":"time","values":[0,0.5],"units":"fs"})",
                                    "float32", cshape, sdata, "md") == 0,
              "C: add_volume_series");
        CHECK(qvf_add_volume_series(cw, grid, R"({"kind":"time","values":[0]})", "float32",
                                    cshape, sdata, nullptr) != 0,
              "C: a series axis with too few values is an error");
        CHECK(qvf_add_volume_series(cw, grid, nullptr, "float32", cshape, sdata, nullptr) != 0,
              "C: a NULL axis is an error");
        const char* path = "c_abi_volumes.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write");
        qvf_destroy(cw);
        CHECK(!has_errors(QvfReader::validate_file(path)), "C: volume archive validates");
        try {
            QvfReader cr = QvfReader::open(path);
            const Json& secs = cr.sections();
            CHECK(secs.size() == 5, "C: four volume sections written, refusals left nothing");
            const std::string vec_id = secs.at(1).at("id").as_string();
            const BinaryMember vm = cr.read_binary_member(vec_id, "data");
            CHECK(vm.shape == std::vector<int64_t>({2, 2, 2, 3}) && vm.dtype == "float32",
                  "C: vector data [2, 2, 2, 3] float32");
            const std::string psi_id = secs.at(2).at("id").as_string();
            CHECK(cr.read_binary_member(psi_id, "data").as_doubles()[3] == -1.5,
                  "C: complex data written from doubles");
            const std::string lap_id = secs.at(3).at("id").as_string();
            const Json lg = cr.read_json_member(lap_id, "grid");
            CHECK(lg.at("includes_endpoint").as_bool() &&
                      lg.at("structure_ref").as_string() == "structure" &&
                      cr.read_binary_member(lap_id, "data").dtype == "float32",
                  "C: grid keys written, default dtype float32");
            CHECK(cr.read_json_member(secs.at(4).at("id").as_string(), "axis")
                          .at("units").as_string() == "fs",
                  "C: series axis written");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading volume archive: ") + e.what());
        }
        std::remove(path);
    }

    // --- property kinds (spec § 5.13) ---
    {
        auto throws = [](auto fn) {
            try {
                fn();
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        // rows x cols matrix with `diag` on the diagonal, as JSON text.
        auto mat = [](int rows, int cols, double diag) {
            std::string s = "[";
            for (int i = 0; i < rows; ++i) {
                s += i ? ",[" : "[";
                for (int j = 0; j < cols; ++j)
                    s += (j ? "," : "") + std::to_string(i == j ? diag : 0.0);
                s += "]";
            }
            return s + "]";
        };
        auto zeros = [](int n) {
            std::string s = "[";
            for (int i = 0; i < n; ++i) s += (i ? ",0" : "0");
            return s + "]";
        };
        const std::string electric =
            "{\"origin\":\"center_of_mass\","
            "\"dipole\":{\"vector\":[0,0,-0.75],\"total\":0.75,\"units\":\"debye\"},"
            "\"quadrupole\":{\"matrix\":" + mat(3, 3, 1.0) +
            ",\"traceless\":false,\"convention\":\"cartesian\"},"
            "\"octupole\":{\"tensor\":" + zeros(27) + "},"
            "\"polarizability\":{\"matrix\":" + mat(3, 3, 9.0) +
            ",\"frequency_ev\":0,\"isotropic\":9},"
            "\"hyperpolarizability\":{\"tensor\":" + zeros(27) +
            ",\"frequencies_ev\":[0,0]},\"units\":{\"dipole\":\"debye\"},"
            "\"x_note\":\"kept\"}";
        const std::string elastic =
            "{\"stiffness\":" + mat(6, 6, 36.0) + ",\"compliance\":" + mat(6, 6, 0.0277) +
            ",\"bulk_modulus_gpa\":18.0,\"poisson_ratio\":0.3,\"average\":\"hill\","
            "\"units\":{\"stiffness\":\"gpa\",\"poisson_ratio\":\"dimensionless\"}}";
        const std::string dielectric =
            "{\"electronic\":" + mat(3, 3, 2.0) + ",\"ionic\":" + mat(3, 3, 1.0) +
            ",\"static\":" + mat(3, 3, 3.0) + ",\"born_charges\":[" + mat(3, 3, -0.8) +
            "," + mat(3, 3, 0.4) + "," + mat(3, 3, 0.4) + "],\"piezoelectric\":" +
            mat(3, 6, 0.0) + ",\"units\":{\"born_charges\":\"e\"}}";

        QvfWriter pw({"props-test", "1.0", "properties"});
        pw.add_structure({{"O", 8, {{0.0, 0.0, 0.117}}},
                          {"H", 1, {{0.0, 0.757, -0.469}}},
                          {"H", 1, {{0.0, -0.757, -0.469}}}});
        Json root = Json::object();
        root.set("total_debye", 0.75).set("vector_debye", Json::from_doubles({0, 0, -0.75}));
        pw.set_dipole_moment(root);
        CHECK(pw.add_electric_properties(Json::parse(electric)) == "electric",
              "add_electric_properties default id");
        CHECK(pw.add_elastic_properties(Json::raw(elastic)) == "elastic",
              "add_elastic_properties default id (a raw fragment)");
        CHECK(pw.add_dielectric_properties(Json::parse(dielectric), "permittivity") ==
                  "permittivity",
              "add_dielectric_properties with an id");

        // Every refusal leaves nothing behind. Each bad payload is one valid
        // payload with one defect.
        const std::string good_dipole = "{\"dipole\":{\"vector\":[0,0,1]}}";
        const std::vector<std::pair<const char*, std::string>> bad_electric = {
            {"empty", "{}"},
            {"origin only", "{\"origin\":\"origin\"}"},
            {"not an object", "[1]"},
            {"unknown key", "{\"dipole\":{\"vector\":[0,0,1]},\"dipol\":1}"},
            {"unknown nested key", "{\"dipole\":{\"vector\":[0,0,1],\"magnitude\":1}}"},
            {"dipole without vector", "{\"dipole\":{\"total\":1}}"},
            {"short vector", "{\"dipole\":{\"vector\":[0,1]}}"},
            {"negative total", "{\"dipole\":{\"vector\":[0,0,1],\"total\":-1}}"},
            {"dipole units", "{\"dipole\":{\"vector\":[0,0,1],\"units\":\"e\"}}"},
            {"origin of two", "{\"dipole\":{\"vector\":[0,0,1]},\"origin\":[0,0]}"},
            {"empty origin", "{\"dipole\":{\"vector\":[0,0,1]},\"origin\":\"\"}"},
            {"quadrupole without convention",
             "{\"quadrupole\":{\"matrix\":" + mat(3, 3, 1.0) + ",\"traceless\":false}}"},
            {"quadrupole convention",
             "{\"quadrupole\":{\"matrix\":" + mat(3, 3, 1.0) +
                 ",\"traceless\":false,\"convention\":\"spherical\"}}"},
            {"buckingham with a trace",
             "{\"quadrupole\":{\"matrix\":" + mat(3, 3, 1.0) +
                 ",\"traceless\":false,\"convention\":\"buckingham\"}}"},
            {"asymmetric quadrupole",
             "{\"quadrupole\":{\"matrix\":[[0,1,0],[2,0,0],[0,0,0]],\"traceless\":true,"
             "\"convention\":\"cartesian\"}}"},
            {"quadrupole of 2 x 2",
             "{\"quadrupole\":{\"matrix\":" + mat(2, 2, 1.0) +
                 ",\"traceless\":false,\"convention\":\"cartesian\"}}"},
            {"octupole of 26", "{\"octupole\":{\"tensor\":" + zeros(26) + "}}"},
            {"negative frequency",
             "{\"polarizability\":{\"matrix\":" + mat(3, 3, 1.0) + ",\"frequency_ev\":-1}}"},
            {"one hyperpolarizability frequency",
             "{\"hyperpolarizability\":{\"tensor\":" + zeros(27) +
                 ",\"frequencies_ev\":[1]}}"},
            {"units with a non-unit", "{\"dipole\":{\"vector\":[0,0,1]},\"units\":{\"dipole\":\"D\"}}"},
            {"units contradicting a fixed unit",
             "{\"dipole\":{\"vector\":[0,0,1]},\"units\":{\"dipole\":\"e\"}}"},
        };
        int n = 0;
        for (const auto& c : bad_electric)
            CHECK(throws([&] {
                      pw.add_electric_properties(Json::raw(c.second),
                                                 "bad_e" + std::to_string(n++));
                  }),
                  std::string("add_electric_properties refuses: ") + c.first);
        const std::vector<std::pair<const char*, std::string>> bad_elastic = {
            {"no stiffness", "{\"average\":\"voigt\"}"},
            {"5 x 6", "{\"stiffness\":" + mat(5, 6, 1.0) + "}"},
            {"asymmetric",
             "{\"stiffness\":[[0,1,0,0,0,0],[2,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],"
             "[0,0,0,0,0,0],[0,0,0,0,0,0]]}"},
            {"asymmetric compliance",
             "{\"stiffness\":" + mat(6, 6, 1.0) +
                 ",\"compliance\":[[0,1,0,0,0,0],[2,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],"
                 "[0,0,0,0,0,0],[0,0,0,0,0,0]]}"},
            {"average", "{\"stiffness\":" + mat(6, 6, 1.0) + ",\"average\":\"mean\"}"},
            {"modulus a string", "{\"stiffness\":" + mat(6, 6, 1.0) + ",\"bulk_modulus_gpa\":\"1\"}"},
            {"unknown key", "{\"stiffness\":" + mat(6, 6, 1.0) + ",\"bulk\":1}"},
            {"units contradicting a fixed unit",
             "{\"stiffness\":" + mat(6, 6, 1.0) + ",\"units\":{\"stiffness\":\"ev\"}}"},
        };
        for (const auto& c : bad_elastic)
            CHECK(throws([&] {
                      pw.add_elastic_properties(Json::raw(c.second),
                                                "bad_c" + std::to_string(n++));
                  }),
                  std::string("add_elastic_properties refuses: ") + c.first);
        const std::vector<std::pair<const char*, std::string>> bad_dielectric = {
            {"no electronic", "{\"ionic\":" + mat(3, 3, 1.0) + "}"},
            {"2 x 3", "{\"electronic\":" + mat(2, 3, 1.0) + "}"},
            {"static is not the sum",
             "{\"electronic\":" + mat(3, 3, 2.0) + ",\"ionic\":" + mat(3, 3, 1.0) +
                 ",\"static\":" + mat(3, 3, 3.5) + "}"},
            {"no Born tensors", "{\"electronic\":" + mat(3, 3, 2.0) + ",\"born_charges\":[]}"},
            {"2 x 3 Born tensor",
             "{\"electronic\":" + mat(3, 3, 2.0) + ",\"born_charges\":[" + mat(2, 3, 1.0) + "]}"},
            {"piezoelectric 3 x 3",
             "{\"electronic\":" + mat(3, 3, 2.0) + ",\"piezoelectric\":" + mat(3, 3, 0.0) + "}"},
            {"unknown key", "{\"electronic\":" + mat(3, 3, 2.0) + ",\"dielectric\":1}"},
            {"units contradicting a fixed unit",
             "{\"electronic\":" + mat(3, 3, 2.0) + ",\"units\":{\"born_charges\":\"hartree\"}}"},
        };
        for (const auto& c : bad_dielectric)
            CHECK(throws([&] {
                      pw.add_dielectric_properties(Json::raw(c.second),
                                                   "bad_d" + std::to_string(n++));
                  }),
                  std::string("add_dielectric_properties refuses: ") + c.first);
        CHECK(throws([&] { pw.add_elastic_properties(Json::parse(elastic), "elastic"); }),
              "a duplicate section id throws");
        CHECK(throws([&] { pw.add_elastic_properties(Json::raw("{")); }),
              "malformed raw JSON throws");

        const std::vector<uint8_t> bytes = pw.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        CHECK(!has_errors(fs), "properties archive validates without errors");
        QvfReader pr = QvfReader::from_bytes(bytes);
        CHECK(pr.section("electric").at("kind").as_string() == "properties.electric",
              "electric kind");
        CHECK(pr.read_json_member("electric", "properties")
                      .at("octupole").at("tensor").size() == 27,
              "octupole tensor read back");
        CHECK(pr.read_json_member("electric", "properties").at("x_note").as_string() == "kept",
              "vendor key kept");
        CHECK(pr.read_json_member("elastic", "properties").at("average").as_string() == "hill",
              "elastic average read back");
        CHECK(pr.read_json_member("permittivity", "properties").at("born_charges").size() == 3,
              "born tensors read back");
        for (int k = 0; k < n; ++k)
            for (const char* prefix : {"bad_e", "bad_c", "bad_d"})
                CHECK(!pr.has_section(prefix + std::to_string(k)),
                      "a refused section left nothing behind");

        // What only the validators judge: the dipole against the root block and
        // the Born tensors against the structure.
        QvfWriter vw({"props-test", "1.0", "inconsistent"});
        vw.add_structure({{"H", 1, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.74}}}});
        Json off = Json::object();
        off.set("total_debye", 1.85).set("vector_debye", Json::from_doubles({0, 0, -1.85}));
        vw.set_dipole_moment(off);
        vw.add_electric_properties(
            Json::parse("{\"dipole\":{\"vector\":[0,0,-1.5],\"total\":1.5}}"));
        vw.add_dielectric_properties(Json::parse(dielectric));   // 3 tensors, 2 atoms
        bool root_flagged = false, born_flagged = false;
        for (const auto& f : QvfReader::validate_bytes(vw.to_bytes())) {
            if (f.code == "E_ROOT_INCONSISTENT" && f.where == "/sections/1/members/properties")
                root_flagged = true;
            if (f.code == "E_MEMBER_MALFORMED" && f.where == "/sections/2/members/properties")
                born_flagged = true;
        }
        CHECK(root_flagged, "a dipole off the root dipole_moment is E_ROOT_INCONSISTENT");
        CHECK(born_flagged, "Born tensors off the atom count are E_MEMBER_MALFORMED");
    }

    // --- C ABI: property kinds ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "properties");
        const char* syms[1] = {"H"};
        int Z[1] = {1};
        double pos[3] = {0, 0, 0};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 1) == 0, "C: add_structure ok");
        CHECK(qvf_add_electric_properties_json(
                  cw, "{\"dipole\":{\"vector\":[0,0,0.5],\"total\":0.5,\"units\":\"debye\"}}") == 0,
              "C: add_electric_properties_json ok");
        CHECK(qvf_add_electric_properties_json(cw, "{\"dipole\":{}}") != 0,
              "C: a dipole without a vector fails");
        CHECK(qvf_add_electric_properties_json(cw, "") != 0, "C: empty electric JSON fails");
        CHECK(qvf_add_electric_properties_json(cw, nullptr) != 0, "C: NULL electric JSON fails");
        const std::string c6 = "[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],"
                               "[0,0,0,0,1,0],[0,0,0,0,0,1]]";
        CHECK(qvf_add_elastic_properties_json(cw, ("{\"stiffness\":" + c6 + "}").c_str()) == 0,
              "C: add_elastic_properties_json ok");
        CHECK(qvf_add_elastic_properties_json(cw, "{\"stiffness\":[[1]]}") != 0,
              "C: a 1 x 1 stiffness fails");
        CHECK(qvf_add_elastic_properties_json(cw, "not json") != 0, "C: invalid JSON fails");
        const std::string c3 = "[[2,0,0],[0,2,0],[0,0,2]]";
        CHECK(qvf_add_dielectric_properties_json(
                  cw, ("{\"electronic\":" + c3 + ",\"born_charges\":[" + c3 + "]}").c_str()) == 0,
              "C: add_dielectric_properties_json ok");
        CHECK(qvf_add_dielectric_properties_json(cw, "{\"ionic\":[[1]]}") != 0,
              "C: a dielectric payload without electronic fails");
        CHECK(qvf_add_dielectric_properties_json(nullptr, "{}") != 0, "C: NULL writer fails");
        const char* path = "c_abi_properties.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: write ok");
        qvf_destroy(cw);
        try {
            QvfReader r = QvfReader::open(path);
            std::vector<std::string> errs;
            CHECK(r.validate(errs), "C: properties archive validates");
            CHECK(r.read_json_member("electric", "properties").at("dipole").at("total")
                          .as_double() == 0.5,
                  "C: electric payload read back");
            CHECK(r.read_json_member("elastic", "properties").at("stiffness").size() == 6,
                  "C: stiffness read back");
            CHECK(r.read_json_member("dielectric", "properties").at("born_charges").size() == 1,
                  "C: born_charges read back");
        } catch (const std::exception& e) {
            CHECK(false, std::string("C: reading properties archive: ") + e.what());
        }
        std::remove(path);
    }

    // --- wavefunction.bloch (spec § 5.3.1, provisional) ---
    {
        auto throws = [](auto fn) {
            try {
                fn();
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        QvfWriter bw({"bloch-test", "1.0", "bloch"});
        bw.add_structure({{"H", 1, {{0.0, 0.0, 0.0}}}, {"H", 1, {{0.8, 0.0, 0.0}}}},
                         {{true, false, false}},
                         {{{2.0, 0.0, 0.0}}, {{0.0, 10.0, 0.0}}, {{0.0, 0.0, 10.0}}});
        WavefunctionBloch wf;
        wf.shells = {Shell{0, 0, {0.4166}, {1.0}, true}, Shell{1, 0, {0.4166}, {1.0}, true}};
        wf.spin = "unrestricted";
        wf.kpoints = {{{0.0, 0.0, 0.0}}, {{0.5, 0.0, 0.0}}};
        wf.weights = {0.5, 0.5};
        wf.labels = {"G", "X"};
        wf.mesh = {2, 1, 1};
        wf.has_shift = true;
        wf.has_time_reversal = true;
        wf.n_mo = 2;
        wf.n_mo_per_k = {2, 2};
        wf.has_fermi_energy = true;
        wf.fermi_energy_ev = -3.5;
        wf.mo_coefficients.resize(2 * 2 * 2 * 2 * 2);
        for (size_t i = 0; i < wf.mo_coefficients.size(); ++i)
            wf.mo_coefficients[i] = 0.01 * static_cast<double>(i);
        wf.energies = std::vector<double>(8, -0.5);
        wf.occupations = {1, 0, 1, 0, 1, 0, 1, 0};
        CHECK(bw.add_wavefunction_bloch(wf) == "wf_bloch", "add_wavefunction_bloch default id");

        auto refused = [&](const char* id, void (*edit)(WavefunctionBloch&)) {
            WavefunctionBloch bad = wf;
            bad.section_id = id;
            edit(bad);
            return throws([&] { bw.add_wavefunction_bloch(bad); });
        };
        CHECK(refused("b1", [](WavefunctionBloch& b) { b.spin = "generalized"; }),
              "bloch: unknown spin throws");
        CHECK(refused("b2", [](WavefunctionBloch& b) { b.kpoints.clear(); }),
              "bloch: no k-points throws");
        CHECK(refused("b3", [](WavefunctionBloch& b) { b.mo_coefficients.pop_back(); }),
              "bloch: short coefficients throw");
        CHECK(refused("b4", [](WavefunctionBloch& b) { b.energies.pop_back(); }),
              "bloch: short energies throw");
        CHECK(refused("b5", [](WavefunctionBloch& b) { b.weights = {0.5, 0.25}; }),
              "bloch: weights not summing to 1 throw");
        CHECK(refused("b6", [](WavefunctionBloch& b) { b.labels = {"G"}; }),
              "bloch: one label for two k-points throws");
        CHECK(refused("b7", [](WavefunctionBloch& b) { b.mesh = {2, 0, 1}; }),
              "bloch: a mesh extent of 0 throws");
        CHECK(refused("b8", [](WavefunctionBloch& b) { b.n_mo_per_k = {2, 3}; }),
              "bloch: n_mo_per_k above n_mo throws");
        CHECK(refused("b9", [](WavefunctionBloch& b) { b.spin = "restricted"; }),
              "bloch: restricted with two spin channels' data throws");

        const std::vector<uint8_t> bytes = bw.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        CHECK(!has_errors(fs), "bloch archive validates without errors");
        QvfReader br = QvfReader::from_bytes(bytes);
        CHECK(br.read_binary_member("wf_bloch", "mo_coefficients").shape ==
                  std::vector<int64_t>({2, 2, 2, 2, 2}),
              "bloch coefficients written as [n_k, n_spin, n_mo, n_ao, 2]");
        CHECK(br.read_binary_member("wf_bloch", "mo_coefficients").as_doubles()[3] == 0.03,
              "bloch coefficients read back");
        CHECK(br.read_binary_member("wf_bloch", "occupations").shape ==
                  std::vector<int64_t>({2, 2, 2}),
              "bloch occupations written as [n_k, n_spin, n_mo]");
        const Json kp = br.read_json_member("wf_bloch", "kpoints");
        CHECK(kp.at("weights").size() == 2 && kp.at("labels").at(1).as_string() == "X" &&
                  kp.at("mesh").size() == 3 && kp.has("shift") && kp.has("time_reversal"),
              "bloch kpoints payload written");
        const Json meta = br.read_json_member("wf_bloch", "mo_metadata");
        CHECK(meta.at("n_spin").as_int() == 2 && meta.at("n_ao").as_int() == 2 &&
                  meta.at("complex").as_bool() &&
                  meta.at("k_point_convention").as_string() == "fractional" &&
                  meta.at("fermi_energy").at("units").as_string() == "ev",
              "bloch mo_metadata written");
        CHECK(br.read_json_member("wf_bloch", "basis").at("n_ao").as_int() == 2,
              "bloch basis written");
        for (const char* id : {"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9"})
            CHECK(!br.has_section(id), std::string("refused bloch section left nothing: ") + id);
    }

    // --- C ABI: wavefunction.bloch ---
    {
        qvf_writer* cw = qvf_create("c-abi", "1.0", "bloch");
        const char* syms[2] = {"H", "H"};
        int Z[2] = {1, 1};
        double pos[6] = {0, 0, 0, 0.8, 0, 0};
        CHECK(qvf_add_structure(cw, syms, Z, pos, 2) == 0, "C: add_structure ok");
        int center[2] = {0, 1}, l[2] = {0, 0}, nprim[2] = {1, 1};
        double exps[2] = {0.4166, 0.4166}, coefs[2] = {1.0, 1.0};
        double kpts[6] = {0, 0, 0, 0.5, 0, 0};
        double weights[2] = {0.5, 0.5}, bad_weights[2] = {0.5, 0.4};
        double c[16], e[4] = {-0.6, 0.4, -0.5, 0.3}, o[4] = {2, 0, 2, 0};
        for (int i = 0; i < 16; ++i) c[i] = 0.125 * i;
        CHECK(qvf_add_wavefunction_bloch(cw, 2, center, l, nprim, exps, coefs, 1, 0, 2, kpts,
                                         bad_weights, 1, 2, 2, c, e, o) != 0,
              "C: bloch weights not summing to 1 fail");
        CHECK(qvf_add_wavefunction_bloch(cw, 2, center, l, nprim, exps, coefs, 1, 0, 2, kpts,
                                         weights, 1, 2, 3, c, e, o) != 0,
              "C: bloch n_ao differing from the shells fails");
        CHECK(qvf_add_wavefunction_bloch(cw, 2, center, l, nprim, exps, coefs, 1, 0, 2, kpts,
                                         weights, 3, 2, 2, c, e, o) != 0,
              "C: bloch n_spin 3 fails");
        CHECK(qvf_add_wavefunction_bloch(cw, 2, center, l, nprim, exps, coefs, 1, 0, 2,
                                         nullptr, weights, 1, 2, 2, c, e, o) != 0,
              "C: bloch NULL k-points fail");
        CHECK(qvf_add_wavefunction_bloch(cw, 2, center, l, nprim, exps, coefs, 1, 0, 2, kpts,
                                         weights, 1, 2, 2, c, e, o) == 0,
              "C: add_wavefunction_bloch ok");
        const char* path = "c_abi_bloch.qvf";
        CHECK(qvf_write(cw, path) == 0, "C: bloch write ok");
        qvf_destroy(cw);
        try {
            // The C ABI writes molecular structures only, so the one finding is
            // the missing periodic axis (spec § 5.3.1).
            const std::vector<Finding> cfs = QvfReader::validate_file(path);
            bool only_periodicity = true;
            for (const auto& f : cfs)
                if (f.severity == Finding::Severity::Error)
                    only_periodicity = only_periodicity && f.code == "E_MEMBER_MALFORMED" &&
                                       f.where == "/sections/1/members/basis";
            CHECK(has_errors(cfs) && only_periodicity,
                  "C: the bloch section's only finding is its molecular structure");
            QvfReader r = QvfReader::open(path);
            CHECK(r.read_binary_member("wf_bloch", "mo_coefficients").shape ==
                      std::vector<int64_t>({2, 1, 2, 2, 2}),
                  "C: bloch coefficients are [2, 1, 2, 2, 2]");
            CHECK(r.read_binary_member("wf_bloch", "mo_coefficients").as_doubles()[5] == 0.625,
                  "C: bloch coefficients read back");
            CHECK(r.read_json_member("wf_bloch", "mo_metadata").at("spin").as_string() ==
                      "restricted",
                  "C: n_spin 1 writes a restricted set");
        } catch (const std::exception& ex) {
            CHECK(false, std::string("C: reading bloch archive: ") + ex.what());
        }
        std::remove(path);
    }

    // --- sketch documents: chemical identity, reaction.scheme, preview, job
    //     sets, opaque vendor members, remove / replace, and updating a
    //     qvf_version 2 archive (spec § 2.3, § 3.1, § 5.1, § 5.8, § 5.9,
    //     § 5.14, § 5.15, § 7.3) ---
    {
        const std::vector<Atom> water = {{"O", 8, {{0, 0, 0.117}}},
                                         {"H", 1, {{0, 0.757, -0.469}}},
                                         {"H", 1, {{0, -0.757, -0.469}}}};
        const std::vector<Atom> hydroxyl = {{"O", 8, {{0, 0, 0}}}, {"H", 1, {{0, 0, 0.97}}}};
        const std::vector<std::array<double, 3>> no_cell;
        const std::array<bool, 3> molecule{{false, false, false}};
        auto bytes_of = [](const std::string& s) {
            return std::vector<uint8_t>(s.begin(), s.end());
        };
        // A PNG / JPEG signature and compressible filler: the validator and the
        // writer judge the signature only.
        const std::vector<uint8_t> png =
            bytes_of(std::string("\x89PNG\r\n\x1a\n", 8) + std::string(4096, '\0'));
        const std::vector<uint8_t> jpeg = bytes_of("\xff\xd8\xff\xe0" + std::string(64, 'j'));
        const std::vector<uint8_t> svg =
            bytes_of(R"(<svg xmlns="http://www.w3.org/2000/svg" width="8" height="8"/>)");
        auto n_sections = [](const QvfWriter& w) {
            return w.build_manifest().at("sections").size();
        };
        auto codes_of = [](const std::vector<Finding>& fs) {
            std::string s;
            for (const auto& f : fs) s += " " + f.code + "@" + f.where;
            return s;
        };
        auto count = [](const std::vector<Finding>& fs, const std::string& code) {
            size_t n = 0;
            for (const auto& f : fs) n += f.code == code;
            return n;
        };

        QvfWriter sw({"sketch", "1.0", "scheme"});
        sw.auto_metadata = false;
        sw.add_structure(water, molecule, no_cell, {}, "s1");
        sw.add_structure(hydroxyl, molecule, no_cell, {}, "s2");
        sw.add_structure(hydroxyl, molecule, no_cell, {}, "s3");
        Json water_id = Json::object();
        water_id.set("charge", 0).set("multiplicity", 1).set("name", "water").set("smiles", "O")
            .set("inchi", "InChI=1S/H2O/h1H2").set("inchikey", "XLYOFNOQVPJJNP-UHFFFAOYSA-N");
        CHECK(sw.set_structure_identity("s1", water_id) == "s1", "identity: water");
        sw.set_section_field("s1", "label", "1");
        sw.set_structure_identity("s2", Json::parse(R"({"charge": -1, "multiplicity": 1})"));
        sw.set_structure_identity("s3", Json::parse(R"({"charge": 0, "multiplicity": 2})"));
        // Refusals change nothing.
        const Json s3_before = sw.build_manifest().at("sections").at(2);
        for (const char* bad : {R"({"multiplicity": 1})", R"({"charge": 10})",
                                R"({"multiplicity": 12})", R"({"charge": 1.5})",
                                R"({"multiplicity": 0})", R"({"spin": 1})",
                                R"({"name": ""})", R"({"smiles": "O O"})",
                                R"({"inchi": "InChl=1S/HO/h1H"})", R"({"inchi": "InChI="})",
                                R"({"inchikey": "xlyofnoqvpjjnp-UHFFFAOYSA-N"})", R"([1])"})
            CHECK(throws([&] { sw.set_structure_identity("s3", Json::parse(bad)); }),
                  std::string("identity refuses ") + bad);
        CHECK(sw.build_manifest().at("sections").at(2).dump(0) == s3_before.dump(0),
              "identity: a refusal changes nothing");
        CHECK(throws([&] { sw.set_structure_identity("nowhere", Json::object()); }) &&
                  throws([&] { sw.set_structure_identity("s1", Json::raw("{")); }),
              "identity: an unknown section or invalid JSON throws");
        sw.set_structure_identity("s3", Json::parse(R"({"name": "hydroxyl radical"})"));
        sw.set_structure_identity("s3", Json::parse(R"({"name": null})"));
        CHECK(!sw.build_manifest().at("sections").at(2).has("name") &&
                  sw.build_manifest().at("sections").at(2).at("multiplicity").as_int() == 2,
              "identity: null removes a field, the others stay");

        // reaction.scheme (spec § 5.14).
        const Json steps = Json::parse(R"([
            {"id": "r1", "label": "a",
             "reactants": ["s1", {"structure_ref": "s2", "coefficient": 2}],
             "products": [{"structure_ref": "s3", "x_role": "radical"}],
             "conditions": ["THF, -78 °C", "2 h"], "yield_percent": 87,
             "yield_text": "87%", "x_note": {"drawn": true}},
            {"id": "r2", "arrow": "equilibrium", "reactants": ["s3"],
             "reagents": ["s1"], "products": ["s2"], "label": null}])");
        const size_t before_scheme = n_sections(sw);
        auto with_step = [&](const char* step) {
            return Json::parse(std::string("[") + step + "]");
        };
        for (const char* bad : {
                 R"({"id": "r1", "reactants": ["nowhere"], "products": ["s2"]})",
                 R"({"id": "r1", "reactants": ["s1"], "products": ["s2", "s2"]})",
                 R"({"id": "r 1", "reactants": ["s1"], "products": ["s2"]})",
                 R"({"id": "r1", "reactants": [], "products": ["s2"]})",
                 R"({"id": "r1", "products": ["s2"]})",
                 R"({"id": "r1", "reactants": ["s1"], "products": [{"structure_ref": ""}]})",
                 R"({"id": "r1", "reactants": [{"structure_ref": "s1", "coefficient": 0}],
                     "products": ["s2"]})",
                 R"({"id": "r1", "reactants": [{"structure_ref": "s1", "role": 1}],
                     "products": ["s2"]})",
                 R"({"id": "r1", "reactants": ["s1"], "products": ["s2"], "arrow": "back"})",
                 R"({"id": "r1", "reactants": ["s1"], "products": ["s2"], "yield_percent": 101})",
                 R"({"id": "r1", "reactants": ["s1"], "products": ["s2"], "conditions": [1]})",
                 R"({"id": "r1", "reactants": ["s1"], "products": ["s2"], "energies": {}})",
                 R"({"id": "r1", "reactants": ["s1"], "products": ["s2"], "label": ""})",
                 R"(7)"})
            CHECK(throws([&] { sw.add_reaction_scheme(with_step(bad)); }),
                  std::string("reaction.scheme refuses ") + bad);
        CHECK(throws([&] { sw.add_reaction_scheme(Json::array()); }) &&
                  throws([&] { sw.add_reaction_scheme(Json::object()); }),
              "reaction.scheme needs a non-empty array of steps");
        CHECK(throws([&] {
                  sw.add_reaction_scheme(with_step(
                      R"({"id": "r1", "reactants": ["s1"], "products": ["s2"]},
                         {"id": "r1", "reactants": ["s2"], "products": ["s1"]})"));
              }),
              "reaction.scheme refuses a repeated step id");
        sw.add_bonds({{0, 1, 1.0}}, "bonds");
        CHECK(throws([&] {
                  sw.add_reaction_scheme(with_step(
                      R"({"id": "r1", "reactants": ["bonds"], "products": ["s2"]})"));
              }),
              "reaction.scheme refuses a participant of another kind");
        CHECK(n_sections(sw) == before_scheme + 1, "reaction.scheme refusals add nothing");
        CHECK(sw.add_reaction_scheme(steps, "Scheme 1") == "reaction_scheme",
              "reaction.scheme added");

        // A job set: two jobs, a run naming one of them (spec § 5.8, § 5.9).
        sw.add_job_spec(Json::parse(R"({"job_type": "molecular", "structure_ref": "s1"})"),
                        "job_s1");
        sw.add_job_spec(Json::parse(R"({"job_type": "molecular", "structure_ref": "s2"})"),
                        "job_s2");
        RunRecord rr;
        rr.program = "vibe-qc";
        rr.log_text = "converged\n";
        rr.job_ref = "job_s1";
        rr.run_status = "converged";
        RunRecord bad_ref = rr, bad_status = rr;
        bad_ref.job_ref = "s1";
        bad_status.run_status = "done";
        CHECK(throws([&] { sw.add_run_record(bad_ref, "run_bad"); }) &&
                  throws([&] { sw.add_run_record(bad_status, "run_bad"); }),
              "run.record refuses a job_ref to another kind and an unknown run_status");
        sw.add_run_record(rr, "run_s1");

        // preview (spec § 5.15).
        CHECK(throws([&] { sw.add_preview(svg); }), "preview: an SVG thumbnail is refused");
        CHECK(throws([&] { sw.add_preview(bytes_of("GIF89a")); }),
              "preview: unknown thumbnail bytes are refused");
        CHECK(throws([&] { sw.add_preview(png, {}, "", "preview", "image/jpeg"); }),
              "preview: a media type that the bytes contradict is refused");
        CHECK(throws([&] { sw.add_preview(png, bytes_of("<html/>")); }),
              "preview: an image that is no SVG, PNG or JPEG is refused");
        CHECK(throws([&] { sw.add_preview(png, svg, "", "preview", "", "image/gif"); }),
              "preview: an image media type outside the list is refused");
        CHECK(throws([&] { sw.add_preview(png, bytes_of("<svg>\xff</svg>")); }),
              "preview: an SVG that is not UTF-8 is refused");
        CHECK(sw.add_preview(png, svg, "the scheme") == "preview", "preview added");
        sw.add_preview(jpeg, {}, "", "preview_jpeg");

        // Opaque members of a vendor section (spec § 4.2, § 7.3).
        Json doc_members = Json::object();
        doc_members.set("document", Json::parse(R"({"format": "ket", "version": 2})"));
        CHECK(throws([&] {
                  sw.add_vendor_section("x_vibesketch.document", doc_members,
                                        {{"drawing", svg, "image"}});
              }) &&
                  throws([&] {
                      sw.add_vendor_section("x_vibesketch.document", doc_members,
                                            {{"document", svg, ""}});
                  }),
              "vendor: a bad media type or a repeated role is refused");
        sw.add_vendor_section("x_vibesketch.document", doc_members,
                              {{"drawing", svg, "image/svg+xml"},
                               {"raw", bytes_of("\x01\x02"), ""}});
        sw.add_vendor_section("x_vibesketch.picture", Json::object(),
                              {{"png", png, "image/png"}});
        CHECK(sw.add_opaque_member("x_vibesketch_document", "thumbnail", jpeg, "image/jpeg") ==
                  "x_vibesketch_document",
              "add_opaque_member onto a vendor section");
        CHECK(throws([&] { sw.add_opaque_member("x_vibesketch_document", "thumbnail", jpeg); }) &&
                  throws([&] { sw.add_opaque_member("s1", "picture", png, "image/png"); }) &&
                  throws([&] { sw.add_opaque_member("x_vibesketch_document", "p", png, "png"); }),
              "add_opaque_member refuses a repeated role, a canonical kind, a bad media type");

        const std::vector<uint8_t> bytes = sw.to_bytes();
        const std::vector<Finding> fs = QvfReader::validate_bytes(bytes);
        CHECK(fs.empty(), "the sketch document validates without findings:" + codes_of(fs));
        QvfReader r = QvfReader::from_bytes(bytes);
        const Json& s1 = r.section("s1");
        CHECK(s1.at("label").as_string() == "1" && s1.at("charge").as_int() == 0 &&
                  s1.at("inchikey").as_string() == "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
              "identity read back");
        const Json scheme = r.read_json_member("reaction_scheme", "scheme");
        const Json& r1 = scheme.at("steps").at(0);
        CHECK(r.section("reaction_scheme").at("label").as_string() == "Scheme 1" &&
                  r1.at("reactants").at(0).at("structure_ref").as_string() == "s1" &&
                  r1.at("reactants").at(1).at("coefficient").as_int() == 2 &&
                  r1.at("products").at(0).at("x_role").as_string() == "radical" &&
                  r1.at("yield_percent").as_int() == 87 && r1.at("x_note").has("drawn") &&
                  !scheme.at("steps").at(1).has("label") &&
                  scheme.at("steps").at(1).at("reagents").at(0).at("structure_ref").as_string() ==
                      "s1",
              "reaction.scheme normalised: participant objects, vendor keys kept, nulls dropped");
        CHECK(r.section("job_s2").at("structure_ref").as_string() == "s2" &&
                  r.section("run_s1").at("job_ref").as_string() == "job_s1" &&
                  r.section("run_s1").at("run_status").as_string() == "converged",
              "job set: section-level structure_ref, job_ref and run_status");
        const Json& preview = r.section("preview");
        CHECK(preview.at("label").as_string() == "the scheme" &&
                  preview.at("members").at("thumbnail").at("path").as_string() ==
                      "preview/thumbnail.png" &&
                  preview.at("members").at("thumbnail").at("media_type").as_string() ==
                      "image/png" &&
                  preview.at("members").at("image").at("path").as_string() ==
                      "preview/image.svg" &&
                  r.section("preview_jpeg").at("members").at("thumbnail").at("path").as_string() ==
                      "preview_jpeg/thumbnail.jpg" &&
                  r.read_member_bytes("preview", "thumbnail") == png,
              "preview members, paths and media types");
        CHECK(r.typed_manifest().section("x_vibesketch_document")->member("drawing")->media_type ==
                      std::string("image/svg+xml") &&
                  r.typed_manifest().section("x_vibesketch_document")->member("drawing")->format ==
                      MemberSpec::Format::Opaque &&
                  !r.section("x_vibesketch_document").at("members").at("raw").has("media_type") &&
                  r.section("x_vibesketch_picture").at("members").at("png").at("path").as_string() ==
                      "x_vibesketch_picture/png.bin",
              "opaque vendor members with their media types");
        {
            container::ZipReader zr(bytes);
            CHECK(zr.find("preview/thumbnail.png")->method == 0 &&
                      zr.find("preview/image.svg")->method == 0 &&
                      zr.find("x_vibesketch_picture/png.bin")->method == 8,
                  "preview members are stored uncompressed");
        }

        // remove / replace in a buffered writer: the members are dropped.
        CHECK(throws([&] { sw.remove_section("nowhere"); }), "remove_section: unknown id");
        const Json removed = sw.remove_section("preview_jpeg");
        CHECK(removed.at("kind").as_string() == "preview" && n_sections(sw) == 11 &&
                  !contains(sw.to_bytes(), "preview_jpeg/thumbnail.jpg"),
              "remove_section returns the section and drops its member");
        sw.remove_section("s3");
        sw.add_structure(water, molecule, no_cell, {}, "s3");
        sw.set_structure_identity("s3", Json::parse(R"({"charge": 0, "multiplicity": 3})"));
        const std::vector<uint8_t> replaced = sw.to_bytes();
        const std::vector<Finding> rfs = QvfReader::validate_bytes(replaced);
        CHECK(rfs.empty(), "a replaced structure validates without findings:" + codes_of(rfs));
        CHECK(QvfReader::from_bytes(replaced).section("s3").at("members").at("structure")
                      .at("path").as_string() == "s3/structure.json" &&
                  QvfReader::from_bytes(replaced).read_json_member("s3", "structure")
                          .at("atoms").size() == 3,
              "a buffered writer reuses the replaced section's paths");
        {   // A dangling reference after a removal is the validator's finding.
            QvfWriter dw({"sketch", "1.0", "dangling"});
            dw.add_structure(water, molecule, no_cell, {}, "s1");
            dw.add_job_spec(Json::parse(R"({"job_type": "molecular", "structure_ref": "s1"})"));
            dw.remove_section("s1");
            const std::vector<Finding> dfs = QvfReader::validate_bytes(dw.to_bytes());
            CHECK(count(dfs, "E_REF_UNRESOLVED") == 2, "removal leaves references dangling");
        }

        // remove / replace in an updating writer: old entries stay, new
        // members get fresh paths (spec § 2.3).
        const std::string upath = "sketch_update.qvf";
        sw.write(upath);
        {
            QvfWriter u = QvfWriter::open_for_update(upath);
            u.remove_section("s3");
            u.add_structure(hydroxyl, molecule, no_cell, {}, "s3");
            u.set_structure_identity("s3", Json::parse(R"({"charge": -1, "multiplicity": 1})"));
            u.add_preview(jpeg, {}, "", "preview_new");
            u.remove_section("preview_new");   // added in this session: dropped
            u.remove_section("preview");
            u.add_preview(png, svg);
            u.write();
        }
        std::vector<Finding> ufs = QvfReader::validate_file(upath);
        CHECK(!has_errors(ufs) && count(ufs, "W_ORPHAN_MEMBER") == 4,
              "an update replacing two sections validates; the old entries and the old "
              "manifest are orphans:" + codes_of(ufs));
        {
            QvfReader ur = QvfReader::open(upath);
            CHECK(ur.section("s3").at("members").at("structure").at("path").as_string() ==
                          "s3/structure.2.json" &&
                      ur.read_json_member("s3", "structure").at("atoms").size() == 2 &&
                      ur.section("preview").at("members").at("thumbnail").at("path").as_string() ==
                          "preview/thumbnail.2.png" &&
                      ur.section("preview").at("members").at("image").at("path").as_string() ==
                          "preview/image.2.svg" &&
                      !ur.has_section("preview_new"),
                  "replaced members get fresh paths; a section removed before write leaves "
                  "nothing");
            container::ZipReader zr = container::ZipReader::open(upath);
            CHECK(!zr.has("preview_new/thumbnail.jpg") && zr.has("s3/structure.json"),
                  "the archive keeps the old entry and never got the dropped one");
        }
        {   // A second replacement in a later update takes the next number.
            QvfWriter u2 = QvfWriter::open_for_update(upath);
            u2.remove_section("s3");
            u2.add_structure(water, molecule, no_cell, {}, "s3");
            u2.write();
            QvfReader ur = QvfReader::open(upath);
            CHECK(ur.section("s3").at("members").at("structure").at("path").as_string() ==
                      "s3/structure.3.json",
                  "a later replacement takes the next free number");
        }
        repack(upath, upath);
        ufs = QvfReader::validate_file(upath);
        CHECK(ufs.empty(), "the repacked update has no orphan:" + codes_of(ufs));
        std::remove(upath.c_str());

        // remove / replace in a streaming writer: the streamed entries are
        // dropped from the archive's directory.
        const std::string spath = "sketch_stream.qvf";
        {
            QvfWriter st({"sketch", "1.0", "stream"}, spath);
            st.add_structure(water, molecule, no_cell, {}, "s");
            st.add_preview(png);
            st.remove_section("s");
            st.remove_section("preview");
            st.add_structure(hydroxyl, molecule, no_cell, {}, "s");
            st.add_preview(jpeg);
            Json blob = Json::object();
            blob.set("note", Json::object());
            st.add_vendor_section("x_sketch.blob", blob);
            MemberSink sink = st.begin_binary("x_sketch_blob", "extra", DType::Int8, {1});
            CHECK(throws([&] { st.remove_section("s"); }),
                  "remove_section refuses while a member sink is open");
            sink.write(std::vector<uint8_t>{1});
            sink.close();
            st.write();
        }
        const std::vector<Finding> sfs = QvfReader::validate_file(spath);
        CHECK(sfs.empty(),
              "a streaming writer's replaced sections leave no orphan:" + codes_of(sfs));
        {
            QvfReader sr = QvfReader::open(spath);
            CHECK(sr.read_json_member("s", "structure").at("atoms").size() == 2 &&
                      sr.section("preview").at("members").at("thumbnail").at("path").as_string() ==
                          "preview/thumbnail.jpg",
                  "the streamed replacements read back");
        }
        std::remove(spath.c_str());

        // An archive stamped with the withdrawn qvf_version 2 is updated as 1
        // (spec § 3.1).
        const std::string vpath = "sketch_v2.qvf";
        {
            const std::string payload =
                R"({"atoms":[{"symbol":"He","position":[0,0,0],"atomic_number":2}],)"
                R"("pbc":[false,false,false],"lattice_vectors":null})";
            Json member = Json::object();
            member.set("path", "structure/structure.json").set("format", "json")
                .set("sha256", detail::Sha256::hex_of(payload));
            Json members = Json::object();
            members.set("structure", member);
            Json section = Json::object();
            section.set("id", "structure").set("kind", "structure").set("members", members);
            Json sections = Json::array();
            sections.push_back(section);
            Json src = Json::object();
            src.set("program", "vibe-qc").set("version", "0.12").set("calculation", "he");
            Json manifest = Json::object();
            manifest.set("qvf_version", 2)
                .set("schema_uri", "https://vibe-qc.org/spec/qvf/2/manifest.schema.json")
                .set("byte_order", "little").set("source", src).set("sections", sections);
            detail::ZipWriter zw;
            zw.add("structure/structure.json", payload);
            zw.add("manifest.json", manifest.dump(), false);
            const std::vector<uint8_t> v2 = zw.finish();
            std::ofstream fh(vpath, std::ios::binary);
            fh.write(reinterpret_cast<const char*>(v2.data()),
                     static_cast<std::streamsize>(v2.size()));
        }
        const std::vector<Finding> v2fs = QvfReader::validate_file(vpath);
        CHECK(!has_errors(v2fs) && count(v2fs, "W_VERSION_WITHDRAWN") == 1,
              "a v2 archive is valid with W_VERSION_WITHDRAWN:" + codes_of(v2fs));
        {
            QvfWriter vu = QvfWriter::open_for_update(vpath);
            vu.set_structure_identity("structure", Json::parse(R"({"charge": 0, "multiplicity": 1})"));
            CHECK(throws([&] {
                      vu.set_structure_identity("structure",
                                                Json::parse(R"({"multiplicity": 2})"));
                  }),
                  "identity is checked against an archived structure's atoms");
            vu.write();
        }
        {
            QvfReader vr = QvfReader::open(vpath);
            CHECK(vr.manifest().at("qvf_version").as_int() == 1 &&
                      vr.manifest().at("schema_uri").as_string() ==
                          "https://vibe-qc.org/spec/qvf/1/manifest.schema.json" &&
                      vr.section("structure").at("charge").as_int() == 0,
                  "the update of a v2 archive writes qvf_version 1 and the v1 schema_uri");
        }
        const std::vector<Finding> v1fs = QvfReader::validate_file(vpath);
        CHECK(!has_errors(v1fs) && count(v1fs, "W_VERSION_WITHDRAWN") == 0 &&
                  count(v1fs, "W_ORPHAN_MEMBER") == 1,
              "the updated archive is version 1; the old manifest is superseded:" +
                  codes_of(v1fs));
        std::remove(vpath.c_str());
    }

    // --- unrestricted NTOs: unequal/empty spin channels and paired weights ---
    {
        auto archive = [](bool empty_beta, bool mismatched) {
            QvfWriter nw({"nto-test", "1", "unrestricted NTO"});
            nw.add_structure({{"H", 1, {{0, 0, 0}}}, {"H", 1, {{0, 0, .74}}}});
            WavefunctionGTO wf;
            wf.shells = {{0, 0, {1.0}, {1.0}, true}, {1, 0, {1.0}, {1.0}, true}};
            wf.unrestricted = true;
            wf.orbital_kind = "nto";
            wf.mo_coefficients_alpha = {2, 2, {1, 0, 0, 1}};
            wf.mo_coefficients_beta = empty_beta ? Matrix{0, 2, {}} : Matrix{1, 2, {1, 0}};
            wf.occupations_alpha = empty_beta ? std::vector<double>{.8, .2} : std::vector<double>{.55, .15};
            wf.occupations_beta = empty_beta ? std::vector<double>{} : std::vector<double>{.3};
            wf.state_ref = "es"; wf.state_index = 1; wf.n_electrons = 2;
            wf.nto_role = "particle"; wf.section_id = "nto_particle";
            nw.add_wavefunction_gto(wf);
            wf.nto_role = "hole"; wf.section_id = "nto_hole";
            if (mismatched) wf.occupations_alpha[0] = .65;
            nw.add_wavefunction_gto(wf);
            nw.add_excited_states(Json::raw(R"([{"index":1,"excitation_energy":4.0,"nto_ref":"nto_hole"}])"), Json(), "", "es");
            return nw.to_bytes();
        };
        CHECK(!has_errors(QvfReader::validate_bytes(archive(false, false))), "unrestricted NTOs validate");
        CHECK(!has_errors(QvfReader::validate_bytes(archive(true, false))), "empty NTO beta channel validates");
        const auto fs = QvfReader::validate_bytes(archive(false, true));
        std::set<std::string> locations;
        for (const auto& f : fs) if (f.code == "E_MEMBER_MALFORMED") locations.insert(f.where);
        CHECK((locations == std::set<std::string>{"/sections/1/members/mo_metadata", "/sections/2/members/mo_metadata"}),
              "NTO pair weight mismatches use shared codes and locations");
        for (int bad = 0; bad < 5; ++bad) {
            QvfWriter nw({"nto-test", "1", "refusal"});
            WavefunctionGTO wf;
            wf.shells = {{0, 0, {1.0}, {1.0}, true}};
            wf.unrestricted = true; wf.orbital_kind = "nto";
            wf.mo_coefficients_alpha = {2, 1, {1, 1}};
            wf.mo_coefficients_beta = {0, 1, {}};
            wf.occupations_alpha = {.8, .2}; wf.state_ref = "es";
            wf.state_index = 1; wf.nto_role = "hole";
            if (bad == 0) wf.occupations_alpha = {.2, .8};
            if (bad == 1) wf.occupations_alpha = {.8, -.2};
            if (bad == 2) wf.occupations_alpha = {.8};
            if (bad == 3) wf.energies_alpha = {0};
            if (bad == 4) wf.state_ref.clear();
            bool rejected = false;
            try { nw.add_wavefunction_gto(wf); } catch (const std::runtime_error&) { rejected = true; }
            CHECK(rejected, "invalid unrestricted NTO writer input rejected");
            CHECK(nw.build_manifest().at("sections").size() == 0, "NTO refusal leaves no section");
        }
    }

    // --- v0.3 writer fields: explicit null and reference basis ---
    {
        QvfWriter w({"field-test", "1", "writer fields"});
        Json job = Json::object();
        job.set("job_type", "molecular").set("basis", Json()).set("functional", Json());
        w.add_job_spec(job, "job");
        WavefunctionGTO wf;
        wf.shells = {{0, 0, {1.0}, {1.0}, true}};
        wf.mo_coefficients = {1, 1, {1}};
        wf.orbital_kind = "iao";
        wf.localization_reference_basis = "minao";
        w.add_structure({{"H", 1, {{0, 0, 0}}}});
        w.add_wavefunction_gto(wf);
        QvfReader r = QvfReader::from_bytes(w.to_bytes());
        CHECK(r.read_json_member("job", "spec").at("basis").is_null(), "explicit basis null is preserved");
        CHECK(r.read_json_member("job", "spec").at("functional").is_null(), "explicit functional null is preserved");
        CHECK(r.read_json_member("wf", "mo_metadata").at("localization_reference_basis").as_string() == "minao",
              "localization reference basis is written");
        CHECK(!has_errors(r.validate_findings()), "v0.3 writer fields validate");
        QvfWriter bad({"field-test", "1", "refusal"});
        bool rejected = false;
        wf.localization_reference_basis = "";
        try { bad.add_wavefunction_gto(wf); } catch (const std::runtime_error&) { rejected = true; }
        CHECK(rejected && bad.build_manifest().at("sections").size() == 0, "empty localization basis rejected atomically");
        rejected = false;
        job.set("basis", 3);
        try { bad.add_job_spec(job); } catch (const std::runtime_error&) { rejected = true; }
        CHECK(rejected && bad.build_manifest().at("sections").size() == 0, "non-string/non-null basis rejected atomically");
    }

    // Duplicate/invalid IDs cannot append orphan entries, including when a
    // writer streams directly to an archive. Rejected calls leave bytes and
    // the manifest unchanged; valid calls still work afterwards.
    {
        auto exercise_ids = [](QvfWriter& w) {
            w.auto_metadata = false;
            w.add_structure({{"H", 1, {{0, 0, 0}}}}, {{false, false, false}}, {}, {}, "taken");
            const auto before = w.build_manifest().dump();
            const std::vector<std::function<void(const std::string&)>> adders = {
                [&](const std::string& id) { w.add_citations("@misc{test}", id); },
                [&](const std::string& id) { w.add_scf_history(Json::array(), id); },
                [&](const std::string& id) { w.add_bond_orders("mayer", Json::array(), id); },
                [&](const std::string& id) { w.add_structure_symmetry(Json::object(), id); },
                [&](const std::string& id) { w.add_job_spec(Json::object().set("job_type", "molecular"), id); },
                [&](const std::string& id) { w.add_spectrum("spectra.ir", Json::object().set("frequencies", Json::array()).set("intensities", Json::array()), "", id); },
                [&](const std::string& id) { w.add_vendor_section("x_test.data", Json::object().set("payload", Json::object()), false, "", id); }
            };
            for (const auto& add : adders) for (const std::string id : {"taken", "camera", "bad/id", "bad\n"}) {
                bool rejected = false;
                try { add(id); } catch (const std::runtime_error&) { rejected = true; }
                CHECK(rejected, "bad section id rejected");
                CHECK(w.build_manifest().dump() == before, "bad section id leaves manifest unchanged");
            }
        };
        QvfWriter w({"id-test", "1", "ID validation"});
        exercise_ids(w);
        QvfReader r = QvfReader::from_bytes(w.to_bytes());
        CHECK(!has_errors(r.validate_findings()), "buffered ID refusal leaves valid archive");
        for (const auto& f : r.validate_findings()) CHECK(f.code != "W_ORPHAN_MEMBER", "buffered ID refusal leaves no orphan member");
        const std::string path = "qvf-id-refusal-test.qvf";
        {
            QvfWriter sw({"id-test", "1", "ID validation"}, path);
            exercise_ids(sw);
            sw.write();
        }
        const auto sf = QvfReader::validate_file(path);
        CHECK(!has_errors(sf), "streamed ID refusal leaves a valid archive");
        for (const auto& f : sf) CHECK(f.code != "W_ORPHAN_MEMBER", "streamed ID refusal leaves no orphan member");
        std::remove(path.c_str());
    }

    if (g_failures == 0) {
        std::cout << "all C++ tests passed\n";
        return 0;
    }
    std::cerr << g_failures << " C++ test(s) failed\n";
    return 1;
}
