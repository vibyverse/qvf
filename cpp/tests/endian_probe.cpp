// endian_probe.cpp: byte-order probe for the C++ writer and reader.
// License: Apache-2.0.
//
//   endian_probe write <out.qvf>   write the probe archive on this host
//   endian_probe check <in.qvf>    validate an archive written on any host and
//                                  compare every member with the expected
//                                  little-endian bytes (spec § 2.2)
//
// The expected bytes are built with shifts, never from host memory, so a
// writer or reader that skips the byte swap on a big-endian host fails here.
// scripts/check_endian_roundtrip.py computes the same values independently.
// Every member has shape {2, 3, 4}; element i has the values of value_u64,
// value_f64 and value_f32 below.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "qvf/c_api.h"
#include "qvf/qvf_reader.hpp"
#include "qvf/reader.hpp"
#include "qvf/writer.hpp"

using namespace qvf;

namespace {

constexpr int kN = 24;
const std::vector<int64_t> kShape{2, 3, 4};

struct DT { const char* name; DType dtype; qvf_dtype c; size_t width; bool flt; };
const DT kTypes[] = {
    {"int8", DType::Int8, QVF_INT8, 1, false},
    {"uint8", DType::UInt8, QVF_UINT8, 1, false},
    {"int16", DType::Int16, QVF_INT16, 2, false},
    {"uint16", DType::UInt16, QVF_UINT16, 2, false},
    {"int32", DType::Int32, QVF_INT32, 4, false},
    {"uint32", DType::UInt32, QVF_UINT32, 4, false},
    {"int64", DType::Int64, QVF_INT64, 8, false},
    {"uint64", DType::UInt64, QVF_UINT64, 8, false},
    {"float32", DType::Float32, QVF_FLOAT32, 4, true},
    {"float64", DType::Float64, QVF_FLOAT64, 8, true},
};

// Distinct bytes in every position; the low `width` bytes are the value.
uint64_t value_u64(int i) { return 0x0102030405060708ULL * static_cast<uint64_t>(i + 1) + 0x8000000000000081ULL * static_cast<uint64_t>(i & 1); }
double value_f64(int i) {
    if (i == 0) return -0.0;
    if (i == 1) return std::ldexp(1.0, -1074);  // smallest subnormal
    const double x = std::ldexp((i + 1) * 0.1, i - 12);
    return (i & 1) ? -x : x;
}
float value_f32(int i) {
    if (i == 0) return -0.0f;
    if (i == 1) return std::ldexp(1.0f, -149);  // smallest subnormal
    return static_cast<float>(value_f64(i));
}

uint64_t bits_of(const DT& t, int i) {
    if (t.dtype == DType::Float64) { uint64_t b; double v = value_f64(i); std::memcpy(&b, &v, 8); return b; }
    if (t.dtype == DType::Float32) { uint32_t b; float v = value_f32(i); std::memcpy(&b, &v, 4); return b; }
    return value_u64(i);
}

// The expected little-endian member bytes, built with shifts.
std::vector<uint8_t> expected_le(const DT& t) {
    std::vector<uint8_t> out;
    for (int i = 0; i < kN; ++i) {
        const uint64_t b = bits_of(t, i);
        for (size_t k = 0; k < t.width; ++k) out.push_back(static_cast<uint8_t>(b >> (8 * k)));
    }
    return out;
}
std::vector<uint8_t> expected_be(const DT& t) {
    std::vector<uint8_t> out;
    for (int i = 0; i < kN; ++i) {
        const uint64_t b = bits_of(t, i);
        for (size_t k = t.width; k-- > 0;) out.push_back(static_cast<uint8_t>(b >> (8 * k)));
    }
    return out;
}
// The values in host memory, as a producer would hand them over.
std::vector<uint8_t> host_native(const DT& t) {
    std::vector<uint8_t> out(t.width * kN);
    for (int i = 0; i < kN; ++i) {
        const uint64_t b = bits_of(t, i);
        uint8_t* p = out.data() + t.width * static_cast<size_t>(i);
        switch (t.width) {
            case 1: { const uint8_t v = static_cast<uint8_t>(b); std::memcpy(p, &v, 1); break; }
            case 2: { const uint16_t v = static_cast<uint16_t>(b); std::memcpy(p, &v, 2); break; }
            case 4: { const uint32_t v = static_cast<uint32_t>(b); std::memcpy(p, &v, 4); break; }
            default: std::memcpy(p, &b, 8);
        }
    }
    return out;
}

// int16 values converted to float64 and float64 converted to float32: the
// TensorView conversion path.
std::vector<uint8_t> le_f64_of_int16() {
    std::vector<uint8_t> out;
    for (int i = 0; i < kN; ++i) {
        const double v = static_cast<int16_t>(static_cast<uint16_t>(value_u64(i)));
        uint64_t b; std::memcpy(&b, &v, 8);
        for (int k = 0; k < 8; ++k) out.push_back(static_cast<uint8_t>(b >> (8 * k)));
    }
    return out;
}
std::vector<uint8_t> le_f32_of_f64() {
    std::vector<uint8_t> out;
    for (int i = 0; i < kN; ++i) {
        const float v = static_cast<float>(value_f64(i));
        uint32_t b; std::memcpy(&b, &v, 4);
        for (int k = 0; k < 4; ++k) out.push_back(static_cast<uint8_t>(b >> (8 * k)));
    }
    return out;
}

bool host_is_little() { const uint16_t one = 1; uint8_t b; std::memcpy(&b, &one, 1); return b == 1; }

int failures = 0;
void check(bool ok, const std::string& what) {
    if (!ok) { ++failures; std::fprintf(stderr, "FAIL %s\n", what.c_str()); }
}

void write_probe(const std::string& path) {
    std::vector<std::vector<uint8_t>> keep;  // views borrow these until finish()
    {
        Writer w({"endian-probe", "1", "byte-order"}, path);
        w.set_auto_metadata(false);
        Section s; s.id = "probe"; s.kind = "x_qvf.endian_probe";
        std::vector<BinaryPayload> members;
        for (const DT& t : kTypes) {
            keep.push_back(host_native(t));
            members.push_back({std::string("native_") + t.name,
                               TensorView(keep.back().data(), t.dtype, kShape, ByteOrder::Native)});
            keep.push_back(expected_be(t));
            members.push_back({std::string("big_") + t.name,
                               TensorView(keep.back().data(), t.dtype, kShape, ByteOrder::Big)});
            keep.push_back(expected_le(t));
            members.push_back({std::string("little_") + t.name,
                               TensorView(keep.back().data(), t.dtype, kShape, ByteOrder::Little)});
        }
        keep.push_back(host_native(kTypes[2]));
        TensorData widened(TensorView(keep.back().data(), DType::Int16, kShape));
        widened.output_dtype = DType::Float64;
        members.push_back({"convert_int16_to_float64", widened});
        keep.push_back(host_native(kTypes[9]));
        TensorData narrowed(TensorView(keep.back().data(), DType::Float64, kShape));
        narrowed.output_dtype = DType::Float32;
        members.push_back({"convert_float64_to_float32", narrowed});
        w.add_payload(s, {}, members);

        // Streaming MemberSink takes caller-encoded little-endian bytes.
        Section st; st.id = "stream"; st.kind = "x_qvf.endian_probe";
        w.add_section(st);
        for (const DT& t : kTypes) {
            auto sink = w.begin_binary("stream", t.name, t.dtype, kShape);
            const auto le = expected_le(t);
            sink.write(le.data(), le.size());
            sink.close();
        }
        w.finish();
    }
}

// C API: host-native views through qvf_writer_add_section, appended with the
// legacy helpers' bytes (tensor_f64/f32/i32/i64) through begin_binary.
void write_probe_c(const std::string& path) {
    qvf_writer* w = nullptr;
    check(qvf_writer_create("endian-probe", "1", "byte-order-c", path.c_str(), &w) == QVF_OK, "C create");
    qvf_writer_set_auto_metadata(w, 0);
    std::vector<std::vector<uint8_t>> keep;
    std::vector<std::string> roles;
    std::vector<qvf_binary_input> bins;
    for (const DT& t : kTypes) { keep.push_back(host_native(t)); roles.push_back(std::string("native_") + t.name); }
    for (size_t k = 0; k < keep.size(); ++k) {
        const DT& t = kTypes[k];
        qvf_binary_input b{};
        b.role = roles[k].c_str();
        b.tensor = {keep[k].data(), keep[k].size(), t.c, kShape.size(), kShape.data(), QVF_NATIVE,
                    QVF_DTYPE_UNKNOWN};
        bins.push_back(b);
    }
    qvf_section_input in{};
    in.section_json = "{\"id\":\"probe_c\",\"kind\":\"x_qvf.endian_probe\"}";
    in.binary = bins.data(); in.binary_count = bins.size();
    check(qvf_writer_add_section(w, &in) == QVF_OK, "C add_section");

    check(qvf_writer_declare_section(w, "{\"id\":\"legacy\",\"kind\":\"x_qvf.endian_probe\"}") == QVF_OK, "C declare");
    std::vector<double> d; std::vector<int32_t> i32; std::vector<int64_t> i64;
    for (int i = 0; i < kN; ++i) {
        d.push_back(value_f64(i));
        i32.push_back(static_cast<int32_t>(static_cast<uint32_t>(value_u64(i))));
        i64.push_back(static_cast<int64_t>(value_u64(i)));
    }
    const std::pair<const char*, Tensor> legacy[] = {
        {"float64", tensor_f64(d, kShape)}, {"float32", tensor_f32(d, kShape)},
        {"int32", tensor_i32(i32, kShape)}, {"int64", tensor_i64(i64, kShape)}};
    for (const auto& kv : legacy) {
        const qvf_dtype c = kv.second.dtype == DType::Float64 ? QVF_FLOAT64
                          : kv.second.dtype == DType::Float32 ? QVF_FLOAT32
                          : kv.second.dtype == DType::Int32 ? QVF_INT32 : QVF_INT64;
        check(qvf_writer_begin_binary(w, "legacy", kv.first, c, kShape.size(), kShape.data()) == QVF_OK, "C begin");
        check(qvf_writer_write(w, kv.second.bytes.data(), kv.second.bytes.size()) == QVF_OK, "C write");
        check(qvf_writer_end_binary(w) == QVF_OK, "C end");
    }
    check(qvf_writer_finish(w) == QVF_OK, "C finish");
    qvf_writer_destroy(w);
}

const DT* type_named(const std::string& n) {
    for (const DT& t : kTypes) if (n == t.name) return &t;
    return nullptr;
}

// The member bytes every role must hold.
std::vector<uint8_t> want_for(const std::string& section, const std::string& role) {
    if (role == "convert_int16_to_float64") return le_f64_of_int16();
    if (role == "convert_float64_to_float32") return le_f32_of_f64();
    std::string n = role;
    const auto us = n.find('_');
    if (section != "stream" && section != "legacy" && us != std::string::npos) n = n.substr(us + 1);
    if (section == "legacy" && n == "float32") return le_f32_of_f64();
    const DT* t = type_named(n);
    return t ? expected_le(*t) : std::vector<uint8_t>{};
}

void check_archive(const std::string& path) {
    auto r = Reader::open(path);
    const auto report = r.validate();
    check(report.ok, "validator accepts " + path);
    for (const auto& f : report.findings) std::fprintf(stderr, "  %s %s\n", f.code.c_str(), f.message.c_str());
    size_t members = 0;
    for (const auto& s : r.manifest().sections.value_or(std::vector<Section>{})) {
        if (!s.members) continue;
        for (const auto& kv : *s.members) {
            const std::string where = *s.id + "/" + kv.first;
            const Tensor t = r.tensor(s, kv.first);
            check(t.bytes == want_for(*s.id, kv.first), where + " bytes");
            check(t.shape == kShape, where + " shape");
            ++members;
            // Read-side conversion re-emits little-endian bytes.
            if (kv.first == "native_int16") {
                std::vector<uint8_t> conv;
                TensorView(t).write_into([&](const uint8_t* p, size_t n) { conv.insert(conv.end(), p, p + n); },
                                         DType::Float64);
                check(conv == le_f64_of_int16(), where + " read-side int16 to float64");
            }
        }
    }
    check(members == 3 * 10 + 2 + 10 || members == 10 + 4, "member count " + std::to_string(members));
    // The compatibility reader decodes float64 into host doubles.
    QvfReader legacy = QvfReader::open(path);
    for (const auto& s : r.manifest().sections.value_or(std::vector<Section>{})) {
        if (!s.members) continue;
        for (const auto& kv : *s.members) {
            if (kv.second.dtype.value_or("") != "float64" || kv.first.rfind("convert", 0) == 0) continue;
            const auto host = legacy.read_binary_member(*s.id, kv.first).as_doubles();
            bool same = host.size() == kN;
            for (int i = 0; same && i < kN; ++i) {
                uint64_t a, b; const double want = value_f64(i);
                std::memcpy(&a, &host[static_cast<size_t>(i)], 8); std::memcpy(&b, &want, 8);
                same = a == b;
            }
            check(same, *s.id + "/" + kv.first + " as_doubles");
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) { std::fprintf(stderr, "usage: endian_probe write|check FILE\n"); return 2; }
    const std::string mode = argv[1], path = argv[2];
    std::printf("host byte order: %s\n", host_is_little() ? "little" : "big");
    try {
        if (mode == "write") {
            write_probe(path);
            write_probe_c(path.substr(0, path.size() - 4) + "_c.qvf");
        } else if (mode == "check") {
            check_archive(path);
        } else {
            std::fprintf(stderr, "unknown mode %s\n", mode.c_str()); return 2;
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL exception: %s\n", e.what()); return 1;
    }
    if (failures) { std::fprintf(stderr, "%d failure(s)\n", failures); return 1; }
    std::printf("%s ok\n", mode.c_str());
    return 0;
}
