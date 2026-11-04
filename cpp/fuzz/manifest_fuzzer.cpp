// manifest_fuzzer.cpp: libFuzzer target over manifest.json and member payloads
// inside a well-formed ZIP, through the JSON parser, the typed model, the
// reader and every validation stage.
//
// zip_reader_fuzzer mutates raw archive bytes, so almost every mutation of a
// member fails its CRC or sha256 before the validator looks at it. This target
// instead takes the archive contents and builds a consistent ZIP from them:
//
//     <flags byte><manifest.json text>
//     \n=qvf=\n<entry name>\n<entry bytes>     (repeated, any number)
//
// Flags bit 0 rewrites every manifest object that has a "path" naming a given
// entry and a "sha256" string with the entry's real digest, so a mutated
// payload still reaches the payload and semantic checks. Bit 1 DEFLATEs the
// entries (STORE otherwise). Bit 2 turns the JSON Schema check off.
// scripts/fuzz_manifest_seeds.py turns .qvf archives into such inputs.
//
// Built only with -DQVF_FUZZ=ON (clang), like zip_reader_fuzzer:
//     manifest_fuzzer -max_total_time=60 SCRATCH_DIR SEED_DIR
#include <cstddef>
#include <cstdint>
#include <exception>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "qvf/container.hpp"
#include "qvf/json.hpp"
#include "qvf/model.hpp"
#include "qvf/qvf_reader.hpp"
#include "qvf/validation.hpp"
#include "../src/qvf_sha256.hpp"

namespace {

const std::string kSep = "\n=qvf=\n";

qvf::Json fix_digests(const qvf::Json& j, const std::map<std::string, std::string>& digests,
                      int depth) {
    if (depth > 64) return j;
    if (j.type() == qvf::Json::Type::Array) {
        qvf::Json out = qvf::Json::array();
        for (const auto& e : j.elements()) out.push_back(fix_digests(e, digests, depth + 1));
        return out;
    }
    if (j.type() != qvf::Json::Type::Object) return j;
    std::string path;
    for (const auto& kv : j.items())
        if (kv.first == "path" && kv.second.type() == qvf::Json::Type::String)
            path = kv.second.as_string();
    const auto hit = path.empty() ? digests.end() : digests.find(path);
    qvf::Json out = qvf::Json::object();
    for (const auto& kv : j.items()) {
        if (hit != digests.end() && kv.first == "sha256" &&
            kv.second.type() == qvf::Json::Type::String)
            out.set(kv.first, qvf::Json(hit->second));
        else
            out.set(kv.first, fix_digests(kv.second, digests, depth + 1));
    }
    return out;
}

}  // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size < 1) return 0;
    const uint8_t flags = data[0];
    const std::string in(reinterpret_cast<const char*>(data) + 1, size - 1);

    qvf::Limits lim;
    lim.max_member_bytes = uint64_t(1) << 24;
    lim.max_total_bytes = uint64_t(1) << 26;
    lim.max_json_bytes = uint64_t(1) << 22;
    qvf::ValidationOptions opt;
    opt.schema_check = (flags & 4) == 0;

    // Split into the manifest text and (name, bytes) entries.
    std::string manifest_text;
    std::vector<std::pair<std::string, std::string>> entries;
    size_t pos = in.find(kSep);
    manifest_text = in.substr(0, pos);
    while (pos != std::string::npos) {
        const size_t start = pos + kSep.size();
        pos = in.find(kSep, start);
        const std::string rec = in.substr(start, pos == std::string::npos ? pos : pos - start);
        const size_t nl = rec.find('\n');
        if (nl == std::string::npos || nl == 0) continue;
        entries.emplace_back(rec.substr(0, nl), rec.substr(nl + 1));
    }

    // The parser and the typed model on their own.
    try {
        qvf::Json j = qvf::Json::parse(manifest_text);
        if (flags & 1) {
            std::map<std::string, std::string> digests;
            for (const auto& e : entries) digests[e.first] = qvf::detail::Sha256::hex_of(e.second);
            j = fix_digests(j, digests, 0);
            manifest_text = j.dump(0);
        }
        if (j.type() == qvf::Json::Type::Object) {
            const qvf::Manifest m = qvf::Manifest::from_json(j);
            (void)m.to_json().dump(0);
        }
    } catch (const std::exception&) {
    }

    // A consistent archive from those contents.
    std::vector<uint8_t> bytes;
    try {
        auto sink = std::make_unique<qvf::container::MemorySink>();
        qvf::container::MemorySink* raw = sink.get();
        qvf::container::ZipWriter w(std::move(sink));
        const auto method = (flags & 2) ? qvf::container::Method::Deflate
                                        : qvf::container::Method::Store;
        w.write("manifest.json", manifest_text, method);
        for (const auto& e : entries) {
            if (w.contains(e.first)) continue;
            w.write(e.first, e.second, method);
        }
        w.finish();
        bytes = raw->take();
    } catch (const std::exception&) {
        return 0;   // a name the writer refuses: nothing to read
    }

    try {
        const qvf::QvfReader r = qvf::QvfReader::from_bytes(bytes, lim);
        (void)r.typed_manifest().to_json().dump(0);
        const qvf::Json& secs = r.sections();
        if (secs.type() == qvf::Json::Type::Array) {
            for (const auto& s : secs.elements()) {
                if (s.type() != qvf::Json::Type::Object || !s.has("id") || !s.has("members"))
                    continue;
                const qvf::Json& id = s.at("id");
                const qvf::Json& members = s.at("members");
                if (id.type() != qvf::Json::Type::String ||
                    members.type() != qvf::Json::Type::Object)
                    continue;
                for (const auto& m : members.items()) {
                    try {
                        (void)r.read_binary_member(id.as_string(), m.first);
                    } catch (const std::exception&) {
                    }
                    try {
                        (void)r.read_json_member(id.as_string(), m.first);
                    } catch (const std::exception&) {
                    }
                }
            }
        }
        (void)r.validate_findings(lim, opt);
    } catch (const std::exception&) {
    }

    try {
        const qvf::container::ZipReader zip(bytes);
        for (const auto level : {qvf::validate::Level::Manifest, qvf::validate::Level::Payloads,
                                 qvf::validate::Level::Semantics, qvf::validate::Level::All})
            (void)qvf::validate::run(zip, level, lim, opt);
    } catch (const std::exception&) {
    }
    return 0;
}
