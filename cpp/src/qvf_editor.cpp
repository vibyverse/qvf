// License: Apache-2.0.
#include "qvf/editor.hpp"
#include <filesystem>
#include <random>
#include <set>
#include <sstream>
namespace qvf {
namespace facade_editor_detail {
std::filesystem::path reserve_directory(const std::string& target) {
    const auto parent = std::filesystem::absolute(target).parent_path();
    std::random_device random;
    for (int attempt = 0; attempt < 64; ++attempt) {
        std::ostringstream suffix; suffix << std::hex << random() << random();
        const auto candidate = parent / (".qvf-update-" + suffix.str());
        std::error_code error;
        if (std::filesystem::create_directory(candidate, error)) return candidate;
        if (error && error != std::errc::file_exists) throw std::filesystem::filesystem_error("reserve QVF update directory", candidate, error);
    }
    throw std::runtime_error("cannot reserve QVF update directory");
}
void cleanup(const std::filesystem::path& p) noexcept {
    if (!p.empty()) { std::error_code ignored; std::filesystem::remove_all(p, ignored); }
}
struct Temporary {
    std::filesystem::path directory;
    explicit Temporary(const std::string& target) : directory(reserve_directory(target)) {}
    ~Temporary() { cleanup(directory); }
    std::string archive() const { return (directory / "archive.qvf").string(); }
};
} // namespace facade_editor_detail
struct Editor::Impl {
    std::string target;
    std::filesystem::path directory;
    std::unique_ptr<Writer> writer;
    bool committed = false;
    ~Impl() { writer.reset(); facade_editor_detail::cleanup(directory); }
};
Editor::Editor(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
Editor::Editor(Editor&&) noexcept = default;
Editor& Editor::operator=(Editor&&) noexcept = default;
Editor::~Editor() = default;
Editor Editor::open(const std::string& path, bool atomic) {
    auto impl = std::make_unique<Impl>(); impl->target = path;
    std::string working = path;
    if (atomic) {
        impl->directory = facade_editor_detail::reserve_directory(path);
        working = (impl->directory / "archive.qvf").string();
        std::filesystem::copy_file(path, working);
    }
    impl->writer.reset(new Writer(QvfWriter::open_for_update(working)));
    return Editor(std::move(impl));
}
Writer& Editor::writer() {
    if (!impl_ || impl_->committed || !impl_->writer)
        throw std::runtime_error("Editor is committed or moved from");
    return *impl_->writer;
}
Section Editor::remove_section(std::string_view id) {
    Writer& w = writer(); w.require_open();
    return Section::from_json(w.writer_->remove_section(std::string(id)));
}
void Editor::set_run_status(RunStatus status, std::string_view id) {
    Writer& w = writer(); w.require_open();
    const Manifest m = w.manifest();
    std::string selected(id);
    if (selected.empty() && m.sections) {
        for (const auto& s : *m.sections) if (s.kind && *s.kind == "run.record") {
            if (!selected.empty()) throw std::invalid_argument("set_run_status needs an id when several runs exist");
            selected = s.id.value_or("");
        }
    }
    const Section* s = m.section(selected);
    if (!s || !s->kind || *s->kind != "run.record") throw std::invalid_argument("set_run_status requires a run.record section");
    const char* value = status == RunStatus::Running ? "running" : status == RunStatus::Converged ? "converged" : "failed";
    w.writer_->set_section_field(selected, "run_status", value);
}
void Editor::commit() {
    if (!impl_) throw std::runtime_error("Editor was moved from");
    if (impl_->committed) return;
    if (impl_->writer) {
        impl_->writer->finish();
        // Release the working copy before replacing the target with it:
        // Windows does not rename a file that is still open.
        impl_->writer.reset();
    }
    if (!impl_->directory.empty()) {
        // Never remove the original first: rename either replaces it
        // atomically or fails with the original still in place.
        std::filesystem::rename(impl_->directory / "archive.qvf", impl_->target);
    }
    impl_->committed = true;
    facade_editor_detail::cleanup(impl_->directory); impl_->directory.clear();
}
void Editor::repack(const std::string& in, const std::string& out, container::Method method) {
    facade_editor_detail::Temporary temporary(out);
    {
        auto zip = container::ZipReader::open(in);
        for (const auto& name : zip.duplicates()) if (name != "manifest.json")
            throw std::runtime_error("E_ZIP_DUPLICATE_ENTRY: " + name);
        const container::Limits limits;
        const auto bytes = zip.read("manifest.json", limits.max_json_bytes);
        const Json manifest = Json::parse(std::string(bytes.begin(), bytes.end()), limits.max_json_depth);
        const Manifest typed = Manifest::from_json(manifest);
        if (!typed.sections) throw std::invalid_argument("manifest sections are absent");
        std::set<std::string> paths;
        for (const auto& s : *typed.sections) if (s.members)
            for (const auto& member : *s.members) {
                if (!member.second.path || *member.second.path == "manifest.json") throw std::invalid_argument("invalid member path");
                paths.insert(*member.second.path);
            }
        container::ZipWriter target(std::make_unique<container::FileSink>(temporary.archive()));
        for (const auto& path : paths) {
            const auto* info = zip.find(path);
            if (!info) throw std::runtime_error("E_MEMBER_MISSING: " + path);
            const bool zip64 = info->zip64 || info->usize >= 0xFFFFFFFFull - info->usize / 4096 - 65536;
            auto entry = target.begin(path, method, zip64);
            zip.read_into(path, [&entry](const uint8_t* p, size_t n) { entry.write(p, n); }, limits);
            entry.close();
        }
        target.write("manifest.json", bytes, container::Method::Store); target.finish();
    }
    std::filesystem::rename(temporary.archive(), out);
}
} // namespace qvf
