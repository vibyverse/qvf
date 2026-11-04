// qvf_container.cpp — the ZIP container layer (qvf/container.hpp): ZIP64,
// streamed entries, data descriptors, append, and capped reads.
// License: Apache-2.0.
#include "qvf/container.hpp"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>

#include "qvf_crc32.hpp"
#include "qvf_deflate.hpp"
#include "qvf_inflate.hpp"
#include "qvf_sha256.hpp"

namespace qvf {
namespace container {

namespace zipfmt {

const uint32_t kSigLocal = 0x04034b50u;
const uint32_t kSigCentral = 0x02014b50u;
const uint32_t kSigEnd = 0x06054b50u;
const uint32_t kSigEnd64 = 0x06064b50u;
const uint32_t kSigLocator64 = 0x07064b50u;
const uint32_t kSigDescriptor = 0x08074b50u;
const uint64_t kMax32 = 0xFFFFFFFFull;
const uint32_t kMax16 = 0xFFFFu;
const uint16_t kFlagEncrypted = 0x0001;
const uint16_t kFlagDescriptor = 0x0008;
const uint16_t kFlagStrongEncryption = 0x0040;
const uint16_t kFlagUtf8 = 0x0800;
const size_t kChunk = size_t(1) << 20;

inline void put16(std::vector<uint8_t>& o, uint32_t v) {
    o.push_back(static_cast<uint8_t>(v & 0xFF));
    o.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
}
inline void put32(std::vector<uint8_t>& o, uint32_t v) {
    for (int i = 0; i < 4; ++i) o.push_back(static_cast<uint8_t>((v >> (8 * i)) & 0xFF));
}
inline void put64(std::vector<uint8_t>& o, uint64_t v) {
    for (int i = 0; i < 8; ++i) o.push_back(static_cast<uint8_t>((v >> (8 * i)) & 0xFF));
}
inline uint16_t get16(const uint8_t* p) {
    return static_cast<uint16_t>(p[0] | (p[1] << 8));
}
inline uint32_t get32(const uint8_t* p) {
    return static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) | (static_cast<uint32_t>(p[3]) << 24);
}
inline uint64_t get64(const uint8_t* p) {
    return static_cast<uint64_t>(get32(p)) | (static_cast<uint64_t>(get32(p + 4)) << 32);
}
inline uint32_t lo32(uint64_t v) { return static_cast<uint32_t>(v & 0xFFFFFFFFull); }

[[noreturn]] inline void fail(const std::string& message) {
    throw std::runtime_error("zip: " + message);
}

inline std::string qname(const std::string& name) { return "'" + name + "'"; }

// General-purpose bit 11 marks a UTF-8 name; an ASCII name needs no flag.
inline uint16_t name_flags(const std::string& name) {
    for (unsigned char c : name)
        if (c >= 0x80) return kFlagUtf8;
    return 0;
}

// Find the ZIP64 extra field (0x0001) in `extra`; returns its body or null.
inline const uint8_t* find_zip64_extra(const uint8_t* extra, size_t n, size_t& body_len) {
    size_t q = 0;
    while (q + 4 <= n) {
        const uint16_t id = get16(extra + q);
        const uint16_t len = get16(extra + q + 2);
        if (q + 4 + len > n) return nullptr;
        if (id == 0x0001) {
            body_len = len;
            return extra + q + 4;
        }
        q += 4u + len;
    }
    return nullptr;
}

// A local file header. With `zip64` it carries the ZIP64 extra field with
// both sizes (APPNOTE 4.5.3); with `descriptor` the CRC and sizes are zero
// and follow the data instead.
inline std::vector<uint8_t> local_header(const std::string& name, uint16_t flags,
                                         uint16_t method, uint32_t crc, uint64_t csize,
                                         uint64_t usize, bool zip64, bool descriptor) {
    std::vector<uint8_t> h;
    h.reserve(30 + name.size() + (zip64 ? 20 : 0));
    put32(h, kSigLocal);
    put16(h, zip64 ? 45 : 20);      // version needed
    put16(h, flags);
    put16(h, method);
    put16(h, 0);                    // mod time (fixed / zeroed)
    put16(h, 0);                    // mod date
    put32(h, descriptor ? 0 : crc);
    if (zip64) {
        // As Python's zipfile: zero with a descriptor, else saturated.
        put32(h, descriptor ? 0 : lo32(kMax32));
        put32(h, descriptor ? 0 : lo32(kMax32));
    } else {
        put32(h, descriptor ? 0 : lo32(csize));
        put32(h, descriptor ? 0 : lo32(usize));
    }
    put16(h, static_cast<uint32_t>(name.size()));
    put16(h, zip64 ? 20 : 0);       // extra length
    h.insert(h.end(), name.begin(), name.end());
    if (zip64) {
        put16(h, 0x0001);
        put16(h, 16);
        put64(h, descriptor ? 0 : usize);
        put64(h, descriptor ? 0 : csize);
    }
    return h;
}

}  // namespace zipfmt

using namespace zipfmt;

// ---------------------------------------------------------------------------
// Sinks
// ---------------------------------------------------------------------------

void Sink::write_at(uint64_t, const uint8_t*, size_t) {
    throw std::logic_error("zip: this sink cannot overwrite earlier bytes");
}

void MemorySink::write(const uint8_t* data, size_t n) {
    if (n) bytes_.insert(bytes_.end(), data, data + n);
}

void MemorySink::write_at(uint64_t offset, const uint8_t* data, size_t n) {
    if (offset > bytes_.size() || n > bytes_.size() - offset)
        throw std::logic_error("zip: write_at past the end of the memory sink");
    if (n) std::memcpy(bytes_.data() + offset, data, n);
}

struct FileSink::Impl {
    std::fstream f;
};

FileSink::FileSink(const std::string& path) : impl_(new Impl), path_(path) {
    impl_->f.open(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!impl_->f) throw std::runtime_error("cannot open for writing: " + path);
}

std::unique_ptr<FileSink> FileSink::open_at(const std::string& path, uint64_t offset) {
    std::unique_ptr<FileSink> s(new FileSink());
    s->impl_.reset(new Impl);
    s->path_ = path;
    s->impl_->f.open(path, std::ios::in | std::ios::out | std::ios::binary);
    if (!s->impl_->f) throw std::runtime_error("cannot open for update: " + path);
    s->impl_->f.seekp(static_cast<std::streamoff>(offset));
    if (!s->impl_->f) throw std::runtime_error("cannot seek in: " + path);
    s->pos_ = offset;
    s->truncate_on_close_ = true;
    return s;
}

FileSink::~FileSink() {
    if (impl_ && impl_->f.is_open()) impl_->f.close();
}

void FileSink::write(const uint8_t* data, size_t n) {
    if (!impl_ || !impl_->f.is_open()) throw std::logic_error("zip: write to a closed file sink");
    if (n == 0) return;
    impl_->f.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(n));
    if (!impl_->f) throw std::runtime_error("write failed: " + path_);
    pos_ += n;
}

void FileSink::write_at(uint64_t offset, const uint8_t* data, size_t n) {
    if (!impl_ || !impl_->f.is_open()) throw std::logic_error("zip: write to a closed file sink");
    if (offset > pos_ || n > pos_ - offset)
        throw std::logic_error("zip: write_at past the end of the file sink");
    impl_->f.seekp(static_cast<std::streamoff>(offset));
    impl_->f.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(n));
    impl_->f.seekp(static_cast<std::streamoff>(pos_));
    if (!impl_->f) throw std::runtime_error("write failed: " + path_);
}

void FileSink::close() {
    if (!impl_ || !impl_->f.is_open()) return;
    impl_->f.flush();
    const bool ok = static_cast<bool>(impl_->f);
    impl_->f.close();
    if (!ok) throw std::runtime_error("write failed: " + path_);
    if (truncate_on_close_) std::filesystem::resize_file(path_, pos_);
}

void CallbackSink::write(const uint8_t* data, size_t n) {
    if (n == 0) return;
    fn_(data, n);
    pos_ += n;
}

// ---------------------------------------------------------------------------
// ZipWriter
// ---------------------------------------------------------------------------

namespace impl {   // internal; named, not anonymous: the amalgamation is a header

struct WriterEntry {
    EntryInfo info;
    bool forced = false;                 // forced ZIP64: every central field saturated
    std::vector<uint8_t> raw_central;    // a record kept verbatim by append()
};

}  // namespace impl
using impl::WriterEntry;

struct ZipWriter::State {
    std::unique_ptr<Sink> sink;
    bool force_zip64 = false;
    bool allow_duplicates = false;       // detail::ZipWriter builds malformed test archives
    std::set<std::string> repeatable;    // names allow_repeated_name() lets repeat
    std::vector<WriterEntry> entries;
    std::map<std::string, size_t> index; // name -> last entry of that name
    bool entry_open = false;
    bool finished = false;
    bool failed = false;

    void reindex() {
        index.clear();
        for (size_t i = 0; i < entries.size(); ++i) index[entries[i].info.name] = i;
    }
    void check_writable() const {
        if (finished) throw std::logic_error("zip: the archive is already finished");
        if (failed) throw std::runtime_error("zip: an earlier write failed; the archive is unusable");
        if (entry_open) throw std::logic_error("zip: an entry is still being streamed");
    }
    void check_name(const std::string& name) const {
        if (name.empty()) throw std::invalid_argument("zip: empty entry name");
        if (name.size() > 0xFFFF) throw std::invalid_argument("zip: entry name too long");
        if (!allow_duplicates && !repeatable.count(name) && index.count(name))
            throw std::invalid_argument("zip: duplicate entry name " + qname(name));
    }
    void emit(const uint8_t* data, size_t n) {
        try {
            sink->write(data, n);
        } catch (...) {
            failed = true;
            throw;
        }
    }
    void emit(const std::vector<uint8_t>& v) { emit(v.data(), v.size()); }
    void patch(uint64_t offset, const std::vector<uint8_t>& v) {
        try {
            sink->write_at(offset, v.data(), v.size());
        } catch (...) {
            failed = true;
            throw;
        }
    }
    void add(WriterEntry e) {
        index[e.info.name] = entries.size();
        entries.push_back(std::move(e));
    }
};

struct EntrySink::State {
    // Shared, so an EntrySink never outlives the writer state it streams into.
    std::shared_ptr<ZipWriter::State> w;
    EntryInfo info;
    bool forced = false;
    bool descriptor = false;
    bool open = true;
    detail::Deflater deflater{detail::Deflater::Framing::Blocks};
    detail::Crc32 crc;
    detail::Sha256 sha;
    std::vector<uint8_t> scratch;
};

ZipWriter::ZipWriter(std::unique_ptr<Sink> out, bool force_zip64) {
    if (!out) throw std::invalid_argument("zip: null sink");
    std::shared_ptr<State> st = std::make_shared<State>();
    st->sink = std::move(out);
    st->force_zip64 = force_zip64;
    s_ = std::move(st);
}

ZipWriter::ZipWriter(ZipWriter&&) noexcept = default;
ZipWriter& ZipWriter::operator=(ZipWriter&&) noexcept = default;
ZipWriter::~ZipWriter() = default;

void ZipWriter::allow_duplicate_names() { s_->allow_duplicates = true; }

void ZipWriter::allow_repeated_name(const std::string& name) { s_->repeatable.insert(name); }

EntrySink ZipWriter::begin(const std::string& name, Method method, bool force_zip64) {
    State& s = *s_;
    s.check_writable();
    s.check_name(name);
    if (method != Method::Store && method != Method::Deflate)
        throw std::invalid_argument("zip: unsupported method");
    std::unique_ptr<EntrySink::State> st(new EntrySink::State);
    st->w = s_;
    st->forced = force_zip64 || s.force_zip64;
    st->descriptor = !s.sink->seekable();
    st->info.name = name;
    st->info.method = method;
    st->info.flags = static_cast<uint16_t>(name_flags(name) |
                                           (st->descriptor ? kFlagDescriptor : 0));
    st->info.zip64 = st->forced;
    st->info.local_offset = s.sink->position();
    s.emit(local_header(name, st->info.flags, static_cast<uint16_t>(method), 0, 0, 0,
                        st->info.zip64, st->descriptor));
    s.entry_open = true;
    return EntrySink(std::move(st));
}

EntryInfo ZipWriter::write(const std::string& name, const uint8_t* data, size_t n,
                           Method method, bool force_zip64) {
    State& s = *s_;
    s.check_writable();
    s.check_name(name);
    WriterEntry e;
    EntryInfo& i = e.info;
    i.name = name;
    i.flags = name_flags(name);
    i.crc32 = detail::crc32_of(data, n);
    i.size = n;
    i.local_offset = s.sink->position();
    // DEFLATE when that is smaller than the raw bytes, else STORE.
    std::vector<uint8_t> deflated;
    const uint8_t* body = data;
    uint64_t csize = n;
    i.method = Method::Store;
    if (method == Method::Deflate && n > 0) {
        deflated = detail::deflate(data, n);
        if (deflated.size() < n) {
            i.method = Method::Deflate;
            body = deflated.data();
            csize = deflated.size();
        }
    }
    i.compressed_size = csize;
    e.forced = force_zip64 || s.force_zip64;
    i.zip64 = e.forced || i.size >= kMax32 || csize >= kMax32;
    s.emit(local_header(name, i.flags, static_cast<uint16_t>(i.method), i.crc32, csize,
                        i.size, i.zip64, false));
    s.emit(body, static_cast<size_t>(csize));
    s.add(e);
    return i;
}

EntryInfo ZipWriter::write(const std::string& name, const std::vector<uint8_t>& bytes,
                           Method method, bool force_zip64) {
    return write(name, bytes.data(), bytes.size(), method, force_zip64);
}

EntryInfo ZipWriter::write(const std::string& name, const std::string& bytes,
                           Method method, bool force_zip64) {
    return write(name, reinterpret_cast<const uint8_t*>(bytes.data()), bytes.size(), method,
                 force_zip64);
}

bool ZipWriter::contains(const std::string& name) const { return s_->index.count(name) != 0; }

bool ZipWriter::drop(const std::string& name) {
    State& s = *s_;
    if (s.finished) throw std::logic_error("zip: the archive is already finished");
    auto it = s.index.find(name);
    if (it == s.index.end()) return false;
    s.entries.erase(s.entries.begin() + static_cast<std::ptrdiff_t>(it->second));
    s.reindex();
    return true;
}

std::vector<EntryInfo> ZipWriter::entries() const {
    std::vector<EntryInfo> out;
    out.reserve(s_->entries.size());
    for (const auto& e : s_->entries) out.push_back(e.info);
    return out;
}

void ZipWriter::set_force_zip64(bool force) { s_->force_zip64 = force; }

bool ZipWriter::finished() const { return s_->finished; }

Sink& ZipWriter::sink() { return *s_->sink; }

// One central-directory record. Saturated fields (forced, or past the 32-bit
// limits) move into a ZIP64 extra field, in APPNOTE order.
static void put_central(std::vector<uint8_t>& o, const WriterEntry& e) {
    if (!e.raw_central.empty()) {
        o.insert(o.end(), e.raw_central.begin(), e.raw_central.end());
        return;
    }
    const EntryInfo& i = e.info;
    const bool sat_u = e.forced || i.size >= kMax32;
    const bool sat_c = e.forced || i.compressed_size >= kMax32;
    const bool sat_o = e.forced || i.local_offset >= kMax32;
    std::vector<uint8_t> extra;
    if (sat_u || sat_c || sat_o) {
        put16(extra, 0x0001);
        put16(extra, static_cast<uint32_t>(8 * (int(sat_u) + int(sat_c) + int(sat_o))));
        if (sat_u) put64(extra, i.size);
        if (sat_c) put64(extra, i.compressed_size);
        if (sat_o) put64(extra, i.local_offset);
    }
    const uint32_t version = (i.zip64 || !extra.empty()) ? 45 : 20;
    put32(o, kSigCentral);
    put16(o, version);              // version made by
    put16(o, version);              // version needed
    put16(o, i.flags);
    put16(o, static_cast<uint16_t>(i.method));
    put16(o, 0);                    // time
    put16(o, 0);                    // date
    put32(o, i.crc32);
    put32(o, sat_c ? lo32(kMax32) : lo32(i.compressed_size));
    put32(o, sat_u ? lo32(kMax32) : lo32(i.size));
    put16(o, static_cast<uint32_t>(i.name.size()));
    put16(o, static_cast<uint32_t>(extra.size()));
    put16(o, 0);                    // comment
    put16(o, 0);                    // disk number
    put16(o, 0);                    // internal attrs
    put32(o, 0);                    // external attrs
    put32(o, sat_o ? lo32(kMax32) : lo32(i.local_offset));
    o.insert(o.end(), i.name.begin(), i.name.end());
    o.insert(o.end(), extra.begin(), extra.end());
}

void ZipWriter::finish() {
    State& s = *s_;
    s.check_writable();
    const uint64_t cd_offset = s.sink->position();
    std::vector<uint8_t> buf;
    uint64_t cd_size = 0;
    for (const auto& e : s.entries) {
        put_central(buf, e);
        if (buf.size() >= kChunk) {
            cd_size += buf.size();
            s.emit(buf);
            buf.clear();
        }
    }
    cd_size += buf.size();
    const uint64_t count = s.entries.size();
    const bool force = s.force_zip64;
    const bool need64 = force || count >= kMax16 || cd_offset >= kMax32 || cd_size >= kMax32;
    if (need64) {
        const uint64_t end64_offset = cd_offset + cd_size;
        put32(buf, kSigEnd64);
        put64(buf, 44);             // size of the rest of the record
        put16(buf, 45);             // version made by
        put16(buf, 45);             // version needed
        put32(buf, 0);              // this disk
        put32(buf, 0);              // disk with the central directory
        put64(buf, count);
        put64(buf, count);
        put64(buf, cd_size);
        put64(buf, cd_offset);
        put32(buf, kSigLocator64);
        put32(buf, 0);              // disk with the ZIP64 end record
        put64(buf, end64_offset);
        put32(buf, 1);              // total disks
    }
    const uint32_t count16 = (force || count >= kMax16) ? kMax16 : static_cast<uint32_t>(count);
    put32(buf, kSigEnd);
    put16(buf, 0);                  // disk numbers
    put16(buf, 0);
    put16(buf, count16);
    put16(buf, count16);
    put32(buf, (force || cd_size >= kMax32) ? lo32(kMax32) : lo32(cd_size));
    put32(buf, (force || cd_offset >= kMax32) ? lo32(kMax32) : lo32(cd_offset));
    put16(buf, 0);                  // comment length
    s.emit(buf);
    s.finished = true;
    try {
        s.sink->close();
    } catch (...) {
        s.failed = true;
        throw;
    }
}

// ---------------------------------------------------------------------------
// EntrySink
// ---------------------------------------------------------------------------

EntrySink::EntrySink(std::unique_ptr<State> state) : s_(std::move(state)) {}
EntrySink::EntrySink(EntrySink&&) noexcept = default;
EntrySink& EntrySink::operator=(EntrySink&& other) noexcept {
    if (this != &other) {
        abandon();
        s_ = std::move(other.s_);
    }
    return *this;
}
EntrySink::~EntrySink() { abandon(); }

void EntrySink::abandon() noexcept {
    if (s_ && s_->open) {
        s_->open = false;
        s_->w->entry_open = false;
    }
}

bool EntrySink::is_open() const { return s_ && s_->open; }
uint64_t EntrySink::size() const { return s_ ? s_->info.size : 0; }

void EntrySink::write(const std::string& bytes) {
    write(reinterpret_cast<const uint8_t*>(bytes.data()), bytes.size());
}

void EntrySink::write(const uint8_t* data, size_t n) {
    if (!is_open()) throw std::logic_error("zip: write to a closed entry");
    State& st = *s_;
    ZipWriter::State& w = *st.w;
    if (w.failed) throw std::runtime_error("zip: an earlier write failed; the archive is unusable");
    while (n > 0) {
        const size_t take = n < kChunk ? n : kChunk;
        st.crc.update(data, take);
        st.sha.update(data, take);
        st.info.size += take;
        if (st.info.method == Method::Store) {
            w.emit(data, take);
            st.info.compressed_size += take;
        } else {
            st.scratch.clear();
            st.deflater.write(data, take, st.scratch);
            if (!st.scratch.empty()) {
                w.emit(st.scratch);
                st.info.compressed_size += st.scratch.size();
            }
        }
        data += take;
        n -= take;
    }
}

EntryInfo EntrySink::close() {
    if (!is_open()) throw std::logic_error("zip: entry already closed");
    State& st = *s_;
    ZipWriter::State& w = *st.w;
    st.open = false;
    w.entry_open = false;
    if (w.failed) throw std::runtime_error("zip: an earlier write failed; the archive is unusable");
    EntryInfo& i = st.info;
    if (i.method == Method::Deflate) {
        st.scratch.clear();
        st.deflater.finish(st.scratch);
        w.emit(st.scratch);
        i.compressed_size += st.scratch.size();
    }
    i.crc32 = st.crc.value();
    if (!i.zip64 && (i.size >= kMax32 || i.compressed_size >= kMax32)) {
        w.failed = true;
        throw std::runtime_error("zip: entry " + qname(i.name) +
                                 " passed 4 GiB but was begun without force_zip64");
    }
    std::vector<uint8_t> fix;
    if (st.descriptor) {
        put32(fix, kSigDescriptor);
        put32(fix, i.crc32);
        if (i.zip64) {
            put64(fix, i.compressed_size);
            put64(fix, i.size);
        } else {
            put32(fix, lo32(i.compressed_size));
            put32(fix, lo32(i.size));
        }
        w.emit(fix);
    } else {
        put32(fix, i.crc32);
        w.patch(i.local_offset + 14, fix);
        fix.clear();
        if (i.zip64) {
            put64(fix, i.size);
            put64(fix, i.compressed_size);
            w.patch(i.local_offset + 30 + i.name.size() + 4, fix);
        } else {
            put32(fix, lo32(i.compressed_size));
            put32(fix, lo32(i.size));
            w.patch(i.local_offset + 18, fix);
        }
    }
    i.sha256 = st.sha.hex();
    WriterEntry e;
    e.info = i;
    e.forced = st.forced;
    w.add(std::move(e));
    return i;
}

// ---------------------------------------------------------------------------
// ZipReader
// ---------------------------------------------------------------------------

namespace impl {

class Source {
public:
    virtual ~Source() = default;
    virtual uint64_t size() const = 0;
    // Copy `n` bytes at `offset`; throws when the range is past the end.
    virtual void read(uint64_t offset, uint8_t* dst, size_t n) const = 0;
    // The whole archive, when it is in memory.
    virtual const uint8_t* contiguous() const { return nullptr; }
};

class MemorySource : public Source {
public:
    explicit MemorySource(std::vector<uint8_t> b) : bytes(std::move(b)) {}
    uint64_t size() const override { return bytes.size(); }
    void read(uint64_t offset, uint8_t* dst, size_t n) const override {
        if (offset > bytes.size() || n > bytes.size() - offset) fail("read past the end");
        if (n) std::memcpy(dst, bytes.data() + offset, n);
    }
    const uint8_t* contiguous() const override { return bytes.data(); }
    std::vector<uint8_t> bytes;
};

class FileSource : public Source {
public:
    explicit FileSource(const std::string& path) {
        in_.open(path, std::ios::in | std::ios::binary);
        if (!in_) throw std::runtime_error("cannot open: " + path);
        in_.seekg(0, std::ios::end);
        const std::streamoff end = in_.tellg();
        if (end < 0) throw std::runtime_error("cannot open: " + path);
        size_ = static_cast<uint64_t>(end);
    }
    uint64_t size() const override { return size_; }
    void read(uint64_t offset, uint8_t* dst, size_t n) const override {
        if (offset > size_ || n > size_ - offset) fail("read past the end");
        if (n == 0) return;
        std::lock_guard<std::mutex> lock(mu_);
        in_.clear();
        in_.seekg(static_cast<std::streamoff>(offset));
        in_.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(n));
        if (static_cast<size_t>(in_.gcount()) != n) fail("unexpected end of file");
    }

private:
    mutable std::ifstream in_;
    mutable std::mutex mu_;
    uint64_t size_ = 0;
};

}  // namespace impl
using impl::Source;
using impl::MemorySource;
using impl::FileSource;

struct ZipReader::Impl {
    std::unique_ptr<Source> src;
    std::vector<ZipEntry> entries;
    std::vector<std::vector<uint8_t>> raw_central;   // each record, verbatim
    std::map<std::string, size_t> index;             // name -> last record
    std::vector<std::string> duplicates;
    uint64_t cd_offset = 0;
    uint64_t total = 0;

    std::vector<uint8_t>& memory() {
        MemorySource* m = dynamic_cast<MemorySource*>(src.get());
        if (!m) throw std::logic_error("zip: not an in-memory archive");
        return m->bytes;
    }

    std::vector<uint8_t> bytes_at(uint64_t offset, size_t n) const {
        std::vector<uint8_t> out(n);
        src->read(offset, out.data(), n);
        return out;
    }

    void parse();
};

void ZipReader::Impl::parse() {
    const uint64_t size = src->size();
    if (size < 22) fail("too small");
    // The end record is in the last 22 + 65535 bytes (comment), and a ZIP64
    // locator, when present, is the 20 bytes before it.
    const uint64_t tail_len = std::min<uint64_t>(size, 22 + 65535 + 20);
    const uint64_t tail_off = size - tail_len;
    const std::vector<uint8_t> tail = bytes_at(tail_off, static_cast<size_t>(tail_len));
    size_t eocd = 0;
    bool found = false;
    const size_t lowest = tail.size() > 22 + 65535 ? tail.size() - 22 - 65535 : 0;
    for (size_t p = tail.size() - 22;; --p) {
        if (get32(&tail[p]) == kSigEnd && p + 22 + get16(&tail[p + 20]) <= tail.size()) {
            eocd = p;
            found = true;
            break;
        }
        if (p == lowest) break;
    }
    if (!found) fail("no end-of-central-directory");
    const uint8_t* e = &tail[eocd];
    uint64_t count = get16(e + 10);
    uint64_t cd_size = get32(e + 12);
    uint64_t cd_off = get32(e + 16);
    const bool saturated = count == kMax16 || cd_size == kMax32 || cd_off == kMax32;
    if (saturated && eocd >= 20 && get32(&tail[eocd - 20]) == kSigLocator64) {
        const uint8_t* loc = &tail[eocd - 20];
        if (get32(loc + 4) != 0 || get32(loc + 16) != 1)
            fail("multi-volume archives are not supported");
        const uint64_t end64 = get64(loc + 8);
        if (end64 > size || size - end64 < 56) fail("ZIP64 end record out of range");
        const std::vector<uint8_t> rec = bytes_at(end64, 56);
        if (get32(rec.data()) != kSigEnd64) fail("bad ZIP64 end record signature");
        if (get32(rec.data() + 16) != 0 || get32(rec.data() + 20) != 0)
            fail("multi-volume archives are not supported");
        if (count == kMax16) count = get64(rec.data() + 32);
        if (cd_size == kMax32) cd_size = get64(rec.data() + 40);
        if (cd_off == kMax32) cd_off = get64(rec.data() + 48);
    } else if (get16(e + 4) != 0 || get16(e + 6) != 0) {
        fail("multi-volume archives are not supported");
    }
    if (cd_off > size || cd_size > size - cd_off) fail("central directory out of range");
    if (count > cd_size / 46) fail("central directory too small for its entry count");
    cd_offset = cd_off;

    const std::vector<uint8_t> cd = bytes_at(cd_off, static_cast<size_t>(cd_size));
    entries.reserve(static_cast<size_t>(count));
    raw_central.reserve(static_cast<size_t>(count));
    size_t p = 0;
    for (uint64_t k = 0; k < count; ++k) {
        if (p + 46 > cd.size() || get32(&cd[p]) != kSigCentral)
            fail("bad central directory entry");
        const uint8_t* r = &cd[p];
        ZipEntry ent;
        ent.flags = get16(r + 8);
        ent.method = get16(r + 10);
        ent.crc32 = get32(r + 16);
        ent.csize = get32(r + 20);
        ent.usize = get32(r + 24);
        const uint16_t namelen = get16(r + 28);
        const uint16_t extralen = get16(r + 30);
        const uint16_t commentlen = get16(r + 32);
        ent.local_offset = get32(r + 42);
        const size_t next = p + 46u + namelen + extralen + commentlen;
        if (next > cd.size()) fail("truncated central directory entry");
        ent.name.assign(reinterpret_cast<const char*>(r + 46), namelen);
        size_t z64_len = 0;
        const uint8_t* z64 = find_zip64_extra(r + 46 + namelen, extralen, z64_len);
        ent.zip64 = z64 != nullptr;
        const bool sat_u = ent.usize == kMax32;
        const bool sat_c = ent.csize == kMax32;
        const bool sat_o = ent.local_offset == kMax32;
        if (sat_u || sat_c || sat_o) {
            // The ZIP64 extra field holds exactly the saturated fields, in
            // the order uncompressed size, compressed size, offset.
            if (!z64)
                fail("entry " + qname(ent.name) +
                     ": a size or offset is saturated but there is no ZIP64 extra field");
            size_t o = 0;
            auto take = [&](uint64_t& field) {
                if (o + 8 > z64_len)
                    fail("entry " + qname(ent.name) + ": truncated ZIP64 extra field");
                field = get64(z64 + o);
                o += 8;
            };
            if (sat_u) take(ent.usize);
            if (sat_c) take(ent.csize);
            if (sat_o) take(ent.local_offset);
        }
        auto it = index.find(ent.name);
        if (it != index.end()) {
            if (std::find(duplicates.begin(), duplicates.end(), ent.name) == duplicates.end())
                duplicates.push_back(ent.name);
            it->second = entries.size();       // last one wins
        } else {
            index[ent.name] = entries.size();
        }
        total = ent.usize > UINT64_MAX - total ? UINT64_MAX : total + ent.usize;
        raw_central.emplace_back(cd.begin() + static_cast<std::ptrdiff_t>(p),
                                 cd.begin() + static_cast<std::ptrdiff_t>(next));
        entries.push_back(std::move(ent));
        p = next;
    }
}

ZipReader ZipReader::open(const std::string& path) {
    ZipReader r;
    r.impl_.reset(new Impl);
    r.impl_->src.reset(new FileSource(path));
    r.impl_->parse();
    return r;
}

ZipReader ZipReader::from_bytes(std::vector<uint8_t> archive) {
    return ZipReader(std::move(archive));
}

ZipReader::ZipReader(std::vector<uint8_t> archive) : impl_(new Impl) {
    impl_->src.reset(new MemorySource(std::move(archive)));
    impl_->parse();
}

ZipReader::ZipReader(ZipReader&&) noexcept = default;
ZipReader& ZipReader::operator=(ZipReader&&) noexcept = default;
ZipReader::~ZipReader() = default;

const std::vector<ZipReader::Entry>& ZipReader::entries() const { return impl_->entries; }

std::vector<std::string> ZipReader::names() const {
    std::vector<std::string> out;
    out.reserve(impl_->index.size());
    std::map<std::string, bool> seen;
    for (const auto& e : impl_->entries)
        if (!seen[e.name]) {
            seen[e.name] = true;
            out.push_back(e.name);
        }
    return out;
}

const std::vector<std::string>& ZipReader::duplicates() const { return impl_->duplicates; }

bool ZipReader::has(const std::string& name) const {
    return impl_->index.find(name) != impl_->index.end();
}

const ZipReader::Entry* ZipReader::find(const std::string& name) const {
    auto it = impl_->index.find(name);
    return it == impl_->index.end() ? nullptr : &impl_->entries[it->second];
}

uint64_t ZipReader::archive_size() const { return impl_->src->size(); }
uint64_t ZipReader::central_directory_offset() const { return impl_->cd_offset; }
uint64_t ZipReader::total_size() const { return impl_->total; }

LocalHeader ZipReader::local_header(const Entry& e) const {
    const Source& src = *impl_->src;
    const uint64_t size = src.size();
    if (e.local_offset > size || size - e.local_offset < 30) fail("truncated local header");
    uint8_t h[30];
    src.read(e.local_offset, h, 30);
    if (get32(h) != kSigLocal) fail("bad local header signature");
    const uint16_t namelen = get16(h + 26);
    const uint16_t extralen = get16(h + 28);
    const uint64_t var_off = e.local_offset + 30;
    if (size - var_off < static_cast<uint64_t>(namelen) + extralen)
        fail("truncated local header");
    const std::vector<uint8_t> var =
        impl_->bytes_at(var_off, static_cast<size_t>(namelen) + extralen);
    // An empty name and extra field leave `var` empty (its data() may be null).
    if (namelen != e.name.size() ||
        (namelen != 0 && std::memcmp(var.data(), e.name.data(), namelen) != 0))
        fail("entry " + qname(e.name) + ": local header name differs from the central directory");
    LocalHeader lh;
    lh.data_offset = var_off + namelen + extralen;
    lh.flags = get16(h + 6);
    lh.crc32 = get32(h + 14);
    lh.csize = get32(h + 18);
    lh.usize = get32(h + 22);
    size_t z64_len = 0;
    const uint8_t* z64 = find_zip64_extra(var.data() + namelen, extralen, z64_len);
    if (z64) {
        lh.zip64 = true;
        // The local ZIP64 field carries both sizes (APPNOTE 4.5.3).
        if (z64_len >= 16) {
            if (lh.usize == kMax32) lh.usize = get64(z64);
            if (lh.csize == kMax32) lh.csize = get64(z64 + 8);
        }
    }
    return lh;
}

void ZipReader::check_limits(const Entry& e, const Limits& lim) const {
    if (e.usize > lim.max_member_bytes)
        fail("entry " + qname(e.name) + " declares " + std::to_string(e.usize) +
             " bytes, over the " + std::to_string(lim.max_member_bytes) +
             "-byte member cap");
    if (impl_->total > lim.max_total_bytes)
        fail("the archive declares " + std::to_string(impl_->total) +
             " bytes in total, over the " + std::to_string(lim.max_total_bytes) +
             "-byte total cap");
    if (lim.max_ratio > 0 && e.method != 0 && e.usize >= lim.ratio_exempt_bytes &&
        (e.csize == 0 ||
         static_cast<double>(e.usize) > lim.max_ratio * static_cast<double>(e.csize)))
        fail("entry " + qname(e.name) + " expands " + std::to_string(e.usize) +
             " bytes from " + std::to_string(e.csize) + ", over the compression-ratio cap");
}

// Decode one entry: to `fn` in chunks, or into `*collect`. `lim` (optional)
// re-checks the ratio cap as bytes are produced.
void ZipReader::stream_entry(const Entry& e, const ChunkFn& fn, const Limits* lim,
                             std::vector<uint8_t>* collect) const {
    if (e.flags & (kFlagEncrypted | kFlagStrongEncryption))
        fail("entry " + qname(e.name) + " is encrypted");
    if (e.method != 0 && e.method != 8)
        fail("entry " + qname(e.name) + ": unsupported compression method " +
             std::to_string(e.method));
    const Source& src = *impl_->src;
    const LocalHeader lh = local_header(e);
    const uint64_t size = src.size();
    if (lh.data_offset > size || e.csize > size - lh.data_offset)
        fail("truncated member data");
    detail::Crc32 crc;
    uint64_t produced = 0;
    auto over_ratio = [&]() {
        return lim && lim->max_ratio > 0 && e.method != 0 && produced >= lim->ratio_exempt_bytes &&
               static_cast<double>(produced) > lim->max_ratio * static_cast<double>(e.csize);
    };
    auto deliver = [&](const uint8_t* p, size_t n) {
        crc.update(p, n);
        produced += n;
        if (over_ratio())
            fail("entry " + qname(e.name) + " expands past the compression-ratio cap");
        if (collect)
            collect->insert(collect->end(), p, p + n);
        else
            fn(p, n);
    };
    const uint8_t* mem = src.contiguous();

    if (e.method == 0) {
        if (e.csize != e.usize) fail("stored member size mismatch: " + e.name);
        if (collect && e.usize <= static_cast<uint64_t>(SIZE_MAX))
            collect->reserve(static_cast<size_t>(e.usize));
        std::vector<uint8_t> chunk;
        uint64_t off = lh.data_offset;
        for (uint64_t left = e.csize; left > 0;) {
            const size_t n = static_cast<size_t>(std::min<uint64_t>(left, kChunk));
            if (mem) {
                deliver(mem + off, n);
            } else {
                chunk.resize(n);
                src.read(off, chunk.data(), n);
                deliver(chunk.data(), n);
            }
            off += n;
            left -= n;
        }
    } else {
        // Input: the entry's compressed bytes, from memory or in chunks.
        std::vector<uint8_t> in_chunk;
        uint64_t in_off = lh.data_offset;
        uint64_t in_left = e.csize;
        detail::InflateBitReader::Refill refill;
        const uint8_t* first = nullptr;
        size_t first_n = 0;
        if (mem) {
            first = mem + lh.data_offset;
            first_n = static_cast<size_t>(e.csize);
        } else {
            refill = [&](const uint8_t*& q) -> size_t {
                if (in_left == 0) return 0;
                const size_t n = static_cast<size_t>(std::min<uint64_t>(in_left, kChunk));
                in_chunk.resize(n);
                src.read(in_off, in_chunk.data(), n);
                in_off += n;
                in_left -= n;
                q = in_chunk.data();
                return n;
            };
        }
        detail::InflateBitReader br(first, first_n, refill);
        if (collect) {
            // Into one buffer; reserve no more than DEFLATE can expand to.
            const uint64_t bound = e.csize > (UINT64_MAX - 64) / 1032 ? UINT64_MAX
                                                                      : e.csize * 1032 + 64;
            const uint64_t reserve = std::min<uint64_t>(std::min<uint64_t>(e.usize, bound),
                                                        static_cast<uint64_t>(SIZE_MAX));
            detail::InflateWindow out(e.usize, detail::InflateWindow::Sink(),
                                      static_cast<size_t>(reserve));
            detail::inflate_stream(br, out);
            if (out.total() != e.usize)
                fail("member " + qname(e.name) + " inflates to a size other than declared");
            std::vector<uint8_t> bytes = out.take();
            crc.update(bytes.data(), bytes.size());
            produced = bytes.size();
            if (over_ratio())
                fail("entry " + qname(e.name) + " expands past the compression-ratio cap");
            *collect = std::move(bytes);
        } else {
            detail::InflateWindow out(e.usize, deliver);
            detail::inflate_stream(br, out);
            if (out.total() != e.usize)
                fail("member " + qname(e.name) + " inflates to a size other than declared");
        }
    }
    if (crc.value() != e.crc32) fail("entry " + qname(e.name) + ": CRC-32 mismatch");
}

std::vector<uint8_t> ZipReader::read(const std::string& name, uint64_t max_bytes) const {
    const Entry* e = find(name);
    if (!e) fail("member not found: " + name);
    if (e->usize > max_bytes)
        fail("member " + qname(name) + " declares " + std::to_string(e->usize) +
             " bytes, over the size cap");
    std::vector<uint8_t> out;
    stream_entry(*e, ChunkFn(), nullptr, &out);
    return out;
}

std::vector<uint8_t> ZipReader::read(const std::string& name, const Limits& limits) const {
    const Entry* e = find(name);
    if (!e) fail("member not found: " + name);
    check_limits(*e, limits);
    std::vector<uint8_t> out;
    stream_entry(*e, ChunkFn(), &limits, &out);
    return out;
}

void ZipReader::read_into(const std::string& name, const ChunkFn& fn,
                          const Limits& limits) const {
    const Entry* e = find(name);
    if (!e) fail("member not found: " + name);
    if (!fn) throw std::invalid_argument("zip: read_into needs a callback");
    check_limits(*e, limits);
    stream_entry(*e, fn, &limits, nullptr);
}

// ---------------------------------------------------------------------------
// ZipWriter::append (after ZipReader::Impl, whose records it adopts)
// ---------------------------------------------------------------------------

ZipWriter ZipWriter::append(const std::string& path) {
    std::unique_ptr<ZipWriter> w;
    {
        ZipReader r = ZipReader::open(path);
        w.reset(new ZipWriter(FileSink::open_at(path, r.central_directory_offset())));
        keep_entries(*w, r);
    }  // the reader's file handle is closed here
    return std::move(*w);
}

ZipWriter ZipWriter::append(std::vector<uint8_t> archive) {
    ZipReader r(std::move(archive));
    std::vector<uint8_t> prefix = std::move(r.impl_->memory());
    prefix.resize(static_cast<size_t>(r.central_directory_offset()));
    ZipWriter w(std::unique_ptr<Sink>(new MemorySink(std::move(prefix))));
    keep_entries(w, r);
    return w;
}

// Adopt the entries of an archive opened for appending, central records
// verbatim (their offsets stay valid: new entries go after them). A repeated
// name is kept as it is; whether to accept it is the caller's decision.
void ZipWriter::keep_entries(ZipWriter& w, const ZipReader& r) {
    const uint64_t cd_offset = r.central_directory_offset();
    const std::vector<ZipEntry>& es = r.entries();
    State& s = *w.s_;
    s.entries.reserve(es.size());
    for (size_t k = 0; k < es.size(); ++k) {
        const ZipEntry& z = es[k];
        if (z.local_offset >= cd_offset)
            throw std::runtime_error("zip: refusing to append: entry " + qname(z.name) +
                                     " lies after the central directory");
        WriterEntry e;
        e.info.name = z.name;
        e.info.method = static_cast<Method>(z.method);
        e.info.flags = z.flags;
        e.info.crc32 = z.crc32;
        e.info.compressed_size = z.csize;
        e.info.size = z.usize;
        e.info.local_offset = z.local_offset;
        e.info.zip64 = z.zip64;
        e.raw_central = r.impl_->raw_central[k];
        s.entries.push_back(std::move(e));
    }
    s.reindex();
}

}  // namespace container
}  // namespace qvf
