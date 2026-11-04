// qvf/tensor_view.hpp — borrowed tensor input with bounded-memory conversion.
// License: Apache-2.0.
#ifndef QVF_TENSOR_VIEW_HPP
#define QVF_TENSOR_VIEW_HPP
#include <functional>
#include <variant>
#include "qvf/qvf.hpp"
namespace qvf {
enum class ByteOrder { Native, Little, Big };
// The caller owns data and keeps it alive through write_into() or Writer::add.
// Elements are contiguous, row-major. shape describes the accessible storage.
struct TensorView {
    const void* data = nullptr;
    DType dtype = DType::Float64;
    std::vector<int64_t> shape;
    ByteOrder byte_order = ByteOrder::Native;
    TensorView() = default;
    TensorView(const void* p, DType d, std::vector<int64_t> s,
               ByteOrder order = ByteOrder::Native)
        : data(p), dtype(d), shape(std::move(s)), byte_order(order) {}
    explicit TensorView(const Tensor& t)
        : data(t.bytes.data()), dtype(t.dtype), shape(t.shape), byte_order(ByteOrder::Little) {
        if (byte_size() != t.bytes.size()) throw std::invalid_argument("Tensor bytes do not match dtype/shape");
    }
    uint64_t element_count() const;
    uint64_t byte_size() const;
    // Emits little-endian chunks of at most 64 KiB. Integer narrowing and
    // non-finite/out-of-range float-to-integer conversion throw. Floating
    // narrowing uses the host IEEE conversion. A callback may already have
    // received bytes when a later conversion fails; discard that member.
    void write_into(const std::function<void(const uint8_t*, size_t)>& sink,
                    DType output_dtype) const;
    void write_into(const std::function<void(const uint8_t*, size_t)>& sink) const {
        write_into(sink, dtype);
    }
};
// Typed sections accept either borrowed producer storage or an owning tensor
// returned by Reader. Copying/moving this record never leaves a dangling view
// into its own Tensor: view() is constructed on demand.
class TensorData {
public:
    TensorData() = default;
    TensorData(Tensor t) : value_(std::move(t)) {}
    TensorData(TensorView t) : value_(std::move(t)) {}
    TensorView view() const {
        if (const auto* t = std::get_if<Tensor>(&value_)) return TensorView(*t);
        return std::get<TensorView>(value_);
    }
    const Tensor* owned() const { return std::get_if<Tensor>(&value_); }
    std::optional<DType> output_dtype; // omitted: preserve the input dtype
private:
    std::variant<Tensor, TensorView> value_ = TensorView{};
};
} // namespace qvf
#endif
