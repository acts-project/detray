// Vecmem include(s)
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/unique_ptr.hpp>

namespace detray {

template <typename T, bool vector_like>
struct vecmem_types {};

template <typename T>
struct vecmem_types<T, true> {
    using host = vecmem::vector<T>;
    using buffer = vecmem::data::vector_buffer<T>;
    using view = vecmem::data::vector_view<T>;
    using const_view = vecmem::data::vector_view<const T>;
    using device = vecmem::device_vector<T>;
    using const_device = vecmem::vector<const T>;
};

template <typename T>
struct vecmem_types<T, false> {
    using value_type = T;

    using host = T;
    using buffer = vecmem::unique_alloc_ptr<T>;
    // should it be pointer instead of reference?
    using view = T*;
    using const_view = const T*;
    using device = T*;
    using const_device = const T*;
};

namespace test {

// This didn't work...
// template <typename T, std::enable_if_t<std::disjunction_v<is_specialization<T, vecmem::data::vector_buffer>, 
//                                                           is_specialization<T, vecmem::vector>, 
//                                                           is_specialization<T, std::vector>>>>
// auto get_data(const T& t) {
//     return vecmem::get_data(t);
// }

template <typename T>
auto get_data(vecmem::data::vector_buffer<T>& t) {
    return vecmem::get_data(t);
}

template <typename T>
auto get_data(vecmem::vector<T>& t) {
    return vecmem::get_data(t);
}

template <typename T>
T* get_data(vecmem::unique_alloc_ptr<T>& buff) {
    return buff.get();
}

template <typename T>
T* get_data(T& host) {
    return &host;
}

} // namespace test

} // namespace detray

