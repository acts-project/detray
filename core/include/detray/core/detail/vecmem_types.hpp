// Vecmem include(s)
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/unique_ptr.hpp>

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
    using host = T;
    using buffer = vecmem::unique_alloc_ptr<T>;
    using view = T&;
    using const_view = const T&;
    using device = T&;
    using const_device = const T&;
};