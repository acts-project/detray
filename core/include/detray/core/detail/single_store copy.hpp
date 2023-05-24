/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/containers/vector.hpp>

#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/// @tparam T The type of the collection data, e.g. transforms
/// @tparam context_t the context with which to retrieve the correct data.
template <typename T, typename context_t>
class single_store_base {
    public:
    using value_type = T;
    using context_type = context_t;

    /// Base vecmem types
    // TO DISCUSS: change naming to builder_type?
    using host_type = vecmem::vector<T>;
    using buffer_type = vecmem::data::vector_buffer<T>;
    using view_type = vecmem::data::vector_view<T>;
    using device_type = vecmem::device_vector<T>;

    using size_type = typename device_type::size_type;
    using link_type = size_type;
    using single_link = link_type;
    using range_link = std::array<link_type, 2>;
};

template <typename T, typename context_t>
DETRAY_HOST
class single_store_host : single_store_base<T, context_t> {

    using typename single_store_base<T, context_t>::host_type;
    using typename single_store_base<T, context_t>::view_type;
    using typename single_store_base<T, context_t>::size_type;
    using typename single_store_base<T, context_t>::link_type;
    
    public:
    
    single_store_host(const size_type size) : _data(size) {};

    // TO DISCUSS: Do we gain much, if anything from having a resource constructor here?

    single_store_host(const size_type size, const T &arg) : _data(size, arg) {};

    T& at(const link_type l) {
        return _data.at(l);
    }

    size_type size() {
        return static_cast<size_type>(_data.size());
    }

    view_type get_data() {
        return vecmem::get_data(_data);
    }

    void reserve(const size_type size) {
        _data.reserve(size);
    }

    void push_back(const T& value) {
        _data.push_back(value);
    }

    private:
    host_type _data;
};

template <typename T, typename context_t>
DETRAY_HOST
class single_store_buffer : single_store_base<T, context_t> {

    using typename single_store_base<T, context_t>::buffer_type;
    using typename single_store_base<T, context_t>::view_type;
    using typename single_store_base<T, context_t>::size_type;

    public:
    single_store_buffer(size_type size, vecmem::memory_resource& resource) 
        : _data(size, resource) {}

    single_store_buffer(vecmem::memory_resource& resource)
        : _data(resource) {}

    single_store_buffer(view_type data, vecmem::memory_resource& resource, vecmem::copy& cpy)
        : _data(data.size(), resource) {
            cpy(data, _data)->wait();
        }

    view_type get_data() {
        return vecmem::get_data(_data);
    }

    private:
    buffer_type _data;
};

template <typename T, typename context_t>
DETRAY_HOST_DEVICE
class single_store_view : single_store_base<T, context_t> {

    using typename single_store_base<T, context_t>::view_type;
    using typename single_store_base<T, context_t>::device_type;
    using typename single_store_base<T, context_t>::size_type;
    using typename single_store_base<T, context_t>::link_type;

    public:
    DETRAY_HOST_DEVICE
    single_store_view(view_type data)
        : _data(data) {}

    DETRAY_HOST_DEVICE
    T& at(link_type l) {
        device_type dev(_data);
        return dev.at(l);
    }

    DETRAY_HOST_DEVICE
    size_type size() const {
        return _data.size();
    }

    private:
    view_type _data;
};

template<typename T, typename context_t = empty_context>
struct single_store_types {
    using host = single_store_host<T, context_t>;
    using buffer = single_store_buffer<T, context_t>;
    using view = single_store_view<T, context_t>;
};

}  // namespace detray
