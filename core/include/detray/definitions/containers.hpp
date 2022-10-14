/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/type_traits.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/containers/vector.hpp>

// Thrust include(s)
#include <thrust/tuple.h>

// System include(s)
#include <array>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

namespace detray {
template <typename value_t, std::size_t kDIM>
using darray = std::array<value_t, kDIM>;

template <typename value_t>
using dvector = vecmem::vector<value_t>;

template <typename value_t>
using djagged_vector = vecmem::jagged_vector<value_t>;

template <typename key_t, typename value_t>
using dmap = std::map<key_t, value_t>;

template <class... types>
using dtuple = std::tuple<types...>;

/// @brief Bundle container type definitions
template <template <typename...> class vector_t = dvector,
          template <typename...> class tuple_t = dtuple,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class jagged_vector_t = djagged_vector,
          template <typename, typename> class map_t = dmap>
struct container_types {
    template <typename T>
    using vector_type = vector_t<T>;

    template <class... T>
    using tuple_type = tuple_t<T...>;

    template <typename T, std::size_t kDIM>
    using array_type = array_t<T, kDIM>;

    template <typename T>
    using jagged_vector_type = jagged_vector_t<T>;

    template <typename K, typename T>
    using map_type = map_t<K, T>;
};

/// Defining some common types
using host_container_types = container_types<>;
using device_container_types =
    container_types<vecmem::device_vector, thrust::tuple, darray,
                    vecmem::jagged_device_vector>;

/// How to obtain views for vecmem types
using vecmem::get_data;

/// Specialized view for @c vecmem::vector containers
template <typename value_t>
using dvector_view = detail::view_wrapper<value_t, vecmem::data::vector_view>;

/// Specialized view for @c vecmem::jagged_vector containers
template <typename value_t>
using djagged_vector_view =
    detail::view_wrapper<value_t, vecmem::data::jagged_vector_view>;

/// Specialization of the view getter for @c vecmem::vector
template <typename T>
struct detail::get_view<vecmem::vector<T>, void> : public std::true_type {
    using type = dvector_view<T>;
};

/// Specialization of the view getter for @c vecmem::vector
template <typename T>
struct detail::get_view<const vecmem::vector<T>, void> : public std::true_type {
    using type = dvector_view<const T>;
};

/// Specialization of the view getter for @c vecmem::jagged_vector
template <typename T>
struct detail::get_view<vecmem::jagged_vector<T>, void>
    : public std::true_type {
    using type = djagged_vector_view<T>;
};

/// Specialization of the view getter for @c vecmem::jagged_vector
template <typename T>
struct detail::get_view<const vecmem::jagged_vector<T>, void>
    : public std::true_type {
    using type = djagged_vector_view<const T>;
};

}  // namespace detray
