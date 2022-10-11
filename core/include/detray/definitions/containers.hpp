/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
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

// Views for types that aggregate containers

/// Empty view type for inheritance template resolution
struct dbase_view {};

/// Tag a vecmem view as @c dbase_view , so that it becomes recognizable in
/// detray as a view type
template <typename value_t, template <typename> class view_t>
struct view_wrapper : public dbase_view {
    /// The vecmem data view
    view_t<value_t> m_view{};

    /// Default constructor
    view_wrapper() = default;

    /// Conversion operator from a view of the same value type
    DETRAY_HOST
    view_wrapper(view_t<value_t>&& view) : m_view{view} {}
};

/// Container view helper that aggregates multiple views and performs
/// compile-time checks.
template <bool /*check value*/, typename... view_ts>
class dmulti_view_helper {};

/// In case the checks fail
template <typename... view_ts>
class dmulti_view_helper<false, view_ts...> {};

/// General view type that aggregates vecmem based view implementations
template <typename... view_ts>
struct dmulti_view_helper<true, view_ts...> : public dbase_view {
    dtuple<view_ts...> m_views;

    dmulti_view_helper() = default;

    /// Tie multiple views together
    DETRAY_HOST
    dmulti_view_helper(view_ts&&... views) { m_views = std::tie(views...); }
};

/// The detray container view exists, if all contained view types also derive
/// from @c dbase_view.
template <typename... view_ts>
using dmulti_view = dmulti_view_helper<
    std::conjunction_v<std::is_base_of<dbase_view, view_ts>...>, view_ts...>;

/// Detray view for vector containers
template <typename value_t>
using dvector_view = view_wrapper<value_t, vecmem::data::vector_view>;

}  // namespace detray
