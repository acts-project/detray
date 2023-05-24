/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/tuple_helpers.hpp"

// Vecmem include(s)
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s)
#include <type_traits>

namespace detray {

namespace detail {

// Views for types that aggregate containers/other bufferable types

/// Empty buffer type for inheritance template resolution
struct dbase_buffer {};

/// Container buffer helper that aggregates multiple vecmem buffers and performs
/// compile-time checks.
template <bool /*check value*/, typename... buffer_ts>
class dmulti_buffer_helper {};

/// In case the checks fail
template <typename... buffer_ts>
class dmulti_buffer_helper<false, buffer_ts...> {};

/// @brief General buffer type that aggregates vecmem based buffer
/// implementations.
///
/// This is for detray types that hold multiple members that all define custom
/// buffer types of their own. The 'sub'-buffers are begin aggregated in this
/// helper and are extracted in the types contructor and then handed down to the
/// member constructors.
template <typename... buffer_ts>
struct dmulti_buffer_helper<true, buffer_ts...> : public dbase_buffer {
    detray::tuple<std::remove_reference_t<std::remove_cv_t<buffer_ts>>...>
        m_buffer;

    dmulti_buffer_helper() = default;

    /// Tie multiple buffers together
    DETRAY_HOST
    dmulti_buffer_helper(buffer_ts&&... buffers) {
        m_buffer = ::detray::tie(buffers...);
    }
};

/// Helper trait to determine if a type can be interpreted as a (composite)
/// vecemem buffer
/// @{
template <typename T, typename = void>
struct is_device_buffer : public std::false_type {};

template <typename T>
struct is_device_buffer<
    T, std::enable_if_t<
           std::is_base_of_v<detray::detail::dbase_buffer,
                             std::remove_reference_t<std::remove_cv_t<T>>>,
           void>> : public std::true_type {};

template <typename T>
inline constexpr bool is_device_buffer_v = is_device_buffer<T>::value;
/// @}

/// Helper trait to check whether a type has a vecmem buffer defined
/// @{
template <class T, typename = void>
struct get_buffer : public std::false_type {
    using type = void;
};

template <class T>
struct get_buffer<
    T, std::enable_if_t<
           detray::detail::is_device_buffer_v<typename T::buffer_type>, void>>
    : public std::true_type {
    using type = typename T::buffer_type;
};

template <class T>
inline constexpr bool is_bufferable_v = get_buffer<T>::value;

template <class T>
using get_buffer_t = typename get_buffer<T>::type;
/// @}

}  // namespace detail

/// The detray container buffer exists, if all contained buffer types also
/// derive from @c dbase_buffer.
template <typename... buffer_ts>
using dmulti_buffer = detray::detail::dmulti_buffer_helper<
    std::conjunction_v<detail::is_device_buffer<buffer_ts>...>, buffer_ts...>;

/// Type trait specializations for vecmem containers
/// @{

/// Specialized buffer for @c vecmem::vector containers
template <typename T>
using dvector_buffer = vecmem::data::vector_buffer<T>;

/// Specialization of 'is buffer' for @c vecmem::data::vector_buffer containers
template <typename T>
struct detail::is_device_buffer<vecmem::data::vector_buffer<T>, void>
    : public std::true_type {};

/// Specialization of 'is buffer' for constant @c vecmem::data::vector_buffer
/// containers
template <typename T>
struct detail::is_device_buffer<const vecmem::data::vector_buffer<T>, void>
    : public std::true_type {};

/// Specialization of the buffer getter for @c vecmem::vector
template <typename T>
struct detail::get_buffer<vecmem::vector<T>, void> : public std::true_type {
    using type = dvector_buffer<T>;
};

/// Specialization of the buffer getter for @c vecmem::vector
template <typename T>
struct detail::get_buffer<const vecmem::vector<T>, void>
    : public std::true_type {
    using type = dvector_buffer<const T>;
};

/// @brief Detray version of 'get_data' - non-const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <class T>
dvector_buffer<T> get_buffer(const dvector_view<T>& vec_view,
                             ::vecmem::memory_resource& mr,
                             ::vecmem::copy& cpy) {
    dvector_buffer<T> buff{vec_view.size(), mr};
    cpy(vec_view, buff);
    return buff;
}

/// @brief Detray version of 'get_data' - const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <class T>
dvector_buffer<const T> get_buffer(const dvector_view<const T>& vec_view,
                                   ::vecmem::memory_resource& mr,
                                   ::vecmem::copy& cpy) {
    dvector_buffer<const T> buff{vec_view.size(), mr};
    cpy(vec_view, buff);
    return buff;
}

/// @brief Detray version of 'get_data' - non-const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <class T,
          std::enable_if_t<detail::is_device_view_v<typename T::view_type>,
                           bool> = true>
dvector_view<T> get_data(dvector_buffer<T>& buff) {
    return vecmem::get_data(buff);
}

/// @brief Detray version of 'get_data' - const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <class T,
          std::enable_if_t<detail::is_device_view_v<typename T::view_type>,
                           bool> = true>
dvector_view<const T> get_data(dvector_buffer<const T>& buff) {
    return vecmem::get_data(buff);
}

/// @brief Detray version of 'get_data' - non-const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <class T,
          std::enable_if_t<detail::is_device_buffer_v<typename T::buffer_type>,
                           bool> = true>
typename T::buffer_type get_buffer(T& bufferable, vecmem::memory_resource& mr,
                                   vecmem::copy& cpy) {
    return detray::get_buffer(bufferable.get_data(), mr, cpy);
}
/// @}

}  // namespace detray
