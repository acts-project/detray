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

// Buffer wrappers for types that aggregate containers/other bufferable types

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
/// helper and are extracted in the types constructor and then handed down to
/// the member constructors.
template <typename... buffer_ts>
struct dmulti_buffer_helper<true, buffer_ts...> : public dbase_buffer {
    std::tuple<std::remove_reference_t<std::remove_cv_t<buffer_ts>>...>
        m_buffer;

    dmulti_buffer_helper() = default;

    /// Tie multiple buffers together
    DETRAY_HOST
    dmulti_buffer_helper(buffer_ts&&... buffers) {
        m_buffer = ::detray::detail::make_tuple<
            std::tuple,
            std::remove_reference_t<std::remove_cv_t<buffer_ts>>...>(
            std::forward<buffer_ts>(buffers)...);
    }
};

/// Helper trait to determine if a type can be interpreted as a (composite)
/// vecemem buffer
/// @{
template <typename T, typename = void>
struct is_device_buffer : public std::false_type {};

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

/// @}

/// The detray container buffer exists, if all contained buffer types also
/// derive from @c dbase_buffer.
template <typename... buffer_ts>
using dmulti_buffer = detray::detail::dmulti_buffer_helper<
    std::conjunction_v<detail::is_device_buffer<buffer_ts>...>, buffer_ts...>;

/// @brief Get the buffer representation of a vecmem vector - non-const
template <class T>
dvector_buffer<T> get_buffer(const dvector_view<T>& vec_view,
                             ::vecmem::memory_resource& mr,
                             ::vecmem::copy& cpy) {
    dvector_buffer<T> buff{vec_view.size(), mr};
    cpy(vec_view, buff);
    return buff;
}

/// @brief Get the buffer representation of a vecmem vector - const
template <class T>
dvector_buffer<const T> get_buffer(const dvector_view<const T>& vec_view,
                                   ::vecmem::memory_resource& mr,
                                   ::vecmem::copy& cpy) {
    dvector_buffer<const T> buff{vec_view.size(), mr};
    cpy(vec_view, buff);
    return buff;
}

/// @brief Unroll the composite view type
///
/// Unwraps the view tpye at compile time and calls @c get_buffer on every view.
/// Then returns the resulting buffer objects and packages them into a
/// @c mutli_buffer to be passed on to the next level.
///
/// @note This does not pick up the vecmem types.
template <class view_t, std::size_t... I,
          std::enable_if_t<detail::is_device_view_v<view_t>, bool> = true>
auto get_buffer(const view_t& data_view, vecmem::memory_resource& mr,
                vecmem::copy& cpy, std::index_sequence<I...> /*seq*/) {
    return dmulti_buffer<decltype(detray::get_buffer(
        detail::get<I>(data_view.m_view), mr, cpy))...>(
        std::move(
            detray::get_buffer(detail::get<I>(data_view.m_view), mr, cpy))...);
}

/// @brief Recursively get the buffer representation of a composite object
///
/// @note This does not pick up the vecmem types.
template <class view_t,
          std::enable_if_t<detail::is_device_view_v<view_t>, bool> = true>
auto get_buffer(const view_t& data_view, vecmem::memory_resource& mr,
                vecmem::copy& cpy) {
    return detray::get_buffer(
        data_view, mr, cpy,
        std::make_index_sequence<
            detail::tuple_size_v<decltype(data_view.m_view)>>{});
}

/// @brief Get the buffer representation of a composite object - non-const
///
/// @note This does not pick up the vecmem types.
template <class T, std::enable_if_t<
                       not detail::is_device_view_v<typename T::buffer_type>,
                       bool> = true>
typename T::buffer_type get_buffer(T& bufferable, vecmem::memory_resource& mr,
                                   vecmem::copy& cpy) {
    return detray::get_buffer(bufferable.get_data(), mr, cpy);
}
/// @}

/// @brief Get the view of a vecmem vector buffer - non-const
template <class T,
          std::enable_if_t<detail::is_device_view_v<typename T::view_type>,
                           bool> = true>
dvector_view<T> get_data(dvector_buffer<T>& buff) {
    return vecmem::get_data(buff);
}

/// @brief Get the view of a vecmem vector buffer - const
template <class T,
          std::enable_if_t<detail::is_device_view_v<typename T::view_type>,
                           bool> = true>
dvector_view<const T> get_data(dvector_buffer<const T>& buff) {
    return vecmem::get_data(buff);
}

/*template <class... Ts, std::size_t... I>
dmulti_view<std::remove_cv_t<std::remove_reference_t<decltype(detray::get_data(std::declval<Ts>()))>>...>
get_data(dmulti_buffer<Ts...>& multi_buff, std::index_sequence<I...>) {
    //using blub = typename
dmulti_view<std::remove_cv_t<std::remove_reference_t<decltype(detray::get_data(std::declval<Ts>()))>>...>::bla;
    return {detray::get_data(detail::get<I>(multi_buff.m_buffer))...};
}*/
template <class... Ts, std::size_t... I>
auto get_data(dmulti_buffer<Ts...>& multi_buff, std::index_sequence<I...>) {
    return dmulti_view<decltype(detray::get_data(std::declval<Ts&>()))...>(
        detray::get_data(detail::get<I>(multi_buff.m_buffer))...);
}

/// @brief Get the view of a @c multi_buffer - const
template <class... Ts>
auto get_data(dmulti_buffer<Ts...>& multi_buff) {
    return detray::get_data(multi_buff,
                            std::make_index_sequence<sizeof...(Ts)>{});
}

}  // namespace detray
