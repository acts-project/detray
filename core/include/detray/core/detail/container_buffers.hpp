/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/tuple_helpers.hpp"

// Vecmem include(s)
#include <vecmem/containers/data/vector_buffer.hpp>
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

/// @brief General buffer type that hierarchically aggregates vecmem based
/// buffer implementations.
///
/// This is for detray classes that hold multiple members that all define custom
/// buffer types of their own. The 'sub'-buffers are being aggregated in this
/// helper and are extracted again in the @c get_data call when building the
/// buffer views. Since they are packaged into the @c dmulti_buffer_helper
/// according to the class member hierarchy, transcribing them into a view type
/// results back in a view type that is compatible with what needs to be passed
/// to the constructors of the original class.
///
/// class -> view -> buffer -> view (passed to kernel) -> device-side class
template <typename... buffer_ts>
struct dmulti_buffer_helper<true, buffer_ts...> : public dbase_buffer {
    dtuple<buffer_ts...> m_buffer;

    dmulti_buffer_helper() = default;

    /// Tie multiple buffers together
    DETRAY_HOST
    explicit dmulti_buffer_helper(buffer_ts&&... buffers)
        : m_buffer(std::forward<buffer_ts>(buffers)...) {}
};

/// Helper trait to determine if a type can be interpreted as a (composite)
/// vecemem buffer. This is the case if it inherits from @c dbase_buffer or if
/// it matches one of the vecmem specializations
/// @{
template <typename T, typename = void>
struct is_buffer : public std::false_type {};

template <typename T>
struct is_buffer<
    T, std::enable_if_t<
           std::is_base_of_v<detray::detail::dbase_buffer,
                             std::remove_reference_t<std::remove_cv_t<T>>>,
           void>> : public std::true_type {};

/// Specialization of @c is_buffer for @c vecmem::data::vector_buffer containers
template <typename T>
struct is_buffer<vecmem::data::vector_buffer<T>, void> : public std::true_type {
};

/// Specialization of @c is_buffer for constant @c vecmem::data::vector_buffer
/// containers
template <typename T>
struct is_buffer<const vecmem::data::vector_buffer<T>, void>
    : public std::true_type {};

template <typename T>
inline constexpr bool is_buffer_v = is_buffer<T>::value;

/// @}

/// Helper trait to check whether a type has a [vecmem] buffer type defined
/// @{
template <class T, typename = void>
struct has_buffer : public std::false_type {
    using type = void;
};

template <class T>
struct has_buffer<
    T, std::enable_if_t<detray::detail::is_buffer_v<typename T::buffer_type>,
                        void>> : public std::true_type {
    using type = typename T::buffer_type;
};

/// Specialization of the buffer getter for @c vecmem::vector
template <typename T>
struct has_buffer<vecmem::vector<T>, void> : public std::true_type {
    using type = vecmem::data::vector_buffer<T>;
};

/// Specialization of the buffer getter for @c vecmem::vector - const
template <typename T>
struct has_buffer<const vecmem::vector<T>, void> : public std::true_type {
    using type = vecmem::data::vector_buffer<const T>;
};

template <class T>
inline constexpr bool has_buffer_v = has_buffer<T>::value;

template <class T>
using get_buffer_t = typename has_buffer<T>::type;
/// @}

}  // namespace detail

/// Specialized buffer for @c vecmem::vector containers
template <typename T>
using dvector_buffer = vecmem::data::vector_buffer<T>;

/// The detray container buffer exists, if all contained buffer types also
/// derive from @c dbase_buffer.
template <typename... buffer_ts>
using dmulti_buffer = detray::detail::dmulti_buffer_helper<
    std::conjunction_v<detail::is_buffer<buffer_ts>...>, buffer_ts...>;

/// @brief Get the buffer representation of a vecmem vector - non-const
template <class T>
dvector_buffer<T> get_buffer(const dvector_view<T>& vec_view,
                             ::vecmem::memory_resource& mr,
                             ::vecmem::copy& cpy) {
    dvector_buffer<T> buff{vec_view.size(), mr};
    cpy(vec_view, buff)->wait();;
    return buff;
}

/// @brief Recursively get the buffer representation of a composite view
///
/// Unwraps the view type at compile time and calls @c get_buffer on every view.
/// Then returns the resulting buffer objects and packages them recursively
/// into a @c multi_buffer .
///
/// @note This does not pick up the vecmem types.
template <
    typename... Ts,
    std::enable_if_t<detail::is_device_view_v<dmulti_view<Ts...>>, bool> = true>
auto get_buffer(const dmulti_view<Ts...>& data_view,
                vecmem::memory_resource& mr, vecmem::copy& cpy) {
    // Evaluate recursive buffer type
    // (e.g. dmulti_view<..., dmulti_view<dvector_view<T>, ...>, ...>
    //       => dmulti_buffer<..., dmulti_buffer<dvector_buffer<T>, ...>, ...>)
    using result_buffer_t = dmulti_buffer<decltype(detray::get_buffer(
        detail::get<Ts>(data_view.m_view), mr, cpy))...>;

    return result_buffer_t{std::move(
        detray::get_buffer(detail::get<Ts>(data_view.m_view), mr, cpy))...};
}

/// @brief Get the buffer representation of a composite object - non-const
template <class T, std::enable_if_t<detail::has_buffer_v<T>, bool> = true>
typename T::buffer_type get_buffer(T& bufferable, vecmem::memory_resource& mr,
                                   vecmem::copy& cpy) {
    return detray::get_buffer(bufferable.get_data(), mr, cpy);
}

/// Get the vecmem view type of a vector buffer
///
/// @note this is needed, because the corresponding vecmem::get_data() function
///       would return a reference to a view.
template <class T>
dvector_view<T> get_data(dvector_buffer<T>& buff) {
    return vecmem::get_data(buff);
}

/// @brief Unroll the composite buffer type
///
/// Unwraps the buffer type at compile time and calls @c get_data on every
/// buffer.
/// Then returns the resulting view objects and packages them into a
/// @c mutli_view to be ultimately passed on to the class constructor.
///
/// @note This does not pick up the vecmem types.
template <class... Ts>
auto get_data(dmulti_buffer<Ts...>& multi_buff) {
    // Evaluate recursive view type (reverse of 'get_buffer(dmulti_view)')
    using result_view_t =
        dmulti_view<decltype(detray::get_data(std::declval<Ts&>()))...>;

    return result_view_t{
        detray::get_data(detail::get<Ts>(multi_buff.m_buffer))...};
}

}  // namespace detray
