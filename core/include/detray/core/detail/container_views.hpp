/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <tuple>
#include <type_traits>

namespace detray {

namespace detail {

/// Helper trait to check whether a type has a vecmem view defined
/// @{
template <class T, typename = void>
struct get_view : public std::false_type {
    using type = void;
};

template <class T>
struct get_view<T, typename T::view_type> : public std::true_type {
    using type = typename T::view_type;
};

template <class T>
inline constexpr bool is_viewable_v = get_view<T>::value;

template <class T>
using get_view_t = typename get_view<T>::type;
/// @}

// Views for types that aggregate containers/other viewable types

/// Empty view type for inheritance template resolution
struct dbase_view {};

/// Tag a vecmem view as @c dbase_view , so that it becomes recognizable in
/// detray as a vecmem view type
template <typename value_t, template <typename> class view_t>
struct view_wrapper : public dbase_view {
    /// The vecmem data view
    view_t<value_t> m_view{};

    /// Default constructor
    view_wrapper() = default;

    /// Conversion operator from a view of the same value type
    DETRAY_HOST
    view_wrapper(view_t<value_t>&& view) : m_view{view} {}

    /// Conversion operator from a view of the same value type
    DETRAY_HOST
    view_wrapper(const view_t<value_t>& view) : m_view{view} {}
};
template <typename value_t, template <typename> class view_t>
view_wrapper(const view_t<value_t>&) -> view_wrapper<const value_t, view_t>;

/// Container view helper that aggregates multiple vecmem views and performs
/// compile-time checks.
template <bool /*check value*/, typename... view_ts>
class dmulti_view_helper {};

/// In case the checks fail
template <typename... view_ts>
class dmulti_view_helper<false, view_ts...> {};

/// @brief General view type that aggregates vecmem based view implementations.
///
/// This is for detray types that hold multiple members that all define custom
/// view types of their own. The 'sub'-views are begin aggregated in this helper
/// and are extracted in the types contructor and then handed down to the
/// member constructors.
template <typename... view_ts>
struct dmulti_view_helper<true, view_ts...> : public dbase_view {
    std::tuple<std::remove_reference_t<std::remove_cv_t<view_ts>>...> m_views;

    dmulti_view_helper() = default;

    /// Tie multiple views together
    DETRAY_HOST
    dmulti_view_helper(view_ts&&... views) { m_views = std::tie(views...); }
};

/// Helper trait to determine if a type can be interpreted as a (composite)
/// vecemem view
/// @{
template <typename T, typename = void>
struct is_device_view : public std::false_type {};

template <typename T>
struct is_device_view<
    T, std::enable_if_t<
           std::is_base_of_v<detray::detail::dbase_view,
                             std::remove_reference_t<std::remove_cv_t<T>>>,
           void>> : public std::true_type {};

template <typename T>
inline constexpr bool is_device_view_v = is_device_view<T>::value;
/// @}

}  // namespace detail

/// The detray container view exists, if all contained view types also derive
/// from @c dbase_view.
template <typename... view_ts>
using dmulti_view = detray::detail::dmulti_view_helper<
    std::conjunction_v<detail::is_device_view<view_ts>...>, view_ts...>;

/// Forward declare a generic 'get_data' function, to which all subsequent
/// definitions will be template specializations.
template <typename T>
typename T::view_type get_data(T&);

}  // namespace detray
