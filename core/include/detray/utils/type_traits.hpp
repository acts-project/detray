/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/accessor.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

/// Helper trait that should emulate std::remove_cvref_t in c++20
template <typename T>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

/// Helper trait for detecting when a type is a non-const version of another
///
/// This comes into play multiple times to enable certain constructors
/// conditionally through SFINAE.
///
/// @{
template <typename CTYPE, typename NCTYPE>
struct is_same_nc {
    static constexpr bool value = false;
};

template <typename TYPE>
struct is_same_nc<const TYPE, TYPE> {
    static constexpr bool value = true;
};
/// @}

/// Extract the value type of e.g. a container.
/// @{
template <typename T, typename = void>
struct get_value_type {
    using type = T;
};

template <typename T>
struct get_value_type<T*, void> {
    using type = T;
};

template <typename container_t>
struct get_value_type<
    container_t,
    std::enable_if_t<
        not std::is_same_v<typename remove_cvref_t<container_t>::value_type,
                           void>,
        void>> {
    using type = typename remove_cvref_t<container_t>::value_type;
};

template <typename T>
using get_value_t = typename get_value_type<T>::type;
/// @}

/// Helper trait that checks if a type models an interval of some value that can
/// be obtained with 'get'.
/// @{
template <typename TYPE, typename = void, typename = void>
struct is_interval : public std::false_type {};

template <typename TYPE>
struct is_interval<
    TYPE,
    std::enable_if_t<not std::is_arithmetic_v<std::remove_reference_t<TYPE>> and
                         std::is_arithmetic_v<std::remove_reference_t<decltype(
                             detray::detail::get<0>(std::declval<TYPE>()))>>,
                     void>,
    std::enable_if_t<not std::is_arithmetic_v<std::remove_reference_t<TYPE>> and
                         std::is_arithmetic_v<std::remove_reference_t<decltype(
                             detray::detail::get<1>(std::declval<TYPE>()))>>,
                     void>> : public std::true_type {};

template <typename TYPE>
inline constexpr bool is_interval_v = is_interval<TYPE>::value;
/// @}

}  // namespace detray::detail