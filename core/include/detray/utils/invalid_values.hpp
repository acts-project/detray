/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <limits>
#include <type_traits>

namespace detray::detail {

/// Invalid value for fundamental types - constexpr
template <typename T,
          typename std::enable_if_t<
              std::is_fundamental<T>::value  || std::is_enum<T>::value , bool> = true>
DETRAY_HOST_DEVICE inline constexpr T invalid_value() noexcept {
    return std::numeric_limits<T>::max();
}

/// Invalid value for types that cannot be constructed constexpr, e.g. Eigen
template <typename T, typename std::enable_if_t<
                          !std::is_fundamental<T>::value  && !std::is_enum<T>::value  &&
                              std::is_default_constructible<T>::value ,
                          bool> = true>
DETRAY_HOST_DEVICE inline T invalid_value() noexcept {
    return T{};
}

template <typename T,
          typename std::enable_if_t<std::is_fundamental<T>::value , bool> = true>
DETRAY_HOST_DEVICE inline constexpr bool is_invalid_value(
    const T value) noexcept {
    if constexpr (std::is_signed<T>::value ) {
        if constexpr (std::is_floating_point<T>::value ) {
            return (math::fabs(value) == detail::invalid_value<T>());
        } else {
            return (math::abs(value) == detail::invalid_value<T>());
        }
    } else {
        return (value == detail::invalid_value<T>());
    }
}

template <typename T,
          typename std::enable_if_t<!std::is_fundamental<T>::value , bool> = true>
DETRAY_HOST_DEVICE inline constexpr bool is_invalid_value(
    const T& value) noexcept {
    return (value == detail::invalid_value<T>());
}

}  // namespace detray::detail
