/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <limits>
#include <type_traits>

namespace detray {

namespace detail {

/// Invalid value for fundamental types - constexpr
template <typename T,
          typename std::enable_if_t<std::is_fundamental_v<T>, bool> = true>
DETRAY_HOST_DEVICE inline constexpr T invalid_value() noexcept {
    return std::numeric_limits<T>::max();
}

/// Invalid value for types that cannot be constructed constexpr, e.g. Eigen
template <typename T,
          typename std::enable_if_t<not std::is_fundamental_v<T> and
                                        std::is_default_constructible_v<T>,
                                    bool> = true>
DETRAY_HOST_DEVICE inline T invalid_value() noexcept {
    return T{};
}

}  // namespace detail

template <typename T,
          typename std::enable_if_t<std::is_fundamental_v<T>, bool> = true>
DETRAY_HOST_DEVICE inline constexpr bool is_invalid_value(
    const T value) noexcept {
    return (value == detail::invalid_value<T>());
}

template <typename T,
          typename std::enable_if_t<!std::is_fundamental_v<T>, bool> = true>
DETRAY_HOST_DEVICE inline constexpr bool is_invalid_value(
    const T& value) noexcept {
    return (value == detail::invalid_value<T>());
}

}  // namespace detray
