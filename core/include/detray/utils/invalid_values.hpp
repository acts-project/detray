/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <limits>
#include <type_traits>

namespace detray::detail {

/// Invalid value for fundamental types - constexpr
template <typename T,
          typename std::enable_if_t<std::is_fundamental_v<T>, bool> = true>
DETRAY_HOST_DEVICE inline constexpr T invalid_value() {
    return std::numeric_limits<T>::max();
}

/// Invalid value for types that cannot be constructed constexpr, e.g. Eigen
template <typename T,
          typename std::enable_if_t<not std::is_fundamental_v<T> and
                                        std::is_default_constructible_v<T>,
                                    bool> = true>
DETRAY_HOST_DEVICE inline T invalid_value() {
    return T{};
}

}  // namespace detray::detail