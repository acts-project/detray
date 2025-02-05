// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <limits>
#include <type_traits>

namespace detray::detail {

/// Invalid value for fundamental types - constexpr
template <typename T>
    requires std::is_fundamental_v<T> ||
    std::is_enum_v<T> DETRAY_HOST_DEVICE constexpr T invalid_value() noexcept {
    return std::numeric_limits<T>::max();
}

/// Invalid value for types that cannot be constructed constexpr, e.g. Eigen
template <typename T>
    requires(!std::is_fundamental_v<T>) &&
    (!std::is_enum_v<T>)&&std::is_default_constructible_v<T> DETRAY_HOST_DEVICE
    inline T invalid_value() noexcept {
    return T{};
}

template <typename T>
requires std::is_fundamental_v<T> DETRAY_HOST_DEVICE constexpr bool
is_invalid_value(const T value) noexcept {
    if constexpr (std::is_signed_v<T>) {
        return (value == detail::invalid_value<T>() ||
                value == -detail::invalid_value<T>());
    } else {
        return (value == detail::invalid_value<T>());
    }
}

template <typename T>
requires(!std::is_fundamental_v<T>) DETRAY_HOST_DEVICE
    constexpr bool is_invalid_value(const T& value) noexcept {
    return (value == detail::invalid_value<T>());
}

}  // namespace detray::detail
