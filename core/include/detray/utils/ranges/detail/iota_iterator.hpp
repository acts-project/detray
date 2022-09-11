/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>

#include "detray/definitions/qualifiers.hpp"

namespace detray::ranges::detail {

/// Nested iterator to generate a range on demand
template <typename value_t>
struct iota_iterator {
    /// Start value of sequence
    value_t i;

    /// Determine whether we reach end of sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const iota_iterator<value_t> &rhs) const -> bool {
        return (i == rhs.i);
    }

    /// Determine whether we reach end of sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(const iota_iterator<value_t> &rhs) const -> bool {
        return (i != rhs.i);
    }

    /// Increment
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> iota_iterator<value_t> & {
        ++i;
        return *this;
    }

    /// Get the current value in the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const -> value_t { return i; }

    /// Advance the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const value_t j) const -> iota_iterator<value_t> {
        return {i + j};
    }
};

}  // namespace detray::ranges::detail
