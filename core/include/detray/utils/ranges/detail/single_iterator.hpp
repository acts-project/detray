/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "detray/definitions/qualifiers.hpp"

namespace detray::ranges::detail {

/// @brief Emulates iterator behaviour for a single value.
template <typename value_t>
struct single_iterator {

    /// @returns true if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const single_iterator<value_t> &rhs) const
        -> bool {
        return *m_value == *rhs.m_value;
    }

    /// @returns false if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(const single_iterator &rhs) const -> bool {
        return *m_value != *rhs.m_value;
    }

    /// Does nothing
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> single_iterator<value_t> {
        ++m_value;
        return *this;
    }

    /// Does nothing
    DETRAY_HOST_DEVICE
    constexpr auto operator--() -> single_iterator<value_t> {
        --m_value;
        return *this;
    }

    /// @returns the single value that the iterator points to - const
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const -> const value_t & { return *m_value; }

    /// @returns the single value that the iterator points to
    DETRAY_HOST_DEVICE
    constexpr auto operator*() -> value_t & { return *m_value; }

    /// Advance the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const value_t j) const
        -> single_iterator<value_t> {
        return {m_value + j};
    }

    value_t *m_value;
};
}  // namespace detray::ranges::detail
