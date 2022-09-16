/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iostream>
#include <tuple>

#include "detray/definitions/qualifiers.hpp"

namespace detray::ranges::detail {

/// @brief Nested iterator to enumerate the elements of a range.
///
/// The enumeration is done by incrementing an index in lockstep with a wrapped
/// iterator of the range. Index and current iterator value are returned
/// using structured binding.
template <typename itr_t, typename value_t>
struct enumerate_iterator {

    /// Determine whether we reach end of range
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(
        const enumerate_iterator<itr_t, value_t> &rhs) const -> bool {
        return (m_iter != rhs.m_iter);
    }

    /// Increment
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> enumerate_iterator<itr_t, value_t> & {
        ++m_i;
        ++m_iter;
        return *this;
    }

    /// Tie them together for returning
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const { return std::tie(m_i, *m_iter); }

    /// Advance the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const value_t j) const
        -> enumerate_iterator<itr_t, value_t> {
        return {m_iter + j, m_i + j};
    }

    /// Start value of index sequence
    itr_t m_iter;
    value_t m_i;
};

}  // namespace detray::ranges::detail
