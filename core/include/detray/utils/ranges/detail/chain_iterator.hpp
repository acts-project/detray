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
template <std::size_t I, typename iterator_coll_t>
struct chain_iterator {

    using iterator_t = std::decay_t<decltype(detray::detail::get<I>(
        std::declval<iterator_coll_t>()))>;
    using value_type = typename std::iterator_traits<iterator_t>::value_type;

    constexpr chain_iterator(iterator_coll_t &&begins, iterator_coll_t &&ends)
        : m_begins(std::forward<iterator_coll_t>(begins)),
          m_ends(std::forward<iterator_coll_t>(ends)),
          m_itr{detray::detail::get<I>(m_begins)} {}

    /// @returns true if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator==(
        const chain_iterator<I, iterator_coll_t> &rhs) const -> bool {
        return m_itr == rhs.m_itr;
    }

    /// @returns false if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(
        const chain_iterator<I, iterator_coll_t> &rhs) const -> bool {
        return m_itr != rhs.m_itr;
    }

    /// Does nothing
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> chain_iterator<I, iterator_coll_t> {
        ++m_itr;

        // Get the next iterator pair
        if (m_itr == detray::detail::get<I>(m_ends)) {
            if constexpr (I + 1 < sizeof(iterator_coll_t)) {
                return chain_iterator<I + 1, iterator_coll_t>(
                    std::move(m_begins, m_ends));
            }
        }
        // keep returning this iterator. If it reached the last iterator pair,
        // this gets compared to the correct global sentinel
        return *this;
    }

    /// @returns the single value that the iterator points to - const
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const -> const value_type & { return *m_itr; }

    /// @returns the single value that the iterator points to
    DETRAY_HOST_DEVICE
    constexpr auto operator*() -> value_type & { return *m_itr; }

    /// Advance the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const std::size_t j) const
        -> chain_iterator<I, iterator_coll_t> {
        return {m_itr + j};
    }

    iterator_coll_t m_begins, m_ends;
    iterator_t &m_itr;
};

}  // namespace detray::ranges::detail
