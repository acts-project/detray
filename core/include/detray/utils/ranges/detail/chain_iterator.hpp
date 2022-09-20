/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/type_traits.hpp"

namespace detray::detail {

template <typename container_t, typename = void>
struct get_value_type {
    using type = void;
};

template <typename container_t>
struct get_value_type<container_t, typename container_t::value_type> {
    using type = typename container_t::value_type;
};

template <typename container_t>
struct get_value_type<
    container_t,
    std::enable_if_t<
        std::is_class_v<std::remove_reference_t<decltype(detray::detail::get<0>(
            std::declval<container_t>()))>>,
        void>> {
    using type = std::decay_t<decltype(detray::detail::get<0>(
        std::declval<container_t>()))>;
};

template <typename T>
using get_value_type_t = typename get_value_type<T>::type;

}  // namespace detray::detail

namespace detray::ranges::detail {

/// @brief Sequntially iterate through multiple ranges.
///
/// Once the sentinel of one range is reached
template <typename iterator_coll_t,
          typename iterator_t =
              detray::detail::get_value_type_t<iterator_coll_t>>
struct chain_iterator {

    /// Default construction
    constexpr chain_iterator(iterator_coll_t &begins, iterator_coll_t &ends)
        : m_begins(begins),
          m_ends(ends),
          m_itr{detray::detail::get<0>(m_begins)},
          m_end{detray::detail::get<0>(m_ends)},
          m_idx{0} {}

    /// Fully parametrized construction
    constexpr chain_iterator(iterator_coll_t &begins, iterator_coll_t &ends,
                             iterator_t &current, iterator_t &end,
                             const std::size_t i)
        : m_begins(begins),
          m_ends(ends),
          m_itr{current},
          m_end{end},
          m_idx{i} {}

    /// Create a new chained iterator from an existing type
    template <std::size_t I>
    static constexpr auto create(iterator_coll_t &begins,
                                 iterator_coll_t &ends) {
        return chain_iterator(begins, ends, detray::detail::get<I>(begins),
                              detray::detail::get<I>(ends), I);
    }

    /// @returns true if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const chain_iterator &rhs) const -> bool {
        return m_itr == rhs.m_itr;
    }

    /// @returns false if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(const chain_iterator &rhs) const -> bool {
        return m_itr != rhs.m_itr;
    }

    /// Does nothing
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> chain_iterator<iterator_coll_t> {
        ++m_itr;

        // Get the next iterator pair
        if (m_itr == m_end and
            m_idx < sizeof(iterator_coll_t) / sizeof(iterator_t) - 1) {
            auto next =
                get_next_itr<sizeof(iterator_coll_t) / sizeof(iterator_t) -
                             1>();
            m_itr = next.m_itr;
            m_end = next.m_end;
            m_idx = next.m_idx;
        }
        // keep returning this iterator. If it reached the last iterator pair,
        // this gets compared to the correct global sentinel
        return *this;
    }

    /// @returns the single value that the iterator points to - const
    DETRAY_HOST_DEVICE
    constexpr const auto &operator*() const { return *m_itr; }

    /// @returns the single value that the iterator points to
    DETRAY_HOST_DEVICE
    constexpr auto &operator*() { return *m_itr; }

    /// Advance the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const std::size_t j) {
        chain_iterator<iterator_coll_t> tmp(*this);
        tmp.m_itr = tmp.m_itr + j;
        return tmp;
    }

    /// Global range collection
    iterator_coll_t &m_begins, m_ends;
    /// Begin and end iterator for the current range
    iterator_t m_itr;
    iterator_t &m_end;
    /// Current index
    std::size_t m_idx;

    private:
    template <std::size_t N, std::size_t I = 0>
    constexpr auto get_next_itr() {
        constexpr bool is_last = (N <= 1);

        if constexpr (!is_last) {
            if ((m_idx != I)) {
                return get_next_itr<N - 1, I + 1>();
            }
        }
        return create<I + 1>(m_begins, m_ends);
    }
};

}  // namespace detray::ranges::detail
