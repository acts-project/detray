/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <iterator>
#include <type_traits>

namespace detray::ranges {

namespace detail {

template <typename T>
struct chain_iterator;

}

/// @brief Implements a subrange by providing start and end iterators on
/// another range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange.
template <std::size_t I, typename range_itr_t>
struct chain_view : public ranges::view_interface<chain_view<I, range_itr_t>> {

    using iterator_coll_t = std::array<range_itr_t, I>;
    using iterator_t = detray::ranges::detail::chain_iterator<iterator_coll_t>;
    using const_iterator_t =
        detray::ranges::detail::chain_iterator<const iterator_coll_t>;

    /// Default constructor
    constexpr chain_view() = default;

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <typename... ranges_t>
    DETRAY_HOST_DEVICE constexpr explicit chain_view(ranges_t &&...ranges)
        : m_begins{detray::ranges::begin(ranges)...},
          m_ends{detray::ranges::end(ranges)...} {}

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <typename... ranges_t>
    DETRAY_HOST_DEVICE constexpr explicit chain_view(const ranges_t &...ranges)
        : m_begins{detray::ranges::cbegin(ranges)...},
          m_ends{detray::ranges::cend(ranges)...} {}

    /// @return start position of range.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator_t {
        return iterator_t(m_begins, m_ends);
    }

    /// @return start position of range.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> const_iterator_t {
        return iterator_t(m_begins, m_ends);
    }

    /// @return sentinel of the range.
    DETRAY_HOST_DEVICE
    constexpr auto end() {
        // Build a chained itr from the last value in the iterator collection
        auto end = iterator_t::template create<I - 1>(m_begins, m_ends);
        end.m_itr = end.m_end;
        return end;
    }

    /// @note For now no 'size()' function
    constexpr typename std::iterator_traits<range_itr_t>::difference_type size()
        const noexcept = delete;

    /// @note For now no 'back()' function
    constexpr typename std::iterator_traits<range_itr_t>::value_type back()
        const noexcept = delete;

    /// Start and end position of the subranges
    iterator_coll_t m_begins, m_ends;
};

namespace views {

/// @brief interface type to construct a @c chain_view with CTAD
template <std::size_t I, typename range_itr_t>
struct chain : public ranges::chain_view<I, range_itr_t> {

    using base_type = ranges::chain_view<I, range_itr_t>;
    constexpr chain() = default;

    template <typename... ranges_t>
    DETRAY_HOST_DEVICE constexpr explicit chain(const ranges_t &...ranges)
        : base_type(ranges...) {}

    template <typename... ranges_t>
    DETRAY_HOST_DEVICE constexpr explicit chain(ranges_t &&...ranges)
        : base_type(std::forward<ranges_t>(ranges)...) {}
};

// deduction guides

template <typename... ranges_t>
DETRAY_HOST_DEVICE chain(const ranges_t &...ranges)
    ->chain<sizeof...(ranges_t),
            typename detray::ranges::const_iterator_t<
                detray::detail::first_t<ranges_t...>>>;

template <typename... ranges_t>
DETRAY_HOST_DEVICE chain(ranges_t &&...ranges)
    ->chain<sizeof...(ranges_t),
            typename detray::ranges::iterator_t<
                detray::detail::first_t<ranges_t...>>>;

}  // namespace views

namespace detail {

/// @brief Sequentially iterate through multiple ranges of the same type.
///
/// Once the sentinel of one range is reached, set the current iterator to the
/// next ranges 'begin' and update the sentinel.
template <typename iterator_coll_t>
struct chain_iterator {

    using iterator_t = detray::detail::get_value_type_t<iterator_coll_t>;

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

        // Time to get the next iterator pair
        constexpr std::size_t max_idx =
            sizeof(iterator_coll_t) / sizeof(iterator_t) - 1;
        if (m_itr == m_end and m_idx < max_idx) {
            // Copy the temporay into the current iterator, which is the only
            // one that is kept being used in the for loop
            auto next = get_next_itr<max_idx>();
            m_itr = next.m_itr;
            m_end = next.m_end;
            m_idx = next.m_idx;
        }
        // keep returning this iterator unchanged
        return *this;
    }

    /// @returns the single value that the iterator points to - const
    DETRAY_HOST_DEVICE
    constexpr const auto &operator*() const { return *m_itr; }

    /// @returns the single value that the iterator points to
    DETRAY_HOST_DEVICE
    constexpr auto &operator*() { return *m_itr; }

    /// Advance the iterator of the current range
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const std::size_t j) {
        chain_iterator<iterator_coll_t> tmp(*this);
        tmp.m_itr = tmp.m_itr + j;
        return tmp;
    }

    /// Global range collection of begin and end iterators
    iterator_coll_t &m_begins, &m_ends;
    /// This is the actual iterator state that will be advanced during iteration
    iterator_t m_itr;
    /// Compare against the corresponding sentinel of the wrapped iterator
    iterator_t &m_end;
    /// Index of the current range in the chain
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

}  // namespace detail

}  // namespace detray::ranges

/// Add the iterator trait specialization for the chain iterator
namespace std {

/// Specialization of std::iterator_traits struct for the sequential iteration
/// of multiple ranges
template <typename T>
struct iterator_traits<detray::ranges::detail::chain_iterator<T>> {
    private:
    using iterator_t = detray::detail::get_value_type_t<T>;

    public:
    using difference_type =
        typename std::iterator_traits<iterator_t>::difference_type;
    using value_type = typename std::iterator_traits<iterator_t>::value_type;
    using pointer = typename std::iterator_traits<iterator_t>::pointer;
    using reference = typename std::iterator_traits<iterator_t>::reference;
    using iterator_category = std::forward_iterator_tag;
};

}  // namespace std
