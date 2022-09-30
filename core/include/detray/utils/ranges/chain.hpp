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
#include <type_traits>

namespace detray::ranges {

namespace detail {

template <typename T>
struct chain_iterator;

}

/// @brief Chain together different ranges of the same type.
///
/// @tparam I the number of ranges in the chain.
/// @tparam range_itr_t the iterator type of the ranges.
///
/// @note Does not take ownership of the ranges it operates on. Their lifetime
/// needs to be guranteed throughout iteration or between iterations with the
/// same chain instance.
/// @note Is not fit for lazy evaluation.
template <std::size_t I, typename range_itr_t>
struct chain_view
    : public detray::ranges::view_interface<chain_view<I, range_itr_t>> {

    using iterator_coll_t = std::array<range_itr_t, I>;
    using iterator_t = detray::ranges::detail::chain_iterator<iterator_coll_t>;
    using value_t = typename std::iterator_traits<iterator_t>::value_type;

    /// Default constructor
    constexpr chain_view() = default;

    /// Construct from a pack of @param ranges.
    template <typename... ranges_t>
    DETRAY_HOST_DEVICE constexpr explicit chain_view(ranges_t &&...ranges)
        : m_begins{detray::ranges::begin(ranges)...},
          m_ends{detray::ranges::end(ranges)...} {}

    /// Construct from a pack of @param ranges - const
    template <typename... ranges_t>
    DETRAY_HOST_DEVICE constexpr explicit chain_view(const ranges_t &...ranges)
        : m_begins{detray::ranges::cbegin(ranges)...},
          m_ends{detray::ranges::cend(ranges)...} {}

    /// Copy assignment operator (otherwise implicitely deleted due to reference
    /// members)
    DETRAY_HOST_DEVICE
    chain_view &operator=(const chain_view &other) {
        m_begins = other.m_begins;
        m_ends = other.m_ends;

        return *this;
    }

    /// @return start position of range - const
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator_t { return {m_begins, m_ends}; }

    /// @return sentinel of the range.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t {
        // Build a chained itr from the last value in the iterator collection
        return {m_begins, m_ends, detray::detail::get<I - 1>(m_ends), I - 1};
    }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - const
    DETRAY_HOST_DEVICE
    constexpr auto data() const noexcept -> const value_t * {
        return &(*(detray::detail::get<0>(m_begins())));
    }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - non-const
    DETRAY_HOST_DEVICE
    constexpr auto data() noexcept -> value_t * {
        return &(*(detray::detail::get<0>(m_begins())));
    }

    /// @returns sum of the number elements of all ranges in the chain
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept ->
        typename std::iterator_traits<range_itr_t>::difference_type {
        std::size_t size{0};
        for (std::size_t i{0}; i < I; ++i) {
            const range_itr_t &begin = m_begins[i];
            const range_itr_t &end = m_ends[i];
            size += detray::ranges::distance(begin, end);
        }
        return size;
    }

    /// Start and end position of the subranges
    iterator_coll_t m_begins{}, m_ends{};
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
///
/// @tparam iterator_coll_t type of iterator collection of ranges to be chained.
///         Can contain const iterators.
///
/// @note The iterator must not be typed on the current range index, so that
/// begin and sentinel type are the same.
/// @todo Add Comparability to fulfill random access iterator traits once
///       needed.
template <typename iterator_coll_t>
struct chain_iterator {

    using iterator_t = detray::detail::get_value_type_t<iterator_coll_t>;

    using difference_type =
        typename std::iterator_traits<iterator_t>::difference_type;
    using value_type = typename std::iterator_traits<iterator_t>::value_type;
    using pointer = typename std::iterator_traits<iterator_t>::pointer;
    using reference = typename std::iterator_traits<iterator_t>::reference;
    using iterator_category =
        typename std::iterator_traits<iterator_t>::iterator_category;

    /// Default constructor required by LegacyIterator trait
    constexpr chain_iterator() = default;

    /// Construct from a collection of @param begin and @param  end positions
    DETRAY_HOST_DEVICE
    constexpr chain_iterator(const iterator_coll_t &begins,
                             const iterator_coll_t &ends)
        : m_begins(&begins),
          m_ends(&ends),
          m_iter{detray::detail::get<0>(*m_begins)},
          m_idx{0} {}

    /// Fully parametrized construction
    DETRAY_HOST_DEVICE
    constexpr chain_iterator(const iterator_coll_t &begins,
                             const iterator_coll_t &ends, iterator_t current,
                             const std::size_t i)
        : m_begins(&begins), m_ends(&ends), m_iter{current}, m_idx{i} {}

    /// @returns true if it points to the same value.
    template <typename T = iterator_category,
              std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, T>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr bool operator==(
        const chain_iterator &rhs) const {
        return m_iter == rhs.m_iter;
    }

    /// @returns false if it points to the same value (usually the global
    /// sentinel of the chain).
    template <typename T = iterator_category,
              std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, T>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr bool operator!=(
        const chain_iterator &rhs) const {
        return m_iter != rhs.m_iter;
    }

    /// Increment current iterator and check for switch between ranges.
    template <typename T = iterator_category,
              std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, T>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator++() -> chain_iterator & {
        ++m_iter;

        // Switch to next range in the collection
        constexpr std::size_t max_idx =
            sizeof(iterator_coll_t) / sizeof(iterator_t) - 1;
        if (m_iter == (*m_ends)[m_idx] and m_idx < max_idx) {
            ++m_idx;
            m_iter = (*m_begins)[m_idx];
        }
        return *this;
    }

    /// Decrement current iterator and check for switch between ranges.
    template <
        typename T = iterator_category,
        std::enable_if_t<std::is_base_of_v<std::bidirectional_iterator_tag, T>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator--() -> chain_iterator & {

        if (m_iter != (*m_begins)[m_idx] and m_idx > 0) {
            // Normal case
            --m_idx;
            --m_iter;
        } else if (m_idx > 0) {
            // Iterator has reached last valid position in this range during the
            // previous decrement. Now go to the end of the next range
            --m_idx;
            m_iter = (*m_ends)[m_idx] - difference_type{1};
        }
        return *this;
    }

    /// @returns the single value that the iterator points to - const
    template <typename T = iterator_category,
              std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, T>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator*() const -> const value_type & {
        return *m_iter;
    }

    /// @returns the single value that the iterator points to.
    template <typename T = iterator_category,
              std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, T>,
                               bool> = true>
    DETRAY_HOST_DEVICE auto operator*() -> value_type & {
        return *m_iter;
    }

    /// @returns an iterator advanced by @param j through the chain.
    template <
        typename T = iterator_category,
        std::enable_if_t<std::is_base_of_v<std::random_access_iterator_tag, T>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator+(const difference_type j) const
        -> chain_iterator {
        chain_iterator<iterator_coll_t> tmp(*this);
        // walk through chain to catch the switch between intermediate ranges
        difference_type i{j};
        if (i >= difference_type{0}) {
            while (i--) {
                ++tmp;
            };
        } else {
            while (i++) {
                --tmp;
            };
        }
        return tmp;
    }

    /// @returns an iterator advanced by @param j through the chain.
    template <
        typename T = iterator_category,
        std::enable_if_t<std::is_base_of_v<std::random_access_iterator_tag, T>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator-(const difference_type j) const
        -> chain_iterator {
        chain_iterator<iterator_coll_t> tmp(*this);

        return tmp + (-j);
    }

    /// @returns the positional difference between two iterators
    template <
        typename T = iterator_category,
        std::enable_if_t<std::is_base_of_v<std::random_access_iterator_tag, T>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator-(
        const chain_iterator &other) const -> difference_type {
        if (m_idx == other.m_idx) {
            return m_iter - other.m_iter;
        }
        if (m_idx < other.m_idx) {
            // Negative distance
            difference_type diff{m_iter - (*m_ends)[m_idx]};
            for (std::size_t i{m_idx + 1}; i < other.m_idx; ++i) {
                diff += (*m_begins)[i] - (*m_ends)[i];
            }
            diff += (*other.m_begins)[m_idx] - other.m_iter;

            return diff;
        } else {
            // Positive distance
            difference_type diff{m_iter - (*m_begins)[m_idx]};
            for (std::size_t i{m_idx - 1}; i > other.m_idx; --i) {
                diff += (*m_ends)[i] - (*m_begins)[i];
            }
            diff += (*other.m_ends)[other.m_idx] - other.m_iter;
            return diff;
        }
    }

    /// @returns advance this iterator state by @param j.
    template <
        typename T = iterator_category,
        std::enable_if_t<std::is_base_of_v<std::random_access_iterator_tag, T>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator+=(const difference_type j)
        -> chain_iterator & {
        // walk through chain to catch the switch between intermediate ranges
        difference_type i{j};
        if (i >= difference_type{0}) {
            while (i--) {
                ++(*this);
            };
        } else {
            while (i++) {
                --(*this);
            };
        }

        return *this;
    }

    /// @returns advance this iterator state by @param j.
    template <
        typename T = iterator_category,
        std::enable_if_t<std::is_base_of_v<std::random_access_iterator_tag, T>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator-=(const difference_type j)
        -> chain_iterator & {
        m_iter += (-j);
        return *this;
    }

    /// Global range collection of begin and end iterators
    const iterator_coll_t *m_begins{nullptr}, *m_ends{nullptr};
    /// This is the actual iterator state that will be advanced during iteration
    iterator_t m_iter{};
    /// Index of the current range in the chain
    std::size_t m_idx{0};
};

}  // namespace detail

}  // namespace detray::ranges