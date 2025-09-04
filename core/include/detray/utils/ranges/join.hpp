/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::ranges {

namespace detail {

template <typename T>
struct join_iterator;

}

/// @brief Range adaptor that joins different ranges of the same type (static)
///
/// @see https://en.cppreference.com/w/cpp/ranges/join_view
///
/// @tparam range_t a range of the ranges that should be joined.
///
/// @note Static implementation: The number of ranges needs to be know at
/// compile time
/// @note Does not take ownership of the ranges it operates on. Their lifetime
/// needs to be guranteed throughout iteration or between iterations with the
/// same join instance.
template <typename range_t>
struct join_view : public detray::ranges::view_interface<join_view<range_t>> {

    /// Iterator over the range of ranges
    using outer_iterator_t = detray::ranges::iterator_t<range_t>;
    // Type of a single range
    using outer_value_t =
        std::conditional_t<std::is_const<range_t>::value ,
                           const detray::ranges::range_value_t<range_t>,
                           detray::ranges::range_value_t<range_t>>;
    // Iterator over a single range
    using inner_iterator_t = detray::ranges::iterator_t<outer_value_t>;

    using iterator_t = detray::ranges::detail::join_iterator<range_t>;
    using value_t = typename std::iterator_traits<iterator_t>::value_type;

    /// Default constructor
    constexpr join_view() = default;

    /// Construct from a range of @param ranges.
    template <typename R>
    DETRAY_HOST_DEVICE constexpr explicit join_view(R &&ranges)
        : m_begin{detray::ranges::begin(std::forward<R>(ranges))},
          m_end{detray::ranges::end(std::forward<R>(ranges))} {}

    /// Copy constructor
    DETRAY_HOST_DEVICE
    constexpr join_view(const join_view &other)
        : m_begin{other.m_begin}, m_end{other.m_end} {}

    /// Default destructor
    ~join_view() = default;

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    join_view &operator=(const join_view &other) {
        m_begin = other.m_begin;
        m_end = other.m_end;
        return *this;
    }

    /// @return start position of range - const
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator_t { return {m_begin, m_end}; }

    /// @return sentinel of the range.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t { return {m_end, m_end}; }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - const
    DETRAY_HOST_DEVICE
    constexpr auto data() const -> const value_t * { return &(*(*m_begin)); }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - non-const
    DETRAY_HOST_DEVICE
    constexpr auto data() -> value_t * { return &(*(*m_begin)); }

    /// @returns sum of the number elements of all ranges in the join
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> std::size_t {
        std::size_t size{0u};
        for (outer_iterator_t itr = m_begin; itr != m_end; ++itr) {
            // subrange
            const auto begin = detray::ranges::begin(*itr);
            const auto end = detray::ranges::end(*itr);
            size +=
                static_cast<std::size_t>(detray::ranges::distance(begin, end));
        }
        return size;
    }

    /// Start and end position of the subranges
    outer_iterator_t m_begin{}, m_end{};
};

namespace views {

/// @brief interface type to construct a @c join_view with CTAD
template <typename range_t>
struct join : public ranges::join_view<range_t> {

    using base_type = ranges::join_view<range_t>;

    constexpr join() = default;

    template <typename deduced_range_t>
    DETRAY_HOST_DEVICE constexpr explicit join(deduced_range_t &&ranges)
        : base_type(std::forward<deduced_range_t>(ranges)) {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    join &operator=(const join &other) {
        base_type::operator=(other);
        return *this;
    }
};

// deduction guides
#ifndef DETRAY_COMPILE_VITIS
template <typename R>
DETRAY_HOST_DEVICE join(R &&ranges)->join<std::remove_reference_t<R>>;
#endif

}  // namespace views


namespace detail {

/// @brief Sequentially iterate through multiple ranges of the same type.
///
/// Once the sentinel of one range is reached, set the current iterator to the
/// next ranges 'begin' (or 'end' if decrementing)
///
/// @tparam range_t a range that contains the ranges to be joined.
template <typename range_t>
struct join_iterator {

    using outer_iterator_t =
        std::conditional_t<std::is_const<range_t>::value ,
                           detray::ranges::const_iterator_t<range_t>,
                           detray::ranges::iterator_t<range_t>>;
    using outer_value_t = decltype(*std::declval<outer_iterator_t>());
    using inner_iterator_t = std::conditional_t<
        std::is_same<outer_iterator_t,
                       detray::ranges::const_iterator_t<range_t>>::value ,
        detray::ranges::const_iterator_t<outer_value_t>,
        detray::ranges::iterator_t<outer_value_t>>;

    using iterator_t = inner_iterator_t;
    using difference_type =
        typename std::iterator_traits<iterator_t>::difference_type;
    using value_type = typename std::iterator_traits<iterator_t>::value_type;
    using pointer = typename std::iterator_traits<iterator_t>::pointer;
    using reference = typename std::iterator_traits<iterator_t>::reference;
    using iterator_category =
        typename std::iterator_traits<iterator_t>::iterator_category;

    /// Default constructor required by LegacyIterator trait
    constexpr join_iterator() = default;

    /// Construct from range of ranges ( @param begin and @param  end )
    DETRAY_HOST_DEVICE
    constexpr join_iterator(outer_iterator_t begin, outer_iterator_t end)
        : m_outer_begin(begin), m_outer_end(end), m_outer_itr(begin) {

        if (m_outer_itr != m_outer_end) {
            m_inner_itr = (*m_outer_itr).begin();
        } else {
            m_inner_itr = {};
        }
        next_inner();
    }

    /// @returns true if it points to the same value.
    DETRAY_HOST_DEVICE constexpr bool operator==(
        const join_iterator &rhs) const {
        return (m_inner_itr == rhs.m_inner_itr);
    }

    /// @returns false if it points to the same value (usually the global
    /// sentinel of the join).
    DETRAY_HOST_DEVICE constexpr bool operator!=(
        const join_iterator &rhs) const {
        return (m_outer_itr != rhs.m_outer_itr);
    }

    /// Increment current iterator and check for switch between ranges.
    DETRAY_HOST_DEVICE auto operator++() -> join_iterator & {
        ++m_inner_itr;
        next_inner();

        return *this;
    }

    /// Decrement current iterator and check for switch between ranges.
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::bidirectional_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator--() -> join_iterator & {
        // If we are calling this at the end of the join iteration, go back into
        // the valid range
        if (m_outer_itr == m_outer_end) {
            --m_outer_itr;
        }

        previous_inner();
        --m_inner_itr;

        return *this;
    }

    /// @returns the single value that the iterator points to - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator*() const { return *m_inner_itr; }

    /// @returns an iterator advanced by @param j through the join.
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator+(const difference_type j) const
        -> join_iterator {
        join_iterator<range_t> tmp(*this);
        // walk through join to catch the switch between intermediate ranges
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

    /// @returns an iterator advanced by @param j through the join.
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator-(const difference_type j) const
        -> join_iterator {
        return *this + (-j);
    }

    /// @returns the positional difference between two iterators
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator-(
        const join_iterator &other) const -> difference_type {
        outer_iterator_t tmp_outer_itr;
        if (tmp_outer_itr == other.m_outer_itr) {
            return m_inner_itr - other.m_inner_itr;
        }
        if ((tmp_outer_itr - other.m_outer_itr) < 0) {
            // Negative distance
            difference_type diff{m_inner_itr -
                                 detray::ranges::end(*m_outer_itr)};
            for (tmp_outer_itr + 1; tmp_outer_itr != other.m_outer_itr;
                 ++tmp_outer_itr) {
                diff += detray::ranges::end(*tmp_outer_itr) -
                        detray::ranges::begin(*tmp_outer_itr);
            }
            diff += other.m_inner_itr - detray::ranges::end(*other.m_outer_itr);
            return diff;
        } else {
            // Positive distance
            difference_type diff{m_inner_itr -
                                 detray::ranges::begin(*tmp_outer_itr)};
            for (tmp_outer_itr - 1; tmp_outer_itr != other.m_outer_itr;
                 --tmp_outer_itr) {
                diff += detray::ranges::end(*tmp_outer_itr) -
                        detray::ranges::begin(*tmp_outer_itr);
            }
            diff +=
                other.m_inner_itr - detray::ranges::begin(*other.m_outer_itr);
            return diff;
        }
    }

    /// @returns advance this iterator state by @param j.
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator+=(const difference_type j)
        -> join_iterator & {
        // walk through join to catch the switch between intermediate ranges
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
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator-=(const difference_type j)
        -> join_iterator & {
        m_inner_itr += (-j);
        return *this;
    }

    /// @returns the value at a given position - const
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator[](const difference_type i) const
        -> const value_type & {
        difference_type offset{
            i - (m_inner_itr - detray::ranges::begin(*m_outer_itr))};
        return *(*this + offset);
    }

    /// @returns the value at a given position - const
    template <typename I = iterator_t,
              std::enable_if_t<detray::ranges::random_access_iterator_v<I>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto operator[](const difference_type i) {
        difference_type offset{
            i - (m_inner_itr - detray::ranges::begin(*m_outer_itr))};
        return *(*this + offset);
    }

    private:
    /// Find the first inner range that is not empty
    constexpr void next_inner() {
        if (m_outer_itr == m_outer_end) {
            return;
        }
        while (m_inner_itr == detray::ranges::end(*m_outer_itr)) {
            ++m_outer_itr;
            if (m_outer_itr != m_outer_end) {
                m_inner_itr = (*m_outer_itr).begin();
            } else {
                break;
            }
        }
    }

    /// Find the last inner range that is not empty
    constexpr void previous_inner() {
        // Get the start of the current inner range
        inner_iterator_t inner_begin = detray::ranges::begin(*m_outer_itr);

        // Iterator has reached last valid position in this range
        // during the previous decrement. Now go to the end of the
        // previous range
        while (m_inner_itr == inner_begin) {
            // No more inner ranges to try
            if (m_outer_itr == m_outer_begin) {
                m_inner_itr = detray::ranges::end(*m_outer_begin);
                return;
            }

            --m_outer_itr;

            inner_begin = detray::ranges::begin(*m_outer_itr);
            m_inner_itr = detray::ranges::end(*m_outer_itr);
        }
    }

    /// Global range collection begin and end (outer iterators)
    outer_iterator_t m_outer_begin{}, m_outer_end{};
    /// Current range
    outer_iterator_t m_outer_itr{};
    /// Current iterators over the inner ranges
    inner_iterator_t m_inner_itr{};
};

}  // namespace detail

}  // namespace detray::ranges
