/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/detail/chain_iterator.hpp"
#include "detray/utils/ranges/ranges.hpp"

namespace detray::ranges {

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
    constexpr auto begin() const -> iterator_t {
        return iterator_t(m_begins, m_ends);
    }

    /// @return sentinel of the range.
    DETRAY_HOST_DEVICE
    constexpr auto end() {
        // Build a chained itr from the last value in the iterator collection
        detray::ranges::detail::chain_iterator<iterator_coll_t> end(
            iterator_t::template create<I - 1>(m_begins, m_ends));
        end.m_itr = end.m_end;
        return end;
    }

    /// @note For now no 'size()' function
    constexpr typename std::iterator_traits<range_itr_t>::difference_type size()
        const noexcept = delete;

    /// @note For now no 'size()' function
    constexpr typename std::iterator_traits<range_itr_t>::value_type back()
        const noexcept = delete;

    /// Start and end position of the subranges
    iterator_coll_t m_begins, m_ends;
};

namespace views {

/// @brief interface type to construct a @c chain_view with CTAD
template <std::size_t I, typename range_t>
struct chain : public ranges::chain_view<I, range_t> {

    using base_type = ranges::chain_view<I, range_t>;
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

}  // namespace detray::ranges
