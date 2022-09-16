/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/detail/enumerate_iterator.hpp"

namespace detray::ranges {

/// @brief Struct that implements a subrange by providing start and end
/// iterators on a range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange
template <typename range_itr_t, typename incrementable_t>
struct enumerate_view : public detray::ranges::view_interface<
                            enumerate_view<range_itr_t, incrementable_t>> {
    /// Always use const iterator of the range as the wrapped iterator
    using iterator_type =
        ranges::detail::enumerate_iterator<range_itr_t, incrementable_t>;
    using const_iterator_type = iterator_type;

    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    enumerate_view() = default;

    template <
        typename range_t,
        std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate_view(range_t &&rng)
        : m_range{detray::ranges::begin(rng), 0},
          m_end{detray::ranges::end(rng),
                static_cast<incrementable_t>(detray::ranges::size(rng))} {}

    template <
        typename range_t,
        std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate_view(
        range_t &&rng, incrementable_t &&start)
        : m_range{detray::ranges::begin(rng), start},
          m_end{detray::ranges::end(rng),
                static_cast<incrementable_t>(detray::ranges::size(rng))} {}

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator_type & { return m_range; }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> const_iterator_type & { return m_range; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() -> const_iterator_type & { return m_end; }

    /// @note Cannot peak at the end of input-iterator based range
    constexpr typename std::iterator_traits<iterator_type>::value_type
    back() noexcept = delete;

    iterator_type m_range, m_end;
};

// deduction guides

template <typename range_t,
          std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
DETRAY_HOST_DEVICE enumerate_view(range_t &&rng)
    ->enumerate_view<detray::ranges::const_iterator_t<range_t>, dindex>;

template <typename range_t, typename incrementable_t,
          std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
DETRAY_HOST_DEVICE enumerate_view(range_t &&rng, incrementable_t &&start)
    ->enumerate_view<detray::ranges::const_iterator_t<range_t>,
                     incrementable_t>;

}  // namespace detray::ranges
