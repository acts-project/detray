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
#include "detray/utils/ranges/ranges.hpp"

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
    using iterator_t =
        ranges::detail::enumerate_iterator<range_itr_t, incrementable_t>;
    using const_iterator_t = iterator_t;

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
    constexpr auto begin() -> iterator_t & { return m_range; }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> const_iterator_t & { return m_range; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() -> const_iterator_t & { return m_end; }

    /// @note Cannot peak at the end of input-iterator based range
    constexpr typename std::iterator_traits<iterator_t>::value_type
    back() noexcept = delete;

    iterator_t m_range, m_end;
};

namespace views {

template <typename range_itr_t, typename incrementable_t>
struct enumerate : public enumerate_view<range_itr_t, incrementable_t> {

    using base_type = enumerate_view<range_itr_t, incrementable_t>;
    using iterator_t = typename base_type::iterator_t;
    using const_iterator_t = typename base_type::const_iterator_t;

    enumerate() = default;

    template <
        typename range_t,
        std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate(range_t &&rng)
        : base_type(std::forward<range_t>(rng)) {}

    template <
        typename range_t,
        std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate(range_t &&rng,
                                                    incrementable_t &&start)
        : base_type(std::forward<range_t>(rng), std::forward<range_t>(start)) {}
};

// deduction guides

template <typename range_t,
          std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
DETRAY_HOST_DEVICE enumerate(range_t &&rng)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, dindex>;

template <typename range_t, typename incrementable_t,
          std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
DETRAY_HOST_DEVICE enumerate(range_t &&rng, incrementable_t &&start)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, incrementable_t>;

}  // namespace views

}  // namespace detray::ranges
