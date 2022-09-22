/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/ranges/subrange.hpp"

// System include(s)
#include <iterator>
#include <tuple>
#include <type_traits>

namespace detray::ranges {

namespace detail {

template <typename U, typename V>
struct enumerate_iterator;

}

/// @brief Struct that implements a subrange by providing start and end
/// iterators on a range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange
template <typename range_itr_t, typename incrementable_t = dindex>
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
        : m_range{detray::ranges::begin(std::forward<range_t>(rng)), 0},
          m_end{detray::ranges::end(std::forward<range_t>(rng)), 0} {}

    template <
        typename range_t,
        std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate_view(range_t &&rng,
                                                         dindex start)
        : m_range{detray::ranges::begin(std::forward<range_t>(rng)), start},
          m_end{detray::ranges::end(std::forward<range_t>(rng)), 0} {}

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

template <typename range_itr_t, typename incrementable_t = dindex>
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
    DETRAY_HOST_DEVICE constexpr enumerate(range_t &&rng, dindex start)
        : base_type(std::forward<range_t>(rng), start) {}

    /// Construct from a range and start/end positions
    ///
    /// @param range container to iterate over
    /// @param vol start and end position for iteration according to the volume
    template <typename deduced_range_t, typename volume_t,
              typename = typename std::remove_reference_t<volume_t>::volume_def>
    DETRAY_HOST_DEVICE enumerate(deduced_range_t &&range, const volume_t &vol)
        : enumerate(
              detray::ranges::subrange(std::forward<deduced_range_t>(range),
                                       vol),
              detray::detail::get<0>(
                  vol.template range<typename detray::ranges::range_value_t<
                      deduced_range_t>>())) {}
};

// deduction guides

template <typename range_t,
          std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
DETRAY_HOST_DEVICE enumerate(range_t &&rng)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, dindex>;

template <typename range_t, typename volume_t,
          typename = typename std::remove_reference_t<volume_t>::volume_def>
DETRAY_HOST_DEVICE enumerate(range_t &&range, const volume_t &vol)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, dindex>;

template <typename range_t,
          std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
DETRAY_HOST_DEVICE enumerate(range_t &&rng, dindex start)
    ->enumerate<detray::ranges::const_iterator_t<range_t>, dindex>;

}  // namespace views

namespace detail {

/// @brief Nested iterator to enumerate the elements of a range.
///
/// The enumeration is done by incrementing an index in lockstep with a wrapped
/// iterator of the range. Index and current iterator value are returned
/// using structured binding.
template <typename iterator_t, typename index_t>
struct enumerate_iterator {

    /// Determine whether we reach end of range
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(const enumerate_iterator &rhs) const -> bool {
        return (m_iter != rhs.m_iter);
    }

    /// Increment
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> enumerate_iterator<iterator_t, index_t> & {
        ++m_i;
        ++m_iter;
        return *this;
    }

    /// Tie them together for returning
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const { return std::tie(m_i, *m_iter); }

    /// Advance the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const index_t j) const
        -> enumerate_iterator<iterator_t, index_t> {
        return {m_iter + j, m_i + j};
    }

    /// Start value of index sequence
    iterator_t m_iter;
    index_t m_i;
};

}  // namespace detail

}  // namespace detray::ranges

namespace std {

/// Specialization of std::iterator_traits struct for the enumrator
template <typename T, typename I>
struct iterator_traits<detray::ranges::detail::enumerate_iterator<T, I>> {
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using value_type = typename std::iterator_traits<T>::value_type;
    using pointer = typename std::iterator_traits<T>::pointer;
    using reference = typename std::iterator_traits<T>::reference;
    using iterator_category = std::input_iterator_tag;
};

}  // namespace std
