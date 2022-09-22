/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <tuple>
#include <type_traits>

#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"

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

namespace detail {

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
