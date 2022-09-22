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

template <typename I>
struct iota_iterator;

}

/// @brief Struct that implements a subrange by providing start and end
/// iterators on a range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange
template <typename incrementable_t>
struct iota_view
    : public detray::ranges::view_interface<iota_view<incrementable_t>> {

    using iterator_t = detail::iota_iterator<incrementable_t>;

    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    iota_view() = default;

    template <typename interval_t,
              std::enable_if_t<detray::detail::is_interval_v<interval_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota_view(interval_t &&interval)
        : m_start{detray::detail::get<0>(interval)},
          m_end{detray::detail::get<1>(interval)} {}

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <typename deduced_incr_t>
    DETRAY_HOST_DEVICE constexpr iota_view(deduced_incr_t &&start,
                                           deduced_incr_t &&end)
        : m_start{start}, m_end{end - 1} {}

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <
        typename deduced_incr_t,
        std::enable_if_t<not detray::detail::is_interval_v<deduced_incr_t>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota_view(deduced_incr_t &&start)
        : m_start{start}, m_end{start + 1} {}

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator_t { return {m_start}; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t { return {m_end}; }

    /// @note Cannot peak at the end of input-iterator based range
    constexpr typename std::iterator_traits<iterator_t>::value_type
    back() noexcept = delete;

    /// Start and end values of the sequence
    incrementable_t m_start;
    incrementable_t m_end;
};

namespace views {

/// @brief interface type to construct a @c iota_view with CTAD
template <typename incrementable_t, typename sentinel_t = void>
struct iota : public detray::ranges::iota_view<incrementable_t> {

    using base_type = detray::ranges::iota_view<incrementable_t>;
    using iterator_t = typename base_type::iterator_t;

    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    iota() = default;

    template <
        typename deduced_interval_t,
        std::enable_if_t<detray::detail::is_interval_v<deduced_interval_t>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota(deduced_interval_t &&interval)
        : base_type(std::forward<deduced_interval_t>(interval)) {}

    template <typename deduced_incr_t>
    DETRAY_HOST_DEVICE constexpr iota(deduced_incr_t &&start,
                                      deduced_incr_t &&end)
        : base_type(std::forward<deduced_incr_t>(start),
                    std::forward<deduced_incr_t>(end)) {}

    template <
        typename deduced_incr_t,
        std::enable_if_t<not detray::detail::is_interval_v<deduced_incr_t>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota(deduced_incr_t &&start)
        : base_type(std::forward<deduced_incr_t>(start)) {}
};

// deduction guides

template <typename deduced_interval_t,
          std::enable_if_t<detray::detail::is_interval_v<deduced_interval_t>,
                           bool> = true>
DETRAY_HOST_DEVICE iota(deduced_interval_t &&interval)
    ->iota<std::remove_cv_t<std::remove_reference_t<
        decltype(std::get<0>(std::declval<deduced_interval_t>()))>>>;

template <typename deduced_incr_t = dindex>
DETRAY_HOST_DEVICE iota(deduced_incr_t &&start, deduced_incr_t &&end)
    ->iota<std::remove_cv_t<std::remove_reference_t<deduced_incr_t>>>;

template <typename deduced_incr_t,
          std::enable_if_t<not detray::detail::is_interval_v<deduced_incr_t>,
                           bool> = true>
DETRAY_HOST_DEVICE iota(deduced_incr_t &&start)
    ->iota<std::remove_cv_t<std::remove_reference_t<deduced_incr_t>>>;

}  // namespace views

namespace detail {

/// Nested iterator to generate a range of values on demand
template <typename value_t>
struct iota_iterator {

    /// @returns true if we reach end of sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const iota_iterator<value_t> &rhs) const -> bool {
        return (i == rhs.i);
    }

    /// @returns true if we reach end of sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(const iota_iterator<value_t> &rhs) const -> bool {
        return (i != rhs.i);
    }

    /// Increment
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> iota_iterator<value_t> & {
        ++i;
        return *this;
    }

    /// @returns the current value in the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const -> const value_t & { return i; }

    /// @returns the current value in the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator*() -> value_t & { return i; }

    /// Advance the sequence by @param j positions
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const value_t j) const -> iota_iterator<value_t> {
        return {i + j};
    }

    /// Current value of sequence
    value_t i;
};

}  // namespace detail

}  // namespace detray::ranges

namespace std {

/// Specialization of std::iterator_traits struct for the iota range factory
template <typename T>
struct iterator_traits<detray::ranges::detail::iota_iterator<T>> {
    using difference_type = T;
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using iterator_category = std::input_iterator_tag;
};

}  // namespace std
