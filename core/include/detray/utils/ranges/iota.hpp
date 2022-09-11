/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <type_traits>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray::ranges {

/// @brief Struct that implements a subrange by providing start and end
/// iterators on a range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange
template <typename incrementable_t, typename sentinel_t = void>
struct iota_view
    : public detray::ranges::view_interface<iota_view<incrementable_t>> {

    /// Nested range
    struct iterator {
        /// Start and end of sequence
        incrementable_t i;

        /// Determine whether we reach end of sequence
        DETRAY_HOST_DEVICE
        bool operator!=(const iterator &rhs) const { return i != rhs.i; }

        /// Increase index and iterator at once
        DETRAY_HOST_DEVICE
        void operator++() { ++i; }

        /// Tie them together for returning
        DETRAY_HOST_DEVICE
        auto operator*() const { return i; }
    };

    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    iota_view() = default;

    template <
        typename interval_t,
        typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
            decltype(std::get<0>(std::declval<interval_t>()))>>>,
        typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
            decltype(std::get<1>(std::declval<interval_t>()))>>>>
    DETRAY_HOST_DEVICE iota_view(interval_t &&interval)
        : m_start{detray::detail::get<0>(interval)},
          m_end{detray::detail::get<1>(interval)} {}

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <typename deduced_incr_t, typename deduced_sentinel_t>
    DETRAY_HOST_DEVICE constexpr iota_view(deduced_incr_t &&start,
                                           deduced_sentinel_t &&end)
        : m_start{start}, m_end{end} {}

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <typename deduced_incr_t>
    DETRAY_HOST_DEVICE constexpr explicit iota_view(deduced_incr_t &&start)
        : m_start{start}, m_end{0} {}

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    inline constexpr auto begin() -> iterator { return iterator{m_start}; }

    /// @return end position of range on co jntainer.
    template <
        std::enable_if_t<not std::is_same_v<sentinel_t, void>, bool> = true>
    DETRAY_HOST_DEVICE inline constexpr auto end() -> iterator {
        return iterator{m_end};
    }

    /// @return end position of range on container.
    template <std::enable_if_t<std::is_same_v<sentinel_t, void>, bool> = true>
    DETRAY_HOST_DEVICE inline constexpr auto end() -> iterator {
        return iterator{m_start + 1};
    }

    /// Start and end position of the subrange
    // TODO: std::conditional for array size
    incrementable_t m_start;
    std::conditional_t<std::is_same_v<sentinel_t, void>, incrementable_t,
                       sentinel_t>
        m_end;
};

// deduction guides

template <
    typename interval_t,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<0>(std::declval<interval_t>()))>>>,
    typename = std::enable_if_t<std::is_integral_v<std::remove_reference_t<
        decltype(std::get<1>(std::declval<interval_t>()))>>>>
DETRAY_HOST_DEVICE iota_view(interval_t &&interval)
    ->iota_view<std::remove_reference_t<
        decltype(std::get<0>(std::declval<interval_t>()))>>;

template <typename deduced_incr_t, typename deduced_sentinel_t>
DETRAY_HOST_DEVICE iota_view(deduced_incr_t &&start, deduced_sentinel_t &&end)
    ->iota_view<std::remove_reference_t<deduced_incr_t>,
                std::remove_reference_t<deduced_sentinel_t>>;

template <typename deduced_incr_t>
DETRAY_HOST_DEVICE iota_view(deduced_incr_t &&start)
    ->iota_view<std::remove_reference_t<deduced_incr_t>>;

}  // namespace detray::ranges
