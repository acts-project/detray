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
#include "detray/utils/ranges/detail/iota_iterator.hpp"
#include "detray/utils/type_traits.hpp"

namespace detray::ranges {

/// @brief Struct that implements a subrange by providing start and end
/// iterators on a range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange
template <typename incrementable_t, typename sentinel_t = void>
struct iota_view : public detray::ranges::view_interface<
                       iota_view<incrementable_t, sentinel_t>> {

    using iterator_type = detail::iota_iterator<incrementable_t>;
    using sentinel_type = detail::iota_iterator<std::conditional_t<
        std::is_same_v<sentinel_t, void>, incrementable_t, sentinel_t>>;

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
    template <typename deduced_incr_t, typename deduced_sentinel_t>
    DETRAY_HOST_DEVICE constexpr iota_view(deduced_incr_t &&start,
                                           deduced_sentinel_t &&end)
        : m_start{start}, m_end{end} {}

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
    constexpr auto begin() const -> iterator_type { return m_start; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE constexpr auto end() const -> sentinel_type {
        return m_end;
    }

    /// Start and end position of the subrange
    // TODO: std::conditional for array size
    iterator_type m_start;
    sentinel_type m_end;
};

// deduction guides

template <
    typename interval_t,
    std::enable_if_t<detray::detail::is_interval_v<interval_t>, bool> = true>
DETRAY_HOST_DEVICE iota_view(interval_t &&interval)
    ->iota_view<std::remove_reference_t<
                    decltype(std::get<0>(std::declval<interval_t>()))>,
                std::remove_reference_t<
                    decltype(std::get<1>(std::declval<interval_t>()))>>;

template <typename deduced_incr_t, typename deduced_sentinel_t>
DETRAY_HOST_DEVICE iota_view(deduced_incr_t &&start, deduced_sentinel_t &&end)
    ->iota_view<std::remove_reference_t<deduced_incr_t>,
                std::remove_reference_t<deduced_sentinel_t>>;

template <typename deduced_incr_t,
          std::enable_if_t<not detray::detail::is_interval_v<deduced_incr_t>,
                           bool> = true>
DETRAY_HOST_DEVICE iota_view(deduced_incr_t &&start)
    ->iota_view<std::remove_reference_t<deduced_incr_t>, void>;

}  // namespace detray::ranges
