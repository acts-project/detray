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
#include "detray/utils/ranges/ranges.hpp"

namespace detray::ranges {

/// @brief Implements a subrange by providing start and end iterators on
/// another range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange.
template <typename range_t>
struct subrange : public ranges::view_interface<subrange<range_t>> {

    using iterator_t = typename detray::ranges::iterator_t<range_t>;
    using const_iterator_t = typename detray::ranges::const_iterator_t<range_t>;
    using range_size_t = typename detray::ranges::range_size_t<range_t>;

    /// Default constructor
    constexpr subrange() = default;

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <typename deduced_range_t>
    DETRAY_HOST_DEVICE constexpr explicit subrange(deduced_range_t &&range)
        : m_start{detray::ranges::begin(std::forward<deduced_range_t>(range))},
          m_end{detray::ranges::end(std::forward<deduced_range_t>(range))} {}

    /// Construct from a range and a single position.
    ///
    /// @param range container to iterate over
    /// @param pos start and end position for iteration
    template <typename deduced_range_t>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range,
                                          range_size_t pos)
        : m_start{detray::ranges::begin(std::forward<deduced_range_t>(range)) +
                  pos},
          m_end{detray::ranges::next(m_start)} {}

    /// Construct from a range and start/end positions
    ///
    /// @param range container to iterate over
    /// @param pos start and end position for iteration
    template <typename deduced_range_t, typename index_range_t>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range,
                                          index_range_t &&pos)
        : m_start{detray::ranges::begin(std::forward<deduced_range_t>(range)) +
                  detray::detail::get<0>(pos)},
          m_end{detray::ranges::begin(std::forward<deduced_range_t>(range)) +
                detray::detail::get<1>(pos)} {}

    /// Construct from a range and start/end positions
    ///
    /// @param range container to iterate over
    /// @param pos start and end position for iteration
    template <typename deduced_range_t, typename volume_t,
              std::enable_if_t<std::is_same_v<volume_t::volume_def, volume_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, volume_t &&vol) {
        dindex_range r = vol.template range<
            typename detray::ranges::range_value_t<deduced_range_t>>();

        auto start =
            detray::ranges::begin(std::forward<deduced_range_t>(range));

        m_start = start + detray::detail::get<0>(r);
        m_end = start + detray::detail::get<1>(r);
    }

    /// Construct from an iterator pair
    ///
    /// @param start container to iterate over
    /// @param end start and end position for iteration
    DETRAY_HOST_DEVICE constexpr subrange(iterator_t &&start, iterator_t &&end)
        : m_start{std::forward<iterator_t>(start)},
          m_end{std::forward<iterator_t>(end)} {}

    /// @return start position of range.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator_t { return m_start; }

    /// @return start position of the range - const
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> const_iterator_t { return m_start; }

    /// @return sentinel of the range.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t { return m_end; }

    /// Start and end position of the subrange
    iterator_t m_start, m_end;
};

template <typename deduced_range_t>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range)->subrange<deduced_range_t>;

template <typename deduced_range_t>
DETRAY_HOST_DEVICE subrange(
    deduced_range_t &&range,
    typename detray::ranges::range_size_t<deduced_range_t> pos)
    ->subrange<deduced_range_t>;

template <typename deduced_range_t, typename index_range_t>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, index_range_t &&pos)
    ->subrange<deduced_range_t>;

template <typename deduced_range_t, typename volume_t,
          std::enable_if_t<std::is_same_v<volume_t::volume_def, volume_t>,
                           bool> = true>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, volume_t &&vol)
    ->subrange<deduced_range_t>;

template <typename deduced_range_t>
DETRAY_HOST_DEVICE subrange(
    typename detray::ranges::iterator_t<deduced_range_t> &&start,
    typename detray::ranges::iterator_t<deduced_range_t> &&end)
    ->subrange<deduced_range_t>;

}  // namespace detray::ranges
