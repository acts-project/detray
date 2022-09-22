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

// System include(s)
#include <type_traits>

namespace detray::ranges {

/// @brief Implements a subrange by providing start and end iterators on
/// another range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange.
template <typename range_t>
class subrange : public ranges::view_interface<subrange<range_t>> {

    public:
    using iterator_t = typename detray::ranges::iterator_t<range_t>;
    using const_iterator_t = typename detray::ranges::const_iterator_t<range_t>;
    using range_size_t = typename detray::ranges::range_size_t<range_t>;

    /// Default constructor
    subrange() = default;

    /// Construct from an @param start and @param end iterator pair.
    DETRAY_HOST_DEVICE constexpr subrange(iterator_t &&start, iterator_t &&end)
        : m_start{std::forward<iterator_t>(start)},
          m_end{std::forward<iterator_t>(end)} {}

    /// Construct from a @param range.
    template <typename deduced_range_t>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range)
        : m_start{detray::ranges::begin(std::forward<deduced_range_t>(range))},
          m_end{detray::ranges::end(std::forward<deduced_range_t>(range))} {}

    /// Construct from a @param range and starting position @param pos. Used
    /// as an overload when only a single position is needed.
    template <typename deduced_range_t>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range,
                                          range_size_t pos)
        : m_start{detray::ranges::begin(std::forward<deduced_range_t>(range)) +
                  pos},
          m_end{detray::ranges::next(m_start)} {}

    /// Construct from a @param range and an index range provided by a volume
    /// @param vol.
    template <typename deduced_range_t, typename volume_t,
              typename = typename std::remove_reference_t<volume_t>::volume_def>
    DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, const volume_t &vol) {
        const dindex_range r = vol.template range<
            typename detray::ranges::range_value_t<deduced_range_t>>();

        auto start =
            detray::ranges::begin(std::forward<deduced_range_t>(range));

        m_start = start + detray::detail::get<0>(r);
        m_end = start + detray::detail::get<1>(r);
    }

    /// Construct from a @param range and an index range @param pos.
    template <typename deduced_range_t, typename index_range_t,
              std::enable_if_t<detray::detail::is_interval_v<index_range_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range,
                                          index_range_t &&pos)
        : m_start{detray::ranges::begin(std::forward<deduced_range_t>(range)) +
                  detray::detail::get<0>(pos)},
          m_end{detray::ranges::begin(std::forward<deduced_range_t>(range)) +
                detray::detail::get<1>(pos)} {}

    /// @return start position of range.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator_t { return m_start; }

    /// @return start position of the range - const
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> const_iterator_t { return m_start; }

    /// @return sentinel of the range.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t { return m_end; }

    private:
    /// Start and end position of the subrange
    iterator_t m_start, m_end;
};

/*template <typename deduced_range_t>
DETRAY_HOST_DEVICE subrange(
    typename detray::ranges::iterator_t<deduced_range_t> &&start,
    typename detray::ranges::iterator_t<deduced_range_t> &&end)
    ->subrange<deduced_range_t>;*/

template <typename deduced_range_t>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range)->subrange<deduced_range_t>;

template <typename deduced_range_t>
DETRAY_HOST_DEVICE subrange(
    deduced_range_t &&range,
    typename detray::ranges::range_size_t<deduced_range_t> pos)
    ->subrange<deduced_range_t>;

template <
    typename deduced_range_t, typename index_range_t,
    std::enable_if_t<detray::detail::is_interval_v<index_range_t>, bool> = true>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, index_range_t &&pos)
    ->subrange<deduced_range_t>;

template <typename deduced_range_t, typename volume_t,
          typename = typename std::remove_reference_t<volume_t>::volume_def>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, const volume_t &vol)
    ->subrange<deduced_range_t>;

}  // namespace detray::ranges
