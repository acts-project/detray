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
#include <type_traits>

namespace detray::ranges {

/// @brief Implements a subrange by providing start and end iterators on
/// another range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange.
template <typename range_t>
class subrange : public detray::ranges::view_interface<subrange<range_t>> {

    public:
    using iterator_t = typename detray::ranges::iterator_t<range_t>;
    using const_iterator_t = typename detray::ranges::const_iterator_t<range_t>;
    using difference_t = typename detray::ranges::range_difference_t<range_t>;

    /// Default constructor
    subrange() = default;

    /// Construct from an @param start and @param end iterator pair.
    template <typename deduced_itr_t,
              std::enable_if_t<std::is_same_v<deduced_itr_t, iterator_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_itr_t &&start,
                                          deduced_itr_t &&end)
        : m_begin{std::forward<deduced_itr_t>(start)},
          m_end{std::forward<deduced_itr_t>(end)} {}

    /// Construct from a @param range.
    template <
        typename deduced_range_t,
        std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range)
        : m_begin{detray::ranges::begin(range)},
          m_end{detray::ranges::end(range)} {}

    /// Construct from a @param range and starting position @param pos. Used
    /// as an overload when only a single position is needed.
    template <
        typename deduced_range_t,
        std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range,
                                          difference_t pos)
        : m_begin{detray::ranges::next(detray::ranges::begin(range), pos)},
          m_end{detray::ranges::next(m_begin)} {}

    /// Construct from a @param range and an index range provided by a volume
    /// @param vol.
    template <
        typename deduced_range_t, typename volume_t,
        typename value_t = detray::ranges::range_value_t<deduced_range_t>,
        std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true,
        typename = typename std::remove_reference_t<volume_t>::volume_def>
    DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, const volume_t &vol)
        : subrange(std::forward<deduced_range_t>(range), vol.get_all()) {}

    /// Construct from a @param range and an index range @param pos.
    template <
        typename deduced_range_t, typename index_range_t,
        std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true,
        std::enable_if_t<detray::detail::is_interval_v<index_range_t>, bool> =
            true>
    DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range,
                                          index_range_t &&pos)
        : m_begin{detray::ranges::next(detray::ranges::begin(range),
                                       detray::detail::get<0>(pos))},
          m_end{detray::ranges::next(detray::ranges::begin(range),
                                     detray::detail::get<1>(pos))} {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    subrange &operator=(const subrange &other) {
        m_begin = other.m_begin;
        m_end = other.m_end;
        return *this;
    };

    /// @return start position of range.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator_t { return m_begin; }

    /// @return start position of the range - const
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator_t { return m_begin; }

    /// @return sentinel of the range.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t { return m_end; }

    private:
    /// Start and end position of the subrange
    iterator_t m_begin, m_end;
};

/*template <typename deduced_range_t,
          std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> =
true> DETRAY_HOST_DEVICE subrange( typename
detray::ranges::iterator_t<deduced_range_t> &&start, typename
detray::ranges::iterator_t<deduced_range_t> &&end)
    ->subrange<deduced_range_t>;*/

template <
    typename deduced_range_t,
    std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range)->subrange<deduced_range_t>;

template <
    typename deduced_range_t,
    std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true>
DETRAY_HOST_DEVICE subrange(
    deduced_range_t &&range,
    typename detray::ranges::range_difference_t<deduced_range_t> pos)
    ->subrange<deduced_range_t>;

template <
    typename deduced_range_t, typename index_range_t,
    std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true,
    std::enable_if_t<detray::detail::is_interval_v<index_range_t>, bool> = true>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, index_range_t &&pos)
    ->subrange<deduced_range_t>;

template <
    typename deduced_range_t, typename volume_t,
    std::enable_if_t<detray::ranges::range_v<deduced_range_t>, bool> = true,
    typename = typename std::remove_reference_t<volume_t>::volume_def>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, const volume_t &vol)
    ->subrange<deduced_range_t>;

/// @see https://en.cppreference.com/w/cpp/ranges/borrowed_iterator_t
template <class R>
using borrowed_subrange_t =
    std::conditional_t<borrowed_range<R>, detray::ranges::subrange<R>,
                       dangling>;

}  // namespace detray::ranges
