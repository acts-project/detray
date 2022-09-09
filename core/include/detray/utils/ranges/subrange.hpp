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

namespace detray {

namespace views {

/// @brief Struct that implements a subrange by providing start and end
/// iterators on a range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange
template <typename range_t>
struct subrange : public view_interface<subrange<range_t>> {

    /// non-owning range
    using range_type = ranges::range<range_t>;

    /// Compliance with std::ranges typedefs
    using iterator_t = typename range_type::iterator_t;
    using const_iterator_t = typename range_type::const_iterator_t;
    using range_size_t = typename range_type::range_size_t;
    using range_difference_t = typename range_type::range_difference_t;
    using range_reference_t = typename range_type::range_reference_t;
    using range_const_reference_t =
        typename range_type::range_const_reference_t;
    using range_rvalue_reference_t =
        typename range_type::range_rvalue_reference_t;

    /// Delete default constructor
    subrange() = delete;

    /// Construct from a range: The subrange spans the entire range
    ///
    /// @param range container to iterate over
    template <typename deduced_range_t>
    DETRAY_HOST_DEVICE subrange(deduced_range_t &range)
        : m_start{ranges::begin(range)}, m_end{ranges::end(range)} {}

    /// Construct from a range and start/end positions
    ///
    /// @param range container to iterate over
    /// @param pos start and end position for iteration
    template <typename deduced_range_t, typename index_range_t>
    DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, index_range_t &&pos)
        : m_start{ranges::begin(range) + detail::get<0>(pos)},
          m_end{ranges::begin(range) + detail::get<1>(pos)} {}

    /// Construct from an iterator pair
    ///
    /// @param start container to iterate over
    /// @param end start and end position for iteration
    /*template<typename deduced_range_t>
    DETRAY_HOST_DEVICE
    subrange(ranges::detail::iterator_traits<deduced_range_t>::iterator_t
    &&start, ranges::detail::iterator_traits<deduced_range_t>::iterator_t &&end)
        : m_start{start},
          m_end{end} {}*/

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr inline auto begin() -> iterator_t & { return m_start; }

    /// @return end position of range on container.
    DETRAY_HOST_DEVICE
    constexpr inline auto end() -> iterator_t & { return m_end; }

    /// Does this describe the same range?
    DETRAY_HOST_DEVICE
    auto operator!=(const subrange<range_t> &rhs) -> bool {
        return m_start != rhs.m_start or m_end != rhs.m_end;
    }

    /// Does this describe the same range?
    DETRAY_HOST_DEVICE
    auto operator++() -> void { ++m_start; }

    /// Start and end position of the subrange
    iterator_t m_start, m_end;
};

template <typename deduced_range_t>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range)->subrange<deduced_range_t>;

template <typename deduced_range_t, typename index_range_t>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, index_range_t &&pos)
    ->subrange<deduced_range_t>;

/*template<typename deduced_range_t>
DETRAY_HOST_DEVICE
subrange(ranges::detail::iterator_traits<deduced_range_t>::iterator_t &&start,
         ranges::detail::iterator_traits<deduced_range_t>::iterator_t &&end)
    -> subrange<deduced_range_t>;*/

/*template<typename deduced_range_t>
DETRAY_HOST_DEVICE
subrange(deduced_range_t &&range) -> subrange<deduced_range_t>;

} // namespace detray::views

/// Get the subrange on an iterable - convenience function
///
/// @tparam iterable_t the type of iterable. Must fulfill
///                    @c detray::ranges::iterable
///
/// @param iterable reference to an iterable range
/// @param range the subrange
///
/// @returns a subrange on the iterable
template <typename iterable_t>
DETRAY_HOST_DEVICE inline constexpr auto range(const iterable_t &iterable,
                                               const dindex_range &range) {

    return detray::views::subrange(iterable, range);
}

/// Overload of the @c range function, that extracts the range from a
/// @param volume.
///
/// @returns an subrange on the @param iterable, according to the volume's range
template <typename iterable_t,
          typename = typename std::remove_reference_t<iterable_t>::value_type,
          typename volume_t,
          typename = typename std::remove_reference_t<volume_t>::volume_def>
DETRAY_HOST_DEVICE inline constexpr auto range(const iterable_t &iterable,
                                               volume_t &&volume) {
    return detray::views::subrange(
        iterable, volume.template range<typename iterable_t::value_type>());
}

/// Overload of the @c range function for a single index
template <typename iterable_t>
DETRAY_HOST_DEVICE inline constexpr auto range(const iterable_t &iterable,
                                               const dindex &i) {

    return detray::views::subrange(iterable, dindex_range{i, i + 1});
}*/

}  // namespace views

}  // namespace detray
