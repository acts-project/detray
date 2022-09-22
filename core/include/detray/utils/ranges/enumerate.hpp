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
#include <tuple>
#include <type_traits>

namespace detray::ranges {

/// @brief Enumerates the elements of a range on the fly.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_itr_t the iterator type of the enumerated range
/// @tparam incr_t a type that can be incremented in lockstep with the
///         iterator 'range_itr_t'.
template <typename range_itr_t, typename incr_t = dindex>
class enumerate_view : public detray::ranges::view_interface<
                           enumerate_view<range_itr_t, incr_t>> {

    private:
    /// @brief Nested iterator to enumerate the elements of a range.
    ///
    /// The enumeration is done by incrementing an index in lockstep with a
    /// wrapped iterator of a range. Index and current iterator value are
    /// returned using structured binding.
    struct iterator {

        using difference_type =
            typename std::iterator_traits<range_itr_t>::difference_type;
        using value_type =
            typename std::iterator_traits<range_itr_t>::value_type;
        using pointer = typename std::iterator_traits<range_itr_t>::pointer;
        using reference = typename std::iterator_traits<range_itr_t>::reference;
        using iterator_category = std::input_iterator_tag;

        /// Determine whether we reach end of range
        DETRAY_HOST_DEVICE
        constexpr auto operator!=(const iterator &rhs) const -> bool {
            return (m_iter != rhs.m_iter);
        }

        /// Increment iterator and index in lockstep
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator & {
            ++m_i;
            ++m_iter;
            return *this;
        }

        /// @return iterator and index together
        DETRAY_HOST_DEVICE
        constexpr auto operator*() const { return std::tie(m_i, *m_iter); }

        /// Advance the iterator and index position by @param j.
        DETRAY_HOST_DEVICE
        constexpr auto operator+(const incr_t j) const -> iterator {
            return {m_iter + j, m_i + j};
        }

        range_itr_t m_iter;
        incr_t m_i;
    };

    iterator m_range, m_end;

    public:
    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    enumerate_view() = default;

    /// Construct from a @param range that will be enumerated beginning at 0
    template <
        typename range_t,
        std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate_view(range_t &&range)
        : m_range{detray::ranges::begin(std::forward<range_t>(range)), 0},
          m_end{detray::ranges::end(std::forward<range_t>(range)), 0} {}

    /// Construct from a @param range that will be enumerated beginning at
    /// @param start.
    template <
        typename range_t,
        std::enable_if_t<detray::ranges::range<range_t>::value, bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit enumerate_view(range_t &&rng,
                                                         dindex start)
        : m_range{detray::ranges::begin(std::forward<range_t>(rng)), start},
          m_end{detray::ranges::end(std::forward<range_t>(rng)), 0} {}

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() -> iterator & { return m_range; }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> const iterator & { return m_range; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> const iterator & { return m_end; }

    /// @note Cannot peek at the end of input-iterator based range
    constexpr typename iterator::value_type back() noexcept = delete;
};

namespace views {

template <typename range_itr_t, typename incr_t = dindex>
struct enumerate : public enumerate_view<range_itr_t, incr_t> {

    using base_type = enumerate_view<range_itr_t, incr_t>;

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

    /// Construct from a @param range and an index range provided by a volume
    /// @param vol.
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

}  // namespace detray::ranges