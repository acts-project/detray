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

/// @brief Range factory that produces a sequence of values.
///
/// @see https://en.cppreference.com/w/cpp/ranges/iota_view
///
/// @tparam incr_t the incrementable type that makes up the sequence
///
/// @note If given single value, does not do infinite iteration, but only jumps
///       to next value.
/// @note Is not fit for lazy evaluation.
template <typename incr_t>
class iota_view : public detray::ranges::view_interface<iota_view<incr_t>> {

    private:
    /// @brief Nested iterator to generate a range of values on demand.
    struct iterator : public std::iterator<detray::ranges::input_iterator_tag,
                                           incr_t, incr_t> {

        /// Default construction only works if incr_t is default constructible
        iterator() = default;

        /// Parametrized Constructor
        DETRAY_HOST_DEVICE
        constexpr explicit iterator(const incr_t i) : m_i{i} {}

        /// @returns true if incremetables are the same
        DETRAY_HOST_DEVICE
        constexpr auto operator==(const iterator &rhs) const -> bool {
            return (m_i == rhs.m_i);
        }

        /// @returns true while incrementables are different
        DETRAY_HOST_DEVICE
        constexpr auto operator!=(const iterator &rhs) const -> bool {
            return (m_i != rhs.m_i);
        }

        /// Increment the index
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator & {
            ++m_i;
            return *this;
        }

        /// @returns the current value in the sequence - copy
        DETRAY_HOST_DEVICE
        constexpr auto operator*() const noexcept -> const incr_t & {
            return m_i;
        }

        /// @returns the current value in the sequence - copy
        DETRAY_HOST_DEVICE
        constexpr auto operator*() noexcept -> incr_t & { return m_i; }

        /// Current value of sequence
        incr_t m_i;
    };

    /// Start and end values of the sequence
    incr_t m_start, m_end;

    public:
    using iterator_t = iterator;

    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    iota_view() = default;

    /// Construct from just a @param start value to represent a single value seq
    template <
        typename deduced_incr_t,
        std::enable_if_t<not detray::detail::is_interval_v<deduced_incr_t>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota_view(deduced_incr_t &&start)
        : m_start{std::forward<deduced_incr_t>(start)},
          m_end{std::forward<deduced_incr_t>(start) + 1} {}

    /// Construct from an @param interval that defines start and end values.
    template <typename interval_t,
              std::enable_if_t<detray::detail::is_interval_v<interval_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota_view(interval_t &&interval)
        : m_start{detray::detail::get<0>(std::forward<interval_t>(interval))},
          m_end{detray::detail::get<1>(std::forward<interval_t>(interval))} {}

    /// Construct from a @param start start and @param end value.
    template <typename deduced_incr_t>
    DETRAY_HOST_DEVICE constexpr iota_view(deduced_incr_t &&start,
                                           deduced_incr_t &&end)
        : m_start{std::forward<deduced_incr_t>(start)},
          m_end{std::forward<deduced_incr_t>(end) - 1} {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    iota_view &operator=(const iota_view &other) {
        m_start = other.m_start;
        m_end = other.m_end;
        return *this;
    }

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator_t { return iterator_t{m_start}; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t { return iterator_t{m_end}; }

    /// @returns the span of the sequence
    DETRAY_HOST_DEVICE
    constexpr auto size() const -> incr_t { return m_end - m_start; }
};

namespace views {

/// @brief interface type to construct a @c iota_view with CTAD
template <typename incr_t, typename sentinel_t = void>
struct iota : public detray::ranges::iota_view<incr_t> {

    using base_type = detray::ranges::iota_view<incr_t>;

    constexpr iota() = default;

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
    ->iota<detray::detail::remove_cvref_t<
        decltype(std::get<0>(std::declval<deduced_interval_t>()))>>;

template <typename deduced_incr_t = dindex>
DETRAY_HOST_DEVICE iota(deduced_incr_t &&start, deduced_incr_t &&end)
    ->iota<detray::detail::remove_cvref_t<deduced_incr_t>>;

template <typename deduced_incr_t,
          std::enable_if_t<not detray::detail::is_interval_v<deduced_incr_t>,
                           bool> = true>
DETRAY_HOST_DEVICE iota(deduced_incr_t &&start)
    ->iota<detray::detail::remove_cvref_t<deduced_incr_t>>;

}  // namespace views

}  // namespace detray::ranges