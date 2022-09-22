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

/// @brief Range factory that produces a sequence of values.
///
/// @see https://en.cppreference.com/w/cpp/ranges/iota
///
/// @tparam incr_t the incrementable type that makes up the sequence
template <typename incr_t>
class iota_view : public detray::ranges::view_interface<iota_view<incr_t>> {

    private:
    /// @brief Nested iterator to generate a range of values on demand.
    struct iterator {

        using difference_type = incr_t;
        using value_type = incr_t;
        using pointer = incr_t *;
        using reference = incr_t &;
        using iterator_category = std::input_iterator_tag;

        /// @returns true if we reach end of sequence
        DETRAY_HOST_DEVICE
        constexpr auto operator==(const iterator &rhs) const -> bool {
            return (i == rhs.i);
        }

        /// @returns true if we reach end of sequence
        DETRAY_HOST_DEVICE
        constexpr auto operator!=(const iterator &rhs) const -> bool {
            return (i != rhs.i);
        }

        /// Increment the index
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator & {
            ++i;
            return *this;
        }

        /// @returns the current value in the sequence
        DETRAY_HOST_DEVICE
        constexpr auto operator*() const -> const incr_t & { return i; }

        /// @returns the current value in the sequence
        DETRAY_HOST_DEVICE
        constexpr auto operator*() -> incr_t & { return i; }

        /// Advance the sequence by @param j positions
        DETRAY_HOST_DEVICE
        constexpr auto operator+(const incr_t j) const -> iterator {
            return {i + j};
        }

        /// Current value of sequence
        incr_t i;
    };

    /// Start and end values of the sequence
    incr_t m_start, m_end;

    public:
    /// Default constructor (only works if @c imrementable_t is default
    /// constructible)
    iota_view() = default;

    template <typename interval_t,
              std::enable_if_t<detray::detail::is_interval_v<interval_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota_view(interval_t &&interval)
        : m_start{detray::detail::get<0>(interval)},
          m_end{detray::detail::get<1>(interval)} {}

    /// Construct from a @param start start and @param end value.
    template <typename deduced_incr_t>
    DETRAY_HOST_DEVICE constexpr iota_view(deduced_incr_t &&start,
                                           deduced_incr_t &&end)
        : m_start{start}, m_end{end - 1} {}

    /// Construct from just a @param start value to represent a single value seq
    template <
        typename deduced_incr_t,
        std::enable_if_t<not detray::detail::is_interval_v<deduced_incr_t>,
                         bool> = true>
    DETRAY_HOST_DEVICE constexpr explicit iota_view(deduced_incr_t &&start)
        : m_start{start}, m_end{start + 1} {}

    /// @return start position of range on container.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator { return {m_start}; }

    /// @return sentinel of a sequence.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator { return {m_end}; }

    /// @note Cannot peak at the end of input-iterator based range
    constexpr typename iterator::value_type back() noexcept = delete;
};

namespace views {

/// @brief interface type to construct a @c iota_view with CTAD
template <typename incr_t, typename sentinel_t = void>
struct iota : public detray::ranges::iota_view<incr_t> {

    using base_type = detray::ranges::iota_view<incr_t>;

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

}  // namespace detray::ranges