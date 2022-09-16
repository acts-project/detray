/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/detail/single_iterator.hpp"
#include "detray/utils/ranges/ranges.hpp"

namespace detray::ranges {

/// @brief Struct that implements a view on a single element.
///
/// @tparam value_t typen of the single element (outside of a container)
template <typename value_t>
struct single_view
    : public detray::ranges::view_interface<single_view<value_t>> {

    using iterator_t = detail::single_iterator<value_t>;

    /// Default constructor
    constexpr single_view() = default;

    /// Construct iterator from a value.
    ///
    /// @param value the single value that this iterator points to
    DETRAY_HOST_DEVICE single_view(value_t &value)
        : m_begin{&value}, m_end{&value + 1} {}

    /// Construct iterator from a value.
    ///
    /// @param value the single value that this iterator points to
    DETRAY_HOST_DEVICE single_view(const value_t &value)
        : m_begin{&value}, m_end{&value + 1} {}

    /// @return start position.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator_t { return m_begin; }

    /// @return sentinel.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator_t { return m_end; }

    /// Start and end values of the sequence
    iterator_t m_begin{nullptr};
    iterator_t m_end{nullptr};
};

namespace views {

/// @brief interface type to construct a @c single_view with CTAD
template <typename value_t>
struct single : public detray::ranges::single_view<value_t> {

    using base_type = detray::ranges::single_view<value_t>;
    using iterator_t = typename base_type::iterator_t;

    /// Default constructor
    constexpr single() = default;

    template <typename deduced_value_t>
    DETRAY_HOST_DEVICE constexpr explicit single(deduced_value_t &value)
        : base_type(value) {}

    template <typename deduced_value_t>
    DETRAY_HOST_DEVICE constexpr explicit single(const deduced_value_t &value)
        : base_type(value) {}
};

// deduction guides

template <typename deduced_value_t>
single(deduced_value_t &value) -> single<deduced_value_t>;

template <typename deduced_value_t>
single(const deduced_value_t &value) -> single<const deduced_value_t>;

}  // namespace views

}  // namespace detray::ranges
