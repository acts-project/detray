/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <iterator>

namespace detray::ranges {

namespace detail {

template <typename I>
struct single_iterator;

}

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

namespace detail {

/// @brief Emulates iterator behaviour for a single value.
template <typename value_t>
struct single_iterator {

    /// @returns true if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const single_iterator<value_t> &rhs) const
        -> bool {
        return *m_value == *rhs.m_value;
    }

    /// @returns false if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    constexpr auto operator!=(const single_iterator &rhs) const -> bool {
        return *m_value != *rhs.m_value;
    }

    /// Does nothing
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> single_iterator<value_t> {
        ++m_value;
        return *this;
    }

    /// Does nothing
    DETRAY_HOST_DEVICE
    constexpr auto operator--() -> single_iterator<value_t> {
        --m_value;
        return *this;
    }

    /// @returns the single value that the iterator points to - const
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const -> const value_t & { return *m_value; }

    /// @returns the single value that the iterator points to
    DETRAY_HOST_DEVICE
    constexpr auto operator*() -> value_t & { return *m_value; }

    /// Advance the sequence
    DETRAY_HOST_DEVICE
    constexpr auto operator+(const value_t j) const
        -> single_iterator<value_t> {
        return {m_value + j};
    }

    value_t *m_value;
};

}  // namespace detail

}  // namespace detray::ranges

namespace std {

/// Specialization of std::iterator_traits struct for the a single value
template <typename T>
struct iterator_traits<detray::ranges::detail::single_iterator<T>> {
    using difference_type = std::size_t;
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using iterator_category = std::bidirectional_iterator_tag;
};

}  // namespace std
