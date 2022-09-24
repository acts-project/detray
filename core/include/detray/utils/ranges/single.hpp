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

namespace detray::ranges {

/// @brief Struct that implements a view on a single element.
///
/// @tparam value_t type of the single element (outside of a container)
template <typename value_t>
class single_view
    : public detray::ranges::view_interface<single_view<value_t>> {

    private:
    /// @brief Emulates iterator behaviour for a single value.
    struct iterator {

        using difference_type = std::size_t;
        using value_type = value_t;
        using pointer = value_t *;
        using reference = value_t &;
        using iterator_category = std::bidirectional_iterator_tag;

        /// @returns true if it points to the same value.
        DETRAY_HOST_DEVICE
        constexpr auto operator==(const iterator &rhs) const -> bool {
            return m_value == rhs.m_value;
        }

        /// @returns false if it points to the same value.
        DETRAY_HOST_DEVICE
        constexpr auto operator!=(const iterator &rhs) const -> bool {
            return m_value != rhs.m_value;
        }

        /// Advance pointer to end immediately
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator & {
            m_value = m_sentinel;
            return *this;
        }

        /// Advance pointer to end immediately
        DETRAY_HOST_DEVICE
        constexpr auto operator--() -> iterator & {
            m_value = m_sentinel;
            return *this;
        }

        /// @returns the single value that the iterator points to - const
        DETRAY_HOST_DEVICE
        constexpr auto operator*() const -> const value_t & { return *m_value; }

        /// @returns the single value that the iterator points to
        DETRAY_HOST_DEVICE
        constexpr auto operator*() -> value_t & { return *m_value; }

        /// Advance the sequence: Immediately points to end
        DETRAY_HOST_DEVICE
        constexpr auto operator+(const difference_type j) const -> iterator {
            if (j == difference_type{0}) {
                return {m_value, m_sentinel};
            }
            return {m_sentinel, m_sentinel};
        }

        value_t *m_value, *m_sentinel;
    };

    /// Start and end values of the sequence
    value_t m_sentinel;
    iterator m_begin{nullptr}, m_end{nullptr};

    public:
    /// Default constructor
    single_view() = default;

    /// Construct iterator from the single @param value this iterator points to.
    DETRAY_HOST_DEVICE constexpr single_view(value_t &value)
        : m_sentinel{value},
          m_begin{&value, &m_sentinel},
          m_end{&m_sentinel, &m_sentinel} {}

    /// Construct iterator from the single @param value this iterator points to
    /// - const
    DETRAY_HOST_DEVICE constexpr single_view(const value_t &value)
        : m_sentinel{value},
          m_begin{&value, &m_sentinel},
          m_end{&m_sentinel, &m_sentinel} {}

    /// @return start position.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const -> iterator { return m_begin; }

    /// @return sentinel.
    DETRAY_HOST_DEVICE
    constexpr auto end() const -> iterator { return m_end; }
};

namespace views {

/// @brief interface type to construct a @c single_view with CTAD
template <typename value_t>
struct single : public detray::ranges::single_view<value_t> {

    using base_type = detray::ranges::single_view<value_t>;

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