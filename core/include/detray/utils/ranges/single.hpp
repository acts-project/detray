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
/// @see https://en.cppreference.com/w/cpp/ranges/single_view
///
/// @tparam value_t type of the single element (outside of a container)
///
/// @note Does not take ownership of the value it operates on. Its lifetime
/// needs to be guaranteed throughout iteration or between iterations with the
/// same view instance.
/// @note Is not fit for lazy evaluation.
template <typename value_t>
class single_view
    : public detray::ranges::view_interface<single_view<value_t>> {

    private:
    /// @brief Emulates iterator behaviour for a single value.
    struct iterator
        : public std::iterator<std::bidirectional_iterator_tag, value_t> {

        iterator() = default;

        DETRAY_HOST_DEVICE
        iterator(value_t *const value, value_t *const sentinel)
            : m_value{value}, m_sentinel{sentinel} {}

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

        value_t *m_value, *m_sentinel;
    };

    /// Start and end values of the sequence
    value_t *m_sentinel{nullptr};
    iterator m_begin{nullptr, nullptr}, m_end{nullptr, nullptr};

    public:
    using iterator_t = iterator;

    /// Default constructor
    single_view() = default;

    /// Construct iterator from the single @param value this iterator points to.
    DETRAY_HOST_DEVICE constexpr explicit single_view(value_t &value)
        : m_sentinel{&value},
          m_begin{&value, m_sentinel},
          m_end{m_sentinel, m_sentinel} {}

    /// Construct iterator from the single @param value this iterator points to
    /// - const
    DETRAY_HOST_DEVICE constexpr explicit single_view(const value_t &value)
        : m_sentinel{&value},
          m_begin{&value, m_sentinel},
          m_end{m_sentinel, m_sentinel} {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    single_view &operator=(const single_view &other) {
        m_sentinel = other.m_sentinel;
        m_begin = other.m_begin;
        m_end = other.m_end;

        return *this;
    };

    /// @return start position.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const noexcept -> iterator { return m_begin; }

    /// @return sentinel.
    DETRAY_HOST_DEVICE
    constexpr auto end() const noexcept -> iterator { return m_end; }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - const
    DETRAY_HOST_DEVICE
    constexpr auto data() const noexcept -> const
        typename iterator_t::value_type * {
        return m_begin.m_value;
    }

    /// @returns a pointer to the beginning of the data of the first underlying
    /// range - non-const
    DETRAY_HOST_DEVICE
    constexpr auto data() noexcept -> typename iterator_t::value_type * {
        return m_begin.m_value;
    }

    /// @returns the size of the single view, which is always 'one'.
    /// @note overwrites the inherited function from the view interface.
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> std::size_t { return 1; }
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