/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <tuple>
#include <type_traits>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/// @brief Iterator-like access to a single value instead of a collection
template<typename value_t>
struct single_iterator {
    using container_type_iter = value_t;

    /// Delete default constructor
    single_iterator() = delete;

    /// Construct iterator from a value.
    ///
    /// @param value the single value that this iterator points to
    DETRAY_HOST_DEVICE single_iterator(const value_t & value)
        : m_value(&value) {}

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    inline auto begin() -> const value_t * { return m_value; }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    inline auto end() -> const value_t * { return m_value + 1; }

    /// @returns true if it points to the same value (not necessarily the same
    /// instance though).
    DETRAY_HOST_DEVICE
    bool operator!=(const single_iterator &rhs) {
        return m_value != rhs.m_value;
    }

    /// Does nothing
    DETRAY_HOST_DEVICE
    inline constexpr auto operator++() -> void {}

    /// @returns the single value that the iterator points to
    DETRAY_HOST_DEVICE
    inline constexpr auto operator*() const -> const value_t {
        return *m_value;
    }

    /// @return element at position i, relative to iterator range. */
    DETRAY_HOST_DEVICE
    inline constexpr auto operator[](const dindex /*i*/) -> value_t {
        return *m_value;
    }

    /// @return element at position i, relative to iterator range - const
    DETRAY_HOST_DEVICE
    inline constexpr auto operator[](const dindex /*i*/) const -> const value_t & {
        return *m_value;
    }

    /** @return the offset of the range start into the container. */
    DETRAY_HOST_DEVICE
    inline constexpr auto offset(const value_t &iterable) ->  { return 0; }

    value_t * value;
};

}  // namespace detray
