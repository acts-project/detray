/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <limits>

namespace detray {

using dindex = unsigned int;
inline constexpr dindex dindex_invalid = detail::invalid_value<dindex>();
using dindex_range = darray<dindex, 2>;
using dindex_sequence = dvector<dindex>;

/// @brief Simple multi-index structure
///
/// @tparam DIM number of indices that are held by this type
/// @tparam index_t type of indices
template <typename index_t = dindex, std::size_t DIM = 3u>
struct dmulti_index {
    using index_type = index_t;

    std::array<index_t, DIM> indices{};

    /// @returns the number of conatined indices
    DETRAY_HOST_DEVICE
    constexpr static auto size() -> std::size_t { return DIM; }

    /// Elementwise access.
    DETRAY_HOST_DEVICE
    auto operator[](const std::size_t i) -> index_t& { return indices[i]; }
    DETRAY_HOST_DEVICE
    auto operator[](const std::size_t i) const -> const index_t& {
        return indices[i];
    }

    /// Equality operator @returns true if all bin indices match.
    DETRAY_HOST_DEVICE
    auto operator==(const dmulti_index<index_t, DIM>& rhs) const -> bool {
        return (indices == rhs.indices);
    }
};

/// @brief Ties an object type and an index into a container together.
///
/// @tparam id_type Represents the type of object that is being indexed
/// @tparam index_type The type of indexing needed for the indexed type's
///         container (e.g. single index, range, multiindex)
template <typename id_t = unsigned int, typename index_t = dindex>
struct dtyped_index {

    using id_type = id_t;
    using index_type = index_t;

    id_type _object_id;
    index_type _index;

    /// @return the type id
    auto id() const -> id_type { return _object_id; }

    /// @return a reference to the index - const
    DETRAY_HOST_DEVICE
    auto index() const -> const index_type& { return _index; }

    /// @return a reference to the index - non-const
    DETRAY_HOST_DEVICE
    auto index() -> index_type& { return _index; }

    /// Equality operator
    DETRAY_HOST_DEVICE
    bool operator==(const dtyped_index<id_type, index_type>& rhs) const {
        return (_object_id == rhs._object_id && _index == rhs._index);
    }

    /// Arithmetic operators
    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type> operator+(
        const dtyped_index<id_type, index_type>& rhs) const {
        return {_object_id, _index + rhs._index};
    }

    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type> operator+(const index_type& index) const {
        return {_object_id, _index + index};
    }

    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type> operator-(
        const dtyped_index<id_type, index_type>& rhs) const {
        return {_object_id, _index - rhs._index};
    }

    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type> operator-(const index_type& index) const {
        return {_object_id, _index - index};
    }

    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type>& operator+=(
        const dtyped_index<id_type, index_type>& rhs) {
        _index += rhs._index;
        return *this;
    }

    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type>& operator+=(const index_type& index) {
        _index += index;
        return *this;
    }

    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type>& operator-=(
        const dtyped_index<id_type, index_type>& rhs) {
        _index -= rhs._index;
        return *this;
    }

    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type>& operator-=(const index_type& index) {
        _index -= index;
        return *this;
    }

    /// Only make the prefix operator available
    DETRAY_HOST_DEVICE
    dtyped_index<id_type, index_type>& operator++() {
        ++_index;
        return *this;
    }
};

namespace detail {

/// Stub function to get a single index
template <std::size_t ID>
dindex get(dindex idx) noexcept {
    return idx;
}

/// Custom get function for the dtyped_index struct. Get the type.
template <std::size_t idx, typename index_type, std::size_t index_size>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const dmulti_index<index_type, index_size>& index) noexcept {
    return index.indices[idx];
}

/// Custom get function for the dtyped_index struct. Get the type.
template <std::size_t idx, typename index_type, std::size_t index_size>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    dmulti_index<index_type, index_size>& index) noexcept {
    return index.indices[idx];
}

/// Custom get function for the dtyped_index struct. Get the type.
template <std::size_t ID, typename id_type, typename index_type,
          std::enable_if_t<ID == 0, bool> = true>
DETRAY_HOST_DEVICE constexpr auto& get(
    const dtyped_index<id_type, index_type>& index) noexcept {
    return index._object_id;
}

/// Custom get function for the dtyped_index struct. Get the index.
template <std::size_t ID, typename id_type, typename index_type,
          std::enable_if_t<ID == 1, bool> = true>
DETRAY_HOST_DEVICE constexpr auto& get(
    const dtyped_index<id_type, index_type>& index) noexcept {
    return index._index;
}

}  // namespace detail

}  // namespace detray
