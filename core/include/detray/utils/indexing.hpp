/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <limits>

namespace detray {
using dindex = unsigned long;
dindex constexpr dindex_invalid = std::numeric_limits<dindex>::max();
using dindex_range = darray<dindex, 2>;
using dindex_sequence = dvector<dindex>;

template <typename id_type = unsigned int, typename index_type = dindex>
struct typed_index {
    id_type _object_id;
    index_type _index;

    /** Equality operator */
    DETRAY_HOST_DEVICE
    bool operator==(const typed_index<id_type, index_type>& rhs) const {
        return (_object_id == rhs._object_id && _index == rhs._index);
    }

    /** Equality operator */
    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type> operator+(
        const typed_index<id_type, index_type>& rhs) const {
        return {_object_id, _index + rhs._index};
    }

    /** Equality operator */
    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type> operator+(const index_type& index) const {
        return {_object_id, _index + index};
    }

    /** Equality operator */
    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type> operator-(
        const typed_index<id_type, index_type>& rhs) const {
        return {_object_id, _index - rhs._index};
    }

    /** Equality operator */
    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type> operator-(const index_type& index) const {
        return {_object_id, _index - index};
    }

    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type>& operator+=(
        const typed_index<id_type, index_type>& rhs) {
        _index += rhs._index;
        return *this;
    }

    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type>& operator+=(const index_type& index) {
        _index += index;
        return *this;
    }

    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type>& operator-=(
        const typed_index<id_type, index_type>& rhs) {
        _index -= rhs._index;
        return *this;
    }

    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type>& operator-=(const index_type& index) {
        _index -= index;
        return *this;
    }

    DETRAY_HOST_DEVICE
    typed_index<id_type, index_type>& operator++() {
        _index++;
        return *this;
    }
};

namespace detail {

using std::get;

template <std::size_t ID, typename id_type, typename index_type,
          std::enable_if_t<ID == 0, bool> = true>
DETRAY_HOST_DEVICE constexpr auto& get(
    const typed_index<id_type, index_type>& index) noexcept {
    return index._object_id;
}

template <std::size_t ID, typename id_type, typename index_type,
          std::enable_if_t<ID == 1, bool> = true>
DETRAY_HOST_DEVICE constexpr auto& get(
    const typed_index<id_type, index_type>& index) noexcept {
    return index._index;
}

}  // namespace detail
}  // namespace detray
