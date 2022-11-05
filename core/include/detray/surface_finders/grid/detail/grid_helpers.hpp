/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

namespace detail {

// TODO: Replace by iterator based approach, along the lines of std::ranges

/// @brief data state of a grid - owning
///
/// The data that is owned and managed by the grid. Does not contain the
/// data of the axes, as that is managed by the multi-axis type directly.
template <typename backend_storage_type>
struct grid_data {
    /// Contains all bin data
    backend_storage_type m_bin_data{};

    /// Default constructor
    grid_data() = default;

    /// Construct containers using a memory resources
    DETRAY_HOST_DEVICE
    grid_data(vecmem::memory_resource &resource) : m_bin_data(&resource) {}

    /// Construct grid data from containers - move
    DETRAY_HOST_DEVICE
    grid_data(backend_storage_type &&bin_data)
        : m_bin_data(std::move(bin_data)) {}

    /// Construct containers from a vecmem view
    DETRAY_HOST_DEVICE
    grid_data(
        const dvector_view<typename backend_storage_type::value_type> &view)
        : m_bin_data(view) {}

    /// @returns pointer to the entire bin data for the grid - const
    DETRAY_HOST_DEVICE
    auto bin_data() const -> const backend_storage_type * {
        return &m_bin_data;
    }

    /// @returns pointer to the entire bin data for the grid - non-const
    DETRAY_HOST_DEVICE
    auto bin_data() -> backend_storage_type * { return &m_bin_data; }

    /// @returns the offset of this grid in a global grid collection
    /// (not needed if the grid owns its data)
    DETRAY_HOST_DEVICE
    constexpr auto offset() const -> dindex { return 0; }
};

/// @brief data state of a grid - non-owning (i.e. container view)
///
/// The grid holds a view onto a global collection of data. This type is NOT
/// used to move @c grid_data from host to device.
/// @todo Use spans
template <typename backend_storage_type>
struct grid_view {
    /// Contains all bin edges for all the axes
    backend_storage_type *m_bin_data{nullptr};
    /// Offset for this grid into the global container
    dindex m_offset{0};

    /// Default constructor
    grid_view() = default;

    /// Construct from explicitly given data pointers
    DETRAY_HOST_DEVICE
    grid_view(backend_storage_type *bin_data_ptr, const dindex offset)
        : m_bin_data(bin_data_ptr), m_offset{offset} {}

    /// @returns pointer to the entire bin edges collection for every axis -
    /// const
    DETRAY_HOST_DEVICE
    auto bin_data() const -> const backend_storage_type * { return m_bin_data; }

    /// @returns pointer to the entire bin edges collection for every axis -
    /// non_const
    DETRAY_HOST_DEVICE
    auto bin_data() -> backend_storage_type * { return m_bin_data; }

    /// @returns offset of this multi-axis in the global data container
    DETRAY_HOST_DEVICE
    auto offset() const -> dindex { return m_offset; }
};

}  // namespace detail

}  // namespace detray
