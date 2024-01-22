/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::detail {

/// @brief bin data state of a grid
///
/// Can be data-owning or not. Does not contain the data of the axes,
/// as that is managed by the multi-axis type directly.
template <bool is_owning, typename bin_t, typename containers>
class bin_storage : public detray::ranges::view_interface<
                        bin_storage<is_owning, bin_t, containers>> {

    template <typename T>
    using vector_t = typename containers::template vector_type<T>;
    using bin_range_t =
        std::conditional_t<is_owning, vector_t<bin_t>,
                           detray::ranges::subrange<vector_t<bin_t>>>;

    public:
    /// Bin type: single or static_array
    using bin_type = bin_t;
    /// Backend storage type for the grid
    using bin_container_type = vector_t<bin_t>;

    // Vecmem based view type
    using view_type = dvector_view<bin_type>;
    using const_view_type = dvector_view<const bin_type>;

    // Vecmem based buffer type
    using buffer_type = dvector_buffer<bin_type>;

    /// Default constructor
    bin_storage() = default;

    /// Construct containers using a memory resources
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST bin_storage(vecmem::memory_resource& resource)
        : m_bin_data(&resource) {}

    /// Construct grid data from containers - move
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST_DEVICE bin_storage(bin_container_type&& bin_data)
        : m_bin_data(std::move(bin_data)) {}

    /// Construct the non-owning type from the @param offset into the global
    /// container @param bin_data and the number of bins @param size
    template <bool owner = is_owning, std::enable_if_t<!owner, bool> = true>
    DETRAY_HOST_DEVICE bin_storage(bin_container_type& bin_data, dindex offset,
                                   dindex size)
        : m_bin_data(bin_data, dindex_range{offset, offset + size}) {}

    /// Construct bin storage from its vecmem view
    template <typename view_t,
              typename std::enable_if_t<detail::is_device_view_v<view_t>,
                                        bool> = true>
    DETRAY_HOST_DEVICE bin_storage(const view_t& view) : m_bin_data(view) {}

    /// begin and end of the bin range
    /// @{
    DETRAY_HOST_DEVICE
    auto begin() { return detray::ranges::begin(m_bin_data); }
    DETRAY_HOST_DEVICE
    auto begin() const { return detray::ranges::cbegin(m_bin_data); }
    DETRAY_HOST_DEVICE
    auto end() { return detray::ranges::end(m_bin_data); }
    DETRAY_HOST_DEVICE
    auto end() const { return detray::ranges::cend(m_bin_data); }
    /// @}

    /// @returns the vecmem view of the bin storage
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() -> view_type {
        return detray::get_data(m_bin_data);
    }

    /// @returns the vecmem view of the bin storage - const
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() const -> const_view_type {
        return detray::get_data(m_bin_data);
    }

    private:
    /// Container that holds all bin data when owning or a view into an
    /// externally owned container
    bin_range_t m_bin_data{};
};

}  // namespace detray::detail
