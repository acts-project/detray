/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/surface_finders/grid/grid.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray {

/// @brief A collection of grids that can be moved to device.
///
/// @tparam grid_t The type of grid in this collection. Must be non-owning, so
///                that the grid collection can manage the underlying memory.
/// @tparam containers Host or device container types
template <typename grid_t, typename = void>
class grid_collection {};

/// Specialization for @c detray::grid
///
/// @todo refactor this, grid_data and grid_view as detray::ranges::grid_view
template <typename multi_axis_t, typename value_t,
          template <std::size_t> class serializer_t, typename populator_t>
class grid_collection<
    detray::grid<multi_axis_t, value_t, serializer_t, populator_t>,
    std::enable_if_t<not detray::grid<multi_axis_t, value_t, serializer_t,
                                      populator_t>::is_owning,
                     void>> {

    public:
    using grid_type =
        detray::grid<multi_axis_t, value_t, serializer_t, populator_t>;
    using const_grid_type = const grid_type;
    using value_type = grid_type;
    using size_type = dindex;

    /// Backend storage type for the grid
    using bin_storage_type = typename grid_type::bin_storage_type;
    /// Data that the axes keep: bin boundary ranges in the edges container
    using axes_storage_type = typename multi_axis_t::boundary_storage_type;
    /// Contains all bin edges for all axes
    using edges_storage_type = typename multi_axis_t::edges_storage_type;
    template <typename T>
    using vector_type = typename multi_axis_t::template vector_type<T>;

    /// Vecmem based grid collection view type
    using view_type = dmulti_view<dvector_view<size_type>,
                                  detail::get_view_t<bin_storage_type>,
                                  detail::get_view_t<axes_storage_type>,
                                  detail::get_view_t<edges_storage_type>>;

    /// Vecmem based grid collection view type
    using const_view_type =
        dmulti_view<dvector_view<const size_type>,
                    detail::get_view_t<const bin_storage_type>,
                    detail::get_view_t<const axes_storage_type>,
                    detail::get_view_t<const edges_storage_type>>;

    /// Make grid default constructible: Empty grid with empty axis
    grid_collection() = default;

    /// Create empty grid with empty axes from specific vecmem memory resource
    DETRAY_HOST
    explicit grid_collection(vecmem::memory_resource &resource)
        : m_offsets(&resource),
          m_bins(&resource),
          m_axes_data(&resource),
          m_bin_edges(&resource) {}

    /// Create grid colection from existing data - move
    DETRAY_HOST_DEVICE
    grid_collection(vector_type<size_type> &&offs, bin_storage_type &&bins,
                    axes_storage_type &&axes_data, edges_storage_type &&edges)
        : m_offsets(std::move(offs)),
          m_bins(std::move(bins)),
          m_axes_data(std::move(axes_data)),
          m_bin_edges(std::move(edges)) {}

    /// Device-side construction from a vecmem based view type
    template <typename coll_view_t,
              typename std::enable_if_t<detail::is_device_view_v<coll_view_t>,
                                        bool> = true>
    DETRAY_HOST_DEVICE grid_collection(coll_view_t &view)
        : m_offsets(detail::get<0>(view.m_views).m_view),
          m_bins(detail::get<1>(view.m_views).m_view),
          m_axes_data(detail::get<2>(view.m_views).m_view),
          m_bin_edges(detail::get<3>(view.m_views).m_view) {}

    /// @returns the number of grids in the collection - const
    DETRAY_HOST_DEVICE
    constexpr auto ngrids() const -> std::size_t { return m_offsets.size(); }

    /// @returns the offsets for the grids in the bin storage - const
    DETRAY_HOST
    auto offsets() const -> const vector_type<size_type> & { return m_offsets; }
    DETRAY_HOST
    auto offsets() -> vector_type<size_type> & { return m_offsets; }

    /// @returns the underlying bin content storage - const
    DETRAY_HOST
    auto bin_storage() const -> const bin_storage_type & { return m_bins; }
    DETRAY_HOST
    auto bin_storage() -> bin_storage_type & { return m_bins; }

    /// @returns the underlying axis boundary storage - const
    DETRAY_HOST
    auto axes_storage() const -> const axes_storage_type & {
        return m_axes_data;
    }
    DETRAY_HOST
    auto axes_storage() -> axes_storage_type & { return m_axes_data; }

    /// @returns the underlying bin edges storage - const
    DETRAY_HOST
    auto bin_edges_storage() const -> const edges_storage_type & {
        return m_bin_edges;
    }
    DETRAY_HOST
    auto bin_edges_storage() -> edges_storage_type & { return m_bin_edges; }

    /// Create grid from container pointers - const
    DETRAY_HOST_DEVICE
    auto operator[](const size_type i) const -> const_grid_type {
        const size_type axes_offset{grid_type::Dim * i};
        return const_grid_type(
            &m_bins, multi_axis_t(&m_axes_data, &m_bin_edges, axes_offset),
            m_offsets[i]);
    }

    /// Create grid from container pointers - non-const
    DETRAY_HOST_DEVICE
    auto operator[](const size_type i) -> grid_type {
        const size_type axes_offset{grid_type::Dim * i};
        return grid_type(&m_bins,
                         multi_axis_t(&m_axes_data, &m_bin_edges, axes_offset),
                         m_offsets[i]);
    }

    private:
    /// Offsets for the respective grids into the bin storage
    vector_type<size_type> m_offsets{};
    /// Contains the bin content for all grids
    bin_storage_type m_bins{};
    /// Contains the axis boundaries/no. bins for all axes of all grids
    axes_storage_type m_axes_data{};
    /// Contains the bin edges for all grids
    edges_storage_type m_bin_edges{};
};

/// @returns view of a grid collection
template <typename grid_t>
DETRAY_HOST inline typename grid_collection<grid_t>::view_type get_data(
    grid_collection<grid_t> &grid_coll) {
    return {vecmem::get_data(grid_coll.offsets()),
            vecmem::get_data(grid_coll.bin_storage()),
            vecmem::get_data(grid_coll.axes_storage()),
            vecmem::get_data(grid_coll.bin_edges_storage())};
}

/// @returns const view of a grid, for every grid that is passed as a const
/// reference
template <typename grid_t>
DETRAY_HOST inline typename grid_collection<grid_t>::const_view_type get_data(
    const grid_collection<grid_t> &grid_coll) {
    return {vecmem::get_data(grid_coll.offsets()),
            vecmem::get_data(grid_coll.bin_storage()),
            vecmem::get_data(grid_coll.axes_storage()),
            vecmem::get_data(grid_coll.bin_edges_storage())};
}

// template<typename grid_t>
// using grid_container = viewable_collection<grid_collection<grid_t>>;
}  // namespace detray
