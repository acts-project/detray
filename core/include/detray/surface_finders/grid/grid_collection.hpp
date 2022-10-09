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

    /// Vecmem based grid view type
    using view_type = dmulti_view<detail::get_view_t<bin_storage_type>,
                                  detail::get_view_t<axes_storage_type>,
                                  detail::get_view_t<edges_storage_type>,
                                  dvector_view<size_type>>;

    /// Make grid default constructible: Empty grid with empty axis
    grid_collection() = default;

    /// Create empty grid with empty axes from specific vecmem memory resource
    DETRAY_HOST
    explicit grid_collection(vecmem::memory_resource &resource)
        : m_offsets(resource),
          m_bins(resource),
          m_axes_data(resource),
          m_bin_edges(resource) {}

    /// Create grid colection from existing data - move
    DETRAY_HOST_DEVICE
    grid_collection(vector_type<size_type> &&offs, bin_storage_type &&bins,
                    axes_storage_type &&axes_data, edges_storage_type &&edges)
        : m_offsets(std::move(offs)),
          m_bins(std::move(bins)),
          m_axes_data(std::move(axes_data)),
          m_bin_edges(std::move(edges)) {}

    /// Device-side construction from a vecmem based view type
    template <typename grid_view_t,
              typename std::enable_if_t<detail::is_device_view_v<grid_view_t>,
                                        bool> = true>
    DETRAY_HOST_DEVICE grid_collection(grid_view_t &view)
        : m_offsets(detail::get<2>(view.m_views)),
          m_bins(detail::get<0>(view.m_views)),
          m_axes_data(detail::get<1>(view.m_views)),
          m_bin_edges(detail::get<2>(view.m_views)) {}

    /// @returns the underlying bin content storage - const
    DETRAY_HOST_DEVICE
    auto bin_data() const -> const bin_storage_type & { return m_bins; }

    /// @returns the underlying bin value storage - non-const for vecmem
    // TODO: Don't do
    DETRAY_HOST_DEVICE
    auto bin_data() -> bin_storage_type & { return m_bins; }

    /// @returns the underlying axis boundary storage - const
    DETRAY_HOST_DEVICE
    auto axes_data() const -> const axes_storage_type & { return m_axes_data; }

    /// @returns the underlying axis boundary storage - non-const for vecmem
    // TODO: Don't do
    DETRAY_HOST_DEVICE
    auto axes_data() -> axes_storage_type & { return m_axes_data; }

    /// @returns the underlying bin edges storage - const
    DETRAY_HOST_DEVICE
    auto bin_edges_data() const -> const edges_storage_type & {
        return m_bin_edges;
    }

    /// @returns the underlying bin edges storage - non-const for vecmem
    // TODO: Don't do
    DETRAY_HOST_DEVICE
    auto bin_edges_data() -> edges_storage_type & { return m_bin_edges; }

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

/// @returns const view of a grid, for every grid that is passed as a const
/// reference
/*template <typename grid_t, typename containers>
inline typename grid_collection<grid_t, containers>::view_type
get_data(grid_collection<grid_t, containers> &g_coll) {
    return {vecmem::get_data(*g.data().bin_data()), detray::get_data(g.axes())};
}

/// @returns const view of a grid, for every grid that is passed as a const
/// reference
template <typename grid_t, typename containers>
inline typename grid_collection<const grid_t, containers>::view_type
get_data(const grid_collection<grid_t, containers> &g_coll) {

    // add const to grid value type
    auto const_g = reinterpret_cast<
        const grid<multi_axis_t, const value_t, serializer_t, populator_t> &>(
        g);

    return {vecmem::get_data(*const_g.data().bin_data()),
            detray::get_data(const_g.axes())};
}*/

// template<typename grid_t, typename containers>
// using grid_container = viewable_collection<grid_collection<grid_t,
// containers>>;
}  // namespace detray
