/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
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
class grid_collection {
    grid_collection() = delete;
    grid_collection(const grid_collection &) = delete;
    grid_collection(grid_collection &&) = delete;
};

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
    using const_grid_type =
        const detray::grid<const multi_axis_t, const value_t, serializer_t,
                           populator_t>;
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

    using buffer_type = dmulti_buffer<dvector_buffer<size_type>,
                                      detail::get_buffer_t<bin_storage_type>,
                                      detail::get_buffer_t<axes_storage_type>,
                                      detail::get_buffer_t<edges_storage_type>>;

    /// Make grid default constructible: Empty grid with empty axis
    grid_collection() = default;

    /// Create empty grid with empty axes from specific vecmem memory resource
    DETRAY_HOST
    explicit grid_collection(vecmem::memory_resource *resource)
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
    template <typename coll_view_t,
              typename std::enable_if_t<detail::is_device_view_v<coll_view_t>,
                                        bool> = true>
    DETRAY_HOST_DEVICE grid_collection(coll_view_t &view)
        : m_offsets(detail::get<0>(view.m_view)),
          m_bins(detail::get<1>(view.m_view)),
          m_axes_data(detail::get<2>(view.m_view)),
          m_bin_edges(detail::get<3>(view.m_view)) {}

    /// @returns the number of grids in the collection - const
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> dindex {
        return static_cast<dindex>(m_offsets.size());
    }

    /// @returns an iterator that points to the first grid
    /// @note Not implemented!
    DETRAY_HOST_DEVICE
    constexpr auto begin() noexcept -> bool { return true; }

    /// @returns an iterator that points to the coll. end
    /// @note Not implemented!
    DETRAY_HOST_DEVICE
    constexpr auto end() noexcept -> bool { return false; }

    /// @returns the number of grids in the collection - const
    DETRAY_HOST_DEVICE
    constexpr auto empty() const noexcept -> bool { return m_offsets.empty(); }

    /// @brief Resize the underlying containers
    /// @note Not defined! The amount of memory can differ for every grid
    DETRAY_HOST_DEVICE
    constexpr void resize(std::size_t) noexcept { /*Not defined*/
    }

    /// @brief Reserve memory
    /// @note Not defined! The amount of memory can differ for every grid
    DETRAY_HOST_DEVICE
    constexpr void reserve(std::size_t) noexcept { /*Not defined*/
    }

    /// Removes all data from the grid collection containers
    DETRAY_HOST_DEVICE
    constexpr void clear() noexcept {
        m_offsets.clear();
        m_bins.clear();
        m_axes_data.clear();
        m_bin_edges.clear();
    }

    /// Insert a number of grids
    /// @note Not defined! There is no grid iterator implementation
    template <typename... Args>
    DETRAY_HOST_DEVICE constexpr void insert(Args &&...) noexcept {
        /*Not defined*/
    }

    /// @returns the offsets for the grids in the bin storage - const
    DETRAY_HOST_DEVICE
    constexpr auto offsets() const -> const vector_type<size_type> & {
        return m_offsets;
    }

    /// @returns the underlying bin content storage - const
    DETRAY_HOST_DEVICE
    constexpr auto bin_storage() const -> const bin_storage_type & {
        return m_bins;
    }

    /// @returns the underlying axis boundary storage - const
    DETRAY_HOST_DEVICE
    constexpr auto axes_storage() const -> const axes_storage_type & {
        return m_axes_data;
    }

    /// @returns the underlying bin edges storage - const
    DETRAY_HOST_DEVICE
    constexpr auto bin_edges_storage() const -> const edges_storage_type & {
        return m_bin_edges;
    }

    /// Create grid from container pointers - const
    DETRAY_HOST_DEVICE
    auto operator[](const size_type i) const -> grid_type {
        const size_type axes_offset{grid_type::Dim * i};
        return grid_type(&m_bins,
                         multi_axis_t(&m_axes_data, &m_bin_edges, axes_offset),
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

    /// @returns a vecmem view on the grid collection data - non-const
    DETRAY_HOST auto get_data() -> view_type {
        return view_type{detray::get_data(m_offsets), detray::get_data(m_bins),
                         detray::get_data(m_axes_data),
                         detray::get_data(m_bin_edges)};
    }

    /// @returns a vecmem view on the grid collection data - const
    DETRAY_HOST
    auto get_data() const -> const_view_type {
        return const_view_type{
            detray::get_data(m_offsets), detray::get_data(m_bins),
            detray::get_data(m_axes_data), detray::get_data(m_bin_edges)};
    }

    /// Add a new grid @param gr to the collection.
    /// @note this takes a data owning grid to transcribe the data from.
    DETRAY_HOST constexpr auto push_back(
        const typename grid_type::template type<true> &gr) noexcept(false)
        -> void {
        // Current offset into the global bin storage for the new grid
        m_offsets.push_back(static_cast<size_type>(m_bins.size()));

        // Add the bins of the new grid to the collection
        const auto *grid_bins = gr.data().bin_data();
        m_bins.insert(m_bins.end(), grid_bins->begin(), grid_bins->end());

        // Add the axes data of the new grid to the collection
        // (how to lookup the axis bin edges)
        const auto *axes_data = gr.axes().data().axes_data();
        m_axes_data.insert(m_axes_data.end(), axes_data->begin(),
                           axes_data->end());

        // Current offset into the global bin edges storage
        dindex bin_edges_offset = static_cast<dindex>(m_bin_edges.size());
        // The binning types of the axes
        const auto binnings = get_binning(gr.axes());

        // Update the bin edges index range for the axes in the grid collection
        for (std::size_t i = m_axes_data.size() - grid_type::Dim;
             i < m_axes_data.size(); ++i) {
            auto &bin_entry_range = m_axes_data[i];
            bin_entry_range[0] += bin_edges_offset;
            // If the axis has irregular binning, the second entry is the index
            // of the last bin edge value instead of the number of bins
            if (binnings[i] == n_axis::binning::e_irregular) {
                bin_entry_range[1] += bin_edges_offset;
            }
        }

        // Add the bin edges of the new grid to the collection
        const auto *bin_edges = gr.axes().data().edges();

        m_bin_edges.insert(m_bin_edges.end(), bin_edges->begin(),
                           bin_edges->end());
    }

    private:
    /// @returns an array that contians the binning enum values of the
    /// corresponding single axes in @param axes
    template <bool ownership, typename local_frame_t, typename... axis_ts>
    static auto get_binning(
        const n_axis::multi_axis<ownership, local_frame_t, axis_ts...> &axes) {

        // Serialize every single axis and construct array from their payloads
        std::array<n_axis::binning, sizeof...(axis_ts)> binning_types{
            get_binning(axes.template get_axis<axis_ts>())...};

        return binning_types;
    }

    /// @returns the binning enum value of a single axis @param axis
    template <typename bounds_t, typename binning_t>
    static auto get_binning(
        const n_axis::single_axis<bounds_t, binning_t> &axis) {
        return axis.binning();
    }

    /// Offsets for the respective grids into the bin storage
    vector_type<size_type> m_offsets{};
    /// Contains the bin content for all grids
    bin_storage_type m_bins{};
    /// Contains the axis boundaries/no. bins for all axes of all grids
    axes_storage_type m_axes_data{};
    /// Contains the bin edges for all grids
    edges_storage_type m_bin_edges{};
};

}  // namespace detray
