/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/detail/grid_bins.hpp"
#include "detray/utils/grid/detail/bin_view.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/grid_collection.hpp"

namespace detray {

/// @brief An N-dimensional spacial grid for geometry object searches.
///
/// @tparam value_t Type of values contained in the grid
template <typename axes_t, typename bin_t,
          template <std::size_t> class serializer_t = simple_serializer>
class accelerator_grid_impl : public grid_impl<axes_t, bin_t, serializer_t> {

    using base_grid = grid_impl<axes_t, bin_t, serializer_t>;

    public:
    /// Adopt type definitions and static variables
    /// @{
    template <typename neighbor_t>
    using neighborhood_type = darray<neighbor_t, base_grid::dim>;

    static constexpr unsigned int dim{base_grid::dim};
    static constexpr bool is_owning{base_grid::is_owning};

    using value_type = typename base_grid::value_type;
    using algebra_type = typename base_grid::algebra_type;
    using scalar_type = typename base_grid::scalar_type;
    using point_type = typename base_grid::point_type;

    using bin_type = typename base_grid::bin_type;
    using glob_bin_index = typename base_grid::glob_bin_index;
    using loc_bin_index = typename base_grid::loc_bin_index;
    using axes_type = typename base_grid::axes_type;
    using bin_storage = typename base_grid::bin_storage;
    using bin_container_type = typename base_grid::bin_container_type;

    using view_type = typename base_grid::view_type;
    using const_view_type = typename base_grid::const_view_type;
    using buffer_type = typename base_grid::buffer_type;

    template <std::size_t DIM>
    using serializer_type = typename base_grid::template serializer_type<DIM>;

    /// Find the corresponding (non-)owning grid type
    template <bool owning>
    using type =
        accelerator_grid_impl<typename axes_t::template type<owning>, bin_t, serializer_t>;
    /// @}

    /// Use all of the grid constructors
    using base_grid::base_grid;

    /// Interface for the navigator
    template <typename detector_t, typename track_t, typename config_t>
    DETRAY_HOST_DEVICE auto search(
        const detector_t &det, const typename detector_t::volume_type &volume,
        const track_t &track, const config_t &cfg,
        const typename detector_t::geometry_context &ctx) const {

        // Track position in grid coordinates
        const auto &trf = det.transform_store().at(volume.transform(), ctx);
        const auto loc_pos = this->project(trf, track.pos(), track.dir());

        // Grid lookup
        return search(loc_pos, cfg.search_window);
    }

    /// Find the value of a single bin - const
    ///
    /// @param p is point in the local (bound) frame
    ///
    /// @return the iterable view of the bin content
    DETRAY_HOST_DEVICE decltype(auto) search(const point_type &p) const {
        return this->bin(p);
    }

    /// Find the value of a single bin
    ///
    /// @param p is point in the local (bound) frame
    ///
    /// @return the iterable view of the bin content
    DETRAY_HOST_DEVICE decltype(auto) search(const point_type &p) {
        return this->bin(p);
    }

    /// @brief Return a neighborhood of values from the grid
    ///
    /// The lookup is done with a search window around the bin
    ///
    /// @param p is point in the local frame
    /// @param win_size size of the binned/scalar search window
    ///
    /// @return the sequence of values
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE auto search(
        const point_type &p, const darray<neighbor_t, 2> &win_size) const {

        // Return iterable over bins in the search window
        auto search_window = this->axes().bin_ranges(p, win_size);
        auto search_area = axis::detail::bin_view(*this, search_window);

        // Join the respective bins to a single iteration
        return detray::views::join(std::move(search_area));
    }
};

/// Type alias for easier construction
template <concepts::algebra algebra_t, typename axes_t, typename bin_t,
          template <std::size_t> class serializer_t = simple_serializer,
          typename containers = host_container_types, bool ownership = true>
using accelerator_grid =
    accelerator_grid_impl<coordinate_axes<axes_t, algebra_t, ownership, containers>, bin_t,
              simple_serializer>;

}  // namespace detray