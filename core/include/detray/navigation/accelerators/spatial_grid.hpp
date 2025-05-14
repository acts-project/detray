/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/navigation/accelerators/search_window.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/detail/bin_view.hpp"
#include "detray/utils/grid/detail/grid_bins.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/grid_collection.hpp"

namespace detray {

/// @brief An N-dimensional spacial grid for geometry object searches.
template <concepts::grid grid_t>
class spatial_grid_impl : public grid_t {

    using base_grid = grid_t;

    public:
    using value_type = typename base_grid::value_type;
    using query_type = typename base_grid::point_type;

    /// Find the corresponding (non-)owning grid type
    template <bool is_owning>
    using type = spatial_grid_impl<typename grid_t::template type<is_owning>>;

    /// Use all of the grid constructors
    using base_grid::base_grid;

    /// Construct from existing grid - move
    DETRAY_HOST_DEVICE
    explicit constexpr spatial_grid_impl(base_grid &&gr)
        : base_grid(std::move(gr)) {}

    /// Construct from existing grid - copy
    DETRAY_HOST_DEVICE
    explicit constexpr spatial_grid_impl(base_grid &gr) : base_grid(gr) {}

    /// Find the value of a single bin - const
    ///
    /// @param p is point in the local (bound) frame
    ///
    /// @return the iterable view of the bin content
    DETRAY_HOST_DEVICE decltype(auto) search(const query_type &p) const {
        return this->bin(p);
    }

    /// Find the value of a single bin
    ///
    /// @param p is point in the local (bound) frame
    ///
    /// @return the iterable view of the bin content
    DETRAY_HOST_DEVICE decltype(auto) search(const query_type &p) {
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
    template <concepts::arithmetic window_size_t>
    DETRAY_HOST_DEVICE auto search(
        const query_type &p,
        const search_window<window_size_t, 2> &win_size) const {

        // Return iterable over bins in the search window
        auto bin_ranges = this->axes().bin_ranges(p, win_size);
        auto search_area = axis::detail::bin_view(*this, bin_ranges);

        // Join the respective bins to a single iteration
        return detray::views::join(std::move(search_area));
    }

    /// Interface for the navigator
    template <typename detector_t, typename track_t,
              concepts::arithmetic window_size_t>
    DETRAY_HOST_DEVICE auto search(
        const detector_t &det, const typename detector_t::volume_type &volume,
        const track_t &track, const search_window<window_size_t, 2> &win_size,
        const typename detector_t::geometry_context &ctx) const {

        // Track position in grid coordinates
        const auto &trf = det.transform_store().at(volume.transform(), ctx);
        const auto loc_pos = this->project(trf, track.pos(), track.dir());

        // Grid lookup
        return search(loc_pos, win_size);
    }
};

template <concepts::algebra algebra_t, typename axes_t, typename bin_t,
          template <std::size_t> class serializer_t = simple_serializer,
          typename containers = host_container_types, bool ownership = true>
using spatial_grid = spatial_grid_impl<
    grid_impl<coordinate_axes<axes_t, algebra_t, ownership, containers>, bin_t,
              simple_serializer>>;

/// Accelerator collection specialization for @c detray::spatial_grid_impl
///
/// @todo add container for masks
template <concepts::grid grid_t>
    requires(!spatial_grid_impl<grid_t>::is_owning)
class grid_collection<spatial_grid_impl<grid_t>>
    : public grid_collection<grid_t> {
    // Use a normal grid collection fro the grid related data
    using base_collection = grid_collection<grid_t>;

    public:
    using size_type = typename base_collection::size_type;
    using value_type = spatial_grid_impl<grid_t>;

    /// Use the same constructors as the grid collection
    using base_collection::base_collection;

    /// Create spatial grid acceleration structure from underlying grid data
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const size_type i) const
        -> spatial_grid_impl<grid_t> {
        return spatial_grid_impl<grid_t>(base_collection::operator[](i));
    }

    /// Create spatial grid acceleration structure from underlying grid data
    DETRAY_HOST_DEVICE
    constexpr auto at(const size_type i) const -> spatial_grid_impl<grid_t> {
        return spatial_grid_impl<grid_t>(base_collection::at(i));
    }
};

}  // namespace detray
