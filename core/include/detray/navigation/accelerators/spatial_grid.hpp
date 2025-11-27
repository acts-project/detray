/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/navigation/accelerators/search_window.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"
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
    using frame_t = typename grid_t::local_frame_type;
    using algebra_t = typename grid_t::algebra_type;
    using mask_t = mask<concentric_cylinder2D, algebra_t, std::uint8_t>;

    public:
    using value_type = typename base_grid::value_type;
    using query_type = typename base_grid::point_type;

    /// Find the corresponding (non-)owning grid type
    template <bool is_owning>
    using type = spatial_grid_impl<typename grid_t::template type<is_owning>>;

    /// Use all of the grid constructors
    using base_grid::base_grid;

    /// Construct from existing grid - move
    ///
    /// @param mask_values additional mask boundary values that cannot be
    /// inferred from the axis spans, e.g. the concentric cylinder radius
    template <typename... Args>
    DETRAY_HOST_DEVICE explicit constexpr spatial_grid_impl(
        base_grid &&gr, Args &&...mask_values)
        : base_grid(std::move(gr)),
          m_mask{get_mask_from_axes(mask_values...)} {}

    /// Construct from existing grid - copy
    template <typename... Args>
    DETRAY_HOST_DEVICE explicit constexpr spatial_grid_impl(
        base_grid &gr, Args &&...mask_values)
        : base_grid(gr), m_mask{get_mask_from_axes(mask_values...)} {}

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

        // Placement of the grid (same as volume)
        const auto &trf = det.transform_store().at(volume.transform(), ctx);

        query_type loc_pos{};
        // For 2-dimensional grids project the track position along the track
        // direction to find the optimal grid bin
        if constexpr (base_grid::dim == 2) {
            // Tangential to the current track parameters
            detray::detail::ray<algebra_t> tangential{track};

            // Intersect the (virtual) reference surface of the grid to find
            // the correct bin
            using intersector_t =
                ray_intersector_impl<frame_t, algebra_t,
                                     intersection::contains_pos>;

            constexpr intersector_t intersector{};
            typename intersector_t::result_type result{};

            if constexpr (concepts::cylindrical<frame_t>) {
                // The cylinder intersector requires the radius from the mask
                result =
                    intersector.point_of_intersection(tangential, trf, m_mask);
            } else {
                result = intersector.point_of_intersection(tangential, trf);
            }

            // Retrieve the closest intersection point
            typename intersector_t::point_type intr_point;
            if constexpr (intersector_t::n_solutions == 1) {
                intr_point = result.point;
            } else {
                // Use the closest intersection
                intr_point = result[0].point;
            }

            // Most intersectors return global positions -> project to grid axes
            if constexpr (std::same_as<typename intersector_t::point_type,
                                       query_type>) {
                loc_pos = intr_point;
            } else if constexpr (std::same_as<frame_t, concentric_cylindrical2D<
                                                           algebra_t>>) {
                // The concentric cylinder intersector only returns the z-pos
                loc_pos = {vector::phi(tangential.pos(result.path)),
                           intr_point[1]};
            } else {
                // Projec the intersection point into the local grid frame
                loc_pos = this->project(trf, intr_point, track.dir());
            }
        } else {
            loc_pos = this->project(trf, track.pos(), track.dir());
        }

        // Grid lookup
        return search(loc_pos, win_size);
    }

    private:
    /// @returns a mask that has boundaries which match the grid axis spans
    template <typename... Args>
    DETRAY_HOST_DEVICE mask_t get_mask_from_axes(Args &&...mask_values) {

        constexpr auto inv_vol_link{
            detray::detail::invalid_value<std::uint8_t>()};

        /// @param mask_values contains the cylinder radius
        if constexpr (concepts::cylindrical<frame_t>) {
            static_assert(sizeof...(Args) == 1,
                          "Spatial cylinder grid: Only the cylinder radius "
                          "needs to be provided externally");

            constexpr auto inf{std::numeric_limits<dscalar<algebra_t>>::max()};
            return mask_t{inv_vol_link, mask_values..., -inf, inf};
        } else {
            /// @TODO: Implement axes to mask conversion for the other shapes
            return {};
        }
    }

    /// Struct that contains the grid's data state
    mask_t m_mask{};
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
    using frame_t = typename grid_t::local_frame_type;

    public:
    using size_type = typename base_collection::size_type;
    using value_type = spatial_grid_impl<grid_t>;

    /// Use the same constructors as the grid collection
    using base_collection::base_collection;

    /// Create spatial grid acceleration structure from underlying grid data
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const size_type i) const
        -> spatial_grid_impl<grid_t> {
        if constexpr (concepts::cylindrical<frame_t>) {
            return spatial_grid_impl<grid_t>(base_collection::operator[](i),
                                             100.f);
        } else {
            return spatial_grid_impl<grid_t>(base_collection::operator[](i));
        }
    }

    /// Create spatial grid acceleration structure from underlying grid data
    DETRAY_HOST_DEVICE
    constexpr auto at(const size_type i) const -> spatial_grid_impl<grid_t> {
        if constexpr (concepts::cylindrical<frame_t>) {
            return spatial_grid_impl<grid_t>(base_collection::at(i), 100.f);
        } else {
            return spatial_grid_impl<grid_t>(base_collection::at(i));
        }
    }
};

}  // namespace detray
