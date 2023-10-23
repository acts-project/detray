/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/surface_finders/grid/detail/axis_helpers.hpp"
#include "detray/tools/bin_association.hpp"

// System include(s)
#include <cassert>
#include <vector>

namespace detray::detail {

/// Fill a surface grid using a local bin index of the grid and a payload.
///
/// @param grid the grid that should be filled
/// @param det the detector from which to get the surface placements
/// @param vol the volume the grid belongs to
/// @param ctx the geometry context
struct fill_by_bin {

    /// Single piece of data for a grid bin
    template <typename grid_t>
    struct bin_data {
        /// Bin index on the grid axes
        n_axis::multi_bin<grid_t::Dim> local_bin_idx;
        /// Single element of the bin content, which can be a collection of grid
        /// values
        typename grid_t::value_type single_element;
    };

    template <typename grid_t, typename volume_t, typename surface_container_t,
              typename mask_container, typename transform_container,
              typename context_t, typename... Args>
    DETRAY_HOST auto operator()(grid_t &grid, const volume_t &,
                                const surface_container_t &,
                                const mask_container &, const context_t,
                                std::vector<bin_data<grid_t>> &bins) const
        -> void {
        for (const bin_data<grid_t> &bd : bins) {
            grid.populate(bd.local_bin_idx, bd.single_element);
        }
    }
};

/// Fill a surface grid using the surface translation.
///
/// @param grid the grid that should be filled
/// @param det the detector from which to get the surface placements
/// @param vol the volume the grid belongs to
/// @param ctx the geometry context
struct fill_by_pos {

    template <typename grid_t, typename volume_t, typename surface_container_t,
              typename mask_container, typename transform_container,
              typename context_t, typename... Args>
    DETRAY_HOST auto operator()(grid_t &grid, const volume_t &vol,
                                const surface_container_t &surfaces,
                                const transform_container &transforms,
                                const mask_container & /*masks*/,
                                const context_t ctx, Args &&...) const -> void {

        // Fill the volumes surfaces into the grid
        for (const auto &sf : surfaces) {

            // no portals in grids allowed
            if (sf.volume() == vol.index() && sf.is_sensitive()) {

                const auto &sf_trf = transforms.at(sf.transform(), ctx);
                const auto &t = sf_trf.translation();

                // transform to axis coordinate system
                const auto loc_pos =
                    grid.global_to_local(vol.transform(), t, t);

                // Populate
                grid.populate(loc_pos, sf);
            }
        }
    }
};

/// Fill a grid surface finder by bin association.
///
/// @param grid the grid that should be filled
/// @param det the detector from which to get the surface placements
/// @param vol the volume the grid belongs to
/// @param ctx the geometry context
struct bin_associator {

    template <typename detector_t, typename volume_type, typename grid_t,
              typename... Args>
    DETRAY_HOST auto operator()(grid_t &grid, detector_t &det,
                                const volume_type &vol,
                                const typename detector_t::geometry_context ctx,
                                Args &&...) const -> void {
        this->operator()(grid, det.surfaces(vol), det.mask_store(),
                         det.transform_store(), ctx);
    }

    template <typename grid_t, typename volume_t, typename surface_container_t,
              typename mask_container, typename transform_container,
              typename context_t, typename... Args>
    DETRAY_HOST auto operator()(grid_t &grid, const volume_t &,
                                const surface_container_t &surfaces,
                                const transform_container &transforms,
                                const mask_container &masks,
                                const context_t ctx, Args &&...) const -> void {
        // Fill the surfaces into the grid by matching their contour onto the
        // grid bins
        bin_association(ctx, surfaces, transforms, masks, grid, {0.1f, 0.1f},
                        false);
    }
};

}  // namespace detray::detail
