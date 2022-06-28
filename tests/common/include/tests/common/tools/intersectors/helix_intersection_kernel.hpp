/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "tests/common/tools/intersectors/helix_cylinder_intersector.hpp"
#include "tests/common/tools/intersectors/helix_plane_intersector.hpp"

namespace detray {

/// A functor to update the closest intersection between helix and
/// surface
struct helix_intersection_update {

    using output_type = line_plane_intersection;

    /// Operator function to update the intersection
    ///
    /// @tparam mask_group_t is the input mask group type found by variadic
    /// unrolling
    /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
    /// @tparam surface_t is the input surface type
    /// @tparam transform_container_t is the input transform store type
    ///
    /// @param mask_group is the input mask group
    /// @param traj is the input trajectory
    /// @param surface is the input surface
    /// @param contextual_transforms is the input transform container
    /// @param edge_tolerance is the tolerance for mask size
    ///
    /// @return the intersection
    ///
    template <typename mask_group_t, typename traj_t, typename surface_t,
              typename transform_container_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const mask_group_t &mask_group, const traj_t &traj,
        const surface_t &surface,
        const transform_container_t &contextual_transforms,
        const scalar edge_tolerance = 0.) const {

        const auto &mask_range = surface.mask_range();
        const auto &ctf = contextual_transforms[surface.transform()];

        // Run over the masks belonged to the surface
        for (const auto &mask : range(mask_group, mask_range)) {

            if constexpr (std::is_same_v<typename mask_group_t::value_type::
                                             intersector_type,
                                         plane_intersector>) {

                auto sfi = std::move(
                    helix_plane_intersector()(traj, mask, ctf, edge_tolerance));

                if (sfi[0].status == intersection::status::e_inside and
                    sfi[0].direction == intersection::direction::e_along) {
                    sfi[0].index = surface.volume();
                    return sfi[0];
                }

            } else if constexpr (std::is_same_v<
                                     typename mask_group_t::value_type::
                                         intersector_type,
                                     cylinder_intersector>) {

                auto sfi = std::move(helix_cylinder_intersector()(
                    traj, mask, ctf, edge_tolerance));

                if (sfi[0].status == intersection::status::e_inside and
                    sfi[0].direction == intersection::direction::e_along) {
                    sfi[0].index = surface.volume();
                    return sfi[0];
                }
            }
        }

        // return null object if the intersection is not valid anymore
        return output_type{};
    }
};

}  // namespace detray