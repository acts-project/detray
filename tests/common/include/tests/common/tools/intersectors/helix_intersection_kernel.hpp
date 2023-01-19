/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include <iostream>

#include "detray/intersection/detail/trajectories.hpp"
#include "detray/utils/ranges.hpp"
#include "tests/common/tools/intersectors/helix_cylinder_intersector.hpp"
#include "tests/common/tools/intersectors/helix_plane_intersector.hpp"

namespace detray {

/// A functor to update the closest intersection between helix and
/// surface
struct helix_intersection_initialize {

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
    /// @param mask_tolerance is the tolerance for mask size
    ///
    /// @return the intersection
    template <typename mask_group_t, typename mask_range_t,
              typename is_container_t, typename traj_t, typename surface_t,
              typename transform_container_t>
    DETRAY_HOST_DEVICE inline void operator()(
        const mask_group_t &mask_group, const mask_range_t &mask_range,
        is_container_t &is_container, const traj_t &traj,
        const surface_t &surface,
        const transform_container_t &contextual_transforms,
        const scalar mask_tolerance = 0.f) const {

        using mask_t = typename mask_group_t::value_type;
        using transform3_t = typename transform_container_t::value_type;

        const auto &ctf = contextual_transforms[surface.transform()];

        // Run over the masks that belong to the surface
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            if (place_in_collection(
                    helix_intersector<transform3_t, mask_t>()(
                        traj, surface, mask, ctf, mask_tolerance),
                    is_container)) {
                return;
            };
        }
    }

    private:
    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        typename is_container_t::value_type &&sfi,
        is_container_t &intersections) const {
        if (sfi.status == intersection::status::e_inside) {
            intersections.push_back(sfi);
            return true;
        } else {
            return false;
        }
    }

    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        std::array<typename is_container_t::value_type, 2> &&solutions,
        is_container_t &intersections) const {
        bool is_valid = false;
        for (auto &sfi : solutions) {
            if (sfi.status == intersection::status::e_inside) {
                intersections.push_back(sfi);
                is_valid = true;
            }
        }
        return is_valid;
    }
};

}  // namespace detray