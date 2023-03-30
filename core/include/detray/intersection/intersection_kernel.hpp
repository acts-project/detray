/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

/** A functor to add all valid intersections between the trajectory and surface
 */
struct intersection_initialize {

    /** Operator function to initalize intersections
     *
     * @tparam mask_group_t is the input mask group type found by variadic
     * unrolling
     * @tparam is_container_t is the intersection container type
     * @tparam traj_t is the input trajectory type (e.g. ray or helix)
     * @tparam surface_t is the input surface type
     * @tparam transform_container_t is the input transform store type
     *
     * @param mask_group is the input mask group
     * @param is_container is the intersection container to be filled
     * @param traj is the input trajectory
     * @param surface is the input surface
     * @param contextual_transforms is the input transform container
     * @param mask_tolerance is the tolerance for mask size
     *
     * @return the number of valid intersections
     */
    template <typename mask_group_t, typename mask_range_t,
              typename is_container_t, typename traj_t, typename surface_t,
              typename transform_container_t>
    DETRAY_HOST_DEVICE inline std::size_t operator()(
        const mask_group_t &mask_group, const mask_range_t &mask_range,
        is_container_t &is_container, const traj_t &traj,
        const surface_t &surface,
        const transform_container_t &contextual_transforms,
        const scalar mask_tolerance = 0.f) const {

        std::size_t count{0u};

        const auto &ctf = contextual_transforms[surface.transform()];

        // Run over the masks belonged to the surface
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            auto sfi = std::move(mask.intersector()(
                traj, mask, ctf, mask_tolerance, traj.overstep_tolerance()));

            for (auto &is : sfi) {
                if (is.status == intersection::status::e_inside &&
                    is.path >= traj.overstep_tolerance()) {
                    is.barcode = surface.barcode();
                    is_container.push_back(is);
                    count++;
                }
            }

            if (count > 0u) {
                return count;
            }
        }

        return count;
    }
};

/** A functor to update the closest intersection between the trajectory and
 * surface
 */
struct intersection_update {

    /** Operator function to update the intersection
     *
     * @tparam mask_group_t is the input mask group type found by variadic
     * unrolling
     * @tparam traj_t is the input trajectory type (e.g. ray or helix)
     * @tparam surface_t is the input surface type
     * @tparam transform_container_t is the input transform store type
     *
     * @param mask_group is the input mask group
     * @param mask_range is the range of masks in the group that belong to the
     *                   surface
     * @param traj is the input trajectory
     * @param surface is the input surface
     * @param contextual_transforms is the input transform container
     * @param mask_tolerance is the tolerance for mask size
     *
     * @return the intersection
     */
    template <typename mask_group_t, typename mask_range_t, typename traj_t,
              typename surface_t, typename transform_container_t>
    DETRAY_HOST_DEVICE inline line_plane_intersection operator()(
        const mask_group_t &mask_group, const mask_range_t &mask_range,
        const traj_t &traj, const surface_t &surface,
        const transform_container_t &contextual_transforms,
        const scalar mask_tolerance = 0.f) const {

        const auto &ctf = contextual_transforms[surface.transform()];

        // Run over the masks belonging to the surface
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            auto sfi = std::move(mask.intersector()(
                traj, mask, ctf, mask_tolerance, traj.overstep_tolerance()));

            if (sfi[0].status == intersection::status::e_inside &&
                sfi[0].path >= traj.overstep_tolerance()) {
                sfi[0].barcode = surface.barcode();
                return sfi[0];
            }
        }

        // return null object if the intersection is not valid anymore
        return {};
    }
};

}  // namespace detray