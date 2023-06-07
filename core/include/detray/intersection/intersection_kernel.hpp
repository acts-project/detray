/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include <iostream>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/utils/ranges.hpp"
namespace detray {

/// A functor to add all valid intersections between the trajectory and surface
struct intersection_initialize {

    /// Operator function to initalize intersections
    ///
    /// @tparam mask_group_t is the input mask group type found by variadic
    /// unrolling
    /// @tparam is_container_t is the intersection container type
    /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
    /// @tparam surface_t is the input surface type
    /// @tparam transform_container_t is the input transform store type
    ///
    /// @param mask_group is the input mask group
    /// @param is_container is the intersection container to be filled
    /// @param traj is the input trajectory
    /// @param surface is the input surface
    /// @param contextual_transforms is the input transform container
    /// @param mask_tolerance is the tolerance for mask size
    ///
    /// @return the number of valid intersections
    template <typename mask_group_t, typename mask_range_t,
              typename is_container_t, typename traj_t, typename surface_t,
              typename transform_container_t>
    DETRAY_HOST_DEVICE inline void operator()(
        const mask_group_t &mask_group, const mask_range_t &mask_range,
        is_container_t &is_container, const traj_t &traj,
        const surface_t &surface,
        const transform_container_t &contextual_transforms,
        const scalar mask_tolerance = 0.f) const {

        using intersection_t = typename is_container_t::value_type;

        const auto &ctf = contextual_transforms[surface.transform()];

        // Run over the masks that belong to the surface (only one can be hit)
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            if (place_in_collection(
                    mask.template intersector<intersection_t>()(
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
        // std::cout << sfi;
        bool is_inside = (sfi.status == intersection::status::e_inside);
        if (is_inside) {
            intersections.push_back(sfi);
        }
        return is_inside;
    }

    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        std::array<typename is_container_t::value_type, 2> &&solutions,
        is_container_t &intersections) const {
        bool is_valid = false;
        for (auto &sfi : solutions) {
            bool is_inside = (sfi.status == intersection::status::e_inside);
            if (is_inside) {
                intersections.push_back(sfi);
            }
            is_valid |= is_inside;
        }
        return is_valid;
    }
};

/// A functor to update the closest intersection between the trajectory and
/// surface
struct intersection_update {

    /// Operator function to update the intersection
    ///
    /// @tparam mask_group_t is the input mask group type found by variadic
    /// unrolling
    /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
    /// @tparam surface_t is the input surface type
    /// @tparam transform_container_t is the input transform store type
    ///
    /// @param mask_group is the input mask group
    /// @param mask_range is the range of masks in the group that belong to the
    ///                   surface
    /// @param traj is the input trajectory
    /// @param surface is the input surface
    /// @param contextual_transforms is the input transform container
    /// @param mask_tolerance is the tolerance for mask size
    ///
    /// @return the intersection
    template <typename mask_group_t, typename mask_range_t, typename traj_t,
              typename intersection_t, typename transform_container_t>
    DETRAY_HOST_DEVICE inline bool operator()(
        const mask_group_t &mask_group, const mask_range_t &mask_range,
        const traj_t &traj, intersection_t &sfi,
        const transform_container_t &contextual_transforms,
        const scalar mask_tolerance = 0.f) const {

        const auto &ctf = contextual_transforms[sfi.surface.transform()];

        // Run over the masks that belong to the surface
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            mask.template intersector<intersection_t>().update(
                traj, sfi, mask, ctf, mask_tolerance);

            if (sfi.status == intersection::status::e_inside) {
                return true;
            }
        }

        return false;
    }
};

}  // namespace detray