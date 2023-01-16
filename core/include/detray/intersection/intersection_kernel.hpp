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
        const scalar mask_tolerance = 0.) const {

        const auto &ctf = contextual_transforms[surface.transform()];

        // Run over the masks that belong to the surface (only one can be hit)
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            if (place_in_collection(mask.intersector()(traj, surface, mask, ctf,
                                                       mask_tolerance),
                                    is_container)) {
                return;
            };
        }
    }

    private:
    template <typename is_container_t>
    bool place_in_collection(
        std::array<typename is_container_t::value_type, 1> &&is,
        is_container_t &intersections) const {
        if (is[0].status == intersection::status::e_inside) {
            intersections.push_back(is[0]);
            return true;
        } else {
            return false;
        }
    }

    template <typename is_container_t>
    bool place_in_collection(
        std::array<typename is_container_t::value_type, 2> &&is,
        is_container_t &intersections) const {
        bool is_valid = false;
        if (is[0].status == intersection::status::e_inside) {
            intersections.push_back(is[0]);
            is_valid = true;
            if (is[1].status == intersection::status::e_inside) {
                intersections.push_back(is[1]);
            }
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
              typename surface_t, typename transform_container_t>
    DETRAY_HOST_DEVICE inline line_plane_intersection<
        surface_t, typename transform_container_t::value_type>
    operator()(const mask_group_t &mask_group, const mask_range_t &mask_range,
               const traj_t &traj, const surface_t &surface,
               const transform_container_t &contextual_transforms,
               const scalar mask_tolerance = 0.) const {

        const auto &ctf = contextual_transforms[surface.transform()];

        // Run over the masks that belong to the surface
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            auto sfi =
                mask.intersector()(traj, surface, mask, ctf, mask_tolerance);
            if (sfi[0].status == intersection::status::e_inside) {
                return sfi[0];
            }
        }

        // return null object if the intersection is not valid anymore
        return {};
    }

    private:
    /*template<typename is_container_t, typename surface_t>
    line_plane_intersection update(std::array<typename
    is_container_t::value_type, 1> &&is, const surface_t &surface) const { if
    (sfi[0].status == intersection::status::e_inside) { sfi[0].surface =
    surface; return sfi[0];
        }
    }

    template<typename is_container_t, typename surface_t>
    line_plane_intersection update(std::array<typename
    is_container_t::value_type, 2> &&is, const surface_t &surface) const { if
    (sfi[0].status == intersection::status::e_inside) { sfi[0].surface =
    surface; return sfi[0];
        }
    }*/
};

}  // namespace detray