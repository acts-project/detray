/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// system include
#include <cmath>
#include <type_traits>

// detray include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/helix_cylinder_intersector.hpp"
#include "detray/intersection/helix_plane_intersector.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/intersection/ray_cylinder_intersector.hpp"
#include "detray/intersection/ray_plane_intersector.hpp"
#include "detray/utils/enumerate.hpp"

namespace detray {

namespace helix {

template <typename surface_t, typename transform_container_t,
          typename mask_container_t>
DETRAY_HOST_DEVICE static inline auto intersect(
    const detail::helix &h, surface_t &surface,
    const transform_container_t &contextual_transforms,
    const mask_container_t &masks, const scalar epsilon);

}

/// @brief struct that holds functionality to shoot a parametrzed particle
/// trajectory through a detector.
///
/// Records intersections with every detector surface along the trajectory.
struct particle_gun {

    using intersection_type = line_plane_intersection;

    /// Intersect all surfaces in a detector with a given ray.
    ///
    /// @param detector the detector.
    /// @param traj the trajectory to be shot through the detector.
    /// @param epsilon numerical precision
    ///
    /// @return a sorted vector of volume indices with the corresponding
    ///         intersections of the surfaces that were encountered.
    template <typename detector_t, typename trajectory_t>
    DETRAY_HOST_DEVICE inline static auto shoot_particle(
        const detector_t &detector, const trajectory_t &traj,
        const scalar epsilon = 1e-5) {

        std::vector<std::pair<dindex, intersection_type>> intersection_record;

        // Loop over all surfaces in the detector
        for (const auto &volume : detector.volumes()) {
            for (const auto &[sf_idx, sf] :
                 enumerate(detector.surfaces(), volume)) {

                // Retrieve candidate from the surface
                intersection_type sfi;
                if constexpr (std::is_same_v<trajectory_t, detail::helix>) {
                    // Call helix specific version instead of the intersection
                    // kernel
                    sfi = helix::intersect(traj, sf, detector.transform_store(),
                                           detector.mask_store(), epsilon);
                } else {
                    // Call detray intersection kernel
                    sfi =
                        detray::intersect(traj, sf, detector.transform_store(),
                                          detector.mask_store());
                }

                // Candidate is invalid if it oversteps too far (this is neg!)
                if (sfi.path < traj.overstep_tolerance()) {
                    continue;
                }
                // Accept if inside
                if (sfi.status == intersection::status::e_inside) {
                    // Volume the candidate belongs to
                    sfi.index = volume.index();
                    intersection_record.emplace_back(sf_idx, sfi);
                }
            }
        }

        // Sort intersections by distance to origin of the trajectory
        auto sort_path = [&](std::pair<dindex, intersection_type> a,
                             std::pair<dindex, intersection_type> b) -> bool {
            return (a.second < b.second);
        };
        std::sort(intersection_record.begin(), intersection_record.end(),
                  sort_path);

        return intersection_record;
    }
};

namespace helix {

/// Helix version of the intersection unrolling. Calls the helix intersectors
/// instead of the mask's native intersectors. This way, the geometry does not
/// need to be recompiled. See documentation of the detray
/// @c intersection_kernel
template <typename mask_defs, typename transform_t, typename mask_container_t,
          typename mask_range_t, unsigned int first_mask_id,
          unsigned int... remaining_mask_ids>
DETRAY_HOST_DEVICE static inline auto unroll_intersect(
    const detail::helix &h, const transform_t &ctf,
    const mask_container_t &masks, const mask_range_t &mask_range,
    const typename mask_container_t::id_type mask_id, dindex volume_index,
    const scalar epsilon,
    std::integer_sequence<unsigned int, first_mask_id, remaining_mask_ids...>
    /*available_ids*/) {

    using intersection_type = line_plane_intersection;

    // Pick the first one for interseciton
    if (mask_id == first_mask_id) {
        // Get the mask id that was found
        constexpr auto id = mask_container_t::to_id(first_mask_id);
        auto &mask_group = masks.template group<id>();

        // Check all masks of this surface for intersection with the helix.
        // In order to not recompile the geometry, use external intersectors
        for (const auto &mask : range(mask_group, mask_range)) {
            using mask_t = typename mask_defs::template get_type<id>::type;

            intersection_type sfi;
            // Make sure helix_intersectors are only called for the correct
            // surface type
            if constexpr (std::is_same_v<typename mask_t::intersector_type,
                                         ray_plane_intersector>) {
                // TODO: Might not work with all mask types
                typename mask_t::mask_tolerance tol{epsilon};
                sfi = helix_plane_intersector::intersect(ctf, h, mask, tol);
            }
            if constexpr (std::is_same_v<typename mask_t::intersector_type,
                                         ray_cylinder_intersector>) {
                typename mask_t::mask_tolerance tol{epsilon, epsilon};
                sfi = helix_cylinder_intersector::intersect(ctf, h, mask, tol);
            }
            // Only keep intersection along the direction, like the ray does
            if (sfi.status == intersection::status::e_inside and
                sfi.direction == intersection::direction::e_along) {
                sfi.index = volume_index;
                return sfi;
            }
        }
    }

    // The reduced integer sequence
    std::integer_sequence<unsigned int, remaining_mask_ids...> remaining;

    // Unroll as long as you have at least 1 entries
    if constexpr (remaining.size() >= 1) {
        return (unroll_intersect<mask_defs>(h, ctf, masks, mask_range, mask_id,
                                            volume_index, epsilon, remaining));
    }

    // No intersection was found
    return intersection_type{};
}

/// Start helix intersection unrolling. See documentation of the detray
/// @c intersection_kernel
template <typename surface_t, typename transform_container_t,
          typename mask_container_t>
DETRAY_HOST_DEVICE static inline auto intersect(
    const detail::helix &h, surface_t &surface,
    const transform_container_t &contextual_transforms,
    const mask_container_t &masks, const scalar epsilon) {
    // Gather all information to perform intersections
    const auto &ctf = contextual_transforms[surface.transform()];
    const auto volume_index = surface.volume();
    const auto mask_id = surface.mask_type();
    const auto &mask_range = surface.mask_range();

    // Unroll the intersection depending on the mask container size
    using mask_defs = typename surface_t::mask_defs;

    return unroll_intersect<mask_defs>(
        h, ctf, masks, mask_range, mask_id, volume_index, epsilon,
        std::make_integer_sequence<unsigned int, mask_defs::n_types>{});
}

}  // namespace helix

}  // namespace detray
