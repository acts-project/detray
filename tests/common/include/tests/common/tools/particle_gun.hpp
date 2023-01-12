/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// system include
#include <cmath>
#include <type_traits>

// detray include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/utils/ranges.hpp"
#include "tests/common/tools/intersectors/helix_intersection_kernel.hpp"

namespace detray {

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
        const detector_t &detector, const trajectory_t &traj) {

        std::vector<std::pair<dindex, intersection_type>> intersection_record;

        using helix_type =
            detail::helix<typename trajectory_t::transform3_type>;

        // Loop over all surfaces in the detector
        const auto &mask_store = detector.mask_store();
        const auto &tf_store = detector.transform_store();

        std::vector<intersection_type> intersections{};

        for (const auto &volume : detector.volumes()) {
            for (const auto &sf : detector.surfaces(volume)) {

                // Retrieve candidate(s) from the surface
                if constexpr (std::is_same_v<trajectory_t, helix_type>) {
                    mask_store.template visit<helix_intersection_initialize>(
                        sf.mask(), intersections, traj, sf, tf_store, 1e-4);
                } else {
                    mask_store.template visit<intersection_initialize>(
                        sf.mask(), intersections, traj, sf, tf_store);
                }
                // Candidate is invalid if it oversteps too far (this is neg!)
                if (intersections.empty() or
                    intersections[0].path < traj.overstep_tolerance()) {
                    continue;
                }
                // Accept if inside
                if (intersections[0].status == intersection::status::e_inside &&
                    intersections[0].direction ==
                        intersection::direction::e_along) {
                    // Volume the candidate belongs to
                    intersections[0].barcode = sf.barcode();
                    intersection_record.emplace_back(volume.index(),
                                                     intersections[0]);
                }
                intersections.clear();
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

}  // namespace detray
