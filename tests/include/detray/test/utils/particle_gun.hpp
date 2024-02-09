/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/test/intersection/helix_intersection_kernel.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/// @brief struct that holds functionality to shoot a parametrzed particle
/// trajectory through a detector.
///
/// Records intersections with every detector surface along the trajectory.
struct particle_gun {

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
        typename detector_t::scalar_type mask_tolerance =
            1.f * unit<typename detector_t::scalar_type>::um) {

        using intersection_t = intersection2D<typename detector_t::surface_type,
                                              typename detector_t::transform3>;

        std::vector<std::pair<dindex, intersection_t>> intersection_record;

        using helix_type =
            detail::helix<typename trajectory_t::transform3_type>;

        using intersection_kernel_t =
            std::conditional_t<std::is_same_v<trajectory_t, helix_type>,
                               helix_intersection_initialize,
                               intersection_initialize>;

        // Loop over all surfaces in the detector
        const auto &trf_store = detector.transform_store();

        std::vector<intersection_t> intersections{};

        for (const auto &sf_desc : detector.surfaces()) {
            // Retrieve candidate(s) from the surface
            const auto sf = surface{detector, sf_desc};
            sf.template visit_mask<intersection_kernel_t>(
                intersections, traj, sf_desc, trf_store,
                sf.is_portal() ? 0.f : mask_tolerance);

            // Candidate is invalid if it lies in the opposite direction
            for (auto &sfi : intersections) {
                if (sfi.direction == intersection::direction::e_along) {
                    sfi.sf_desc = sf_desc;
                    // Volume the candidate belongs to
                    intersection_record.emplace_back(sf.volume(), sfi);
                }
            }
            intersections.clear();
        }

        // Sort intersections by distance to origin of the trajectory
        auto sort_path = [&](std::pair<dindex, intersection_t> a,
                             std::pair<dindex, intersection_t> b) -> bool {
            return (a.second < b.second);
        };
        std::stable_sort(intersection_record.begin(), intersection_record.end(),
                         sort_path);

        // Make sure the intersection record terminates at world portals
        auto is_world_exit = [](const std::pair<dindex, intersection_t> &r) {
            return r.second.volume_link ==
                   detray::detail::invalid_value<decltype(
                       r.second.volume_link)>();
        };

        if (auto it = std::find_if(intersection_record.begin(),
                                   intersection_record.end(), is_world_exit);
            it != intersection_record.end()) {
            auto n{static_cast<std::size_t>(it - intersection_record.begin())};
            intersection_record.resize(n + 1u);
        }

        return intersection_record;
    }
};

}  // namespace detray
