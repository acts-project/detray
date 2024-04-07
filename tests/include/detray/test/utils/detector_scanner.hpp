/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection_kernel.hpp"
#include "detray/navigation/intersector.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <algorithm>
#include <cmath>
#include <iostream>
#include <type_traits>

namespace detray {

/// Record of a surface intersection along a track
template <typename detector_t>
struct intersection_record {
    using intersection_type = intersection2D<typename detector_t::surface_type,
                                             typename detector_t::algebra_type>;

    /// Current global track parameters
    free_track_parameters<typename detector_t::algebra_type> track_param;
    /// Index of the volume the intersection was found in
    dindex vol_idx;
    /// The intersection result
    intersection_type intersection;
};

/// @brief struct that holds functionality to shoot a parametrzed particle
/// trajectory through a detector.
///
/// Records intersections with every detector surface along the trajectory.
template <typename trajectory_t>
struct brute_force_scan {

    template <typename D>
    using intersection_trace_type = std::vector<intersection_record<D>>;

    template <typename detector_t>
    inline auto operator()(const detector_t &detector, const trajectory_t &traj,
                           typename detector_t::scalar_type mask_tolerance =
                               1.f * unit<typename detector_t::scalar_type>::um,
                           const typename detector_t::scalar_type p =
                               1.f *
                               unit<typename detector_t::scalar_type>::GeV) {

        using sf_desc_t = typename detector_t::surface_type;
        using nav_link_t = typename detector_t::surface_type::navigation_link;

        using intersection_t =
            intersection2D<sf_desc_t, typename detector_t::algebra_type>;

        using intersection_kernel_t = intersection_initialize<intersector>;

        intersection_trace_type<detector_t> intersection_trace;

        // Loop over all surfaces in the detector
        const auto &trf_store = detector.transform_store();

        std::vector<intersection_t> intersections{};

        for (const sf_desc_t &sf_desc : detector.surfaces()) {
            // Retrieve candidate(s) from the surface
            const auto sf = surface{detector, sf_desc};
            sf.template visit_mask<intersection_kernel_t>(
                intersections, traj, sf_desc, trf_store,
                sf.is_portal() ? 0.f : mask_tolerance);

            // Candidate is invalid if it lies in the opposite direction
            for (auto &sfi : intersections) {
                if (sfi.direction == intersection::direction::e_along) {
                    sfi.sf_desc = sf_desc;
                    // Record the intersection
                    intersection_trace.push_back(
                        {{traj.pos(sfi.path), 0.f, p * traj.dir(sfi.path),
                          -1.f},
                         sf.volume(),
                         sfi});
                }
            }
            intersections.clear();
        }

        // Save initial track position as dummy intersection record
        const auto &first_record = intersection_trace.front();
        intersection_t start_intersection{};
        start_intersection.sf_desc = first_record.intersection.sf_desc;
        start_intersection.sf_desc.set_id(surface_id::e_passive);
        start_intersection.path = 0.f;
        start_intersection.volume_link =
            static_cast<nav_link_t>(first_record.vol_idx);

        intersection_trace.insert(intersection_trace.begin(),
                                  intersection_record<detector_t>{
                                      {traj.pos(), 0.f, p * traj.dir(), -1.f},
                                      first_record.vol_idx,
                                      start_intersection});

        return intersection_trace;
    }
};

template <typename algebra_t>
using ray_scan = brute_force_scan<detail::ray<algebra_t>>;

template <typename algebra_t>
using helix_scan = brute_force_scan<detail::helix<algebra_t>>;

/// Run a scan on detector object by shooting test particles through it
namespace detector_scanner {

template <template <typename> class scan_type, typename detector_t,
          typename trajectory_t, typename... Args>
inline auto run(const detector_t &detector, const trajectory_t &traj,
                Args &&...args) {

    using algebra_t = typename detector_t::algebra_type;

    auto intersection_record =
        scan_type<algebra_t>{}(detector, traj, std::forward<Args>(args)...);

    using record_t = typename decltype(intersection_record)::value_type;

    // Sort intersections by distance to origin of the trajectory
    auto sort_path = [&](const record_t &a, const record_t &b) -> bool {
        return (a.intersection < b.intersection);
    };
    std::stable_sort(intersection_record.begin(), intersection_record.end(),
                     sort_path);

    // Make sure the intersection record terminates at world portals
    auto is_world_exit = [](const record_t &r) {
        return r.intersection.volume_link ==
               detray::detail::invalid_value<
                   decltype(r.intersection.volume_link)>();
    };

    if (auto it = std::find_if(intersection_record.begin(),
                               intersection_record.end(), is_world_exit);
        it != intersection_record.end()) {
        auto n{static_cast<std::size_t>(it - intersection_record.begin())};
        intersection_record.resize(n + 1u);
    }

    return intersection_record;
}

}  // namespace detector_scanner

}  // namespace detray
