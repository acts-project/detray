/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/tools/intersection_kernel.hpp"
#include "detray/utils/enumerate.hpp"

namespace detray {

/** ray object **/
struct ray {
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    point3 _pos = {0., 0., 0.};
    vector3 _dir = {0., 0., 0.};
    scalar _overstep_tolerance = -1e-4;

    DETRAY_HOST_DEVICE point3 pos() const { return _pos; }

    DETRAY_HOST_DEVICE vector3 dir() const { return _dir; }

    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }
};

/** Intersect all portals in a detector with a given ray.
 *
 * @param detector the detector
 * @param ray the ray to be shot through the detector
 *
 * @return a sorted vector of volume indices with the corresponding
 *         intersections of the portals that were encountered
 */
template <typename detector_t>
inline auto shoot_ray(const detector_t &detector, const ray &r) {

    std::vector<std::pair<dindex, intersection>> intersection_record;

    // Loop over volumes
    for (const auto &volume : detector.volumes()) {
        for (const auto [sf_idx, sf] : enumerate(detector.surfaces(), volume)) {
            // Retrieve candidate from the object
            auto sfi = intersect(r, sf, detector.transform_store(),
                                 detector.mask_store());

            // Candidate is invalid if it oversteps too far (this is neg!)
            if (sfi.path < r.overstep_tolerance()) {
                continue;
            }
            // Accept if inside
            if (sfi.status == e_inside) {
                // object the candidate belongs to
                sfi.index = volume.index();
                intersection_record.emplace_back(sf_idx, sfi);
            }
        }
    }

    // Sort intersections by distance to origin of the ray
    auto sort_path = [&](std::pair<dindex, intersection> a,
                         std::pair<dindex, intersection> b) -> bool {
        return (a.second < b.second);
    };
    std::sort(intersection_record.begin(), intersection_record.end(),
              sort_path);

    return intersection_record;
}

/** Intersect all portals in a detector with a given ray.
 *
 * @param detector the detector
 * @param origin the origin of the ray in global coordinates
 * @param direction the direction of the ray in global coordinater
 *
 * @return a sorted vector of volume indices with the corresponding
 *         intersections of the portals that were encountered
 */
template <typename detector_t>
inline auto shoot_ray(const detector_t &detector, const point3 &origin,
                      const point3 &direction) {
    ray r{origin, direction};
    return shoot_ray(detector, r);
}

}  // namespace detray
