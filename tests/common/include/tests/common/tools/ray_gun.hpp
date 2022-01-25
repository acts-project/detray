/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if defined(__CUDACC__)
#include <thrust/sort.h>
#endif

#include "detray/tools/intersection_kernel.hpp"
#include "detray/utils/enumerate.hpp"

namespace detray {

struct intersection_record {
    dindex sf_idx;
    intersection sfi;

    DETRAY_HOST_DEVICE
    bool operator<(const intersection_record &other) const {
        return (sfi < other.sfi);
    }
};

/** Intersect all portals in a detector with a given ray.
 *
 * @param detector the detector
 * @param ray the ray to be shot through the detector
 *
 * @return a sorted vector of volume indices with the corresponding
 *         intersections of the portals that were encountered
 */
template <typename detector_t, typename ray_t,
          template <typename...> class vector_t>
DETRAY_HOST_DEVICE inline auto shoot_ray(
    const detector_t &detector, const ray_t &ray,
    vector_t<intersection_record> &intersection_records) {

    // Loop over volumes
    for (const auto &volume : detector.volumes()) {
        for (const auto [sf_idx, sf] : enumerate(detector.surfaces(), volume)) {
            // Retrieve candidate from the object
            auto sfi =
                intersect(ray, sf, detector.transforms(), detector.masks());

            // Candidate is invalid if it oversteps too far (this is neg!)
            if (sfi.path < ray.overstep_tolerance) {
                continue;
            }
            // Accept if inside
            if (sfi.status == e_inside) {
                // object the candidate belongs to
                sfi.index = volume.index();
                // the next volume if we encounter the candidate
                sfi.link = std::get<0>(sf.edge());
                // NOTE: emplace_back seems not defined for device
                // intersection_record.emplace_back(sf_idx, sfi);
                intersection_records.push_back({sf_idx, sfi});
            }
        }
    }

#if defined(__CUDACC__)
    thrust::sort(thrust::seq, intersection_records.begin(),
                 intersection_records.end());
#else
    std::sort(intersection_records.begin(), intersection_records.end());
#endif
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

    using detray_context = typename detector_t::context;
    // detray_context default_context;

    track<detray_context> ray = {.pos = origin, .dir = direction};

    std::vector<intersection_record> intersection_records;

    shoot_ray(detector, ray, intersection_records);

    return intersection_records;
};

}  // namespace detray
