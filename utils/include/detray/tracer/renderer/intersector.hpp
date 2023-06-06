/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
//#include "detray/intersection/intersection.hpp"
#include "detray/propagator/base_actor.hpp"

// System include(s)
#include <iostream>
#include <limits>

namespace detray {

/// Calculates the color of a pixel. Starting point of the shader pipeline
/*template<typename intersection_impl = detail::intersection_visitor>
struct intersector : detray::actor {

    // dummy output type: will be gone with future update
    using output_type = bool;

    struct state {

        DETRAY_HOST_DEVICE
        state(const point3D& ori, const vector3D& dir, scalar min, scalar max) :
m_ray{ori, 0.f, dir, 0.f}, m_interval{min, max} {}

        /// The input ray into the scene
        detail::ray<transform3D> m_ray;
        /// Interval in which ray intersections are being considered
        std::array<scalar, 2> m_interval{0.f,
std::numeric_limits<scalar>::infinity()};
        /// Resulting intersection
        std::array<line_plane_intersection, 2> m_intersections;
    };

    /// Intersect the ray with the geometry
    template <typename scene_handle_t>
    DETRAY_HOST_DEVICE
    void operator()(state &st, const scene_handle_t & geo) const {

        // Loop over volumes
        for (const auto &v : geo.volumes()) {
            // Loop over all surfaces in volume
            for (const auto &sf :
                    detray::ranges::subrange(geo.surfaces(), v)) {

                const transform3D& ctf = geo.transform_store()[sf.transform()];
                geo.mask_store().template call<intersection_impl>(
                        sf.mask(), st.m_intersections, st.m_ray, sf, ctf,
st.m_interval, 1e-4f);
            }
        }
    }
};*/

}  // namespace detray
