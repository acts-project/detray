/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/** A functor to find intersections between trajectory and planar mask
 */
struct plane_intersector {

    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 1>;
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /** Operator function to find intersections between ray and planar mask
     *
     * @tparam mask_t is the input mask type
     * @tparam transform_t is the input transform type
     *
     * @param ray is the input ray trajectory
     * @param mask is the input mask
     * @param trf is the transform
     * @param mask_tolerance is the tolerance for mask edges
     * @param overstep_tolerance is the tolerance for track overstepping
     *
     * @return the intersection
     */
    template <typename mask_t, typename transform_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const detail::ray &ray, const mask_t &mask, const transform_t &trf,
        const scalar mask_tolerance = 0,
        const scalar overstep_tolerance = 0.) const {

        output_type ret;

        using local_frame = typename mask_t::local_type;

        // Retrieve the surface normal & translation (context resolved)
        const auto &sm = trf.matrix();
        const vector3 sn = getter::vector<3>(sm, 0, 2);
        const vector3 st = getter::vector<3>(sm, 0, 3);

        // Intersection code
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const scalar denom = vector::dot(rd, sn);
        if (denom != scalar{0.}) {
            intersection_type &is = ret[0];
            is.path = vector::dot(sn, st - ro) / denom;
            is.p3 = ro + is.path * rd;
            constexpr local_frame local_converter{};
            is.p2 = local_converter(trf, is.p3);
            is.status =
                mask.template is_inside<local_frame>(is.p2, mask_tolerance);
            is.direction = is.path > overstep_tolerance
                               ? intersection::direction::e_along
                               : intersection::direction::e_opposite;
            is.link = mask.volume_link();

            // Get incidene angle
            is.cos_incidence_angle = std::abs(denom);
        }
        return ret;
    }
};

}  // namespace detray