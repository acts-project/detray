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

/** A functor to find intersections between trajectory and line mask
 */
struct line_intersector {

    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 1>;
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /** Operator function to find intersections between ray and line mask
     *
     * @tparam mask_t is the input mask type
     * @tparam transform_t is the input transform type
     *
     * @param ray is the input ray trajectory
     * @param mask is the input mask
     * @param trf is the transform
     * @param edge_tolerance is the tolerance for mask edges
     * @param overstep_tolerance is the tolerance for track overstepping
     *
     * @return the intersection
     */
    template <typename mask_t, typename transform_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const detail::ray &ray, const mask_t &mask, const transform_t &trf,
        const scalar edge_tolerance = 0,
        const scalar overstep_tolerance = 0.) const {

        output_type ret;

        // line direction
        const vector3 _z = getter::vector<3>(trf.matrix(), 0, 2);

        // line center
        const point3 _t = trf.translation();

        // track direction
        const vector3 _d = ray.dir();

        // track position
        const point3 _p = ray.pos();

        // Projection of line to track direction
        const scalar zd = vector::dot(_z, _d);

        // Case for wire is parallel to track
        if (1 - static_cast<detray::scalar>(std::abs(zd)) < 1e-5) {
            return ret;
        }

        const scalar denom = 1 - (zd * zd);

        // vector from track position to line center
        const vector3 t2l = _t - _p;

        // t2l projection on line direction
        const scalar t2l_on_line = vector::dot(t2l, _z);

        // t2l projection on track direction
        const scalar t2l_on_track = vector::dot(t2l, _d);

        // path length to the point of closest approach on the track
        const scalar A = 1. / denom * (t2l_on_track - t2l_on_line * zd);

        // distance to the point of closest approarch on the
        // line from line center
        const scalar B = zd * A - t2l_on_line;

        // point of closest approach on the track
        const vector3 m = _p + _d * A;

        intersection_type &is = ret[0];
        is.path = A;
        is.p3 = m;

        // global to local transform in cartesian coordinate
        auto loc = trf.point_to_local(is.p3);
        is.p2[0] = getter::perp(loc);
        is.p2[1] = B;

        is.status = mask.is_inside(loc, edge_tolerance);

        is.direction = is.path > overstep_tolerance
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        return ret;
    }
};

}  // namespace detray