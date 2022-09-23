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
template <typename transform3_t>
struct line_intersector {

    /// Transformation matching this struct
    using scalar_type = typename transform3_t::scalar_type;
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using ray_type = detail::ray<transform3_t>;
    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 1>;

    /** Operator function to find intersections between ray and line mask
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
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const ray_type &ray, const mask_t &mask, const transform3_t &trf,
        const scalar_type mask_tolerance = 0,
        const scalar_type overstep_tolerance = 0.) const {

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

        const scalar_type denom = scalar_type{1.} - (zd * zd);

        // Case for wire is parallel to track
        if (denom < scalar_type{1e-5}) {
            return ret;
        }

        // vector from track position to line center
        const vector3 t2l = _t - _p;

        // t2l projection on line direction
        const scalar_type t2l_on_line = vector::dot(t2l, _z);

        // t2l projection on track direction
        const scalar_type t2l_on_track = vector::dot(t2l, _d);

        // path length to the point of closest approach on the track
        const scalar_type A =
            scalar_type{1.} / denom * (t2l_on_track - t2l_on_line * zd);

        // distance to the point of closest approarch on the
        // line from line center
        const scalar_type B = zd * A - t2l_on_line;

        // point of closest approach on the track
        const vector3 m = _p + _d * A;

        intersection_type &is = ret[0];
        is.path = A;
        is.p3 = m;

        // global to local transform in cartesian coordinate
        auto loc = trf.point_to_local(is.p3);

        // assign the sign depending on the position w.r.t line surface
        // Right: -1
        // Left: 1
        const auto r = vector::cross(_z, _d);
        const scalar_type sign = vector::dot(r, t2l) > 0. ? -1. : 1.;

        is.p2[0] = sign * getter::perp(loc);
        is.p2[1] = B;

        is.status = mask.is_inside(loc, mask_tolerance);

        is.direction = is.path > overstep_tolerance
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        // Get incidence angle
        is.cos_incidence_angle = std::abs(vector::dot(_d, _z));

        return ret;
    }
};

}  // namespace detray