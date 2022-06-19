/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/ray_plane_intersector.hpp"
#include "detray/propagator/track.hpp"

namespace detray {

/** This is an intersector struct for line surface
 * The intersectino point is obtained by calculating the point of closest
 * approach between line and track direction
 */
struct ray_line_intersector {

    using intersection_type = line_plane_intersection;
    using transform3 = __plugin::transform3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    template <
        typename transform_t, typename track_t, typename mask_t,
        std::enable_if_t<std::is_class_v<typename mask_t::local_type>, bool> =
            true,
        std::enable_if_t<not std::is_same_v<track_t, detail::ray>, bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t &trf, const track_t &track, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) const {

        return intersect(trf, detail::ray(track), mask, tolerance,
                         track.overstep_tolerance());
    }

    /** Intersection method for line surfaces
     *
     * @tparam track_t The type of the track (which carries the context
     *         object)
     * @tparam mask_t The mask type applied to the local frame
     * @tparam local_frame The local frame type to be intersected
     *
     * Contextual part:
     * @param trf the surface to be intersected
     * @param track the track information
     *
     * Non-contextual part:
     * @param mask the local mask
     * @param tolerance is the mask specific tolerance
     *
     * @return the intersection with optional parameters
     **/
    template <typename transform_t, typename mask_t,
              std::enable_if_t<std::is_class_v<typename mask_t::local_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t &trf, const detail::ray ray, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon,
        const scalar overstep_tolerance = 0.) const {

        using local_frame = typename mask_t::local_type;

        // line direction
        const auto _z = getter::vector<3>(trf.matrix(), 0, 2);

        // line center
        const auto _t = trf.translation();

        // track direction
        const auto _d = ray.dir();

        // track position
        const auto _p = ray.pos();

        const scalar zd = vector::dot(_z, _d);
        // Case for wire is parallel to track
        if (1 - std::abs(zd) < 1e-5) {
            return intersection_type{};
        }

        const scalar td = vector::dot(_t, _d);
        const scalar tz = vector::dot(_t, _z);
        const scalar pz = vector::dot(_p, _z);
        const scalar pd = vector::dot(_p, _d);

        scalar A = td + zd * pz - zd * tz - pd;
        A *= 1. / 1 - (zd * zd);
        const scalar B = pz + zd * A - tz;

        // m is the intersection point on track
        const vector3 m = _p + _d * A;

        // n is the corresponding wire position
        // const vector3 n = _t + _z * B;

        // Vector of closest approach
        // const vector3 u = m - n;

        intersection_type is;
        is.path = A;
        is.p3 = m;

        auto loc = trf.point_to_local(is.p3);
        is.p2[0] = getter::perp(loc);
        is.p2[1] = B;

        is.status = mask.template is_inside<local_frame>(loc, tolerance);

        is.direction = is.path > overstep_tolerance
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        return is;
    }
};

}  // namespace detray
