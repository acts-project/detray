/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/quadratic_equation.hpp"
#include "detray/propagator/track.hpp"

namespace detray {

struct unbound;

/** This is an intersector struct for generic cylinder surface
 */
struct cylinder_intersector {

    using intersection_t = line_plane_intersection;

    using transform3 = __plugin::transform3<detray::scalar>;
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;
    using point2 = __plugin::point2<detray::scalar>;
    using cylindrical2 = __plugin::cylindrical2<detray::scalar>;

    /** Intersection method for cylindrical surfaces
     *
     * @tparam track_t The type of the track carrying the context object
     * @tparam mask_t The mask type applied to the local frame
     *
     * Contextual part:
     * @param trf the transform of the surface to be intersected
     * @param track the track information
     *
     * Non-contextual part:
     * @param mask the local mask
     * @param tolerance is the mask specific tolerance
     *
     * @return the intersection with optional parameters
     **/
    template <
        typename track_t, typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detray::unbound>,
            bool> = true>
    DETRAY_HOST_DEVICE inline intersection_t intersect(
        const transform3 &trf, const track_t &track, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) const {
        return intersect(trf, track.pos(), track.dir(), mask, tolerance,
                         track.overstep_tolerance());
    }

    /** Intersection method for cylindrical surfaces
     *
     * @tparam mask_t The mask type applied to the local frame
     *
     * Contextual part:
     * @param trf the transform of the surface to be intersected
     * @param ro the origin of the ray
     * @param rd the direction of the ray
     * @param local to the local local frame
     *
     * Non-contextual part:
     * @param mask the local mask
     * @param tolerance is the mask specific tolerance
     * @param overstep_tolerance  is the stepping specific tolerance
     *
     * @return the intersection with optional parameters
     **/
    template <
        typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detray::unbound>,
            bool> = true>
    DETRAY_HOST_DEVICE inline intersection_t intersect(
        const transform3 &trf, const point3 &ro, const vector3 &rd,
        const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon,
        scalar overstep_tolerance = 0.) const {

        using local_frame = typename mask_t::local_type;

        scalar r = mask[0];
        const auto &m = trf.matrix();
        auto sz = getter::vector<3>(m, 0, 2);
        auto sc = getter::vector<3>(m, 0, 3);

        vector3 pc_cross_sz = vector::cross(ro - sc, sz);
        vector3 rd_cross_sz = vector::cross(rd, sz);
        scalar a = vector::dot(rd_cross_sz, rd_cross_sz);
        scalar b = 2. * vector::dot(rd_cross_sz, pc_cross_sz);
        scalar c = vector::dot(pc_cross_sz, pc_cross_sz) - (r * r);

        quadratic_equation<scalar> qe = {a, b, c};
        auto qe_solution = qe();

        if (std::get<0>(qe_solution) > 0) {
            auto t01 = std::get<1>(qe_solution);
            scalar t = (t01[0] > overstep_tolerance) ? t01[0] : t01[1];
            if (t > overstep_tolerance) {
                intersection_t is;
                is.path = t;
                is.p3 = ro + is.path * rd;
                constexpr local_frame local_converter{};
                is.p2 = local_converter(trf, is.p3);
                auto local3 = trf.point_to_local(is.p3);
                is.status =
                    mask.template is_inside<local_frame>(local3, tolerance);
                scalar rdr = getter::perp(local3 + 0.1 * rd);
                is.direction = rdr > r ? intersection::direction::e_along
                                       : intersection::direction::e_opposite;
                is.link = mask.volume_link();
                return is;
            }
        }
        return intersection_t{};
    }
};

}  // namespace detray
