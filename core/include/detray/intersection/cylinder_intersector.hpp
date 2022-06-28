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
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/** A functor to find intersections between trajectory and cylinder mask
 */
struct cylinder_intersector {

    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 2>;
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /** Operator function to find intersections between ray and cylinder mask
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

        using local_frame = typename mask_t::local_type;

        const scalar r{mask[0]};
        const auto &m = trf.matrix();
        const vector3 sz = getter::vector<3>(m, 0, 2);
        const vector3 sc = getter::vector<3>(m, 0, 3);

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        const vector3 pc_cross_sz = vector::cross(ro - sc, sz);
        const vector3 rd_cross_sz = vector::cross(rd, sz);
        const scalar a = vector::dot(rd_cross_sz, rd_cross_sz);
        const scalar b = scalar{2.} * vector::dot(rd_cross_sz, pc_cross_sz);
        const scalar c = vector::dot(pc_cross_sz, pc_cross_sz) - (r * r);

        quadratic_equation<scalar> qe = {a, b, c};
        auto qe_solution = qe();

        if (std::get<0>(qe_solution) > scalar{0.}) {
            const auto t01 = std::get<1>(qe_solution);
            const scalar t{(t01[0] > overstep_tolerance) ? t01[0] : t01[1]};

            if (t > overstep_tolerance) {
                intersection_type &is = ret[0];
                is.path = t;
                is.p3 = ro + is.path * rd;
                constexpr local_frame local_converter{};
                is.p2 = local_converter(trf, is.p3);
                const auto local3 = trf.point_to_local(is.p3);
                is.status = mask.template is_inside<local_frame>(
                    local3, edge_tolerance);
                const scalar rdr = getter::perp(local3 + scalar{0.1} * rd);
                is.direction = rdr > r ? intersection::direction::e_along
                                       : intersection::direction::e_opposite;
                is.link = mask.volume_link();
            }
        }

        return ret;
    }
};

}  // namespace detray