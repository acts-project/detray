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

/** A functor to find intersections between trajectory and concentric cylinder
 * mask
 */
struct concentric_cylinder_intersector {

    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 1>;
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;
    using point2 = __plugin::point2<detray::scalar>;

    /** Operator function to find intersections between ray and concentric
     * cylinder mask
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

        const scalar r{mask[0]};

        // Two points on the line, thes are in the cylinder frame
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const point3 &l0 = ro;
        const point3 l1 = ro + rd;

        // swap coorinates x/y for numerical stability
        const bool swap_x_y = std::abs(rd[0]) < scalar{1e-3};

        std::size_t _x = swap_x_y ? 1 : 0;
        std::size_t _y = swap_x_y ? 0 : 1;
        const scalar k{(l0[_y] - l1[_y]) / (l0[_x] - l1[_x])};
        const scalar d{l1[_y] - k * l1[_x]};

        quadratic_equation<scalar> qe = {(1 + k * k), 2 * k * d, d * d - r * r};
        auto qe_solution = qe();

        if (std::get<0>(qe_solution) > overstep_tolerance) {
            std::array<point3, 2> candidates;
            const auto u01 = std::get<1>(qe_solution);
            std::array<scalar, 2> t01 = {0., 0.};

            candidates[0][_x] = u01[0];
            candidates[0][_y] = k * u01[0] + d;
            t01[0] = (candidates[0][_x] - ro[_x]) / rd[_x];
            candidates[0][2] = ro[2] + t01[0] * rd[2];

            candidates[1][_x] = u01[1];
            candidates[1][_y] = k * u01[1] + d;
            t01[1] = (candidates[1][_x] - ro[_x]) / rd[_x];
            candidates[1][2] = ro[2] + t01[1] * rd[2];

            // Chose the index, take the smaller positive one
            const std::size_t cindex =
                (t01[0] < t01[1] and t01[0] > overstep_tolerance)
                    ? 0
                    : (t01[0] < overstep_tolerance and
                               t01[1] > overstep_tolerance
                           ? 1
                           : 0);
            if (t01[0] > overstep_tolerance or t01[1] > overstep_tolerance) {
                intersection_type &is = ret[0];
                is.p3 = candidates[cindex];
                is.path = t01[cindex];

                // Global to local transform in cartesian coordinate
                const auto loc = trf.point_to_local(is.p3);
                is.status = mask.is_inside(loc, mask_tolerance);

                // Get intersection in local coordinate
                is.p2 = typename mask_t::local_type()(loc);

                const scalar rdir{getter::perp(is.p3 + scalar{0.1} * rd)};
                is.direction = rdir > r ? intersection::direction::e_along
                                        : intersection::direction::e_opposite;
                is.link = mask.volume_link();
            }
        }
        return ret;
    }
};

}  // namespace detray