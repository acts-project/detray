/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/definitions/math.hpp"
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
template <typename transform3_t>
struct concentric_cylinder_intersector {

    /// Transformation matching this struct
    using scalar_type = typename transform3_t::scalar_type;
    using point3 = typename transform3_t::point3;
    using point2 = typename transform3_t::point2;
    using vector3 = typename transform3_t::vector3;
    using ray_type = detail::ray<transform3_t>;
    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 1>;

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
    template <
        typename mask_t,
        std::enable_if_t<std::is_same_v<typename mask_t::measurement_frame_type,
                                        cylindrical2<transform3_t>>,
                         bool> = true>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const ray_type &ray, const mask_t &mask, const transform3_t & /*trf*/,
        const scalar_type mask_tolerance = 0.,
        const scalar_type overstep_tolerance = 0.) const {

        output_type ret;

        const scalar_type r{mask[0]};
        // Two points on the line, thes are in the cylinder frame
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const point3 &l0 = ro;
        const point3 l1 = ro + rd;

        // swap coorinates x/y for numerical stability
        const bool swap_x_y = std::abs(rd[0]) < scalar_type{1e-3};

        std::size_t _x = swap_x_y ? 1 : 0;
        std::size_t _y = swap_x_y ? 0 : 1;
        const scalar_type k{(l0[_y] - l1[_y]) / (l0[_x] - l1[_x])};
        const scalar_type d{l1[_y] - k * l1[_x]};

        quadratic_equation<scalar_type> qe = {(1 + k * k), 2 * k * d,
                                              d * d - r * r};
        auto qe_solution = qe();

        if (std::get<0>(qe_solution) > overstep_tolerance) {
            std::array<point3, 2> candidates;
            const auto u01 = std::get<1>(qe_solution);
            std::array<scalar_type, 2> t01 = {0., 0.};

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

                is.p2 = point2{r * getter::phi(is.p3), is.p3[2]};
                // In this case, the point has to be in cylinder3 coordinates
                // for the r-check
                if constexpr (mask_t::shape::check_radius) {
                    point3 loc3D = {r, getter::phi(is.p3), is.p3[2]};
                    is.status = mask.is_inside(is.p3, mask_tolerance);
                } else {
                    is.status = mask.is_inside(is.p2, mask_tolerance);
                }
                is.direction = vector::dot(is.p3, rd) > scalar_type{0.}
                                   ? intersection::direction::e_along
                                   : intersection::direction::e_opposite;
                is.link = mask.volume_link();

                // Get incidence angle
                const scalar_type phi{is.p2[0] / mask[mask_t::shape::e_r]};
                const vector3 normal = {math_ns::cos(phi), std::sin(phi), 0};
                is.cos_incidence_angle = vector::dot(rd, normal);
            }
        }
        return ret;
    }
};

}  // namespace detray