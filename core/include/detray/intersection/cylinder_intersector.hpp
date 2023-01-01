/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2.hpp"
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
template <typename transform3_t>
struct cylinder_intersector {

    /// Transformation matching this struct
    using scalar_type = typename transform3_t::scalar_type;
    using point2 = typename transform3_t::point2;
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using ray_type = detail::ray<transform3_t>;
    using intersection_type = line_plane_intersection;

    /** Operator function to find intersections between ray and cylinder mask
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
    DETRAY_HOST_DEVICE inline std::array<intersection_type, 2> operator()(
        const ray_type &ray, const mask_t &mask, const transform3_t &trf,
        const scalar_type mask_tolerance = 0,
        const scalar_type overstep_tolerance = 0.) const {

        std::array<intersection_type, 2> ret;

        const scalar_type r{mask[0]};
        const auto &m = trf.matrix();
        const vector3 sz = getter::vector<3>(m, 0, 2);
        const vector3 sc = getter::vector<3>(m, 0, 3);

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        const vector3 pc_cross_sz = vector::cross(ro - sc, sz);
        const vector3 rd_cross_sz = vector::cross(rd, sz);
        const scalar_type a{vector::dot(rd_cross_sz, rd_cross_sz)};
        const scalar_type b{scalar_type{2.} *
                            vector::dot(rd_cross_sz, pc_cross_sz)};
        const scalar_type c{vector::dot(pc_cross_sz, pc_cross_sz) - (r * r)};

        quadratic_equation<scalar_type> qe = {a, b, c};
        auto qe_solution = qe();

        if (std::get<0>(qe_solution) > scalar_type{0.}) {
            const auto t01 = std::get<1>(qe_solution);
            const scalar_type t{(t01[0] > overstep_tolerance) ? t01[0]
                                                              : t01[1]};

            if (t > overstep_tolerance) {
                intersection_type &is = ret[0];
                is.path = t;
                is.p3 = ro + is.path * rd;
                // In this case, the point has to be in cylinder3 coordinates
                // for the r-check
                if constexpr (mask_t::shape::check_radius) {
                    auto loc3D = mask.to_local_frame(trf, is.p3);
                    is.status = mask.is_inside(loc3D, mask_tolerance);
                    is.p2 = {loc3D[0] * loc3D[1], loc3D[2]};
                } else {
                    is.p2 = mask.to_local_frame(trf, is.p3);
                    is.status = mask.is_inside(is.p2, mask_tolerance);
                }
                is.direction = vector::dot(is.p3, rd) > scalar_type{0.}
                                   ? intersection::direction::e_along
                                   : intersection::direction::e_opposite;
                is.volume_link = mask.volume_link();

                // Get incidence angle
                const scalar_type phi{is.p2[0] / mask[mask_t::shape::e_r]};
                const vector3 normal = {std::cos(phi), std::sin(phi), 0};
                is.cos_incidence_angle = vector::dot(rd, normal);
            }
        }

        return ret;
    }
};

}  // namespace detray