/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
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

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <typename transform3_t>
struct cylinder_portal_intersector {

    /// linear algebra types
    /// @{
    using scalar_type = typename transform3_t::scalar_type;
    using point3 = typename transform3_t::point3;
    using point2 = typename transform3_t::point2;
    using vector3 = typename transform3_t::vector3;
    /// @}
    using ray_type = detail::ray<transform3_t>;
    using intersection_type = line_plane_intersection;

    /// Operator function to find intersections between ray and cylinder mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam transform_t is the input transform type
    ///
    /// @param ray is the input ray trajectory
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the closest intersection
    template <
        typename mask_t,
        std::enable_if_t<std::is_same_v<typename mask_t::measurement_frame_type,
                                        cylindrical2<transform3_t>>,
                         bool> = true>
    DETRAY_HOST_DEVICE inline std::array<intersection_type, 1> operator()(
        const ray_type &ray, const mask_t &mask, const transform3_t &trf,
        const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tolerance = 0.f) const {

        std::array<intersection_type, 1> ret;

        const scalar_type r{mask[mask_t::shape::e_r]};
        const auto &m = trf.matrix();
        const vector3 sz = getter::vector<3>(m, 0, 2);
        const vector3 sc = getter::vector<3>(m, 0, 3);

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        const auto pc_cross_sz = vector::cross(ro - sc, sz);
        const auto rd_cross_sz = vector::cross(rd, sz);
        const scalar_type a{vector::dot(rd_cross_sz, rd_cross_sz)};
        const scalar_type b{2.f * vector::dot(rd_cross_sz, pc_cross_sz)};
        const scalar_type c{vector::dot(pc_cross_sz, pc_cross_sz) - (r * r)};

        detail::quadratic_equation<scalar_type> qe{a, b, c};

        if (qe.solutions > 0) {
            // Find the first valid intersection
            const scalar_type t{(qe.smaller() > overstep_tolerance)
                                    ? qe.smaller()
                                    : qe.larger()};

            if (t > overstep_tolerance) {
                intersection_type &is = ret[0];
                is.path = t;
                is.p3 = ro + is.path * rd;
                // In this case, the point has to be in cylinder3 coordinates
                // for the r-check
                if constexpr (mask_t::shape::check_radius) {
                    const auto loc3D = mask.to_local_frame(trf, is.p3);
                    is.status = mask.is_inside(loc3D, mask_tolerance);
                    is.p2 = {loc3D[0] * loc3D[1], loc3D[2]};
                } else {
                    is.p2 = mask.to_local_frame(trf, is.p3);
                    is.status = mask.is_inside(is.p2, mask_tolerance);
                }
                // prepare some additional information in case the intersection
                // is valid
                if (is.status == intersection::status::e_inside) {
                    is.direction = vector::dot(is.p3, rd) > 0.f
                                       ? intersection::direction::e_along
                                       : intersection::direction::e_opposite;
                    is.volume_link = mask.volume_link();

                    // Get incidence angle
                    const scalar_type phi{is.p2[0] / r};
                    const vector3 normal = {std::cos(phi), std::sin(phi), 0.f};
                    is.cos_incidence_angle = vector::dot(rd, normal);
                }
            }
        }

        return ret;
    }
};

}  // namespace detray