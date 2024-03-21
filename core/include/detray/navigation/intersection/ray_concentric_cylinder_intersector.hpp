/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// A functor to find intersections between trajectory and concentric cylinder
/// mask
template <typename algebra_t>
struct ray_concentric_cylinder_intersector {

    /// linear algebra types
    /// @{
    using transform3_type = algebra_t;
    using scalar_type = typename transform3_type::scalar_type;
    using point3 = typename transform3_type::point3;
    using point2 = typename transform3_type::point2;
    using vector3 = typename transform3_type::vector3;
    /// @}

    template <typename surface_descr_t>
    using intersection_type = intersection2D<surface_descr_t, algebra_t>;
    using ray_type = detail::ray<transform3_type>;

    /// Operator function to find intersections between ray and concentric
    /// cylinder mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_descr_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    ///
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type & /*trf*/, const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        const scalar_type r{mask[mask_t::shape::e_r]};
        // Two points on the line, these are in the cylinder frame
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const point3 &l0 = ro;
        const point3 l1 = ro + rd;

        // swap coorinates x/y for numerical stability
        const bool swap_x_y = math::abs(rd[0]) < 1e-3f;

        unsigned int _x = swap_x_y ? 1u : 0u;
        unsigned int _y = swap_x_y ? 0u : 1u;
        const scalar_type k{(l0[_y] - l1[_y]) / (l0[_x] - l1[_x])};
        const scalar_type d{l1[_y] - k * l1[_x]};

        detail::quadratic_equation<scalar_type> qe{(1.f + k * k), 2.f * k * d,
                                                   d * d - r * r};

        if (qe.solutions() > 0) {
            const scalar_type overstep_tolerance{overstep_tol};
            std::array<point3, 2> candidates;
            std::array<scalar_type, 2> t01 = {0.f, 0.f};

            candidates[0][_x] = qe.smaller();
            candidates[0][_y] = k * qe.smaller() + d;
            t01[0] = (candidates[0][_x] - ro[_x]) / rd[_x];
            candidates[0][2] = ro[2] + t01[0] * rd[2];

            candidates[1][_x] = qe.larger();
            candidates[1][_y] = k * qe.larger() + d;
            t01[1] = (candidates[1][_x] - ro[_x]) / rd[_x];
            candidates[1][2] = ro[2] + t01[1] * rd[2];

            // Chose the index, take the smaller positive one
            const unsigned int cindex =
                (t01[0] < t01[1] and t01[0] > overstep_tolerance)
                    ? 0u
                    : (t01[0] < overstep_tolerance and
                               t01[1] > overstep_tolerance
                           ? 1u
                           : 0u);
            if (t01[0] > overstep_tolerance or t01[1] > overstep_tolerance) {

                const point3 p3 = candidates[cindex];
                const scalar_type phi{getter::phi(p3)};
                is.local = {r * phi, p3[2], r};

                is.path = t01[cindex];
                // In this case, the point has to be in cylinder3 coordinates
                // for the r-check
                is.status = mask.is_inside(is.local, mask_tolerance);

                // prepare some additional information in case the intersection
                // is valid
                if (is.status == intersection::status::e_inside) {
                    is.sf_desc = sf;
                    is.direction = detail::signbit(is.path)
                                       ? intersection::direction::e_opposite
                                       : intersection::direction::e_along;
                    is.volume_link = mask.volume_link();

                    // Get incidence angle
                    const vector3 normal = {math::cos(phi), math::sin(phi),
                                            0.f};
                    is.cos_incidence_angle = vector::dot(rd, normal);
                }
            }
        }
        return is;
    }

    /// Operator function to find intersections between a ray and a 2D cylinder
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_descr_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_type<surface_descr_t> &sfi,
        const mask_t &mask, const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tol = 0.f) const {
        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance,
                               overstep_tol)[0];
    }
};

}  // namespace detray
