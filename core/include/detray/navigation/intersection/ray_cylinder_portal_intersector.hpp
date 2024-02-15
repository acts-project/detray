/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/concentric_cylindrical2.hpp"
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, typename algebra_t>
struct ray_intersector_impl;

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <typename algebra_t>
struct ray_intersector_impl<concentric_cylindrical2<algebra_t>, algebra_t>
    : public ray_intersector_impl<cylindrical2<algebra_t>, algebra_t> {

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

    /// Operator function to find intersections between ray and cylinder mask
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
    /// @return the closest intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type &trf, const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        // Intersecting the cylinder from the inside yield one intersection
        // along the direction of the track and one behind it
        const auto qe = this->solve_intersection(ray, mask, trf);

        // Find the closest valid intersection
        if (qe.solutions() > 0 and qe.larger() > overstep_tol) {
            // Only the closest intersection that is outside the overstepping
            // tolerance is needed
            const scalar_type t{(qe.smaller() > overstep_tol) ? qe.smaller()
                                                              : qe.larger()};
            is = this->template build_candidate<surface_descr_t>(
                ray, mask, trf, t, mask_tolerance, overstep_tol);
            is.sf_desc = sf;
        } else {
            is.status = intersection::status::e_missed;
        }

        return is;
    }

    /// Operator function to find intersections between a ray and a 2D cylinder
    ///
    /// @tparam mask_t is the input mask type
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
                               overstep_tol);
    }
};

}  // namespace detray
