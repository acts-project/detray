/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2D.hpp"
#include "detray/definitions/boolean.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/utils/soa/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray::soa {

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <typename intersection_t>
struct cylinder_portal_intersector
    : public cylinder_intersector<intersection_t> {

    /// Linear algebra types
    /// @{
    using algebra = typename intersection_t::algebra;
    using transform3_type = typename intersection_t::transform3D;
    using value_type = typename intersection_t::value_t;
    using scalar_type = typename intersection_t::scalar_t;
    using point3 = typename intersection_t::point3D;
    using point2 = typename intersection_t::point2D;
    using vector3 = typename intersection_t::vector3D;
    /// @}

    using intersection_type = intersection_t;
    using ray_type =
        detray::detail::ray<dtransform3D<ALGEBRA_PLUGIN<value_type>>>;

    /// Operator function to find intersections between ray and cylinder mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return the closest intersection
    template <typename mask_t, typename surface_t,
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              cylindrical2D<algebra>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline intersection_t operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_t is;

        // Intersecting the cylinder from the inside yield one intersection
        // along the direction of the track and one behind it
        const auto qe = this->solve_intersection(ray, mask, trf);

        // None of the cylinders has a valid intersection
        if (detray::detail::all_of(qe.solutions() <= 0) or
            detray::detail::all_of(qe.larger() <= ray.overstep_tolerance())) {
            is.status = decltype(is.status)(false);
            return is;
        }

        // Only the closest intersection that is outside the overstepping
        // tolerance is needed
        const auto valid_smaller = (qe.smaller() > ray.overstep_tolerance());
        scalar_type t = 0.f;
        t(valid_smaller) = qe.smaller();
        t(!valid_smaller) = qe.larger();

        is = this->template build_candidate(ray, mask, trf, t, mask_tolerance);
        is.sf_desc = sf;

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
    template <typename mask_t,
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              cylindrical2D<algebra>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance);
    }
};

}  // namespace detray::soa
