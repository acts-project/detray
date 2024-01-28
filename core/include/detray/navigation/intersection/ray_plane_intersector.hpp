/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cartesian2.hpp"
#include "detray/coordinates/polar2.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/navigation/intersection/intersection.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename algebra_t, typename fame_t>
struct ray_intersector_impl;

/// A functor to find intersections between straight line and planar surface
template <typename algebra_t>
struct ray_intersector_impl<algebra_t, cartesian2<algebra_t>> {

    /// linear algebra types
    /// @{
    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point3 = typename algebra_t::point3;
    using point2 = typename algebra_t::point2;
    using vector3 = typename algebra_t::vector3;
    /// @}

    template <typename surface_descr_t>
    using intersection_type = intersection2D<surface_descr_t, algebra_t>;
    using ray_type = detail::ray<algebra_t>;

    /// Operator function to find intersections between ray and planar mask
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
        const transform3_type &trf, const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        // Retrieve the surface normal & translation (context resolved)
        const auto &sm = trf.matrix();
        const vector3 sn = getter::vector<3>(sm, 0u, 2u);
        const vector3 st = getter::vector<3>(sm, 0u, 3u);

        // Intersection code
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const scalar_type denom = vector::dot(rd, sn);
        // this is dangerous
        if (denom != 0.f) {
            is.path = vector::dot(sn, st - ro) / denom;

            // Intersection is valid for navigation - continue
            if (is.path >= overstep_tol) {

                const point3 p3 = ro + is.path * rd;
                is.local = mask.to_local_frame(trf, p3, ray.dir());
                is.status = mask.is_inside(is.local, mask_tolerance);

                // prepare some additional information in case the intersection
                // is valid
                if (is.status == intersection::status::e_inside) {
                    is.sf_desc = sf;

                    is.direction = detail::signbit(is.path)
                                       ? intersection::direction::e_opposite
                                       : intersection::direction::e_along;
                    is.volume_link = mask.volume_link();

                    // Get incidene angle
                    is.cos_incidence_angle = math::abs(denom);
                }
            }
        } else {
            is.status = intersection::status::e_missed;
        }

        return is;
    }

    /// Operator function to updtae an intersections between a ray and a planar
    /// surface.
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

template <typename algebra_t>
struct ray_intersector_impl<algebra_t, polar2<algebra_t>>
    : public ray_intersector_impl<algebra_t, cartesian2<algebra_t>> {};

}  // namespace detray
