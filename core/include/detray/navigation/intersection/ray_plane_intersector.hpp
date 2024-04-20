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
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, typename algebra_t>
struct ray_intersector_impl;

/// A functor to find intersections between straight line and planar surface
template <typename algebra_t>
struct ray_intersector_impl<cartesian2D<algebra_t>, algebra_t> {

    /// linear algebra types
    /// @{
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
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
        const transform3_type &trf,
        const std::array<scalar_type, 2u> mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type overstep_tol = 0.f) const {

        assert((mask_tolerance[0] <= mask_tolerance[1]) &&
               "Minimal mask tolerance needs to be smaller or equal maximal "
               "mask tolerance");

        intersection_type<surface_descr_t> is;

        // Retrieve the surface normal & translation (context resolved)
        const auto &sm = trf.matrix();
        const vector3_type sn = getter::vector<3>(sm, 0u, 2u);
        const vector3_type st = getter::vector<3>(sm, 0u, 3u);

        // Intersection code
        const point3_type &ro = ray.pos();
        const vector3_type &rd = ray.dir();
        const scalar_type denom = vector::dot(rd, sn);
        // this is dangerous
        if (denom != 0.f) {
            is.path = vector::dot(sn, st - ro) / denom;

            // Intersection is valid for navigation - continue
            if (is.path >= overstep_tol) {

                const point3_type p3 = ro + is.path * rd;
                is.local = mask.to_local_frame(trf, p3, ray.dir());
                // Tolerance: per mille of the distance
                is.status = mask.is_inside(
                    is.local, math::max(mask_tolerance[0],
                                        math::min(mask_tolerance[1],
                                                  1e-3f * math::abs(is.path))));

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

    /// Interface to use fixed mask tolerance
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type &trf, const scalar_type mask_tolerance,
        const scalar_type overstep_tol = 0.f) const {
        return this->operator()(ray, sf, mask, trf, {mask_tolerance, 0.f},
                                overstep_tol);
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
        const std::array<scalar_type, 2u> &mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type overstep_tol = 0.f) const {
        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance,
                               overstep_tol);
    }
};

template <typename algebra_t>
struct ray_intersector_impl<polar2D<algebra_t>, algebra_t>
    : public ray_intersector_impl<cartesian2D<algebra_t>, algebra_t> {};

}  // namespace detray
