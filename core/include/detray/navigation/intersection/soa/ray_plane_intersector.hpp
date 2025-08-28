/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool do_debug>
struct ray_intersector_impl;

/// A functor to find intersections between straight line and planar surface
template <algebra::concepts::soa algebra_t, bool do_debug>
struct ray_intersector_impl<cartesian2D<algebra_t>, algebra_t, do_debug> {

    /// Linear algebra types
    /// @{
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    /// @}

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, do_debug>;

    /// Operator function to find intersections between ray and planar mask
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
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t,
              typename other_algebra_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const detail::ray<other_algebra_t> &ray, const surface_descr_t &sf,
        const mask_t &mask, const transform3_type &trf,
        const darray<scalar_type, 2u> &mask_tolerance = {0.f, 1.f},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        // Retrieve the surface normal & translation (context resolved)
        const vector3_type &sn = trf.z();
        const vector3_type &st = trf.translation();

        // Broadcast ray data
        const auto &pos = ray.pos();
        const auto &dir = ray.dir();
        const vector3_type ro{pos[0], pos[1], pos[2]};
        const vector3_type rd{dir[0], dir[1], dir[2]};

        const scalar_type denom = vector::dot(rd, sn);
        const vector3_type diff = st - ro;
        const scalar_type s = vector::dot(sn, diff) / denom;

        // Check if we divided by zero
        const auto check_sum = s.sum();
        if (!std::isnan(check_sum) && !std::isinf(check_sum)) {

            const point3_type glob_pos = ro + s * rd;

            build_intersection(ray, is, glob_pos, s, sf, mask, trf,
                               mask_tolerance, mask_tol_scalor, overstep_tol);
        } else {
            is.status = decltype(is.status)(false);
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
    template <typename surface_descr_t, typename mask_t,
              typename other_algebra_t>
    DETRAY_HOST_DEVICE inline void update(
        const detail::ray<other_algebra_t> &ray,
        intersection_type<surface_descr_t> &sfi, const mask_t &mask,
        const transform3_type &trf,
        const darray<scalar_type, 2u> &mask_tolerance = {0.f, 1.f},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {
        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance,
                               mask_tol_scalor, overstep_tol);
    }
};

template <algebra::concepts::soa algebra_t, bool do_debug>
struct ray_intersector_impl<polar2D<algebra_t>, algebra_t, do_debug>
    : public ray_intersector_impl<cartesian2D<algebra_t>, algebra_t, do_debug> {
};

}  // namespace detray
