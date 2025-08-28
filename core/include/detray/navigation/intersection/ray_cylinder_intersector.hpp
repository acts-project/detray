/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool do_debug>
struct ray_intersector_impl;

/// A functor to find intersections between a ray and a 2D cylinder mask
template <algebra::concepts::aos algebra_t, bool do_debug>
struct ray_intersector_impl<cylindrical2D<algebra_t>, algebra_t, do_debug> {

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
    using ray_type = detail::ray<algebra_t>;

    /// Operator function to find intersections between a ray and a 2D cylinder
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
    /// @return the intersections.
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline darray<intersection_type<surface_descr_t>, 2>
    operator()(const ray_type &ray, const surface_descr_t &sf,
               const mask_t &mask, const transform3_type &trf,
               const darray<scalar_type, 2u> mask_tolerance =
                   {0.f, 100.f * unit<scalar_type>::um},
               const scalar_type mask_tol_scalor = 0.f,
               const scalar_type overstep_tol = 0.f) const {

        // One or both of these solutions might be invalid
        const auto qe = solve_intersection(ray, mask, trf);

        darray<intersection_type<surface_descr_t>, 2> ret;

        switch (qe.solutions()) {
            case 2: {
                point3_type glob_pos = ray.pos() + qe.larger() * ray.dir();
                build_intersection(ray, ret[1], glob_pos, qe.larger(), sf, mask,
                                   trf, mask_tolerance, mask_tol_scalor,
                                   overstep_tol);
                // If there are two solutions, reuse the case for a single
                // solution to setup the intersection with the smaller path
                // in ret[0]
                [[fallthrough]];
            }
            case 1: {
                point3_type glob_pos = ray.pos() + qe.smaller() * ray.dir();
                build_intersection(ray, ret[0], glob_pos, qe.smaller(), sf,
                                   mask, trf, mask_tolerance, mask_tol_scalor,
                                   overstep_tol);
                break;
            }
            case 0: {
                ret[0].status = false;
                ret[1].status = false;
            }
        }

        // Even if there are two geometrically valid solutions, the smaller one
        // might not be passed on if it is below the overstepping tolerance:
        // see 'build_candidate'
        return ret;
    }

    /// Interface to use fixed mask tolerance
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline darray<intersection_type<surface_descr_t>, 2>
    operator()(const ray_type &ray, const surface_descr_t &sf,
               const mask_t &mask, const transform3_type &trf,
               const scalar_type mask_tolerance,
               const scalar_type overstep_tol = 0.f) const {
        return this->operator()(ray, sf, mask, trf, {mask_tolerance, 0.f}, 0.f,
                                overstep_tol);
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
        const darray<scalar_type, 2u> mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        // One or both of these solutions might be invalid
        const auto qe = solve_intersection(ray, mask, trf);

        switch (qe.solutions()) {
            case 1: {
                const point3_type glob_pos =
                    ray.pos() + qe.smaller() * ray.dir();
                build_intersection(ray, sfi, glob_pos, qe.smaller(),
                                   sfi.sf_desc, mask, trf, mask_tolerance,
                                   mask_tol_scalor, overstep_tol);
                break;
            }
            case 0: {
                sfi.status = false;
            }
        }
    }

    protected:
    /// Calculates the distance to the (two) intersection points on the
    /// cylinder in global coordinates.
    ///
    /// @returns a quadratic equation object that contains the solution(s).
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline detail::quadratic_equation<scalar_type>
    solve_intersection(const ray_type &ray, const mask_t &mask,
                       const transform3_type &trf) const {
        const scalar_type r{mask[mask_t::shape::e_r]};
        const vector3_type &sz = trf.z();
        const vector3_type &sc = trf.translation();

        const point3_type &ro = ray.pos();
        const vector3_type &rd = ray.dir();

        const auto pc_cross_sz = vector::cross(ro - sc, sz);
        const auto rd_cross_sz = vector::cross(rd, sz);
        const scalar_type a{vector::dot(rd_cross_sz, rd_cross_sz)};
        const scalar_type b{2.f * vector::dot(rd_cross_sz, pc_cross_sz)};
        const scalar_type c{vector::dot(pc_cross_sz, pc_cross_sz) - (r * r)};

        return detail::quadratic_equation<scalar_type>{a, b, c};
    }
};

}  // namespace detray
