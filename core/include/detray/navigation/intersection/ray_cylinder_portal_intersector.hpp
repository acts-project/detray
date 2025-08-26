/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool do_debug>
struct ray_intersector_impl;

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <algebra::concepts::aos algebra_t, bool do_debug>
struct ray_intersector_impl<concentric_cylindrical2D<algebra_t>, algebra_t,
                            do_debug> {

    /// linear algebra types
    /// @{
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    /// @}

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, do_debug>;
    using ray_type = detail::ray<algebra_t>;

    /// Operator function to find intersections between ray and cylinder mask
    ///
    /// Intersecting the cylinder from the inside yields one intersection
    /// along the direction of the track and one behind it. These intersections
    /// can be calculated in a simplified way, since the cylinder cannot be
    /// shifted or rotated. Solve: perp(ro + t * rd) = r_cyl
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
        const transform3_type &trf,
        const darray<scalar_type, 2u> mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;
        const scalar_type r{mask[mask_t::shape::e_r]};

        const point3_type &ro = ray.pos();
        const vector3_type &rd = ray.dir();
        scalar_type path{0.f};

        const scalar_type rd_perp_2{rd[0] * rd[0] + rd[1] * rd[1]};

        // The ray is parallel to the cylinder axis (z-axis)...
        if (rd_perp_2 < std::numeric_limits<scalar_type>::epsilon())
            [[unlikely]] {
            is.status = false;
            return is;
        }

        // ...otherwise, two solutions should exist, if the descriminator is
        // greater than zero
        const scalar_type ro_perp_2{ro[0] * ro[0] + ro[1] * ro[1]};
        const scalar_type rad_diff{r * r - ro_perp_2};

        // Only calculate the path, when not already on surface
        if (math::fabs(rad_diff) > 1.f * unit<scalar_type>::um) {

            const scalar_type rd_perp_inv_2{1.f / rd_perp_2};
            const scalar_type k{-rd_perp_inv_2 *
                                (ro[0] * rd[0] + ro[1] * rd[1])};
            const scalar_type discr{rd_perp_inv_2 * rad_diff + k * k};

            // No intersection found
            if (discr < 0.f) [[unlikely]] {
                is.status = false;
                return is;
            }

            const scalar_type sqrt_discr{math::sqrt(discr)};
            const scalar_type s1{k + sqrt_discr};
            const scalar_type s2{k - sqrt_discr};

            // Take the nearest solution
            path = (s1 < s2) ? s1 : s2;

            // If the near solution is outside the overstepping tolerance, take
            // the far solution (if it exists)
            if (path < overstep_tol) {
                path = (s1 >= s2) ? s1 : s2;
                // Ray is outside of cylinder or second solution does not exist:
                // Should not happen for portal
                if (path < overstep_tol || discr == 0.f) {
                    is.status = false;
                    return is;
                }
            }
        }

        if constexpr (intersection_type<surface_descr_t>::is_debug()) {
            is.local = mask_t::to_local_frame(trf, ro + path * rd);
        }

        // Tolerance: per mille of the distance, scaled with distance
        const auto base_tol =
            math::max(mask_tolerance[0],
                      math::min(mask_tolerance[1],
                                mask_tol_scalor * math::fabs(is.path)));
        // Portal cylinder is closed in phi, only need to check z
        is.status =
            mask.is_inside(point2_type{0.f, ro[2] + path * rd[2]}, base_tol);
        is.direction = !detail::signbit(path);
        is.volume_link = mask.volume_link();

        is.path = path;
        is.sf_desc = sf;

        return is;
    }

    /// Interface to use fixed mask tolerance
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type &trf, const scalar_type mask_tolerance,
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
        const darray<scalar_type, 2u> &mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {
        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance,
                               mask_tol_scalor, overstep_tol);
    }
};

}  // namespace detray
