/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/soa/ray_cylinder_intersector.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <algebra::concepts::soa algebra_t, bool do_debug>
struct ray_intersector_impl<concentric_cylindrical2D<algebra_t>, algebra_t,
                            do_debug> {

    /// Linear algebra types
    /// @{
    using value_type = dvalue<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    /// @}

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, do_debug>;

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
    template <typename surface_descr_t, typename mask_t,
              typename other_algebra_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const detail::ray<other_algebra_t> &ray, const surface_descr_t &sf,
        const mask_t &mask, const transform3_type &trf,
        const darray<scalar_type, 2u> &mask_tolerance = {0.f, 1.f},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        const scalar_type r = mask[mask_t::shape::e_r];

        const auto &pos = ray.pos();
        const auto &dir = ray.dir();

        const scalar_type rd_perp_2 = dir[0] * dir[0] + dir[1] * dir[1];

        // The ray is parallel to all/any cylinder axes (z-axis)...
        if (detray::detail::all_of(rd_perp_2 <
                                   std::numeric_limits<scalar_type>::epsilon()))
            [[unlikely]] {
            is.status = decltype(is.status)(false);
            return is;
        }

        // ...otherwise, two solutions should exist, if the descriminator is
        // greater than zero
        const scalar_type ro_perp_2 = pos[0] * pos[0] + pos[1] * pos[1];
        const scalar_type rad_diff = r * r - ro_perp_2;

        const scalar_type rd_perp_inv_2 = 1.f / rd_perp_2;
        const scalar_type k =
            -rd_perp_inv_2 * (pos[0] * dir[0] + pos[1] * dir[1]);
        const scalar_type discr = rd_perp_inv_2 * rad_diff + k * k;

        // No intersection found for any cylinder
        if (detray::detail::all_of(discr < 0.f)) [[unlikely]] {
            is.status = decltype(is.status)(false);
            return is;
        }

        const scalar_type sqrt_discr = math::sqrt(discr);
        const scalar_type s1 = k + sqrt_discr;
        const scalar_type s2 = k - sqrt_discr;

        // Take the nearest solution in every lane
        auto is_smaller_sol = s1 < s2;

        scalar_type path = 0.f;
        path(is_smaller_sol) = s1;
        path(!is_smaller_sol) = s2;

        // If any of the the near solutions is outside the overstepping
        // tolerance, take the far solution (if it exists)
        const auto outside_overstep_tol = path < overstep_tol;
        if (detray::detail::any_of(outside_overstep_tol)) {
            is_smaller_sol = (s1 >= s2) && outside_overstep_tol;
            path(is_smaller_sol) = s1;
            path(!is_smaller_sol) = s2;
            // Ray is outside of all cylinders or second solutions do not
            // exist: Should not happen for portal
            if (detray::detail::all_of(path < overstep_tol) ||
                detray::detail::all_of(discr == 0.f)) {
                is.status = decltype(is.status)(false);
                return is;
            }
        }

        if constexpr (intersection_type<surface_descr_t>::is_debug()) {
            is.local = mask_t::to_local_frame(trf, pos + is.path * dir);
        }

        const auto base_tol =
            math::max(mask_tolerance[0],
                      math::min(mask_tolerance[1],
                                mask_tol_scalor * math::fabs(is.path)));
        point2_type loc;
        loc[0] = 0.f;
        loc[1] = pos[2] + path * dir[2];
        is.status = mask.is_inside(loc, base_tol);

        is.direction = !math::signbit(is.path);
        is.volume_link = mask.volume_link();

        // Mask the values where the overstepping tolerance was not met
        is.status &= (is.path >= overstep_tol);
        is.path = path;
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

}  // namespace detray
