/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"

// System include(s).
#include <type_traits>

namespace detray {

/// A functor to find intersections between trajectory and planar surface
template <typename intersection_t>
struct plane_intersector {

    /// linear algebra types
    /// @{
    using transform3_type = typename intersection_t::transform3_type;
    using scalar_type = typename transform3_type::scalar_type;
    using point3 = typename transform3_type::point3;
    using point2 = typename transform3_type::point2;
    using vector3 = typename transform3_type::vector3;
    /// @}

    using intersection_type = intersection_t;
    using ray_type = detail::ray<transform3_type>;
    using helix_type = detail::helix<transform3_type>;

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
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline intersection_t operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_t is;

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
            if (is.path >= ray.overstep_tolerance()) {

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
                    is.cos_incidence_angle = std::abs(denom);
                }
            }
        } else {
            is.status = intersection::status::e_missed;
        }

        return is;
    }

    /// Operator function to find intersections between helix and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam transform_t is the input transform type
    ///
    /// @param h is the input helix trajectory
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the intersection
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline intersection_t operator()(
        const helix_type &h, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_t sfi;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{1000u};
        // Tolerance for convergence
        constexpr scalar_type tol{1e-3f};

        // Get the surface info
        const auto &sm = trf.matrix();
        // Surface normal
        const vector3 sn = getter::vector<3>(sm, 0u, 2u);
        // Surface translation
        const point3 st = getter::vector<3>(sm, 0u, 3u);

        // Starting point on the helix for the Newton iteration
        scalar_type s{getter::norm(st - h.pos(0.f))};
        scalar_type s_prev{0.f};

        // f(s) = sn * (h.pos(s) - st) == 0
        // Run the iteration on s
        std::size_t n_tries{0u};
        while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {
            // f'(s) = sn * h.dir(s)
            const scalar_type denom{vector::dot(sn, h.dir(s))};
            // No intersection can be found if dividing by zero
            if (denom == 0.f) {
                return sfi;
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= vector::dot(sn, h.pos(s) - st) / denom;
            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return sfi;
        }

        // Build intersection struct from helix parameters
        sfi.path = s;
        sfi.local = mask.to_local_frame(trf, h.pos(s), h.dir(s));
        sfi.status = mask.is_inside(sfi.local, mask_tolerance);

        // Compute some additional information if the intersection is valid
        if (sfi.status == intersection::status::e_inside) {
            sfi.sf_desc = sf;
            sfi.direction = std::signbit(s)
                                ? intersection::direction::e_opposite
                                : intersection::direction::e_along;
            sfi.volume_link = mask.volume_link();
        }

        return sfi;
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
    template <typename traj_t, typename mask_t>
    DETRAY_HOST_DEVICE inline void update(
        const traj_t &traj, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(traj, sfi.sf_desc, mask, trf, mask_tolerance);
    }
};

}  // namespace detray
