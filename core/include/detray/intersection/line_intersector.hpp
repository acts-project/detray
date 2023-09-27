/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/line2.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// A functor to find intersections between trajectory and line mask
template <typename intersection_t>
struct line_intersector {

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

    /// Operator function to find intersections between ray and line mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    //
    /// @return the intersection
    template <typename mask_t, typename surface_t,
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              line2<transform3_type>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline intersection_t operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_type is;

        // line direction
        const vector3 _z = getter::vector<3>(trf.matrix(), 0u, 2u);

        // line center
        const point3 _t = trf.translation();

        // track direction
        const vector3 _d = ray.dir();

        // track position
        const point3 _p = ray.pos();

        // Projection of line to track direction
        const scalar_type zd{vector::dot(_z, _d)};

        const scalar_type denom{1.f - (zd * zd)};

        // Case for wire is parallel to track
        if (denom < 1e-5f) {
            is.status = intersection::status::e_missed;
            return is;
        }

        // vector from track position to line center
        const auto t2l = _t - _p;

        // t2l projection on line direction
        const scalar_type t2l_on_line{vector::dot(t2l, _z)};

        // t2l projection on track direction
        const scalar_type t2l_on_track{vector::dot(t2l, _d)};

        // path length to the point of closest approach on the track
        const scalar_type A{1.f / denom * (t2l_on_track - t2l_on_line * zd)};

        is.path = A;
        // Intersection is not valid for navigation - return early
        if (is.path >= ray.overstep_tolerance()) {

            // point of closest approach on the track
            const point3 m = _p + _d * A;

            is.local = mask.to_local_frame(trf, m, _d);

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
                is.cos_incidence_angle = std::abs(zd);
            }
        }
        return is;
    }

    /// Operator function to find intersections between helix and line mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam transform_t is the input transform type
    ///
    /// @param h is the input helix trajectory
    /// @param sf is the input surface
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
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

        // line axis direction
        const vector3 l = getter::vector<3>(trf.matrix(), 0u, 2u);

        // line center
        const point3 c = trf.translation();

        // initial track direction
        const vector3 t0 = h.dir(0.f);

        // initial track position
        const point3 r0 = h.pos(0.f);

        // Projection of line to track direction
        const scalar_type lt0{vector::dot(l, t0)};

        const scalar_type denom{1.f - (lt0 * lt0)};

        // Case for wire is parallel to track
        // @NOTE We might not have to call this which is meant to be for ray
        // intersection...
        if (denom < 1e-5f) {
            sfi.status = intersection::status::e_missed;
            return sfi;
        }

        // vector from track position to line center
        const vector3 D = c - r0;

        // D projection on line direction
        const scalar_type P{vector::dot(D, l)};

        // D projection on track direction
        const scalar_type Q{vector::dot(D, t0)};

        // Path length to the point of closest approach on the track
        // @NOTE Ray intersection algorithm is used for the initial guess on the
        // path length
        scalar_type s{1.f / denom * (Q - P * lt0)};
        scalar_type s_prev{0.f};

        // Run the iteration on s
        std::size_t n_tries{0u};
        while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {

            // track direction
            const vector3 t = h.dir(s);

            // track position
            const point3 r = h.pos(s);

            // Projection of (track position - center) to the line
            const scalar_type A = vector::dot(r - c, l);

            // Vector orthogonal to the line and passing the track position
            // w = r - (c + ((r - c) * l)l)
            const vector3 w = r - (c + A * l);

            // f(s) = t * w = 0
            const scalar_type f = vector::dot(t, w);

            // dtds = d^2r/ds^2 = qop * (t X b_field)
            const vector3 dtds = h.qop() * vector::cross(t, *h._mag_field);
            // dwds = t - (t * l)l
            const vector3 dwds = t - vector::dot(t, l) * l;

            // f'(s) = dtds * w + t * dwds
            const scalar_type dfds =
                vector::dot(dtds, w) + vector::dot(t, dwds);

            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= f / dfds;

            ++n_tries;
        }

        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return sfi;
        }

        // Build intersection struct from helix parameters
        sfi.path = s;
        sfi.local = mask.to_local_frame(trf, h.pos(s), h.dir(s));

        const point3 local = mask.to_local_frame(trf, h.pos(s), h.dir(s));
        sfi.status = mask.is_inside(local, mask_tolerance);

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

    /// Operator function to find intersections between a ray and a line.
    ///
    /// @tparam mask_t is the input mask type
    ///
    /// @param ray is the input ray trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    template <typename traj_t, typename mask_t,
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              line2<transform3_type>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline void update(
        const traj_t &traj, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(traj, sfi.sf_desc, mask, trf, mask_tolerance);
    }
};

}  // namespace detray
