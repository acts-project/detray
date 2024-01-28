/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/line2.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/navigation/intersection/intersection.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename algebra_t, typename fame_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for helical trajectories with line
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <typename algebra_t>
struct helix_intersector_impl<algebra_t, line2<algebra_t>> {

    using transform3_type = algebra_t;
    using scalar_type = typename transform3_type::scalar_type;
    using point3 = typename transform3_type::point3;
    using vector3 = typename transform3_type::vector3;

    template <typename surface_descr_t>
    using intersection_type = intersection2D<surface_descr_t, algebra_t>;
    using helix_type = detail::helix<transform3_type>;

    /// Operator function to find intersections between helix and line mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_desc_t is the input surface descriptor type
    ///
    /// @param h is the input helix trajectory
    /// @param sf_desc is the surface descriptor
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const helix_type &h, const surface_descr_t &sf_desc, const mask_t &mask,
        const transform3_type &trf, const scalar_type mask_tolerance = 0.f,
        const scalar_type = 0.f) const {

        intersection_type<surface_descr_t> sfi;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{1000u};

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
        while (math::abs(s - s_prev) > convergence_tolerance and
               n_tries < max_n_tries) {

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
            sfi.sf_desc = sf_desc;
            sfi.direction = math::signbit(s)
                                ? intersection::direction::e_opposite
                                : intersection::direction::e_along;
            sfi.volume_link = mask.volume_link();
        }

        return sfi;
    }

    /// Tolerance for convergence
    scalar_type convergence_tolerance{1e-3f};
};

}  // namespace detray
