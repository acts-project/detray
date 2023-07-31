/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/line_intersector.hpp"
#include "tests/common/tools/intersectors/helix_intersector.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

namespace detail {

/// @brief Intersection implementation for helical trajectories with line
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <typename intersection_t>
struct helix_line_intersector {

    using transform3_type = typename intersection_t::transform3_type;
    using scalar_type = typename transform3_type::scalar_type;
    using matrix_operator = typename transform3_type::matrix_actor;
    using point3 = typename transform3_type::point3;
    using vector3 = typename transform3_type::vector3;
    using helix_type = detail::helix<transform3_type>;

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
};

}  // namespace detail

/// Specialization of the @c helix_intersector for line surfaces
template <typename intersection_t, typename mask_t>
struct helix_intersector<
    intersection_t, mask_t,
    std::enable_if_t<
        std::is_same_v<
            typename mask_t::shape::template intersector_type<intersection_t>,
            line_intersector<intersection_t>>,
        void>> : public detail::helix_line_intersector<intersection_t> {

    using intersector_impl = detail::helix_line_intersector<intersection_t>;

    using scalar_type = typename intersector_impl::scalar_type;
    using matrix_operator = typename intersector_impl::matrix_operator;
    using point3 = typename intersector_impl::point3;
    using vector3 = typename intersector_impl::vector3;
    using helix_type = typename intersector_impl::helix_type;
};

}  // namespace detray
