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
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/helix.hpp"
#include "detray/utils/root_finding.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for helical trajectories with planar
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <algebra::concepts::aos algebra_t>
struct helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {

    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, intersection::contains_pos>;

    template <typename other_algebra_t>
    using trajectory_type = detail::helix<other_algebra_t>;

    // Maximum number of solutions this intersector can produce
    static constexpr std::uint8_t n_solutions{1u};

    using result_type = intersection_point_err<algebra_t>;

    /// Operator function to find intersections between ray and planar mask
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    ///
    /// @return the intersection
    DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
        const trajectory_type<algebra_t> &h, const transform3_type &trf,
        const scalar_type = 0.f) const {

        if (!run_rtsafe) {
            // Get the surface info
            const auto &sm = trf.matrix();
            // Surface normal
            const vector3_type sn = getter::vector<3>(sm, 0u, 2u);
            // Surface translation
            const point3_type st = getter::vector<3>(sm, 0u, 3u);

            // Starting point on the helix for the Newton iteration
            const vector3_type dist{st - h.pos(0.f)};
            scalar_type denom{vector::dot(sn, h.dir(0.f))};

            scalar_type s;
            if (denom == 0.f) {
                s = vector::norm(dist);
            }
            s = math::fabs(vector::dot(sn, dist) / denom);

            scalar_type s_prev{0.f};

            // f(s) = sn * (h.pos(s) - st) == 0
            // Run the iteration on s
            std::size_t n_tries{0u};
            while (math::fabs(s - s_prev) > convergence_tolerance &&
                   n_tries < max_n_tries) {
                // f'(s) = sn * h.dir(s)
                denom = vector::dot(sn, h.dir(s));
                // No intersection can be found if dividing by zero
                if (denom == 0.f) {
                    return {};
                }
                // x_n+1 = x_n - f(s) / f'(s)
                s_prev = s;
                s -= vector::dot(sn, h.pos(s) - st) / denom;
                ++n_tries;
            }
            // No intersection found within max number of trials
            if (n_tries == max_n_tries) {
                return {};
            }

            return {s, h.pos(s), s - s_prev};
        } else {
            // Surface normal
            const vector3_type sn = trf.z();
            // Surface translation
            const point3_type st = trf.translation();

            // Starting point on the helix for the Newton iteration
            const vector3_type dist{st - h.pos(0.f)};
            scalar_type denom{
                vector::dot(sn, h.dir(0.5f * vector::norm(dist)))};
            scalar_type s_ini;
            if (denom == 0.f) {
                s_ini = vector::norm(dist);
            } else {
                s_ini = vector::dot(sn, dist) / denom;
            }

            /// Evaluate the function and its derivative at the point @param x
            auto plane_inters_func = [&h, &st, &sn](const scalar_type x) {
                // f(s) = sn * (h.pos(s) - st) == 0
                const scalar_type f_s{vector::dot(sn, (h.pos(x) - st))};
                // f'(s) = sn * h.dir(s)
                const scalar_type df_s{vector::dot(sn, h.dir(x))};

                return std::make_tuple(f_s, df_s);
            };

            // Run the root finding algorithm
            const auto [s, ds] = newton_raphson_safe(plane_inters_func, s_ini,
                                                     convergence_tolerance,
                                                     max_n_tries, max_path);

            return {s, h.pos(s), ds};
        }
    }

    /// Tolerance for convergence
    scalar_type convergence_tolerance{1.f * unit<scalar_type>::um};
    // Guard against inifinite loops
    std::size_t max_n_tries{1000u};
    // Early exit, if the intersection is too far away
    scalar_type max_path{5.f * unit<scalar_type>::m};
    // Complement the Newton algorithm with Bisection steps
    bool run_rtsafe{true};
};

template <algebra::concepts::aos algebra_t>
struct helix_intersector_impl<polar2D<algebra_t>, algebra_t>
    : public helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {};

}  // namespace detray
