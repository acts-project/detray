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
#include "detray/navigation/detail/helix.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/utils/root_finding.hpp"

// System include(s)
#include <iostream>
#include <tuple>
#include <type_traits>

namespace detray {

template <typename frame_t, typename algebra_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for helical trajectories with planar
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <typename algebra_t>
struct helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {

    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    template <typename surface_descr_t>
    using intersection_type = intersection2D<surface_descr_t, algebra_t>;
    using helix_type = detail::helix<algebra_t>;

    /// Operator function to find intersections between helix and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_desc_t is the input surface descriptor type
    ///
    /// @param h is the input helix trajectory
    /// @param sf_desc is the surface descriptor
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const helix_type &h, const surface_descr_t &sf_desc, const mask_t &mask,
        const transform3_type &trf,
        const std::array<scalar_type, 2u> mask_tolerance =
            {detail::invalid_value<scalar_type>(),
             detail::invalid_value<scalar_type>()},
        const scalar_type = 0.f, const scalar_type = 0.f) const {

        assert((mask_tolerance[0] == mask_tolerance[1]) &&
               "Helix intersectors use only one mask tolerance value");

        intersection_type<surface_descr_t> sfi;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{1000u};
        // Early exit, if the intersection is too far away
        constexpr auto max_path{50.f * unit<scalar_type>::m};

        // Get the surface info
        const auto &sm = trf.matrix();
        // Surface normal
        const vector3_type sn = getter::vector<3>(sm, 0u, 2u);
        // Surface translation
        const point3_type st = getter::vector<3>(sm, 0u, 3u);

        // Starting point on the helix for the Newton iteration
        const vector3_type dist{trf.point_to_global(mask.centroid()) -
                                h.pos(0.f)};
        scalar_type denom{vector::dot(sn, h.dir(0.5f * getter::norm(dist)))};
        scalar_type s;
        if (denom == 0.f) {
#ifdef DEBUG
            std::cout
                << "WARNING: Helix plane intersector encountered invalid value!"
                << std::endl;
#endif
            s = getter::norm(dist);
        } else {
            s = vector::dot(sn, dist) / denom;
        }

        // f(s) = sn * (h.pos(s) - st) == 0
        auto f = [&](const scalar_type x) {
            return vector::dot(sn, (h.pos(x) - st));
        };

        // f'(s) = sn * h.dir(s)
        auto df = [&](const scalar_type x) {
            return vector::dot(sn, h.dir(x));
        };

        auto evaluate_f = [&](const scalar_type x) {
            return std::make_tuple(f(x), df(x));
        };

        // Try to bracket a root

        // Initial bracket
        scalar_type a{0.9f * s};
        scalar_type b{1.1f * s};
        const std::array<scalar_type, 2> br = bracket(a, b, f);

        // Check bracket
        scalar_type f_a{f(br[0])};
        scalar_type f_b{f(br[1])};
        bool is_bracketed{std::signbit(f_a * f_b)};
        // Root is not in the detector
        bool bracket_outside_detector{
            s > max_path && ((br[0] < -max_path && br[1] < -max_path) ||
                             (br[0] > max_path && br[1] > max_path))};
        if (bracket_outside_detector) {
#ifdef DEBUG
            std::cout << "ERROR: Root outside maximum search area" << std::endl;
#endif
            return sfi;
        }

        // Run iteration only if root was not already found
        scalar_type s_prev{s};
        if (f_a != 0.f && f_b != 0.f) {

            // Update initial guess on the root after bracketing
            bool is_lower_a{f_a < 0.f};
            a = br[is_lower_a ? 0u : 1u];
            b = br[is_lower_a ? 1u : 0u];
            s = std::clamp(s, br[0], br[1]);
            // std::cout << a << ", " << b << std::endl;

            // Run the iteration on s
            s_prev = 0.f;
            std::size_t n_tries{0u};

            // f(s) = sn * (h.pos(s) - st)
            auto [f_s, df_s] = evaluate_f(s);
            while (math::fabs(s - s_prev) > convergence_tolerance &&
                   n_tries < max_n_tries) {

                // Does Newton step escape bracket?
                const scalar_type s_newton{s - f_s / df_s};
                const bool bracket_escape{
                    std::signbit((s_newton - a) * (b - s_newton))};
                // Is the convergence of Newton too slow (possibly oscillating)?
                const bool slow_convergence{
                    math::fabs(2.f * f_s) >
                    math::fabs(math::fabs(s_prev - s) * df_s)};

                // Run bisection if Newton-Raphson would be poor
                if (is_bracketed &&
                    (bracket_escape || slow_convergence || df_s == 0.f)) {
                    // Test the function sign in the middle of the interval
                    s_prev = s;
                    s = 0.5f * (a + b);
                    // std::cout << s << " Bisection" << std::endl;
                } else {
                    // No intersection can be found if dividing by zero
                    if (!is_bracketed && df_s == 0.f) {
                        return sfi;
                    }

                    // x_n+1 = x_n - f(s) / f'(s)
                    s_prev = s;
                    s = s_newton;
                    // std::cout << s << " Newton" << std::endl;
                }

                // Going out of the detector
                if (math::fabs(s) > max_path && math::fabs(s_prev) > max_path) {
#ifdef DEBUG
                    std::cout << "WARNING: Helix plane intersector: Root "
                                 "finding diverges: s = "
                              << s << std::endl;
#endif
                    return sfi;
                }

                // Update function and bracket
                std::tie(f_s, df_s) = evaluate_f(s);
                if (std::signbit(f_s)) {
                    a = s;
                } else {
                    b = s;
                }

                ++n_tries;
            }
            // No intersection found within max number of trials
            if (n_tries == max_n_tries) {
                // Should have found the root
                if (is_bracketed) {
                    std::cout << "ERROR: Helix plane intersector did not "
                                 "converge to root in ["
                              << a << ", " << b << "]" << std::endl;
                }
#ifdef DEBUG
                std::cout << "WARNING: Helix plane intersector did not "
                             "converge after "
                          << n_tries << " steps!" << std::endl;
#endif
                return sfi;
            }
        }
        // std::cout << "Converged" << std::endl;

        // Build intersection struct from helix parameters
        sfi.path = s;
        sfi.local = mask.to_local_frame(trf, h.pos(s), h.dir(s));
        sfi.cos_incidence_angle =
            vector::dot(mask.local_frame().normal(trf, sfi.local), h.dir(s));

        scalar_type tol{mask_tolerance[1]};
        if (detail::is_invalid_value(tol)) {
            // Due to floating point errors this can be negative if cos ~ 1
            const scalar_type sin_inc2{math::abs(
                1.f - sfi.cos_incidence_angle * sfi.cos_incidence_angle)};

            tol = math::abs((s - s_prev) * math::sqrt(sin_inc2));
        }
        sfi.status = mask.is_inside(sfi.local, tol);

        // Compute some additional information if the intersection is valid
        if (sfi.status) {
            sfi.sf_desc = sf_desc;
            sfi.direction = !math::signbit(s);
            sfi.volume_link = mask.volume_link();
        }

        return sfi;
    }

    /// Interface to use fixed mask tolerance
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const helix_type &h, const surface_descr_t &sf_desc, const mask_t &mask,
        const transform3_type &trf, const scalar_type mask_tolerance,
        const scalar_type = 0.f, const scalar_type = 0.f) const {
        return this->operator()(h, sf_desc, mask, trf,
                                {mask_tolerance, mask_tolerance}, 0.f);
    }

    /// Tolerance for convergence
    scalar_type convergence_tolerance{1.f * unit<scalar_type>::um};
};

template <typename algebra_t>
struct helix_intersector_impl<polar2D<algebra_t>, algebra_t>
    : public helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {};

}  // namespace detray
