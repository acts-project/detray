/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s).
#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

namespace detray {

/// @brief Try to find a bracket around a root
///
/// @param a lower initial boundary
/// @param b upper initial boundary
/// @param f function for which to find the root
/// @param k scale factor with which to widen the bracket at every step
///
/// @see Numerical Recepies pp. 445
///
/// @return bracket around root
template <typename scalar_t, typename function_t>
DETRAY_HOST_DEVICE inline std::array<scalar_t, 2> expand_bracket(
    const scalar_t a, const scalar_t b, function_t &f,
    const scalar_t k = 0.1f) {

    if (a == b) {
        throw std::invalid_argument(
            "Root bracketing: Not a valid start interval [" +
            std::to_string(a) + ", " + std::to_string(b) + "]");
    }

    scalar_t lower{a > b ? b : a};
    scalar_t upper{a > b ? a : b};

    // Sample function points at interval
    scalar_t f_l{f(lower)};
    scalar_t f_u{f(upper)};
    std::size_t n_tries{0u};

    // If there is no sign change in interval, we don't know if there is a root
    while (!math::signbit(f_l * f_u)) {
        // No interval could be found to bracket the root
        // Might be correct, if there is not root
        if ((n_tries == 1000u) || !std::isfinite(f_l) || !std::isfinite(f_u)) {
#ifdef DEBUG
            std::cout << "WARNING: Could not bracket a root" << std::endl;
#endif
            return {a, b};
        }
        scalar_t d{k * (upper - lower)};
        // Make interval larger in the direction where the function is smaller
        if (math::fabs(f_l) < math::fabs(f_u)) {
            lower -= d;
            f_l = f(lower);
        } else {
            upper += d;
            f_u = f(upper);
        }
        ++n_tries;
    }

    return {lower, upper};
}

/// @brief Find a root using the Newton-Raphson and Bisection algorithms
///
/// @param s initial guess for the root
/// @param evaluate_func evaluate the function and its derivative
/// @param max_path don't consider root if it is too far away
///
/// @see Numerical Recepies pp. 445
///
/// @return pathlength to root and the last step size
template <typename scalar_t, typename function_t>
DETRAY_HOST_DEVICE inline std::pair<scalar_t, scalar_t> newton_raphson_safe(
    function_t &evaluate_func, scalar_t s,
    const scalar_t convergence_tolerance = 1.f * unit<scalar_t>::um,
    const std::size_t max_n_tries = 1000u,
    const scalar_t max_path = 5.f * unit<scalar_t>::m) {

    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    constexpr scalar_t epsilon{std::numeric_limits<scalar_t>::epsilon()};

    // First, try to bracket a root
    auto f = [&evaluate_func](const scalar_t x) {
        auto [f_x, df_x] = evaluate_func(x);

        return f_x;
    };

    // Initial bracket
    scalar_t a{s == 0.f ? -0.01f : 0.99f * s};
    scalar_t b{s == 0.f ? 0.01f : 1.09f * s};
    const std::array<scalar_t, 2> br = expand_bracket(a, b, f);

    // Check bracket
    [[maybe_unused]] auto [f_a, df_a] = evaluate_func(br[0]);
    [[maybe_unused]] auto [f_b, df_b] = evaluate_func(br[1]);
    bool is_bracketed{math::signbit(f_a * f_b)};

    // Root is not within the maximal pathlength
    bool bracket_outside_tol{s > max_path &&
                             ((br[0] < -max_path && br[1] < -max_path) ||
                              (br[0] > max_path && br[1] > max_path))};
    if (bracket_outside_tol) {
#ifdef DEBUG
        std::cout << "INFO: Root outside maximum search area - skipping"
                  << std::endl;
#endif
        return std::make_pair(inv, inv);
    }

    // Root already found?
    if (math::fabs(f_a) < convergence_tolerance) {
        return std::make_pair(a, epsilon);
    }
    if (math::fabs(f_b) < convergence_tolerance) {
        return std::make_pair(b, epsilon);
    }

    // Update initial guess on the root after bracketing
    // Did the original guess already contain the root?
    s = math::fabs(b - a) < math::fabs(br[1] - br[0]) ? 0.5f * (br[1] - br[0])
                                                      : s;
    // Make 'a' the boundary for the negative function value -> easier to update
    bool is_lower_a{f_a < 0.f};
    a = br[is_lower_a ? 0u : 1u];
    b = br[is_lower_a ? 1u : 0u];

    // Run the iteration on s
    scalar_t s_prev{0.f};
    std::size_t n_tries{0u};
    auto [f_s, df_s] = evaluate_func(s);

    while (math::fabs(s - s_prev) > convergence_tolerance) {

        // Allow more Newton steps for faster convergence

        // Does Newton step escape bracket?
        /*bool bracket_escape{true};
        scalar_t s_newton{0.f};
        if (math::fabs(df_s) != 0.f) {
            s_newton = s - f_s / df_s;
            bracket_escape = math::signbit((s_newton - a) * (b - s_newton));
        }*/

        // This criterion from Numerical Recipes seems to work, but why?
        /*const bool slow_convergence{math::fabs(2.f * f_s) >
                                    math::fabs((s_prev - s) * df_s)};*/

        // Take a bisection step if it converges faster than Newton
        // |f(next_newton_s)| > |f(next_bisection_s)|
        /*bool slow_convergence{true};
        // If not converged far enough, take bisection step
        if (math::fabs(s - s_prev) < 0.1f * unit<scalar_t>::mm) {
            const scalar_t ds_bisection{0.5f * (a + b) - s};
            slow_convergence = (2.f * math::fabs(f_s) > math::fabs(df_s *
        ds_bisection + f_s));
        }*/

        s_prev = s;

        // Run bisection if Newton-Raphson would be poor
        if (is_bracketed/* &&
            (bracket_escape || slow_convergence || df_s == 0.f)*/) {
            // Test the function sign in the middle of the interval
            s = 0.5f * (a + b);
        } else {
            // No intersection can be found if dividing by zero
            if (!is_bracketed && df_s == 0.f) {
                std::cout << "WARNING: Encountered invalid derivative "
                          << std::endl;

                return std::make_pair(inv, inv);
            }

            s = s_newton;
        }

        // Update function and bracket
        std::tie(f_s, df_s) = evaluate_func(s);
        if (is_bracketed && math::signbit(f_s)) {
            a = s;
        } else {
            b = s;
        }

        // Converges to a point outside the search space
        if (math::fabs(s) > max_path && math::fabs(s_prev) > max_path &&
            ((a < -max_path && b < -max_path) ||
             (a > max_path && b > max_path))) {
#ifdef DEBUG
            std::cout << "WARNING: Root finding left the search space"
                      << std::endl;
#endif
            return std::make_pair(inv, inv);
        }

        ++n_tries;
        // No intersection found within max number of trials
        if (n_tries >= max_n_tries) {

            // Should have found the root
            if (is_bracketed) {
                throw std::runtime_error(
                    "ERROR: Helix intersector did not "
                    "find root for s=" +
                    std::to_string(s) + " in [" + std::to_string(a) + ", " +
                    std::to_string(b) + "]");
            } else {
#ifdef DEBUG
                std::cout << "WARNING: Helix intersector did not "
                             "converge after "
                          << n_tries << " steps unbracketed search"
                          << std::endl;
#endif
            }
            return std::make_pair(inv, inv);
        }
    }

    // Final pathlengt to root
    return std::make_pair(s, math::fabs(s - s_prev));
}

/// @brief Fill an intersection with the result of the root finding
///
/// @param [in] traj the test trajectory that intersects the surface
/// @param [out] sfi the surface intersection
/// @param [in] s path length to the root
/// @param [in] ds approximation error for the root
/// @param [in] mask the mask of the surface
/// @param [in] trf the transform of the surface
/// @param [in] mask_tolerance minimal and maximal mask tolerance
template <typename scalar_t, typename intersection_t, typename surface_descr_t,
          typename mask_t, typename trajectory_t, typename transform3_t>
DETRAY_HOST_DEVICE inline void build_intersection(
    const trajectory_t &traj, intersection_t &sfi, const scalar_t s,
    const scalar_t ds, const surface_descr_t sf_desc, const mask_t &mask,
    const transform3_t &trf, const std::array<scalar_t, 2> &mask_tolerance) {

    // Build intersection struct from test trajectory, if the distance is valid
    if (!detail::is_invalid_value(s)) {
        sfi.path = s;
        sfi.local = mask.to_local_frame(trf, traj.pos(s), traj.dir(s));
        sfi.cos_incidence_angle =
            vector::dot(mask.local_frame().normal(trf, sfi.local), traj.dir(s));

        scalar_t tol{mask_tolerance[1]};
        if (detail::is_invalid_value(tol)) {
            // Due to floating point errors this can be negative if cos ~ 1
            const scalar_t sin_inc2{math::fabs(
                1.f - sfi.cos_incidence_angle * sfi.cos_incidence_angle)};

            tol = math::fabs(ds * math::sqrt(sin_inc2));
        }
        sfi.status = mask.is_inside(sfi.local, tol);

        // Compute some additional information if the intersection is valid
        if (sfi.status) {
            sfi.sf_desc = sf_desc;
            sfi.direction = !math::signbit(s);
            sfi.volume_link = mask.volume_link();
        }
    } else {
        // Not a valid intersection
        sfi.status = false;
    }
}

}  // namespace detray
