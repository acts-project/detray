/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/boolean.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray::soa::detail {

/// Class to solve a quadratic equation of type a * x^2 + b * x + c = 0
///
/// @note If there are no real solutions, the result is undefined
/// @note The solutions are sorted by default. If there is only one
/// solution, the larger value is undefined.
template <typename scalar_t>
class quadratic_equation {
    public:
    quadratic_equation() = delete;

    DETRAY_HOST_DEVICE
    constexpr quadratic_equation(const scalar_t &a, const scalar_t &b,
                                 const scalar_t &c,
                                 const scalar_t &tolerance = 1e-6f) {
        // Linear case
        auto one_sol = (math_ns::abs(a) <= tolerance);
        m_solutions(one_sol) = 1.f;
        m_values[0] = -c / b;

        // Early exit
        if (detray::detail::all_of(one_sol)) {
            return;
        }

        const scalar_t discriminant = b * b - (4.f * a) * c;

        const auto two_sol = (discriminant > tolerance);
        one_sol = !two_sol && (discriminant >= 0.f);

        // If there is more than one solution, then a != 0 and q != 0
        if (detray::detail::any_of(two_sol)) {
            m_solutions = 2.f;
            m_solutions.setZeroInverted(two_sol);

            const scalar_t q =
                -0.5f * (b + math_ns::copysign(math_ns::sqrt(discriminant), b));

            scalar_t first = q / a;
            scalar_t second = c / q;
            first.setZeroInverted(two_sol);
            second.setZeroInverted(two_sol);

            // Sort the solutions
            const auto do_swap = (second < first);
            if (detray::detail::all_of(do_swap)) {
                m_values = {second, first};
            } else if (detray::detail::none_of(do_swap)) {
                m_values = {first, second};
            } else {
                const auto tmp = second;
                second(do_swap) = first;
                first(do_swap) = tmp;
                m_values = {first, second};
            }
        }

        // Only one solution and a != 0
        if (detray::detail::any_of(one_sol)) {
            scalar_t sol = 1.f;
            scalar_t result = -0.5f * b / a;
            sol.setZeroInverted(one_sol);
            result.setZeroInverted(one_sol);

            m_solutions += sol;
            m_values[0] += result;
        }
        // discriminant < 0 is not allowed, since all solutions should
        // be real
    }

    /// Getters for the solution(s)
    /// @{
    constexpr const auto &solutions() const { return m_solutions; }
    constexpr const scalar_t &smaller() const { return m_values[0]; }
    constexpr const scalar_t &larger() const { return m_values[1]; }
    /// @}

    private:
    /// Number of solutions of the equation (needs to be floating point to
    /// apply the masks correctly)
    scalar_t m_solutions = 0.f;
    /// The solutions
    darray<scalar_t, 2> m_values{scalar_t(0.f), scalar_t(0.f)};
};

template <typename S>
quadratic_equation(const S a, const S &b, const S &c, const S &tolerance)
    -> quadratic_equation<S>;

}  // namespace detray::soa::detail
