/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <cmath>
#include <limits>

namespace detray::detail {

/// Class to solve a quadratic equation of type a * x^2 + b * x + c = 0
///
/// @note If there are no real solutions, the result is undefined
/// @note The solutions are sorted by default. If there is only one solution,
/// the larger value is undefined.
template <typename scalar_t>
class quadratic_equation {
    public:
    quadratic_equation() = delete;

    /// Solve the quadratic equation with the coefficients @param a, @param b
    /// and @param c
    ///
    /// @param tolerance threshhold to compare the discrimant against to decide
    ///                  if we have two separate solutions.
    DETRAY_HOST_DEVICE
    constexpr quadratic_equation(
        const scalar_t a, const scalar_t b, const scalar_t c,
        const scalar_t tolerance = std::numeric_limits<scalar_t>::epsilon()) {
        // linear case
        if (std::abs(a) <= tolerance) {
            m_solutions = 1;
            m_values[0] = -c / b;
        } else {
            const scalar_t discriminant{b * b - 4.f * a * c};
            // If there is more than one solution, then a != 0 and q != 0
            if (discriminant > tolerance) {
                m_solutions = 2;
                const scalar_t q{
                    -0.5f * (b + std::copysign(std::sqrt(discriminant), b))};
                m_values = {q / a, c / q};
                // Sort the two solutions
                if (m_values[0] > m_values[1]) {
                    m_values = {m_values[1], m_values[0]};
                }
            }
            // Only one solution and a != 0
            else if (discriminant >= 0.f) {
                m_solutions = 1;
                m_values[0] = -0.5f * b / a;
            }
            // discriminant < 0 is not allowed, since all solutions should be
            // real
        }
    }

    /// Getters for the solution(s)
    /// @{
    constexpr int solutions() const { return m_solutions; }
    constexpr scalar_t smaller() const { return m_values[0]; }
    constexpr scalar_t larger() const { return m_values[1]; }
    /// @}

    private:
    /// Number of solutions of the equation
    int m_solutions{0};
    /// The solutions
    std::array<scalar_t, 2> m_values{std::numeric_limits<scalar_t>::infinity(),
                                     std::numeric_limits<scalar_t>::infinity()};
};

}  // namespace detray::detail
