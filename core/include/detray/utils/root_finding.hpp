/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <array>
#include <iostream>
#include <stdexcept>

namespace detray {

/// @brief Try to find a bracket around a root
///
/// @param a lower initial boundary
/// @param b upper initial boundary
/// @param f function for which to find the root
/// @param k scale factor with which to widen the bracket at every step
///
/// @see Numerical Recepies pp. 445
template <typename scalar_t, typename function_t>
DETRAY_HOST_DEVICE inline std::array<scalar_t, 2> bracket(
    scalar_t a, scalar_t b, function_t& f, const scalar_t k = 1.1f) {

    assert((k > 1.f) &&
           "Root bracketing: scale factor has to be larger than one");

    if (a == b) {
        throw std::invalid_argument(
            "Root bracketing: Not a valid start interval [" +
            std::to_string(a) + ", " + std::to_string(b) + "]");
    }

    // Sample function points at interval
    scalar_t f_a{f(a)};
    scalar_t f_b{f(b)};
    std::size_t n_tries{0u};

    // If there is no sign change in interval, we don't know if there is a root
    while (!math::signbit(f_a * f_b)) {
        if (math::fabs(f_a) < math::fabs(f_b)) {
            a += k * (a - b);
            f_a = f(a);
        } else {
            b += k * (b - a);
            f_b = f(b);
        }
        ++n_tries;
        // No interval could be found to bracket the root
        // Might be correct, if there is not root
        if (n_tries == 100u) {
#ifdef DEBUG
            std::cout << "WARNING: Could not bracket a root" << std::endl;
#endif
            break;
        }
    }

    return {a, b};
}

}  // namespace detray
