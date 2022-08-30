/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"

// System include(s).
#include <cmath>

namespace detray::detail {

template <typename scalar_t>
struct trigonometrics {
    scalar_t cos_theta;
    scalar_t sin_theta;
    scalar_t cos_phi;
    scalar_t sin_phi;
};

}  // namespace detray::detail