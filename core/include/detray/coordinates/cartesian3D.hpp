/** Detray library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/cartesian2D.hpp"

namespace detray {

/// @brief Frame projection into a 3D cartesian coordinate frame
template <typename algebra_t>
struct cartesian3D : cartesian2D<algebra_t> {
    // Import all constructors
    using cartesian2D<algebra_t>::cartesian2D;
};

}  // namespace detray
