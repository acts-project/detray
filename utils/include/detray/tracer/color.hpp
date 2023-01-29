/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/definitions/math.hpp"

// System include(s)
#include <cmath>

namespace detray {

/// @brief holds rgb and alpha values for color shading
struct color : public point3D {
    scalar alpha{0.f};
};

}  // namespace detray
