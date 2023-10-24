/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

enum class shape_types {
    annulus2D,
    cuboid3D,
    cylinder2D,
    cylinder3D,
    line,
    rectangle2D,
    ring2D,
    single3D,
    trapezoid2D,
    unbounded,
    unmasked,
};
}  // namespace detray