/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

using real_io = double;

/// The following enums are defined per detector in the detector metadata
namespace io::detail {

/// Enumerate the shape primitives globally
enum class mask_shape : unsigned int {
    annulus2 = 0u,
    cuboid3 = 1u,
    cylinder2 = 2u,
    cylinder3 = 3u,
    line = 4u,
    rectangle2 = 5u,
    ring2 = 6u,
    single3 = 7u,
    trapezoid2 = 8u,
    square2 = 9u,
    n_shapes = 10u,
    unknown = n_shapes
};

/// Enumerate the different material types
enum class material_type : unsigned int {
    // Material texture (grid) shapes
    annulus2 = 0u,
    cuboid3 = 1u,
    cylinder2 = 2u,
    cylinder3 = 3u,
    line = 4u,
    rectangle2 = 5u,
    ring2 = 6u,
    single3 = 7u,
    trapezoid2 = 8u,
    unknown = 9u,
    // Simple materials
    slab = 11u,
    rod = 12u
};

/// Enumerate the different acceleration data structures
enum class acc_type : unsigned int {
    cyl_grid = 0u,
    disc_grid = 1u,
    unknown = 2u,
};

}  // namespace io::detail

}  // namespace detray
