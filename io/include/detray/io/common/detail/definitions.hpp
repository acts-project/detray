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
    portal_cylinder2 = 4u,
    rectangle2 = 5u,
    ring2 = 6u,
    trapezoid2 = 7u,
    cell_wire = 8u,
    straw_wire = 9u,
    single1 = 10u,
    single2 = 11u,
    single3 = 12u,
    n_shapes = 13u,
    unknown = n_shapes
};

/// Enumerate the different material types
enum class material_type : unsigned int {
    // Material texture (grid) shapes
    annulus2 = 0u,
    cuboid3 = 1u,
    cylinder2 = 2u,
    cylinder3 = 3u,
    rectangle2 = 4u,
    ring2 = 5u,
    trapezoid2 = 6u,
    cell_wire = 7u,
    straw_wire = 8u,
    // Homogeneous materials
    slab = 9u,
    rod = 10u,
    unknown = 11u
};

/// Enumerate the different acceleration data structures
enum class acc_type : unsigned int {
    brute_force = 0u,      // try all
    cartesian2_grid = 1u,  // rectangle, trapezoid, (triangle) grids
    cuboid3_grid = 2u,     // cuboid grid
    polar2_grid = 3u,      // ring/disc, annulus grids
    cylinder2_grid = 4u,   // 2D cylinder grid
    cylinder3_grid = 5u,   // 3D cylinder grid
    n_accel = 6u,
    unknown = n_accel
};

}  // namespace io::detail

}  // namespace detray
