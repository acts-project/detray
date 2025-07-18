/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <detray/definitions/detail/qualifiers.hpp>
#include <ostream>
namespace detray {

using real_io = double;

/// The following enums are defined per detector in the detector metadata
namespace io {

enum class format { json = 0u };

/// Enumerate the shape primitives globally
enum class shape_id : unsigned int {
    annulus2 = 0u,
    cuboid3 = 1u,
    cylinder2 = 2u,
    cylinder3 = 3u,
    portal_cylinder2 = 4u,
    rectangle2 = 5u,
    ring2 = 6u,
    trapezoid2 = 7u,
    drift_cell = 8u,
    straw_tube = 9u,
    single1 = 10u,
    single2 = 11u,
    single3 = 12u,
    n_shapes = 13u,
    unknown = n_shapes
};

/// Enumerate the different material types
enum class material_id : unsigned int {
    // Material texture (grid) shapes
    annulus2_map = 0u,
    rectangle2_map = 1u,
    cuboid3_map = 2u,
    concentric_cylinder2_map = 3u,
    cylinder2_map = 4u,
    cylinder3_map = 5u,
    ring2_map = 0u,
    trapezoid2_map = 1u,
    // Homogeneous materials
    slab = 6u,
    rod = 7u,
    raw_material = 8u,
    n_mats = 9u,
    unknown = n_mats
};

/// Enumerate the different acceleration data structures
enum class accel_id : unsigned int {
    brute_force = 0u,                // try all
    cartesian2_grid = 1u,            // rectangle, trapezoid, (triangle) grids
    cuboid3_grid = 2u,               // cuboid grid
    polar2_grid = 3u,                // ring/disc, annulus grids
    concentric_cylinder2_grid = 4u,  // 2D concentric cylinder grid
    cylinder2_grid = 5u,             // 2D cylinder grid
    cylinder3_grid = 6u,             // 3D cylinder grid
    n_accel = 7u,
    unknown = n_accel
};

#define _enum_print(x) \
    case x:            \
        os << #x;      \
        break

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, shape_id sid) {
    switch (sid) {
        using enum shape_id;
        _enum_print(annulus2);
        _enum_print(cuboid3);
        _enum_print(cylinder2);
        _enum_print(cylinder3);
        _enum_print(portal_cylinder2);
        _enum_print(rectangle2);
        _enum_print(ring2);
        _enum_print(trapezoid2);
        _enum_print(drift_cell);
        _enum_print(straw_tube);
        _enum_print(single1);
        _enum_print(single2);
        _enum_print(single3);
        case unknown:
            // n_shapes has same value (13u)
            os << "unknown/n_shapes";
            break;
    }
    return os;
}

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, material_id mid) {
    switch (mid) {
        using enum material_id;
        case annulus2_map:
            // ring2_map has same value (0u)
            os << "annulus2_map/ring2_map";
            break;
        case rectangle2_map:
            // trapezoid2_map has same value (1u)
            os << "rectangle2_map/trapezoid2_map";
            break;
            _enum_print(cuboid3_map);
            _enum_print(concentric_cylinder2_map);
            _enum_print(cylinder2_map);
            _enum_print(cylinder3_map);
            _enum_print(slab);
            _enum_print(rod);
            _enum_print(raw_material);
        case unknown:
            // n_mats has same value (9u)
            os << "unknown/n_mats";
            break;
    }
    return os;
}

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, accel_id aid) {
    switch (aid) {
        using enum accel_id;
        _enum_print(brute_force);
        _enum_print(cartesian2_grid);
        _enum_print(cuboid3_grid);
        _enum_print(polar2_grid);
        _enum_print(concentric_cylinder2_grid);
        _enum_print(cylinder2_grid);
        _enum_print(cylinder3_grid);
        _enum_print(unknown);
    }
    return os;
}

#undef _enum_print

}  // namespace io

}  // namespace detray
