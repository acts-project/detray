/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/coordinates.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/masks/annulus2D.hpp"
#include "detray/masks/cuboid3D.hpp"
#include "detray/masks/cylinder2D.hpp"
#include "detray/masks/cylinder3D.hpp"
#include "detray/masks/line.hpp"
#include "detray/masks/rectangle2D.hpp"
#include "detray/masks/ring2D.hpp"
#include "detray/masks/single3D.hpp"
#include "detray/masks/trapezoid2D.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/type_registry.hpp"

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

/// Infer the IO shape id from the shape type
template <typename shape_t>
constexpr mask_shape get_shape_id() {

    /// Register the mask shapes to the @c mask_shape enum
    using shape_registry =
        type_registry<mask_shape, annulus2D<>, cuboid3D<>, cylinder2D<>,
                      cylinder3D,
                      cylinder2D<false, cylinder_portal_intersector>,
                      rectangle2D<>, ring2D<>, trapezoid2D<>, line<true>,
                      line<false>, single3D<0>, single3D<1>, single3D<2>>;

    // Find the correct shape IO id;
    if constexpr (shape_registry::is_defined(shape_t{})) {
        return shape_registry::get_id(shape_t{});
    } else {
        return mask_shape::unknown;
    }
}

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

/// Infer the IO material id from the material type
template <typename material_t>
constexpr material_type get_material_id() {
    using scalar_t = typename material_t::scalar_type;

    /// Register the material types to the @c material_type enum
    using mat_registry =
        type_registry<material_type, annulus2D<>, cuboid3D<>, cylinder2D<>,
                      cylinder3D, rectangle2D<>, ring2D<>, trapezoid2D<>,
                      line<true>, line<false>, material_slab<scalar_t>,
                      material_rod<scalar_t>>;

    // Find the correct material IO id;
    if constexpr (mat_registry::is_defined(material_t{})) {
        return mat_registry::get_id(material_t{});
    } else {
        return material_type::unknown;
    }
}

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

/// Infer the grid id from its coordinate system
template <typename grid_t>
constexpr acc_type get_grid_id() {

    using frame_t = typename grid_t::local_frame_type;
    using algebra_t = typename frame_t::algebra;

    /// Register the grid shapes to the @c acc_type enum
    /// @note the first type corresponds to a non-grid type in the enum
    /// (brute force)
    using frame_registry =
        type_registry<acc_type, void, cartesian2D<algebra_t>,
                      cartesian3D<algebra_t>, polar2D<algebra_t>,
                      cylindrical2D<algebra_t>, cylindrical3D<algebra_t>>;

    // Find the correct grid shape IO id;
    if constexpr (frame_registry::is_defined(frame_t{})) {
        return frame_registry::get_id(frame_t{});
    } else {
        return acc_type::unknown;
    }
}

}  // namespace io::detail

}  // namespace detray
