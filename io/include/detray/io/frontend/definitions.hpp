/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
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
namespace io {

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
    cell_wire = 8u,
    straw_wire = 9u,
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
    cylinder2_map = 3u,
    cylinder3_map = 4u,
    ring2_map = 0u,
    trapezoid2_map = 1u,
    // Homogeneous materials
    slab = 5u,
    rod = 6u,
    n_mats = 7u,
    unknown = n_mats
};

/// Enumerate the different acceleration data structures
enum class accel_id : unsigned int {
    brute_force = 0u,      // try all
    cartesian2_grid = 1u,  // rectangle, trapezoid, (triangle) grids
    cuboid3_grid = 2u,     // cuboid grid
    polar2_grid = 3u,      // ring/disc, annulus grids
    cylinder2_grid = 4u,   // 2D cylinder grid
    cylinder3_grid = 5u,   // 3D cylinder grid
    n_accel = 6u,
    unknown = n_accel
};

namespace detail {

/// Infer the IO shape id from the shape type
template <
    typename shape_t,
    std::enable_if_t<std::is_enum_v<typename shape_t::boundaries>, bool> = true>
constexpr io::shape_id get_id() {

    /// Register the mask shapes to the @c shape_id enum
    using shape_registry =
        type_registry<io::shape_id, annulus2D<>, cuboid3D<>, cylinder2D<>,
                      cylinder3D,
                      cylinder2D<false, cylinder_portal_intersector>,
                      rectangle2D<>, ring2D<>, trapezoid2D<>, line<true>,
                      line<false>, single3D<0>, single3D<1>, single3D<2>>;

    // Find the correct shape IO id;
    if constexpr (shape_registry::is_defined(shape_t{})) {
        return shape_registry::get_id(shape_t{});
    } else {
        return io::shape_id::unknown;
    }
}

/// Infer the IO material id from the material type - homogeneous material
template <
    typename material_t,
    std::enable_if_t<
        std::is_base_of_v<detray::detail::homogeneous_material_tag, material_t>,
        bool> = true>
constexpr io::material_id get_id() {
    using scalar_t = typename material_t::scalar_type;

    /// Register the material types to the @c material_id enum
    using mat_registry =
        type_registry<io::material_id, void, void, void, void, void,
                      material_slab<scalar_t>, material_rod<scalar_t>>;

    // Find the correct material IO id;
    if constexpr (mat_registry::is_defined(material_t{})) {
        return mat_registry::get_id(material_t{});
    } else {
        return io::material_id::unknown;
    }
}

/// Infer the IO material id from the material type - material maps
template <typename material_t,
          std::enable_if_t<
              std::is_same_v<typename material_t::value_type,
                             material_slab<typename material_t::scalar_type>>,
              bool> = true>
constexpr io::material_id get_id() {

    using map_frame_t = typename material_t::local_frame_type;
    using algebra_t = typename map_frame_t::transform3_type;

    /// Register the material types to the @c material_id enum
    using mat_registry =
        type_registry<io::material_id, polar2<algebra_t>, cartesian2<algebra_t>,
                      cartesian3<algebra_t>, cylindrical2<algebra_t>,
                      cylindrical3<algebra_t>, void, void>;

    // Find the correct material IO id;
    if constexpr (mat_registry::is_defined(map_frame_t{})) {
        return mat_registry::get_id(map_frame_t{});
    } else {
        return io::material_id::unknown;
    }
}

/// Infer the grid id from its coordinate system
template <typename grid_t,
          std::enable_if_t<
              !std::is_same_v<typename grid_t::value_type,
                              material_slab<typename grid_t::scalar_type>>,
              bool> = true>
constexpr io::accel_id get_id() {

    using frame_t = typename grid_t::local_frame_type;
    using algebra_t = typename frame_t::transform3_type;

    /// Register the grid shapes to the @c accel_id enum
    /// @note the first type corresponds to a non-grid type in the enum
    /// (brute force)
    using frame_registry =
        type_registry<io::accel_id, void, cartesian2<algebra_t>,
                      cartesian3<algebra_t>, polar2<algebra_t>,
                      cylindrical2<algebra_t>, cylindrical3<algebra_t>>;

    // Find the correct grid shape IO id;
    if constexpr (frame_registry::is_defined(frame_t{})) {
        return frame_registry::get_id(frame_t{});
    } else {
        return io::accel_id::unknown;
    }
}

}  // namespace detail

}  // namespace io

}  // namespace detray
