/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/coordinates/coordinates.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/materials/detail/concepts.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/navigation/accelerators/concepts.hpp"
#include "detray/utils/tuple_helpers.hpp"
#include "detray/utils/type_registry.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::io::detail {

/// Infer the IO shape id from the shape type
template <typename shape_t>
    requires std::is_enum_v<typename shape_t::boundaries>
constexpr io::shape_id get_id() {

    /// Register the mask shapes to the @c shape_id enum
    using shape_registry =
        types::registry<io::shape_id, annulus2D, cuboid3D, cylinder2D,
                        cylinder3D, concentric_cylinder2D, rectangle2D, ring2D,
                        trapezoid2D, line_square, line_circular, single3D<0>,
                        single3D<1>, single3D<2>>;

    // Find the correct shape IO id;
    if constexpr (types::contains<shape_registry, shape_t>) {
        return types::id<shape_registry, shape_t>;
    } else {
        return io::shape_id::unknown;
    }
}

/// Infer the IO material id from the material type - homogeneous material
template <detray::concepts::homogeneous_material material_t>
constexpr io::material_id get_id() {
    using scalar_t = typename material_t::scalar_type;

    /// Register the material types to the @c material_id enum
    using mat_registry =
        types::registry<io::material_id, void, void, void, void, void, void,
                        material_slab<scalar_t>, material_rod<scalar_t>,
                        material<scalar_t>>;

    // Find the correct material IO id;
    if constexpr (types::contains<mat_registry, material_t>) {
        return types::id<mat_registry, material_t>;
    } else {
        return io::material_id::unknown;
    }
}

/// Infer the IO material id from the material type - material maps
template <detray::concepts::material_map material_t>
constexpr io::material_id get_id() {

    using map_frame_t = typename material_t::local_frame_type;
    using algebra_t = typename map_frame_t::algebra_type;

    /// Register the material types to the @c material_id enum
    using mat_registry = types::registry<
        io::material_id, polar2D<algebra_t>, cartesian2D<algebra_t>,
        cartesian3D<algebra_t>, concentric_cylindrical2D<algebra_t>,
        cylindrical2D<algebra_t>, cylindrical3D<algebra_t>, void, void>;

    // Find the correct material IO id;
    if constexpr (types::contains<mat_registry, map_frame_t>) {
        return types::id<mat_registry, map_frame_t>;
    } else {
        return io::material_id::unknown;
    }
}

/// Infer the grid id from its coordinate system
template <detray::concepts::surface_grid grid_t>
constexpr io::accel_id get_id() {

    using frame_t = typename grid_t::local_frame_type;
    using algebra_t = typename frame_t::algebra_type;

    /// Register the grid shapes to the @c accel_id enum
    /// @note the first type corresponds to a non-grid type in the enum
    /// (brute force)
    using frame_registry =
        types::registry<io::accel_id, void, cartesian2D<algebra_t>,
                        cartesian3D<algebra_t>, polar2D<algebra_t>,
                        concentric_cylindrical2D<algebra_t>,
                        cylindrical2D<algebra_t>, cylindrical3D<algebra_t>>;

    // Find the correct grid shape IO id;
    if constexpr (types::contains<frame_registry, frame_t>) {
        return types::id<frame_registry, frame_t>;
    } else {
        return io::accel_id::unknown;
    }
}

/// Infer the grid id from its coordinate system
template <detray::concepts::volume_grid grid_t>
constexpr io::accel_id get_id() {

    using frame_t = typename grid_t::local_frame_type;
    using algebra_t = typename frame_t::algebra_type;

    /// Register the grid shapes to the @c accel_id enum
    /// @note the first type corresponds to a non-grid type in the enum
    /// (brute force)
    using frame_registry =
        types::registry<io::accel_id, void, cartesian2D<algebra_t>,
                        cartesian3D<algebra_t>, polar2D<algebra_t>,
                        concentric_cylindrical2D<algebra_t>,
                        cylindrical2D<algebra_t>, cylindrical3D<algebra_t>>;

    // Find the correct grid shape IO id;
    if constexpr (types::contains<frame_registry, frame_t>) {
        return types::id<frame_registry, frame_t>;
    } else {
        return io::accel_id::unknown;
    }
}

/// Determine the type and id of a shape of a mask without triggering a compiler
/// error (sfinae) if the detector does not know the type / enum entry
/// @{
/// Mask shape unknown by detector
template <io::shape_id shape, typename detector_t>
struct mask_info {
    using shape_id = typename detector_t::masks::id;
    using type = void;
    static constexpr shape_id value{detray::detail::invalid_value<shape_id>()};
};

/// Check for a stereo annulus shape
template <typename detector_t>
    requires(types::contains<typename detector_t::masks,
                             mask<annulus2D, typename detector_t::algebra_type,
                                  std::uint_least16_t>>)
struct mask_info<io::shape_id::annulus2, detector_t> {
    using type = annulus2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_annulus2D};
};

/// Check for a 2D cylinder shape
template <typename detector_t>
    requires(types::contains<typename detector_t::masks,
                             mask<cylinder2D, typename detector_t::algebra_type,
                                  std::uint_least16_t>>)
struct mask_info<io::shape_id::cylinder2, detector_t> {
    using type = cylinder2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_cylinder2D};
};

/// Check for a 2D cylinder portal shape
template <typename detector_t>
    requires(types::contains<
             typename detector_t::masks,
             mask<concentric_cylinder2D, typename detector_t::algebra_type,
                  std::uint_least16_t>>)
struct mask_info<io::shape_id::portal_cylinder2, detector_t> {
    using type = concentric_cylinder2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_concentric_cylinder2D};
};

/// Check for a cell wire line shape
template <typename detector_t>
    requires(
        types::contains<typename detector_t::masks,
                        mask<line_square, typename detector_t::algebra_type,
                             std::uint_least16_t>>)
struct mask_info<io::shape_id::drift_cell, detector_t> {
    using type = line_square;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_drift_cell};
};

/// Check for a straw wire line shape
template <typename detector_t>
    requires(
        types::contains<typename detector_t::masks,
                        mask<line_circular, typename detector_t::algebra_type,
                             std::uint_least16_t>>)
struct mask_info<io::shape_id::straw_tube, detector_t> {
    using type = line_circular;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_straw_tube};
};

/// Check for a rectangle shape
template <typename detector_t>
    requires(
        types::contains<typename detector_t::masks,
                        mask<rectangle2D, typename detector_t::algebra_type,
                             std::uint_least16_t>>)
struct mask_info<io::shape_id::rectangle2, detector_t> {
    using type = rectangle2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_rectangle2D};
};

/// Check for a ring/disc shape
template <typename detector_t>
    requires(types::contains<typename detector_t::masks,
                             mask<ring2D, typename detector_t::algebra_type,
                                  std::uint_least16_t>>)
struct mask_info<io::shape_id::ring2, detector_t> {
    using type = ring2D;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_ring2D};
};

/// Check for a single masked value (1st value is checked)
template <typename detector_t>
    requires(
        types::contains<typename detector_t::masks,
                        mask<single3D<0>, typename detector_t::algebra_type,
                             std::uint_least16_t>>)
struct mask_info<io::shape_id::single1, detector_t> {
    using type = single3D<0>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single1D};
};

/// Check for a single masked value (2nd value is checked)
template <typename detector_t>
    requires(
        types::contains<typename detector_t::masks,
                        mask<single3D<1>, typename detector_t::algebra_type,
                             std::uint_least16_t>>)
struct mask_info<io::shape_id::single2, detector_t> {
    using type = single3D<1>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single2D};
};

/// Check for a single masked value (3rd value is checked)
template <typename detector_t>
    requires(
        types::contains<typename detector_t::masks,
                        mask<single3D<2>, typename detector_t::algebra_type,
                             std::uint_least16_t>>)
struct mask_info<io::shape_id::single3, detector_t> {
    using type = single3D<2>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single3D};
};

/// Check for a trapezoid shape
template <typename detector_t>
    requires(
        types::contains<typename detector_t::masks,
                        mask<trapezoid2D, typename detector_t::algebra_type,
                             std::uint_least16_t>>)
struct mask_info<io::shape_id::trapezoid2, detector_t> {
    using type = trapezoid2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_trapezoid2D};
};
/// @}

/// Determine the type and id of a material map without triggering a compiler
/// error (sfinae) if the detector does not know the type / enum entry
/// @{
/// Material map unknown by detector
template <io::material_id mat, typename detector_t>
struct mat_map_info {
    using material_id = typename detector_t::materials::id;
    using type = void;
    static constexpr material_id value{
        detray::detail::invalid_value<material_id>()};
};

/// Check for a 2D disc material map
template <typename detector_t>
    requires(types::contains<
             typename detector_t::materials,
             material_map<typename detector_t::algebra_type, ring2D>>)
struct mat_map_info<io::material_id::ring2_map, detector_t> {
    using type = material_map<typename detector_t::algebra_type, ring2D>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_ring2D_map};
};

/// Check for a 2D cartesian material map
template <typename detector_t>
    requires(types::contains<
             typename detector_t::materials,
             material_map<typename detector_t::algebra_type, rectangle2D>>)
struct mat_map_info<io::material_id::rectangle2_map, detector_t> {
    using type = material_map<typename detector_t::algebra_type, rectangle2D>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_rectangle2D_map};
};

/// Check for a 3D cuboid volume material map
template <typename detector_t>
    requires(types::contains<
             typename detector_t::materials,
             material_map<typename detector_t::algebra_type, cuboid3D>>)
struct mat_map_info<io::material_id::cuboid3_map, detector_t> {
    using type = material_map<typename detector_t::algebra_type, cuboid3D>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cuboid3D_map};
};

/// Check for a 2D cylindrical material map
template <typename detector_t>
    requires(types::contains<
             typename detector_t::materials,
             material_map<typename detector_t::algebra_type, cylinder2D>>)
struct mat_map_info<io::material_id::cylinder2_map, detector_t> {
    using type = material_map<typename detector_t::algebra_type, cylinder2D>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cylinder2D_map};
};

/// Check for a 2D concentric cylindrical material map
template <typename detector_t>
    requires(types::contains<typename detector_t::materials,
                             material_map<typename detector_t::algebra_type,
                                          concentric_cylinder2D>>)
struct mat_map_info<io::material_id::concentric_cylinder2_map, detector_t> {
    using type =
        material_map<typename detector_t::algebra_type, concentric_cylinder2D>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_concentric_cylinder2D_map};
};

/// Check for a 3D cylindrical volume material map
template <typename detector_t>
    requires(types::contains<
             typename detector_t::materials,
             material_map<typename detector_t::algebra_type, cylinder3D>>)
struct mat_map_info<io::material_id::cylinder3_map, detector_t> {
    using type = material_map<typename detector_t::algebra_type, cylinder3D>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cylinder3D_map};
};
/// @}

}  // namespace detray::io::detail
