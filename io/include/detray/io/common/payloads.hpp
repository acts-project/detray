/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/grid_axis.hpp"

// System include(s)
#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

/// Raw indices (std::size_t) denote links between data components in different
/// files, while links used in detray detector objects are modelled as e.g.
/// @c single_link_payload
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
    unknown = 9u
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

/// @brief A payload for a single object link
struct single_link_payload {
    std::size_t link;
};

/// Geometry payloads
/// @{

/// @brief a payload for the geometry file header
struct geo_header_payload {
    std::string version, detector, tag, date;
    std::size_t n_volumes, n_surfaces;
};

/// @brief A payload for an affine transformation in homogeneous coordinates
struct transform_payload {
    std::array<real_io, 3u> tr;
    // Column major
    std::array<real_io, 9u> rot;
};

/// @brief A payload object for surface masks
struct mask_payload {
    using mask_shape = io::detail::mask_shape;
    mask_shape shape = mask_shape::unknown;
    single_link_payload volume_link;
    std::vector<real_io> boundaries;
};

/// @brief A payload object to link a surface to its material
struct material_link_payload {
    using material_type = io::detail::material_type;
    material_type type = material_type::unknown;
    std::size_t index;
};

/// @brief  A payload for surfaces
struct surface_payload {
    transform_payload transform;
    mask_payload mask;
    std::optional<material_link_payload> material;
    single_link_payload source;
    // Write the surface barcode as an additional information
    std::uint64_t barcode;
    detray::surface_id type = detray::surface_id::e_sensitive;
};

/// @brief A payload for volume bounds
struct volume_bounds_payload {
    detray::volume_id type = detray::volume_id::e_cylinder;
    std::vector<real_io> values = {};
};

/// @brief A payload object to link a volume to its acceleration data structures
struct acc_links_payload {
    using acc_type = io::detail::acc_type;
    acc_type type;
    std::size_t index;
};

/// @brief A payload for volumes
struct volume_payload {
    std::string name = "";
    volume_bounds_payload bounds;
    transform_payload transform;
    std::vector<surface_payload> surfaces;
    // Index of the volume in the detector volume container
    single_link_payload index;
    // Optional accelerator data structures
    std::optional<std::vector<acc_links_payload>> acc_links;
};

/// @}

/// Material payloads
/// @{

/// @brief a payload for the simple material file header
struct homogeneous_material_header_payload {
    std::string version, detector, tag, date;
    std::size_t n_slabs, n_rods;
};

/// @brief A payload object for a material parametrization
struct material_payload {
    std::array<real_io, 7u> params;
};

/// @brief A payload object for a material slab
struct material_slab_payload {
    using material_type = io::detail::material_type;
    material_type type = material_type::unknown;
    std::size_t index;
    real_io thickness;
    material_payload mat;
};

/// @brief A payload for a simple detector material description
struct detector_homogeneous_material_payload {
    std::vector<material_slab_payload> mat_slabs = {};
    std::optional<std::vector<material_slab_payload>> mat_rods;
};

/// @}

/// Payloads for a uniform grid
/// @{

/// @brief axis definition and bin edges
struct axis_payload {
    /// axis lookup type
    n_axis::binning binning = n_axis::binning::e_regular;
    n_axis::bounds bounds = n_axis::bounds::e_closed;
    n_axis::label label = n_axis::label::e_r;

    std::vector<real_io> edges = {};
    std::size_t bins = 0u;
};

/// @brief A payload for a grid definition
struct grid_payload {
    std::vector<axis_payload> axes = {};
    std::vector<std::vector<unsigned int>> entries = {};
};

/// @brief A payload for objects within a grid
struct grid_objects_payload {
    grid_payload grid;
    std::optional<transform_payload> transform;
};

/// @brief navigation links definition
struct links_payload {
    std::vector<single_link_payload> single_links;
    std::optional<grid_objects_payload> grid_links;
};

/// @}

/// @brief A payload for a detector
struct detector_payload {
    std::vector<volume_payload> volumes = {};
    grid_objects_payload volume_grid;
};

}  // namespace detray
