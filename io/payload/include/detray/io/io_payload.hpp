/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <vector>

namespace detray {

using real_io = float;

/// @brief  A small pa
struct transform_payload {
    std::array<real_io, 3u> tr;
    std::array<real_io, 9u> rot;
};

/// @brief axis definition
struct axis_payload {
    /// axis lookup type
    enum class axis_label : unsigned int {
        x = 0u,
        y = 1u,
        z = 2u,
        r = 3u,
        phi = 4u
    };
    /// How the axis is done
    enum class axis_binning : unsigned int { equidistant = 0u, variable = 1u };
    /// How the axis is bound
    enum class axis_bounds : unsigned int { closed = 0u, circular = 1u };

    axis_binning binning = axis_binning::equidistant;
    axis_bounds bounds = axis_bounds::closed;
    axis_label label = axis_label::r;

    std::vector<real_io> edges = {};
    std::size_t bins = 0u;
};

/// @brief A payload for a grid definition
struct grid_payload {
    std::vector<axis_payload> axes = {};
    std::vector<std::vector<unsigned int>> entries = {};
};

/// @brief A payload object for material
struct material_slab_payload {
    std::array<real_io, 5u> slab;
};

/// @brief A payload for a single object link
struct single_object_payload {
    unsigned int link;
};

/// @brief A payload for objects within a grid
struct grid_objects_payload {
    grid_payload grid;
    std::optional<transform_payload> transform;
};

/// @brief navigation links definition
struct links_payload {
    std::vector<single_object_payload> single_links;
    std::optional<grid_objects_payload> grid_links;
};

/// @brief A payload object for surface masks
struct mask_payload {
    enum class mask_shape : unsigned int {
        annulus2 = 0u,
        cuboid3 = 1u,
        cylinder2 = 2u,
        cylinder3 = 3u,
        line = 4u,
        rectangle2 = 5u,
        ring2 = 6u,
        single3 = 7u,
        trapezoid2 = 8u
    };
    mask_shape shape = mask_shape::ring2;
    std::vector<real_io> boundaries;
};

/// @brief  A payload for surfaces
struct surface_payload {
    transform_payload transform;
    mask_payload mask;
    material_slab_payload material;
    std::size_t gid;
};

/// @brief A payload for portals
struct portal_payload {
    surface_payload surface;
    links_payload volume_links;
};

/// @brief A payload for volume bounds
struct volume_bounds_payload {
    enum class volume_bounds_type {
        cuboid = 0u,
        cylindrical = 1u,
        generic_cuboid = 2u,
        trapezoid = 3u
    };
    volume_bounds_type type = volume_bounds_type::cylindrical;
    std::vector<real_io> values = {};
};

/// @brief A payload for volumes
struct volume_payload {
    std::string name = "";
    transform_payload transform;
    volume_bounds_payload volume_bounds;
    std::vector<portal_payload> portals;
    std::vector<surface_payload> surfaces;
    links_payload surface_links;
};

/// @brief A payload for a detector
struct detector_payload {
    std::string name = "";
    std::vector<volume_payload> volumes = {};
    grid_objects_payload volume_grid;
};

}  // namespace detray
