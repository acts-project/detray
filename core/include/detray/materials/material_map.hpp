/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/containers.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/tools/grid_factory.hpp"

namespace detray {

/// Definition of binned material
template <typename axes_shape, typename scalar_t,
          typename container_t = host_container_types, bool owning = false>
using material_map = grid<coordinate_axes<axes_shape, owning, container_t>,
                          material_slab<scalar_t>, simple_serializer, replacer>;

/// How to build material maps of various shapes
// TODO: Move to material_map_builder once available
template <typename scalar_t = detray::scalar>
using material_map_factory =
    grid_factory<material_slab<scalar_t>, simple_serializer, replacer>;

}  // namespace detray
