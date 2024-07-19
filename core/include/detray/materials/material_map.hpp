/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/builders/grid_factory.hpp"
#include "detray/definitions/detail/containers.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/detail/grid_bins.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/serializers.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// Definition of binned material
template <typename shape, typename scalar_t,
          typename container_t = host_container_types, bool owning = false>
using material_map = grid<axes<shape>, bins::single<material_slab<scalar_t>>,
                          simple_serializer, container_t, owning>;

/// How to build material maps of various shapes
// TODO: Move to material_map_builder once available
template <typename scalar_t = detray::scalar>
using material_grid_factory =
    grid_factory<bins::single<material_slab<scalar_t>>, simple_serializer>;

namespace detail {

// TODO: Define concepts

template <class material_t>
struct is_material_map<
    material_t,
    std::enable_if_t<
        is_grid_v<material_t> &&
            std::is_same_v<
                typename material_t::value_type,
                material_slab<typename material_t::value_type::scalar_type>>,
        void>> : public std::true_type {};

// Pick the 2D material map types up for surface material maps
template <typename material_t>
struct is_surface_material<
    material_t,
    std::enable_if_t<is_material_map_v<material_t> && (material_t::dim == 2),
                     void>> : public std::true_type {};

// Pick the 3D material map types up for volume material maps
template <typename material_t>
struct is_volume_material<
    material_t,
    std::enable_if_t<is_material_map_v<material_t> && (material_t::dim == 3),
                     void>> : public std::true_type {};

}  // namespace detail

}  // namespace detray
