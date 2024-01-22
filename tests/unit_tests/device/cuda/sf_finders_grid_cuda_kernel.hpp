/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray core
#include "detray/definitions/algebra.hpp"
#include "detray/masks/masks.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/grid_collection.hpp"
#include "detray/surface_finders/grid/populators.hpp"
#include "detray/surface_finders/grid/serializers.hpp"
#include "detray/tools/grid_builder.hpp"

namespace {

// type definitions
using size_type = __plugin::size_type;
using vector3 = __plugin::vector3<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;

}  // namespace

namespace detray {

using namespace n_axis;

static constexpr scalar tol{1e-7f};
static constexpr std::size_t n_points{3u};
static constexpr bool is_owning = true;

// Coordinate system definitions (all data-owning)

// 3D cartesian coordinate axes with closed bin boundaries and regular binning
template <typename containers = host_container_types>
using cartesian_3D = coordinate_axes<cuboid3D<>::axes<>, is_owning, containers>;

// 2D polar coordinate axes with closed bin boundaries, irregular binning in r
// and regular binning in phi
template <typename containers = host_container_types>
using polar_ir =
    coordinate_axes<ring2D<>::axes<bounds::e_closed, irregular, regular>,
                    is_owning, containers>;

// 2D polar coordinate axes with closed bin boundaries and regular binning
template <typename containers = host_container_types>
using polar = coordinate_axes<ring2D<>::axes<>, is_owning, containers>;

// 3D cylindrical coordinate axes with open bin boundaries and regular binning
// non-owning
template <typename containers = host_container_types>
using cylindrical_3D = coordinate_axes<cylinder3D::axes<bounds::e_open>,
                                       not is_owning, containers>;

// host and device grid definitions

// replacer
using host_grid3_single = grid<cartesian_3D<>, bins::single<point3>>;

using device_grid3_single =
    grid<cartesian_3D<device_container_types>, bins::single<point3>>;

using host_grid2_single_ci = grid<polar_ir<>, bins::single<point3>>;

using device_grid2_single_ci =
    grid<polar_ir<device_container_types>, bins::single<point3>>;

// completer/attacher
using host_grid2_array = grid<polar<>, bins::static_array<point3, n_points>>;

using device_grid2_array =
    grid<polar<device_container_types>, bins::static_array<point3, n_points>>;

// grid collection
using n_own_host_grid2_array =
    grid<cylindrical_3D<>, bins::static_array<dindex, n_points>>;

using n_own_device_grid2_array = grid<cylindrical_3D<device_container_types>,
                                      bins::static_array<dindex, n_points>>;

/// test function for replace populator
void grid_replace_test(host_grid3_single::view_type grid_view,
                       std::size_t dim_x, std::size_t dim_y, std::size_t dim_z);

/// test function for replace populator with circular and irregular axis
void grid_replace_ci_test(host_grid2_single_ci::view_type grid_view,
                          std::size_t dim_x, std::size_t dim_y);

// test function for complete populator
void grid_complete_test(host_grid2_array::view_type grid_view,
                        std::size_t dim_x, std::size_t dim_y);

// read test function for attach populator
void grid_attach_test(host_grid2_array::view_type grid_view, std::size_t dim_x,
                      std::size_t dim_y);

// read test function for attach populator
template <typename device_grid_t, typename view_t, typename... I>
void print_grid(view_t grid_view, I... dims);

// test function for a collection of grids
void grid_collection_test(
    grid_collection<n_own_host_grid2_array>::view_type grid_collection_view,
    vecmem::data::vector_view<dindex> n_bins_view,
    vecmem::data::vector_view<std::array<dindex, 3>> result_bins_view,
    std::size_t n_grids, std::size_t dim_x, std::size_t dim_y,
    std::size_t dim_z);

}  // namespace detray
