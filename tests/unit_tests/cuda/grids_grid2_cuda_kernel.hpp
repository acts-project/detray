/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/utils/indexing.hpp"

#pragma once

namespace detray {

static constexpr int n_points = 3;

using host_grid2_replace =
    grid2<replace_populator, axis::regular<>, axis::regular<>, serializer2,
          vecmem::vector, vecmem::jagged_vector, darray, std::tuple,
          test::point3>;

using device_grid2_replace =
    grid2<replace_populator, axis::regular<>, axis::regular<>, serializer2,
          vecmem::device_vector, vecmem::jagged_device_vector, darray,
          std::tuple, test::point3>;

using host_grid2_complete =
    grid2<complete_populator, axis::regular<>, axis::regular<>, serializer2,
          vecmem::vector, vecmem::jagged_vector, darray, std::tuple,
          test::point3, n_points, false>;

using device_grid2_complete =
    grid2<complete_populator, axis::regular<>, axis::regular<>, serializer2,
          vecmem::device_vector, vecmem::jagged_device_vector, darray,
          std::tuple, test::point3, n_points, false>;

using host_grid2_attach =
    grid2<attach_populator, axis::circular<>, axis::regular<>, serializer2,
          vecmem::vector, vecmem::jagged_vector, darray, std::tuple,
          test::point3, false>;

using device_grid2_attach =
    grid2<attach_populator, axis::circular<>, axis::regular<>, serializer2,
          vecmem::device_vector, vecmem::jagged_device_vector, darray,
          std::tuple, test::point3, false>;

// test function for replace populator
void grid_replace_test(grid2_view<host_grid2_replace> grid_view);

// test function for complete populator
void grid_complete_test(grid2_view<host_grid2_complete> grid_view);

// read test function for grid with attach populator
void grid_attach_read_test(grid2_view<host_grid2_attach> grid_view);

// fill test function for grid buffer with attach populator
void grid_attach_fill_test(grid2_view<host_grid2_attach> grid_view);

}  // namespace detray
