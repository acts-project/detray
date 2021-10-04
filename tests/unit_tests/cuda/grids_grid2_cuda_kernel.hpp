/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/populator.hpp"
#include "grids/serializer2.hpp"
#include "utils/indexing.hpp"

#pragma once

namespace detray {

static constexpr int n_points = 3;

using host_grid2_replace = grid2<host_replace_populator<test::point3>,
                                 axis::regular<>, axis::regular<>, serializer2>;

using device_grid2_replace =
    grid2<device_replace_populator<test::point3>, axis::regular<>,
          axis::regular<>, serializer2>;

using host_grid2_complete =
    grid2<host_complete_populator<n_points, false, test::point3>,
          axis::regular<>, axis::regular<>, serializer2>;

using device_grid2_complete =
    grid2<device_complete_populator<n_points, false, test::point3>,
          axis::regular<>, axis::regular<>, serializer2>;

using host_grid2_attach = grid2<host_attach_populator<false, test::point3>,
                                axis::circular<>, axis::regular<>, serializer2>;

using device_grid2_attach =
    grid2<device_attach_populator<false, test::point3>, axis::circular<>,
          axis::regular<>, serializer2>;

// test function for replace populator
void grid_replace_test(grid2_view<host_grid2_replace>& grid_view);

// test function for complete populator
void grid_complete_test(grid2_view<host_grid2_complete>& grid_view);

// read test function for grid with attach populator
void grid_attach_read_test(grid2_view<host_grid2_attach>& grid_view);

// fill test function for grid buffer with attach populator
void grid_attach_fill_test(grid2_view<host_grid2_attach>& grid_view);

}  // namespace detray
