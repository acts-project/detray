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

// test function for replace populator
template <typename grid2_data_t>
void grid_replace_test(grid2_data_t& grid_data);

// test function for complete populator
template <typename grid2_data_t>
void grid_complete_test(grid2_data_t& grid_data);

// read test function for grid with attach populator
template <typename grid2_data_t>
void grid_attach_read_test(grid2_data_t& grid_data);

// fill test function for grid buffer with attach populator
template <typename grid2_data_t>
void grid_attach_fill_test(grid2_data_t& grid_data);

}  // namespace detray
