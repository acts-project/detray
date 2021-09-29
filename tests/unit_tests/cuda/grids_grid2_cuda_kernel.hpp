/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "grids_types.hpp"

#pragma once

namespace detray {

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
