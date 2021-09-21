/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "grids_types.hpp"

#pragma once

namespace detray{

    // test function for replace populator
    template <typename grid2_data_t>    
    void grid_test1(grid2_data_t& grid_data);    
    
    // test function for complete and attach populator
    template <typename grid2_data_t>    
    void grid_test2(grid2_data_t& grid_data);    
    
} // namespace detray
