/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "grid_types.hpp"

#pragma once

namespace detray{
    
    template< typename value_type,
	      typename axis_p0_type,
	      typename axis_p1_type>
    void grid_test1(vecmem::data::vector_view<value_type> data_view,
		    const axis_p0_type axis0,
		    const axis_p1_type axis1);
    
    template <typename populator_type,
              typename axis_p0_type,
              typename axis_p1_type,
              typename serializer_type>    
    void grid_test2(grid2_data<populator_type, axis_p0_type, axis_p1_type, serializer_type >& grid_data);

    template <typename grid2_data_t>    
    void grid_test3(grid2_data_t& grid_data);
    
    
} // namespace detray
