/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "utils/indexing.hpp"

#pragma once

namespace detray{

    using populator_t = complete_populator<3, false, test::point3>;
    
    using grid2r = grid2<populator_t, axis::regular<>, axis::regular<>, serializer2>;

    using grid2r_data = grid2_data<populator_t, axis::regular<>, axis::regular<>, serializer2>;
    
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
    

    
} // namespace detray
