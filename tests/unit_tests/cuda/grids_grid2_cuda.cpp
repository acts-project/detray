/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "grids_grid2_cuda_kernel.cuh"

#include <iostream>
#include <gtest/gtest.h>

#include <climits>

using namespace detray;

TEST(grids, grid2_complete_populator)
{
    // axis
    axis::regular<> xaxis{7, -1., 6.};
    axis::regular<> yaxis{3, 0., 3.};
    
    // managed memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // declare grid
    grid2r g2(std::move(xaxis), std::move(yaxis), test::point3{0,0,0}, &mng_mr);

    // get grid_data
    grid2_data g2_data(g2);

    // pre-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++){
	for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++){
	    
	    auto& data = g2.bin(i_x,i_y);
	    
	    for (auto pt: data){
		EXPECT_EQ(pt, g2.populator().kInvalid);
	    }	    
	}
    }
    
    // test grid_data in cuda
    grid_test(g2_data._data_serialized, g2_data._axis_p0, g2_data._axis_p1);

    auto x_interval = (xaxis.max - xaxis.min)/xaxis.n_bins;
    auto y_interval = (yaxis.max - yaxis.min)/yaxis.n_bins;	
    
    // post-check
    for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++){
	for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++){
	    
	    auto& data = g2.bin(i_x,i_y);
	    
	    for (int i_p = 0 ; i_p<data.size(); i_p++){
		auto& pt = data[i_p];

		auto bin_id = i_x + i_y * xaxis.bins();
		auto gid = i_p + bin_id * data.size() ;
		
		test::point3 tp({xaxis.min + gid*x_interval,
				 yaxis.min + gid*y_interval,
				 0.5});
		EXPECT_EQ(pt, tp);
	    }	    
	}
    }    
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
