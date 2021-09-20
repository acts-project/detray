/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "definitions/cuda_defs.hpp"
#include "grids_grid2_cuda_kernel.cuh"
#include <vecmem/containers/device_vector.hpp>

namespace detray{

    // test kernel declaration
    template <typename grid_data_t>
    __global__ void grid_test_kernel(grid_data_t grid_data);

    // test instantiation for complete populator
    template void grid_test<grid2r_complete_data>(
        grid2r_complete_data& grid_data);

    // test instantiation for attach populator
    template void grid_test<grid2r_attach_data>(
        grid2r_attach_data& grid_data);
    
    // test function implementation
    template <typename grid2_data_t>    
    void grid_test(grid2_data_t& grid_data){

	auto& data_view = grid_data._data_serialized;
	const auto& axis0 = grid_data._axis_p0;
	const auto& axis1 = grid_data._axis_p1;
	
	int num_blocks = axis0.bins() * axis1.bins();
	
	int num_threads = data_view.ptr()[0].size();
	
	grid_test_kernel<<< num_blocks, num_threads >>>(grid_data);

	// cuda error check
	DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
	DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
	
    }

    // test kernel implementation
    template <typename grid_data_t>
    __global__ void grid_test_kernel(grid_data_t grid_data){
	vecmem::device_vector<typename grid_data_t::store_value_t> data_device(grid_data._data_serialized);
	const auto& axis0 = grid_data._axis_p0;
	const auto& axis1 = grid_data._axis_p1;

	auto& bin_data = data_device[blockIdx.x];
	auto& pt = bin_data[threadIdx.x];

	auto x_interval = (axis0.max - axis0.min)/axis0.n_bins;
	auto y_interval = (axis1.max - axis1.min)/axis1.n_bins;	

	auto gid = threadIdx.x + blockIdx.x * blockDim.x;

	pt = test::point3{axis0.min + gid*x_interval, axis1.min + gid*y_interval, 0.5};		
    }
    
} // namespace
