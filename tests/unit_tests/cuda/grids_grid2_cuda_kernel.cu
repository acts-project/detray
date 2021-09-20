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

    /*--------------------------------------------------
      test function for the members of grid_data object
      --------------------------------------------------*/
    
    // test1 kernel declaration
    template <typename value_type, typename axis_p0_type, typename axis_p1_type>
    __global__ void grid_test1_kernel(
	  vecmem::data::vector_view<value_type> data_view,
	  const axis_p0_type axis0,
	  const axis_p1_type axis1);

    // test1 instantiation
    template void grid_test1<typename populator_t::store_value,
			     axis::regular<>,
			     axis::regular<>>(
	  vecmem::data::vector_view<populator_t::store_value> data_view,
	  const axis::regular<> axis0,
	  const axis::regular<> axis1);

    // test1 function implementation
    template< typename value_type,
	      typename axis_p0_type,
	      typename axis_p1_type>
    void grid_test1(vecmem::data::vector_view<value_type> data_view,
		    const axis_p0_type axis0,
		    const axis_p1_type axis1){

	int num_blocks = axis0.bins() * axis1.bins();

	int array_size = data_view.ptr()[0].size();
	
	int num_threads = array_size;
	grid_test1_kernel<<< num_blocks, num_threads >>>(data_view, axis0, axis1);

	// cuda error check
	DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
	DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
    }

    // test1 kernel implementation
    template <typename value_type, typename axis_p0_type, typename axis_p1_type>
    __global__ void grid_test1_kernel(
	  vecmem::data::vector_view<value_type> data_view,
	  const axis_p0_type axis0,
	  const axis_p1_type axis1){

	vecmem::device_vector<value_type> data_device(data_view);

	auto& pt = data_device[blockIdx.x][threadIdx.x];

	auto x_interval = (axis0.max - axis0.min)/axis0.n_bins;
	auto y_interval = (axis1.max - axis1.min)/axis1.n_bins;	

	auto gid = threadIdx.x + blockIdx.x * blockDim.x;

	pt = test::point3{axis0.min + gid*x_interval, axis1.min + gid*y_interval, 0.5};
    }

    /*------------------------------------
      test function for grid_data object
      ------------------------------------*/
    
    // test2 kernel declaration
    template <typename populator_type,
              typename axis_p0_type,
              typename axis_p1_type,
              typename serializer_type>
    __global__ void grid_test2_kernel(grid2_data<populator_type, axis_p0_type, axis_p1_type, serializer_type > grid_data);
    
    // test2 instantiation
    template void grid_test2<populator_t,
			     axis::regular<>,
			     axis::regular<>,
			     serializer2>(
        grid2_data<populator_t, axis::regular<>, axis::regular<>, serializer2>& grid_data);
    
    // test2 function implementation
    template <typename populator_type,
              typename axis_p0_type,
              typename axis_p1_type,
              typename serializer_type>    
    void grid_test2(grid2_data<populator_type, axis_p0_type, axis_p1_type, serializer_type >& grid_data){

	auto& data_view = grid_data._data_serialized;
	const auto& axis0 = grid_data._axis_p0;
	const auto& axis1 = grid_data._axis_p1;
	
	int num_blocks = axis0.bins() * axis1.bins();
	
	int num_threads = data_view.ptr()[0].size();
	
	grid_test2_kernel<<< num_blocks, num_threads >>>(grid_data);

	// cuda error check
	DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
	DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
    }

    // test2 kernel declaration
    template <typename populator_type,
              typename axis_p0_type,
              typename axis_p1_type,
              typename serializer_type>
    __global__ void grid_test2_kernel(grid2_data<populator_type, axis_p0_type, axis_p1_type, serializer_type > grid_data){

	vecmem::device_vector<typename populator_t::store_value> data_device(grid_data._data_serialized);
	const auto& axis0 = grid_data._axis_p0;
	const auto& axis1 = grid_data._axis_p1;

	auto& pt = data_device[blockIdx.x][threadIdx.x];

	auto x_interval = (axis0.max - axis0.min)/axis0.n_bins;
	auto y_interval = (axis1.max - axis1.min)/axis1.n_bins;	

	auto gid = threadIdx.x + blockIdx.x * blockDim.x;

	pt = test::point3{axis0.min + gid*x_interval, axis1.min + gid*y_interval, 0.5};	
    }
    
    /*--------------------------------------------------
      test function for the grid type template
      --------------------------------------------------*/

    // test3 kernel declaration
    template <typename grid_data_t>
    __global__ void grid_test3_kernel(grid_data_t grid_data);

    // test2 instantiation
    template void grid_test3<grid2r_complete_data>(
        grid2r_complete_data& grid_data);
    
    template <typename grid2_data_t>    
    void grid_test3(grid2_data_t& grid_data){

	auto& data_view = grid_data._data_serialized;
	const auto& axis0 = grid_data._axis_p0;
	const auto& axis1 = grid_data._axis_p1;
	
	int num_blocks = axis0.bins() * axis1.bins();
	
	int num_threads = data_view.ptr()[0].size();
	
	grid_test2_kernel<<< num_blocks, num_threads >>>(grid_data);

	// cuda error check
	DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
	DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
	
    }

    template <typename grid_data_t>
    __global__ void grid_test3_kernel(grid_data_t grid_data){
	vecmem::device_vector<typename grid_data_t::populator_t::store_value> data_device(grid_data._data_serialized);
	const auto& axis0 = grid_data._axis_p0;
	const auto& axis1 = grid_data._axis_p1;

	auto& pt = data_device[blockIdx.x][threadIdx.x];

	auto x_interval = (axis0.max - axis0.min)/axis0.n_bins;
	auto y_interval = (axis1.max - axis1.min)/axis1.n_bins;	

	auto gid = threadIdx.x + blockIdx.x * blockDim.x;

	pt = test::point3{axis0.min + gid*x_interval, axis1.min + gid*y_interval, 0.5};		
    }
    
} // namespace
