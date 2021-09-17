/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "definitions/cuda_defs.hpp"
#include "grids_grid2_cuda_kernel.cuh"

namespace detray{

    template <typename value_type, typename axis_p0_type, typename axis_p1_type>
    __global__ void grid_test_kernel(
	  vecmem::data::vector_view<value_type> data_view,
	  const axis_p0_type axis0,
	  const axis_p1_type axis1);

    // instantiation
    template void grid_test<typename populator_t::store_value,
			    axis::regular<>,
			    axis::regular<>>(
	  vecmem::data::vector_view<populator_t::store_value> data_view,
	  const axis::regular<> axis0,
	  const axis::regular<> axis1);
    
    template< typename value_type,
	      typename axis_p0_type,
	      typename axis_p1_type>
    void grid_test(vecmem::data::vector_view<value_type> data_view,
		   const axis_p0_type axis0,
		   const axis_p1_type axis1){

	int num_blocks = axis0.bins() * axis1.bins();

	int array_size = data_view.ptr()[0].size();
	
	int num_threads = array_size;
	grid_test_kernel<<< num_blocks, num_threads >>>(data_view, axis0, axis1);

	// cuda error check
	DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
	DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    }


    template <typename value_type, typename axis_p0_type, typename axis_p1_type>
    __global__ void grid_test_kernel(
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
    
    template void test1<std::vector, int>(std::vector<int>& vec);
    
    template < template < typename > class vector_t, typename T>
    void test1(vector_t<T>& vec){
	vector_t<T> a;
	printf("hi \n");
    }
    
    
} // namespace
