/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/containers/device_vector.hpp>

#include "definitions/cuda_defs.hpp"
#include "grids_grid2_cuda_kernel.cuh"

namespace detray {

/*------------------------------------
  test function for replace populator
  ------------------------------------*/

// test1 kernel declaration
template <typename grid_data_t>
__global__ void grid_test1_kernel(grid_data_t grid_data);

// test1 instantiation for replace populator
template void grid_test1<grid2r_replace_data>(grid2r_replace_data& grid_data);

// test2 function implementation
template <typename grid2_data_t>
void grid_test1(grid2_data_t& grid_data) {

    // auto& data_view = grid_data._data_serialized;
    const auto& axis0 = grid_data._axis_p0;
    const auto& axis1 = grid_data._axis_p1;

    int num_blocks = 1;
    int num_threads = axis0.bins() * axis1.bins();

    // run the kernel
    grid_test1_kernel<<<num_blocks, num_threads>>>(grid_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// test2 kernel implementation
template <typename grid_data_t>
__global__ void grid_test1_kernel(grid_data_t grid_data) {

    typename grid_data_t::populator_t::device_vector_t data_device(
        grid_data._data_serialized);
    const auto& axis0 = grid_data._axis_p0;
    const auto& axis1 = grid_data._axis_p1;

    auto& pt = data_device[threadIdx.x];

    auto x_interval = (axis0.max - axis0.min) / axis0.n_bins;
    auto y_interval = (axis1.max - axis1.min) / axis1.n_bins;

    pt = test::point3{axis0.min + threadIdx.x * x_interval,
                      axis1.min + threadIdx.x * y_interval, 0.5};
}

/*------------------------------------------------
  test function for complete and attach populator
  ------------------------------------------------*/

// test2 kernel declaration
template <typename grid_data_t>
__global__ void grid_test2_kernel(grid_data_t grid_data);

// test2 instantiation for complete populator
template void grid_test2<grid2r_complete_data>(grid2r_complete_data& grid_data);

// test2 instantiation for attach populator
template void grid_test2<grid2r_attach_data>(grid2r_attach_data& grid_data);

// test2 function implementation
template <typename grid2_data_t>
void grid_test2(grid2_data_t& grid_data) {

    // auto& data_view = grid_data._data_serialized;
    const auto& axis0 = grid_data._axis_p0;
    const auto& axis1 = grid_data._axis_p1;

    int num_blocks = axis0.bins() * axis1.bins();
    int num_threads = n_points;

    // run the kernel
    grid_test2_kernel<<<num_blocks, num_threads>>>(grid_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// test2 kernel implementation
template <typename grid_data_t>
__global__ void grid_test2_kernel(grid_data_t grid_data) {

    typename grid_data_t::populator_t::device_vector_t data_device(
        grid_data._data_serialized);
    const auto& axis0 = grid_data._axis_p0;
    const auto& axis1 = grid_data._axis_p1;

    auto& bin_data = data_device[blockIdx.x][threadIdx.x];
    auto& pt = bin_data;

    auto x_interval = (axis0.max - axis0.min) / axis0.n_bins;
    auto y_interval = (axis1.max - axis1.min) / axis1.n_bins;

    auto gid = threadIdx.x + blockIdx.x * blockDim.x;

    pt = test::point3{axis0.min + gid * x_interval,
                      axis1.min + gid * y_interval, 0.5};
}


__global__ void mem_test_kernel(vecmem::data::vector_view<int> vec_view);
    
void mem_test(vecmem::vector<int>& vec){
    vecmem::data::vector_view<int> vec_view = vecmem::get_data(vec);

    int num_blocks = 1;
    int num_threads = 1;

    // run the kernel
    mem_test_kernel<<<num_blocks, num_threads>>>(vec_view);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());    
}


__global__ void mem_test_kernel(vecmem::data::vector_view<int> vec_view){

    vecmem::device_vector<int> dev_vec(vec_view);

    printf("%f \n", dev_vec[0]);
}
    
}  // namespace detray
