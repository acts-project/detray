/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/containers/device_vector.hpp>

#include "definitions/cuda_defs.hpp"
#include "grids_grid2_cuda_kernel.hpp"

namespace detray {

/*---------------------------------------------------
  test function for grid data with replace populator
  ---------------------------------------------------*/

// test1 kernel declaration
template <typename grid_data_t>
__global__ void grid_replace_test_kernel(grid_data_t grid_data);

// test1 instantiation for replace populator
template void grid_replace_test<
    grid2_data<device_replace, axis::regular<>, axis::regular<>, serializer2>>(
    grid2_data<device_replace, axis::regular<>, axis::regular<>, serializer2>&
        grid_data);

// test2 function implementation
template <typename grid2_data_t>
void grid_replace_test(grid2_data_t& grid_data) {

    const auto& axis0 = grid_data._axis_p0;
    const auto& axis1 = grid_data._axis_p1;

    int block_dim = 1;
    dim3 thread_dim(axis0.bins(), axis1.bins());

    // run the kernel
    grid_replace_test_kernel<<<block_dim, thread_dim>>>(grid_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// test2 kernel implementation
template <typename grid_data_t>
__global__ void grid_replace_test_kernel(grid_data_t grid_data) {

    using grid2_device_t =
        grid2<typename grid_data_t::populator_t,
              typename grid_data_t::axis_p0_t, typename grid_data_t::axis_p1_t,
              typename grid_data_t::serializer_t>;

    // Let's try building the grid object
    grid2_device_t g2_device(grid_data, test::point3{0, 0, 0});

    const auto& axis0 = g2_device.axis_p0();
    const auto& axis1 = g2_device.axis_p1();

    auto x_interval = (axis0.max - axis0.min) / axis0.n_bins;
    auto y_interval = (axis1.max - axis1.min) / axis1.n_bins;

    auto gid = threadIdx.x + threadIdx.y * blockDim.x;
    auto pt = test::point3{axis0.min + gid * x_interval,
                           axis1.min + gid * y_interval, 0.5};

    // replace the bin elements
    g2_device.populate(threadIdx.x, threadIdx.y, std::move(pt));
}

/*---------------------------------------------------------------
  test function for grid data with complete populator
  ---------------------------------------------------------------*/

// test2 kernel declaration
template <typename grid_data_t>
__global__ void grid_complete_kernel(grid_data_t grid_data);

// test2 instantiation for complete populator
template void grid_complete_test<
    grid2_data<device_complete, axis::regular<>, axis::regular<>, serializer2>>(
    grid2_data<device_complete, axis::regular<>, axis::regular<>, serializer2>&
        grid_data);

// test2 function implementation
template <typename grid2_data_t>
void grid_complete_test(grid2_data_t& grid_data) {

    const auto& axis0 = grid_data._axis_p0;
    const auto& axis1 = grid_data._axis_p1;

    int block_dim = 1;
    dim3 thread_dim(axis0.bins(), axis1.bins());

    // run the kernel
    grid_complete_kernel<<<block_dim, thread_dim>>>(grid_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// test2 kernel implementation
template <typename grid_data_t>
__global__ void grid_complete_kernel(grid_data_t grid_data) {

    using grid2_device_t =
        grid2<typename grid_data_t::populator_t,
              typename grid_data_t::axis_p0_t, typename grid_data_t::axis_p1_t,
              typename grid_data_t::serializer_t>;

    // Let's try building the grid object
    grid2_device_t g2_device(grid_data, test::point3{0, 0, 0});

    const auto& axis0 = g2_device.axis_p0();
    const auto& axis1 = g2_device.axis_p1();

    auto x_interval = (axis0.max - axis0.min) / axis0.n_bins;
    auto y_interval = (axis1.max - axis1.min) / axis1.n_bins;

    auto bin_id = threadIdx.x + threadIdx.y * blockDim.x;

    for (int i_p = 0; i_p < n_points; i_p++) {
        auto gid = i_p + bin_id * n_points;
        auto pt = test::point3{axis0.min + gid * x_interval,
                               axis1.min + gid * y_interval, 0.5};
        // printf("%f %f %f \n", pt[0], pt[1], pt[2]);
        g2_device.populate(threadIdx.x, threadIdx.y, std::move(pt));
    }
}

/*----------------------------------------------------
  test function for grid buffer with attach populator
  ----------------------------------------------------*/

// buffer_test kernel declaration
template <typename grid_data_t>
__global__ void grid_attach_test_kernel(grid_data_t grid_data);

// buffer_test instantiation for attach populator
template void
grid_attach_test<grid2_data<device_attach_populator<false, test::point3>,
                            axis::regular<>, axis::regular<>, serializer2>>(
    grid2_data<device_attach_populator<false, test::point3>, axis::regular<>,
               axis::regular<>, serializer2>& grid_data);

template <typename grid2_data_t>
void grid_attach_test(grid2_data_t& grid_data) {

    const auto& axis0 = grid_data._axis_p0;
    const auto& axis1 = grid_data._axis_p1;

    dim3 block_dim(axis0.bins(), axis1.bins());
    int thread_dim = 100;

    // run the kernel
    grid_attach_test_kernel<<<block_dim, thread_dim>>>(grid_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// buffer_test kernel declaration
template <typename grid_data_t>
__global__ void grid_attach_test_kernel(grid_data_t grid_data) {

    using grid2_device_t =
        grid2<typename grid_data_t::populator_t,
              typename grid_data_t::axis_p0_t, typename grid_data_t::axis_p1_t,
              typename grid_data_t::serializer_t>;

    // Let's try building the grid object
    grid2_device_t g2_device(grid_data, test::point3{0, 0, 0});

    // Fill with 100 points
    auto pt =
        test::point3{1. * threadIdx.x, 1. * threadIdx.x, 1. * threadIdx.x};
    g2_device.populate(blockIdx.x, blockIdx.y, std::move(pt));

    __syncthreads();

    if (threadIdx.x == 0 && blockIdx.x == 0 && blockIdx.y == 0) {
        auto pts = g2_device.bin(blockIdx.x, blockIdx.y);
        for (int i = 0; i < 100; i++) {
            // printf("%f %f %f \n", pts[i][0], pts[i][1], pts[i][2]);
        }
    }
}

}  // namespace detray
