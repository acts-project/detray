/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/containers/device_vector.hpp>

#include "definitions/cuda_defs.hpp"
#include "index_geometry_cuda_kernel.hpp"

namespace detray {

__global__ void index_geometry_test_kernel(
    index_geometry_data<geometry> geometry_data) {}

void index_geometry_test(index_geometry_data<geometry>& geometry_data) {

    int block_dim = 1;
    int thread_dim = 1;

    // run the kernel
    index_geometry_test_kernel<<<block_dim, thread_dim>>>(geometry_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
