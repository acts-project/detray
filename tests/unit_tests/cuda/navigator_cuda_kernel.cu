/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_defs.hpp"
#include "navigator_cuda_kernel.hpp"

namespace detray {

__global__ void navigator_test_kernel(
    navigator_view<navigator_host_t> n_data,
    vecmem::data::vector_view<intersection> candidates_data) {

    navigator_device_t n(n_data);
}

/// test function for navigator
void navigator_test(navigator_view<navigator_host_t> n_data,
                    vecmem::data::vector_view<intersection>& candidates_data) {

    constexpr int block_dim = 1;
    constexpr int thread_dim = 1;

    // run the test kernel
    navigator_test_kernel<<<block_dim, thread_dim>>>(n_data, candidates_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray