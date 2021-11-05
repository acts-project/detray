/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_defs.hpp"
#include "tuple_test_cuda_kernel.hpp"

namespace detray {

__global__ void tuple_test_kernel(vec_tuple_data<int, float, double> data) {

    vec_tuple<vecmem::device_vector, int, float, double> input(data);

    auto a = thrust::get<0>(input._tuple);
    printf("%f \n", a[0]);
}

void tuple_test(vec_tuple_data<int, float, double>& data) {

    int thread_dim = 1;
    int block_dim = 1;

    // run the test kernel
    tuple_test_kernel<<<block_dim, thread_dim>>>(data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray