/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_defs.hpp"
#include "tuple_test_cuda_kernel.hpp"

namespace detray {

__global__ void tuple_copy_kernel(
    vec_tuple_data<int, float, double> data,
    vecmem::data::vector_view<int> output1_data,
    vecmem::data::vector_view<float> output2_data,
    vecmem::data::vector_view<double> output3_data) {

    // Get vec_tuple with vecmem::device_vector
    vec_tuple<vecmem::device_vector, int, float, double> input(data);
    vecmem::device_vector<int> output1(output1_data);
    vecmem::device_vector<float> output2(output2_data);
    vecmem::device_vector<double> output3(output3_data);

    output1 = thrust::get<0>(input._tuple);
    output2 = thrust::get<1>(input._tuple);
    output3 = thrust::get<2>(input._tuple);
}

void tuple_copy(vec_tuple_data<int, float, double>& data,
                vecmem::data::vector_view<int>& output1_data,
                vecmem::data::vector_view<float>& output2_data,
                vecmem::data::vector_view<double>& output3_data) {

    int thread_dim = 1;
    int block_dim = 1;

    // run the test kernel
    tuple_copy_kernel<<<block_dim, thread_dim>>>(data, output1_data,
                                                 output2_data, output3_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray