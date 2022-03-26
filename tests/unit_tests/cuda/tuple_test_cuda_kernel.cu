/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "tuple_test_cuda_kernel.hpp"

namespace detray {

template <template <typename...> class tuple_type>
__global__ void tuple_copy_kernel(
    vec_tuple_data<vec_tuple<tuple_type, dvector, int, float, double>> data,
    vecmem::data::vector_view<int> output1_data,
    vecmem::data::vector_view<float> output2_data,
    vecmem::data::vector_view<double> output3_data) {

    // Get vec_tuple with vecmem::device_vector
    vec_tuple<thrust::tuple, vecmem::device_vector, int, float, double> input(
        data);
    vecmem::device_vector<int> output1(output1_data);
    vecmem::device_vector<float> output2(output2_data);
    vecmem::device_vector<double> output3(output3_data);

    const auto& input1 = detail::get<0>(input._tuple);
    const auto& input2 = detail::get<1>(input._tuple);
    const auto& input3 = detail::get<2>(input._tuple);

    for (unsigned int i = 0; i < input1.size(); i++) {
        output1[i] = input1[i];
    }

    for (unsigned int i = 0; i < input2.size(); i++) {
        output2[i] = input2[i];
    }

    for (unsigned int i = 0; i < input3.size(); i++) {
        output3[i] = input3[i];
    }
}

template void tuple_copy<thrust::tuple>(
    vec_tuple_data<vec_tuple<thrust::tuple, dvector, int, float, double>>& data,
    vecmem::data::vector_view<int>& output1_data,
    vecmem::data::vector_view<float>& output2_data,
    vecmem::data::vector_view<double>& output3_data);

/* currently does not work
template void tuple_copy<std::tuple>(
    vec_tuple_data<std::tuple, int, float, double>& data,
    vecmem::data::vector_view<int>& output1_data,
    vecmem::data::vector_view<float>& output2_data,
    vecmem::data::vector_view<double>& output3_data);
*/

template <template <typename...> class tuple_type>
void tuple_copy(
    vec_tuple_data<vec_tuple<tuple_type, dvector, int, float, double>>& data,
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