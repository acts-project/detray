/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "container_cuda_kernel.hpp"
#include "detray/definitions/cuda_definitions.hpp"

namespace detray {

__global__ void get_sum_kernel(tuple_vector_container_data_type container_data,
                               vecmem::data::vector_view<double> sum_data) {

    tuple_vector_container<thrust::tuple, vecmem::device_vector, std::size_t,
                           int, float, double>
        container(container_data);
    vecmem::device_vector<double> sum(sum_data);

    const auto& g0 = container.group<0>();
    const auto& g1 = container.group<1>();
    const auto& g2 = container.group<2>();

    for (auto e : g0) {
        sum[0] += e;
    }
    for (auto e : g1) {
        sum[0] += e;
    }
    for (auto e : g2) {
        sum[0] += e;
    }
}

void get_sum(tuple_vector_container_data_type& container_data,
             vecmem::data::vector_view<double>& sum_data) {

    // run the test kernel
    get_sum_kernel<<<1, 1>>>(container_data, sum_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray