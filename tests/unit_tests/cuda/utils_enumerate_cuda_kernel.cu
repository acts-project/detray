/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_defs.hpp"
#include "utils_enumerate_cuda_kernel.hpp"

namespace detray {

__global__ void sequence_single_kernel(
    vecmem::data::vector_view<dindex> check_data,
    vecmem::data::vector_view<dindex> single_data) {

    vecmem::device_vector<dindex> check(check_data);
    vecmem::device_vector<dindex> single(single_data);

    for (auto i : sequence(single[0])) {
        check[0] += i;
    }
}

// test function for enumeration with single integer
void sequence_single(vecmem::data::vector_view<dindex>& check_data,
                     vecmem::data::vector_view<dindex>& single_data) {

    // run the kernel
    sequence_single_kernel<<<1, 1>>>(check_data, single_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

__global__ void sequence_range_kernel(
    const darray<dindex, 2> range,
    vecmem::data::vector_view<dindex> check_data) {

    vecmem::device_vector<dindex> check(check_data);

    for (auto i : sequence(range)) {
        check.push_back(i);
    }
}

// test function for enumeration with range
void sequence_range(const darray<dindex, 2> range,
                    vecmem::data::vector_view<dindex>& check_data) {

    // run the kernel
    sequence_range_kernel<<<1, 1>>>(range, check_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
