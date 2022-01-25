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

void sequence_range(const darray<dindex, 2> range,
                    vecmem::data::vector_view<dindex>& check_data) {

    // run the kernel
    sequence_range_kernel<<<1, 1>>>(range, check_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

__global__ void enumerate_sequence_kernel(
    vecmem::data::vector_view<dindex> idx_data,
    vecmem::data::vector_view<unsigned int> uint_data,
    vecmem::data::vector_view<uint_holder> seq_data) {

    vecmem::device_vector<dindex> idx_vec(idx_data);
    vecmem::device_vector<unsigned int> uint_vec(uint_data);
    vecmem::device_vector<uint_holder> seq(seq_data);

    for (auto [i, v] : enumerate(seq)) {
        idx_vec.push_back(i);
        uint_vec.push_back(v.ui);
    }
}

void enumerate_sequence(vecmem::data::vector_view<dindex>& idx_data,
                        vecmem::data::vector_view<unsigned int>& uint_data,
                        vecmem::data::vector_view<uint_holder>& seq_data) {

    // run the kernel
    enumerate_sequence_kernel<<<1, 1>>>(idx_data, uint_data, seq_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

__global__ void iterate_range_kernel(vecmem::data::vector_view<int> check_data,
                                     vecmem::data::vector_view<int> seq_data,
                                     const size_t begin, const size_t end) {

    vecmem::device_vector<int> check(check_data);
    vecmem::device_vector<int> seq(seq_data);

    for (const auto& v :
         iterator_range(seq, std::array<size_t, 2>{begin, end})) {
        check.push_back(v);
    }
}

void iterate_range(vecmem::data::vector_view<int>& check_data,
                   vecmem::data::vector_view<int>& seq_data,
                   const size_t& begin, const size_t& end) {

    // run the kernel
    iterate_range_kernel<<<1, 1>>>(check_data, seq_data, begin, end);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
