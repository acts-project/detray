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

__global__ void test_single_store_kernel(
    single_store_t::view_type store_view,
    vecmem::data::vector_view<double> sum_data) {

    single_store_dev_t store(store_view);
    vecmem::device_vector<double> sum(sum_data);

    for (const double elem : store) {
        sum[0] += elem;
    }
}

__global__ void test_tuple_cont_kernel(
    typename tuple_cont_t::view_type store_view,
    vecmem::data::vector_view<double> sum_data) {

    tuple_cont_dev_t store(store_view);
    vecmem::device_vector<double> sum(sum_data);

    const auto& g0 = store.get<0>();
    const auto& g1 = store.get<1>();

    for (auto e : g0) {
        sum[0] += e;
    }
    for (auto e : g1) {
        sum[0] += e;
    }
}

void test_single_store(single_store_t::view_type store_view,
                       vecmem::data::vector_view<double> sum_data) {

    // run the test kernel
    test_single_store_kernel<<<1, 1>>>(store_view, sum_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

void test_tuple_container(typename tuple_cont_t::view_type store_view,
                          vecmem::data::vector_view<double> sum_data) {

    // run the test kernel
    test_tuple_cont_kernel<<<1, 1>>>(store_view, sum_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/*void test_multi_store(typename multi_store_t::view_type store_view,
             vecmem::data::vector_view<double> sum_data) {

    // run the test kernel
    get_sum_kernel<<<1, 1>>>(store_view, sum_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

void test_reg_multi_store(typename reg_multi_store_t::view_type store_view,
             vecmem::data::vector_view<double> sum_data) {

    // run the test kernel
    get_sum_kernel<<<1, 1>>>(store_view, sum_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}*/

}  // namespace detray
