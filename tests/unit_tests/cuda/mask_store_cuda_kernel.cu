/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "definitions/cuda_defs.hpp"
#include "mask_store_cuda_kernel.hpp"

namespace detray {

__global__ void mask_test_kernel(
    mask_store_data<rectangle, trapezoid, ring, cylinder, single, annulus>
        store_data) {

    mask_store<thrust::tuple, vecmem::device_vector, rectangle, trapezoid, ring,
               cylinder, single, annulus>
        store(store_data);

    // const int bid = blockIdx.x;
    // const int tid = threadIdx.x;

    const auto& masks = store.group<0>();
    const auto& mask = masks[0];

    const auto& values = mask.values();

    printf("%f %f \n", values[0], values[1]);
}

void mask_test(mask_store_data<rectangle, trapezoid, ring, cylinder, single,
                               annulus>& store_data) {

    /// block dim = number of groups
    int block_dim = thrust::tuple_size<decltype(store_data._data)>::value;

    /// thread dim = max({number of masks})
    std::vector<size_t> n_masks(0);
    n_masks.push_back(store_data.group<0>().size());
    n_masks.push_back(store_data.group<1>().size());
    n_masks.push_back(store_data.group<2>().size());
    n_masks.push_back(store_data.group<3>().size());
    n_masks.push_back(store_data.group<4>().size());
    n_masks.push_back(store_data.group<5>().size());

    int thread_dim = *max_element(n_masks.begin(), n_masks.end());

    mask_test_kernel<<<block_dim, thread_dim>>>(store_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
