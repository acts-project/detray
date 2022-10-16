/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "detray/utils/ranges.hpp"
#include "transform_store_cuda_kernel.hpp"

namespace detray {

__global__ void transform_test_kernel(
    vecmem::data::vector_view<point3<detray::scalar>> input_data,
    typename host_transform_store_t::view_type store_data,
    vecmem::data::vector_view<point3<detray::scalar>> output_data) {

    typename device_transform_store_t::context_type ctx0;
    device_transform_store_t store(store_data);

    vecmem::device_vector<point3<detray::scalar>> input(input_data);
    vecmem::device_vector<point3<detray::scalar>> output(output_data);

    auto range = detray::ranges::subrange(store.get(ctx0),
                                          dindex_range{0, store.size(ctx0)});
    output[threadIdx.x] =
        range[threadIdx.x].point_to_global(input[threadIdx.x]);
}

void transform_test(
    vecmem::data::vector_view<point3<detray::scalar>> input_data,
    typename host_transform_store_t::view_type store_data,
    vecmem::data::vector_view<point3<detray::scalar>> output_data,
    std::size_t n_transforms) {

    int block_dim = 1;
    int thread_dim(n_transforms);

    // run the kernel
    transform_test_kernel<<<block_dim, thread_dim>>>(input_data, store_data,
                                                     output_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}
}  // namespace detray
