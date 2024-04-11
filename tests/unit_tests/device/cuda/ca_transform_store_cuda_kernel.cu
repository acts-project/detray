/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "ca_transform_store_cuda_kernel.hpp"
#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

__global__ void ca_transform_test_kernel(
    vecmem::data::vector_view<detray::point3> input_data,
    typename host_ca_transform_store_t::view_type store_data,
    vecmem::data::vector_view<detray::point3> output_data,
    detray::detail::data_context context) {

    device_ca_transform_store_t store(store_data);

    vecmem::device_vector<detray::point3> input(input_data);
    vecmem::device_vector<detray::point3> output(output_data);

    auto range = detray::ranges::subrange(
        store.get().at(context.get()), dindex_range{0u, store.size(context)});
    output[threadIdx.x] = range[threadIdx.x].point_to_global(
        input[context.get() + 2 * threadIdx.x]);
}

void ca_transform_test(
    vecmem::data::vector_view<detray::point3> input_data,
    typename host_ca_transform_store_t::view_type ca_store_data,
    vecmem::data::vector_view<detray::point3> output_data,
    std::size_t n_transforms, detray::detail::data_context context) {
    int block_dim = 1;
    int thread_dim(n_transforms);

    // run the kernel
    ca_transform_test_kernel<<<block_dim, thread_dim>>>(
        input_data, ca_store_data, output_data, context);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}
}  // namespace detray
