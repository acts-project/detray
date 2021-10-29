/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/containers/device_vector.hpp>

#include "definitions/cuda_defs.hpp"
#include "transform_store_cuda_kernel.hpp"

namespace detray {

__global__ void transform_test_kernel(
    vecmem::data::vector_view<point3> input_data,
    static_transform_store_data store_data,
    vecmem::data::vector_view<point3> output_data) {
    static_transform_store<vecmem::device_vector>::context ctx0;
    static_transform_store<vecmem::device_vector> store(store_data);

    vecmem::device_vector<point3> input(input_data);
    vecmem::device_vector<point3> output(output_data);

    auto range = store.range(0, store.size(ctx0), ctx0);

    output[threadIdx.x] =
        range[threadIdx.x].point_to_global(input[threadIdx.x]);
}

void transform_test(vecmem::data::vector_view<point3>& input_data,
                    static_transform_store_data& store_data,
                    vecmem::data::vector_view<point3>& output_data) {

    int block_dim = 1;
    int thread_dim(store_data._data.size());

    // run the kernel
    transform_test_kernel<<<block_dim, thread_dim>>>(input_data, store_data,
                                                     output_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}
}  // namespace detray
