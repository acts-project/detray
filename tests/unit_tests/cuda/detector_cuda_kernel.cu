/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detector_cuda_kernel.hpp"
#include "detray/definitions/cuda_defs.hpp"

namespace detray {

__global__ void detector_test_kernel(
    detector_data<detector_t> det_data,
    vecmem::data::vector_view<volume_t> volumes_data,
    vecmem::data::vector_view<surface_t> surfaces_data,
    static_transform_store_data<transform_store_t> transforms_data,
    vecmem::data::vector_view<rectangle_t> rectangles_data,
    vecmem::data::vector_view<disc_t> discs_data,
    vecmem::data::vector_view<cylinder_t> cylinders_data) {

    detector<darray, thrust::tuple, vecmem::device_vector,
             vecmem::jagged_device_vector>
        det_device(det_data);
}

void detector_test(
    detector_data<detector_t>& det_data,
    vecmem::data::vector_view<volume_t>& volumes_data,
    vecmem::data::vector_view<surface_t>& surfaces_data,
    static_transform_store_data<transform_store_t>& transforms_data,
    vecmem::data::vector_view<rectangle_t>& rectangles_data,
    vecmem::data::vector_view<disc_t>& discs_data,
    vecmem::data::vector_view<cylinder_t>& cylinders_data) {

    int block_dim = 1;
    int thread_dim = 1;

    // run the test kernel
    detector_test_kernel<<<block_dim, thread_dim>>>(
        det_data, volumes_data, surfaces_data, transforms_data, rectangles_data,
        discs_data, cylinders_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray