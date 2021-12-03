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

    // convert toy detector_data into detector w/ device vectors
    detector<darray, thrust::tuple, vecmem::device_vector,
             vecmem::jagged_device_vector>
        det_device(det_data);

    // convert subdetector data objects into objects w/ device vectors
    vecmem::device_vector<volume_t> volumes_device(volumes_data);
    vecmem::device_vector<surface_t> surfaces_device(surfaces_data);
    static_transform_store<vecmem::device_vector> transforms_device(
        transforms_data);
    vecmem::device_vector<rectangle_t> rectangles_device(rectangles_data);
    vecmem::device_vector<disc_t> discs_device(discs_data);
    vecmem::device_vector<cylinder_t> cylinders_device(cylinders_data);

    // copy objects - volume
    for (unsigned int i = 0; i < det_device.volumes().size(); i++) {
        // printf("%d \n", i);
        // auto vol = det_device.volumes_test()[i];
        // printf("%d \n", vol.index());

        // volumes_device[i] = det_device.volumes()[i];
        // surfaces_device[i] = det_device.surfaces()[i];
    }

    // copy objects - mask
    auto& masks = det_device.masks();
    auto& rectangles = masks.template group<detector_t::e_rectangle2>();
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles_device[i] = rectangles[i];
        // rectangles_device[i][0] = rectangles[i][0];
        // rectangles_device[i][1] = rectangles[i][1];
    }
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