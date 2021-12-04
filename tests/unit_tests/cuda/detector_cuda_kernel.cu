/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detector_cuda_kernel.hpp"
#include "detray/definitions/cuda_defs.hpp"

namespace detray {

// cuda kernel to copy sub-detector objects
__global__ void detector_test_kernel(
    detector_data<detector_t> det_data,
    vecmem::data::vector_view<volume_t> volumes_data,
    vecmem::data::vector_view<surface_t> surfaces_data,
    static_transform_store_data<transform_store_t> transforms_data,
    vecmem::data::vector_view<rectangle_t> rectangles_data,
    vecmem::data::vector_view<disc_t> discs_data,
    vecmem::data::vector_view<cylinder_t> cylinders_data) {

    // convert toy detector_data into detector w/ device vectors
    detector_device_t det_device(det_data);

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
        volumes_device[i] = det_device.volume_by_index(i);
    }

    // copy objects - surfaces
    for (unsigned int i = 0; i < det_device.surfaces().size(); i++) {
        surfaces_device[i] = det_device.surfaces()[i];
    }

    // copy objects - transforms
    auto& trfs = det_device.transforms();
    for (unsigned int i = 0; i < trfs.size(typename detector_t::context());
         i++) {
        transforms_device.data()[i] = trfs.data()[i];
    }

    // copy objects - masks
    auto& masks = det_device.masks();
    auto& rectangles = masks.template group<detector_t::e_rectangle2>();
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles_device[i] = rectangles[i];
    }

    auto& discs = masks.template group<detector_t::e_portal_ring2>();
    for (unsigned int i = 0; i < discs.size(); i++) {
        discs_device[i] = discs[i];
    }

    auto& cylinders = masks.template group<detector_t::e_portal_cylinder3>();
    for (unsigned int i = 0; i < cylinders.size(); i++) {
        cylinders_device[i] = cylinders[i];
    }
}

/// implementation of the test function for detector
void detector_test(
    detector_data<detector_t>& det_data,
    vecmem::data::vector_view<volume_t>& volumes_data,
    vecmem::data::vector_view<surface_t>& surfaces_data,
    static_transform_store_data<transform_store_t>& transforms_data,
    vecmem::data::vector_view<rectangle_t>& rectangles_data,
    vecmem::data::vector_view<disc_t>& discs_data,
    vecmem::data::vector_view<cylinder_t>& cylinders_data) {

    constexpr int block_dim = 1;
    constexpr int thread_dim = 1;

    // run the test kernel
    detector_test_kernel<<<block_dim, thread_dim>>>(
        det_data, volumes_data, surfaces_data, transforms_data, rectangles_data,
        discs_data, cylinders_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray