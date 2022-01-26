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
    detector_view<detector_host_t> det_data,
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
    for (unsigned int i = 0; i < trfs.size(typename detector_host_t::context());
         i++) {
        transforms_device.data()[i] = trfs.data()[i];
    }

    // copy objects - masks
    auto& masks = det_device.masks();
    auto& rectangles = masks.template group<detector_host_t::e_rectangle2>();
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles_device[i] = rectangles[i];
    }

    auto& discs = masks.template group<detector_host_t::e_portal_ring2>();
    for (unsigned int i = 0; i < discs.size(); i++) {
        discs_device[i] = discs[i];
    }

    auto& cylinders =
        masks.template group<detector_host_t::e_portal_cylinder3>();
    for (unsigned int i = 0; i < cylinders.size(); i++) {
        cylinders_device[i] = cylinders[i];
    }

    // print output test for surface finder
    auto& sf_finder_device = det_device.get_surfaces_finder();
    for (unsigned int i_s = 0; i_s < sf_finder_device.effective_size(); i_s++) {
        auto& grid = sf_finder_device[i_s];
        for (unsigned int i = 0; i < grid.axis_p0().bins(); i++) {
            for (unsigned int j = 0; j < grid.axis_p1().bins(); j++) {
                const auto& bin = grid.bin(i, j);
                for (auto& id : bin) {
                    // printf("%d \n", id);
                }
            }
        }
    }
}

/// implementation of the test function for detector
void detector_test(
    detector_view<detector_host_t> det_data,
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

// cuda kernel to test enumeration
__global__ void enumerate_test_kernel(
    detector_view<detector_host_t> det_data,
    vecmem::data::jagged_vector_view<surface_t> surfaces_data) {

    // Get detector in device
    detector_device_t detector(det_data);

    // Get surface vector per volume
    vecmem::device_vector<surface_t> surfaces_device(
        surfaces_data.m_ptr[threadIdx.x]);

    // Get volume
    auto& vol = detector.volume_by_index(threadIdx.x);

    // Push_back surfaces to the surface vector
    for (const auto [obj_idx, obj] : enumerate(detector.surfaces(), vol)) {
        surfaces_device.push_back(obj);
    }
}

// implementation of a test function for surface enumeration
void enumerate_test(detector_view<detector_host_t> det_data,
                    vecmem::data::jagged_vector_view<surface_t> surfaces_data) {

    constexpr int block_dim = 1;

    // number of threads = number of volumes
    int thread_dim = surfaces_data.m_size;

    // run the test kernel
    enumerate_test_kernel<<<block_dim, thread_dim>>>(det_data, surfaces_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray