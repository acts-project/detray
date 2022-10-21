/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detector_cuda_kernel.hpp"
#include "detray/definitions/cuda_definitions.hpp"

namespace detray {

// cuda kernel to copy sub-detector objects
__global__ void detector_test_kernel(
    detector_view<detector_host_t> det_data,
    vecmem::data::vector_view<volume_t> volumes_data,
    vecmem::data::vector_view<surface_t> surfaces_data,
    vecmem::data::vector_view<transform_t> transforms_data,
    vecmem::data::vector_view<rectangle_t> rectangles_data,
    vecmem::data::vector_view<disc_t> discs_data,
    vecmem::data::vector_view<cylinder_t> cylinders_data) {

    // convert toy detector_data into detector w/ device vectors
    detector_device_t det_device(det_data);

    // convert subdetector data objects into objects w/ device vectors
    vecmem::device_vector<volume_t> volumes_device(volumes_data);
    vecmem::device_vector<surface_t> surfaces_device(surfaces_data);
    vecmem::device_vector<transform_t> transforms_device(transforms_data);
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
    auto& trfs = det_device.transform_store();
    auto ctx = typename detector_host_t::geometry_context{};
    for (unsigned int i = 0; i < trfs.size(ctx); i++) {
        transforms_device[i] = trfs.at(i, ctx);
    }

    // copy objects - masks
    auto& masks = det_device.mask_store();
    auto& rectangles =
        masks.template get<detector_host_t::masks::id::e_rectangle2>();
    for (unsigned int i = 0; i < rectangles.size(); i++) {
        rectangles_device[i] = rectangles[i];
    }

    auto& discs =
        masks.template get<detector_host_t::masks::id::e_portal_ring2>();
    for (unsigned int i = 0; i < discs.size(); i++) {
        discs_device[i] = discs[i];
    }

    auto& cylinders =
        masks.template get<detector_host_t::masks::id::e_portal_cylinder2>();
    for (unsigned int i = 0; i < cylinders.size(); i++) {
        cylinders_device[i] = cylinders[i];
    }

    // print output test for surface finder
    /*auto& sf_finder_device = det_device.sf_finders_store();
    for (unsigned int i_s = 0; i_s < sf_finder_device.size(); i_s++) {
        auto& grid = sf_finder_device[i_s];
        for (unsigned int i = 0; i < grid.axis_p0().bins(); i++) {
            for (unsigned int j = 0; j < grid.axis_p1().bins(); j++) {
                const auto& bin = grid.bin(i, j);
                for (auto& id : bin) {
                    // printf("%d \n", id);
                }
            }
        }
    }*/
}

/// implementation of the test function for detector
void detector_test(detector_view<detector_host_t> det_data,
                   vecmem::data::vector_view<volume_t> volumes_data,
                   vecmem::data::vector_view<surface_t> surfaces_data,
                   vecmem::data::vector_view<transform_t> transforms_data,
                   vecmem::data::vector_view<rectangle_t> rectangles_data,
                   vecmem::data::vector_view<disc_t> discs_data,
                   vecmem::data::vector_view<cylinder_t> cylinders_data) {

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

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    // Get detector in device
    detector_device_t detector(det_data);

    // Get surface vector per volume
    vecmem::jagged_device_vector<surface_t> all_surfaces(surfaces_data);

    if (gid >= all_surfaces.size()) {
        return;
    }

    vecmem::device_vector<surface_t> surfaces = all_surfaces.at(gid);

    // Get volume
    auto& vol = detector.volume_by_index(gid);

    // Push_back surfaces to the surface vector
    for (const auto [obj_idx, obj] :
         detray::views::enumerate(detector.surfaces(), vol)) {
        surfaces.push_back(obj);
    }
}

// implementation of a test function for surface enumeration
void enumerate_test(detector_view<detector_host_t> det_data,
                    vecmem::data::jagged_vector_view<surface_t> surfaces_data) {

    constexpr int thread_dim = WARP_SIZE * 2;

    int block_dim = surfaces_data.m_size / thread_dim + 1;

    // run the test kernel
    enumerate_test_kernel<<<block_dim, thread_dim>>>(det_data, surfaces_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray