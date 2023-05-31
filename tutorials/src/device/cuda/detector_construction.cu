/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detector_construction.hpp"
#include "detray/definitions/cuda_definitions.hpp"

namespace detray::tutorial {

/// Kernel that runs the entire propagation loop
__global__ void print_kernel(
    typename detray::tutorial::detector_host_t::detector_view_type det_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    if (gid > 0) {
        return;
    }

    // Setup of the device-side detector
    detray::tutorial::detector_device_t det(det_data);

    printf("Number of volumes: %d\n", det.volumes().size());
    printf("Number of transforms: %d\n", det.transform_store().size());
    printf("Number of rectangles: %d\n",
           det.mask_store().get<mask_id::e_rectangle2>().size());
    printf("Number of trapezoids: %d\n",
           det.mask_store().get<mask_id::e_trapezoid2>().size());
    printf("Number of portal discs: %d\n",
           det.mask_store().get<mask_id::e_portal_ring2>().size());
    printf("Number of portal cylinders: %d\n",
           det.mask_store().get<mask_id::e_portal_cylinder2>().size());
}

void print(
    typename detray::tutorial::detector_host_t::detector_view_type det_data) {

    // run the tutorial kernel
    print_kernel<<<1, 1>>>(det_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray::tutorial
