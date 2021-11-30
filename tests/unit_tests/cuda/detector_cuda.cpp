/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <iostream>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "detector_cuda_kernel.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

TEST(detector_cuda, detector) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // create toy geometry
    auto toy_det = create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                                       vecmem::jagged_vector>(mng_mr);

    // host objects
    auto volumes_host = toy_det.volumes();
    auto surfaces_host = toy_det.surfaces();
    auto transforms_host = toy_det.transforms();
    auto masks_host = toy_det.masks();
    auto& discs_host = masks_host.group<detector_t::e_portal_ring2>();
    auto& cylinders_host = masks_host.group<detector_t::e_portal_cylinder3>();
    auto& rectangles_host = masks_host.group<detector_t::e_rectangle2>();

    // copied outpus from device side
    vecmem::vector<volume_t> volumes_device(volumes_host.size(), &mng_mr);
    vecmem::vector<surface_t> surfaces_device(surfaces_host.size(), &mng_mr);
    transform_store_t transforms_device(mng_mr);
    auto& trfs = transforms_device.data();
    trfs.resize(transforms_host.size(typename detector_t::context()));
    vecmem::vector<rectangle_t> rectangles_device(discs_host.size(), &mng_mr);
    vecmem::vector<disc_t> discs_device(cylinders_host.size(), &mng_mr);
    vecmem::vector<cylinder_t> cylinders_device(rectangles_host.size(),
                                                &mng_mr);

    // get data object for toy detector
    auto toy_det_data = get_data(toy_det);

    // get data object for device outputs
    auto volumes_data = vecmem::get_data(volumes_device);
    auto surfaces_data = vecmem::get_data(surfaces_device);
    auto transforms_data = get_data(transforms_device);
    auto rectangles_data = vecmem::get_data(rectangles_device);
    auto discs_data = vecmem::get_data(discs_device);
    auto cylinders_data = vecmem::get_data(cylinders_device);

    // run the test code
    detector_test(toy_det_data, volumes_data, surfaces_data, transforms_data,
                  rectangles_data, discs_data, cylinders_data);
}