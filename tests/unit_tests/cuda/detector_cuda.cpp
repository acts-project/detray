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
    detector_host_t toy_det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr);

    auto ctx0 = typename detector_host_t::context();

    // host objects
    auto& volumes_host = toy_det.volumes();
    auto& surfaces_host = toy_det.surfaces();
    auto& transforms_host = toy_det.transforms();
    auto& masks_host = toy_det.masks();
    auto& discs_host = masks_host.group<detector_host_t::e_portal_ring2>();
    auto& cylinders_host =
        masks_host.group<detector_host_t::e_portal_cylinder3>();
    auto& rectangles_host = masks_host.group<detector_host_t::e_rectangle2>();

    // copied outpus from device side
    vecmem::vector<volume_t> volumes_device(volumes_host.size(), &mng_mr);
    vecmem::vector<surface_t> surfaces_device(surfaces_host.size(), &mng_mr);
    transform_store_t transforms_device(mng_mr);
    auto& trfs = transforms_device.data();
    trfs.resize(transforms_host.size(ctx0));
    vecmem::vector<rectangle_t> rectangles_device(rectangles_host.size(),
                                                  &mng_mr);
    vecmem::vector<disc_t> discs_device(discs_host.size(), &mng_mr);
    vecmem::vector<cylinder_t> cylinders_device(cylinders_host.size(), &mng_mr);

    // get data object for toy detector
    auto toy_det_data = get_data(toy_det);

    // get data object for device outputs
    auto volumes_data = vecmem::get_data(volumes_device);
    auto surfaces_data = vecmem::get_data(surfaces_device);
    auto transforms_data = get_data(transforms_device);
    auto rectangles_data = vecmem::get_data(rectangles_device);
    auto discs_data = vecmem::get_data(discs_device);
    auto cylinders_data = vecmem::get_data(cylinders_device);

    // run the test code to copy the objects
    detector_test(toy_det_data, volumes_data, surfaces_data, transforms_data,
                  rectangles_data, discs_data, cylinders_data);

    // check if the same volume objects are copied
    for (unsigned int i = 0; i < volumes_host.size(); i++) {
        EXPECT_EQ(volumes_host[i] == volumes_device[i], true);
    }

    // check if the same surface objects are copied
    for (unsigned int i = 0; i < surfaces_host.size(); i++) {
        EXPECT_EQ(surfaces_host[i] == surfaces_device[i], true);
    }

    // check if the same transform objects are copied
    for (unsigned int i = 0; i < transforms_host.size(ctx0); i++) {
        EXPECT_EQ(transforms_host.contextual_transform(ctx0, i) ==
                      transforms_device.contextual_transform(ctx0, i),
                  true);
    }

    // check if the same masks are copied
    for (unsigned int i = 0; i < rectangles_host.size(); i++) {
        EXPECT_EQ(rectangles_host[i] == rectangles_device[i], true);
    }

    for (unsigned int i = 0; i < discs_host.size(); i++) {
        EXPECT_EQ(discs_host[i] == discs_device[i], true);
    }

    for (unsigned int i = 0; i < cylinders_host.size(); i++) {
        EXPECT_EQ(cylinders_host[i] == cylinders_device[i], true);
    }
}

// test volume enumeration
TEST(detector_cuda, enumerate) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // create toy geometry
    detector_host_t toy_det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr);

    auto ctx0 = typename detector_host_t::context();

    // get data object for toy detector
    auto toy_det_data = get_data(toy_det);

    // run the test code to test enumerate
    enumerate_test(toy_det_data);
}