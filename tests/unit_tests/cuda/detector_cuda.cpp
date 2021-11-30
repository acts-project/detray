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
    detector<darray, thrust::tuple, vecmem::vector, vecmem::jagged_vector>
        toy_det = create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                                      vecmem::jagged_vector>(mng_mr);

    // get data object for toy detector
    auto toy_det_data = get_data(toy_det);

    // run the test code
    // detector_test(toy_det_data);
}