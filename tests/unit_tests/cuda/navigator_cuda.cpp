/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "navigator_cuda_kernel.hpp"

TEST(navigator_cuda, navigator) {

    /** Tolerance for tests */
    constexpr double tol = 0.01;

    using namespace detray;
    vecmem::cuda::managed_memory_resource mng_mr;

    std::size_t n_brl_layers = 4;
    std::size_t n_edc_layers = 3;

    detector_host_t toy_det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                                   n_edc_layers);

    using nav_context = decltype(toy_det)::context;

    navigator_host_t n(toy_det);

    auto n_data = get_data(n);

    navigator_test(n_data);
}