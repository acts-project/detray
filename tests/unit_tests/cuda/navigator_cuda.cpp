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
#include "vecmem/utils/cuda/copy.hpp"

TEST(navigator_cuda, navigator) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // vecmem managed memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    using namespace detray;

    /** Tolerance for tests */
    constexpr double tol = 0.01;

    // create toy geometry
    std::size_t n_brl_layers = 4;
    std::size_t n_edc_layers = 3;
    detector_host_t det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                                   n_edc_layers);

    // create navigator
    navigator_host_t n(det);

    // navigator state
    navigator_host_t::state state(mng_mr);

    // create navigator state buffer
    navigator_host_t::state_buffer state_buffer(
        det.get_n_max_objects_per_volume(), mng_mr);
    copy.setup(state_buffer._candidates_buffer);

    auto n_data = get_data(n);

    // run navigator test
    navigator_test(n_data, state_buffer);

    copy(state_buffer._candidates_buffer, state.candidates());
}