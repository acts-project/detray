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

    // Create detector
    detector_host_t det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                                   n_edc_layers);

    // Create navigator
    navigator_host_t n(det, mng_mr);

    // Get navigator data
    auto n_data = get_data(n);

    // Create Navigator state
    navigator_host_t::state state(mng_mr);

    // Create navigator candidates buffer
    vecmem::data::vector_buffer<intersection> candidates_buffer(
        det.get_n_max_objects_per_volume(), 0, mng_mr);
    copy.setup(candidates_buffer);

    // Create a track
    const track<nav_context> traj{.pos = {0., 0., 0.},
                                  .dir = vector::normalize(vector3{1., 1., 0.}),
                                  .momentum = 100,
                                  .ctx = nav_context{},
                                  .overstep_tolerance = -1e-4};

    vecmem::vector<track<nav_context> > tracks(&mng_mr);
    tracks.push_back(traj);

    auto tracks_data = vecmem::get_data(tracks);

    // Run navigator test
    navigator_test(n_data, candidates_buffer, tracks_data);

    // Copy candidates buffer into state candidates
    copy(candidates_buffer, state.candidates());
}