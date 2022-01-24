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

    detector_host_t det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                                   n_edc_layers);

    // create navigator
    navigator_host_t n(det);

    // navigator state
    navigator_host_t::state state(mng_mr);

    // create navigator candidates buffer
    vecmem::data::vector_buffer<intersection> candidates_buffer(
        det.get_n_max_objects_per_volume(), 0, mng_mr);

    copy.setup(candidates_buffer);

    auto n_data = get_data(n);

    // run navigator test
    navigator_test(n_data, candidates_buffer);

    // copy candidates buffer into state candidates
    copy(candidates_buffer, state.candidates());
}

TEST(geometry_navigation_cuda, navigator) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // vecmem managed memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    detector_host_t det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                                   n_edc_layers);

    // create navigator
    navigator_host_t n(det);

    // create the vector of initial track parameters
    vecmem::vector<track<nav_context>> tracks(&mng_mr);

    const point3 ori{0., 0., 0.};

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.001 + itheta * (M_PI - 0.001) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            const point3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                             cos_theta};

            track<nav_context> ray = {.pos = ori, .dir = dir};

            tracks.push_back(ray);
        }
    }

    // navigator data
    auto n_data = get_data(n);

    // create navigator candidates buffer
    const auto n_max_objects = det.get_n_max_objects_per_volume();
    vecmem::data::jagged_vector_buffer<intersection> candidates_buffer(
        std::vector<std::size_t>(theta_steps * phi_steps, 0),
        std::vector<std::size_t>(theta_steps * phi_steps, n_max_objects),
        mng_mr);
    copy.setup(candidates_buffer);

    // track data
    auto tracks_data = vecmem::get_data(tracks);

    // create intersection buffer
    const auto n_volumes = det.volumes().size();
    vecmem::data::jagged_vector_buffer<std::pair<dindex, intersection>>
        intersection_record_buffer(
            std::vector<std::size_t>(theta_steps * phi_steps, 0),
            std::vector<std::size_t>(theta_steps * phi_steps, n_volumes),
            mng_mr);
    copy.setup(intersection_record_buffer);

    // run navigator test
    geometry_navigation_test(n_data, candidates_buffer, tracks_data,
                             intersection_record_buffer);

    // copy candidates buffer into state candidates
    // copy(candidates_buffer, state.candidates());
}