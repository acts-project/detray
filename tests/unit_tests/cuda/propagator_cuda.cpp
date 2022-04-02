/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "propagator_cuda_kernel.hpp"
#include "vecmem/utils/cuda/copy.hpp"

class CudaPropagatorWithRkStepper
    : public ::testing::TestWithParam<__plugin::vector3<scalar>> {};

TEST_P(CudaPropagatorWithRkStepper, propagator) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    // Create the toy geometry
    detector_host_type det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                                   n_edc_layers);

    // Create the vector of initial track parameters
    vecmem::vector<free_track_parameters> tracks_host(&mng_mr);
    vecmem::vector<free_track_parameters> tracks_device(&mng_mr);

    // Set origin position of tracks
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
            vector3 mom{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};
            mom = 10. * mom;

            // intialize a track
            free_track_parameters traj(ori, 0, mom, -1);

            // Put it into vector of trajectories
            tracks_host.push_back(traj);
            tracks_device.push_back(traj);
        }
    }

    /**
     * Host propagation
     */

    // Set the magnetic field
    vector3 B = GetParam();
    field_type B_field(B);

    // Create RK stepper
    rk_stepper_type s(B_field);

    // Create navigator
    navigator_host_type n(det);

    // Create propagator
    propagator_host_type p(std::move(s), std::move(n));

    auto actor_states =
        detail::make_tuple(track_inspector_host_type::state(mng_mr));

    // Create vector for track recording
    vecmem::jagged_vector<intersection_t> host_intersection_records(&mng_mr);

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {

        // Create the propagator state
        propagator_host_type::state<decltype(actor_states)> state(
            tracks_host[i], actor_states);

        // Run propagation
        p.propagate(state);

        // push back the intersection record
        host_intersection_records.push_back(
            detail::get<0>(actor_states)._intersections);
    }

    /**
     * Device propagation
     */

    // Get detector data
    auto det_data = get_data(det);

    // Get tracks data
    auto tracks_data = vecmem::get_data(tracks_device);

    // Create navigator candidates buffer
    auto candidates_buffer =
        create_candidates_buffer(det, theta_steps * phi_steps, dev_mr);
    copy.setup(candidates_buffer);

    // Create vector for track recording
    vecmem::jagged_vector<intersection_t> device_intersection_records(&mng_mr);

    // Create vector buffer for track recording
    std::vector<std::size_t> sizes(theta_steps * phi_steps, 0);
    std::vector<std::size_t> capacities;
    for (auto& r : host_intersection_records) {
        capacities.push_back(r.size());
    }

    vecmem::data::jagged_vector_buffer<intersection_t> intersections_buffer(
        sizes, capacities, dev_mr, &mng_mr);
    copy.setup(intersections_buffer);

    // Run the propagator test for GPU device
    propagator_test(det_data, tracks_data, candidates_buffer,
                    intersections_buffer);

    // copy back intersection record
    copy(intersections_buffer, device_intersection_records);

    for (unsigned int i = 0; i < host_intersection_records.size(); i++) {
        for (unsigned int j = 0; j < host_intersection_records[i].size(); j++) {
            auto& host_intersection = host_intersection_records[i][j];
            auto& device_intersection = device_intersection_records[i][j];

            EXPECT_EQ(host_intersection.link, device_intersection.link);
            EXPECT_NEAR(host_intersection.p3[0], device_intersection.p3[0],
                        pos_diff_tolerance);
            EXPECT_NEAR(host_intersection.p3[1], device_intersection.p3[1],
                        pos_diff_tolerance);
            EXPECT_NEAR(host_intersection.p3[2], device_intersection.p3[2],
                        pos_diff_tolerance);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation1, CudaPropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             0. * unit_constants::T, 0. * unit_constants::T,
                             2. * unit_constants::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation2, CudaPropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             0. * unit_constants::T, 1. * unit_constants::T,
                             1. * unit_constants::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation3, CudaPropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             1. * unit_constants::T, 0. * unit_constants::T,
                             1. * unit_constants::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation4, CudaPropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             1. * unit_constants::T, 1. * unit_constants::T,
                             1. * unit_constants::T}));
