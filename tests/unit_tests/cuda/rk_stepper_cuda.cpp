/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#include <gtest/gtest.h>

#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

#include "rk_stepper_cuda_kernel.hpp"

TEST(rk_stepper_cuda, rk_stepper) {

    // VecMem memory resource(s)
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;

    // Create the vector of initial track parameters
    vecmem::vector<free_track_parameters> tracks_host(&mng_mr);
    vecmem::vector<free_track_parameters> tracks_device(&mng_mr);

    // Create the vector of accumulated path lengths
    vecmem::vector<scalar> path_lengths(&host_mr);

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};

    // Set the magnetic field
    const vector3 B{0, 0, 2 * unit_constants::T};

    // Define RK stepper
    rk_stepper_t rk_stepper(B);
    crk_stepper_t crk_stepper(B);
    nav_state n_state{};

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
            const vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                              cos_theta};

            // intialize a track
            free_track_parameters traj(ori, 0, dir, -1);

            tracks_host.push_back(traj);
            tracks_device.push_back(traj);
        }
    }

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {

        auto& traj = tracks_host[i];
        free_track_parameters c_traj(traj);

        // Forward direction
        rk_stepper_t::state rk_state(traj);
        crk_stepper_t::state crk_state(c_traj);

        crk_state.template set_constraint<step::constraint::e_user>(
            0.5 * unit_constants::mm);
        n_state._step_size = 1. * unit_constants::mm;
        ASSERT_NEAR(crk_state.constraints().template size<>(),
                    0.5 * unit_constants::mm, epsilon);
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk_stepper.step(rk_state, n_state);
            crk_stepper.step(crk_state, n_state);
            crk_stepper.step(crk_state, n_state);
        }

        // check constrained steps
        EXPECT_NEAR(rk_state.path_length(), crk_state.path_length(), epsilon);

        // Backward direction
        // Roll the same track back to the origin
        // Use the same path length, since there is no overstepping
        scalar path_length = rk_state.path_length();
        n_state._step_size *= -1. * unit_constants::mm;
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk_stepper.step(rk_state, n_state);
            crk_stepper.step(crk_state, n_state);
            crk_stepper.step(crk_state, n_state);
        }

        EXPECT_NEAR(rk_state.path_length(), crk_state.path_length(), epsilon);

        path_lengths.push_back(2 * path_length);
    }

    // Get tracks data
    auto tracks_data = vecmem::get_data(tracks_device);

    // Run RK stepper cuda kernel
    rk_stepper_test(tracks_data, B);

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {
        auto host_pos = tracks_host[i].pos();
        auto device_pos = tracks_device[i].pos();

        auto host_relative_error = 1. / path_lengths[i] * (host_pos - ori);
        auto device_relative_error = 1. / path_lengths[i] * (device_pos - ori);

        EXPECT_NEAR(getter::norm(host_relative_error), 0, epsilon);
        EXPECT_NEAR(getter::norm(device_relative_error), 0, epsilon);
    }
}