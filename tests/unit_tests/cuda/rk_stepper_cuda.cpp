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
    rk_stepper_type rk(B);

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
            vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

            // intialize a track
            free_track_parameters traj(ori, 0, dir, -1);

            tracks_host.push_back(traj);
            tracks_device.push_back(traj);
        }
    }

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {

        auto& traj = tracks_host[i];

        // Forward direction
        rk_stepper_type::state forward_state(traj);

        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk.step(forward_state, path_limit);
        }

        // Backward direction
        traj.flip();
        rk_stepper_type::state backward_state(traj);
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk.step(backward_state, path_limit);
        }

        path_lengths.push_back(forward_state._path_accumulated +
                               backward_state._path_accumulated);
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