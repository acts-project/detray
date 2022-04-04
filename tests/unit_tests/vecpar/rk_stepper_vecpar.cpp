/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "rk_stepper_vecpar.hpp"

#include <gtest/gtest.h>
#include <chrono>

#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

#include "vecpar/all/main.hpp"

TEST(rk_stepper_algo_vecpar, rk_stepper) {

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
            rk.step(forward_state, n_state);
        }

        // Backward direction
        traj.flip();
        rk_stepper_type::state backward_state(traj);
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk.step(backward_state, n_state);
        }

        path_lengths.push_back(2 * path_limit -
                               forward_state.dist_to_path_limit() -
                               backward_state.dist_to_path_limit());
    }

    // Run RK stepper
    rk_stepper_algorithm rk_stepper_algo;
    vecpar::parallel_map(rk_stepper_algo, mng_mr, tracks_device, B);

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {
        auto host_pos = tracks_host[i].pos();
        auto device_pos = tracks_device[i].pos();

        auto host_relative_error = 1. / path_lengths[i] * (host_pos - ori);
        auto device_relative_error = 1. / path_lengths[i] * (device_pos - ori);

        EXPECT_NEAR(getter::norm(host_relative_error), 0, epsilon);
        EXPECT_NEAR(getter::norm(device_relative_error), 0, epsilon);
    }
}

/*
TEST(rk_stepper_algo_vecpar, rk_stepper_timed) {

    // VecMem memory resource(s)
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;

    // Create the vector of initial track parameters
    vecmem::vector<free_track_parameters> tracks_host(&mng_mr);
    vecmem::vector<free_track_parameters> tracks_device(&mng_mr);
    vecmem::vector<free_track_parameters> tracks_benchmark(&mng_mr);

    // Create the vector of accumulated path lengths
    vecmem::vector<scalar> path_lengths(&host_mr);

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};

    // Set the magnetic field
    const vector3 B{0, 0, 2 * unit_constants::T};

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
            tracks_benchmark.push_back(traj);
        }
    }

    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time;

#if defined(_OPENMP)
    start_time = std::chrono::high_resolution_clock::now();

    // Define RK stepper
    rk_stepper_type rk(B);
    nav_state n_state{};

    #pragma omp parallel for
    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {

        auto& traj = tracks_host[i];

        // Forward direction
        rk_stepper_type::state forward_state(traj);

        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk.step(forward_state, n_state);
        }

        // Backward direction
        traj.flip();
        rk_stepper_type::state backward_state(traj);
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk.step(backward_state, n_state);
        }
    }

    end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_cpu = end_time - start_time;
    printf("CPU OMP time  = %f s\n", time_cpu.count());
#else
    start_time = std::chrono::high_resolution_clock::now();

   // Define RK stepper
   rk_stepper_type rk(B);
   nav_state n_state{};

   for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {

       auto& traj = tracks_host[i];

       // Forward direction
       rk_stepper_type::state forward_state(traj);

       for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
           rk.step(forward_state, n_state);
       }

       // Backward direction
       traj.flip();
       rk_stepper_type::state backward_state(traj);
       for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
           rk.step(backward_state, n_state);
       }
   }

   end_time = std::chrono::high_resolution_clock::now();

   std::chrono::duration<double> time_cpu = end_time - start_time;
   printf("CPU seq time  = %f s\n", time_cpu.count());
#endif

    /// get data for test bench

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {

        auto& traj = tracks_benchmark[i];

        // Forward direction
        rk_stepper_type::state forward_state(traj);

        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk.step(forward_state, n_state);
        }

        // Backward direction
        traj.flip();
        rk_stepper_type::state backward_state(traj);
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk.step(backward_state, n_state);
        }

        path_lengths.push_back(2 * path_limit -
                               forward_state.dist_to_path_limit() -
                               backward_state.dist_to_path_limit());

    }

    // Run RK stepper in parallel on CPU/GPU
    start_time = std::chrono::high_resolution_clock::now();
    rk_stepper_algorithm rk_stepper_algo;
    vecpar::parallel_map(rk_stepper_algo, mng_mr, tracks_device, B);
    end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_par = end_time - start_time;
    printf("CPU/GPU_vecpar_clang time  = %f s\n", time_par.count());

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {
        auto host_pos = tracks_host[i].pos();
        auto device_pos = tracks_device[i].pos();

        auto host_relative_error = 1. / path_lengths[i] * (host_pos - ori);
        auto device_relative_error = 1. / path_lengths[i] * (device_pos - ori);

        EXPECT_NEAR(getter::norm(host_relative_error), 0, epsilon);
        EXPECT_NEAR(getter::norm(device_relative_error), 0, epsilon);
    }
}
*/