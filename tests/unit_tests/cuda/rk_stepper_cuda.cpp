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

        auto &traj = tracks_host[i];
        free_track_parameters c_traj(traj);

        // RK Stepping into forward direction
        prop_state<rk_stepper_t::state, nav_state> propagation{
            rk_stepper_t::state{traj}, nav_state{}};
        prop_state<crk_stepper_t::state, nav_state> c_propagation{
            crk_stepper_t::state{c_traj}, nav_state{}};

        rk_stepper_t::state &rk_state = propagation._stepping;
        crk_stepper_t::state &crk_state = c_propagation._stepping;

        nav_state &n_state = propagation._navigation;
        nav_state &cn_state = c_propagation._navigation;

        crk_state.template set_constraint<step::constraint::e_user>(
            0.5 * unit_constants::mm);
        n_state._step_size = 1. * unit_constants::mm;
        cn_state._step_size = 1. * unit_constants::mm;
        ASSERT_NEAR(crk_state.constraints().template size<>(),
                    0.5 * unit_constants::mm, epsilon);
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        // check constrained steps
        EXPECT_NEAR(rk_state.path_length(), crk_state.path_length(), epsilon);

        // Backward direction
        // Roll the same track back to the origin
        // Use the same path length, since there is no overstepping
        scalar path_length = rk_state.path_length();
        n_state._step_size *= -1. * unit_constants::mm;
        cn_state._step_size *= -1. * unit_constants::mm;
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
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

TEST(rk_stepper_cuda, bound_state) {

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;

    // test surface
    const vector3 u{0, 1, 0};
    const vector3 w{1, 0, 0};
    const vector3 t{0, 0, 0};
    const transform3 trf(t, w, u);

    // Generate track starting point
    vector3 local{2, 3, 0};
    vector3 mom{0.02, 0., 0.};
    scalar time = 0.;
    scalar q = -1.;

    // bound vector
    typename bound_track_parameters::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0) = local[1];
    getter::element(bound_vector, e_bound_phi, 0) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0) = time;

    // bound covariance
    typename bound_track_parameters::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.;

    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.;

    // bound track parameter
    const bound_track_parameters in_param(0, bound_vector, bound_cov);
    const vector3 B{0, 0, 1. * unit_constants::T};

    /**
     * Get CPU bound parameter after one turn
     */

    mag_field_t mag_field(B);
    prop_state<crk_stepper_t::state, nav_state> propagation{
        crk_stepper_t::state(in_param, trf), nav_state{}};
    crk_stepper_t::state &crk_state = propagation._stepping;
    nav_state &n_state = propagation._navigation;

    // Decrease tolerance down to 1e-8
    crk_state.set_tolerance(rk_tolerance);

    // RK stepper and its state
    crk_stepper_t crk_stepper(mag_field);

    // Path length per turn
    scalar S = 2. * std::fabs(1. / in_param.qop()) / getter::norm(B) * M_PI;

    // Run stepper for half turn
    unsigned int max_steps = 1e4;

    for (unsigned int i = 0; i < max_steps; i++) {

        crk_state.set_constraint(S - crk_state.path_length());

        n_state._step_size = S;

        crk_stepper.step(propagation);

        if (std::abs(S - crk_state.path_length()) < 1e-6) {
            break;
        }

        // Make sure that we didn't reach the end of for loop
        ASSERT_TRUE(i < max_steps - 1);
    }

    // Bound state after one turn propagation
    const auto out_param_cpu = crk_stepper.bound_state(propagation, trf);

    /**
     * Get CUDA bound parameter after one turn
     */
    vecmem::vector<bound_track_parameters> out_param_cuda(1, &mng_mr);

    bound_state_test(vecmem::get_data(out_param_cuda), in_param, B, trf);

    /**
     * Compare CPU and CUDA
     */
    const auto bvec_cpu = out_param_cpu.vector();
    const auto bcov_cpu = out_param_cpu.covariance();

    const auto bvec_cuda = out_param_cuda[0].vector();
    const auto bcov_cuda = out_param_cuda[0].covariance();

    for (size_type i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bvec_cpu, i, 0),
                    matrix_operator().element(bvec_cuda, i, 0), epsilon);
    }

    for (size_type i = 0; i < e_bound_size; i++) {
        for (size_type j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bcov_cpu, i, j),
                        matrix_operator().element(bcov_cuda, i, j), epsilon);
        }
    }
}