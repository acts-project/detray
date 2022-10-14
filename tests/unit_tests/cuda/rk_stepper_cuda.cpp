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
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0) = local[1];
    getter::element(bound_vector, e_bound_phi, 0) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.;

    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.;

    // bound track parameter
    const bound_track_parameters<transform3> in_param(0, bound_vector,
                                                      bound_cov);
    const vector3 B{0, 0, 1. * unit_constants::T};

    /**
     * Get CPU bound parameter after one turn
     */

    mag_field_t mag_field(
        mag_field_t::backend_t::configuration_t{B[0], B[1], B[2]});
    prop_state<crk_stepper_t::state, nav_state> propagation{
        crk_stepper_t::state(in_param, trf, mag_field), nav_state{}};
    crk_stepper_t::state &crk_state = propagation._stepping;
    nav_state &n_state = propagation._navigation;

    // Decrease tolerance down to 1e-8
    crk_state.set_tolerance(rk_tolerance);

    // RK stepper and its state
    crk_stepper_t crk_stepper;

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
    vecmem::vector<bound_track_parameters<transform3>> out_param_cuda(1,
                                                                      &mng_mr);

    bound_state_test(vecmem::get_data(out_param_cuda), in_param, B, trf);

    /**
     * Compare CPU and CUDA
     */
    const auto bvec_cpu = out_param_cpu.vector();
    const auto bcov_cpu = out_param_cpu.covariance();

    const auto bvec_cuda = out_param_cuda[0].vector();
    const auto bcov_cuda = out_param_cuda[0].covariance();

    for (std::size_t i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bvec_cpu, i, 0),
                    matrix_operator().element(bvec_cuda, i, 0), epsilon);
    }

    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bcov_cpu, i, j),
                        matrix_operator().element(bcov_cuda, i, j), epsilon);
        }
    }
}
