/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "material_interaction_cuda_kernel.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

TEST(material_interaction, telescope_geometry) {

    // Helper object for performing memory copies.
    vecmem::copy copy;

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;

    using ln_stepper_t = line_stepper<transform3>;

    typename ln_stepper_t::free_track_parameters_type default_trk{
        {0, 0, 0}, 0, {0, 0, 1}, -1};

    // Build from given module positions
    std::vector<scalar> positions = {0.,  10., 20., 30., 40., 50.,
                                     60., 70,  80,  90., 100.};

    // Use unbounded surfaces
    constexpr bool unbounded = true;

    // Create the telescope geometry
    detector_host_type det =
        create_telescope_detector<unbounded, ln_stepper_t, darray,
                                  thrust::tuple, vecmem::vector,
                                  vecmem::jagged_vector>(
            mng_mr, positions, ln_stepper_t(),
            typename ln_stepper_t::state{default_trk}, 20. * unit_constants::mm,
            20. * unit_constants::mm, silicon_tml<scalar>(),
            0.17 * unit_constants::cm);

    const scalar q = -1.;
    const scalar iniP = 10 * unit_constants::GeV;
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = 0.;
    getter::element(bound_vector, e_bound_loc1, 0) = 0.;
    getter::element(bound_vector, e_bound_phi, 0) = 0;
    getter::element(bound_vector, e_bound_theta, 0) = M_PI_4;
    getter::element(bound_vector, e_bound_qoverp, 0) = q / iniP;
    getter::element(bound_vector, e_bound_time, 0) = 0.;
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template identity<e_bound_size, e_bound_size>();

    // Field
    field_type B_field(B);

    std::vector<bound_track_parameters<transform3>> host_initial_states;
    std::vector<bound_track_parameters<transform3>> host_final_states;
    vecmem::vector<bound_track_parameters<transform3>> device_initial_states(
        &mng_mr);
    vecmem::vector<bound_track_parameters<transform3>> device_final_states(
        &mng_mr);

    for (int i = 0; i < n_tracks; i++) {
        // Propagator is built from the stepper and navigator
        propagator_host_type p({}, {});

        const scalar phi_step = 2 * M_PI / n_tracks;
        const scalar phi = -M_PI + phi_step * i;
        getter::element(bound_vector, e_bound_phi, 0) = phi;

        // bound track parameter
        const bound_track_parameters<transform3> bound_params(0, bound_vector,
                                                              bound_cov);

        host_initial_states.push_back(bound_params);
        device_initial_states.push_back(bound_params);
    }

    for (int i = 0; i < n_tracks; i++) {

        // Propagator is built from the stepper and navigator
        propagator_host_type p({}, {});

        // Actor states
        pathlimit_aborter::state aborter_state{};
        bound_to_bound_updater<transform3>::state bound_updater{};
        pointwise_material_interactor<transform3>::state interactor_state{};
        resetter<transform3>::state resetter_state{};

        // Create actor states tuples
        actor_chain_t::state actor_states = thrust::tie(
            aborter_state, bound_updater, interactor_state, resetter_state);

        propagator_host_type::state state(host_initial_states[i], B_field, det,
                                          actor_states);

        state._stepping().set_overstep_tolerance(overstep_tolerance);

        // Propagate the entire detector
        p.propagate(state);

        host_final_states.push_back(state._stepping._bound_params);
    }

    // Get detector data
    auto det_data = get_data(det);

    // Create navigator candidates buffer
    auto candidates_buffer = create_candidates_buffer(det, n_tracks, mng_mr);
    copy.setup(candidates_buffer);

    material_interaction_test(det_data, B, candidates_buffer,
                              vecmem::get_data(device_initial_states),
                              vecmem::get_data(device_final_states));
}