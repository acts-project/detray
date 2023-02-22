/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "tests/common/test_base/propagator_test.hpp"

namespace detray {

/// test function for propagator on the host
template<typename detector_t, typename propagator_t>
auto run_propagation_host(
    vecmem::memory_resource *mr,
    const detector_t &det, propagator_t &&p,
    const dvector<free_track_parameters<transform3>> &tracks) -> std::tuple<vecmem::jagged_vector<scalar>, vecmem::jagged_vector<vector3>, vecmem::jagged_vector<free_matrix>> {

    // Create vector for track recording
    vecmem::jagged_vector<scalar> host_path_lengths(mr);
    vecmem::jagged_vector<vector3> host_positions(mr);
    vecmem::jagged_vector<free_matrix> host_jac_transports(mr);

    for (const auto& trk : tracks) {

        // Create the propagator state
        inspector_host_t::state insp_state{*mr};
        pathlimit_aborter::state pathlimit_state{path_limit};
        parameter_transporter<transform3>::state transporter_state{};
        pointwise_material_interactor<transform3>::state interactor_state{};
        parameter_resetter<transform3>::state resetter_state{};
        auto actor_states =
            thrust::tie(insp_state, pathlimit_state, transporter_state,
                        interactor_state, resetter_state);

        typename propagator_t::state state(trk, det.get_bfield(), det);

        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            constrainted_step_size);

        state._stepping.set_tolerance(rk_tolerance);

        // Run propagation
        p.propagate(state, actor_states);

        // Record the step information
        host_path_lengths.push_back(insp_state._path_lengths);
        host_positions.push_back(insp_state._positions);
        host_jac_transports.push_back(insp_state._jac_transports);
    }

    return std::make_tuple(std::move(host_path_lengths), std::move(host_positions), std::move(host_jac_transports));
}

/// Launch the propagation test kernel
template<typename bfield_bknd_t>
extern void propagator_test(
    typename detector_host_t<bfield_bknd_t>::detector_view_type,
    vecmem::data::vector_view<free_track_parameters<transform3>>&,
    vecmem::data::jagged_vector_view<intersection_t>&,
    vecmem::data::jagged_vector_view<scalar>&,
    vecmem::data::jagged_vector_view<vector3>&,
    vecmem::data::jagged_vector_view<free_matrix>&);

/// test function for propagator on the device
/*template<typename bfield_bknd_t, typename detector_t>
auto run_propagation_device(
    vecmem::memory_resource *mr,
    detector_t &det, 
    dvector<free_track_parameters<transform3>> &tracks,
    const vecmem::jagged_vector<vector3> &host_positions) -> std::tuple<vecmem::jagged_vector<scalar>, vecmem::jagged_vector<vector3>, vecmem::jagged_vector<free_matrix>> {

    // Helper object for performing memory copies.
    vecmem::copy copy;

    // Get detector data
    auto det_data = get_data(det);

    // Get tracks data
    auto tracks_data = vecmem::get_data(tracks);

    // Create navigator candidates buffer
    auto candidates_buffer =
        create_candidates_buffer(det, theta_steps * phi_steps, *mr);
    copy.setup(candidates_buffer);

    // Create vector buffer for track recording
    std::vector<std::size_t> sizes(theta_steps * phi_steps, 0);
    std::vector<std::size_t> capacities;
    for (auto& r : host_positions) {
        capacities.push_back(r.size());
    }

    vecmem::data::jagged_vector_buffer<scalar> path_lengths_buffer(
        sizes, capacities, *mr);
    vecmem::data::jagged_vector_buffer<vector3> positions_buffer(
        sizes, capacities, *mr);
    vecmem::data::jagged_vector_buffer<free_matrix> jac_transports_buffer(
        sizes, capacities, *mr);

    copy.setup(path_lengths_buffer);
    copy.setup(positions_buffer);
    copy.setup(jac_transports_buffer);

    // Run the propagator test for GPU device
    propagator_test<bfield_bknd_t>(det_data, tracks_data, candidates_buffer,
                    path_lengths_buffer, positions_buffer,
                    jac_transports_buffer);

    vecmem::jagged_vector<scalar> device_path_lengths(mr);
    vecmem::jagged_vector<vector3> device_positions(mr);
    vecmem::jagged_vector<free_matrix> device_jac_transports(mr);

    copy(path_lengths_buffer, device_path_lengths);
    copy(positions_buffer, device_positions);
    copy(jac_transports_buffer, device_jac_transports);

    return std::make_tuple(std::move(device_path_lengths), std::move(device_positions), std::move(device_jac_transports));
}

template
std::tuple<vecmem::jagged_vector<scalar>, vecmem::jagged_vector<vector3>, vecmem::jagged_vector<free_matrix>> run_propagation_device<const_bfield_bknd_t, detector_host_t<const_bfield_bknd_t>>(
    vecmem::memory_resource *,
    detector_host_t<const_bfield_bknd_t> &,
    dvector<free_track_parameters<transform3>> &,
    const vecmem::jagged_vector<vector3> &);*/

}  // namespace detray
