/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/test/propagator_test.hpp"
#include "queue_wrapper.hpp"

// Vecmem include(s)
#include "vecmem/utils/sycl/copy.hpp"

namespace detray {

/// Launch the propagation test kernel
template <typename bfield_bknd_t, typename detector_t>
void propagator_test(
    typename detector_t::view_type, covfie::field_view<bfield_bknd_t>,
    vecmem::data::vector_view<track_t> &,
    vecmem::data::jagged_vector_view<intersection_t<detector_t>> &,
    vecmem::data::jagged_vector_view<scalar> &,
    vecmem::data::jagged_vector_view<vector3_t> &,
    vecmem::data::jagged_vector_view<free_matrix> &, sycl::queue_wrapper);

/// test function for propagator on the device
template <typename bfield_bknd_t, typename detector_t>
inline auto run_propagation_device(
    vecmem::memory_resource *mr, detector_t &det,
    typename detector_t::view_type det_view,
    covfie::field_view<bfield_bknd_t> field_data, sycl::queue_wrapper queue,
    dvector<track_t> &tracks,
    const vecmem::jagged_vector<vector3_t> &host_positions)
    -> std::tuple<vecmem::jagged_vector<scalar>,
                  vecmem::jagged_vector<vector3_t>,
                  vecmem::jagged_vector<free_matrix>> {

    // Helper object for performing memory copies.
    vecmem::copy copy;

    // Get tracks data
    auto tracks_data = vecmem::get_data(tracks);

    // Create navigator candidates buffer
    auto candidates_buffer =
        create_candidates_buffer(det, theta_steps * phi_steps, *mr);
    copy.setup(candidates_buffer);

    // Create vector buffer for track recording
    std::vector<std::size_t> sizes(theta_steps * phi_steps, 0);
    std::vector<std::size_t> capacities;
    for (auto &r : host_positions) {
        capacities.push_back(r.size());
    }

    vecmem::data::jagged_vector_buffer<scalar> path_lengths_buffer(
        capacities, *mr, nullptr, vecmem::data::buffer_type::resizable);
    vecmem::data::jagged_vector_buffer<vector3_t> positions_buffer(
        capacities, *mr, nullptr, vecmem::data::buffer_type::resizable);
    vecmem::data::jagged_vector_buffer<free_matrix> jac_transports_buffer(
        capacities, *mr, nullptr, vecmem::data::buffer_type::resizable);

    copy.setup(path_lengths_buffer);
    copy.setup(positions_buffer);
    copy.setup(jac_transports_buffer);

    // Run the propagator test for GPU device
    propagator_test<bfield_bknd_t, detector_t>(
        det_view, field_data, tracks_data, candidates_buffer,
        path_lengths_buffer, positions_buffer, jac_transports_buffer, queue);

    vecmem::jagged_vector<scalar> device_path_lengths(mr);
    vecmem::jagged_vector<vector3_t> device_positions(mr);
    vecmem::jagged_vector<free_matrix> device_jac_transports(mr);

    copy(path_lengths_buffer, device_path_lengths);
    copy(positions_buffer, device_positions);
    copy(jac_transports_buffer, device_jac_transports);

    return std::make_tuple(std::move(device_path_lengths),
                           std::move(device_positions),
                           std::move(device_jac_transports));
}

/// Test chain for the propagator
template <typename device_bfield_bknd_t, typename host_bfield_bknd_t,
          typename detector_t>
inline auto run_propagation_test(vecmem::memory_resource *mr, ::sycl::queue *q,
                                 detector_t &det,
                                 typename detector_t::view_type det_view,
                                 covfie::field<host_bfield_bknd_t> &field) {

    // Create the vector of initial track parameterizations
    auto tracks_host = generate_tracks(mr);
    vecmem::vector<track_t> tracks_device(tracks_host, mr);

    // Host propagation
    auto &&[host_path_lengths, host_positions, host_jac_transports] =
        run_propagation_host(mr, det, field, tracks_host);

    // Device propagation
    detray::sycl::queue_wrapper queue(q);

    covfie::field<device_bfield_bknd_t> device_field(field);
    auto &&[device_path_lengths, device_positions, device_jac_transports] =
        run_propagation_device<device_bfield_bknd_t>(
            mr, det, det_view, device_field, queue, tracks_device,
            host_positions);

    // Check the results
    compare_propagation_results(host_positions, device_positions,
                                host_path_lengths, device_path_lengths,
                                host_jac_transports, device_jac_transports);
}

}  // namespace detray
