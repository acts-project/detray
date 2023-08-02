/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/test/device/propagator_test.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Covfie include(s)
#include <covfie/cuda/backend/primitive/cuda_device_array.hpp>

namespace detray {

namespace bfield::cuda {

// Inhomogeneous field (cuda)
using inhom_bknd_t = covfie::backend::affine<covfie::backend::linear<
    covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                             covfie::backend::cuda_device_array<
                                 covfie::vector::vector_d<scalar, 3>>>>>;

}  // namespace bfield::cuda

/// Launch the propagation test kernel
template <typename bfield_bknd_t, typename detector_t>
void propagator_test(
    typename detector_t::view_type, const propagation::config &,
    covfie::field_view<bfield_bknd_t>, vecmem::data::vector_view<track_t> &,
    vecmem::data::jagged_vector_view<scalar> &,
    vecmem::data::jagged_vector_view<point3_t> &,
    vecmem::data::jagged_vector_view<free_matrix_t> &);

/// Test function for propagator on the device
template <typename bfield_bknd_t, typename detector_t>
inline auto run_propagation_device(
    vecmem::memory_resource *mr, detector_t &det,
    const propagation::config &cfg, typename detector_t::view_type det_view,
    covfie::field_view<bfield_bknd_t> field_data, dvector<track_t> &tracks,
    const vecmem::jagged_vector<point3_t> &host_positions)
    -> std::tuple<vecmem::jagged_vector<scalar>,
                  vecmem::jagged_vector<point3_t>,
                  vecmem::jagged_vector<free_matrix_t>> {

    // Helper object for performing memory copies.
    vecmem::copy copy;

    // Get tracks data
    auto tracks_data = vecmem::get_data(tracks);

    // Create vector buffer for track recording
    std::vector<std::size_t> sizes(tracks.size(), 0);
    std::vector<std::size_t> capacities;
    for (auto &r : host_positions) {
        capacities.push_back(r.size());
    }

    vecmem::data::jagged_vector_buffer<scalar> path_lengths_buffer(
        capacities, *mr, nullptr, vecmem::data::buffer_type::resizable);
    vecmem::data::jagged_vector_buffer<point3_t> positions_buffer(
        capacities, *mr, nullptr, vecmem::data::buffer_type::resizable);
    vecmem::data::jagged_vector_buffer<free_matrix_t> jac_transports_buffer(
        capacities, *mr, nullptr, vecmem::data::buffer_type::resizable);

    copy.setup(path_lengths_buffer);
    copy.setup(positions_buffer);
    copy.setup(jac_transports_buffer);

    // Run the propagator test for GPU device
    propagator_test<bfield_bknd_t, detector_t>(
        det_view, cfg, field_data, tracks_data,
        path_lengths_buffer, positions_buffer, jac_transports_buffer);

    vecmem::jagged_vector<scalar> device_path_lengths(mr);
    vecmem::jagged_vector<point3_t> device_positions(mr);
    vecmem::jagged_vector<free_matrix_t> device_jac_transports(mr);

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
inline auto run_propagation_test(vecmem::memory_resource *mr, detector_t &det,
                                 const propagator_test_config &cfg,
                                 typename detector_t::view_type det_view,
                                 covfie::field<host_bfield_bknd_t> &&field) {

    // Create the vector of initial track parameterizations
    auto tracks_host = generate_tracks<generator_t>(mr, cfg.track_generator);
    vecmem::vector<track_t> tracks_device(tracks_host, mr);

    // Host propagation
    auto &&[host_path_lengths, host_positions, host_jac_transports] =
        run_propagation_host(mr, det, cfg.propagation, field, tracks_host);

    // Device propagation (device backend specific implementation)
    covfie::field<device_bfield_bknd_t> device_field(field);
    auto &&[device_path_lengths, device_positions, device_jac_transports] =
        run_propagation_device<device_bfield_bknd_t>(
            mr, det, cfg.propagation, det_view, device_field, tracks_device,
            host_positions);

    // Check the results
    compare_propagation_results(host_positions, device_positions,
                                host_path_lengths, device_path_lengths,
                                host_jac_transports, device_jac_transports);
}

}  // namespace detray
