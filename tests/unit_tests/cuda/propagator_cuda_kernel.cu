/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "propagator_cuda_kernel.hpp"

namespace detray {

__global__ void propagator_test_kernel(
    detector_view<detector_host_type> det_data, const vector3 B,
    vecmem::data::vector_view<free_track_parameters<transform3>> tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> candidates_data,
    vecmem::data::jagged_vector_view<scalar> path_lengths_data,
    vecmem::data::jagged_vector_view<vector3> positions_data,
    vecmem::data::jagged_vector_view<free_matrix> jac_transports_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_type det(det_data);
    vecmem::device_vector<free_track_parameters<transform3>> tracks(
        tracks_data);
    vecmem::jagged_device_vector<intersection_t> candidates(candidates_data);
    vecmem::jagged_device_vector<scalar> path_lengths(path_lengths_data);
    vecmem::jagged_device_vector<vector3> positions(positions_data);
    vecmem::jagged_device_vector<free_matrix> jac_transports(
        jac_transports_data);

    if (gid >= tracks.size()) {
        return;
    }

    // Set the magnetic field
    field_type B_field(B);

    // Create RK stepper
    rk_stepper_type s;

    // Create navigator
    navigator_device_type n;

    // Create propagator
    propagator_device_type p(std::move(s), std::move(n));

    // Create actor states
    inspector_device_t::state insp_state(
        path_lengths.at(gid), positions.at(gid), jac_transports.at(gid));
    pathlimit_aborter::state aborter_state{path_limit};

    // Create the propagator state
    propagator_device_type::state state(tracks[gid], B_field, det,
                                        thrust::tie(insp_state, aborter_state),
                                        candidates.at(gid));

    state._stepping.set_tolerance(rk_tolerance);

    state._stepping.template set_constraint<step::constraint::e_accuracy>(
        constrainted_step_size);

    // Run propagation
    p.propagate(state);
}

void propagator_test(
    detector_view<detector_host_type> det_data, const vector3 B,
    vecmem::data::vector_view<free_track_parameters<transform3>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    vecmem::data::jagged_vector_view<scalar>& path_lengths_data,
    vecmem::data::jagged_vector_view<vector3>& positions_data,
    vecmem::data::jagged_vector_view<free_matrix>& jac_transports_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

    // run the test kernel
    propagator_test_kernel<<<block_dim, thread_dim>>>(
        det_data, B, tracks_data, candidates_data, path_lengths_data,
        positions_data, jac_transports_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray