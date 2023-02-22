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
    typename detector_host_t<const_bfield_bknd_t>::detector_view_type det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>> tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> candidates_data,
    vecmem::data::jagged_vector_view<scalar> path_lengths_data,
    vecmem::data::jagged_vector_view<vector3> positions_data,
    vecmem::data::jagged_vector_view<free_matrix> jac_transports_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_t<const_bfield_bknd_t> det(det_data);
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

    detector_device_t<const_bfield_bknd_t>::bfield_type B_field = det.get_bfield();


    /*using propagator_device_t = 
        propagator_t<detector_device_t<const_bfield_bknd_t>, actor_chain_device_t>;*/

    // Create RK stepper
    /*rk_stepper_type s;

    // Create navigator
    navigator_device_type n;

    // Create propagator
    propagator_device_type p(std::move(s), std::move(n));

    // Create actor states
    inspector_device_t::state insp_state(
        path_lengths.at(gid), positions.at(gid), jac_transports.at(gid));
    pathlimit_aborter::state aborter_state{path_limit};
    parameter_transporter<transform3>::state transporter_state{};
    pointwise_material_interactor<transform3>::state interactor_state{};
    parameter_resetter<transform3>::state resetter_state{};

    // Create the actor states
    auto actor_states =
        thrust::tie(insp_state, aborter_state, transporter_state,
                    interactor_state, resetter_state);
    // Create the propagator state
    propagator_device_type::state state(tracks[gid], B_field, det,
                                        candidates.at(gid));

    state._stepping.set_tolerance(rk_tolerance);

    state._stepping.template set_constraint<step::constraint::e_accuracy>(
        constrainted_step_size);

    // Run propagation
    p.propagate(state, actor_states);*/
}

template<typename bfield_bknd_t>
void propagator_test(
    typename detector_host_t<bfield_bknd_t>::detector_view_type det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    vecmem::data::jagged_vector_view<scalar>& path_lengths_data,
    vecmem::data::jagged_vector_view<vector3>& positions_data,
    vecmem::data::jagged_vector_view<free_matrix>& jac_transports_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

    // run the test kernel
    propagator_test_kernel<<<block_dim, thread_dim>>>(
        det_data, tracks_data, candidates_data, path_lengths_data,
        positions_data, jac_transports_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}
using queue = decltype(propagator_test<const_bfield_bknd_t>)::bla;

/// Explicit instantiation for constant bfield
template void propagator_test<const_bfield_bknd_t>(
    typename detector_host_t<const_bfield_bknd_t>::detector_view_type,
    vecmem::data::vector_view<free_track_parameters<transform3>>&,
    vecmem::data::jagged_vector_view<intersection_t>&,
    vecmem::data::jagged_vector_view<scalar>&,
    vecmem::data::jagged_vector_view<vector3>&,
    vecmem::data::jagged_vector_view<free_matrix>&);

}  // namespace detray
