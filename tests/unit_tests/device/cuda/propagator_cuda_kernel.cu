/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "propagator_cuda_kernel.hpp"

namespace detray {

template <typename detector_t, typename bfield_bknd_t>
__global__ void propagator_test_kernel(
    detector_view<covfie::field<bfield_bknd_t>, detector_t> det_data,
    vecmem::data::vector_view<track_t> tracks_data,
    vecmem::data::jagged_vector_view<intersection_t<detector_t>>
        candidates_data,
    vecmem::data::jagged_vector_view<scalar> path_lengths_data,
    vecmem::data::jagged_vector_view<vector3_t> positions_data,
    vecmem::data::jagged_vector_view<free_matrix> jac_transports_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_t<bfield_bknd_t> det(det_data);
    vecmem::device_vector<track_t> tracks(tracks_data);
    vecmem::jagged_device_vector<intersection_t<detector_t>> candidates(
        candidates_data);
    vecmem::jagged_device_vector<scalar> path_lengths(path_lengths_data);
    vecmem::jagged_device_vector<vector3_t> positions(positions_data);
    vecmem::jagged_device_vector<free_matrix> jac_transports(
        jac_transports_data);

    if (gid >= tracks.size()) {
        return;
    }

    auto stepr =
        rk_stepper_t<typename detector_device_t<bfield_bknd_t>::bfield_type>{};
    auto nav = navigator_t<detector_device_t<bfield_bknd_t>>{};

    // Create propagator
    using propagator_device_t =
        propagator<decltype(stepr), decltype(nav), actor_chain_device_t>;
    propagator_device_t p(std::move(stepr), std::move(nav));

    // Create actor states
    inspector_device_t::state insp_state(
        path_lengths.at(gid), positions.at(gid), jac_transports.at(gid));
    pathlimit_aborter::state aborter_state{path_limit};
    parameter_transporter<transform3>::state transporter_state{};
    pointwise_material_interactor<transform3>::state interactor_state{};
    parameter_resetter<transform3>::state resetter_state{};

    // Create the actor states
    auto actor_states =
        ::detray::tie(insp_state, aborter_state, transporter_state,
                      interactor_state, resetter_state);
    // Create the propagator state
    typename propagator_device_t::state state(tracks[gid], det.get_bfield(),
                                              det, candidates.at(gid));

    state._stepping.set_tolerance(rk_tolerance);

    state._stepping.template set_constraint<step::constraint::e_accuracy>(
        constrainted_step_size);

    // Run propagation
    p.propagate(state, actor_states);
}

/// Launch the device kernel
template <typename detector_t, typename bfield_bknd_t>
void propagator_test(
    detector_view<covfie::field<bfield_bknd_t>, detector_t> det_view,
    vecmem::data::vector_view<track_t>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t<detector_t>>&
        candidates_data,
    vecmem::data::jagged_vector_view<scalar>& path_lengths_data,
    vecmem::data::jagged_vector_view<vector3_t>& positions_data,
    vecmem::data::jagged_vector_view<free_matrix>& jac_transports_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

    // run the test kernel
    propagator_test_kernel<<<block_dim, thread_dim>>>(
        det_view, tracks_data, candidates_data, path_lengths_data,
        positions_data, jac_transports_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Explicit instantiation for a constant magnetic field
template void propagator_test<detector_host_t<const_backend_t>, const_backend_t>(
    detector_view<covfie::field<const_backend_t>, detector_host_t<const_backend_t>> det_view,
    vecmem::data::vector_view<track_t>&,
    vecmem::data::jagged_vector_view<intersection_t<detector_host_t<const_backend_t>>>&,
    vecmem::data::jagged_vector_view<scalar>&,
    vecmem::data::jagged_vector_view<vector3_t>&,
    vecmem::data::jagged_vector_view<free_matrix>&);

/// Explicit instantiation for an inhomogeneous magnetic field
template void propagator_test<detector_host_t<inhom_bfield_bknd_t>, inhom_bfield_bknd_cuda_t>(
    detector_view<covfie::field<inhom_bfield_bknd_cuda_t>, detector_host_t<inhom_bfield_bknd_t>> det_view,
    vecmem::data::vector_view<track_t>&,
    vecmem::data::jagged_vector_view<intersection_t<detector_host_t<inhom_bfield_bknd_t>>>&,
    vecmem::data::jagged_vector_view<scalar>&,
    vecmem::data::jagged_vector_view<vector3_t>&,
    vecmem::data::jagged_vector_view<free_matrix>&);

}  // namespace detray
