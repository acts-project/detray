/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "propagator_cuda_kernel.hpp"

namespace detray {

template <typename bfield_bknd_t, typename detector_t>
__global__ void propagator_test_kernel(
    typename detector_t::view_type det_data,
    covfie::field_view<bfield_bknd_t> field_data,
    vecmem::data::vector_view<track_t> tracks_data,
    vecmem::data::jagged_vector_view<intersection_t<detector_t>>
        candidates_data,
    vecmem::data::jagged_vector_view<scalar_t> path_lengths_data,
    vecmem::data::jagged_vector_view<vector3_t> positions_data,
    vecmem::data::jagged_vector_view<free_matrix_t> jac_transports_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    using detector_device_t =
        detector<typename detector_t::metadata, device_container_types>;

    static_assert(std::is_same_v<typename detector_t::view_type,
                                 typename detector_device_t::view_type>,
                  "Host and device detector views do not match");

    detector_device_t det(det_data);
    vecmem::device_vector<track_t> tracks(tracks_data);
    vecmem::jagged_device_vector<intersection_t<detector_t>> candidates(
        candidates_data);
    vecmem::jagged_device_vector<scalar> path_lengths(path_lengths_data);
    vecmem::jagged_device_vector<vector3_t> positions(positions_data);
    vecmem::jagged_device_vector<free_matrix_t> jac_transports(
        jac_transports_data);

    if (gid >= tracks.size()) {
        return;
    }

    auto stepr = rk_stepper_t<covfie::field_view<bfield_bknd_t>>{};
    auto nav = navigator_t<detector_device_t>{};

    // Create propagator
    using propagator_device_t =
        propagator<decltype(stepr), decltype(nav), actor_chain_device_t>;

    propagation::config<scalar> cfg{};
    cfg.navigation.search_window = {3u, 3u};
    cfg.stepping.rk_error_tol = rk_tolerance;
    propagator_device_t p{cfg};

    // Create actor states
    inspector_device_t::state insp_state(
        path_lengths.at(gid), positions.at(gid), jac_transports.at(gid));
    pathlimit_aborter::state aborter_state{path_limit};
    parameter_transporter<algebra_t>::state transporter_state{};
    pointwise_material_interactor<algebra_t>::state interactor_state{};
    parameter_resetter<algebra_t>::state resetter_state{};

    // Create the actor states
    auto actor_states =
        ::detray::tie(insp_state, aborter_state, transporter_state,
                      interactor_state, resetter_state);
    // Create the propagator state
    typename propagator_device_t::state state(tracks[gid], field_data, det,
                                              candidates.at(gid));

    state._stepping.template set_constraint<step::constraint::e_accuracy>(
        constrainted_step_size);

    // Run propagation
    p.propagate(state, actor_states);
}

/// Launch the device kernel
template <typename bfield_bknd_t, typename detector_t>
void propagator_test(
    typename detector_t::view_type det_view,
    covfie::field_view<bfield_bknd_t> field_data,
    vecmem::data::vector_view<track_t>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t<detector_t>>&
        candidates_data,
    vecmem::data::jagged_vector_view<scalar_t>& path_lengths_data,
    vecmem::data::jagged_vector_view<vector3_t>& positions_data,
    vecmem::data::jagged_vector_view<free_matrix_t>& jac_transports_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

    // run the test kernel
    propagator_test_kernel<bfield_bknd_t, detector_t>
        <<<block_dim, thread_dim>>>(det_view, field_data, tracks_data,
                                    candidates_data, path_lengths_data,
                                    positions_data, jac_transports_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Explicit instantiation for a constant magnetic field
template void propagator_test<bfield::const_bknd_t,
                              detector<toy_metadata, host_container_types>>(
    detector<toy_metadata, host_container_types>::view_type,
    covfie::field_view<bfield::const_bknd_t>,
    vecmem::data::vector_view<track_t>&,
    vecmem::data::jagged_vector_view<
        intersection_t<detector<toy_metadata, host_container_types>>>&,
    vecmem::data::jagged_vector_view<scalar_t>&,
    vecmem::data::jagged_vector_view<vector3_t>&,
    vecmem::data::jagged_vector_view<free_matrix_t>&);

/// Explicit instantiation for an inhomogeneous magnetic field
template void propagator_test<bfield::cuda::inhom_bknd_t,
                              detector<toy_metadata, host_container_types>>(
    detector<toy_metadata, host_container_types>::view_type,
    covfie::field_view<bfield::cuda::inhom_bknd_t>,
    vecmem::data::vector_view<track_t>&,
    vecmem::data::jagged_vector_view<
        intersection_t<detector<toy_metadata, host_container_types>>>&,
    vecmem::data::jagged_vector_view<scalar_t>&,
    vecmem::data::jagged_vector_view<vector3_t>&,
    vecmem::data::jagged_vector_view<free_matrix_t>&);

}  // namespace detray
