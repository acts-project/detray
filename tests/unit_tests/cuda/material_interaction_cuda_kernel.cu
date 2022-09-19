/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "material_interaction_cuda_kernel.hpp"

namespace detray {

__global__ void material_interaction_test_kernel(
    detector_view<detector_host_type> det_data, const vector3 B,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    vecmem::data::vector_view<bound_track_parameters<transform3>>
        initial_states,
    vecmem::data::vector_view<bound_track_parameters<transform3>>
        final_states) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_type det(det_data);
    vecmem::jagged_device_vector<intersection_t> candidates(candidates_data);
    vecmem::device_vector<bound_track_parameters<transform3>>
        device_initial_states(initial_states);
    vecmem::device_vector<bound_track_parameters<transform3>>
        device_final_states(final_states);

    if (gid >= device_final_states.size()) {
        return;
    }

    // Set the magnetic field
    field_type B_field(B);

    // Propagator is built from the stepper and navigator
    propagator_device_type p({}, {});

    // Actor states
    pathlimit_aborter::state aborter_state{};
    bound_to_bound_updater<transform3>::state bound_updater{};
    pointwise_material_interactor<transform3>::state interactor_state{};
    resetter<transform3>::state resetter_state{};

    // Create actor states tuples
    actor_chain_t::state actor_states = thrust::tie(
        aborter_state, bound_updater, interactor_state, resetter_state);

    propagator_device_type::state state(device_initial_states.at(gid), B_field,
                                        det, actor_states, candidates.at(gid));

    state._stepping().set_overstep_tolerance(overstep_tolerance);

    // Run propagation
    p.propagate(state);

    device_final_states[gid] = state._stepping._bound_params;
}

void material_interaction_test(
    detector_view<detector_host_type> det_data, const vector3 B,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    vecmem::data::vector_view<bound_track_parameters<transform3>>
        initial_states,
    vecmem::data::vector_view<bound_track_parameters<transform3>>
        final_states) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = n_tracks / thread_dim + 1;

    // run the test kernel
    material_interaction_test_kernel<<<block_dim, thread_dim>>>(
        det_data, B, candidates_data, initial_states, final_states);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray