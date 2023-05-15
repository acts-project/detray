/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "propagation.hpp"

namespace detray::tutorial {

// Propagation configurations
inline constexpr detray::scalar path_limit{2.f *
                                           detray::unit<detray::scalar>::m};

/// Kernel that runs the entire propagation loop
__global__ void propagation_kernel(
    typename detray::tutorial::detector_host_t::detector_view_type det_data,
    const vecmem::data::vector_view<
        detray::free_track_parameters<detray::tutorial::transform3>>
        tracks_data,
    vecmem::data::jagged_vector_view<detray::tutorial::intersection_t>
        candidates_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    // Setup device-side track collection
    vecmem::device_vector<
        detray::free_track_parameters<detray::tutorial::transform3>>
        tracks(tracks_data);

    if (gid >= tracks.size()) {
        return;
    }

    // Setup of the device-side detector
    detray::tutorial::detector_device_t det(det_data);
    // Setup of the avigator cache
    vecmem::jagged_device_vector<detray::tutorial::intersection_t> candidates(
        candidates_data);
    // Setup of the device b-field
    detray::tutorial::detector_device_t::bfield_type B_field = det.get_bfield();

    // Create propagator from a stepper and a navigator
    detray::tutorial::propagator_t p(detray::tutorial::stepper_t{},
                                     detray::tutorial::navigator_t{});

    // Create actor states
    detray::pathlimit_aborter::state aborter_state{path_limit};
    detray::parameter_transporter<detray::tutorial::transform3>::state
        transporter_state{};
    detray::pointwise_material_interactor<detray::tutorial::transform3>::state
        interactor_state{};
    detray::parameter_resetter<detray::tutorial::transform3>::state
        resetter_state{};

    auto actor_states = detray::tie(aborter_state, transporter_state,
                                    interactor_state, resetter_state);

    // Create the propagator state for the track
    detray::tutorial::propagator_t::state state(tracks[gid], B_field, det,
                                                candidates.at(gid));

    // Run propagation
    p.propagate(state, actor_states);
}

void propagation(
    typename detray::tutorial::detector_host_t::detector_view_type det_data,
    const vecmem::data::vector_view<
        detray::free_track_parameters<detray::tutorial::transform3>>
        tracks_data,
    vecmem::data::jagged_vector_view<detray::tutorial::intersection_t>
        candidates_data) {

    int thread_dim = 2 * WARP_SIZE;
    int block_dim = tracks_data.size() / thread_dim + 1;

    // run the tutorial kernel
    propagation_kernel<<<block_dim, thread_dim>>>(det_data, tracks_data,
                                                  candidates_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray::tutorial
