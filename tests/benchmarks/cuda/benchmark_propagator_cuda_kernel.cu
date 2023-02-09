/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "benchmark_propagator_cuda_kernel.hpp"
#include "detray/definitions/cuda_definitions.hpp"

namespace detray {

__global__ void propagator_benchmark_kernel(
    typename detector_host_type::detector_view_type det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>> tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> candidates_data,
    const propagate_option opt) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_type det(det_data);
    vecmem::device_vector<free_track_parameters<transform3>> tracks(
        tracks_data);
    vecmem::jagged_device_vector<intersection_t> candidates(candidates_data);

    if (gid >= tracks.size()) {
        return;
    }

    // Create RK stepper
    rk_stepper_type s;

    // Create navigator
    navigator_device_type n;

    // Create propagator
    propagator_device_type p(std::move(s), std::move(n));

    parameter_transporter<transform3>::state transporter_state{};
    pointwise_material_interactor<transform3>::state interactor_state{};
    parameter_resetter<transform3>::state resetter_state{};

    // Create the actor states
    auto actor_states =
        thrust::tie(transporter_state, interactor_state, resetter_state);
    // Create the propagator state
    propagator_device_type::state p_state(tracks.at(gid), det.get_bfield(), det,
                                          candidates.at(gid));

    // Run propagation
    if (opt == propagate_option::e_unsync) {
        p.propagate(p_state, actor_states);
    } else if (opt == propagate_option::e_sync) {
        p.propagate_sync(p_state, actor_states);
    }
}

void propagator_benchmark(
    typename detector_host_type::detector_view_type det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    const propagate_option opt) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    int block_dim = tracks_data.size() / thread_dim + 1;

    // run the test kernel
    propagator_benchmark_kernel<<<block_dim, thread_dim>>>(
        det_data, tracks_data, candidates_data, opt);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
