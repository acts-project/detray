/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_defs.hpp"
#include "navigator_cuda_kernel.hpp"

namespace detray {

__global__ void navigator_test_kernel(
    navigator_view<navigator_host_t> n_data,
    vecmem::data::vector_view<track<nav_context>> tracks_data,
    vecmem::data::jagged_vector_view<intersection> candidates_data,
    vecmem::data::jagged_vector_view<dindex> volume_records_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    navigator_device_t n(n_data);

    vecmem::device_vector<track<nav_context>> tracks(tracks_data);

    auto& traj = tracks[gid];

    navigator_device_t::state state(candidates_data.m_ptr[gid]);

    vecmem::jagged_device_vector<dindex> volume_records(volume_records_data);

    if (gid >= tracks.size()) {
        return;
    }

    // Set initial volume
    state.set_volume(0u);

    // Start propagation and record volume IDs
    bool heartbeat = n.status(state, traj);

    while (heartbeat) {
        heartbeat = n.target(state, traj);

        traj.pos = traj.pos + state() * traj.dir;

        heartbeat = n.status(state, traj);

        // Record volume
        volume_records[gid].push_back(state.volume());
    }
}

void navigator_test(
    navigator_view<navigator_host_t> n_data,
    vecmem::data::vector_view<track<nav_context>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection>& candidates_data,
    vecmem::data::jagged_vector_view<dindex>& volume_records_data) {

    constexpr int block_dim = theta_steps;
    constexpr int thread_dim = phi_steps;

    // run the test kernel
    navigator_test_kernel<<<block_dim, thread_dim>>>(
        n_data, tracks_data, candidates_data, volume_records_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray