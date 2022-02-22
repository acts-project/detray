/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/containers/device_vector.hpp>

#include "detray/definitions/cuda_defs.hpp"
#include "rk_stepper_cuda_kernel.hpp"

namespace detray {

__global__ void rk_stepper_test_kernel(
    vecmem::data::vector_view<free_track_parameters> tracks_data,
    const vector3 B) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    vecmem::device_vector<free_track_parameters> tracks(tracks_data);

    // Prevent overflow
    if (gid >= tracks.size()) {
        return;
    }

    // Define RK stepper
    rk_stepper_type rk(B);

    // Get a track
    auto& traj = tracks.at(gid);

    // Forward direction
    rk_stepper_type::state forward_state(traj);
    for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
        rk.step(forward_state, path_limit);
    }

    // Backward direction
    traj.flip();
    rk_stepper_type::state backward_state(traj);
    for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
        rk.step(backward_state, path_limit);
    }
}

void rk_stepper_test(
    vecmem::data::vector_view<free_track_parameters>& tracks_data,
    const vector3 B) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

    // run the test kernel
    rk_stepper_test_kernel<<<block_dim, thread_dim>>>(tracks_data, B);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray