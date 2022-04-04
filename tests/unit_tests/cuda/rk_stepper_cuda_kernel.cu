/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/containers/device_vector.hpp>

#include "detray/definitions/cuda_definitions.hpp"
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

    // Define RK [constrained] stepper
    rk_stepper_t rk_stepper(B);
    crk_stepper_t crk_stepper(B);
    nav_state n_state{};

    // Get a track
    auto& traj = tracks.at(gid);

    rk_stepper_t::state rk_state(traj);
    crk_stepper_t::state crk_state(traj);

    // Forward direction
    crk_state.template set_constraint<constraint::e_user>(0.5 *
                                                          unit_constants::mm);
    n_state._step_size = 1. * unit_constants::mm;
    for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
        rk_stepper.step(rk_state, n_state);
        crk_stepper.step(crk_state, n_state);
        crk_stepper.step(crk_state, n_state);
    }

    // Backward direction
    // Roll the same track back to the origin
    scalar path_length = rk_state.path_length();
    n_state._step_size *= -1. * unit_constants::mm;
    for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
        rk_stepper.step(rk_state, n_state);
        crk_stepper.step(crk_state, n_state);
        crk_stepper.step(crk_state, n_state);
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