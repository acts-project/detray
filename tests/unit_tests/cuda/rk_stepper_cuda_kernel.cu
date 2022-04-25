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

    // Get a track
    auto &traj = tracks.at(gid);
    free_track_parameters c_traj(traj);

    // RK Stepping into forward direction
    prop_state<rk_stepper_t::state> propagation{rk_stepper_t::state{traj},
                                                nav_state{}};
    prop_state<crk_stepper_t::state> c_propagation{crk_stepper_t::state{c_traj},
                                                   nav_state{}};

    crk_stepper_t::state &crk_state = c_propagation._stepping;

    nav_state &n_state = propagation._navigation;
    nav_state &cn_state = c_propagation._navigation;

    // Forward direction
    crk_state.template set_constraint<step::constraint::e_user>(
        0.5 * unit_constants::mm);
    n_state._step_size = 1. * unit_constants::mm;
    cn_state._step_size = 1. * unit_constants::mm;
    for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
        rk_stepper.step(propagation);
        crk_stepper.step(c_propagation);
        crk_stepper.step(c_propagation);
    }

    // Backward direction
    // Roll the same track back to the origin
    n_state._step_size *= -1. * unit_constants::mm;
    cn_state._step_size = 1. * unit_constants::mm;
    for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
        rk_stepper.step(propagation);
        crk_stepper.step(c_propagation);
        crk_stepper.step(c_propagation);
    }
}

void rk_stepper_test(
    vecmem::data::vector_view<free_track_parameters> &tracks_data,
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