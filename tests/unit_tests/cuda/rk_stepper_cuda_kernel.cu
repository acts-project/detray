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

__global__ void bound_state_test_kernel(
    vecmem::data::vector_view<bound_track_parameters<transform3>> out_param,
    const bound_track_parameters<transform3> in_param, const vector3 B,
    const transform3 trf) {
    /*
    vecmem::device_vector<bound_track_parameters<transform3>> out_param_cuda(
        out_param);

    mag_field_t mag_field(B);
    prop_state<crk_stepper_t::state, nav_state> propagation{
        crk_stepper_t::state(in_param, mag_field), nav_state{}};
    crk_stepper_t::state &crk_state = propagation._stepping;
    nav_state &n_state = propagation._navigation;

    // Decrease tolerance down to 1e-8
    crk_state.set_tolerance(rk_tolerance);

    // RK stepper and its state
    crk_stepper_t crk_stepper;

    // Path length per turn
    scalar S = 2. * std::fabs(1. / in_param.qop()) / getter::norm(B) * M_PI;

    // Run stepper for one turn
    unsigned int max_steps = 1e4;
    for (unsigned int i = 0; i < max_steps; i++) {

        crk_state.set_constraint(S - crk_state.path_length());

        n_state._step_size = S;

        crk_stepper.step(propagation);

        if (std::abs(S - crk_state.path_length()) < 1e-6) {
            break;
        }
    }

    // Bound state after one turn propagation
    out_param_cuda[0] = crk_stepper.bound_state(propagation, trf);
    */
}

void bound_state_test(
    vecmem::data::vector_view<bound_track_parameters<transform3>> out_param,
    const bound_track_parameters<transform3> in_param, const vector3 B,
    const transform3 trf) {

    constexpr int thread_dim = 1;
    constexpr int block_dim = 1;

    // run the test kernel
    bound_state_test_kernel<<<block_dim, thread_dim>>>(out_param, in_param, B,
                                                       trf);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray