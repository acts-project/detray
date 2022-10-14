/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/vector.hpp>
#include <vecmem/containers/device_vector.hpp>

#include "detray/definitions/cuda_definitions.hpp"
#include "rk_stepper_cuda_kernel.hpp"

namespace {
using field_t = covfie::field<
    covfie::backend::constant<covfie::vector::vector_d<detray::scalar, 3>,
                              covfie::vector::vector_d<detray::scalar, 3>>>;
}

namespace detray {

__global__ void bound_state_test_kernel(
    vecmem::data::vector_view<bound_track_parameters<transform3>> out_param,
    const bound_track_parameters<transform3> in_param, const field_t::view_t B,
    const transform3 trf) {

    vecmem::device_vector<bound_track_parameters<transform3>> out_param_cuda(
        out_param);

    prop_state<crk_stepper_t::state, nav_state> propagation{
        crk_stepper_t::state(in_param, trf, B), nav_state{}};
    crk_stepper_t::state &crk_state = propagation._stepping;
    nav_state &n_state = propagation._navigation;

    // Decrease tolerance down to 1e-8
    crk_state.set_tolerance(rk_tolerance);

    // RK stepper and its state
    crk_stepper_t crk_stepper;

    vector3 Bv{B.at(0.f, 0.f, 0.f)[0], B.at(0.f, 0.f, 0.f)[1],
               B.at(0.f, 0.f, 0.f)[2]};

    // Path length per turn
    scalar S = 2. * std::fabs(1. / in_param.qop()) / getter::norm(Bv) * M_PI;

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
}

void bound_state_test(
    vecmem::data::vector_view<bound_track_parameters<transform3>> out_param,
    const bound_track_parameters<transform3> in_param, const vector3 B,
    const transform3 trf) {

    constexpr int thread_dim = 1;
    constexpr int block_dim = 1;

    field_t f(field_t::backend_t::configuration_t{B[0], B[1], B[2]});

    // run the test kernel
    bound_state_test_kernel<<<block_dim, thread_dim>>>(out_param, in_param, f,
                                                       trf);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
