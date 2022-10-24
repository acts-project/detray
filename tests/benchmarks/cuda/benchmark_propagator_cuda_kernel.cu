/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "benchmark_propagator_cuda_kernel.hpp"
#include "detray/definitions/cuda_definitions.hpp"

namespace detray {

template <typename stepper_policy_t>
__global__ void propagator_benchmark_kernel(
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>> tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> candidates_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_type det(det_data);
    vecmem::device_vector<free_track_parameters<transform3>> tracks(
        tracks_data);
    vecmem::jagged_device_vector<intersection_t> candidates(candidates_data);

    if (gid >= tracks.size()) {
        return;
    }

    // Create RK stepper
    rk_stepper_type<stepper_policy_t> s;

    // Create navigator
    navigator_device_type n;

    // Create propagator
    propagator_device_type<stepper_policy_t> p(std::move(s), std::move(n));

    parameter_transporter<transform3>::state transporter_state{};
    pointwise_material_interactor<transform3>::state interactor_state{};
    parameter_resetter<transform3>::state resetter_state{};

    // Create the propagator state
    typename propagator_device_type<stepper_policy_t>::state p_state(
        tracks.at(gid), det.get_bfield(), det,
        thrust::tie(transporter_state, interactor_state, resetter_state),
        candidates.at(gid));

    // Run propagation
    p.propagate(p_state);
}

template <typename stepper_policy_t>
void propagator_benchmark(
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    int block_dim = tracks_data.size() / thread_dim + 1;

    // run the test kernel
    propagator_benchmark_kernel<stepper_policy_t>
        <<<block_dim, thread_dim>>>(det_data, tracks_data, candidates_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

template void propagator_benchmark<always_init>(
    detector_view<detector_host_type>,
    vecmem::data::vector_view<free_track_parameters<transform3>>&,
    vecmem::data::jagged_vector_view<intersection_t>&);

template void propagator_benchmark<stepper_default_policy>(
    detector_view<detector_host_type>,
    vecmem::data::vector_view<free_track_parameters<transform3>>&,
    vecmem::data::jagged_vector_view<intersection_t>&);

template void propagator_benchmark<stepper_rk_policy>(
    detector_view<detector_host_type>,
    vecmem::data::vector_view<free_track_parameters<transform3>>&,
    vecmem::data::jagged_vector_view<intersection_t>&);

}  // namespace detray
