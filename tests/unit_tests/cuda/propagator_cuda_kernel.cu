/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "propagator_cuda_kernel.hpp"

namespace detray {

__global__ void propagator_test_kernel(
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters> tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> candidates_data,
    vecmem::data::jagged_vector_view<intersection_t> intersections_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_type det(det_data);
    vecmem::device_vector<free_track_parameters> tracks(tracks_data);
    vecmem::jagged_device_vector<intersection_t> candidates(candidates_data);
    vecmem::jagged_device_vector<intersection_t> intersections(
        intersections_data);

    if (gid >= tracks.size()) {
        return;
    }

    // Set the magnetic field
    const vector3 B{0, 0, 2 * unit_constants::T};
    field_type B_field(B);

    // Create RK stepper
    rk_stepper_type s(B_field);

    // Create navigator
    navigator_device_type n(det);

    // Create propagator
    propagator_device_type p(std::move(s), std::move(n));

    // Create track inspector
    auto actor_states = thrust::make_tuple(
        inspector_t::state<vecmem::device_vector>(intersections.at(gid)));

    // Create the propagator state
    propagator_device_type::state<decltype(actor_states)> state(
        tracks[gid], std::move(actor_states), candidates.at(gid));

    // Run propagation
    p.propagate(state);
}

void propagator_test(
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    vecmem::data::jagged_vector_view<intersection_t>& intersections_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

    // run the test kernel
    propagator_test_kernel<<<block_dim, thread_dim>>>(
        det_data, tracks_data, candidates_data, intersections_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray