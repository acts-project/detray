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
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters> tracks_data,
    vecmem::data::jagged_vector_view<intersection> candidates_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_type det(det_data);
    vecmem::device_vector<free_track_parameters> tracks(tracks_data);
    vecmem::jagged_device_vector<intersection> candidates(candidates_data);

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

    propagation::void_inspector vi;

    // Create the propagator state
    propagator_device_type::state p_state(tracks.at(gid), candidates.at(gid));

    // Run propagation
    p.propagate(p_state, vi);
}

void propagator_benchmark(
    detector_view<detector_host_type> det_data,
    vecmem::data::vector_view<free_track_parameters>& tracks_data,
    vecmem::data::jagged_vector_view<intersection>& candidates_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    int block_dim = tracks_data.size() / thread_dim + 1;

    // run the test kernel
    propagator_benchmark_kernel<<<block_dim, thread_dim>>>(
        det_data, tracks_data, candidates_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray