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
    vecmem::data::vector_view<intersection> candidates_data) {

    navigator_device_t n(n_data);
    navigator_device_t::state state(candidates_data);
}

void navigator_test(navigator_view<navigator_host_t> n_data,
                    vecmem::data::vector_view<intersection>& candidates_data) {

    constexpr int block_dim = 1;
    constexpr int thread_dim = 1;

    // run the test kernel
    navigator_test_kernel<<<block_dim, thread_dim>>>(n_data, candidates_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

__global__ void geometry_navigation_kernel(
    navigator_view<navigator_host_t> n_data,
    vecmem::data::jagged_vector_view<intersection> candidates_data,
    vecmem::data::vector_view<track<nav_context>> tracks_data,
    vecmem::data::jagged_vector_view<std::pair<dindex, intersection>>
        intersection_record_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    if (gid >= candidates_data.m_size) {
        return;
    }

    navigator_device_t n(n_data);

    detector_device_t& det = n.get_detector();

    navigator_device_t::state state(candidates_data.m_ptr[gid]);

    vecmem::device_vector<track<nav_context>> tracks(tracks_data);

    vecmem::device_vector<std::pair<dindex, intersection>> intersection_trace(
        intersection_record_data.m_ptr[gid]);

    /*
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        printf("%d \n",
               thrust::tuple_size<
                   detector_device_t::mask_container::mask_tuple>::value);
    }
    */

    // shoot_ray(n.get_detector(), tracks[gid], intersection_trace);

    /*
    state.candidates().push_back(intersection{});

    printf("%lu %d \n", static_cast<long unsigned
    int>(state.candidates()[0].index), state.candidates().size());
    */
}

void geometry_navigation_test(
    navigator_view<navigator_host_t> n_data,
    vecmem::data::jagged_vector_view<intersection>& candidates_data,
    vecmem::data::vector_view<track<nav_context>>& tracks_data,
    vecmem::data::jagged_vector_view<std::pair<dindex, intersection>>&
        intersection_record_data) {

    constexpr int block_dim = theta_steps;
    constexpr int thread_dim = phi_steps;

    // run the test kernel
    geometry_navigation_kernel<<<block_dim, thread_dim>>>(
        n_data, candidates_data, tracks_data, intersection_record_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray