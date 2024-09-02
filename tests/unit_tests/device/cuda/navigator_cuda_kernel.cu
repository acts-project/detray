/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "navigator_cuda_kernel.hpp"

namespace detray {

__global__ void navigator_test_kernel(
    typename detector_host_t::view_type det_data, navigation::config cfg,
    vecmem::data::vector_view<free_track_parameters<algebra_t>> tracks_data,
    vecmem::data::jagged_vector_view<dindex> volume_records_data,
    vecmem::data::jagged_vector_view<point3> position_records_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    detector_device_t det(det_data);
    vecmem::device_vector<free_track_parameters<algebra_t>> tracks(tracks_data);
    vecmem::jagged_device_vector<dindex> volume_records(volume_records_data);
    vecmem::jagged_device_vector<point3> position_records(
        position_records_data);

    if (gid >= tracks.size()) {
        return;
    }

    navigator_device_t nav;

    auto& traj = tracks.at(gid);
    stepper_t stepper;

    prop_state<navigator_device_t::state> propagation{
        stepper_t::state{traj}, navigator_device_t::state(det)};

    navigator_device_t::state& navigation = propagation._navigation;
    stepper_t::state& stepping = propagation._stepping;

    // Set initial volume
    navigation.set_volume(0u);

    // Start propagation and record volume IDs
    bool heartbeat = nav.init(propagation, cfg);
    while (heartbeat) {

        heartbeat &= stepper.step(propagation);

        navigation.set_high_trust();

        heartbeat = nav.update(propagation, cfg);

        // Record volume
        volume_records[gid].push_back(navigation.volume());
        position_records[gid].push_back(stepping().pos());
    }
}

void navigator_test(
    typename detector_host_t::view_type det_data, navigation::config& cfg,
    vecmem::data::vector_view<free_track_parameters<algebra_t>>& tracks_data,
    vecmem::data::jagged_vector_view<dindex>& volume_records_data,
    vecmem::data::jagged_vector_view<point3>& position_records_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

    // run the test kernel
    navigator_test_kernel<<<block_dim, thread_dim>>>(
        det_data, cfg, tracks_data, volume_records_data, position_records_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
