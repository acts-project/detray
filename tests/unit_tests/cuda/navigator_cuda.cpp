/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "detray/definitions/units.hpp"
#include "navigator_cuda_kernel.hpp"
#include "tests/common/tools/track_generators.hpp"
#include "vecmem/utils/cuda/copy.hpp"

TEST(navigator_cuda, navigator) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    // Create detector
    detector_host_t det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                                   n_edc_layers);

    // Create navigator
    navigator_host_t nav(det);

    // Create the vector of initial track parameters
    vecmem::vector<free_track_parameters> tracks_host(&mng_mr);
    vecmem::vector<free_track_parameters> tracks_device(&mng_mr);

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};
    const scalar mom_mag = 10. * unit_constants::GeV;

    // Iterate through uniformly distributed momentum directions
    for (auto traj : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps, ori, mom_mag)) {

        tracks_host.push_back(ray);
        tracks_device.push_back(ray);
    }

    /**
     * Host Volume Record
     */
    vecmem::jagged_vector<dindex> volume_records_host(theta_steps * phi_steps,
                                                      &mng_mr);
    vecmem::jagged_vector<point3> position_records_host(theta_steps * phi_steps,
                                                        &mng_mr);

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {

        auto& traj = tracks_host[i];
        stepper_t stepper;

        prop_state<navigator_host_t::state> propagation{
            stepper_t::state{traj}, navigator_host_t::state{mng_mr}};

        navigator_host_t::state& navigation = propagation._navigation;
        stepper_t::state& stepping = propagation._stepping;

        // Start propagation and record volume IDs
        bool heartbeat = nav.init(propagation);
        while (heartbeat) {

            heartbeat &= stepper.step(propagation);

            navigation.set_high_trust();

            heartbeat &= nav.update(propagation);

            // Record volume
            volume_records_host[i].push_back(navigation.volume());
            position_records_host[i].push_back(stepping().pos());
        }
    }

    /**
     * Device Volume Record
     */

    vecmem::jagged_vector<dindex> volume_records_device(&mng_mr);
    vecmem::jagged_vector<point3> position_records_device(&mng_mr);

    // Create size and capacity vectors for volume record buffer
    std::vector<size_t> sizes;
    std::vector<size_t> capacities;

    for (unsigned int i = 0; i < theta_steps * phi_steps; i++) {
        sizes.push_back(0);
        capacities.push_back(volume_records_host[i].size());
    }

    vecmem::data::jagged_vector_buffer<dindex> volume_records_buffer(
        sizes, capacities, dev_mr, &mng_mr);
    copy.setup(volume_records_buffer);

    vecmem::data::jagged_vector_buffer<point3> position_records_buffer(
        sizes, capacities, dev_mr, &mng_mr);
    copy.setup(position_records_buffer);

    // Get detector data
    auto det_data = get_data(det);

    // Get tracks data
    auto tracks_data = vecmem::get_data(tracks_device);

    // Create navigator candidates buffer
    auto candidates_buffer =
        create_candidates_buffer(det, theta_steps * phi_steps, dev_mr);
    copy.setup(candidates_buffer);

    // Run navigator test
    navigator_test(det_data, tracks_data, candidates_buffer,
                   volume_records_buffer, position_records_buffer);

    // Copy volume record buffer into volume & position records device
    copy(volume_records_buffer, volume_records_device);
    copy(position_records_buffer, position_records_device);

    for (unsigned int i = 0; i < volume_records_host.size(); i++) {

        EXPECT_EQ(volume_records_host[i].size(),
                  volume_records_device[i].size());

        for (unsigned int j = 0; j < volume_records_host[i].size(); j++) {

            EXPECT_EQ(volume_records_host[i][j], volume_records_device[i][j]);

            auto& pos_host = position_records_host[i][j];
            auto& pos_device = position_records_device[i][j];

            EXPECT_NEAR(pos_host[0], pos_device[0], pos_diff_tolerance);
            EXPECT_NEAR(pos_host[1], pos_device[1], pos_diff_tolerance);
            EXPECT_NEAR(pos_host[2], pos_device[2], pos_diff_tolerance);
        }
    }
}