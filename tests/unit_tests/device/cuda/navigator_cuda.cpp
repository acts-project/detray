/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/build_toy_detector.hpp"
#include "navigator_cuda_kernel.hpp"

// vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s)
#include <gtest/gtest.h>

TEST(navigator_cuda, navigator) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    // Create detector
    auto [det, names] = build_toy_detector(mng_mr);

    // Create navigator
    navigator_host_t nav;
    navigation::config cfg{};
    cfg.search_window = {1u, 1u};

    // Create the vector of initial track parameters
    vecmem::vector<free_track_parameters<algebra_t>> tracks_host(&mng_mr);
    vecmem::vector<free_track_parameters<algebra_t>> tracks_device(&mng_mr);

    // Magnitude of total momentum of tracks
    const scalar_t p_mag{10.f * unit<scalar_t>::GeV};

    // Iterate through uniformly distributed momentum directions
    for (auto track : uniform_track_generator<free_track_parameters<algebra_t>>(
             phi_steps, theta_steps, p_mag)) {
        tracks_host.push_back(track);
        tracks_device.push_back(track);
    }

    /**
     * Host Volume Record
     */
    vecmem::jagged_vector<dindex> volume_records_host(theta_steps * phi_steps,
                                                      &mng_mr);
    vecmem::jagged_vector<point3> position_records_host(theta_steps * phi_steps,
                                                        &mng_mr);

    for (unsigned int i = 0u; i < theta_steps * phi_steps; i++) {

        auto& track = tracks_host[i];
        stepper_t stepper;

        // Propagator is built from the stepper and navigator
        prop_state<navigator_host_t::state> propagation{
            stepper_t::state{track}, navigator_host_t::state(det)};

        navigator_host_t::state& navigation = propagation._navigation;
        stepper_t::state& stepping = propagation._stepping;

        // Start propagation and record volume IDs
        bool heartbeat = nav.init(propagation, cfg);
        while (heartbeat) {

            heartbeat &= stepper.step(propagation);

            navigation.set_high_trust();

            heartbeat &= nav.update(propagation, cfg);

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
    std::vector<size_t> capacities;

    for (unsigned int i = 0u; i < theta_steps * phi_steps; i++) {
        capacities.push_back(volume_records_host[i].size());
    }

    vecmem::data::jagged_vector_buffer<dindex> volume_records_buffer(
        capacities, dev_mr, &mng_mr, vecmem::data::buffer_type::resizable);
    copy.setup(volume_records_buffer);

    vecmem::data::jagged_vector_buffer<point3> position_records_buffer(
        capacities, dev_mr, &mng_mr, vecmem::data::buffer_type::resizable);
    copy.setup(position_records_buffer);

    // Get detector data
    auto det_data = detray::get_data(det);

    // Get tracks data
    auto tracks_data = vecmem::get_data(tracks_device);

    // Run navigator test
    navigator_test(det_data, cfg, tracks_data, volume_records_buffer,
                   position_records_buffer);

    // Copy volume record buffer into volume & position records device
    copy(volume_records_buffer, volume_records_device);
    copy(position_records_buffer, position_records_device);

    for (unsigned int i = 0u; i < volume_records_host.size(); i++) {

        EXPECT_EQ(volume_records_host[i].size(),
                  volume_records_device[i].size());

        for (unsigned int j = 0u; j < volume_records_host[i].size(); j++) {

            EXPECT_EQ(volume_records_host[i][j], volume_records_device[i][j]);

            auto& pos_host = position_records_host[i][j];
            auto& pos_device = position_records_device[i][j];

            EXPECT_NEAR(pos_host[0], pos_device[0], pos_diff_tolerance);
            EXPECT_NEAR(pos_host[1], pos_device[1], pos_diff_tolerance);
            EXPECT_NEAR(pos_host[2], pos_device[2], pos_diff_tolerance);
        }
    }
}
