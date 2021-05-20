/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "core/track.hpp"
#include "io/csv_io.hpp"
#include "tools/navigator.hpp"

#include <iostream>
#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction and general methods of the navigator
TEST(__plugin, navigator)
{
    using namespace detray;

    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr)
    {
        throw std::ios_base::failure("Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string surface_file = data_directory + std::string("tml.csv");
    std::string surface_grid_file = data_directory + std::string("tml-surface-grids.csv");
    std::string layer_volume_file = data_directory + std::string("tml-layer-volumes.csv");

    auto d = detector_from_csv<static_transform_store>("tml", surface_file, surface_grid_file, layer_volume_file);

    // Create the navigator
    using detray_navigator = navigator<decltype(d)>;

    detray_navigator n(std::move(d));

    track<static_transform_store::context> traj;
    traj.pos = { 0., 0., 0. };
    traj.dir = vector::normalize(vector3{ 1., 1., 0. });
    traj.ctx = static_transform_store::context{};
    traj.momentum = 100.;


    detray_navigator::navigation_state state;

    // Check that the state is unitialized
    // Volume is invalid
    ASSERT_EQ(state.volume_index, dindex_invalid);
    // No surface candidates
    ASSERT_EQ(state.surface_kernel.candidates.size(), 0u);
    // No portal candidates
    ASSERT_EQ(state.portal_kernel.candidates.size(), 0u);
    // You can not trust the state
    ASSERT_EQ(state.trust_level, detray_navigator::navigation_trust_level::e_not);

    // Initial status call
    n.status(state, traj);
    // Now the volume, surfaces are defined and are trust worthy
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 1u);
    ASSERT_EQ(state.trust_level, detray_navigator::navigation_trust_level::e_full);
    ASSERT_TRUE(std::abs(state.distance_to_next-19.) < 0.01);
    // Still no portals defined
    ASSERT_EQ(state.portal_kernel.candidates.size(), 0u);

    // Let's immediately target, nothing should change, as there is full trust
    n.target(state, traj);
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 1u);
    ASSERT_EQ(state.trust_level, detray_navigator::navigation_trust_level::e_full);
    ASSERT_TRUE(std::abs(state.distance_to_next-19.) < 0.01);

    // Let's make half the step towards the surface
    traj.pos = 0.5 * state.distance_to_next * traj.dir;
    // In this case, this would set the trust level to 'high'
    state.trust_level = detray_navigator::navigation_trust_level::e_high;

}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
