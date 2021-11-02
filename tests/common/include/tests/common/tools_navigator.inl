/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/detector.hpp"
#include "detray/core/track.hpp"
#include "detray/core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "tests/common/read_geometry.hpp"
#include "detray/tools/navigator.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, navigator) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;

    auto [d, name_map] = read_from_csv(tml_files, host_mr);

    // Create the navigator
    using detray_navigator = navigator<decltype(d)>;
    using detray_context = decltype(d)::transform_store::context;

    detray_navigator n(std::move(d));

    track<detray_context> traj;
    traj.pos = {0., 0., 0.};
    traj.dir = vector::normalize(vector3{1., 1., 0.});
    traj.ctx = detray_context{};
    traj.momentum = 100.;
    traj.overstep_tolerance = -1e-4;

    detray_navigator::state state;

    // Check that the state is unitialized
    // Volume is invalid
    ASSERT_EQ(state.volume_index, dindex_invalid);
    // No surface candidates
    ASSERT_EQ(state.surface_kernel.candidates.size(), 0u);
    // No portal candidates
    ASSERT_EQ(state.portal_kernel.candidates.size(), 0u);
    // You can not trust the state
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_no_trust);
    // The status is unkown
    ASSERT_EQ(state.status, detray_navigator::navigation_status::e_unknown);

    // Initial status call
    bool heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is towards surface
    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_surface);
    // Now the volume, surfaces are defined and are trust worthy
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 1u);
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_full_trust);
    ASSERT_TRUE(std::abs(state() - 19.) < 0.01);
    // Still no portals defined
    ASSERT_EQ(state.portal_kernel.candidates.size(), 0u);

    // Let's immediately target, nothing should change, as there is full trust
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status remains: towards surface
    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_surface);
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 1u);
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_full_trust);
    ASSERT_TRUE(std::abs(state() - 19.) < 0.01);

    // Let's make half the step towards the surface
    traj.pos = traj.pos + 0.5 * state() * traj.dir;
    // In this case, this would set the trust level to 'high'
    state.trust_level = detray_navigator::navigation_trust_level::e_high_trust;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status remains: towards surface
    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_surface);
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 1u);
    ASSERT_TRUE(std::abs(state() - 9.5) < 0.01);
    // Trust level is restored
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_full_trust);

    // Let's immediately target, nothing should change, as there is full trust
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status remains: towards surface
    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_surface);
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 1u);
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_full_trust);
    ASSERT_TRUE(std::abs(state() - 9.5) < 0.01);

    // Now step onto the surface
    traj.pos = traj.pos + state() * traj.dir;
    state.trust_level = detray_navigator::navigation_trust_level::e_high_trust;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    // The status is: on surface
    ASSERT_EQ(state.status, detray_navigator::navigation_status::e_on_surface);
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_TRUE(std::abs(state()) < state.on_surface_tolerance);
    ASSERT_EQ(state.status, detray_navigator::navigation_status::e_on_surface);
    // Kernel is exhaused, and trust level is gone
    ASSERT_EQ(state.surface_kernel.next, state.surface_kernel.candidates.end());
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_high_trust);

    // New target call
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    // The status is: towards portal
    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_portal);
    ASSERT_EQ(state.volume_index, 0u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 0u);
    ASSERT_EQ(state.portal_kernel.candidates.size(), 1u);
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_full_trust);

    // Now step towards the portal
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    // The status is: on portal - points towards volume 16
    ASSERT_EQ(state.status, detray_navigator::navigation_status::e_on_portal);
    ASSERT_EQ(state.volume_index, 16u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 0u);
    ASSERT_EQ(state.portal_kernel.candidates.size(), 0u);
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_no_trust);

    // Let's target now - new volume should be volume 16 and is empty
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_portal);
    ASSERT_EQ(state.volume_index, 16u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 0u);

    // Jump to the next portal
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    // The status is: on portal - points towards volume 17
    ASSERT_EQ(state.status, detray_navigator::navigation_status::e_on_portal);
    ASSERT_EQ(state.volume_index, 17u);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 0u);
    ASSERT_EQ(state.portal_kernel.candidates.size(), 0u);
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_no_trust);

    // Let's target now - new volume should be volume 17 and should not be empty
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    ASSERT_EQ(state.volume_index, 17u);

    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_surface);
    ASSERT_EQ(state.surface_kernel.candidates.size(), 4u);
    ASSERT_EQ(std::distance(state.surface_kernel.next,
                            state.surface_kernel.candidates.end()),
              4u);
    ASSERT_EQ(state.trust_level,
              detray_navigator::navigation_trust_level::e_full_trust);

    // Intersect the remaining ones
    for (unsigned int is = 0; is < 3; ++is) {
        // Step towards the surface
        traj.pos = traj.pos + state() * traj.dir;
        heartbeat = n.status(state, traj);
        // Test that the navigator has a heartbeat
        ASSERT_TRUE(heartbeat);

        // The status is: on portal - points towards volume 17
        ASSERT_EQ(state.status,
                  detray_navigator::navigation_status::e_on_surface);
        ASSERT_EQ(state.volume_index, 17u);
        // We should have switched by one
        ASSERT_EQ(std::distance(state.surface_kernel.next,
                                state.surface_kernel.candidates.end()),
                  3u - is);
        ASSERT_EQ(state.trust_level,
                  detray_navigator::navigation_trust_level::e_high_trust);

        heartbeat = n.target(state, traj);
        // Test that the navigator has a heartbeat
        ASSERT_TRUE(heartbeat);

        ASSERT_EQ(state.status,
                  detray_navigator::navigation_status::e_towards_surface);
        ASSERT_EQ(state.volume_index, 17u);
        ASSERT_EQ(state.surface_kernel.candidates.size(), 4u);
        ASSERT_EQ(std::distance(state.surface_kernel.next,
                                state.surface_kernel.candidates.end()),
                  3u - is);
        ASSERT_EQ(state.trust_level,
                  detray_navigator::navigation_trust_level::e_full_trust);
    }

    // Surface kernel is now exhausted, status call should invalidate
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    ASSERT_EQ(state.status, detray_navigator::navigation_status::e_on_surface);
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);

    ASSERT_EQ(state.status,
              detray_navigator::navigation_status::e_towards_portal);

    // Let's try to see if the heartbeat dies off at the end of world
    traj.pos = point3{1011, 0., 1355};
    detray_navigator::state late_state;
    heartbeat = n.status(late_state, traj);
    // Heartbeat should be alive
    ASSERT_TRUE(heartbeat);
    // We should be in volume 96
    ASSERT_EQ(late_state.volume_index, 96u);
    ASSERT_EQ(late_state.surface_kernel.candidates.size(), 0u);
    heartbeat = n.target(late_state, traj);
    // Heartbeat should be alive
    ASSERT_EQ(late_state.volume_index, 96u);
    ASSERT_EQ(late_state.surface_kernel.candidates.size(), 0u);
    ASSERT_EQ(late_state.portal_kernel.candidates.size(), 1u);
    // Step to the last portal
    traj.pos = traj.pos + late_state() * traj.dir;
    heartbeat = n.status(late_state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_EQ(late_state.status,
              detray_navigator::navigation_status::e_on_portal);
    ASSERT_TRUE(heartbeat);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
