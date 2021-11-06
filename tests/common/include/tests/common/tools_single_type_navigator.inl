/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/mask_store.hpp"
#include "detray/core/track.hpp"
#include "detray/tools/single_type_navigator.hpp"
#include "tests/common/read_geometry.hpp"

/// @note __plugin has to be defined with a preprocessor command

auto [volumes, surfaces, transforms, discs, cylinders, rectangles] =
    toy_geometry();

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, single_type_navigator) {
    using namespace detray;

    /** Tolerance for tests */
    constexpr double tol = 0.01;

    /** Empty context type struct */
    struct empty_context {};

    mask_store<dtuple, dvector, decltype(discs)::value_type,
               decltype(cylinders)::value_type,
               decltype(rectangles)::value_type>
        masks;
    // populate mask store
    masks.add_masks(discs);
    masks.add_masks(cylinders);
    masks.add_masks(rectangles);

    single_type_navigator n(volumes, surfaces, transforms, masks);
    using toy_navigator = decltype(n);

    // test track
    track<empty_context> traj;
    traj.pos = {0., 0., 0.};
    traj.dir = vector::normalize(vector3{1., 1., 0.});
    traj.ctx = empty_context{};
    traj.momentum = 100.;
    traj.overstep_tolerance = -1e-4;

    toy_navigator::state state;
    state.set_volume(0u);

    // Check that the state is unitialized
    // Volume is invalid
    ASSERT_EQ(state.volume(), 0u);
    // No surface candidates
    ASSERT_EQ(state.candidates().size(), 0u);
    // You can not trust the state
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_no_trust);
    // The status is unkown
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_unknown);

    //
    // beampipe
    //

    // Initial status call
    bool heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is towards portal (beampipe)
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    // Now the volume, surfaces are defined and are trustworthy
    ASSERT_EQ(state.volume(), 0u);
    // Only the cylinder portal
    ASSERT_EQ(state.candidates().size(), 1u);
    ASSERT_EQ(state.next()->index, 2u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);
    ASSERT_NEAR(state(), 27., tol);

    // Let's immediately target, nothing should change, as there is full trust
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status remains: towards portal
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 0u);
    ASSERT_EQ(state.candidates().size(), 1u);
    ASSERT_EQ(state.next()->index, 2u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);
    ASSERT_NEAR(state(), 27., tol);

    // Let's make half the step towards the portal
    traj.pos = traj.pos + 0.5 * state() * traj.dir;
    // Could be externally set by actor (in the future)
    state.set_trust_level(toy_navigator::navigation_trust_level::e_high_trust);
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status remains: towards portal
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 0u);
    ASSERT_EQ(state.candidates().size(), 1u);
    ASSERT_EQ(state.next()->index, 2u);
    ASSERT_NEAR(state(), 13.5, tol);
    // Trust level is restored
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Let's immediately target, nothing should change, as there is full trust
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status remains: towards surface
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 0u);
    ASSERT_EQ(state.candidates().size(), 1u);
    ASSERT_EQ(state.next()->index, 2u);
    ASSERT_NEAR(state(), 13.5, tol);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Now step onto the portal (2)
    traj.pos = traj.pos + state() * traj.dir;
    state.set_trust_level(toy_navigator::navigation_trust_level::e_high_trust);
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on portal
    ASSERT_TRUE(std::abs(state()) < state.tolerance());
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    // Switched to next volume
    ASSERT_EQ(state.volume(), 1u);
    // Kernel is exhaused, and trust level is gone
    ASSERT_EQ(state.next(), state.candidates().end());
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_no_trust);

    //
    // layer 1
    //

    // New target call will initialize volume 1
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on adjacent portal in volume 1, towards next candidate
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 1u);
    // This includes overlapping modules and the adjacent portal we are
    // already on
    ASSERT_EQ(state.candidates().size(), 6u);
    // We are already on this portal, so switch to next candidate which must
    // be a surface
    ASSERT_EQ(state.next()->index, 128u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Now step onto the surface in volume 1 (128)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 128
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 1u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // points to the next surface now
    ASSERT_EQ(state.next()->index, 129u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target - update distance 129
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 1u);
    // Should be on our way to the next ovelapping module
    ASSERT_EQ(state.candidates().size(), 6u);
    // this is still the next surface, since we did not step
    ASSERT_EQ(state.next()->index, 129u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Jump to the next surface in volume 1 (129)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 129
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 1u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // points to the next surface now
    ASSERT_EQ(state.next()->index, 112u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target - update distance to 112
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 1u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // this is still the next surface, since we did not step
    ASSERT_EQ(state.next()->index, 112u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Jump to the next surface in volume 1 (112)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 112
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 1u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // points to the next surface now
    ASSERT_EQ(state.next()->index, 113u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target - update distance to 113
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 1u);
    // Should be on our way to the next ovelapping module
    ASSERT_EQ(state.candidates().size(), 6u);
    // this is still the next surface, since we did not step
    ASSERT_EQ(state.next()->index, 113u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Jump to the next surface (113)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 113
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 1u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // points to the portal towards the gap volume now
    ASSERT_EQ(state.next()->index, 6u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target again - should go towards portal 6 next
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 1u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // the portal is still the next object, since we did not step
    ASSERT_EQ(state.next()->index, 6u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Step onto the portal in volume 1
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is on portal 6
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    // Switch to volume 2
    ASSERT_EQ(state.volume(), 2u);
    // Kernel is exhaused, and trust level is gone
    ASSERT_EQ(state.next(), state.candidates().end());
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_no_trust);

    //
    // gap volume
    //

    // With the new target call all surfaces of vol.2 should be initialized
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.volume(), 2u);
    // The status is: on adjacent portal in volume 2, towards next candidate,
    // which is portal 234
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    // This includes the adjacent portal we are already on
    ASSERT_EQ(state.candidates().size(), 2u);
    ASSERT_EQ(state.next()->index, 234u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Step onto the portal 234
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on portal
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    // Switch to volume 3
    ASSERT_EQ(state.volume(), 3u);
    // Kernel is exhaused, and trust level is gone
    ASSERT_EQ(state.next(), state.candidates().end());
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_no_trust);

    //
    // layer 2
    //

    // With the new target call all objects of vol.3 should be initialized
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.volume(), 3u);
    // The status is: on adjacent portal in volume 3, towards next candidate,
    // which should be a module surface
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    // This includes the adjacent portal we are already on
    ASSERT_EQ(state.candidates().size(), 6u);
    // We are already on this portal, so switch to next candidate which must
    // be a surface
    ASSERT_EQ(state.next()->index, 482u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Now step onto the surface (482)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 482
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 3u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // Next module surface from 482
    ASSERT_EQ(state.next()->index, 450u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target - update distance to 450
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 3u);
    // Should be on our way to the next ovelapping module
    ASSERT_EQ(state.candidates().size(), 6u);
    ASSERT_EQ(state.next()->index, 450u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Jump to the next surface (450)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 450
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 3u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // Next module surface from 450
    ASSERT_EQ(state.next()->index, 483u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target - update distance to 483
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 3u);
    // Should be on our way to the next ovelapping module
    ASSERT_EQ(state.candidates().size(), 6u);
    // Next module surface from 450
    ASSERT_EQ(state.next()->index, 483u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Jump to the next surface (483)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 483
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 3u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // Next module surface from 483
    ASSERT_EQ(state.next()->index, 451u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target - update distance to 451
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), 3u);
    // Should be on our way to the next ovelapping module
    ASSERT_EQ(state.candidates().size(), 6u);
    // Next module surface from 483
    ASSERT_EQ(state.next()->index, 451u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Jump to the next surface (451)
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on surface 451
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_object);
    ASSERT_EQ(state.volume(), 3u);
    ASSERT_EQ(state.candidates().size(), 6u);
    // Next is the portal that leaves the detector world (238)
    ASSERT_EQ(state.next()->index, 238u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_high_trust);

    // Let's target again - should go towards portal 238 next
    heartbeat = n.target(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    ASSERT_EQ(state.volume(), 3u);
    ASSERT_EQ(state.status(),
              toy_navigator::navigation_status::e_towards_object);
    ASSERT_EQ(state.candidates().size(), 6u);
    ASSERT_EQ(state.next()->index, 238u);
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Step onto the portal
    traj.pos = traj.pos + state() * traj.dir;
    heartbeat = n.status(state, traj);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is: on portal
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_on_target);
    // Switch to next volume leads out of the detector world -> exit
    ASSERT_EQ(state.volume(), dindex_invalid);
    // We know we went out of the detector
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
