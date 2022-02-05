/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/track.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/tools/line_stepper.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the base functionality of the line stepper
TEST(ALGEBRA_PLUGIN, line_stepper) {
    using namespace detray;
    using namespace __plugin;

    using detray_track = track<static_transform_store<>::context>;

    detray_track traj;
    traj.pos = {0., 0., 0.};
    traj.dir = vector::normalize(vector3<detray::scalar>{1., 1., 0.});
    traj.ctx = static_transform_store<>::context{};
    traj.momentum = 100.;
    traj.overstep_tolerance = -1e-5;

    line_stepper<detray_track>::state lstate(traj);

    line_stepper<detray_track> lstepper;
    bool heartbeat = lstepper.step(lstate, 10.);
    ASSERT_TRUE(heartbeat);

    ASSERT_FLOAT_EQ(traj.pos[0], 10. / sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos[1], 10. / sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos[2], 0.);

    // Step with limit
    lstate.set_limit(20.);
    heartbeat = lstepper.step(lstate, 10.);
    ASSERT_TRUE(heartbeat);

    heartbeat = lstepper.step(lstate, 10.);
    ASSERT_FALSE(heartbeat);
}
