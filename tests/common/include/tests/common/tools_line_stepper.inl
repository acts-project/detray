/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/track.hpp"
#include "core/transform_store.hpp"
#include "tools/line_stepper.hpp"

#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command

// This tests the base functionality of the line stepper
TEST(__plugin, line_stepper)
{
    using namespace detray;
    using namespace __plugin;

    using detray_track = track<static_transform_store::context>;

    detray_track traj;
    traj.pos = {0., 0., 0.};
    traj.dir = vector::normalize(vector3{1., 1., 0.});
    traj.ctx = static_transform_store::context{};
    traj.momentum = 100.;
    traj.overstep_tolerance = -1e-5;

    line_stepper<detray_track>::state lstate(traj);

    line_stepper<detray_track> lstepper;
    bool heartbeat = lstepper.step(lstate, 10.);
    ASSERT_TRUE(heartbeat);

    ASSERT_FLOAT_EQ(traj.pos[0],  10. / sqrt(2) );
    ASSERT_FLOAT_EQ(traj.pos[1],  10. / sqrt(2) );
    ASSERT_FLOAT_EQ(traj.pos[2], 0.);

    // Step with limit
    lstate.set_limit(20.);
    heartbeat = lstepper.step(lstate, 10.);
    ASSERT_TRUE(heartbeat);

    heartbeat = lstepper.step(lstate, 10.);
    ASSERT_FALSE(heartbeat);


}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

