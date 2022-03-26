/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/transform_store.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/track.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the base functionality of the line stepper
TEST(ALGEBRA_PLUGIN, line_stepper) {
    using namespace detray;
    using namespace __plugin;

    using detray_track = free_track_parameters;

    // dummy navigation struct
    struct nav_state {
        scalar operator()() const { return 1. * unit_constants::mm; }
        inline void set_full_trust() {}
        inline void set_high_trust() {}
        inline void set_fair_trust() {}
        inline void set_no_trust() {}
        inline bool abort() { return false; }
    };

    point3<scalar> pos{0., 0., 0.};
    vector3<scalar> mom{1., 1., 0.};
    detray_track traj(pos, 0, mom, -1);

    line_stepper<detray_track>::state lstate(traj);
    nav_state n_state{};

    line_stepper<detray_track> lstepper;
    ASSERT_TRUE(lstepper.step(lstate, n_state));

    ASSERT_FLOAT_EQ(traj.pos()[0], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[1], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[2], 0.);

    // Only step 1.5 mm further, then break
    lstate.set_path_limit(1.5);
    ASSERT_TRUE(lstepper.step(lstate, n_state));

    ASSERT_FLOAT_EQ(traj.pos()[0], std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[1], std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[2], 0.);

    // breaks
    ASSERT_FALSE(lstepper.step(lstate, n_state));
}
