/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/tools/rk_stepper.hpp"
#include "detray/tools/track.hpp"
#include "tests/common/tools/helix_gun.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using namespace __plugin;

constexpr scalar epsilon = 1e-4;

// This tests the base functionality of the Runge-Kutta stepper
TEST(ALGEBRA_PLUGIN, rk_stepper) {

    // type definitions
    using vector3 = vector3<scalar>;
    using point3 = point3<scalar>;

    using mag_field_type = constant_magnetic_field<>;
    using rk_stepper_type = rk_stepper<mag_field_type, free_track_parameters>;

    // RK stepper configurations
    constexpr unsigned int theta_steps = 100;
    constexpr unsigned int phi_steps = 100;
    constexpr unsigned int rk_steps = 100;

    // Constant magnetic field
    vector3 B{1 * unit_constants::T, 1 * unit_constants::T,
              1 * unit_constants::T};
    mag_field_type mag_field(B);

    // RK stepper
    rk_stepper_type rk_stepper(mag_field);

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};
    scalar p_mag = 10;

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {

            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            const vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                              cos_theta};

            vector3 mom = p_mag * dir;
            free_track_parameters traj(ori, 0, mom, -1);

            // helix gun
            helix_gun helix(traj, B);

            // RK Stepping into forward direction
            rk_stepper_type::state forward_state(traj);
            for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
                rk_stepper.step(forward_state);
            }

            // get relative error by dividing error with path length
            auto path_accumulated = forward_state._path_accumulated;
            auto helix_pos = helix(path_accumulated);
            auto forward_pos = forward_state().pos();
            auto forward_relative_error =
                1. / path_accumulated * (forward_pos - helix_pos);

            // Make sure that relative error is smaller than epsion
            EXPECT_NEAR(getter::norm(forward_relative_error), 0, epsilon);

            // RK Stepping into backward direction
            traj.flip();
            rk_stepper_type::state backward_state(traj);
            for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
                rk_stepper.step(backward_state);
            }

            // get relative error by dividing error with path length
            path_accumulated += backward_state._path_accumulated;
            auto backward_pos = backward_state().pos();
            auto backward_relative_error =
                1. / path_accumulated * (backward_pos - ori);

            // Make sure that relative error is smaller than epsion
            EXPECT_NEAR(getter::norm(backward_relative_error), 0, epsilon);
        }
    }
}
