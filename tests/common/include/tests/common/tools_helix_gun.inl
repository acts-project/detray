/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/tools/track.hpp"
#include "tests/common/tools/helix_gun.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using namespace __plugin;

constexpr const float epsilon = 1e-6;

// This tests the base functionality of the Helix Gun
TEST(tools, helix_gun) {

    using vector3 = vector3<scalar>;
    using point3 = point3<scalar>;

    point3 pos{0., 0., 0.};
    scalar time = 0.;
    vector3 mom{1., 0., 1.};
    scalar q = -1.;

    // vertex
    free_track_parameters vertex(pos, time, mom, q);

    // magnetic field
    vector3 B{0, 0, 1 * unit_constants::T};

    scalar p_mag = getter::norm(mom);
    scalar B_mag = getter::norm(B);
    scalar pz_along = vector::dot(mom, vector::normalize(B));
    scalar pt = std::sqrt(std::pow(p_mag, 2) - std::pow(pz_along, 2));

    // helix gun
    helix_gun helix(vertex, B);

    // radius of helix
    scalar R = helix.radius();
    EXPECT_FLOAT_EQ(R, pt / B_mag);

    // After half turn
    point3 half_turn = helix(p_mag / B_mag * M_PI);

    EXPECT_NEAR(half_turn[0], 0., R * epsilon);
    EXPECT_NEAR(half_turn[1], 2 * R, R * epsilon);
    EXPECT_NEAR(half_turn[2], pz_along / B_mag * M_PI, R * epsilon);

    // After one full turn
    point3 full_turn = helix(2 * p_mag / B_mag * M_PI);

    EXPECT_NEAR(full_turn[0], 0., R * epsilon);
    EXPECT_NEAR(full_turn[1], 0., R * epsilon);
    EXPECT_NEAR(full_turn[2], 2 * pz_along / B_mag * M_PI, R * epsilon);
}

TEST(tools, helix_gun_small_pT) {

    using vector3 = vector3<scalar>;
    using point3 = point3<scalar>;

    point3 pos{0., 0., 0.};
    scalar time = 0.;
    vector3 mom{0., 0., 1.};
    scalar q = -1.;

    // vertex
    free_track_parameters vertex(pos, time, mom, q);

    // magnetic field
    vector3 B{0, 0, 1 * unit_constants::T};

    // helix gun
    helix_gun helix(vertex, B);

    // After 10 mm
    scalar path_length = 10;
    point3 helix_pos = helix(path_length);
    point3 true_pos = pos + path_length * vector::normalize(mom);

    EXPECT_EQ(true_pos, helix_pos);
}