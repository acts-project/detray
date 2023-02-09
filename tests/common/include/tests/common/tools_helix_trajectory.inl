/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/tracks/tracks.hpp"

// System include(s)
#include <cmath>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

constexpr const scalar tol{1e-5f};

// This tests the base functionality of the Helix Gun
TEST(tools, helix_trajectory) {

    using vector3 = __plugin::vector3<scalar>;
    using point3 = __plugin::point3<scalar>;
    using transform3_type = __plugin::transform3<scalar>;

    const point3 pos{0.f, 0.f, 0.f};
    const scalar time{0.f};
    const vector3 mom{1.f, 0.f, 1.f * unit<scalar>::GeV};
    const scalar q{-1.f * unit<scalar>::e};

    // vertex
    free_track_parameters<transform3_type> vertex(pos, time, mom, q);

    // magnetic field
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};

    const scalar p_mag{getter::norm(mom)};
    const scalar B_mag{getter::norm(B)};
    const scalar pz_along{vector::dot(mom, vector::normalize(B))};
    const scalar pt{std::sqrt(p_mag * p_mag - pz_along * pz_along)};

    // helix trajectory
    detail::helix helix_traj(vertex, &B);

    // radius of helix
    scalar R{helix_traj.radius()};
    EXPECT_NEAR(R, pt / B_mag, tol);

    // After half turn
    point3 half_turn = helix_traj(p_mag / B_mag * constant<scalar>::pi);

    EXPECT_NEAR(half_turn[0], 0.f, R * tol);
    EXPECT_NEAR(half_turn[1], 2.f * R, R * tol);
    EXPECT_NEAR(half_turn[2], pz_along / B_mag * constant<scalar>::pi, R * tol);

    // After one full turn
    point3 full_turn = helix_traj(2.f * p_mag / B_mag * constant<scalar>::pi);

    EXPECT_NEAR(full_turn[0], 0.f, R * tol);
    EXPECT_NEAR(full_turn[1], 0.f, R * tol);
    EXPECT_NEAR(full_turn[2], 2.f * pz_along / B_mag * constant<scalar>::pi,
                R * tol);
}

TEST(tools, helix_trajectory_small_pT) {

    using vector3 = __plugin::vector3<scalar>;
    using point3 = __plugin::point3<scalar>;
    using transform3_type = __plugin::transform3<scalar>;

    const point3 pos{0.f, 0.f, 0.f};
    const scalar time{0.f};
    const vector3 mom{0.f, 0.f, 1.f * unit<scalar>::GeV};
    const scalar q{-1. * unit<scalar>::e};

    // vertex
    free_track_parameters<transform3_type> vertex(pos, time, mom, q);

    // magnetic field
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};

    // helix trajectory
    detail::helix helix_traj(vertex, &B);

    // After 10 mm
    const scalar path_length{10.f * unit<scalar>::mm};
    const point3 helix_pos = helix_traj(path_length);
    const point3 true_pos = pos + path_length * vector::normalize(mom);

    EXPECT_NEAR(true_pos[0], helix_pos[0], tol);
    EXPECT_NEAR(true_pos[1], helix_pos[1], tol);
    EXPECT_NEAR(true_pos[2], helix_pos[2], tol);
}