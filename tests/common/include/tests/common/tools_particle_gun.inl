/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/particle_gun.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using namespace __plugin;

constexpr const float epsilon = 1e-5;

// This tests the base functionality of the Helix Gun
TEST(tools, helix_trajectory) {

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

    // helix trajectory
    helix helix_traj(vertex, &B);

    // radius of helix
    scalar R = helix_traj.radius();
    EXPECT_FLOAT_EQ(R, pt / B_mag);

    // After half turn
    point3 half_turn = helix_traj(p_mag / B_mag * M_PI);

    EXPECT_NEAR(half_turn[0], 0., R * epsilon);
    EXPECT_NEAR(half_turn[1], 2 * R, R * epsilon);
    EXPECT_NEAR(half_turn[2], pz_along / B_mag * M_PI, R * epsilon);

    // After one full turn
    point3 full_turn = helix_traj(2 * p_mag / B_mag * M_PI);

    EXPECT_NEAR(full_turn[0], 0., R * epsilon);
    EXPECT_NEAR(full_turn[1], 0., R * epsilon);
    EXPECT_NEAR(full_turn[2], 2 * pz_along / B_mag * M_PI, R * epsilon);
}

TEST(tools, helix_trajectory_small_pT) {

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

    // helix trajectory
    helix helix_traj(vertex, &B);

    // After 10 mm
    scalar path_length = 10;
    point3 helix_pos = helix_traj(path_length);
    point3 true_pos = pos + path_length * vector::normalize(mom);

    EXPECT_FLOAT_EQ(true_pos[0], helix_pos[0]);
    EXPECT_FLOAT_EQ(true_pos[1], helix_pos[1]);
    EXPECT_FLOAT_EQ(true_pos[2], helix_pos[2]);
}

TEST(tools, uniform_track_generator) {

    using vector3 = __plugin::vector3<scalar>;

    constexpr std::size_t phi_steps = 5;
    constexpr std::size_t theta_steps = 5;

    std::array<vector3, phi_steps * theta_steps> momenta{};

    // Loops of theta values ]0,pi[
    for (std::size_t itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.01 + itheta * (M_PI - 0.01) / theta_steps;

        // Loops of phi values [-pi, pi]
        for (std::size_t iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2. * M_PI) / phi_steps;

            // intialize a track
            vector3 mom{std::cos(phi) * std::sin(theta),
                        std::sin(phi) * std::sin(theta), std::cos(theta)};
            vector::normalize(mom);
            free_track_parameters traj({0., 0., 0.}, 0, mom, -1);

            momenta[itheta * phi_steps + iphi] = traj.mom();
        }
    }

    // Now run the track generator and compare
    std::size_t n_tracks = 0;
    for (const auto track : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps)) {
        vector3 &expected = momenta[n_tracks];
        vector3 result = track.mom();

        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Track: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);

    // Genrate rays
    n_tracks = 0;
    for (const auto r : uniform_track_generator<ray>(theta_steps, phi_steps)) {
        vector3 &expected = momenta[n_tracks];
        vector3 result = r.dir();

        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Ray: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);

    // Generate helix trajectories
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    2. * unit_constants::T};
    n_tracks = 0;
    for (const auto track : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps)) {
        helix helix_traj(track, &B);
        vector3 &expected = momenta[n_tracks];
        vector3 result = helix_traj.dir(0.);

        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Helix: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);
}