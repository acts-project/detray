/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/track.hpp"
#include "detray/utils/enumerate.hpp"
#include "tests/common/tools/track_generators.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

constexpr const float epsilon = 1e-5;

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

        // Compare with for loop
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
    for (const auto r :
         uniform_track_generator<detail::ray>(theta_steps, phi_steps)) {
        vector3 &expected = momenta[n_tracks];
        vector3 result = r.dir();

        // Compare with for loop
        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Ray: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);

    // Generate helical trajectories
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    2. * unit_constants::T};
    n_tracks = 0;
    for (const auto track : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps)) {
        detail::helix helix_traj(track, &B);
        vector3 &expected = momenta[n_tracks];
        vector3 result = helix_traj.dir(0.);

        // Compare with for loop
        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Helix: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);
}