/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/simulation/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/statistics.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3 = __plugin::transform3<detray::scalar>;

TEST(tools, uniform_track_generator) {

    constexpr const scalar tol{1e-5f};

    constexpr std::size_t phi_steps{50u};
    constexpr std::size_t theta_steps{50u};

    std::array<vector3, phi_steps * theta_steps> momenta{};

    // Loops of theta values ]0,pi[
    for (std::size_t itheta{0u}; itheta < theta_steps; ++itheta) {
        const scalar theta{0.01f + static_cast<scalar>(itheta) *
                                       (constant<scalar>::pi - 0.01f) /
                                       static_cast<scalar>(theta_steps)};

        // Loops of phi values [-pi, pi]
        for (std::size_t iphi{0u}; iphi < phi_steps; ++iphi) {
            // The direction
            const scalar phi{-constant<scalar>::pi +
                             static_cast<scalar>(iphi) *
                                 (2.f * constant<scalar>::pi) /
                                 static_cast<scalar>(phi_steps)};

            // intialize a track
            vector3 mom{std::cos(phi) * std::sin(theta),
                        std::sin(phi) * std::sin(theta), std::cos(theta)};
            vector::normalize(mom);
            free_track_parameters<transform3> traj({0.f, 0.f, 0.f}, 0.f, mom,
                                                   -1.f);

            momenta[itheta * phi_steps + iphi] = traj.mom();
        }
    }

    // Now run the track generator and compare
    std::size_t n_tracks{0u};
    for (const auto track :
         uniform_track_generator<free_track_parameters<transform3>>(
             theta_steps, phi_steps)) {
        vector3 &expected = momenta[n_tracks];
        vector3 result = track.mom();

        // Compare with for loop
        EXPECT_NEAR(getter::norm(expected - result), 0.f, tol)
            << "Track: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);

    // Genrate rays
    n_tracks = 0u;
    for (const auto r : uniform_track_generator<detail::ray<transform3>>(
             theta_steps, phi_steps)) {
        vector3 &expected = momenta[n_tracks];
        vector3 result = r.dir();

        // Compare with for loop
        EXPECT_NEAR(getter::norm(expected - result), 0.f, tol)
            << "Ray: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);

    // Generate helical trajectories
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    2.f * unit<scalar>::T};
    n_tracks = 0u;
    for (const auto track :
         uniform_track_generator<free_track_parameters<transform3>>(
             theta_steps, phi_steps)) {
        detail::helix<transform3> helix_traj(track, &B);
        vector3 &expected = momenta[n_tracks];
        vector3 result = helix_traj.dir(0.f);

        // Compare with for loop
        EXPECT_NEAR(getter::norm(expected - result), 0.f, tol)
            << "Helix: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);
}

TEST(tools, uniform_track_generator_with_range) {

    constexpr const scalar tol{1e-5f};

    constexpr std::size_t theta_steps{2u};
    constexpr std::size_t phi_steps{4u};

    std::vector<std::array<scalar, 2>> theta_phi;

    for (const auto track :
         uniform_track_generator<free_track_parameters<transform3>>(
             theta_steps, phi_steps, {0.f, 0.f, 0.f}, 1.f * unit<scalar>::GeV,
             {1.f, 2.f}, {-2.f, 2.f})) {
        const auto dir = track.dir();
        theta_phi.push_back({getter::theta(dir), getter::phi(dir)});
    }

    EXPECT_EQ(theta_phi.size(), 8u);
    EXPECT_NEAR(theta_phi[0][0], 1.f, tol);
    EXPECT_NEAR(theta_phi[0][1], -2.f, tol);
    EXPECT_NEAR(theta_phi[1][0], 1.f, tol);
    EXPECT_NEAR(theta_phi[1][1], -1.f, tol);
    EXPECT_NEAR(theta_phi[2][0], 1.f, tol);
    EXPECT_NEAR(theta_phi[2][1], 0.f, tol);
    EXPECT_NEAR(theta_phi[3][0], 1.f, tol);
    EXPECT_NEAR(theta_phi[3][1], 1.f, tol);
    EXPECT_NEAR(theta_phi[4][0], 1.5f, tol);
    EXPECT_NEAR(theta_phi[4][1], -2.f, tol);
    EXPECT_NEAR(theta_phi[5][0], 1.5f, tol);
    EXPECT_NEAR(theta_phi[5][1], -1.f, tol);
    EXPECT_NEAR(theta_phi[6][0], 1.5f, tol);
    EXPECT_NEAR(theta_phi[6][1], 0.f, tol);
    EXPECT_NEAR(theta_phi[7][0], 1.5f, tol);
    EXPECT_NEAR(theta_phi[7][1], 1.f, tol);
}

/// Tests a random number based track state generator - uniform distribution
TEST(tools, random_track_generator_uniform) {

    // Use deterministic random number generator for testing
    using uniform_gen_t =
        random_numbers<scalar, std::uniform_real_distribution<scalar>,
                       std::seed_seq>;

    // Tolerance depends on sample size
    constexpr scalar tol{0.05f};

    // Track counter
    std::size_t n_tracks{0u};
    constexpr std::size_t n_gen_tracks{10000u};

    // Other params
    point3 ori = {0.f, 0.f, 0.f};
    scalar mom{2.f * unit<scalar>::GeV};
    std::array<scalar, 2> theta_range = {0.01f, constant<scalar>::pi};
    std::array<scalar, 2> phi_range = {-0.9f * constant<scalar>::pi,
                                       0.8f * constant<scalar>::pi};

    // Catch the results
    std::array<scalar, n_gen_tracks> phi{};
    std::array<scalar, n_gen_tracks> theta{};

    for (const auto track :
         random_track_generator<free_track_parameters<transform3>,
                                uniform_gen_t>(n_gen_tracks, ori, mom,
                                               theta_range, phi_range)) {
        phi[n_tracks] = getter::phi(track.dir());
        theta[n_tracks] = getter::theta(track.dir());

        ++n_tracks;
    }

    ASSERT_EQ(n_gen_tracks, n_tracks);

    // Check uniform distrubution

    // Mean
    EXPECT_NEAR(statistics::mean(phi), 0.5f * (phi_range[0] + phi_range[1]),
                tol);
    EXPECT_NEAR(statistics::mean(theta),
                0.5f * (theta_range[0] + theta_range[1]), tol);

    // variance
    EXPECT_NEAR(statistics::variance(phi),
                1.0f / 12.0f * std::pow(phi_range[1] - phi_range[0], 2.f), tol);
    EXPECT_NEAR(statistics::variance(theta),
                1.0f / 12.0f * std::pow(theta_range[1] - theta_range[0], 2.f),
                tol);
}

/// Tests a random number based track state generator - normal distribution
TEST(tools, random_track_generator_normal) {

    // Use deterministic random number generator for testing
    using normal_gen_t =
        random_numbers<scalar, std::normal_distribution<scalar>, std::seed_seq>;

    // Tolerance depends on sample size
    constexpr scalar tol{0.05f};

    // Track counter
    std::size_t n_tracks{0u};
    constexpr std::size_t n_gen_tracks{10000u};

    // Other params
    point3 ori = {0.f, 0.f, 0.f};
    scalar mom{2.f * unit<scalar>::GeV};
    std::array<scalar, 2> theta_range = {0.01f, constant<scalar>::pi};
    std::array<scalar, 2> phi_range = {-0.9f * constant<scalar>::pi,
                                       0.8f * constant<scalar>::pi};

    // Catch the results
    std::array<scalar, n_gen_tracks> phi{};
    std::array<scalar, n_gen_tracks> theta{};

    for (const auto track :
         random_track_generator<free_track_parameters<transform3>,
                                normal_gen_t>(n_gen_tracks, ori, mom,
                                              theta_range, phi_range)) {
        phi[n_tracks] = getter::phi(track.dir());
        theta[n_tracks] = getter::theta(track.dir());

        ++n_tracks;
    }

    ASSERT_EQ(n_gen_tracks, n_tracks);

    // check gaussian distribution - values are clamped to phi/theta range

    // Mean
    EXPECT_NEAR(statistics::mean(phi),
                phi_range[0] + 0.5f * (phi_range[1] - phi_range[0]), tol);
    EXPECT_NEAR(statistics::mean(theta),
                theta_range[0] + 0.5f * (theta_range[1] - theta_range[0]), tol);

    // variance
    EXPECT_NEAR(statistics::variance(phi),
                std::pow(0.5f / 3.0f * (phi_range[1] - phi_range[0]), 2.f),
                tol);
    EXPECT_NEAR(statistics::variance(theta),
                std::pow(0.5f / 3.0f * (theta_range[1] - theta_range[0]), 2.f),
                tol);
}