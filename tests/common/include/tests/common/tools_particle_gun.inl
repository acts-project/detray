/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

using transform3_type = __plugin::transform3<scalar>;

constexpr const scalar epsilon = 1e-3;

/// Brute force test: Intersect toy geometry and compare between ray and helix
/// without B-field
TEST(tools, particle_gun) {

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto toy_det = create_toy_geometry(host_mr);

    unsigned int theta_steps = 50;
    unsigned int phi_steps = 50;
    const point3 ori{0., 0., 0.};

    // Record ray tracing
    std::vector<std::vector<std::pair<dindex, particle_gun::intersection_type>>>
        expected;
    //  Iterate through uniformly distributed momentum directions with ray
    for (const auto test_ray :
         uniform_track_generator<detail::ray<transform3_type>>(
             theta_steps, phi_steps, ori)) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            particle_gun::shoot_particle(toy_det, test_ray);

        expected.push_back(intersection_record);
    }

    // Simulate straight line track
    const vector3 B{0. * unit<scalar>::T, 0. * unit<scalar>::T,
                    epsilon * unit<scalar>::T};
    // Iterate through uniformly distributed momentum directions with helix
    std::size_t n_tracks{0};
    for (const auto track :
         uniform_track_generator<free_track_parameters<transform3_type>>(
             theta_steps, phi_steps, ori)) {
        const detail::helix test_helix(track, &B);

        // Record all intersections and objects along the ray
        const auto intersection_trace =
            particle_gun::shoot_particle(toy_det, test_helix);

        // Should have encountered the same number of tracks (vulnerable to
        // floating point errors)
        EXPECT_EQ(expected[n_tracks].size(), intersection_trace.size());

        // Check every single recorded intersection
        for (std::size_t i = 0; i < intersection_trace.size(); ++i) {
            if (expected[n_tracks][i].first != intersection_trace[i].first) {
                // Intersection record at portal bound might be flipped
                // (the portals overlap completely)
                if (expected[n_tracks][i].first ==
                        intersection_trace[i + 1].first and
                    expected[n_tracks][i + 1].first ==
                        intersection_trace[i].first) {
                    // Have already checked the next record
                    ++i;
                    continue;
                }
            }
            EXPECT_EQ(expected[n_tracks][i].first, intersection_trace[i].first);
        }

        ++n_tracks;
    }
}