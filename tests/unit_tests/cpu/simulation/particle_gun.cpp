/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/test/utils/particle_gun.hpp"

#include "detray/definitions/units.hpp"
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using algebra_t = test::algebra;
using vector3 = test::vector3;

constexpr const scalar tol{1e-3f};

/// Brute force test: Intersect toy geometry and compare between ray and helix
/// without B-field
GTEST_TEST(detray_simulation, particle_gun) {

    // Simulate straight line track
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    tol * unit<scalar>::T};

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto [toy_det, names] = build_toy_detector(host_mr);

    unsigned int theta_steps{50u};
    unsigned int phi_steps{50u};

    // Record ray tracing
    using detector_t = decltype(toy_det);
    using algebra_type = typename detector_t::algebra_type;
    using intersection_t =
        intersection2D<typename detector_t::surface_type, algebra_type>;
    std::vector<std::vector<std::pair<dindex, intersection_t>>> expected;
    //  Iterate through uniformly distributed momentum directions with ray
    for (const auto test_ray : uniform_track_generator<detail::ray<algebra_t>>(
             phi_steps, theta_steps)) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            particle_gun::shoot_particle(toy_det, test_ray);

        expected.push_back(intersection_record);
    }

    // Iterate through uniformly distributed momentum directions with helix
    std::size_t n_tracks{0u};
    for (const auto track :
         uniform_track_generator<free_track_parameters<algebra_t>>(
             phi_steps, theta_steps)) {
        const detail::helix test_helix(track, &B);

        // Record all intersections and objects along the ray
        const auto intersection_trace =
            particle_gun::shoot_particle(toy_det, test_helix);

        // Should have encountered the same number of tracks (vulnerable to
        // floating point errors)
        EXPECT_EQ(expected[n_tracks].size(), intersection_trace.size());

        // Check every single recorded intersection
        for (std::size_t i = 0u; i < intersection_trace.size(); ++i) {
            if (expected[n_tracks][i].first != intersection_trace[i].first) {
                // Intersection record at portal bound might be flipped
                // (the portals overlap completely)
                if (expected[n_tracks][i].first ==
                        intersection_trace[i + 1u].first and
                    expected[n_tracks][i + 1u].first ==
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
