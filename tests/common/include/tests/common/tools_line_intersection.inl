/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/ray_line_intersector.hpp"
#include "detray/masks/line.hpp"
#include "detray/propagator/track.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <climits>
#include <cmath>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using cartesian = __plugin::cartesian2<detray::scalar>;
using transform3 = __plugin::transform3<detray::scalar>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using point2 = __plugin::point2<scalar>;
constexpr scalar tolerance = 1e-5;

// This defines the local frame test suite
TEST(tools, line_intersector) {

    // tf3 with Identity rotation and no translation
    vector3 x{1, 0, 0};
    vector3 y{0, 1, 0};
    vector3 z{0, 0, 1};
    vector3 t{0, 0, 0};

    const transform3 tf{t, x, y, z};

    /// Create a track
    std::vector<free_track_parameters> trks(2);
    trks[0] = free_track_parameters({1, -1, 0}, 0, {0, 1, 0}, -1);
    trks[1] = free_track_parameters({-1, -1, 0}, 0, {0, 1, 0}, -1);

    // Infinite wire with 10 mm radial cell size
    line<> ln{10., std::numeric_limits<scalar>::infinity(), 0u};

    // Test intersect
    std::vector<ray_line_intersector::intersection_type> is(2);
    is[0] = ray_line_intersector().intersect(tf, trks[0], ln);
    is[1] = ray_line_intersector().intersect(tf, trks[1], ln);

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_EQ(is[0].path, 1);
    EXPECT_EQ(is[0].p3, point3({1, 0, 0}));
    EXPECT_EQ(is[0].p2, point2({1, 0}));

    EXPECT_EQ(is[1].status, intersection::status::e_inside);
    EXPECT_EQ(is[1].path, 1);
    EXPECT_EQ(is[1].p3, point3({-1, 0, 0}));
    EXPECT_EQ(is[1].p2, point2({1, 0}));
}

TEST(tools, line_intersector_square_scope) {

    // tf3 with Identity rotation and no translation
    vector3 x{1, 0, 0};
    vector3 y{0, 1, 0};
    vector3 z{0, 0, 1};
    vector3 t{0, 0, 0};

    const transform3 tf{t, x, y, z};

    /// Create a track
    std::vector<free_track_parameters> trks;
    trks.push_back({{2, 0, 0}, 0, {-1, 1, 0}, -1});
    trks.push_back({{1.9, 0, 0}, 0, {-1, 1, 0}, -1});
    trks.push_back({{2.1, 0, 0}, 0, {-1, 1, 0}, -1});

    trks.push_back({{-2, 0, 0}, 0, {1, 1, 0}, -1});
    trks.push_back({{-1.9, 0, 0}, 0, {1, 1, 0}, -1});
    trks.push_back({{-2.1, 0, 0}, 0, {1, 1, 0}, -1});

    trks.push_back({{0, -2, 0}, 0, {-1, 1, 0}, -1});
    trks.push_back({{0, -1.9, 0}, 0, {-1, 1, 0}, -1});
    trks.push_back({{0, -2.1, 0}, 0, {-1, 1, 0}, -1});

    trks.push_back({{0, -2, 0}, 0, {1, 1, 0}, -1});
    trks.push_back({{0, -1.9, 0}, 0, {1, 1, 0}, -1});
    trks.push_back({{0, -2.1, 0}, 0, {1, 1, 0}, -1});

    // Infinite wire with 1 mm square cell size
    line<ray_line_intersector, cartesian, dindex, true> ln{
        1., std::numeric_limits<scalar>::infinity(), 0u};

    // Test intersect
    std::vector<ray_line_intersector::intersection_type> is;
    for (const auto& trk : trks) {
        is.push_back(
            ray_line_intersector().intersect(tf, trk, ln, {1e-5, 1e-5}));
    }

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_NEAR(is[0].path, sqrt(2), tolerance);
    EXPECT_NEAR(is[0].p3[0], 1, tolerance);
    EXPECT_NEAR(is[0].p3[1], 1, tolerance);
    EXPECT_NEAR(is[0].p3[2], 0, tolerance);
    EXPECT_NEAR(is[0].p2[0], sqrt(2), tolerance);
    EXPECT_NEAR(is[0].p2[1], 0, tolerance);

    EXPECT_EQ(is[1].status, intersection::status::e_inside);
    EXPECT_EQ(is[2].status, intersection::status::e_outside);

    EXPECT_EQ(is[3].status, intersection::status::e_inside);
    EXPECT_EQ(is[4].status, intersection::status::e_inside);
    EXPECT_EQ(is[5].status, intersection::status::e_outside);

    EXPECT_EQ(is[6].status, intersection::status::e_inside);
    EXPECT_EQ(is[7].status, intersection::status::e_inside);
    EXPECT_EQ(is[8].status, intersection::status::e_outside);

    EXPECT_EQ(is[9].status, intersection::status::e_inside);
    EXPECT_EQ(is[10].status, intersection::status::e_inside);
    EXPECT_EQ(is[11].status, intersection::status::e_outside);
}

TEST(tools, line_intersector_stereo_wire) {

    // tf3 with skewed axis
    vector3 x{1, 0, -1};
    vector3 z{1, 0, 1};
    vector3 t{1, 1, 1};
    const transform3 tf{t, vector::normalize(z), vector::normalize(x)};

    // Create a track
    point3 pos{1, -1, 0};
    vector3 dir{0, 1, 0};
    free_track_parameters trk(pos, 0, dir, -1);

    // Infinite wire with 10 mm radial cell size
    line<> ln{10., std::numeric_limits<scalar>::infinity(), 0u};

    // Test intersect
    ray_line_intersector::intersection_type is =
        ray_line_intersector().intersect(tf, trk, ln);

    EXPECT_EQ(is.status, intersection::status::e_inside);
    EXPECT_FLOAT_EQ(is.path, 2.);
    EXPECT_NEAR(is.p3[0], 1., tolerance);
    EXPECT_NEAR(is.p3[1], 1., tolerance);
    EXPECT_NEAR(is.p3[2], 0., tolerance);
    EXPECT_NEAR(is.p2[0], sqrt(0.5 * 0.5 + 0.5 * 0.5), tolerance);
    EXPECT_NEAR(is.p2[1], -1 * sqrt(0.5 * 0.5 + 0.5 * 0.5), tolerance);
}