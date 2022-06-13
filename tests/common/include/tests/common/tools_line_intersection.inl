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
using transform3 = __plugin::transform3<detray::scalar>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using point2 = __plugin::point2<scalar>;
constexpr scalar tolerance = 1e-5;

// This defines the local frame test suite
TEST(tools, line_intersector_straight_wire) {

    // tf3 with Identity rotation and no translation
    vector3 x{1, 0, 0};
    vector3 y{0, 1, 0};
    vector3 z{0, 0, 1};
    vector3 t{0, 0, 0};

    const transform3 tf{t, x, y, z};

    /// Create a track
    point3 pos{1, -1, 0};
    vector3 dir{0, 1, 0};
    free_track_parameters trk(pos, 0, dir, -1);

    // Infinite wire with 10 mm and 0.1 mm radial scope
    std::vector<line<>> lines(2);
    lines[0] = line<>{std::numeric_limits<scalar>::infinity(), 10., 0u};
    lines[1] = line<>{std::numeric_limits<scalar>::infinity(), 0.1, 0u};

    // Test line to planar transform
    const auto new_tf =
        ray_line_intersector().line_to_planar_transform3(tf, trk);

    const vector3 new_x = getter::vector<3>(new_tf.matrix(), 0, 0);
    const vector3 new_y = getter::vector<3>(new_tf.matrix(), 0, 1);
    const vector3 new_z = getter::vector<3>(new_tf.matrix(), 0, 2);
    const vector3 new_tsl = new_tf.translation();

    EXPECT_EQ(new_x, vector3({-1, 0, 0}));
    EXPECT_EQ(new_y, vector3({0, 0, 1}));
    EXPECT_EQ(new_z, vector3({0, 1, 0}));
    EXPECT_EQ(new_tsl, vector3({0, 0, 0}));

    // Test intersect
    std::vector<ray_line_intersector::intersection_type> is(2);
    is[0] = ray_line_intersector().intersect(tf, trk, lines[0]);
    is[1] = ray_line_intersector().intersect(tf, trk, lines[1]);

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_EQ(is[0].p3, point3({1, 0, 0}));
    EXPECT_EQ(is[0].p2, point2({-1, 0}));

    EXPECT_EQ(is[1].status, intersection::status::e_outside);
    EXPECT_EQ(is[1].p3, point3({1, 0, 0}));
    EXPECT_EQ(is[1].p2, point2({-1, 0}));
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

    // Infinite wire with 10 mm radial scope length
    line<> ln{std::numeric_limits<scalar>::infinity(), 10., 0u};

    // Test intersect
    std::vector<ray_line_intersector::intersection_type> is(2);
    is[0] = ray_line_intersector().intersect(tf, trk, ln);

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_NEAR(is[0].p3[0], 1., tolerance);
    EXPECT_NEAR(is[0].p3[1], 1., tolerance);
    EXPECT_NEAR(is[0].p3[2], 0., tolerance);

    EXPECT_NEAR(is[0].p2[0], -0.5 * sqrt(2), tolerance);
    EXPECT_NEAR(is[0].p2[1], 0, tolerance);
}