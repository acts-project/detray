/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/line_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tracks/tracks.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <cmath>
#include <limits>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using transform3 = __plugin::transform3<detray::scalar>;
using cartesian = cartesian2<transform3>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using point2 = __plugin::point2<scalar>;
using line_intersector_type = line_intersector<transform3>;
using intersection_t = line_plane_intersection<dindex, transform3>;

constexpr scalar tolerance = 1e-5;
constexpr dindex sf_handle = std::numeric_limits<dindex>::max();

// Test simplest case
TEST(tools, line_intersector_case1) {

    // tf3 with Identity rotation and no translation
    const vector3 x{1, 0, 0};
    const vector3 y{0, 1, 0};
    const vector3 z{0, 0, 1};
    const vector3 t{0, 0, 0};

    const transform3 tf{t, x, y, z};

    // Create a track
    std::vector<free_track_parameters<transform3>> trks;
    trks.emplace_back(point3{1, -1, 0}, 0, vector3{0, 1, 0}, -1);
    trks.emplace_back(point3{-1, -1, 0}, 0, vector3{0, 1, 0}, -1);
    trks.emplace_back(point3{1, -1, 2}, 0, vector3{0, 1, -1}, -1);

    // Infinite wire with 10 mm radial cell size
    const mask<line<>> ln{0UL, 10.f, std::numeric_limits<scalar>::infinity()};

    // Test intersect
    std::vector<intersection_t> is(3);
    is[0] = line_intersector_type()(detail::ray(trks[0]), sf_handle, ln, tf)[0];
    is[1] = line_intersector_type()(detail::ray(trks[1]), sf_handle, ln, tf)[0];
    is[2] = line_intersector_type()(detail::ray(trks[2]), sf_handle, ln, tf)[0];

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_EQ(is[0].path, 1);
    EXPECT_EQ(is[0].p3, point3({1, 0, 0}));
    EXPECT_EQ(is[0].p2, point2({-1, 0}));  // right
    EXPECT_NEAR(is[0].cos_incidence_angle, 0, tolerance);

    EXPECT_EQ(is[1].status, intersection::status::e_inside);
    EXPECT_EQ(is[1].path, 1);
    EXPECT_EQ(is[1].p3, point3({-1, 0, 0}));
    EXPECT_EQ(is[1].p2, point2({1, 0}));  // left
    EXPECT_NEAR(is[1].cos_incidence_angle, 0, tolerance);

    EXPECT_EQ(is[2].status, intersection::status::e_inside);
    EXPECT_NEAR(is[2].path, std::sqrt(2), tolerance);
    EXPECT_NEAR(is[2].p3[0], 1., tolerance);
    EXPECT_NEAR(is[2].p3[1], 0., tolerance);
    EXPECT_NEAR(is[2].p3[2], 1., tolerance);
    EXPECT_NEAR(is[2].p2[0], -1., tolerance);  // right
    EXPECT_NEAR(is[2].p2[1], 1., tolerance);
    EXPECT_NEAR(is[2].cos_incidence_angle, 1. / std::sqrt(2), tolerance);
}

// Test inclined wire
TEST(tools, line_intersector_case2) {
    // tf3 with skewed axis
    const vector3 x{1, 0, -1};
    const vector3 z{1, 0, 1};
    const vector3 t{1, 1, 1};
    const transform3 tf{t, vector::normalize(z), vector::normalize(x)};

    // Create a track
    const point3 pos{1, -1, 0};
    const vector3 dir{0, 1, 0};
    const free_track_parameters<transform3> trk(pos, 0, dir, -1);

    // Infinite wire with 10 mm
    // radial cell size
    const mask<line<>> ln{0UL, 10.f, std::numeric_limits<scalar>::infinity()};

    // Test intersect
    const intersection_t is = line_intersector_type()(
        detail::ray<transform3>(trk), sf_handle, ln, tf)[0];

    EXPECT_EQ(is.status, intersection::status::e_inside);
    EXPECT_NEAR(is.path, 2., tolerance);
    EXPECT_NEAR(is.p3[0], 1., tolerance);
    EXPECT_NEAR(is.p3[1], 1., tolerance);
    EXPECT_NEAR(is.p3[2], 0., tolerance);
    EXPECT_NEAR(is.p2[0], -1. / std::sqrt(2), tolerance);  // right
    EXPECT_NEAR(is.p2[1], -1. / std::sqrt(2), tolerance);
}

TEST(tools, line_intersector_square_scope) {

    // tf3 with Identity rotation and no translation
    const vector3 x{1, 0, 0};
    const vector3 y{0, 1, 0};
    const vector3 z{0, 0, 1};
    const vector3 t{0, 0, 0};

    const transform3 tf{t, x, y, z};

    /// Create a track
    std::vector<free_track_parameters<transform3>> trks;
    trks.emplace_back(point3{2, 0, 0}, 0, vector3{-1, 1, 0}, -1);
    trks.emplace_back(point3{1.9, 0, 0}, 0, vector3{-1, 1, 0}, -1);
    trks.emplace_back(point3{2.1, 0, 0}, 0, vector3{-1, 1, 0}, -1);

    trks.emplace_back(point3{-2, 0, 0}, 0, vector3{1, 1, 0}, -1);
    trks.emplace_back(point3{-1.9, 0, 0}, 0, vector3{1, 1, 0}, -1);
    trks.emplace_back(point3{-2.1, 0, 0}, 0, vector3{1, 1, 0}, -1);

    trks.emplace_back(point3{0, -2, 0}, 0, vector3{-1, 1, 0}, -1);
    trks.emplace_back(point3{0, -1.9, 0}, 0, vector3{-1, 1, 0}, -1);
    trks.emplace_back(point3{0, -2.1, 0}, 0, vector3{-1, 1, 0}, -1);

    trks.emplace_back(point3{0, -2, 0}, 0, vector3{1, 1, 0}, -1);
    trks.emplace_back(point3{0, -1.9, 0}, 0, vector3{1, 1, 0}, -1);
    trks.emplace_back(point3{0, -2.1, 0}, 0, vector3{1, 1, 0}, -1);

    // Infinite wire with 1 mm square cell size
    mask<line<true, line_intersector>, dindex, transform3> ln{
        0UL, 1.f, std::numeric_limits<scalar>::infinity()};

    // Test intersect
    std::vector<intersection_t> is;
    for (const auto& trk : trks) {
        is.push_back(line_intersector_type()(detail::ray<transform3>(trk),
                                             sf_handle, ln, tf, 1e-5)[0]);
    }

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_NEAR(is[0].path, std::sqrt(2), tolerance);
    EXPECT_NEAR(is[0].p3[0], 1, tolerance);
    EXPECT_NEAR(is[0].p3[1], 1, tolerance);
    EXPECT_NEAR(is[0].p3[2], 0, tolerance);
    EXPECT_NEAR(is[0].p2[0], -1 * std::sqrt(2), tolerance);
    EXPECT_NEAR(is[0].p2[1], 0, tolerance);

    EXPECT_EQ(is[1].status, intersection::status::e_inside);
    EXPECT_TRUE(is[1].p2[0] < 0);
    EXPECT_EQ(is[2].status, intersection::status::e_outside);
    EXPECT_TRUE(is[2].p2[0] < 0);

    EXPECT_EQ(is[3].status, intersection::status::e_inside);
    EXPECT_TRUE(is[3].p2[0] > 0);
    EXPECT_EQ(is[4].status, intersection::status::e_inside);
    EXPECT_TRUE(is[4].p2[0] > 0);
    EXPECT_EQ(is[5].status, intersection::status::e_outside);
    EXPECT_TRUE(is[5].p2[0] > 0);

    EXPECT_EQ(is[6].status, intersection::status::e_inside);
    EXPECT_TRUE(is[6].p2[0] > 0);
    EXPECT_EQ(is[7].status, intersection::status::e_inside);
    EXPECT_TRUE(is[7].p2[0] > 0);
    EXPECT_EQ(is[8].status, intersection::status::e_outside);
    EXPECT_TRUE(is[8].p2[0] > 0);

    EXPECT_EQ(is[9].status, intersection::status::e_inside);
    EXPECT_TRUE(is[9].p2[0] < 0);
    EXPECT_EQ(is[10].status, intersection::status::e_inside);
    EXPECT_TRUE(is[10].p2[0] < 0);
    EXPECT_EQ(is[11].status, intersection::status::e_outside);
    EXPECT_TRUE(is[11].p2[0] < 0);
}
