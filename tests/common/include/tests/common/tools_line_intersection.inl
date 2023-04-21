/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/geometry/surface.hpp"
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
using intersection_t = intersection2D_point<surface<>, transform3>;
using line_intersector_type = line_intersector<intersection_t>;

constexpr scalar tol{1e-5f};

// Test simplest case
TEST(tools, line_intersector_case1) {

    // tf3 with Identity rotation and no translation
    const transform3 tf{};

    // Create a track
    std::vector<free_track_parameters<transform3>> trks;
    trks.emplace_back(point3{1.f, -1.f, 0.f}, 0.f, vector3{0.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{-1.f, -1.f, 0.f}, 0.f, vector3{0.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{1.f, -1.f, 2.f}, 0.f, vector3{0.f, 1.f, -1.f},
                      -1.f);

    // Infinite wire with 10 mm radial cell size
    const mask<line<>> ln{0u, 10.f, std::numeric_limits<scalar>::infinity()};

    // Test intersect
    std::vector<intersection_t> is(3u);
    is[0] = line_intersector_type()(detail::ray(trks[0]), surface<>{}, ln, tf);
    is[1] = line_intersector_type()(detail::ray(trks[1]), surface<>{}, ln, tf);
    is[2] = line_intersector_type()(detail::ray(trks[2]), surface<>{}, ln, tf);

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_EQ(is[0].path, 1.f);
    EXPECT_EQ(is[0].p3, point3({1.f, 0.f, 0.f}));
    EXPECT_EQ(is[0].p2, point2({-1.f, 0.f}));  // right
    EXPECT_NEAR(is[0].cos_incidence_angle, 0.f, tol);

    EXPECT_EQ(is[1].status, intersection::status::e_inside);
    EXPECT_EQ(is[1].path, 1.f);
    EXPECT_EQ(is[1].p3, point3({-1.f, 0.f, 0.f}));
    EXPECT_EQ(is[1].p2, point2({1.f, 0.f}));  // left
    EXPECT_NEAR(is[1].cos_incidence_angle, 0.f, tol);

    EXPECT_EQ(is[2].status, intersection::status::e_inside);
    EXPECT_NEAR(is[2].path, constant<scalar>::sqrt2, tol);
    EXPECT_NEAR(is[2].p3[0], 1.f, tol);
    EXPECT_NEAR(is[2].p3[1], 0.f, tol);
    EXPECT_NEAR(is[2].p3[2], 1.f, tol);
    EXPECT_NEAR(is[2].p2[0], -1.f, tol);  // right
    EXPECT_NEAR(is[2].p2[1], 1.f, tol);
    EXPECT_NEAR(is[2].cos_incidence_angle, constant<scalar>::inv_sqrt2, tol);
}

// Test inclined wire
TEST(tools, line_intersector_case2) {
    // tf3 with skewed axis
    const vector3 x{1.f, 0.f, -1.f};
    const vector3 z{1.f, 0.f, 1.f};
    const vector3 t{1.f, 1.f, 1.f};
    const transform3 tf{t, vector::normalize(z), vector::normalize(x)};

    // Create a track
    const point3 pos{1.f, -1.f, 0.f};
    const vector3 dir{0.f, 1.f, 0.f};
    const free_track_parameters<transform3> trk(pos, 0.f, dir, -1.f);

    // Infinite wire with 10 mm
    // radial cell size
    const mask<line<>> ln{0u, 10.f, std::numeric_limits<scalar>::infinity()};

    // Test intersect
    const intersection_t is = line_intersector_type()(
        detail::ray<transform3>(trk), surface<>{}, ln, tf);

    EXPECT_EQ(is.status, intersection::status::e_inside);
    EXPECT_NEAR(is.path, 2.f, tol);
    EXPECT_NEAR(is.p3[0], 1.f, tol);
    EXPECT_NEAR(is.p3[1], 1.f, tol);
    EXPECT_NEAR(is.p3[2], 0.f, tol);
    EXPECT_NEAR(is.p2[0], -constant<scalar>::inv_sqrt2, tol);  // right
    EXPECT_NEAR(is.p2[1], -constant<scalar>::inv_sqrt2, tol);
}

TEST(tools, line_intersector_square_scope) {

    // tf3 with Identity rotation and no translation
    const transform3 tf{};

    /// Create a track
    std::vector<free_track_parameters<transform3>> trks;
    trks.emplace_back(point3{2.f, 0.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{1.9f, 0.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{2.1f, 0.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                      -1.f);

    trks.emplace_back(point3{-2.f, 0.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{-1.9f, 0.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{-2.1f, 0.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f},
                      -1.f);

    trks.emplace_back(point3{0.f, -2.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{0.f, -1.9f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{0.f, -2.1f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                      -1.f);

    trks.emplace_back(point3{0.f, -2.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{0.f, -1.9f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f},
                      -1.f);
    trks.emplace_back(point3{0.f, -2.1f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f},
                      -1.f);

    // Infinite wire with 1 mm square cell size
    mask<line<true, line_intersector>, std::uint_least16_t, transform3> ln{
        0u, 1.f, std::numeric_limits<scalar>::infinity()};

    // Test intersect
    std::vector<intersection_t> is;
    for (const auto& trk : trks) {
        is.push_back(line_intersector_type()(detail::ray<transform3>(trk),
                                             surface<>{}, ln, tf, 1e-5f));
    }

    EXPECT_EQ(is[0].status, intersection::status::e_inside);
    EXPECT_NEAR(is[0].path, constant<scalar>::sqrt2, tol);
    EXPECT_NEAR(is[0].p3[0], 1.f, tol);
    EXPECT_NEAR(is[0].p3[1], 1.f, tol);
    EXPECT_NEAR(is[0].p3[2], 0.f, tol);
    EXPECT_NEAR(is[0].p2[0], -constant<scalar>::sqrt2, tol);
    EXPECT_NEAR(is[0].p2[1], 0.f, tol);

    EXPECT_EQ(is[1].status, intersection::status::e_inside);
    EXPECT_TRUE(std::signbit(is[1].p2[0]));
    EXPECT_EQ(is[2].status, intersection::status::e_outside);
    EXPECT_TRUE(std::signbit(is[2].p2[0]));

    EXPECT_EQ(is[3].status, intersection::status::e_inside);
    EXPECT_FALSE(std::signbit(is[3].p2[0]));
    EXPECT_EQ(is[4].status, intersection::status::e_inside);
    EXPECT_FALSE(std::signbit(is[4].p2[0]));
    EXPECT_EQ(is[5].status, intersection::status::e_outside);
    EXPECT_FALSE(std::signbit(is[5].p2[0]));

    EXPECT_EQ(is[6].status, intersection::status::e_inside);
    EXPECT_FALSE(std::signbit(is[6].p2[0]));
    EXPECT_EQ(is[7].status, intersection::status::e_inside);
    EXPECT_FALSE(std::signbit(is[7].p2[0]));
    EXPECT_EQ(is[8].status, intersection::status::e_outside);
    EXPECT_FALSE(std::signbit(is[8].p2[0]));

    EXPECT_EQ(is[9].status, intersection::status::e_inside);
    EXPECT_TRUE(std::signbit(is[9].p2[0]));
    EXPECT_EQ(is[10].status, intersection::status::e_inside);
    EXPECT_TRUE(std::signbit(is[10].p2[0]));
    EXPECT_EQ(is[11].status, intersection::status::e_outside);
    EXPECT_TRUE(std::signbit(is[11].p2[0]));
}
