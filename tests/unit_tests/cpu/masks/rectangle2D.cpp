/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"
#include "detray/test/types.hpp"

// GTest include
#include <gtest/gtest.h>

using namespace detray;
using point3_t = test::point3;

constexpr scalar tol{1e-7f};

constexpr scalar hx{1.f * unit<scalar>::mm};
constexpr scalar hy{9.3f * unit<scalar>::mm};
constexpr scalar hz{0.5f * unit<scalar>::mm};

/// This tests the basic functionality of a rectangle
GTEST_TEST(detray_masks, rectangle2D) {
    using point_t = point3_t;

    point_t p2_in = {0.5f, -9.f, 0.f};
    point_t p2_edge = {1.f, 9.3f, 0.f};
    point_t p2_out = {1.5f, -9.f, 0.f};

    mask<rectangle2D<>> r2{0u, hx, hy};

    ASSERT_NEAR(r2[rectangle2D<>::e_half_x], hx, tol);
    ASSERT_NEAR(r2[rectangle2D<>::e_half_y], hy, tol);

    ASSERT_TRUE(r2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside(p2_out, 1.f) == intersection::status::e_inside);

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = r2.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);

    const auto centroid = r2.centroid();
    ASSERT_NEAR(centroid[0], 0.f, tol);
    ASSERT_NEAR(centroid[1], 0.f, tol);
    ASSERT_NEAR(centroid[2], 0.f, tol);
}

/// This tests the basic functionality of a cuboid3D
GTEST_TEST(detray_masks, cuboid3D) {
    using point_t = point3_t;

    point_t p2_in = {0.5f, 8.0f, -0.4f};
    point_t p2_edge = {1.f, 9.3f, 0.5f};
    point_t p2_out = {1.5f, -9.f, 0.55f};

    mask<cuboid3D<>> c3{0u, -hx, -hy, -hz, hx, hy, hz};

    ASSERT_NEAR(c3[cuboid3D<>::e_min_x], -hx, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_min_y], -hy, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_min_z], -hz, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_max_x], hx, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_max_y], hy, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_max_z], hz, tol);

    ASSERT_TRUE(c3.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(c3.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(c3.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c3.is_inside(p2_out, 1.f) == intersection::status::e_inside);

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = c3.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -(hz + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], (hz + envelope), tol);
}
