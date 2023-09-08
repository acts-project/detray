/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"
#include "detray/test/types.hpp"

using namespace detray;
using point3_t = test::point3;
using transform3_t = test::transform3;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of a trapezoid
GTEST_TEST(detray_masks, trapezoid2D) {
    using point_t = point3_t;

    point_t p2_in = {1.f, -0.5f, 0.f};
    point_t p2_edge = {2.5f, 1.f, 0.f};
    point_t p2_out = {3.f, 1.5f, 0.f};

    constexpr scalar hx_miny{1.f * unit<scalar>::mm};
    constexpr scalar hx_maxy{3.f * unit<scalar>::mm};
    constexpr scalar hy{2.f * unit<scalar>::mm};
    constexpr scalar divisor{1.f / (2.f * hy)};

    mask<trapezoid2D<>> t2{0u, hx_miny, hx_maxy, hy, divisor};

    ASSERT_NEAR(t2[trapezoid2D<>::e_half_length_0], hx_miny, tol);
    ASSERT_NEAR(t2[trapezoid2D<>::e_half_length_1], hx_maxy, tol);
    ASSERT_NEAR(t2[trapezoid2D<>::e_half_length_2], hy, tol);
    ASSERT_NEAR(t2[trapezoid2D<>::e_divisor], divisor, tol);

    ASSERT_TRUE(t2.is_inside(p2_in));
    ASSERT_TRUE(t2.is_inside(p2_edge));
    ASSERT_FALSE(t2.is_inside(p2_out));
    // Move outside point inside using a tolerance
    ASSERT_TRUE(t2.is_inside(p2_out, 1.));

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = t2.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(hx_maxy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (hx_maxy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}
