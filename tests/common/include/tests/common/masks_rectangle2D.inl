/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>

#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

/// This tests the basic functionality of a rectangle
TEST(mask, rectangle2D) {
    using point_t = typename mask<rectangle2D<>>::loc_point_t;

    point_t p2_in = {0.5, -9.};
    point_t p2_edge = {1., 9.3};
    point_t p2_out = {1.5, -9.};

    scalar hx = 1.;
    scalar hy = 9.3;

    mask<rectangle2D<>> r2{0UL, hx, hy};

    ASSERT_FLOAT_EQ(r2[0], hx);
    ASSERT_FLOAT_EQ(r2[1], hy);

    ASSERT_TRUE(r2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside(p2_out, 1.) == intersection::status::e_inside);
}
