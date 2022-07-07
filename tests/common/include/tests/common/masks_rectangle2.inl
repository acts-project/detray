/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>

#include "detray/masks/rectangle2.hpp"

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, rectangle2) {
    using point2 = __plugin::point2<detray::scalar>;

    point2 p2_in = {0.5, -9.};
    point2 p2_edge = {1., 9.3};
    point2 p2_out = {1.5, -9.};

    scalar hx = 1.;
    scalar hy = 9.3;

    rectangle2<> r2{hx, hy, 0u};

    ASSERT_EQ(r2[0], hx);
    ASSERT_EQ(r2[1], hy);

    ASSERT_TRUE(r2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside(p2_out, 1.) == intersection::status::e_inside);
}
