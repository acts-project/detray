/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

/// This tests the basic functionality of a trapezoid
TEST(mask, trapezoid2D) {
    using point_t = typename mask<trapezoid2D<>>::loc_point_t;

    point_t p2_in = {1., -0.5};
    point_t p2_edge = {2.5, 1.};
    point_t p2_out = {3., 1.5};

    scalar hx_miny = 1.;
    scalar hx_maxy = 3.;
    scalar hy = 2.;
    scalar divisor = 1. / (2. * hy);

    mask<trapezoid2D<>> t2{0UL, hx_miny, hx_maxy, hy, divisor};

    ASSERT_EQ(t2[0], hx_miny);
    ASSERT_EQ(t2[1], hx_maxy);
    ASSERT_EQ(t2[2], hy);
    ASSERT_EQ(t2[3], divisor);

    ASSERT_TRUE(t2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(t2.is_inside(p2_out, 1.) == intersection::status::e_inside);
}
