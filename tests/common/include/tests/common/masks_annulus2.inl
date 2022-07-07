/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/annulus2.hpp"

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, annulus2) {
    using point2 = __plugin::point2<detray::scalar>;

    scalar minR = 7.2;
    scalar maxR = 12.0;
    scalar minPhi = 0.74195;
    scalar maxPhi = 1.33970;
    point2 offset = {-2., 2.};

    // points in cartesian module frame
    point2 p2_in = {7., 7.};
    point2 p2_out1 = {5., 5.};
    point2 p2_out2 = {10., 3.};
    point2 p2_out3 = {10., 10.};
    point2 p2_out4 = {4., 10.};

    annulus2<> ann2{minR, maxR, minPhi, maxPhi, offset[0], offset[1], 0., 0u};

    ASSERT_EQ(ann2[0], static_cast<scalar>(7.2));
    ASSERT_EQ(ann2[1], static_cast<scalar>(12.0));
    ASSERT_EQ(ann2[2], static_cast<scalar>(0.74195));
    ASSERT_EQ(ann2[3], static_cast<scalar>(1.33970));
    ASSERT_EQ(ann2[4], static_cast<scalar>(-2.0));
    ASSERT_EQ(ann2[5], static_cast<scalar>(2.0));
    ASSERT_EQ(ann2[6], static_cast<scalar>(0.));

    ASSERT_TRUE(ann2.is_inside(point2(p2_in + offset)) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ann2.is_inside(point2(p2_out1 + offset)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(point2(p2_out2 + offset)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(point2(p2_out3 + offset)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(point2(p2_out4 + offset)) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(ann2.is_inside(point2(p2_out1 + offset), 1.3) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ann2.is_inside(point2(p2_out4 + offset), 0.07) ==
                intersection::status::e_inside);
}
