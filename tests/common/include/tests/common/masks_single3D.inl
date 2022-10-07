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

/// This tests the basic functionality of a single value mask (index 0)
TEST(mask, single3_0) {
    using point_t = typename mask<single3D<>>::loc_point_t;

    point_t p3_in = {0.5, -9., 0.};
    point_t p3_edge = {1., 9.3, 2.};
    point_t p3_out = {1.5, -9.8, 8.};

    const scalar h0{1.};
    mask<single3D<>> m1_0{0UL, -h0, h0};

    ASSERT_FLOAT_EQ(m1_0[0], -h0);
    ASSERT_FLOAT_EQ(m1_0[1], h0);

    ASSERT_TRUE(m1_0.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t0 not t1
    ASSERT_TRUE(m1_0.is_inside(p3_out, 0.6) == intersection::status::e_inside);
}

/// This tests the basic functionality of a single value mask (index 1)
TEST(mask, single3_1) {
    using point_t = typename mask<single3D<1>>::loc_point_t;

    point_t p3_in = {0.5, -9., 0.};
    point_t p3_edge = {1., 9.3, 2.};
    point_t p3_out = {1.5, -9.8, 8.};

    const scalar h1{9.3};
    mask<single3D<1>> m1_1{0UL, -h1, h1};

    ASSERT_FLOAT_EQ(m1_1[0], -h1);
    ASSERT_FLOAT_EQ(m1_1[1], h1);

    ASSERT_TRUE(m1_1.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_1.is_inside(p3_out, 0.6) == intersection::status::e_inside);
}

/// This tests the basic functionality of a single value mask (index 2)
TEST(mask, single3_2) {
    using point_t = typename mask<single3D<2>>::loc_point_t;

    point_t p3_in = {0.5, -9., 0.};
    point_t p3_edge = {1., 9.3, 2.};
    point_t p3_out = {1.5, -9.8, 8.};

    const scalar h2{2.};
    mask<single3D<2>> m1_2{0UL, -h2, h2};

    ASSERT_FLOAT_EQ(m1_2[0], -h2);
    ASSERT_FLOAT_EQ(m1_2[1], h2);

    ASSERT_TRUE(m1_2.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_2.is_inside(p3_out, 6.1) == intersection::status::e_inside);
}
