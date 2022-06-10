/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>

#include "detray/masks/single3.hpp"

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, single3_0) {
    using local_type = __plugin::transform3<detray::scalar>;
    using point3 = __plugin::point3<detray::scalar>;

    point3 p3_in = {0.5, -9., 0.};
    point3 p3_edge = {1., 9.3, 2.};
    point3 p3_out = {1.5, -9.8, 8.};

    scalar h0 = 1.;
    single3<0> m1_0{-h0, h0, 0u};

    ASSERT_EQ(m1_0[0], -h0);
    ASSERT_EQ(m1_0[1], h0);

    ASSERT_TRUE(m1_0.is_inside<local_type>(p3_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside<local_type>(p3_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside<local_type>(p3_out) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t0 not t1
    ASSERT_TRUE(m1_0.is_inside<local_type>(p3_out, 0.6) ==
                intersection::status::e_inside);
}

// This tests the basic function of a rectangle
TEST(mask, single3_1) {
    using local_type = __plugin::transform3<detray::scalar>;
    using point3 = __plugin::point3<detray::scalar>;

    point3 p3_in = {0.5, -9., 0.};
    point3 p3_edge = {1., 9.3, 2.};
    point3 p3_out = {1.5, -9.8, 8.};

    scalar h1 = 9.3;
    single3<1> m1_1 = {-h1, h1, 0u};

    ASSERT_EQ(m1_1[0], -h1);
    ASSERT_EQ(m1_1[1], h1);

    ASSERT_TRUE(m1_1.is_inside<local_type>(p3_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside<local_type>(p3_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside<local_type>(p3_out) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_1.is_inside<local_type>(p3_out, 0.6) ==
                intersection::status::e_inside);
}

// This tests the basic function of a rectangle
TEST(mask, single3_2) {
    using local_type = __plugin::transform3<detray::scalar>;
    using point3 = __plugin::point3<detray::scalar>;

    point3 p3_in = {0.5, -9., 0.};
    point3 p3_edge = {1., 9.3, 2.};
    point3 p3_out = {1.5, -9.8, 8.};

    scalar h2 = 2.;
    single3<2> m1_2 = {-h2, h2, 0u};

    ASSERT_EQ(m1_2[0], -h2);
    ASSERT_EQ(m1_2[1], h2);

    ASSERT_TRUE(m1_2.is_inside<local_type>(p3_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside<local_type>(p3_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside<local_type>(p3_out) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_2.is_inside<local_type>(p3_out, 6.1) ==
                intersection::status::e_inside);
}
