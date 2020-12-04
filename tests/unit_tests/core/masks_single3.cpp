/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "masks/single3.hpp"

#include <climits>

#include <gtest/gtest.h>

using namespace detray;
using namespace plugin;

// This tests the basic function of a rectangle
TEST(mask, single3_0)
{
    point3 p3_in = {0.5, -9., 0.};
    point3 p3_edge = {1., 9.3, 2.};
    point3 p3_out = {1.5, -9.8, 8.};

    scalar h0 = 1.;
    single3<scalar,0> m1_0 = {h0};

    ASSERT_EQ(m1_0[0], h0);

    ASSERT_TRUE(m1_0(p3_in) == intersection_status::e_inside);
    ASSERT_TRUE(m1_0(p3_edge) == intersection_status::e_inside);
    ASSERT_TRUE(m1_0(p3_out) == intersection_status::e_outside);
    // Move outside point inside using a tolerance - take t0 not t1
    ASSERT_TRUE(m1_0(p3_out, 0.6) == intersection_status::e_inside);

}

// This tests the basic function of a rectangle
TEST(mask, single3_1)
{
    point3 p3_in = {0.5, -9., 0.};
    point3 p3_edge = {1., 9.3, 2.};
    point3 p3_out = {1.5, -9.8, 8.};

    scalar h1 = 9.3;
    single3<scalar,1> m1_1 = {h1};

    ASSERT_EQ(m1_1[0], h1);

    ASSERT_TRUE(m1_1(p3_in) == intersection_status::e_inside);
    ASSERT_TRUE(m1_1(p3_edge) == intersection_status::e_inside);
    ASSERT_TRUE(m1_1(p3_out) == intersection_status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_1(p3_out, 0.6) == intersection_status::e_inside);

}

// This tests the basic function of a rectangle
TEST(mask, single3_2)
{
    point3 p3_in = {0.5, -9., 0.};
    point3 p3_edge = {1., 9.3, 2.};
    point3 p3_out = {1.5, -9.8, 8.};

    scalar h2 = 2.;
    single3<scalar,2> m1_2 = {h2};

    ASSERT_EQ(m1_2[0], h2);

    ASSERT_TRUE(m1_2(p3_in) == intersection_status::e_inside);
    ASSERT_TRUE(m1_2(p3_edge) == intersection_status::e_inside);
    ASSERT_TRUE(m1_2(p3_out) == intersection_status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_2(p3_out, 6.1) == intersection_status::e_inside);

}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
