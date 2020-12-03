/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "masks/cylinder3.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace plugin;

// This tests the basic function of a rectangle
TEST(mask, cylinder3)
{
    scalar r = 3.;
    scalar hz = 4.;

    point3 p3_in = {r, 0., -1.};
    point3 p3_edge = {0., r, hz};
    point3 p3_out = {static_cast<scalar>(r / std::sqrt(2.)), static_cast<scalar>(r / std::sqrt(2.)), 4.5};
    point3 p3_off = {1., 1., -9.};

    cylinder3<scalar> c = {r, hz};

    ASSERT_EQ(c[0], r);
    ASSERT_EQ(c[1], hz);

    ASSERT_TRUE(c(p3_in) == intersection_status::e_inside);
    ASSERT_TRUE(c(p3_edge) == intersection_status::e_inside);
    ASSERT_TRUE(c(p3_out) == intersection_status::e_outside);
    ASSERT_TRUE(c(p3_off) == intersection_status::e_missed);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c(p3_out, 0., 0.6) == intersection_status::e_inside);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
