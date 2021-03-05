/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/rectangle2.hpp"

#include <climits>
#include <gtest/gtest.h>

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, rectangle2)
{
    using local_type = __plugin::cartesian2;
    using point2 = local_type::point2;

    point2 p2_in = {0.5, -9.};
    point2 p2_edge = {1., 9.3};
    point2 p2_out = {1.5, -9.};

    scalar hx = 1.;
    scalar hy = 9.3;

    rectangle2<> r2 = {hx, hy};

    ASSERT_EQ(r2[0], hx);
    ASSERT_EQ(r2[1], hy);

    ASSERT_TRUE(r2.is_inside<local_type>(p2_in) == intersection_status::e_inside);
    ASSERT_TRUE(r2.is_inside<local_type>(p2_edge) == intersection_status::e_inside);
    ASSERT_TRUE(r2.is_inside<local_type>(p2_out) == intersection_status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside<local_type>(p2_out, {1., 0.5}) == intersection_status::e_inside);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
