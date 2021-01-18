/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/trapezoid2.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace __plugin;

// This tests the basic function of a trapezoid
TEST(mask, trapezoid2)
{
    using local_type = __plugin::cartesian2;
    using point2 = local_type::point2;

    point2 p2_in = {1., -0.5};
    point2 p2_edge = {2.5, 1.};
    point2 p2_out = {3., 1.5};

    scalar hx_miny = 1.;
    scalar hx_maxy = 3.;
    scalar hy = 2.;

    trapezoid2<scalar> t2 = {hx_miny, hx_maxy, hy};

    ASSERT_EQ(t2[0], hx_miny);
    ASSERT_EQ(t2[1], hx_maxy);
    ASSERT_EQ(t2[2], hy);

    ASSERT_TRUE(t2.is_inside<local_type>(p2_in) == intersection_status::e_inside);
    ASSERT_TRUE(t2.is_inside<local_type>(p2_edge) == intersection_status::e_inside);
    ASSERT_TRUE(t2.is_inside<local_type>(p2_out) == intersection_status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(t2.is_inside<local_type>(p2_out, 1., 0.5) == intersection_status::e_inside);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
