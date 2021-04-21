/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/cylinder3.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, cylinder3)
{
    using local_type = __plugin::transform3;
    using point3 = local_type::point3;

    scalar r = 3.;
    scalar hz = 4.;

    point3 p3_in = {r, 0., -1.};
    point3 p3_edge = {0., r, hz};
    point3 p3_out = {static_cast<scalar>(r / std::sqrt(2.)), static_cast<scalar>(r / std::sqrt(2.)), 4.5};
    point3 p3_off = {1., 1., -9.};

    cylinder3<> c = {r, -hz, hz};

    ASSERT_EQ(c[0], r);
    ASSERT_EQ(c[1], -hz);
    ASSERT_EQ(c[2], hz);

    ASSERT_TRUE(c.is_inside<local_type>(p3_in) == intersection_status::e_inside);
    ASSERT_TRUE(c.is_inside<local_type>(p3_edge) == intersection_status::e_inside);
    ASSERT_TRUE(c.is_inside<local_type>(p3_out) == intersection_status::e_outside);
    ASSERT_TRUE(c.is_inside<local_type>(p3_off) == intersection_status::e_missed);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c.is_inside<local_type>(p3_out, {0., 0.6}) == intersection_status::e_inside);
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
