/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/annulus2.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, annulus2)
{
    using polar     = __plugin::polar2;
    using cartesian = __plugin::cartesian2;
    using point2_pl = polar::point2;
    using point2_c  = cartesian::point2;

    scalar minR = 7.2;
    scalar maxR = 12.0;
    scalar minPhi = 0.74195;
    scalar maxPhi = 1.33970;
    point2_c offset = {-2., 2.};

    // points in cartesian module frame
    point2_c p2_in   = {7., 7.};
    point2_c p2_out1 = {5., 5.};
    point2_c p2_out2 = {10., 3.};
    point2_c p2_out3 = {10., 10.};
    point2_c p2_out4 = {4., 10.};

    auto toStripFrame = [&](const point2_pl& xy) -> point2_pl {
        auto shifted = xy + offset;
        scalar r   = getter::perp(shifted);
        scalar phi = getter::phi(shifted);
        return point2_pl{r, phi};
    };

    annulus2<scalar> ann2 = {minR, maxR, minPhi, maxPhi, offset[0], offset[1]};

    ASSERT_EQ(ann2[0], static_cast<scalar>(7.2));
    ASSERT_EQ(ann2[1], static_cast<scalar>(12.0));
    ASSERT_EQ(ann2[2], static_cast<scalar>(0.74195));
    ASSERT_EQ(ann2[3], static_cast<scalar>(1.33970));
    ASSERT_EQ(ann2[4], static_cast<scalar>(-2.0));
    ASSERT_EQ(ann2[5], static_cast<scalar>(2.0));
    ASSERT_EQ(ann2[6], static_cast<scalar>(0.));

    ASSERT_TRUE(ann2.is_inside<cartesian>(p2_in + offset) == intersection_status::e_inside);
    ASSERT_TRUE(ann2.is_inside<cartesian>(p2_out1 + offset) == intersection_status::e_outside);
    ASSERT_TRUE(ann2.is_inside<cartesian>(p2_out2 + offset) == intersection_status::e_outside);
    ASSERT_TRUE(ann2.is_inside<cartesian>(p2_out3 + offset) == intersection_status::e_outside);
    ASSERT_TRUE(ann2.is_inside<cartesian>(p2_out4 + offset) == intersection_status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(ann2.is_inside<cartesian>(p2_out1 + offset, {1.3, 0.}) == intersection_status::e_inside);
    ASSERT_TRUE(ann2.is_inside<cartesian>(p2_out4 + offset, {0., 0.07}) == intersection_status::e_inside);

    ASSERT_TRUE(ann2.is_inside<polar>(toStripFrame(p2_in)) == intersection_status::e_inside);
    ASSERT_TRUE(ann2.is_inside<polar>(toStripFrame(p2_out1)) == intersection_status::e_outside);
    ASSERT_TRUE(ann2.is_inside<polar>(toStripFrame(p2_out2)) == intersection_status::e_outside);
    ASSERT_TRUE(ann2.is_inside<polar>(toStripFrame(p2_out3)) == intersection_status::e_outside);
    ASSERT_TRUE(ann2.is_inside<polar>(toStripFrame(p2_out4)) == intersection_status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(ann2.is_inside<polar>(toStripFrame(p2_out1), {1.3, 0.}) == intersection_status::e_inside);
    ASSERT_TRUE(ann2.is_inside<polar>(toStripFrame(p2_out4), {0., 0.07}) == intersection_status::e_inside);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
